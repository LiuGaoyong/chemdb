#!/usr/bin/env python3
import os, hashlib
import time, uuid
import dpdata, json
import sqlalchemy
import pandas as pd
from typing import List
from pathlib import Path
from sqlalchemy_utils import database_exists, create_database
from multiprocessing.dummy import Pool

from chemdb.utils import get_all_dirs as all_dirs
from chemdb.utils import FileContents


class ChemDB():

    TB_INFO = {
        "main": { 'id'              : sqlalchemy.types.Text(),             
                'software'          : sqlalchemy.types.Text(),
                'files_hash_json'   : sqlalchemy.types.Text(),
                'type_of_elements'  : sqlalchemy.types.Text(),
                'type_of_systems'   : sqlalchemy.types.Text(),
                'num_of_point'      : sqlalchemy.types.Integer() },
        "files": {'id'              : sqlalchemy.types.Text(),
                'binary_contents'   : sqlalchemy.types.LargeBinary() }
    }
    TB_NMAES = set(TB_INFO.keys())


    def __init__(self, db_url:str=None, parallel:bool=False):
        self.parallel = parallel
        # setup SQL engine
        if db_url is None: db_url = "sqlite:///:memory:"
        if not database_exists(db_url):
            create_database(db_url)
        self.engine = sqlalchemy.create_engine(db_url)#,
        #    echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
        #    pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
        #    pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
        #)
        self.inspct = sqlalchemy.inspect(self.engine)
        self.is_empty  = self.__is_empty()
        self.is_chemdb = self.__url_is_for_chemdb()
        if not (self.is_empty or self.is_chemdb):
            raise KeyError("{} isn't a chemdb-typed url".format(db_url))

    def __is_empty(self):
        """判断是不是空数据库"""
        return len(self.inspct.get_table_names()) == 0

    def __url_is_for_chemdb(self):
        """判断是不是chemdb类型的数据库"""
        result = True
        table_names = self.inspct.get_table_names()
        if len(table_names) == 0: 
            # new database
            return True
        if not set(table_names) == set(self.TB_NMAES):
            return False
        md = sqlalchemy.MetaData()
        for key in self.TB_INFO.keys():
            tb  = sqlalchemy.Table(key,  md, 
                autoload=True, autoload_with=self.engine)
            result = (result and set(tb.c.keys()) == \
                set(self.TB_INFO[key].keys()))
            print(tb.primary_key)
        return result

    @property
    def abstract(self):
        if self.is_empty:
            print("main:\t{}\nfiles:\t{}".format(0, 0))
        if self.is_chemdb:

            print("main:\t{}\nfiles:\t{}".format(0, 0))

    def append_data(self, p:Path):
        # setup data
        try:
            assert os.path.isdir(p)
        except AssertionError as e:
            raise KeyError("'{}' don't exists or isn't a dir.".format(p))
        # generate dirs list
        dirs_list = all_dirs(p)
        if len(dirs_list) > 200: dirs_list = dirs_list[:200] #for test 
        # parse dirs_list ---> pd.DataFrame
        data = self.__analyze_data(dirs_list)
        # pandas to sql
        for key in self.TB_NMAES:
            data[key].to_sql(key, self.engine, 
                dtype=self.TB_INFO[key], index=False, 
                if_exists="replace", chunksize=100 )


    def __analyze_data(self, data_list:List[Path]) -> dict:
        data_main = []
        data_files = []
        def process_vasp(i):
            outcar = os.path.join(i, "OUTCAR")
            try:
                outcar_dpdata = dpdata.LabeledSystem(outcar, fmt='vasp/outcar')
            except Exception as e:
                # cannot be parsed, return None
                return None
            v_num_of_point = len(outcar_dpdata)
            if v_num_of_point == 0:
                # parsed successfully, but there's not available data, return None
                return None
            _tmp = outcar_dpdata[0].data
            v_type_of_ele = '-'.join(sorted(_tmp['atom_names']))
            v_type_of_sys = sorted(list(zip(_tmp['atom_names'], 
                    [str(i) for i in _tmp['atom_numbs']])))
            v_type_of_sys = '-'.join(['_'.join(i) for i in v_type_of_sys])
            # generate v_hash_json
            outcar = FileContents(outcar)
            outcar_md5 = outcar.md5
            data_files.append([outcar_md5, outcar.to_compress_binary()])
            v_hash_json = {
                'OUTCAR':   outcar_md5,
                'vasprun':  None    }
            _tmp = os.path.join(i, "vasprun.xml")
            if os.path.exists(_tmp):
                vasprun = FileContents(_tmp)
                vasprun_md5 = vasprun.md5
                v_hash_json['vasprun'] = vasprun_md5
            v_hash_json = json.dumps(v_hash_json)
            return v_hash_json, v_type_of_ele, v_type_of_sys, v_num_of_point

        def process(i):
                if os.path.exists(os.path.join(i, "OUTCAR")):
                    v_software = "VASP"
                    values = process_vasp(i)
                    if values is None:
                        return None
                    #assert len(values) == 4
                    v_hash_json,    v_type_of_ele  = values[0:2]
                    v_type_of_sys,  v_num_of_point = values[2:4]
                    v_id   = str(uuid.uuid5(uuid.NAMESPACE_DNS, v_hash_json))
                    data_main.append([
                        v_id,           v_software,
                        v_hash_json,    v_type_of_ele,
                        v_type_of_sys,  v_num_of_point  ]) 
                elif os.path.exists(os.path.join(i, "type.raw")):
                    v_software = "deepmd"
                    
        if self.parallel:
            pool = Pool()
            pool.map(process, data_list)
            pool.close()
            pool.join()
        else:
            for i in data_list:
                process(i)
        data_main = pd.DataFrame(data_main, columns=[
            'id',               'software',
            'files_hash_json',  'type_of_elements',
            'type_of_systems',  'num_of_point'           ])
        data_main.drop_duplicates(['id'])
        data_files = pd.DataFrame(data_files, columns=[
            'id',        'binary_contents'])
        data_files.drop_duplicates(['id'])
        return {"main": data_main, "files": data_files}

    

    

    

            


if __name__ == "__main__":
    engine_url = 'sqlite:///./test.db'
    data = "/home/lgy/cal/pd-mo-cluster-cal/"
    #"/home/lgy/cal/pd_mo_o-dpmd_pot/",
    
    test = ChemDB(engine_url, parallel=True)
    test.append_data(data)


