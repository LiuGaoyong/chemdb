#!/usr/bin/env python3
import os, hashlib
import time, uuid
import dpdata, json
import pandas as pd
from typing import List, Union
from pathlib import Path
from sqlalchemy import create_engine
from sqlalchemy.engine.base import Engine
from multiprocessing.dummy import Pool

from chemdb.utils import get_all_dirs as all_dirs
from chemdb.utils import FileContents


class ChemDB():
    def __init__(self, data:List[Path], engine:Union[str, Engine], parallel:bool=False):
        self.parallel = parallel
        # setup SQL engine
        if not isinstance(engine, Engine):
            engine = create_engine(engine)
        self.engine = engine
        # setup data
        try:
            assert len(data) != 0
        except AssertionError as e:
            raise KeyError("the length of data is 0.")
        for i in range(len(data)): 
            try:
                assert os.path.isdir(data[i])
            except AssertionError as e:
                raise KeyError("the {}th data: '{}' ".format(i+1, data[i]) + \
                            "don't exists or isn't a dir.")
        # generate dirs list
        dirs_list = []
        if self.parallel:
            def process(i):
                dirs_list.extend(all_dirs(i))
            pool = Pool()
            pool.map(process, data)
            pool.close()
            pool.join()
        else:
            for i in data:
                dirs_list.extend(all_dirs(i))
        if len(dirs_list) > 200: dirs_list = dirs_list[:200]
        self.df_main, self.df_files = self.__analyze_data(dirs_list)

    def __analyze_data(self, data_list:List[Path]):
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
                v_id   = uuid.uuid5(uuid.NAMESPACE_DNS, str(time.time()))
                if os.path.exists(os.path.join(i, "OUTCAR")):
                    v_software = "VASP"
                    values = process_vasp(i)
                    if values is None:
                        return None
                    #assert len(values) == 4
                    v_hash_json,    v_type_of_ele  = values[0:2]
                    v_type_of_sys,  v_num_of_point = values[2:4]
                    data_main.append([
                        v_id,           v_software,
                        v_hash_json,    v_type_of_ele,
                        v_type_of_sys,  v_num_of_point  ])  
                    print("{:10} {}".format(v_software, i)) 
                elif os.path.exists(os.path.join(i, "type.raw")):
                    v_software = "deepmd"
                    print("{:10} {}".format(v_software, i))
                    
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
        data_files = pd.DataFrame(data_files, columns=[
            'hash',     'binary_contents'])
        return data_main, data_files

    def show_database_abstract(self):
        pass

    

            


if __name__ == "__main__":
    engine = 'sqlite:///./test.db'
    data = [
       #"/home/lgy/cal/pd_mo_o-dpmd_pot/",
       "/home/lgy/cal/pd-mo-cluster-cal/",
    ]
    test = ChemDB(data, engine, parallel=True)


