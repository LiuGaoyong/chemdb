#!/usr/bin/env python3
import pickle, zlib, os
from bson.binary import Binary
from uuid import NAMESPACE_DNS, uuid5
from ase.io.sdf import read_sdf
from ase.io import read
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import Table, Column, create_engine, MetaData, inspect
from sqlalchemy import Integer, LargeBinary, UnicodeText, Boolean
from sqlalchemy.sql import select, insert
from pymongo import MongoClient

class ChemCache():
    metadata = MetaData()
    def __init__(self, url:str):
        if not "{}" in url: raise KeyError("init fail!")
        self.db_url = url.format("chem_cache")
        if not database_exists(self.db_url):
            create_database(self.db_url)
        self.engine = create_engine(self.db_url,
            echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
            pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
            pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
        )
        self.inspect = inspect(self.engine)
        self.__init_table()
        
    
    def __init_table(self):
        tablenames = [
            "nci_open",
            "roadmap",
            "cod",
            "amcsd",
            "alloy_slab",
            "alloy_ads"]
        self.table = dict()
        for name in tablenames:
            if name not in self.inspect.get_table_names():
                self.table[name] = Table(name, self.metadata,
                    Column("id", Integer, primary_key=True),
                    Column("uuid", UnicodeText, unique=True, index=True),
                    Column("contents", LargeBinary, unique=True),
                    Column("can_by_ase", Boolean, index=True)  )
                self.metadata.create_all(self.engine)
            else:
                self.metadata = MetaData(self.engine)
                self.table[name] = Table(name, self.metadata, autoload=True)

    def _scan_sdf2sql(self, fname, tablename):
        assert os.path.exists(fname)
        assert tablename in self.table.keys()
        conn = self.engine.connect()
        with open(fname) as f:
            tmp, clock, insert_list = '', 0, []
            while True:
                line  = f.readline()
                tmp += line
                if '$$$$' in line:
                    clock += 1
                    # 此时tmp恰好存储了一个Atoms的sdf信息
                    tmp_sdf_name = "/tmp/tmp.sdf"
                    with open(tmp_sdf_name, 'w') as f_sdf:
                        f_sdf.write(tmp)
                    try:
                        can_by_ase = True
                        atoms = read_sdf(tmp_sdf_name)
                        print(atoms)
                    except Exception as e:
                        can_by_ase = False
                    uuid = uuid5(NAMESPACE_DNS, tmp)
                    contents = Binary(zlib.compress(pickle.dumps(tmp)))
                    insert_list.append({"uuid": uuid, 
                        "contents":contents, "can_by_ase":can_by_ase})
                    if len(insert_list) == 100:
                        conn.execute(self.table[tablename].insert(), insert_list)
                        insert_list = list()
                    tmp = ''
                if line == '':
                    conn.execute(self.table[tablename].insert(), insert_list)
                    break

    def _scan_cif2sql(self, tablename):
        mongo = MongoClient("mongodb://root:LGY1990907u-mongo@114.212.167.205:27017")
        col = mongo["cache"][tablename]
        conn = self.engine.connect()
        clock, insert_list = 0, list()
        for i in col.find():
            contents = i['result']
            contents = zlib.decompress(contents)
            contents = pickle.loads(contents)
            uuid = uuid5(NAMESPACE_DNS, contents)
            with open("/tmp/CIF.cif", "w") as f:
                f.write(contents)
            try:
                can_by_ase = True
                atoms = read("/tmp/CIF.cif", format="cif")
                print(atoms)
                clock += 1
            except Exception as e:
                can_by_ase = True
                continue
            try:
                conn.execute(self.table[tablename].insert(), {"uuid": str(uuid), 
                    "contents":i['result'], "can_by_ase":can_by_ase})
            except Exception as e:
                print(e)
                continue
            if clock > 10:
                break
        conn.execute(self.table[tablename].insert(), insert_list)

    def get_one(self, tablename, key):
        table = self.table[tablename]
        conn = self.engine.connect()
        if isinstance(key, int):
            rows = conn.execute(select([table]).where(table.c.id == key))
        elif isinstance(key, str):
            rows = conn.execute(select([table]).where(table.c.uuid == key))
        else:
            raise TypeError("key must be int|str typed")
        result = [i for i in rows]
        assert len(result) == 1
        return result[0]

    def get_all_id_can_by_ase(self, tablename):
        tb = self.table[tablename]
        conn = self.engine.connect()
        rows = conn.execute(select(tb.c.id).where(tb.c.can_by_ase==1))
        return [i for i, in rows.fetchall()]

    




            


if __name__ == "__main__":
    from urllib.parse import quote_plus
    passwd = quote_plus("chem123@G505")
    db_url = "mysql+pymysql://chem:{}".format(passwd) + \
             "@114.212.167.205:13306/{}?charset=utf8"
    test = ChemCache(db_url)
    all_can_by_ase = test.get_all_id_can_by_ase('nci_open')
    print(all_can_by_ase)
    for i in all_can_by_ase[:3]:
        print(test.get_one("nci_open", i))











