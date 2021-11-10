#!/usr/bin/env python3
import pickle, zlib, os
from bson.binary import Binary
from uuid import NAMESPACE_DNS, uuid5
from ase.io.sdf import read_sdf
from sqlalchemy_utils import database_exists, create_database
from sqlalchemy import Table, Column, create_engine, MetaData
from sqlalchemy import Integer, Float, LargeBinary, UnicodeText, Boolean

class ChemCache():
    metadata = MetaData()
    def __init__(self, url:str):
        if not "{}" in url: raise KeyError("init fail!")
        self.db_url = url.format("chem_4web")
        if not database_exists(self.db_url):
            create_database(self.db_url)
        self.engine = create_engine(self.db_url,
            echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
            pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
            pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
        )
        self.__init_table()
        
    
    def __init_table(self):
        tablenames = [
            "molcule",
            "crystal",
            "alloy",
            "nanocages",
            "molecular_sieve",
            "metal_oxide"]
        self.table = dict()
        for name in tablenames:
            if name not in self.engine.table_names():
                self.table[name] = Table(name, self.metadata,
                    Column("id",      Integer, primary_key=True),
                    Column("img",     UnicodeText, index=True, unique=True),
                    Column("name",    UnicodeText, index=True),
                    Column("formula", UnicodeText, index=True),
                    Column("weight",  UnicodeText),
                    Column("cas",     UnicodeText, index=True),
                    Column("inchi",   UnicodeText, index=True),
                    Column("inchikey",  UnicodeText, index=True),
                    Column("downloads", UnicodeText),
                    Column("elements",  UnicodeText, index=True, nullable=False),
                    Column("typed", UnicodeText, index=True),
                    Column("spacegroup", UnicodeText, index=True),
                    Column("a", Float),
                    Column("b", Float),
                    Column("c", Float),
                    Column("alpha", Float),
                    Column("beta",  Float),
                    Column("gamma", Float),
                    Column("energy", Float),
                    Column("efermi", Float),
                    Column("s_band", Float),
                    Column("p_band", Float),
                    Column("d_band", Float),
                    Column("attach", UnicodeText, index=True) )
                self.metadata.create_all(self.engine)
            else:
                self.metadata = MetaData(self.engine)
                self.table[name] = Table(name, self.metadata, autoload=True)

    

if __name__ == "__main__":
    from urllib.parse import quote_plus
    passwd = quote_plus("chem123@G505")
    db_url = "mysql+pymysql://chem:{}".format(passwd) + \
             "@114.212.167.205:13306/{}?charset=utf8"
    test = ChemCache(db_url)
    print(test.table)
    print()

    








