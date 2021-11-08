#!/usr/bin/env python3
import pickle, zlib, os
from ase import Atoms
from bson.binary import Binary
from sqlalchemy import create_engine, Column
from sqlalchemy import Integer, LargeBinary, String, Text
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
from urllib.parse import quote_plus

class Many_SDF():
    def __init__(self, fname:str) -> None:
        assert os.path.exists(fname)
        self.fname = os.path.abspath(fname)

    @property
    def generation(self):
        with open(self.fname) as f:
            tmp = ""
            while True:
                line = f.readline()
                tmp = tmp + line
                if "$$$$" in line:
                    yield tmp
                    tmp = ""
    
    @property
    def generation_to_ase(self):
        for i in self.generation:
            try:
                with open("/tmp/tmp.sdf", "w") as f:
                    f.write(i)
                atoms = read_sdf("/tmp/tmp.sdf")
                yield i, atoms
            except KeyError as e:
                #print("{:20} {}".format("KeyError", e))
                pass
            except ValueError as e:
                #print("{:20} {}".format("ValueError", e))
                pass
            except Exception as e:
                #print("{:20} {}".format("General Exception", e))
                pass



passwd = quote_plus("chem123@G505")
test_engine = create_engine("sqlite:///:memory:", echo=True)
chem_cache_engine = create_engine( "mysql+pymysql://chem:{}".format(passwd) +
                        "@114.212.167.205:13306/chem_cache?charset=utf8",
    echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
    pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
    pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
    )
Base = declarative_base()
class NCI_Open(Base):
    __tablename__ = "NCI_Open"
    id       = Column(Integer, primary_key=True, index=True)
    formula  = Column(String, index=True)
    smiles   = Column(String, index=True)

    contents = Column(LargeBinary, unique=True)
    def __init__(self, contents:str, atoms:Atoms):
        
        contents = pickle.dumps(contents)
        contents = zlib.compress(contents)
        self.contents = Binary(contents)
Base.metadata.create_all(chem_cache_engine) 
print("fasdf")

if True:
    sdf = "/home/data/sync/hp-ChemDB_Files/sdf/" + \
        "NCI-Open_2012-05-01.sdf"
    sdf = Many_SDF(sdf)
    sess = sessionmaker(bind=chem_cache_engine)
    sess = sess()
    print("fasdfassda")
    num = 0
    print("fasdfasd")

    for contents in sdf.generation:
        print("fasdf")
        num += 1
        
        sess.add(NCI_Open(contents))
        if num > 10:
            break
    sess.commit()
    
    

    