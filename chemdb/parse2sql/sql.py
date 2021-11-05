#!/usr/bin/env python3
from sqlalchemy import create_engine
from sqlalchemy import Column
from sqlalchemy import Integer
from sqlalchemy import String
from urllib.parse import quote_plus


class ChemCOD():
    def __init__(self, engine):
        pass


if __name__ == "__main__":
    passwd = quote_plus("chem123@G505")
    engine = create_engine( "mysql+pymysql://chem:{}".format(passwd) +
                            "@114.212.167.205:13306/chem_cod?charset=utf8",
            echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
            pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
            pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
            )

    from sqlalchemy.ext.declarative import declarative_base
    Base = declarative_base()
    class Users(Base):
        __tablename__ = "users"

        id = Column(Integer, primary_key=True)
        name = Column(String(64), unique=True)
        email = Column(String(64))

        def __init__(self, name, email):
            self.name = name
            self.email = email 
    Base.metadata.create_all(engine) 
    
    
    