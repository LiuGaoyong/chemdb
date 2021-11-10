#!/usr/bin/env python3
import pickle, zlib, os
from urllib.parse import quote_plus
from ase import Atoms
from ase.io import read, write
from ase.io.sdf import read_sdf
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions
from bson.binary import Binary
from sqlalchemy import create_engine, Column
from sqlalchemy import Integer, LargeBinary, String, Text, UnicodeText
from sqlalchemy.orm import sessionmaker
from sqlalchemy.ext.declarative import declarative_base
Base = declarative_base()

def atoms2smi(atoms, fname:str="/tmp/tmp_frag") -> str :
    assert isinstance(atoms, Atoms)
    write(fname+".xyz", atoms, format="xyz", append=False)
    a = list(pybel.readfile("xyz",fname+".xyz"))[0] #pybel.mol
    a.write("can", fname+".can", overwrite=True)
    with open(fname+".can") as f:
        result = f.readlines()[0].split()[0]
    return str(result)        
    

def nci_open_generation():
    fname = "/home/data/sync/hp-ChemDB_Files/sdf/" + \
        "roadmap-2011-09-23-1.sdf"
        #"NCI-Open_2012-05-01.sdf"
    with open(fname) as f:
        tmp = ""
        while True:
            line = f.readline()
            tmp = tmp + line
            if "$$$$" in line:
                try:
                    """if "99.99" in tmp and " R " in tmp:
                        _tmp = list()
                        for line in tmp.split("\n"):
                            if not("99.99" in line and " R " in line):
                                _tmp.append(line)
                        tmp = "\n".join(_tmp)"""
                    with open("/tmp/tmp.sdf", "w") as _f:
                        _f.write(tmp)
                    #atoms = read_sdf("/tmp/tmp.sdf")
                    mol = Chem.MolFromMolFile("/tmp/tmp.sdf")
                    print(Chem.MolToSmiles(mol))
                    #print()
                except Exception as e:
                    print("================= ", e)
                    #print(tmp)
                    continue
                #yield tmp, atoms
                tmp = ""

class SDF():
    def __init__(self, file_key):
        if os.path.exists(file_key):
            self.fname = file_key
        else:
            self.fname = "/tmp/SDF.sdf"
            with open(self.fname, "w") as f:
                f.write(file_key)
        self.atoms = read(self.fname, format="sdf")

    @property
    def content(self):
        with open(self.fname) as f:
            return f.read()

    @property
    def content_parted(self):
        with open(self.fname) as f:
            return f.read().split("\n\n")

    @property
    def NAME(self):
        for i in self.content_parted:
            if "CAS_NAME" in i.upper():
                return i.upper().split("\n")[-1]
        return "None"
    
    @property
    def InChI(self):
        for i in self.content_parted:
            if "INCHI" in i.upper():
                return i.upper().split("\n")[-1]
        return "None"

    @property
    def InChIKey(self):
        for i in self.content_parted:
            if "INCHI" in i.upper() and "KEY" in i.upper():
                return i.upper().split("\n")[-1]
        return "None"

    @property
    def weight(self):
        for i in self.content_parted:
            if "WEIGHT" in i.upper():
                return i.upper().split("\n")[-1]
        return "None"

    @property
    def formula(self):
        return self.atoms.symbols.formula.format("metal")

    def draw_structure(self, fname:str):
        _fname = "/tmp/tmp_frag"
        write(_fname+".xyz", self.atoms, format="xyz", append=False)
        a = list(pybel.readfile("xyz",_fname+".xyz"))[0] #pybel.mol
        a.write("can", _fname+".can", overwrite=True)
        with open(_fname+".can") as f:
            smiles = f.readlines()[0].split()[0]
        mol = Chem.MolFromSmiles(smiles)
        opts = DrawingOptions()
        opts.includeAtomNumbers = True
        draw = Draw.MolToImage(mol, options=opts)
        result = "{}.png".format(fname)
        draw.save(result)
        return result

class NCI_Open(Base):
    __tablename__ = "NCI_Open"
    id       = Column(Integer, primary_key=True, index=True)
    uuid     = Column(UnicodeText, index=True, unique=True, nullable=False)
    formula  = Column(UnicodeText, index=True)
    name     = Column(UnicodeText, index=True)
    smiles   = Column(UnicodeText, index=True)
    cas_num  = Column(UnicodeText, index=True, unique=True)
    inchi    = Column(UnicodeText, index=True, unique=True)
    inchikey = Column(UnicodeText, index=True, unique=True)
    contents = Column(LargeBinary, unique=True, nullable=False)
    def __init__(self, contents:str):
        self.uuid = hash(contents)
        sdf = SDF(contents)
        self.formula = sdf.formula
        self.name = sdf.NAME
        self.smiles = atoms2smi(atoms)
        self.cas_num = None
        self.inchi = sdf.InChI
        self.inchikey = sdf.InChIKey
        contents = pickle.dumps(contents)
        contents = zlib.compress(contents)
        self.contents = Binary(contents)

if __name__ == "__main__":
    """# setup for database
    passwd = quote_plus("chem123@G505")
    test_engine = create_engine("sqlite:///:memory:", echo=True)
    chem_cache_engine = create_engine( "mysql+pymysql://chem:{}".format(passwd) +
                            "@114.212.167.205:13306/chem_cache?charset=utf8",
        echo=True,          # 当设置为True时会将orm语句转化为sql语句打印，一般debug的时候可用
        pool_size=8,        # 连接池的大小，默认为5个，设置为0时表示连接无限制
        pool_recycle=60*30  # 设置时间以限制数据库多久没连接自动断开
    )
    Base.metadata.create_all(chem_cache_engine)

    sess = sessionmaker(bind=chem_cache_engine)
    # start handle
    num = 0
    with sess() as ss:
        for i, atoms in nci_open_generation():
            num += 1
            print(atoms)
            ss.add(NCI_Open(i))
            if int(num%100) == 0:
                ss.commit()
            if num > 1000:
                break"""

    nci_open_generation()



