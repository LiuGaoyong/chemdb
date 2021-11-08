

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
        print(sdf.atoms)
        self.formula = sdf.formula
        self.name = sdf.NAME
        self.smiles = atoms2smi(atoms)
        self.cas_num = None
        self.inchi = sdf.InChI
        self.inchikey = sdf.InChIKey
        contents = pickle.dumps(contents)
        contents = zlib.compress(contents)
        self.contents = Binary(contents)