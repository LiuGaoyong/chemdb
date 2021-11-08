#!/usr/bin/env python3
import os
from ase.io import read, write
from openbabel import pybel
from rdkit import Chem
from rdkit.Chem import Draw
from rdkit.Chem.Draw.MolDrawing import DrawingOptions

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
    def smiles(self):
        write("/tmp/tmp.xyz", self.atoms, format="xyz", append=False)
        a = list(pybel.readfile("xyz","/tmp/tmp.xyz"))[0] #pybel.mol
        a.write("can", "/tmp/tmp.can", overwrite=True)
        with open("/tmp/tmp.can") as f:
            result = f.readlines()[0].split()[0]
        return str(result)  

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
        
    
