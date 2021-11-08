#!/usr/bin/env python3
import os
from ase.io import read,write
from ase.spacegroup import get_spacegroup



class CIF():
    def __init__(self, file_key):
        if os.path.exists(file_key):
            self.fname = file_key
        else:
            self.fname = "/tmp/CIF.cif"
            with open(self.fname, "w") as f:
                f.write(file_key)
        self.atoms = read(self.fname, format="cif")

    @property
    def content(self):
        with open(self.fname) as f:
            return f.read()

    @property
    def spacegroup(self) -> str:
        sg = get_spacegroup(self.atoms)
        return "{}({})".format(str(sg.symbol
            ).replace(" ", ""), sg.no)


    @property
    def formula(self):
        return self.atoms.symbols.formula.format("metal")

    def draw_structure(self, fname:str):
        result = "{}.png".format(fname)
        write(result, self.atoms, rotation='10z,-80x')
        return result
        
    
        


if __name__ == "__main__":
    from mongo import Mongo
    test = Mongo()
    
