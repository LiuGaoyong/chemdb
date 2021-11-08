#!/usr/bin/env python3
import os, random, numpy as np, pandas as pd
from pymatgen.io.vasp.outputs import Vasprun
from pymatgen.electronic_structure.core import OrbitalType, Spin
from ase.io.vasp import read_vasp
from ase.io import write


def get_all_folder_has(top:str=".", key:str="/") -> list:
    assert os.path.exists(top)
    result = []
    for root,dirs,files in os.walk(top):
        for name in files:
            fname = os.path.join(root, name)
            if key in fname:
                result.append(os.path.dirname(fname))
    return result

class VASP():
    def __init__(self, folder:str):
        assert os.path.exists(folder)
        self.folder = folder
        if os.path.exists(os.path.join(folder, 'vasprun.xml')):
            self.vasprun = Vasprun(os.path.join(folder, 'vasprun.xml'), parse_potcar_file=False)
        else:
            raise KeyError('DIR: '+folder+" don't have a vasprun.xml file, please CHECK!")

    @property
    def atoms(self):
        if os.path.exists(os.path.join(self.folder, 'CONTCAR')):
            return read_vasp(os.path.join(self.folder, 'CONTCAR'))
        else:
            raise KeyError('DIR: '+self.folder+" don't have a vasprun.xml file, please CHECK!")

    @property
    def energy(self) -> float :
        return float(self.vasprun.final_energy)

    @property
    def efermi(self) -> float :
        return float(self.vasprun.efermi)

    @property
    def formula(self):
        return self.atoms.symbols.formula.format("metal")

    def draw_structure(self, fname:str):
        result = "{}.png".format(fname)
        write(result, self.atoms, rotation='10z,-80x')
        return result

    def band_center(self, spd='') -> float :
        assert spd in ("s", "p", "d")
        #判断spd轨道类型
        if spd == 'd':
            type_spd = OrbitalType.d
        elif spd == 'p':
            type_spd = OrbitalType.p
        elif spd == 's':
            type_spd = OrbitalType.s
        else:
            raise KeyError()('parameter "spd" error! only "s" "p" "d" can be identified.')
            return
        #计算band center
        data = self.vasprun.complete_dos
        energies = data.energies
        delta_E = np.average(energies[1:]-energies[:-1])
        dos = data.get_spd_dos()[type_spd]
        if not self.vasprun.is_spin:
            densities = dos.densities
        else:
            densities = 0.5*(dos.densities[Spin.up] + dos.densities[Spin.down])
        dos_de   = densities * delta_E
        e_dos_de = energies * dos_de
        return np.float32(np.sum(e_dos_de)/np.sum(dos_de))


if __name__ == "__main__":
    test = get_all_folder_has("/home/lgy/Documents/0-MGE/1/", "OUTCAR")
    test = random.sample(test, 15)
    result, num = [],0
    for i in test:
        num += 1
        _id = "{:03d}".format(num)
        i = VASP(i)
        fo = i.formula
        name = i.formula
        eng = i.energy
        ef = i.efermi
        s = i.band_center("s")
        p = i.band_center("p")
        d = i.band_center("d")
        tmp = [_id, fo, name, eng, ef, s,p,d]
        tmp.extend(i.atoms.cell.cellpar())
        img = i.draw_structure("{:03d}".format(num))
        elements = "-".join(sorted(list(set(i.atoms.symbols[:]))))
        tmp.extend([img, elements])
        print(tmp)
        result.append(tmp)
    result = pd.DataFrame(result, columns=[
        "id", "formula","name","energy", "efermi", 
        "s band center", "p band center", "d band center",
        "a", "b", "c", "alpha", "beta", "gamma", 
        "image_path", "elements"])
    result.to_csv("data.csv")
    print(result)

        
