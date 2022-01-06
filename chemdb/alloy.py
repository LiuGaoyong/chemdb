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

    def downloads(self, fname:str):
        result = "{}.xyz".format(fname)
        write(result, self.atoms)
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
    if not os.path.exists("alloy_gen"):
        os.mkdir("alloy_gen")
    os.chdir("alloy_gen")
    test = get_all_folder_has("/media/lgy/Elements/23alloy_slab_adsorption_z_diff/1", "OUTCAR")
    test = random.sample(test, 15)
    result, num = [],0
    for i in test:
        _id  = str("{:08d}".format(num))
        img  = None
        name = str(i.split("/")[-1])
        formula = None
        weight  = None
        cas     = None
        inchi     = None
        inchikey  = None
        downloads = None
        elements  = None
        typed     = None
        spacegroup = None
        a,b,c = None, None, None
        alpha, beta, gamma = None, None, None
        energy = None
        efermi = None
        s_band, p_band, d_band = None, None, None
        attach = None
        try:
            num += 1
            _id = "{:07d}".format(num)
            i = VASP(i)
        except Exception as e:
            print("fdgasdf")
            num -= 1
            continue
        formula = i.formula
        name = i.formula
        energy = i.energy
        efermi = i.efermi
        typed  = "alloy"
        s_band = i.band_center("s")
        p_band = i.band_center("p")
        d_band = i.band_center("d")

        if not os.path.exists("imgs"): os.makedirs("imgs")
        if not os.path.exists("files"): os.makedirs("files")
        img = "/"+i.draw_structure("img/alloy-{}".format(_id))
        downloads = "/"+i.downloads("img/alloy-{}".format(_id))
        elements = "-".join([str(i) for i in sorted(list(set(i.atoms.numbers)))])
        a,b,c,alpha,beta,gamma = i.atoms.cell.cellpar()

        tmp = [_id, img, name, formula, weight, cas, inchi, inchikey,
                downloads, elements, typed, spacegroup, a, b, c,
                alpha, beta, gamma, energy, efermi, 
                s_band, p_band, d_band, attach]
        result.append(tmp)
    result = pd.DataFrame(result, columns=['id', 'img', 'name', 
                'formula', 'weight', 'cas', 'inchi', 'inchikey',
                'downloads', 'elements', 'typed', 'spacegroup', 'a', 'b', 'c',
                'alpha', 'beta', 'gamma', 'energy', 'efermi', 
                's_band', 'p_band', 'd_band', 'attach'])
    result.to_csv("data.csv")
    print(result)
    os.chdir("..")
