#!/usr/bin/env python3

import os,sys,random
import pandas as pd
sys.path.append("/home/lgy/Documents/GitHub/chemdb/chemdb/parse")
from mongo import Mongo
from cif import CIF


test = Mongo()
if not os.path.exists("sdf_gen"):
    os.mkdir("sdf_gen")
os.chdir("sdf_gen")
num,result = 0, []
for i in ("cod",
          "amcsd"):
    a = test.all_id_from(i,1000)
    a = list(random.sample(a,30))
    _num = 0
    for key in a:
        _id  = None
        img  = None
        name = None
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
        num += 1
        _num += 1
        if _num == 10:
            break
        try:
            sdf  = CIF(test.get_cache(i, key))
        except Exception as e:
            print("sdfasd")
            num -= 1
            _num -= 1
            continue
        if True:
            formula  = sdf.formula
            spacegroup = sdf.spacegroup
            name = formula
            _id = "{:07d}".format(num)
            if not os.path.exists("img"): os.makedirs("img")
            img = sdf.draw_structure("img/crystal-{:03d}".format(num))
            elements = "-".join([str(i) for i in sorted(list(set(sdf.atoms.numbers)))])
            a,b,c,alpha,beta,gamma = sdf.atoms.cell.cellpar()
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
