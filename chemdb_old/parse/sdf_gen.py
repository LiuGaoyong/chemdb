#!/usr/bin/env python3

import os,sys,random
import pandas as pd
sys.path.append("/home/lgy/Documents/GitHub/chemdb/chemdb/parse")
from mongo import Mongo
from sdf import SDF


test = Mongo()
if not os.path.exists("sdf_gen"):
    os.mkdir("sdf_gen")
os.chdir("sdf_gen")
num,result = 0, []
for i in ("sdf_roadmap-2011-09-23-1",
          "sdf_NCI-Open_2012-05-01"):
    a = test.all_id_from(i,1000)
    a = list(random.sample(a,10))
    print(a)
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
        try:
            num += 1
            sdf = SDF(test.get_cache(i, key))
            formula  = sdf.formula
            _id = "{:07d}".format(num)
            name = sdf.NAME
            weight = sdf.weight
            cas  = "None"
            inchi = sdf.InChI
            inchikey = sdf.InChIKey
            if not os.path.exists("img"): os.makedirs("img")
            img = sdf.draw_structure("img/mol-{:03d}".format(num))
        except Exception as e:
            print("fdgasdf")
            num -= 1
            continue
        elements = "-".join([str(i) for i in 
            sorted(list(set(sdf.atoms.numbers)))])
        print("\n",elements, "\n")
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
