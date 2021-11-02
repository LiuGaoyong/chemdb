#!/usr/bin/env python3

import os,sys,random
import pandas as pd
sys.path.append("/home/lgy/Documents/GitHub/chemdb/chemdb/parse")
from mongo import Mongo
from sdf import SDF


test = Mongo()

num,result = 0, []
for i in ("sdf_roadmap-2011-09-23-1",
          "sdf_NCI-Open_2012-05-01"):
    a = test.all_id_from(i,1000)
    a = list(random.sample(a,10))
    for key in a:
        num += 1
        sdf = SDF(test.get_cache(i, key))
        fo  = sdf.formula
        _id = "{:02d}_{}".format(num, fo)
        name = sdf.NAME
        wght = sdf.weight
        cas  = "None"
        inchi = sdf.InChI
        inchikey = sdf.InChIKey
        try:
            img = sdf.draw_structure("{:03d}".format(num))
            print("{}\t".format(_id) +
              "{}\t".format(fo) +
              "{}\t".format(name) +
              "{}\t".format(wght) +
              "{}\t".format(cas) +
              "{}\t".format(inchi) +
              "{}\t".format(inchikey) +
              "{}\t".format(img) )
            elements = "-".join(set(sdf.atoms.symbols[:]))
            tmp = [_id, fo, name, wght, cas, inchi, inchikey, img, elements]
            result.append(tmp)
        except Exception as e:
            num -= 1
            continue
result = pd.DataFrame(result, columns=["id", "formula", "name", "weight",
    "cas", "inchi", "inchikey", "image_path", "elements"])
result.to_csv("data.csv")
print(result)
