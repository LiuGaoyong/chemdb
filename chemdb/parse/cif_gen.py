#!/usr/bin/env python3

import os,sys,random
import pandas as pd
sys.path.append("/home/lgy/Documents/GitHub/chemdb/chemdb/parse")
from mongo import Mongo
from cif import CIF


test = Mongo()

num,result = 0, []
for i in ("cod",
          "amcsd"):
    a = test.all_id_from(i,1000)
    a = list(random.sample(a,20))
    _num = 0
    for key in a:
        num += 1
        _num += 1
        if _num == 10:
            break
        try:
            sdf = CIF(test.get_cache(i, key))
            fo  = sdf.formula
            sg = sdf.spacegroup
            name = fo
            _id = "{:02d}_{}".format(num, fo)
            img = sdf.draw_structure("{:03d}".format(num))
            elements = "-".join(set(sdf.atoms.symbols[:]))
            latt_para = sdf.atoms.cell.cellpar()
            assert len(latt_para) == 6
            tmp = [_id, fo, name, sg ]
            tmp.extend(latt_para)
            tmp.extend([img, elements])
            print(tmp)
            result.append(tmp)
        except Exception as e:
            num -= 1
            _num -= 1
            continue
result = pd.DataFrame(result, columns=["id", "formula","name","spacegroup", "a", "b",
    "c", "alpha", "beta", "gamma", "image_path", "elements"])
result.to_csv("data.csv")
print(result)
