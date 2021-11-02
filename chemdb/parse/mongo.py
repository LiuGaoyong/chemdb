#!/usr/bin/env python3
import pickle, zlib
from bson.binary import Binary
from pymongo import MongoClient

class Mongo():
    def __init__(self) -> None:
        self.client = MongoClient("mongodb://root:LGY19909" + \
            "07u-mongo@114.212.167.205:27017")

    def show(self):
        mongo_client = self.client
        for db in mongo_client.list_database_names():
            print('MongoDB DataBase: {:^18}'.format(db))
            for col in mongo_client[db].list_collection_names():
                print('      Collection: {n1:28} is {n2:^9}'.format(n1=col, 
                        n2= mongo_client[db][col].estimated_document_count() ))
            print()

    def all_id_from(self, col:str, cutoff:int=100):
        result, num = [],0
        col = self.client["cache"][col]
        for i in col.find():
            result.append(i['_id'])
            num = num + 1
            if num > cutoff:
                break
        return result

    def get_cache(self, col:str, key):
        key = str(key)
        col = self.client["cache"][col]
        record = col.find_one({"_id": key})
        if not (record is None):
            return pickle.loads(zlib.decompress(
                    record['result']))
        else:
            raise KeyError("{} dose not exist.".format(key))



if __name__ == "__main__":
    test = Mongo()

    #print(test.get_cache("amcsd", 1))
    print(test.get_cache("sdf_roadmap-2011-09-23-1", 6603517))

    print(test.all_id_from("sdf_roadmap-2011-09-23-1", 1000))
