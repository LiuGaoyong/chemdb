#!/usr/bin/env python3
import os, hashlib
import zlib, pickle
from bson.binary import Binary

def get_all_files(dirname:str=".") -> list:
    assert os.path.exists(dirname)
    dirname = os.path.abspath(dirname)
    files_list = []
    for root, dirs, files in os.walk(dirname):
        for file in files:
            files_list.append(os.path.join(root, file))
    return files_list

def get_all_dirs(dirname:str=".") -> list:
    assert os.path.exists(dirname)
    dirname = os.path.abspath(dirname)
    dirs_list = []
    for root, dirs, files in os.walk(dirname):
        for _dir in dirs:
            dirs_list.append(os.path.join(root, _dir))
    return dirs_list



class FileContents():
    def __init__(self, contents, is_filename=False):
        if isinstance(contents, bytes):
            self.contents = contents.decode('utf-8')
        elif isinstance(contents, str):
            self.contents = contents
        else:
            raise KeyError("contents must be str|bytes.")
    
    @property
    def md5(self) -> str:
        md5 = hashlib.md5(self.to_bytes())
        return str(md5.hexdigest())

    def to_bytes(self) -> bytes:
        return self.contents.encode('utf-8')

    def to_compress_binary(self) -> Binary:
        data = pickle.dumps(self.contents)
        return Binary(zlib.compress(data))


    @staticmethod
    def from_filename(fname):
        assert os.path.exists(fname)
        with open(fname) as f:
            data = f.read()
        return FileContents(data)

    @staticmethod
    def from_binary(contents):
        data = zlib.decompress(contents)
        return FileContents(pickle.loads(data))
