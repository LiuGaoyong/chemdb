#!/usr/bin/env python3
from setuptools import setup, find_packages

setup(
    name="chemdb",
    version="0.0.1",
    author="Liu Gaoyong",
    author_email="liugaoyon_88@163.com",
    description="The chemdb is package for save chemical data    "   + \
            "everywhere. It will take the ways of online(i.e. use " + \
            "network spider) and offline(i.e. local database file) combination.",

    # 项目主页
    url="https://github.com/LiuGaoyong/chemdb",

    # 你要安装的包，通过 setuptools.find_packages 找到当前目录下有哪些包
    packages=find_packages(
        where='.',
        exclude=['test/*'], )
)