#!/usr/bin/env python3

# 处理模板
import json
with open('nju_vasp_catalyst.json') as f:
    nju_vasp_catalyst = json.load(f)
template = nju_vasp_catalyst['template']
# 模版组织结构：
# 
# ```
# root(dict)---------
#     |----template(dict)
#     |        |----...
#     |        |---- 不用动
#     |        |----...
#     |----data(list)
#     |        |----1(dict)
#     |        |    |----meta(dict)
#     |        |    |----content(dict)
#     |        |    
#     |        |----2(dict)
# ```



import os, sys, shutil
import numpy as np 
import pandas as pd 
import matplotlib.pyplot as plt
from ase.io.vasp import read_vasp
from ase.io import write
from pymatgen.io.vasp.inputs import Potcar
from pymatgen.io.vasp.outputs import Vasprun, Outcar
from pymatgen.electronic_structure.core import Spin, Orbital, OrbitalType
from pymatgen.electronic_structure.dos import DOS, CompleteDos
from pymatgen.core.structure import Site, Structure
from pymatgen.core.periodic_table import Element

class VaspParse:
    def __init__(self, origin_dir='.', target_dir='.'):
        self.origin_dir = origin_dir
        self.target_dir = target_dir
        if os.path.exists(os.path.joi, parse_potcar_file=Falsen(origin_dir, 'vasprun.xml')):
            self.vasprun = Vasprun(os.path.join(origin_dir, 'vasprun.xml'), parse_potcar_file=False)
        else:
            print('DIR: '+origin_dir+" don't have a vasprun.xml file, please CHECK!")
            self.vasprun = None
        
    @property
    def label(self):
        result = self.origin_dir.split(r'/')[-1]
        return result
    
    @property
    def incar(self):
        incar = self.vasprun.parameters
        result = {'ENCUT':  float(self.vasprun.incar['ENCUT']),   #ENCUT  (数值型)'
                  'EDIFF':  float(incar['EDIFF']),                #EDIFF  (数值型)'
                  'EDIFFG': float(incar['EDIFFG']),               #EDIFFG (数值型)', 
                  'LDAU':   str(incar['LDAU']),                   #LDAU   (字符串型)', 
                  'NSW':    int(incar['NSW']),                    #NSW    (数值型)', 
                  'NELM':   int(incar['NELM']),                   #NELM   (数值型)', 
                  'ISIF':   int(incar['ISIF']),                   #ISIF   (数值型)', 
                  'ISMEAR': int(incar['ISMEAR']),                 #ISMEAR (数值型)', 
                  'SIGMA':  float(incar['SIGMA']),                #SIGMA  (数值型)', 
                  'IBRION': int(incar['IBRION']),                 #IBRION (数值型)', 
                  'PREC':   str(self.vasprun.incar['PREC']),      #PREC   (字符串型)', 
                  'ISYM':   int(incar['ISYM']),                   #ISYM   (数值型)', 
                  'ISPIN':  int(incar['ISPIN']),                  #ISPIN  (数值型)', 
                  'ALGO':   str(self.vasprun.incar['ALGO']),      #ALOG   (字符串型)', 
                  'LREAL':  str(incar['LREAL'])}                  #LREAL  (字符串型)'
        return result
    
    @property
    def kpoints(self):
        kpoints = self.vasprun.kpoints.as_dict()
        result = {'format': 'Auto' if kpoints['nkpoints'] <= 0 
                                   else str(kpoints['nkpoints']), 
                  'Grids-Type': str(kpoints['generation_style']),  #Grids-Type (字符串型) 
                  'K-mesh': {'N1': int(kpoints['kpoints'][0][0]),   #N1 (数值型)
                             'N2': int(kpoints['kpoints'][0][1]),   #N2 (数值型) 
                             'N3': int(kpoints['kpoints'][0][2]),}, #N3 (数值型) 
                  'Option shift': {'S1': float(kpoints['shift'][0]),   #S1 (数值型) 
                                   'S2': float(kpoints['shift'][1]),   #S2 (数值型) 
                                   'S3': float(kpoints['shift'][2]),}} #S3 (数值型)
        return result
    
    @property
    def poscar(self):
        #plot POSCAR image
        dir_path = os.path.join(self.target_dir, self.label)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        filename = os.path.join(dir_path, 'POSCAR.png')
        poscar = read_vasp(os.path.join(self.origin_dir, 'POSCAR'))
        write(filename, poscar*(3,3,1), format='png', rotation='10z,-80x')
        #parse data
        poscar = self.vasprun.structures[0]
        result = {'Elements': ''.join(poscar.symbol_set),  #Elements (字符串型)
                  '真空层':    round(float(poscar.lattice.c - max([i.z for i in poscar.sites])), 2),  #真空层 (数值型)
                  '简图':     [r'/'.join(filename.split(r'/')[-2:]), ]}
        return result
    
    @property
    def potcar(self):
        with open(os.path.join(self.origin_dir, 'POTCAR')) as f:
            data = f.read().split('\n')
        tmp = []
        for line in data:
            if 'TITEL' in line:
                tmp.append(line)
        return {'TITLE': '\n'.join(tmp)} #TITLE (字符串型)
    
    
    def band_center(self, spd=''):
        #判断spd轨道类型
        if spd == 'd':
            type_spd = OrbitalType.d
        elif spd == 'p':
            type_spd = OrbitalType.p
        elif spd == 's':
            type_spd = OrbitalType.s
        else:
            print('parameter "spd" error! only "s" "p" "d" can be identified.')
            return
        #计算band center
        out = self.vasprun
        data = out.complete_dos
        energies = data.energies
        delta_E = np.average(energies[1:]-energies[:-1])
        dos = data.get_spd_dos()[type_spd]
        if not out.is_spin:
            densities = dos.densities
        else:
            densities = 0.5*(dos.densities[Spin.up] + dos.densities[Spin.down])
        df = pd.DataFrame()
        df['E'] = energies
        df['DOS'] = densities
        df['DOS.dE'] = densities * delta_E
        df['E.DOS.dE'] = df['E'] * df['DOS.dE']
        return np.around(np.sum(df['E.DOS.dE']) / np.sum(df['DOS.dE']), 2)
        
    def plot_dos(self, types=''):
        #判断DOS类型
        if types == 'total':
            return str( self._plot_tdos() )    #总DOS
        elif types == 'spd':
            return str( self._plot_spddos() )  #spd分DOS
        elif types == 'element':
            return str( self._plot_edos() )    #元素类型分DOS
        else:
            print('parameter "types" error! only "total" "spd" "element" can be identified.')
            return ''
        
    def _plot_tdos(self):
        dir_path = os.path.join(self.target_dir, self.label)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        filename = os.path.join(dir_path, 'DOS_total.png')
        data = self.vasprun.tdos
        plt.figure(figsize=(3,7), dpi=150)
        plt.plot(data.densities[Spin.up], data.energies, label='total dos')
        plt.axhline(self.vasprun.efermi, color='black', linestyle='-.', label='fermi level')
        plt.axvline(0, color='gray', linewidth=0.8, linestyle='--')
        plt.legend()
        plt.xlabel('DOS')
        plt.ylabel(r'energies / eV')
        plt.title('total DOS')
        plt.savefig(filename, format='png')
        return r'/'.join(filename.split(r'/')[-2:])
    
    def _plot_spddos(self):
        dir_path = os.path.join(self.target_dir, self.label)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        filename = os.path.join(dir_path, 'DOS_spd.png')
        data = self.vasprun.complete_dos.get_spd_dos()
        plt.figure(figsize=(3,7), dpi=150)
        labels = ['s', 'p', 'd']
        keys   = list(data.keys())
        for i in range(len(keys)):
            tmp = data[keys[i]]
            plt.plot(tmp.densities[Spin.up], tmp.energies, label=labels[i])
        plt.axhline(test_vasprun.efermi, color='black', linestyle='-.', label='fermi level')
        plt.axvline(0, color='gray', linewidth=0.8, linestyle='--')
        plt.legend()
        plt.xlabel('DOS')
        plt.ylabel(r'energies / eV')
        plt.title('spd DOS')
        plt.savefig(filename, format='png')
        return r'/'.join(filename.split(r'/')[-2:])

    def _plot_edos(self):
        dir_path = os.path.join(self.target_dir, self.label)
        if not os.path.exists(dir_path):
            os.makedirs(dir_path)
        filename = os.path.join(dir_path, 'DOS_element.png')
        data = self.vasprun.complete_dos.get_element_dos()
        plt.figure(figsize=(3,7), dpi=150)
        keys   = list(data.keys())
        for i in range(len(keys)):
            tmp = data[keys[i]]
            plt.plot(tmp.densities[Spin.up], tmp.energies, label=keys[i])
        plt.axhline(test_vasprun.efermi, color='black', linestyle='-.', label='fermi level')
        plt.axvline(0, color='gray', linewidth=0.8, linestyle='--')
        plt.legend()
        plt.xlabel('DOS')
        plt.ylabel(r'energies / eV')
        plt.title('spd DOS')
        plt.savefig(filename, format='png')
        return r'/'.join(filename.split(r'/')[-2:])
        
    
    @property
    def dos_image(self):
        result = {'d band center': self.band_center('d'),   # d band center (数值型)
                  'p band center': self.band_center('p'),   # p band center (数值型)
                  's band center': self.band_center('s'),   # s band center (数值型)
                  'tDOS图':        [ self.plot_dos('total')   ],
                  'spd-分DOS图':   [ self.plot_dos('spd')     ],
                  '元素-分DOS图':  [ self.plot_dos('element') ]}
        return result
    
    @property
    def outcar(self):
        lattice = self.vasprun.structures[-1].lattice
        result = {'Volume': np.around(lattice.volume, 3),                     #Volume (数值型)
                  'Fermi energy':  np.around(self.vasprun.efermi, 3),         #Fermi energy (数值型)
                  'total energy ': np.around(self.vasprun.final_energy, 3),   #total energy  (数值型)
                  'Latice Parameter': {'a': np.around(lattice.a, 3),          # a (数值型)
                                       'b': np.around(lattice.b, 3),          # b (数值型)
                                       'c': np.around(lattice.c, 3),          # c (数值型)
                                       'alpha': np.around(lattice.alpha, 3),  # alpha (数值型)
                                       'beta':  np.around(lattice.beta,  3),  # beta  (数值型)
                                       'gamma': np.around(lattice.gamma, 3)}, # gamma (数值型)
                  'Area': np.around(lattice.volume/lattice.c, 3)}             # Area  (数值型)
        return result
    
    @property
    def cal_files(self):
        def check_file(name):
            #路径不存在则创建
            dir_path = os.path.join(self.target_dir, self.label)
            if not os.path.exists(dir_path):
                os.makedirs(dir_path)
            #文件不存在则复制
            filename = os.path.join(dir_path, name)
            if not os.path.exists(filename):
                shutil.copy(os.path.join(self.origin_dir, name), filename)
            return r'/'.join(filename.split(r'/')[-2:])
        result = {'输入文件': { 'INCAR':   [ check_file('INCAR')   ],
                               'POSCAR':   [ check_file('POSCAR')  ],
                               'POTCAR':   [ check_file('POTCAR')  ],
                               'KPOINTS':  [ check_file('KPOINTS') ]},
                  'CONTCAR': [ check_file('CONTCAR') ],
                  'DOSCAR':  [ check_file('DOSCAR')  ],
                  'OUTCAR':  [ check_file('OUTCAR')  ]}
        return result

    @property
    def all_result_in_moban_nju_vasp_catalyst(self):
        return {
            'meta' : {  
                  '数据 ID':  self.label,
                  '标题':     self.label,
                  'DOI':      '',
                  '数据摘要': self.label,
                  '关键词':   self.label,
                  '来源':     'MGE-SOURCE_HEADER v1 1000 10 #',
                  '引用':     '',
                  '其他信息':        'project: 2017YFB0702800；subject: 2017YFB0702801',
                  '数据生产机构':    '南京大学化学化工学院',
                  '数据生产者':      '史涛涛、刘高勇、陈兆旭',
                  '公开时间': '0',
                  '公开范围': '0'},
            'content':{
                   'Lable':   self.label,   #Lable (字符串型)',
                   'INCAR':   self.incar,
                   'KPOINTS': self.kpoints,
                   'POSCAR':  self.poscar,
                   'POTCAR':  self.potcar,
                   'DOS图像':  self.dos_image,
                   'OUTCAR':  self.outcar,
                   '计算文件':  self.cal_files,
                   '数据产生及校对': {
                            '模型构建人员': '史涛涛',
                            '模型计算人员': '史涛涛',
                            '数据校对人员': '陈兆旭，史涛涛，孙宏亮，刘高勇'
                   }
            }
        }
        
def all_vasp_cal_dir_list(dirname):
    if os.path.exists(os.path.join(dirname, 'vasprun.xml')):
        return dirname




if __name__ == '__main__':
    origin_dir = sys.argv[1]
    # 生成解析数据列表，并将必要的文件复制到必要的目录
    target_dir = './data'
    data = []
    for i in all_vasp_cal_dir_list(origin_dir):
        parse = VaspParse('/home/lgy/VerySync/hp-File-NJU-cal/Ag_Cu_Ag' , target_dir)
        data.append(parse.all_result_in_moban_nju_vasp_catalyst)
    # 将数据组合起来，并写入json文件中
    result = {'template': template, 
              'data'    : data    }
    if not os.path.exists(target_dir):
        os.makedirs(target_dir)
    with open(os.path.join(target_dir, 'data.json'), 'w') as f:
        json.dump(result, f, indent=4, separators=(',', ': '))






