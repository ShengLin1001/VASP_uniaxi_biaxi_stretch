#!/usr/bin/env python
# coding: utf-8

# # Postprocess for uniaxial and biaxial stretch
# 

# In[110]:


# pre-function
from re import findall, search, IGNORECASE
from os import listdir
from numpy import sqrt, array
from inspect import getargvalues, currentframe, stack

dir_path = './y_dir/'
stretch_list = [0.997, 0.998, 0.999, 1.000, 1.001, 1.002, 1.003]

POSTTIME_CONFIG = [
    # PostTime class items
    {'char': "Iteration",                 'basic_pattern': '', 'match_type': r'\d+\(\s*\d+\)',                   'switch': 1, 'basic_content1':'Iteration:','myindex1':0    },
    {'char': 'Elapsed time',              'basic_pattern': '', 'match_type': 'float',                         'switch': 1, 'basic_content1':'Elapsed time ','myindex1':0 },
    {'char': 'NCORES_PER_BAND',           'basic_pattern': '', 'match_type': r'(\d+\s+cores,\s+\d+\s+groups)',   'switch': 1, 'basic_content1':'','myindex1':0  },
    {'char': "Maximum memory used (kb):", 'basic_pattern': '', 'match_type': 'float',                           'switch': 1, 'basic_content1':'memory:','myindex1':0  }, 
    {'char': 'LOOP:  cpu time',           'basic_pattern': r'real\stime\s*', 'match_type': 'float',         'switch': 1, 'basic_content1':'LOOP:','myindex1':0  },
]
POSTDATA_CONFIG = [
    # PostData class items
    {'char': "energy(sigma->0)", 'basic_pattern': r'energy\(sigma->0\) =\s*', 'match_type': 'float', 'switch': 1, 'basic_content1':'','myindex1':0  , 'search_once' : False},
    {'char': "EENTRO",           'basic_pattern': r'EENTRO\s*=\s*',           'match_type': 'float', 'switch': 1, 'basic_content1':'','myindex1':0  , 'search_once' : False},
    {'char': "in kB",            'basic_pattern': '',                         'match_type': 'float', 'switch': 6, 'basic_content1':'in kB ', 'search_once' : False},
]
POSTDATA2_CONFIG = [
    # PostData2 class items
    {'char': "volume of cell :",       'basic_pattern': '', 'match_type': 'float', 'switch': 1, 'basic_content1':'','myindex1':0   },
    {'char': 'external pressure',      'basic_pattern': '', 'match_type': 'float', 'switch': 1, 'basic_content1':'','myindex1':0   },
#    {'char': 'FORCES acting on ions',                      'match_type': 'scientific notation', 'length': 12},
    {'char': 'TOTAL-FORCE (eV/Angst)',                      'match_type': 'float', 'length': 6},
]          
POSTPARAM_CONFIG = [                                                                           
    # PostParam class items
    {'char': "ENCUT",    'basic_pattern': '',                   'match_type': 'float', 'switch': 1, 'basic_content1':'ENCUT=','myindex1':0   },
    {'char': "ISMEAR",   'basic_pattern': '',                   'match_type': 'float', 'switch': 2, 'basic_content1':'ISMEAR=','myindex1':0, 'basic_content2':';SIGMA=','myindex2':1    },
    {'char': 'EDIFF  =', 'basic_pattern': r'EDIFF(?!G)\s+=\s+', 'match_type': 'scientific notation', 'basic_content1':'EDIFF=','myindex1':0   },
    {'char': 'EDIFFG =', 'basic_pattern': r'EDIFFG\s+=\s+',     'match_type': 'scientific notation', 'basic_content1':'EDIFFG=','myindex1':0   },
    {'char': 'LREAL',    'basic_pattern': r'LREAL\s*=\s*',      'match_type': 'bool',   'switch': 1, 'basic_content1':'LREAL=','myindex1':0 },
    {'char': 'NKPTS',    'basic_pattern': '',                   'match_type': 'int',   'switch': 1, 'basic_content1':'NKPTS=','myindex1':0   },
    {'char': 'NBANDS',   'basic_pattern': r'NBANDS=\s*',        'match_type': 'int',   'switch': 1, 'basic_content1':'EBANDS=','myindex1':0   }
]
POSTPARAM2_CONFIG = [                                                                           
    # PostParam2 class items
    {'char': "ISTART", 'basic_pattern': '',                   'match_type': 'int',  'switch': 1, 'basic_content1':'ISTART=','myindex1':0   },
    {'char': "ICHARG", 'basic_pattern': '',                   'match_type': 'int',  'switch': 1, 'basic_content1':'ICHARG=','myindex1':0   },
    {'char': 'PREC',   'basic_pattern': r'PREC\s*=\s*',       'match_type': 'char', 'switch': 1,'basic_content1':'PREC=','myindex1':0   },
    {'char': 'ISPIN',  'basic_pattern': '',                   'match_type': 'int', 'switch': 1, 'basic_content1':'ISPIN=','myindex1':0   },
    {'char': 'VOSKOWN','basic_pattern': '',                   'match_type': 'int',  'switch': 1,'basic_content1':'VOSKOWN=','myindex1':0   },
    {'char': 'RWIGS',  'basic_pattern': r'RWIGS\s+=\s+',      'match_type': 'float', 'end_pattern':r'\s*;', 'switch': 1, 'basic_content1':'RWIGS=','myindex1':0  , 'search_once' : True},
    {'char': 'IALGO',  'basic_pattern': '',                   'match_type': 'int', 'switch': 1,'basic_content1':'IALGO=','myindex1':0   },
    {'char': 'NGX',    'basic_pattern': '',                   'match_type': 'int', 'switch': 1,'basic_content1':'NGX=','myindex1':0, 'specified_num': [1, 2, 3]   },
]
POSTPARAMSTA_CONFIG = [                                                                           
    # PostSta class items
    {'char': "ENCUT",    'basic_pattern': '',                   'match_type': 'float', 'switch': 1, 'basic_content1':'ENCUT:','myindex1':0   },
    {'char': "ISMEAR",   'basic_pattern': '',                   'match_type': 'float', 'switch': 2, 'basic_content1':'ISMEAR:','myindex1':0, 'basic_content2':'   SIGMA:','myindex2':1    },
    {'char': 'EDIFF  =', 'basic_pattern': r'EDIFF(?!G)\s+=\s+', 'match_type': 'scientific notation', 'basic_content1':'EDIFF:','myindex1':0   },
    {'char': 'EDIFFG =', 'basic_pattern': r'EDIFFG\s+=\s+',     'match_type': 'scientific notation', 'basic_content1':'EDIFFG:','myindex1':0   },
    {'char': 'ISIF',    'basic_pattern': r'ISIF\s+=\s+',        'match_type': 'int', 'basic_content1':'ISIF:','myindex1':0   },
    {'char': 'LREAL',    'basic_pattern': '',                   'match_type': 'bool',   'switch': 1, 'basic_content1':'LREAL:','myindex1':0},
    {'char': 'NKPTS',    'basic_pattern': '',                   'match_type': 'int',   'switch': 1, 'basic_content1':'NKPTS:','myindex1':0   },
    {'char': "ISTART", 'basic_pattern': '',                    'match_type': 'int',  'switch': 1, 'basic_content1':'ISTART:','myindex1':0   },
    {'char': "ICHARG", 'basic_pattern': '',                    'match_type': 'int',  'switch': 1, 'basic_content1':'ICHARG:','myindex1':0   },
    {'char': 'PREC',   'basic_pattern': r'PREC\s*=\s*',        'match_type': 'char', 'switch': 1,'basic_content1':'PREC:','myindex1':0   },
    {'char': 'ISPIN',  'basic_pattern': '',                    'match_type': 'int', 'switch': 1, 'basic_content1':'ISPIN:','myindex1':0   },
    {'char': 'VOSKOWN','basic_pattern': '',                    'match_type': 'int',  'switch': 1,'basic_content1':'VOSKOWN:','myindex1':0   },
    {'char': 'RWIGS',  'basic_pattern': r'RWIGS\s+=\s+',       'match_type': 'float', 'end_pattern':r'\s*;', 'switch': 1, 'basic_content1':'RWIGS:','myindex1':0  , 'search_once' : True},
    {'char': 'IALGO',  'basic_pattern': '',                    'match_type': 'int', 'switch': 1,'basic_content1':'IALGO:','myindex1':0   },
    {'char': 'support grid',  'basic_pattern': '',             'match_type': '\s+(\w+\s+\w+)\s+NGXF', 'switch': 1,'basic_content1':'ADDGRID:','myindex1':0   },
]
POSTPARAMSTA_CHECK = [
    {'ENCUT':550.0,   'ISMEAR':1, 'SIGMA':0.20,   'EDIFF': 0.1E-07,   'EDIFFG':-.5E-03,   'ISIF':3,   'LREAL':'F',   'NKPTS':3268,
    'ISTART':0,   'ICHARG':1,   'PREC':'accura',   'ISPIN':1,   'VOSKOWN':1 ,  'RWIGS':2.840,   'IALGO':38,   'ADDGRID':'supportgrid'  },
]
POSTWARNING_CONFIG = [                                                                           
    # PostWarning class items
    {'char': "warning",   'match_type':'full line',  'move_list_down_up':[-10,5] },
    {'char': "ADVICE",   'match_type':'full line',  'move_list_down_up':[-2,7] }
]


class PostTime:
    def __init__(self, my_path=dir_path, post_path='./y_post_time.txt', name = 'Post Time', config = POSTTIME_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1 job state: relaxed? time? CPUs? memory?'):
        count_finish = 0
        count_relax = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            tag_convergence = 0
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            #print(filename)
            try:
                with open(filename, 'r') as file:
                    count_finish = count_finish + 1
                    tag_convergence = check_convergence(file)
                    count_relax = count_relax + tag_convergence
                    self.post(file, tag_convergence, formatted_i)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_finish, stretch_factor,special_name, self.name)
        self.end(total_job = len(stretch_factor), unfinished = len(stretch_factor)-count_finish, unrelaxed = len(stretch_factor)-count_relax)

    def post(self, file, tag_convergence=1, my_stretch_factor=0.001):
        mycontent = [f'{my_stretch_factor}']
        write_convergence(tag_convergence, mycontent)
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'  
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)

    def end(self, total_job=7, unfinished=0, unrelaxed=0):
        with open(self.post_path, 'a') as post_file:
            post_file.write('-------------------------\n')
            post_file.write(f'| total number of jobs: {total_job}\n')
            post_file.write(f'|     un-finished jobs: {unfinished}\n')
            post_file.write(f'|      un-relaxed jobs: {unrelaxed}\n')
            post_file.write('-------------------------\n')

class PostData:
    def __init__(self, my_path=dir_path, post_path='./y_post_data.txt', name = 'Post Data', config = POSTDATA_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1   energy(sigma->0)(eV)  EENTRO(eV)  -stress(kB)'):
        count_read = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            #print(filename)
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    self.post(file, formatted_i)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
    
    def post(self, file, my_stretch_factor=0.001):
        mycontent = [f'{my_stretch_factor}']
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
            #print(mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)

class PostData2:
    def __init__(self, my_path=dir_path , post_path='./y_post_data_2.txt', name = 'Post Data2', config = POSTDATA2_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name 
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1   volume  pressure(kB)  Fmax'):
        count_read = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            #print(filename)
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    self.post(file, formatted_i)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
    
    def post(self, file, my_stretch_factor=0.001):
        mycontent = [f'{my_stretch_factor}']
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)

class PostDiff:
    def __init__(self, my_path=dir_path  , post_path='./y_post_diff.txt', name = 'Post Diff'):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name

    def read_diff(self, stretch_factor= stretch_list, title = 'i1   diff  INCAR  KPOINTS  POSCAR  POTCAR  sub.vasp', specified_file = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'sub.vasp', 'Y_CONSTR_CELL.IN'] ):
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            directory_name = f'{self.my_path}{formatted_i}/'
            # print(directory_name)
            self.post(directory_name, formatted_i, specified_file)
        print(f'Read all diff   - {self.name}')
    
    def post(self, directory_name, my_stretch_factor=0.001, specified_file = ['INCAR', 'KPOINTS', 'POSCAR', 'POTCAR', 'sub.vasp', 'Y_CONSTR_CELL.IN']):
        content = []
        list_dir = listdir(directory_name)
        for file in specified_file:
            if file in list_dir:
                content.append(file)
        self.output_list(content, my_stretch_factor)
        
    def output_list(self, content, my_stretch_factor = 0.001):
        with open(self.post_path, 'a') as post_file:
            post_file.write('\n' + '='*15 + '\n')
            post_file.write(f'{my_stretch_factor}\n')
            post_file.write('='*15+'\n'+ '\n')
            for file in content:
                post_file.write('\n' + f'{file}\n' + '\n')
            
class PostParam:
    def __init__(self, my_path=dir_path , post_path='./y_post_param.txt', name = 'POST Param',  config = POSTPARAM_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1   input parameters'):
        count_read = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            # print(filename)
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    # print(f'read {formatted_i} {special_name}')
                    self.post(file, formatted_i)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
    
    def post(self, file, my_stretch_factor=0.001):
        mycontent = [f'{my_stretch_factor}']
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        #print(mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'
        #print(mycontent_char)
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)

class PostParam2:
    def __init__(self, my_path=dir_path , post_path='./y_post_param_2.txt', name = 'POST Param 2',  config = POSTPARAM2_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1   input parameters 2'):
        count_read = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    self.post(file, formatted_i)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
    
    def post(self, file, my_stretch_factor=0.001):
        mycontent = [f'{my_stretch_factor}']
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'
        #print(mycontent)
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)

class PostParamSta:
    def __init__(self, my_path=dir_path , post_path='./y_post_param_statistics.txt', name = 'POST Param Statistics', config = POSTPARAMSTA_CONFIG,  config2 = POSTPARAMSTA_CHECK):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config
        self.config2= config2

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = '# VASP param statistics. OK = the same in all jobs. ',
                    tag_split = '   ', tag_end = '\n\n', tag_begin='\n', left = False, num=13):
        count_read = 0
        all_is_same_list = []
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            # print(filename)
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    # print(f'read {formatted_i} {special_name}')
                    self.post(file, formatted_i, all_is_same_list)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
        write_check(all_is_same_list, self.post_path, self.config2, tag_split, tag_end, tag_begin, left, num)

    def post(self, file, my_stretch_factor=0.001, all_is_same_list=[]):
        #print('1')
        #print(all_is_same_list)
        # mycontent = [f'{my_stretch_factor}']
        mycontent = []
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        #print(mycontent)
        mycontent_char = list_to_char(mycontent)
        #print(mycontent_char)
        mycontent_dic = char_to_dic(mycontent_char)
        #print(mycontent_dic)
        all_is_same_list = check_same(my_stretch_factor, mycontent_dic, self.config2, all_is_same_list)
        #print('2')
        #print(all_is_same_list)
        # all_is_same_char = list_of_dic_to_char(all_is_same_list)
        #write_content_to_file(all_is_same_char, self.post_path, 'a')
        file.seek(0)

class PostWarning:
    def __init__(self, my_path=dir_path , post_path='./y_post_warning.txt', name = 'POST Warning',  config = POSTWARNING_CONFIG):
        self.my_path = my_path
        self.post_path = post_path
        self.name = name
        self.config = config

    def read_OUTCAR(self, stretch_factor= stretch_list, special_name='OUTCAR', title = 'i1   WARNING', move_list_down_up =[-10,5]):
        count_read = 0
        write_content_to_file(title + '\n', self.post_path, 'w')
        for i in stretch_factor:
            formatted_i = f'{i:.3f}'
            filename = f'{self.my_path}{formatted_i}/{special_name}'
            # print(filename)
            try:
                with open(filename, 'r') as file:
                    count_read = count_read + 1
                    # print(f'read {formatted_i} {special_name}')
                    self.post(file, formatted_i, move_list_down_up)
            except FileNotFoundError:
                print(f"In {formatted_i} directory, The file {special_name} was not found! - {self.name}")
        print_after_read(count_read, stretch_factor,special_name, self.name)
    
    def post(self, file, my_stretch_factor=0.001, move_list_down_up = [-10,5]):
        mycontent = ['-'*20 + f'\n{my_stretch_factor}\n' + '-'*20 +'\n']
        for config in self.config:
            myfindall(file, **config, mycontent=mycontent)
        mycontent_char = list_to_char(mycontent) + '\n'
        write_content_to_file(mycontent_char, self.post_path, 'a')
        file.seek(0)


def myfindall(file, char, basic_pattern = r'', match_type = 'int', end_pattern = r'', mycontent=[], switch = 1, basic_content1='', myindex1 = 0,     
            basic_content2='', myindex2 = 1, search_once = False, specified_num = [],
            move_list_down_up = [10,5],
            length = 6):
    # parameters: first,second - universal, third - fmax, 4-th - warning     
    # functions:  myfindall_universal (extract other parameters), myfindall_warning (extract warning and advice), myfindall_fmax (extract fmax)
    # functions:  myfindall_universal (extract other parameters), myfindall_warning (extract warning and advice), myfindall_fmax (extract fmax)
    #             my_switch (in universal function to simplify the code to change its type, and add more functions)
    #             find_fmax (in fmax function to extract the force component, and sqrt() it)

    def myfindall_universal (file, char, basic_pattern = r'', match_type = 'int', end_pattern = r'', mycontent=[], switch = 1, basic_content1='', myindex1 = 0,
                basic_content2='', myindex2 = 1, search_once = False, specified_num = [1,2,3,4,5]):

        def my_switch(extr_content = [], mycontent1 = [],switch = 1, basic_content1='', myindex1 = 0, basic_content2='', myindex2 = 1, specified_num = []):
            content = ''
            mycontent2 = []
            if switch == 1:
                if extr_content and myindex1 < len(extr_content):
                    content = basic_content1 + f'{extr_content[myindex1]}'
            elif switch == 2:
                if extr_content and myindex1 < len(extr_content) and myindex2 < len(extr_content):
                    #print(extr_content)
                    content = basic_content1 + f'{extr_content[myindex1]}' + basic_content2 + f'{extr_content[myindex2]}'
            elif switch == 6:
                if extr_content and myindex1 < len(extr_content) and myindex2 < len(extr_content):
                    content = basic_content1 + f'{extr_content[0]} {extr_content[1]} {extr_content[2]} {extr_content[3]} {extr_content[4]} {extr_content[5]}'
            if specified_num:
                mycontent1.append(content)
            else:
                mycontent2.append(content)
                #print(content)
                return content

        # the last parameter {specified_num} from 1, not 0. 1-row, 2-row ... etc
        extr_content = []
        count = 0
        mypattern = my_pattern_trans(match_type, char)
        full_pattern = basic_pattern + mypattern + end_pattern

        mycontent1 = []
        for line in file:      
            if char in line:
                count = count + 1
                extr_content = findall(full_pattern, line)
                #print(extr_content)
                rm_blank(extr_content)
                #print(count)
                #print(line)
                if specified_num:
                    if count in specified_num:
                        #print(line)
                        my_switch(extr_content, mycontent1, switch, basic_content1, myindex1, basic_content2, myindex2, specified_num)
                else:
                    content = my_switch(extr_content, [], switch, basic_content1, myindex1, basic_content2, myindex2)
                    if count == 1 and search_once:
                        break
        if specified_num:
            # mycontent1 = ['item1', 'item2', 'item3']
            # joined_string = 'item1 item2 item3'
            joined_string = ' '.join(mycontent1)
            mycontent.append(joined_string)
        else:
            mycontent.append(content)
        #print(mycontent)
        file.seek(0)

    def myfindall_warning (file, char, match_type, move_list_down_up = [10,5], mycontent = []):
        mypattern = my_pattern_trans(match_type, char)
        start_line = find_line_position(file, mypattern)
        if start_line:
            extract_lines = extract_line_at_position(file, start_line, move_list_down_up)
            mycontent.append(extract_lines)
        else:
            print("myfindall_warning didn't find the line")

    def myfindall_fmax (file, char, match_type = 'int', length = 6, mycontent = []):

        def find_fmax(file, char, match_type, length = 6):
            file.seek(0)
            force = []
            start_line_list = find_line_position(file, char)
            start_line = start_line_list[0]
            mypattern = my_pattern_trans(match_type, char)
            force = []
            for line_number, line in enumerate(file, start=1):
                #print(line)
                if line_number > start_line:
                    #print(line)
                    match = findall(mypattern, line)
                    #print(match)
                    if length == 6 and len(match) == length:
                        #print(line)
                        #print(match)
                        forcex = float(match[3])
                        forcey = float(match[4])
                        forcez = float(match[5])
                        temp = sqrt(forcex*forcex+forcey*forcey+forcez*forcez)
                        force.append(temp)
                        #print(match)
                    elif length == 12 and len(match) == length:
                        #print(line)
                        #print(match)
                        forcex = 0
                        forcey = 0
                        forcez = 0
                        for i in range(4):
                            #print(i)
                            forcex = forcex + float(match[i*3])
                            forcey = forcey + float(match[i*3+1])
                            forcez = forcez + float(match[i*3+2])
                        temp = sqrt(forcex*forcex+forcey*forcey+forcez*forcez)
                        force.append(temp)
            if force:
                if length == 6:
                    fmax = max(force)
                elif length == 12:
                    force.sort()
                    fmax = force[-2]
                formatted_fmax = f'{fmax:.10f}'
                file.seek(0)
                return formatted_fmax
            else:
                return None

        fmax = find_fmax(file,  char, match_type, length)
        if fmax:
            mycontent.append(f'{fmax}')
        else:
            print("myfindall_fmax didn't calculate the fmax")

    if search('force', char, IGNORECASE):
        myfindall_fmax(file, char, match_type, length, mycontent)
    elif search('warning|advice', char, IGNORECASE):
        myfindall_warning (file, char, match_type, move_list_down_up, mycontent)
    else:
        myfindall_universal (file, char, basic_pattern, match_type, end_pattern, mycontent, switch, basic_content1, myindex1,
              basic_content2, myindex2, search_once, specified_num)

##############################################################
############# change type
############# adjust variable 
############# construct the content
##############################################################

def list_to_char(mycontent = [], list_type = "char", char_tag = '   ',
                 dic_tag = '\n', key_value_interval = ': ', item_interval = ', '):
    # consisted of two function
    # list_of_char_to_char to output the list of char to file
    # list_of_dic_to_char to output the list of dic to file

    # mycontent1 = ['item1', 'item2', 'item3']
    # joined_string = 'item1 item2 item3'
    # joined_string = ' '.join(mycontent1)

    def list_of_char_to_char(mycontent = [], tag = '   '):
        # for 1-D list of char, change it to char to optput to file
        mycontent_char = ''
        for content in mycontent:
            mycontent_char = mycontent_char  + f'{content}' + tag
        return mycontent_char
    
    def list_of_dic_to_char(list_of_dicts, tag = '\n', key_value_interval = ': ', item_interval = ', '):
        # for 1-D list of dic, change it to char
        #          [{"name":"1","gender":"female"}, {"name":"2","gender":"male"}] 
        # =>>      "name: 1,gender: female
        #           name: 2, gender: male"
        formatted_dicts = [item_interval.join([f"{key}{key_value_interval}{value}" for key, value in d.items()]) + tag for d in list_of_dicts]
        return ''.join(formatted_dicts)
    
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)

    temp = ''
    if list_type == "char":
        temp = list_of_char_to_char(mycontent, char_tag )
    elif list_type == "dic":
        temp = list_of_dic_to_char(mycontent, dic_tag, key_value_interval, item_interval)
    else:
        print(f"the list_type wanted to transformed is not supported, please fix it! - {calling_function}")
    return temp

def char_to_dic(mycontent_char):
    # "name:n1 gender:male" =>> {"name":"n1", "gender":"male"}
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mycontent_dic = {}
    # Split the input line by spaces to separate parameters
    mycontent_char_copy = mycontent_char.split()
    # Create a dictionary to store the parsed parameters
    # Iterate through the parameters and parse them into key-value pairs
    for char in mycontent_char_copy:
        dic_name, dic_value = char.split(':')
        mycontent_dic[dic_name] = dic_value
    return mycontent_dic

def dic_to_char(mycontent_dic={}, tag = '   '):
    #  f"name:n1 {tag} gender:male" <<= {"name":"n1", "gender":"male"}
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    # Initialize an empty string to store the formatted content
    mycontent_char = ''
    if mycontent_dic:
        # Iterate over key-value pairs in the dictionary
        for key, value in mycontent_dic.items():
            # Append each key-value pair to the formatted string
            mycontent_char = mycontent_char + f'{key}:{value}' + tag
        # Remove the trailing space at the end of the string
        mycontent_char += '\n'
    else:
        print_after_blank('mycontent_dic', calling_function, [])
    return mycontent_char

def rm_blank(extr_content, remove_string=' ', removed_string=''):
    # for 1-D list, remove blank
    # like ["reached required accuracy","reached required accuracy"]  to 
    #      [ "reachedrequiredaccuracy" , "reachedrequiredaccuracy" ]
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    for num in range(len(extr_content)):
        extr_content[num] = extr_content[num].replace(remove_string, removed_string)
        # print(extr_content)
        # print('\n')

def my_add_list(list1 = [], list2 = []):
    # combine two list, and sort it
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    my_list = []
    for item in list1 + list2:
        if item not in my_list:
            my_list.append(item)
    my_list.sort()
    return my_list

def normalize_float_int(value):
    # Helper function to normalize values
    # Convert to float and remove trailing zeros, if it's a float
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    try:
        float_value = float(value)
        return str(float_value).rstrip('0').rstrip('.')  # Remove trailing zeros and dot if it's an int
    except ValueError:
        return value  # Leave other types as they are

def myjust(left = True, num = 10, output = ''):
    # left just or right just of output
    # num is digit of variable
    # like 'dig' =>> 'dig       '  =>> '       dig'
    # return char_out
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    char_out = ''
    if left:
        char_out = str(output).ljust(num)
    else:
        char_out = str(output).rjust(num)
    return char_out

def create_content(list = [], tag_split = '   ', tag_end = '', tag_begin = '', left = True, num = 10):
    # have a list of char, using myjust() to formatted them.
    # output a char to write in a specified file
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mycontent_char = ''
    mycontent = []
    mycontent.append(tag_begin)
    for i in list:
        formatted_i = myjust(left, num, i)
        content_i = f'{formatted_i}' + tag_split
        mycontent.append(content_i)
    mycontent.append(tag_end)
    mycontent_char = list_to_char(mycontent, char_tag = '')
    return mycontent_char

def my_down_up(list, move_list_down_up = [0, 0]):
    # enlarge list
    # like      [1,4], move_list = [-2,3]
    #      =>>  [-1,0,1,2,3,4,5,6,7]  (-1, 7)
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    my_list = []
    lower_move, upper_move = move_list_down_up
    for item in list:
        #my_item = [item + i for i in range(lower_move, upper_move + 1)]
        my_item = [item + i for i in range(lower_move, upper_move + 1)]
        my_list = my_add_list(my_list, my_item)
    return my_list

def my_pattern_trans(match_type, char = ''):
    # generate the regular expression
    # char to regular expression
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    mypattern = ''
    if match_type:
        if match_type == 'int':
            mypattern = r'([+-]?\d+)'
        elif match_type == 'float':
            mypattern = r'([+-]?\d*\.?\d+|[+-]?\d+\.)'
        elif match_type == 'scientific notation':
            mypattern = r'([+-]?\d*\.?\d+[eE][-+]?\d+)'
        elif match_type == 'bool':
            mypattern = r'\b(true|false|True|False|1|0|yes|no|Yes|No|T|F|t|f)\b'
        elif match_type == 'char':
            mypattern = r'(\w+)\s+'
        elif match_type == 'full line':
            mypattern = f'{char}'
        else:
            mypattern = match_type
    else:
        print_after_blank('match_type', calling_function, [])
    return mypattern

##############################################################
############# about convergence information 
############# note after read
############# check input 
############# check same
############# input and output 
############# find and extract the matched lines (up-down could be extract)
##############################################################

def check_input(args = [], check_type = 'none',
                char_out = '', specified_blank = ''):
    # check_input_blank function is not correct, I think it's useless
    # first line  - check_input_none
    # second line - check_input_blank

    def check_input_blank(specified_blank = '', calling_function=''):
        print(f"the value is blanked! - {calling_function}")
        return specified_blank

    def check_input_none(args = [], calling_function = ''):
        # print(calling_function)
        for arg_name, arg_value in args.items():
            if arg_value is None:
                raise ValueError(f"In '{calling_function}': Argument '{arg_name}' is missing or None")
        
    # Get the name of the calling function
    calling_function = stack()[1].function
    if check_type == 'blank':
        char_out = check_input_blank(char_out, specified_blank, calling_function)
    elif check_type == 'none':
        check_input_none(args, calling_function )
    else:
        print(f"this check type {check_type} is not supported")

def print_after_read(count_read, stretch_factor = [], special_name = '', name = ''):
    # note if read all {special_name} - OUTCAR
    # {name} is class name {self.name}
    # nead a count_read variable
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    if count_read:
        if count_read == len(stretch_factor):
            print(f'Read all {special_name} - {name}')
        else:
            print(f'Not all {special_name} are readed - {name}')
    else:
        print_after_blank('count_read', calling_function, [])

def print_after_blank(char_name = '', calling_function = '', specified_blank = ''):
    # combined two lines
    #         print(f"the value is None! - {calling_function}")
    #         char_out = ''
    print(f"the || {char_name} || is blanked! - {calling_function}")
    return specified_blank

def print_after_not_supported(char_name = '', calling_function = '', specified_blank = ''):
    # combined two lines
    #         print(f"the value is None! - {calling_function}")
    #         char_out = ''
    print(f"the || {char_name} || is None! - {calling_function}")
    return specified_blank

def print_after_cant_read(char_name = '', calling_function = '', specified_blank = ''):
    # combined two lines
    #         print(f"the value is None! - {calling_function}")
    #         char_out = ''
    print(f"the || {char_name} || is None! - {calling_function}")
    return specified_blank

def check_convergence(file, special_char = "reached required accuracy"):
    # check if {special_char} in file
    # return tag
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    tag = 0
    if file:
        for line in file:
            if special_char in line:
                #print(f'{formatted_i} {special_char}')
                tag = 1
        file.seek(0)
    else:
        print_after_blank('file', calling_function, [])     
    return tag

def write_convergence(tag_convergence, mycontent):
    # add convergence to {mycontent}
    # append char
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    if tag_convergence == 1:
        mycontent.append('reachedrequiredaccuracy')
    elif tag_convergence == 0:
        mycontent.append('didnot reached accuracy')
    else:
        print(f"tag_convergence is not right! - {calling_function}")

def check_same(my_stretch_factor=0.001, mycontent_dic = [{'ENCUT': 300, 'ISMEAR': 0.20}], mycontent_check = POSTPARAMSTA_CHECK, all_is_same_list = []):
    # Now we compare the parameter in my_stretch_factor with the mycontent_check
    # we have two dictionaries, mycontent_dic (dictionary), mycontent_check (list of dictionary)
    # do a loop for mycontent_check, if the value of mycontent_dic is not same as those of my_content_check
    # save a flase value for is_same_dic[key]
    # combined the results with all_is_same_list (a list of a dictionary)  [{'ENCUT': 300, 'ISMEAR': 0.20}]
    # after search on every 7 my_stretch_factor, we have a all_is_same_list ( a list of 7 dictionary)
    # Note, for folat or int, normalize_float_int to keep 0.2 is same as 0.20
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    #print(mycontent_dic)
    #print(mycontent_check)
    # Initialize is_same as a dictionary
    if mycontent_dic and mycontent_check:
        is_same_dic = {}
        is_same_dic[str('NAME')] = my_stretch_factor
        # Iterate over each parameter dictionary in mycontent_check
        for param_check in mycontent_check:
            # all_match = True  # Flag to track if all parameters in param_check match
            for key, value in param_check.items():
                single_match = True
                #print()
                value_dic = normalize_float_int(mycontent_dic[key])
                value_check = normalize_float_int(value)
                if value_dic != value_check:
                    # all_match = False
                    single_match = False
                    #break  # No need to check further, exit inner loop
                # Add the result to is_same with the parameter dictionary as the key
                is_same_dic[str(key)] = single_match
        all_is_same_list.append(is_same_dic)
        #print(is_same_list)
    else:
        print_after_blank('mycontent_dic or mycontent_check', calling_function, [])     
    return all_is_same_list

def write_check(all_is_same_list, post_path = './y_dir/y_post_path.txt', mycontent_check = POSTPARAMSTA_CHECK, tag_split = '   ', tag_end = '\n\n', tag_begin='\n', left = False, num=13):
    # After loop, we have a all_is_same_list ( a list of 7 dictionary)
    # and_opera_check is a all True dictionary
    # is_common_dic saved combined results
    # create_content return a char "\n       Encut          500           OK\n", and append to mycontent
    # list_to_char to change to char, and write it
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    if mycontent_check and all_is_same_list:
        # sum all_is_same_list to is_common_dic
        is_common_dic = {}
        and_opera_check = {}
        mycontent = []
        # create and_opera_check all True
        for param_check in mycontent_check:
            for key, value in param_check.items():
                and_opera_check[key] = True
        # AND operation
        for key, value in and_opera_check.items():
            is_common = all(d.get(key) == value for d in all_is_same_list)
            if is_common:
                is_common_dic[key] = 'OK'
            else:
                is_common_dic[key] = 'NO'
        #print(is_common_dic)
        for param_check in mycontent_check:
            for key, value in param_check.items():
                content = create_content([key, value, is_common_dic.get(key)], tag_split, tag_end, tag_begin, left, num)
                mycontent.append(content)
        #print(mycontent)
        mycontent_char = list_to_char(mycontent, char_tag = '')
        #print(mycontent_char)
        write_content_to_file(mycontent_char, post_path, 'a')
    else:
        print_after_blank('mycontent_check or all_is_same_list', calling_function, [])    

def write_content_to_file(content_char, path, write_type='a'):
    # write content_char to file controled by write_type
    # for example, a-append, w-write overlap
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    try:
        with open(path, write_type) as post_file:
            post_file.write(content_char)
            #print(f"Content written to {path}")
    except FileNotFoundError:
        print(f"i can't open the path: '{path}'! - {calling_function}")

def find_line_position(file, search_string):
    # find the serch_string in file
    # return a list of position of search_string
    # noted, from 1, not 0. it's controled by start=1 or start=0
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)

    file.seek(0)
    line_number_list = []
    for line_number, line in enumerate(file, start=1):
        #print(line)
        if search_string in line:
            #print(line)
            line_number_list.append(line_number)
    file.seek(0)
    return line_number_list

def extract_line_at_position(file, line_number_list, move_list_down_up = [-10,5]):
    # input is a position list generated by find_line_position()
    # using my_down_up to enlarge the list controled by move_list_down_up (sorted)
    # extract the line in position of list
    calling_function = stack()[1].function
    check_input(getargvalues(currentframe()).locals)
    extract_lines = ''
    my_line_number_list = my_down_up(line_number_list, move_list_down_up)
    #print(my_line_number_list)
    for current_line_number, line in enumerate(file, start=1):
        if current_line_number in my_line_number_list:
            extract_lines = extract_lines + line
    file.seek(0)
    return extract_lines






# In[111]:


job = 0

if job:
    posttime = PostTime()
    posttime.read_OUTCAR()

    postdata = PostData()
    postdata.read_OUTCAR()

    postdata2 = PostData2()
    postdata2.read_OUTCAR()

    postdiff = PostDiff()
    postdiff.read_diff()

    postparam = PostParam()
    postparam.read_OUTCAR()

    postparam2 = PostParam2()
    postparam2.read_OUTCAR()
    
    postparamsta = PostParamSta()
    postparamsta.read_OUTCAR()

    postwarning = PostWarning()
    postwarning.read_OUTCAR()
else:
    postparamsta = PostParamSta()
    postparamsta.read_OUTCAR()

    postwarning = PostWarning()
    postwarning.read_OUTCAR()


# 

# 
