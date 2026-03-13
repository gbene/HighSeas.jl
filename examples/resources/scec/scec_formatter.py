"""
Script used to decently format files from scec site. WTF is wrong with these people?
"""
import pandas as pd
import glob
from pathlib import Path
import sys 

input_folder = input("Input folder: ")

if not Path(input_folder).is_dir():
    print(f"{input_folder} does not exist")
    sys.exit() 
    

start_line = int(input("Data start line: "))-1
mode = int(input('Header type (1: onfault, 2: global, 3: rupture): '))

if mode == 1:
    cols= ["t","slip_2" ,"slip_3","slip_rate_2","slip_rate_3","shear_stress_2", "shear_stress_3", "state"]
    ext = ".txt"
elif mode == 2:
    cols= ["t", "max_slip_rate", "moment_rate"]
    ext = ".glb"
elif mode == 3:
    cols = ["x2", "x3", "t"]
    ext = ".rup"
file_names = glob.glob(f"{input_folder}/*{ext}")



for path in file_names:
    stupid_file_name = Path(path).stem

    badly_used_sep = " " # separator in the csv that i not properly separating the columns

    with open(path,'r') as f:
          lines = f.readlines()[start_line:]
          split = [l for l in lines[0].strip().split(badly_used_sep) if l !='']

          
          pd_dataframe = pd.DataFrame(columns=cols)
          for i, line in enumerate(lines):
                print(i)
                line = line.replace('\n','')
                split = [l for l in line.split(badly_used_sep) if l !='']
                if len(split) !=0:
                    pd_dataframe.loc[i] = pd.to_numeric(split) 



    pd_dataframe.to_csv(f'{input_folder}/{stupid_file_name}.csv',sep=' ', index=True)
