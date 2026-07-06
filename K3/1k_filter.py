# to get rid of admix groups. 
# will create a file for BCFtools to use to filter out these samples 

import sys
import pandas as pd
from tabulate import tabulate

# prior used  to get sample list 
file_samples="/path/to/data/K3/samples.txt"
pop_file="/path/to/shakira/resources/1K_pops.txt"
exclude="/path/to/data/K3/noAsia.txt"

#ASW excluded due to admixture
#ACB excluded due to admixture
#BEB - Asian
#CDX - Asian
#FIN - included but I want to circle back due to known asian Ancestry in Finish populations 
#GIH -Asian
#CHB - Asian
#CHS asian
#ITU SAS
# JPT
# KHV
# MXL - admixture
# PJL - asian
#STU - Asian 


pops_to_exclude=['MXL','ACB','ASW', 'BEB', 'CDX', 'GIH', 'CHB', 'CHS', "ITU", "JPT", 'KHV', 'PJL', 'STU']

samples=open(file_samples,"r")

final=open(exclude,"w")

table_pop=pd.read_csv(pop_file,'\s+',header='infer')

print(table_pop['SAMPLE_NAME'])
for line in samples:
    if line.strip("\n") in table_pop['SAMPLE_NAME'].tolist(): #to do save this list later as var prior very bad rn
#        print("sample " + line)
        pop_row= table_pop[table_pop['SAMPLE_NAME']==line.strip('\n')]
        pop=pop_row['POPULATION'].iloc[0]
        if pop in pops_to_exclude:
            final.write("0 \t" + line)