import numpy as np
import pandas as pd
import argparse 

parser = argparse.ArgumentParser()
parser.add_argument('-file_name', action='store', dest='file_name',default='output')
args = parser.parse_args()

check_cat = pd.read_csv('output_check.csv.bz2')
cat    = pd.read_csv(args.file_name+'.bz2')

print(check_cat.shape)
print(cat.shape)

for j in cat.columns: 
    m = np.abs(check_cat[j] > 0.) 
    print(j,np.max(((cat[j]/check_cat[j])[m])))