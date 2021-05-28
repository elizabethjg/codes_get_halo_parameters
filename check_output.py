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

for j in range(len(check_cat.columns)):
    col_check = check_cat.columns[j]
    col_cat   = cat.columns[j]
    m = np.abs(check_cat[col_check] > 0.) 
    print(col_cat,np.max(((cat[col_cat]/check_cat[col_check])[m])))

'''

cat = check_cat 

rc = np.sqrt((cat.xc_fof-cat.xc_rc)**2 + (cat.yc_fof-cat.yc_rc)**2 + (cat.zc_fof-cat.zc_rc)**2)

xold = cat.xc_fof - cat.xc_rc
yold = cat.yc_fof - cat.yc_rc
zold = cat.zc_fof - cat.zc_rc

lM = cat['log10(mass)']
idhalo = cat['Halo number']


part0 = np.loadtxt('ind_halos/particles_halo'+str(int(idpru))).T

f, ax = plt.subplots(1, 3, figsize=(12,4),sharex=True,sharey=True)
f.subplots_adjust(hspace=0,wspace=0)

ax[0].plot(part0[3],part0[4],'C7.')
ax[1].plot(part0[3],part0[5],'C7.')
ax[2].plot(part0[4],part0[5],'C7.')

ax[0].plot(part0[3],part0[4],'C7.')
ax[1].plot(part0[3],part0[5],'C7.')
ax[2].plot(part0[4],part0[5],'C7.')

ax[0].plot(0,0,'C3o')
ax[1].plot(0,0,'C3o')
ax[2].plot(0,0,'C3o')

ax[0].plot(xold[idhalo == idpru],yold[idhalo == idpru],'C1o')
ax[1].plot(xold[idhalo == idpru],zold[idhalo == idpru],'C1o')
ax[2].plot(yold[idhalo == idpru],zold[idhalo == idpru],'C1o')
'''
