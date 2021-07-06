import numpy as np
import pandas as pd
import argparse 

part0 = np.loadtxt('/mnt/projects/lensing/HALO_SHAPE/MICEv1.0/catalogs/ind_halos/coords568416').T  
# part0 = np.loadtxt('/mnt/projects/lensing/HALO_SHAPE/MICEv1.0/catalogs/ind_halos/coords5707').T  

nrings = 10

xr,yr,zr = part0[3],part0[4],part0[5]
xp,yp = part0[8],part0[9]

rin = 0
step = pro[1]/float(nrings)

# s = 1.
# q = 1.
# q2 = 1.

rhop = np.zeros(nrings)
Sp = np.zeros(nrings)
rp = np.zeros(nrings)

Ntot = 0
for ring in range(nrings):
    
    abin_in = rin/(q*s)**(1./3.)
    bbin_in = abin_in*q
    cbin_in = abin_in*s

    abin_out = (rin+step)/(q*s)**(1./3.)
    bbin_out = abin_out*q
    cbin_out = abin_out*s
    
    # print(abin_out)
    
    rp[ring] = (rin + 0.5*step)/1.e3
    
    # rpart_E = (abc*(xr**2/abin**2 + yr**2/bbin**2 + zr**2/cbin**2))**(1./3.)/1.e3 
    rpart_E_in = (xr**2/abin_in**2 + yr**2/bbin_in**2 + zr**2/cbin_in**2)
    rpart_E_out = (xr**2/abin_out**2 + yr**2/bbin_out**2 + zr**2/cbin_out**2)
    
    V    = (4./3.)*np.pi*(((rin+step)/1.e3)**3 - (rin/1.e3)**3)
    mask = (rpart_E_in >= 1)*(rpart_E_out < 1)
    rhop[ring] = (mask.sum()*mp)/V

    # print(mask.sum())

    abin_in = rin/np.sqrt(q2) 
    bbin_in = abin_in*q2

    abin_out = (rin+step)/np.sqrt(q2) 
    bbin_out = abin_out*q2

    rpart_E_in = (xp**2/abin_in**2 + yp**2/bbin_in**2)
    rpart_E_out = (xp**2/abin_out**2 + yp**2/bbin_out**2)
    
    A    = np.pi*(((rin+step)/1.e3)**2 - (rin/1.e3)**2)
    mask = (rpart_E_in >= 1)*(rpart_E_out < 1)
    Sp[ring] = (mask.sum()*mp)/A

    rin += step
    
    


