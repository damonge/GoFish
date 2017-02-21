import numpy as np
import py_cosmo_mad as csm

KMAX=0.2
pcs=csm.PcsPar()
pcs.set_verbosity(1)
pcs.background_set(0.315,0.685,0.049,-1.,0.,0.67,2.725)

def reprocess_file(fname,forgetlmax) :
    data=np.loadtxt(fname,unpack=True)
    z0=data[0]
    zf=data[1]
    sz=data[2]
    marg=data[3]
    zm=(z0+zf)/2
    if forgetlmax :
        lmax=2000*np.ones_like(sz)
    else :
        lmax=np.array([int(pcs.radial_comoving_distance(1./(1+z))*KMAX) for z in zm])
    np.savetxt(fname,np.transpose([z0,zf,sz,marg,lmax]),fmt='%lf %lf %lf %d %d',
               header='[0]z0 [1]zf [2]sigma_z [3]marg_sz [4]lmax')

reprocess_file('bins_red.txt',False)
reprocess_file('bins_blue.txt',False)
reprocess_file('bins_gold.txt',True)
