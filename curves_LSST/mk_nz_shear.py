import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad

case="optimistic"

def mk_nz(case) :
    if case=="optimistic" :
        neff=31.
        alpha=1.23
        beta=1.05
        z0=0.59
    elif case=="fiducial" :
        neff=26.
        alpha=1.24
        beta=1.01
        z0=0.51
    elif case=="pessimistic" :
        neff=18.
        alpha=1.28
        beta=0.97
        z0=0.41
        
    numz=256

    #N(z)
    zarr=4.0*np.arange(numz)/(numz-1.)
    pzarr=zarr**alpha*np.exp(-(zarr/z0)**beta)
    pzi=interp1d(zarr,pzarr)
    norm=quad(pzi,0,3)[0]
    nzarr=pzarr*neff/norm
    nzi=interp1d(zarr,nzarr)
    print quad(nzi,0,3)[0]

    #FR(z)
    zred,nzred=np.loadtxt("nz_red.txt",unpack=True)
    nzri=interp1d(zred,nzred)
    def nzrp(z) :
        if z<zred[0] :
            return 0
        elif z>=zred[-1] :
            return 0
        else :
            return nzri(z)
    fred=np.array([nzrp(z)/nzi(z) for z in zarr])
    fred[0]=0
    plt.title(case); plt.plot(zarr,fred); plt.show()

    np.savetxt("nz_shear_"+case+".txt",np.transpose([zarr,nzarr]))
    np.savetxt("rf_shear_"+case+".txt",np.transpose([zarr,fred,np.zeros(len(fred))]),fmt='%E %E %d')

mk_nz("optimistic")
mk_nz("fiducial")
mk_nz("pessimistic")
