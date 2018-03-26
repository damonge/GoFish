import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

larr=np.logspace(np.log10(2),np.log10(2000),10).astype(int)

def nzpar(z,z0,a) :
    return z**2*np.exp(-(z/z0)**a)

ndens_c=20.; z0_c=0.13; a_c=0.78;
zarr_c=3*np.arange(256)/255.
nzarr_c=ndens_c*nzpar(zarr_c,z0_c,a_c)/quad(nzpar,0,3,args=(z0_c,a_c))[0]
zbins_c=np.array([[0.2,0.3],[0.4,0.5]])

ndens_s=10.; z0_s=0.13; a_s=0.78;
zarr_s=3*np.arange(256)/255.
nzarr_s=ndens_s*nzpar(zarr_s,z0_s,a_s)/quad(nzpar,0,3,args=(z0_s,a_s))[0]
zbins_s=np.array([[0.2,0.3],[0.4,0.5]])

zarr_b=np.array([0.0,0.3,0.5,0.65,0.9,1.6])
bbias_c=(1+zarr_b)**0.8
sbias_c=0.4*np.ones_like(zarr_b)
ebias_c=np.zeros_like(zarr_b)
abias_s=0.1*(1+zarr_b)
marg_bias_c=np.zeros_like(zarr_b); marg_bias_c[2]=1;
marg_bias_s=np.zeros_like(zarr_b); marg_bias_s[3]=1;
zarr_rf=1.6*np.arange(64)/63.
rf_arr=0.2/(1+(zarr_rf/0.8)**6)


np.savetxt("nz_c.txt",np.transpose([zarr_c,nzarr_c]),fmt='%lE %lE')
np.savetxt("nz_s.txt",np.transpose([zarr_s,nzarr_s]),fmt='%lE %lE')
np.savetxt("bz_c.txt",np.transpose([zarr_b,bbias_c,marg_bias_c]),fmt='%lE %lE %d -1')
np.savetxt("bz_s.txt",np.transpose([zarr_b,bbias_c,marg_bias_s]),fmt='%lE %lE %d -1')
np.savetxt("sz_c.txt",np.transpose([zarr_b,sbias_c,marg_bias_c]),fmt='%lE %lE %d -1')
np.savetxt("sz_s.txt",np.transpose([zarr_b,sbias_c,marg_bias_s]),fmt='%lE %lE %d -1')
np.savetxt("ez_c.txt",np.transpose([zarr_b,ebias_c,marg_bias_c]),fmt='%lE %lE %d -1')
np.savetxt("ez_s.txt",np.transpose([zarr_b,ebias_c,marg_bias_s]),fmt='%lE %lE %d -1')
np.savetxt("az_s.txt",np.transpose([zarr_b,abias_s,marg_bias_s]),fmt='%lE %lE %d -1')
np.savetxt("rf_s.txt",np.transpose([zarr_rf,rf_arr]),fmt='%lE %lE 0 -1')
np.savetxt("bins_c.txt",zbins_c,fmt='%.3lf %.3lf 0.01 0 0 2000',header='[0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax\n')
np.savetxt("bins_s.txt",zbins_s,fmt='%.3lf %.3lf 0.01 0 0 2000',header='[0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax\n')
np.savetxt("larr.txt",np.transpose([larr[:-1],larr[1:]]),fmt='%d %d')

plt.figure()
plt.plot(zarr_c,nzarr_c,'r-');
plt.plot(zarr_s,nzarr_s,'b-');
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$dN(z)/(dz\,d\\Omega)\\,[{\\rm amin}^{-2}]$',fontsize=16)
plt.xlim([zarr_c[0],zarr_s[-1]])

plt.figure()
plt.plot(zarr_b,bbias_c,'r-');
plt.plot(zarr_b,sbias_c,'g-');
plt.plot(zarr_b,ebias_c,'b-');
plt.plot(zarr_b,abias_s,'y-');
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$b_c(z)$',fontsize=16)
plt.xlim([0,zarr_s[-1]])

plt.figure()
plt.plot(zarr_b,abias_s,'y-');
plt.plot(zarr_rf,rf_arr,'k-');
plt.xlabel('$z$',fontsize=16)
plt.ylabel('$b_s(z)$',fontsize=16)
plt.xlim([0,zarr_s[-1]])

plt.show()
