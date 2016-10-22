import numpy as np
import matplotlib.pyplot as plt

numz=256

nsigma=6

ngals_c=4*np.pi*1E6
z0_c=0.5
dz_c=0.045
zarr_c=z0_c-nsigma*dz_c+2*nsigma*dz_c*np.arange(numz)/(numz-1.)
nzarr_c=np.exp(-0.5*((zarr_c-z0_c)/dz_c)**2)/np.sqrt(2*np.pi*dz_c**2)*ngals_c/(4*np.pi*(180*60/np.pi)**2)

ngals_s=4*np.pi*1E6
z0_s=0.65
dz_s=0.05
zarr_s=z0_s-nsigma*dz_s+2*nsigma*dz_s*np.arange(numz)/(numz-1.)
nzarr_s=np.exp(-0.5*((zarr_s-z0_s)/dz_s)**2)/np.sqrt(2*np.pi*dz_s**2)*ngals_s/(4*np.pi*(180*60/np.pi)**2)

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

stout ="# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]lmax\n"
stout+="%.3lf "%(z0_c-(nsigma-1)*dz_c)
stout+="%.3lf 0.001 0 500\n"%(z0_c+(nsigma-1)*dz_c)
f=open("bins_c.txt","w")
f.write(stout)
f.close()

stout ="# [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]lmax\n"
stout+="%.3lf "%(z0_s-(nsigma-1)*dz_s)
stout+="%.3lf 0.001 1 2000\n"%(z0_s+(nsigma-1)*dz_s)
f=open("bins_s.txt","w")
f.write(stout)
f.close()

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
