import numpy as np
import pyccl as ccl
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.optimize import minimize

#Generate array of multipoles
n_ell=25
d_ln_ell=np.log(5400./20.)/n_ell
ln_ells_0=np.log(20.)+d_ln_ell*np.arange(n_ell)
ln_ells_1=np.log(20.)+d_ln_ell*(np.arange(n_ell)+1)
ells_0=np.floor(np.exp(ln_ells_0))
ells_1=np.floor(np.exp(ln_ells_1))
np.savetxt("larr.txt",np.transpose([ells_0,ells_1]),fmt='%d %d')


cosmo=ccl.Cosmology(Omega_c=0.3156-0.0492,Omega_b=0.0492,h=0.6727,sigma8=0.831,n_s=0.9645)
z_arr=3.5*np.arange(512)/511.
pz_c_y10=z_arr**2*np.exp(-(z_arr/0.28)**0.9)
norm_c_y10=48./(np.sum(pz_c_y10)*(z_arr[1]-z_arr[0]))
nz_c_y10=norm_c_y10*pz_c_y10
pz_s_y10=z_arr**2*np.exp(-(z_arr/0.11)**0.68)
norm_s_y10=27./(np.sum(pz_s_y10)*(z_arr[1]-z_arr[0]))
nz_s_y10=norm_s_y10*pz_s_y10
np.savetxt("nz_clustering_y10.txt",np.transpose([z_arr,nz_c_y10]))
np.savetxt("nz_shear_y10.txt",np.transpose([z_arr,nz_s_y10]))

pz_c_y1=z_arr**2*np.exp(-(z_arr/0.26)**0.94)
norm_c_y1=18./(np.sum(pz_c_y1)*(z_arr[1]-z_arr[0]))
nz_c_y1=norm_c_y1*pz_c_y1
pz_s_y1=z_arr**2*np.exp(-(z_arr/0.13)**0.78)
norm_s_y1=10./(np.sum(pz_s_y1)*(z_arr[1]-z_arr[0]))
nz_s_y1=norm_s_y1*pz_s_y1
np.savetxt("nz_clustering_y1.txt",np.transpose([z_arr,nz_c_y1]))
np.savetxt("nz_shear_y1.txt",np.transpose([z_arr,nz_s_y1]))

zb_c_y10_0=0.2+0.1*np.arange(10)
zb_c_y10_1=0.2+0.1*(np.arange(10)+1)
sz_c_y10=0.03*(1+0.5*(zb_c_y10_0+zb_c_y10_1))
zbz_y10=[zb_c_y10_0[0]]
bz_y10=[0.95/ccl.growth_factor(cosmo,1./(1+zb_c_y10_0[0]))]
mask_y10=[0]
lmax_c_y10=[]
for z in 0.5*(zb_c_y10_0+zb_c_y10_1) :
    bz_y10.append(0.95/ccl.growth_factor(cosmo,1./(1+z)))
    zbz_y10.append(z)
    mask_y10.append(0)
    #mask_y10.append(1)
    lmax_c_y10.append(np.fmin(0.75*0.6727*ccl.comoving_radial_distance(cosmo,1./(1+z)),5400))
bz_y10.append(0.95/ccl.growth_factor(cosmo,1./(1+3.4)))
zbz_y10.append(3.4)
mask_y10.append(0)
bz_y10=np.array(bz_y10)
zbz_y10=np.array(zbz_y10)
mask_y10=np.array(mask_y10)
lmax_c_y10=np.array(lmax_c_y10)
np.savetxt("bz_clustering_y10.txt",np.transpose([zbz_y10,bz_y10,mask_y10]),fmt='%lf %lf %d -1')
np.savetxt("bins_clustering_y10.txt",np.transpose([zb_c_y10_0,zb_c_y10_1,sz_c_y10,lmax_c_y10]),
           fmt='%.5lf %.5lf %.5lf 0 0 %d',
           header=' [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax')

zb_c_y1_0=0.2+0.2*np.arange(5)
zb_c_y1_1=0.2+0.2*(np.arange(5)+1)
sz_c_y1=0.03*(1+0.5*(zb_c_y1_0+zb_c_y1_1))
zbz_y1=[zb_c_y1_0[0]]
bz_y1=[1.05/ccl.growth_factor(cosmo,1./(1+zb_c_y1_0[0]))]
mask_y1=[0]
lmax_c_y1=[]
for z in 0.5*(zb_c_y1_0+zb_c_y1_1) :
    bz_y1.append(1.05/ccl.growth_factor(cosmo,1./(1+z)))
    zbz_y1.append(z)
    mask_y1.append(0)
    #mask_y1.append(1)
    lmax_c_y1.append(np.fmin(0.75*0.6727*ccl.comoving_radial_distance(cosmo,1./(1+z)),5400))
bz_y1.append(1.05/ccl.growth_factor(cosmo,1./(1+3.4)))
zbz_y1.append(3.4)
mask_y1.append(0)
bz_y1=np.array(bz_y1)
zbz_y1=np.array(zbz_y1)
mask_y1=np.array(mask_y1)
lmax_c_y1=np.array(lmax_c_y1)
np.savetxt("bz_clustering_y1.txt",np.transpose([zbz_y1,bz_y1,mask_y1]),fmt='%lf %lf %d -1')
np.savetxt("bins_clustering_y1.txt",np.transpose([zb_c_y1_0,zb_c_y1_1,sz_c_y1,lmax_c_y1]),
           fmt='%.5lf %.5lf %.5lf 0 0 %d',
           header=' [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax')

nzf_y10=interp1d(z_arr,nz_s_y10,fill_value=0,bounds_error=False)
znzf_y10=interp1d(z_arr,z_arr*nz_s_y10,fill_value=0,bounds_error=False)
nzt_y10=quad(nzf_y10,0,3.5)[0]
zb_s_y10_0=np.zeros(10); zb_s_y10_1=np.zeros(10); sz_s_y10=np.zeros(10); lmax_s_y10=np.zeros(10)
for i in np.arange(10) :
    def func(z) :
        nz=quad(nzf_y10,zb_s_y10_0[i],z)[0]
        return (nz-nzt_y10/10.)**2
    res=minimize(func,zb_s_y10_0[i]+0.05,method='Powell')
    if i<9 :
        zb_s_y10_0[i+1]=np.fmin(res.x,3.4)
    zb_s_y10_1[i]=np.fmin(res.x,3.4)
    zm=quad(znzf_y10,zb_s_y10_0[i],zb_s_y10_1[i])[0]/quad(nzf_y10,zb_s_y10_0[i],zb_s_y10_1[i])[0]
    sz_s_y10[i]=0.05*(1+zm)
    lmax_s_y10[i]=np.fmin(1.5*0.6727*ccl.comoving_radial_distance(cosmo,1./(1+zm)),5400)
    print i,zb_s_y10_0[i],zb_s_y10_1[i],zm,func(zb_s_y10_1[i])
np.savetxt("bins_shear_y10.txt",np.transpose([zb_s_y10_0,zb_s_y10_1,sz_s_y10,lmax_s_y10]),
           fmt='%.5lf %.5lf %.5lf 0 0 %d',
           header=' [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax')

nzf_y1=interp1d(z_arr,nz_s_y1,fill_value=0,bounds_error=False)
znzf_y1=interp1d(z_arr,z_arr*nz_s_y1,fill_value=0,bounds_error=False)
nzt_y1=quad(nzf_y1,0,3.5)[0]
zb_s_y1_0=np.zeros(5); zb_s_y1_1=np.zeros(5); sz_s_y1=np.zeros(5); lmax_s_y1=np.zeros(5)
for i in np.arange(5) :
    def func(z) :
        nz=quad(nzf_y1,zb_s_y1_0[i],z)[0]
        return (nz-nzt_y1/5.)**2
    res=minimize(func,zb_s_y1_0[i]+0.05,method='Powell')
    if i<4 :
        zb_s_y1_0[i+1]=np.fmin(res.x,3.4)
    zb_s_y1_1[i]=np.fmin(res.x,3.4)
    zm=quad(znzf_y1,zb_s_y1_0[i],zb_s_y1_1[i])[0]/quad(nzf_y1,zb_s_y1_0[i],zb_s_y1_1[i])[0]
    sz_s_y1[i]=0.05*(1+zm)
    lmax_s_y1[i]=np.fmin(1.5*0.6727*ccl.comoving_radial_distance(cosmo,1./(1+zm)),5400)
    print i,zb_s_y1_0[i],zb_s_y1_1[i],zm,func(zb_s_y1_1[i])
np.savetxt("bins_shear_y1.txt",np.transpose([zb_s_y1_0,zb_s_y1_1,sz_s_y1,lmax_s_y1]),
           fmt='%.5lf %.5lf %.5lf 0 0 %d',
           header=' [0]z0 [1]zf [2]sigma_z [3]marg_sz [4]marg_bz [5]lmax')
