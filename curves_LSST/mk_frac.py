import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import os as os

def reform_bins(fname_in,fname_out) :
    zm,dz,sz=np.loadtxt(fname_in,unpack=True)
    z0=zm-dz; zf=zm+dz;
    np.savetxt(fname_out,np.transpose([z0,zf,sz]))

reform_bins("/home/damonge/Science/PastProjects/2015_XlargeXcorr/run_LSST_redVSblue/bins/bins_gold_nw3.txt","bins_gold.txt")
reform_bins("/home/damonge/Science/PastProjects/2015_XlargeXcorr/run_LSST_redVSblue/bins/bins_blue_nw3.txt","bins_blue.txt")
reform_bins("/home/damonge/Science/PastProjects/2015_XlargeXcorr/run_LSST_redVSblue/bins/bins_red_nw3.txt" ,"bins_red.txt")


z_blue,nz_blue=np.loadtxt("nz_blue.txt",unpack=True)
nz_blue_i=interp1d(z_blue,nz_blue)
def nz_blue_f(z) :
    if z<=z_blue[0] :
        return 0
    elif z>=z_blue[-1] :
        return 0
    else :
        return nz_blue_i(z)

z_red,nz_red=np.loadtxt("nz_red.txt",unpack=True)
nz_red_i=interp1d(z_red,nz_red)
def nz_red_f(z) :
    if z<=z_red[0] :
        return 0
    elif z>=z_red[-1] :
        return 0
    else :
        return nz_red_i(z)

z_gold,nz_gold=np.loadtxt("nz_gold.txt",unpack=True)
nz_gold_i=interp1d(z_gold,nz_gold)
def nz_gold_f(z) :
    if z<=z_gold[0] :
        return 0
    elif z>=z_gold[-1] :
        return 0
    else :
        return nz_gold_i(z)

z_al,al_al=np.loadtxt("alignment.dat",unpack=True)
al_i=interp1d(z_al,al_al)
def al_f(z) :
    if z<=z_al[0] :
        return al_al[0]
    elif z>=z_al[-1] :
        return al_al[-1]
    else :
        return al_i(z)

def rfrac_f(z) :
    n_gold=nz_gold_f(z)
    n_red=nz_red_f(z)
    
    if n_gold<=0 :
        return 0
    else :
        return n_red/n_gold

z_nodes_al=np.array([0.0,0.5,1.0,1.5,2.0,3.0,4.0])
v_nodes_al=np.array([0,1,1,1,0,0,0])
al_nodes_al=np.array([al_f(z) for z in z_nodes_al])
abias=np.array([al_f(z) for z in z_gold])
plt.plot(z_nodes_al,al_nodes_al)
plt.plot(z_nodes_al,al_nodes_al,'ro')
plt.plot(z_gold,abias)
plt.show()
np.savetxt("dum.txt",np.transpose([z_gold,abias]))
os.system("awk '{print $1\" \"$2\" 0\"}' < dum.txt > az_gold.txt")
np.savetxt("dum.txt",np.transpose([z_nodes_al,al_nodes_al]))
os.system("awk '{print $1\" \"$2\" 0\"}' < dum.txt > az_gold.txt")

rfrac=np.array([rfrac_f(z) for z in z_gold])
plt.plot(z_gold,rfrac)
plt.show()
np.savetxt("dum.txt",np.transpose([z_gold,rfrac]))
os.system("awk '{print $1\" \"$2\" 0\"}' < dum.txt > rf_gold.txt")
