import sys, os
sys.path.append("..")
import numpy as np
# import matplotlib.pyplot as plt
import common as com
import sys as sys
import os as os
import fisher_plot as fsh
import argparse

parser = argparse.ArgumentParser()

parser.add_argument("--red",
                  action="store", dest="red", default=1,type=int)
parser.add_argument("--vary_baryons",
                  action="store", dest="vary_baryons", default=1,type=int)
parser.add_argument("--m_bias",
                  action="store", dest="m_bias", default=1,type=int)
parser.add_argument("--photoz",
                  action="store", dest="photoz", default=1,type=int)
parser.add_argument("--shear_lmax5000",
                  action="store", dest="shear_lmax5000", default=1,type=int)

results = parser.parse_args()
red=bool(results.red)
vary_baryons=bool(results.vary_baryons)
m_bias=bool(results.m_bias)
photoz=bool(results.photoz)
shear_lmax5000=bool(results.shear_lmax5000)

plots_dir = "/group/hepheno/smsharma/GoFish/plots/"

param_och2 = """
#Cosmological parameters
[och2]
#Fiducial value
x= 0.1197
#Increment used for numerical derivatives
dx= 0.001
#Set to 'no' if this parameter should be keep fixed
is_free= yes
onesided=0

"""

param_obh2 = """
[obh2]
x= 0.02222
dx= 0.0001
is_free= yes
onesided=0

"""

param_ok = """
[ok]
x= 0.0
dx= 0.015
is_free= yes
onesided=0

"""

param_hh = """
[hh]
x= 0.69
dx= 0.01
is_free= yes
onesided=0

"""

param_w0 = """
[w0]
x= -1.0
dx= 0.01
is_free= yes
onesided=0

"""

param_wa = """
[wa]
x=0.0
dx=0.05
is_free= yes
onesided=0

"""

param_As = """
[A_s]
x= 2.1955
dx= 0.01
is_free= yes
onesided=0

"""

param_ns = """
[ns]
x= 0.9655
dx= 0.005
is_free= yes
onesided=0

"""

param_tau = """
[tau]
x=0.06
dx=0.02
is_free= yes
onesided = 0

"""

param_rt = """
[rt]
x=0.00
dx=0.005
is_free= no
onesided = 1

"""

param_mnu = """
[mnu]
x=60.
dx=10.
is_free= yes
onesided = 0

"""

param_pan = """
[pan]
x=0.
dx=0.02
is_free= no
onesided = 1

"""

param_lmcb = """
[lmcb]
x=14.08
dx=0.3
is_free= yes
onesided = 0

"""

param_etab = """
[etab]
x=0.5
dx=0.1
is_free= yes
onesided = 0

"""

tracer_cmb_primary_nos4 = """
[Tracer 1]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 31.13 31.13
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 200.94 56.9402
#Beam size in units of arcmin in the two ell regimes
beam_amin= 10.0 7.0 
#Ell value that marks the transition
l_transition= 50
#Minimum ell
lmin= 2
#Maximum ell
lmax= 3000 
use_tracer= yes

"""

tracer_cmb_primary_cv = """
[Tracer 1]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 31.13 1.0
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 4.0 1.4
#Beam size in units of arcmin in the two ell regimes
beam_amin= 1.6 3.0 
#Ell value that marks the transition
l_transition= 50
#Minimum ell
lmin= 2
#Maximum ell
lmax= 3000 
use_tracer= yes

"""

tracer_cmb_primary = """
[Tracer 1]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 31.13 1.0
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 150.0 1.4
#Beam size in units of arcmin in the two ell regimes
beam_amin= 10.0 3.0 
#Ell value that marks the transition
l_transition= 50
#Minimum ell
lmin= 2
#Maximum ell
lmax= 3000 
use_tracer= yes

"""

tracer_cmb_lensing = """
[Tracer 2]
tracer_name= CMB_exp_lensing
tracer_type= cmb_lensing
#Noise level in temperature in units of uK-arcmin (implicitly scaled by sqrt(2) in polarization)
sigma_t= 1.0
#Beam size in units of arcmin
beam_amin= 3.0
lmin= 30
lmax= 3000
use_tracer= yes

"""

tracer_lensing = """
[Tracer 3]
tracer_name= LSST_gold_sh
tracer_type= gal_shear
bins_file= ../curves_LSST/bins_gold.txt
nz_file= ../curves_LSST/nz_shear_fiducial.txt
abias_file= ../curves_LSST/az_gold.txt
rfrac_file= ../curves_LSST/rf_gold.txt
sigma_gamma= 0.28
include_m_bias= yes
m_step = 0.005
use_tracer= yes

"""

tracer_clustering_blue = """
[Tracer 4]
tracer_name= LSST_blue_cl
tracer_type= gal_clustering
bins_file= ../curves_LSST/bins_blue.txt
nz_file= ../curves_LSST/nz_blue.txt
bias_file= ../curves_LSST/bz_blue.txt
sbias_file= ../curves_LSST/sz_blue.txt
ebias_file= ../curves_LSST/ez_blue.txt
use_tracer= yes

"""

tracer_clustering_red = """
[Tracer 5]
tracer_name= LSST_red_cl
tracer_type= gal_clustering
bins_file= ../curves_LSST/bins_red.txt
nz_file= ../curves_LSST/nz_red.txt
bias_file= ../curves_LSST/bz_red.txt
sbias_file= ../curves_LSST/sz_red.txt
ebias_file= ../curves_LSST/ez_red.txt
use_tracer= yes

"""

tracer_bao = """
[BAO 1]
# Absolute errors on dV, dA, H
#fname_dv=../curves_DESI/DESI_sigma_dV.txt
fname_da=../curves_DESI/DESI_sigma_dA.txt
fname_hh=../curves_DESI/DESI_sigma_H.txt
#use_relative_errors = yes

"""

class_params = """
[CLASS parameters]
lmax_cmb= 5000
lmax_lss= 5000
lmin_limber= 200
include_alignment= yes
include_rsd= yes
include_magnification= no
include_gr_vel= no
include_gr_pot= no
exec_path=./run_class.sh
use_nonlinear= yes
use_baryons= yes
f_sky= 0.4

"""

output_params = """
[Output parameters]
output_dir= ../outputs_LSST_DESI_lmax5000
output_spectra= run
output_fisher= Fisher

[Behaviour parameters]
model= wCDM
save_cl_files= yes
save_param_files= yes

"""

def get_ini(setup = "fid", desi=False, shear=True, red=False, blue=True, ok=True, w=True, m_bias=True, vary_baryons=True, photoz=True, shear_lmax5000=True):
    """ Get param.ini file for a given spec as defined in the args
    """
    
    param_ini = (param_och2 + param_obh2 + param_ok +
        param_hh + param_w0 + param_wa + 
        param_As + param_ns + param_tau + 
        param_rt + param_mnu + param_pan 
        )
    
    param_ini += (param_lmcb + param_etab).replace("is_free= yes", "is_free= yes" if vary_baryons else "is_free= no")

    if setup == "fid":
        tracers_ini = tracer_cmb_lensing + tracer_cmb_primary
    elif setup == "nos4":
        tracers_ini = tracer_cmb_lensing.replace("is_free= yes", "is_free= no") + tracer_cmb_primary_nos4
    elif setup == "cv":
        tracers_ini = tracer_cmb_lensing + tracer_cmb_primary_cv


    
    if desi: tracers_ini += tracer_bao
    tracers_ini += tracer_lensing.replace("use_tracer= yes", "use_tracer= yes" if shear else "use_tracer= no").replace("include_m_bias= yes", "include_m_bias= yes" if m_bias else "include_m_bias= no").replace("bins_gold.txt", "bins_gold.txt" if shear_lmax5000 else "bins_gold_lmax2000.txt")
    tracers_ini += tracer_clustering_red.replace("use_tracer= yes", "use_tracer= yes" if red else "use_tracer= no").replace("bins_blue.txt", "bins_blue.txt" if photoz else "bins_blue_nomarg.txt")
    tracers_ini += tracer_clustering_blue.replace("use_tracer= yes", "use_tracer= yes" if blue else "use_tracer= no").replace("bins_red.txt", "bins_red.txt" if photoz else "bins_red_nomarg.txt")
    tracers_ini += param_ok.replace("is_free= yes", "is_free= yes" if ok else "is_free= no")
    tracers_ini += (param_w0 + param_wa).replace("is_free= yes", "is_free= yes" if w else "is_free= no")
    
    ini_file = param_ini + tracers_ini + class_params + output_params.replace("Fisher", "FisherDella"+str(setup)+str(desi)+str(shear)+str(red)+str(blue)+str(ok)+str(w)+str(m_bias)+str(vary_baryons)+str(photoz) +str(shear_lmax5000))
    ini_name = "param"+str(setup)+str(desi)+str(shear)+str(red)+str(blue)+str(ok)+str(w)+str(m_bias)+str(vary_baryons) +str(photoz) +str(shear_lmax5000)+".ini"

    return ini_file, ini_name

parlist = []
parnames = []

fstr, fname = get_ini(desi=False,shear=True,red=red,blue=True,ok=True,w=True,m_bias=m_bias, vary_baryons=vary_baryons, shear_lmax5000=shear_lmax5000,photoz=photoz)
parlist.append(fname)
parnames.append("Marginalize over baryons")
f=open(fname, "w")
f.write(fstr)
f.close()

pars=[com.ParamRun(parname) for parname in parlist]

for i in range(len(parlist)):
#     pars[i] = com.ParamRun(parlist[i])
    print " "
    if (not os.path.isfile(pars[i].output_dir+"/"+pars[i].output_fisher+"/fisher_raw.npz")) :
        print "<> Computing/reading relevant signal power spectra"
        pars[i].get_cls_all()

        if pars[i].just_run_cls==False :
            print "<> Computing relevant noise power spectra"
            pars[i].get_cls_noise()
            print " "
#             pars[i].plot_cls_params()
    print " "
    pars[i].get_fisher_cls()
    pars[i].get_fisher_bao()
    pars[i].join_fishers()
    pars[i].plot_fisher()
