get_fisher_cls_vec.py                                                                               0000644 0001750 0001750 00000012363 13267425643 015074  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                    def get_fisher_cls_vec(self) :
        """ Compute Fisher matrix from numerical derivatives """

        if self.n_tracers<=0 :
            return

        fname_save=self.output_dir+"/"+self.output_fisher+"/fisher_raw.npz"
        if os.path.isfile(fname_save) : # Fisher matrices have already been computed
            self.fshr_l=np.load(fname_save)['fisher_l'] # no computation is performed,
            self.fshr_cls=np.load(fname_save)['fisher_cls'] # fshr is taken from the computed .npz file
        else : # computation is necessary (no .npz file found)

            lmax_arr=np.zeros(self.nbins_total) # vectors with length equal number of bins,
            lmin_arr=np.zeros(self.nbins_total) # will contain l-intervals
            nbt=0 # number(bins)*number(tracers)
            for tr in self.tracers : # iterate over the tracers
                if tr.consider_tracer :
                    zarr=None
                    if ((tr.tracer_type=='gal_clustering') or
                        (tr.tracer_type=='intensity_mapping') or
                        (tr.tracer_type=='gal_shear')) :
                        data=np.loadtxt(tr.bins_file,unpack=True)
                        zarr=(data[0]+data[1])/2 # redshifts are taken as the averages of the given intervals
                    for ib in np.arange(tr.nbins)  : # iterate over bins

                        if zarr is not None : # What is this???
                            lmn=tr.lmin
                        else :
                            lmn=tr.lmin

                        lmax_arr[nbt]=min(tr.lmax_bins[ib],tr.lmax)
                        lmin_arr[nbt]=lmn
                        nbt+=1



            # Number of bins for each tracer
            nbins_lft = 6
            # define a weighting vector, currently contains only 0 and 1
            weight_lft = np.ones(nbins_lft*(2*nbins_lft+1))

            # find the startpoints of the individiual diagonals
            startpoints_lft = np.zeros(nbins_lft*2, dtype=int)
            for ii in np.arange(1, nbins_lft*2):
                startpoints_lft[ii] = startpoints_lft[ii-1] + (nbins_lft*2 - ii + 1)
                        
            # (1) Cut the dndn part 
            dndn_scheme_lft = 'A'
            # A ... take only leading diagonal
            # B ... take all
            if dndn_scheme_lft == 'A':
                for ii in np.arange(1, nbins_lft):
                    weight_lft[ startpoints_lft[ii] : startpoints_lft[ii] + nbins_lft - ii ] = 0.

            # (2) Cut the e-dn part
            edn_scheme_lft = 'a'
            nrows_edn_lft = 2
            # a ... take all
            # b ... remove lower triangle
            # c ... remove b + main diagonal
            # d ... remove c + nrows_edn_lft
            if edn_scheme_lft == 'b' or edn_scheme_lft == 'c' or edn_scheme_lft == 'd':
                for ii in np.arange(1, nbins_lft - 1):
                    weight_lft[ startpoints_lft[ii] + nbins_lft - ii : startpoints_lft[ii] + nbins_lft   ] = 0.
            if edn_scheme_lft == 'c' or edn_scheme_lft == 'd':
                weight_lft[ startpoints_lft[nbins_lft] : startpoints_lft[nbins_lft + 1] ] = 0.
            if edn_scheme_lft == 'd':
                for ii in np.arange(0, nrows_edn_lft):
                weight_lft[ startpoints_lft[nbins_lft + ii + 2] - (nrows_edn_lft - ii) : startpoints_lft[nbins_lft + ii + 2]  ] = 0.



            for il,l in enumerate(self.larr) :
                dl_bpw=self.d_larr[il]
                if l==0 :
                    continue
                ell=float(l)
                fish_pre=dl_bpw*self.fsky*(2*ell+1.) # normalization prefactor
                indices=np.where((lmin_arr<=l) & (lmax_arr>=l))[0]
                gst=GaussSt(len(indices))
                cl_fid_mat=(self.cl_fid_arr[il,indices,:][:,indices]+
                            self.cl_noise_arr[il,indices,:][:,indices])
                covar_lft = gst.gaussian_covariance(cl_fid_mat,1./fish_pre)

                # Modify the covariance matrix: wherever weight_lft[i] = 0., set ith column and row to zero, but retain unity on diagonal
                for ii in np.arange(0, len(weight_lft)):
                    if weight_lft[ii] < 0.5:
                        covar_lft[ii,:] = 0.
                        covar_lft[:,ii] = 0.
                        covar_flt[ii, ii] = 1.

                i_covar=np.linalg.inv(covar_lft)
                for i in np.arange(self.npar_vary+self.npar_mbias) :
                    dcl1=gst.unwrap_matrix(self.dcl_arr[i,il,indices,:][:,indices])
		            # convert matrix of the derivatives of the cl's into a vector
                    for j in np.arange(self.npar_vary-i+self.npar_mbias)+i :
                        dcl2=gst.unwrap_matrix(self.dcl_arr[j,il,indices,:][:,indices])
                        

                        # Multiply the datavectors dcl1, dcl2 with the weights
                        for ii in np.arange(len(dcl1)):
                            dcl1[ii] *= weight_lft[ii]
                            dcl2[ii] *= weight_lft[ii]

                        self.fshr_l[i,j,il]=np.dot(dcl1,np.dot(i_covar,dcl2))
                        # eq 28 of Duncan et al 2013
                        if i!=j :
                            self.fshr_l[j,i,il]=self.fshr_l[i,j,il]

            self.fshr_cls[:,:]=np.sum(self.fshr_l,axis=2)
                                                                                                                                                                                                                                                                             param_may14_magnification_no.ini                                                                    0000644 0001750 0001750 00000011646 13276365456 017115  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                #Cosmological parameters
[och2]
#Fiducial value
x= 0.1197
#Increment used for numerical derivatives
dx= 0.001
#Set to 'no' if this parameter should be keep fixed
is_free= no
onesided=0

[obh2]
x= 0.02222
dx= 0.0001
is_free= no
onesided=0

#Note om, ob and s8 are only used if use_cmb_params= no
[om]
x= 0.3156
dx= 0.002
is_free= yes
onesided=0

# Added this - LFT
[ok]
x= 0.0
dx= 0.002
is_free= yes
onesided=0

[ob]
x= 0.0492
dx= 0.0002
is_free= yes
onesided=0

[hh]
x= 0.69
dx= 0.01
is_free= yes
onesided=0

[w0]
x= -1.0
dx= 0.01
is_free= yes
onesided=0

[wa]
x=0.0
dx=0.05
is_free = no
onesided=0

[s8]
x= 0.831
dx= 0.01
is_free= yes
onesided=0

[A_s]
x= 2.1955
dx= 0.01
is_free= no
onesided=0

[ns]
x= 0.9655
dx= 0.005
is_free= yes
onesided=0

[tau]
x=0.06
dx=0.01
is_free = no
onesided = 0

[rt]
x=0.00
dx=0.005
is_free = no
onesided = 1

[mnu]
x=60.
dx=10.
is_free = no
onesided = 0

[pan]
x=0.
dx=0.02
is_free = no
onesided = 1

[Tracer 1]
tracer_name= Galaxy_survey
tracer_type= gal_clustering
#File describing the redshift bins for this tracer
bins_file= scripts_23Nov/bins_c.txt
#Overall true-redshift distribution
nz_file= scripts_23Nov/nz_c.txt
#Clustering bias
bias_file= scripts_23Nov/bz_c.txt
#Magnification bias
sbias_file= scripts_23Nov/sz_c.txt
#Evolution bias
ebias_file= irrelevant
use_tracer= yes

[Tracer 2]
tracer_name= Shear_survey
tracer_type= gal_shear
bins_file= scripts_23Nov/bins_s.txt
nz_file= scripts_23Nov/nz_s.txt
#Alignment bias (only needed if including IA)
abias_file= irrelevant 
#Red (aligned) fraction (only needed if including IA)
rfrac_file= irrelevant
#Intrinsic ellipticity dispertion
sigma_gamma= 0.28
use_tracer= yes

[Tracer 3]
tracer_name= Shear_survey
tracer_type= gal_clustering
bins_file= scripts_23Nov/bins_s.txt
nz_file= scripts_23Nov/nz_s.txt
bias_file= scripts_23Nov/bz_s.txt
sbias_file= scripts_23Nov/sz_s.txt
ebias_file= irrelevant
use_tracer= yes

[Tracer 0]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 1.0 1.0
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 1.4 1.4
#Beam size in units of arcmin in the two ell regimes
beam_amin= 3.0 3.0 
#Ell value that marks the transition
l_transition= 2 
#Minimum ell
lmin= 30 
#Maximum ell
lmax= 5000 
use_tracer= yes

[Tracer 0]
tracer_name= CMB_exp_lensing
tracer_type= cmb_lensing
#Noise level in temperature in units of uK-arcmin (implicitly scaled by sqrt(2) in polarization)
sigma_t= 1.0
#Beam size in units of arcmin
beam_amin= 3.0
lmin= 30
lmax= 3000
use_tracer= yes

[CLASS parameters]
#File with bandpowers. Two columns corresponding to the l_min and l_max(inclusive) of each bandpower
bandpowers_file= none
#Maximum multipole for which the CMB Cls will be computed
#(only relevant if CMB primary or CMB lensing are included)
lmax_cmb= 5000
#Maximum multipole for which the clustering and lensing Cls will be computed
#(only relevant if galaxy clustering or shear are included)
lmax_lss= 5000
#Minimum multipole from which Limber will be used
lmin_limber= 0
#Include intrinsic alignments in the shear power spectra?
include_alignment= no
#Include RSDs in the clustering power spectra?
include_rsd= no
#Include lensing magnification in the clustering power spectra?
include_magnification= no
#Include relativistic effects in the clustering power spectra?
include_gr_vel= no
include_gr_pot= no
#Command used to execute classt (it should take the CLASS param file as an argument)
exec_path= ./class_mod
#Use non-linear matter transfer function (HALOFit)?
use_nonlinear= yes
#Set to "yes" if you want to include baryonic effects in the power spectrum
use_baryons= no
#Sky fraction (20,000 sq deg for LSST)
f_sky= 0.37

[Cutting parameters]
# Number of bins for each tracer (should be identical)
nbins_lft= 6
# Cutting scheme for the dn-dn part
# A: take only leading diagonal
# B: take all
dndn_scheme_lft= A
# Cutting scheme for e-dn part
# a: take all
# b: remove lower triangle
# c: b + remove leading diagonal
# d: c + remove nrows_edn_lft
edn_scheme_lft= a
nrows_edn_lft= 0



[Output parameters]
#Directory where all the data will be output
output_dir= outputs_may14_magnification_no
#Prefix of all power-spectrum-related output files
output_spectra= run
#Directory where the Fisher information and plots will be stored
output_fisher= Fisher

[Behaviour parameters]
#Cosmological  model
model= wCDM
#Use [Omega_c*h^2,Omega_b*h^2,A_s]? If no, the code will use [Omega_m,Omega_b,sigma_8]
use_cmb_params= no
#Do you wanna keep the Cl files?
save_cl_files= yes
#Do you wanna keep the CLASS param files?
save_param_files= yes
#Do you want to just generate the power spectra and not the Fisher matrix?
just_run_cls= no
#File containing a "measured" power spectrum to estimate the associated parameter bias
bias_file= outputs_may14_magnification_no/run_fidcl.dat
#Do you want to write the Fisher matrix in terms of power spectra?
do_pspec_fisher= yes
                                                                                          param_may14_magnification_yes.ini                                                                   0000644 0001750 0001750 00000011600 13276364571 017264  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                #Cosmological parameters
[och2]
#Fiducial value
x= 0.1197
#Increment used for numerical derivatives
dx= 0.001
#Set to 'no' if this parameter should be keep fixed
is_free= no
onesided=0

[obh2]
x= 0.02222
dx= 0.0001
is_free= no
onesided=0

#Note om, ob and s8 are only used if use_cmb_params= no
[om]
x= 0.3156
dx= 0.002
is_free= yes
onesided=0

# Added this - LFT
[ok]
x= 0.0
dx= 0.002
is_free= yes
onesided=0

[ob]
x= 0.0492
dx= 0.0002
is_free= yes
onesided=0

[hh]
x= 0.69
dx= 0.01
is_free= yes
onesided=0

[w0]
x= -1.0
dx= 0.01
is_free= yes
onesided=0

[wa]
x=0.0
dx=0.05
is_free = no
onesided=0

[s8]
x= 0.831
dx= 0.01
is_free= yes
onesided=0

[A_s]
x= 2.1955
dx= 0.01
is_free= no
onesided=0

[ns]
x= 0.9655
dx= 0.005
is_free= yes
onesided=0

[tau]
x=0.06
dx=0.01
is_free = no
onesided = 0

[rt]
x=0.00
dx=0.005
is_free = no
onesided = 1

[mnu]
x=60.
dx=10.
is_free = no
onesided = 0

[pan]
x=0.
dx=0.02
is_free = no
onesided = 1

[Tracer 1]
tracer_name= Galaxy_survey
tracer_type= gal_clustering
#File describing the redshift bins for this tracer
bins_file= scripts_23Nov/bins_c.txt
#Overall true-redshift distribution
nz_file= scripts_23Nov/nz_c.txt
#Clustering bias
bias_file= scripts_23Nov/bz_c.txt
#Magnification bias
sbias_file= scripts_23Nov/sz_c.txt
#Evolution bias
ebias_file= irrelevant
use_tracer= yes

[Tracer 2]
tracer_name= Shear_survey
tracer_type= gal_shear
bins_file= scripts_23Nov/bins_s.txt
nz_file= scripts_23Nov/nz_s.txt
#Alignment bias (only needed if including IA)
abias_file= irrelevant 
#Red (aligned) fraction (only needed if including IA)
rfrac_file= irrelevant
#Intrinsic ellipticity dispertion
sigma_gamma= 0.28
use_tracer= yes

[Tracer 3]
tracer_name= Shear_survey
tracer_type= gal_clustering
bins_file= scripts_23Nov/bins_s.txt
nz_file= scripts_23Nov/nz_s.txt
bias_file= scripts_23Nov/bz_s.txt
sbias_file= scripts_23Nov/sz_s.txt
ebias_file= irrelevant
use_tracer= yes

[Tracer 0]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 1.0 1.0
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 1.4 1.4
#Beam size in units of arcmin in the two ell regimes
beam_amin= 3.0 3.0 
#Ell value that marks the transition
l_transition= 2 
#Minimum ell
lmin= 30 
#Maximum ell
lmax= 5000 
use_tracer= yes

[Tracer 0]
tracer_name= CMB_exp_lensing
tracer_type= cmb_lensing
#Noise level in temperature in units of uK-arcmin (implicitly scaled by sqrt(2) in polarization)
sigma_t= 1.0
#Beam size in units of arcmin
beam_amin= 3.0
lmin= 30
lmax= 3000
use_tracer= yes

[CLASS parameters]
#File with bandpowers. Two columns corresponding to the l_min and l_max(inclusive) of each bandpower
bandpowers_file= none
#Maximum multipole for which the CMB Cls will be computed
#(only relevant if CMB primary or CMB lensing are included)
lmax_cmb= 5000
#Maximum multipole for which the clustering and lensing Cls will be computed
#(only relevant if galaxy clustering or shear are included)
lmax_lss= 5000
#Minimum multipole from which Limber will be used
lmin_limber= 0
#Include intrinsic alignments in the shear power spectra?
include_alignment= no
#Include RSDs in the clustering power spectra?
include_rsd= no
#Include lensing magnification in the clustering power spectra?
include_magnification= yes
#Include relativistic effects in the clustering power spectra?
include_gr_vel= no
include_gr_pot= no
#Command used to execute classt (it should take the CLASS param file as an argument)
exec_path= ./class_mod
#Use non-linear matter transfer function (HALOFit)?
use_nonlinear= yes
#Set to "yes" if you want to include baryonic effects in the power spectrum
use_baryons= no
#Sky fraction (20,000 sq deg for LSST)
f_sky= 0.37

[Cutting parameters]
# Number of bins for each tracer (should be identical)
nbins_lft= 6
# Cutting scheme for the dn-dn part
# A: take only leading diagonal
# B: take all
dndn_scheme_lft= A
# Cutting scheme for e-dn part
# a: take all
# b: remove lower triangle
# c: b + remove leading diagonal
# d: c + remove nrows_edn_lft
edn_scheme_lft= a
nrows_edn_lft= 0



[Output parameters]
#Directory where all the data will be output
output_dir= outputs_may14_magnification_yes
#Prefix of all power-spectrum-related output files
output_spectra= run
#Directory where the Fisher information and plots will be stored
output_fisher= Fisher

[Behaviour parameters]
#Cosmological  model
model= wCDM
#Use [Omega_c*h^2,Omega_b*h^2,A_s]? If no, the code will use [Omega_m,Omega_b,sigma_8]
use_cmb_params= no
#Do you wanna keep the Cl files?
save_cl_files= yes
#Do you wanna keep the CLASS param files?
save_param_files= yes
#Do you want to just generate the power spectra and not the Fisher matrix?
just_run_cls= no
#File containing a "measured" power spectrum to estimate the associated parameter bias
bias_file= none
#Do you want to write the Fisher matrix in terms of power spectra?
do_pspec_fisher= yes
                                                                                                                                param_mods_srd_sv.ini                                                                               0000644 0001750 0001750 00000010622 13276363447 015107  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                #Cosmological parameters
[och2]
#Fiducial value
x= 0.1197
#Increment used for numerical derivatives
dx= 0.001
#Set to 'no' if this parameter should be keep fixed
is_free= yes
onesided=0

[obh2]
x= 0.02222
dx= 0.0001
is_free= yes
onesided=0

#Note om, ob and s8 are only used if use_cmb_params= no
[om]
x= 0.3156
dx= 0.002
is_free= yes
onesided=0

[ob]
x= 0.0492
dx= 0.0002
is_free= yes
onesided=0

[hh]
x= 0.69
dx= 0.01
is_free= no
onesided=0

[w0]
x= -1.0
dx= 0.01
is_free= no
onesided=0

[wa]
x=0.0
dx=0.05
is_free = no
onesided=0

[s8]
x= 0.831
dx= 0.01
is_free= yes
onesided=0

[A_s]
x= 2.1955
dx= 0.01
is_free= yes
onesided=0

[ns]
x= 0.9655
dx= 0.005
is_free= no
onesided=0

[tau]
x=0.06
dx=0.01
is_free = no
onesided = 0

[rt]
x=0.00
dx=0.005
is_free = no
onesided = 1

[mnu]
x=60.
dx=10.
is_free = no
onesided = 0

[pan]
x=0.
dx=0.02
is_free = no
onesided = 1

[Tracer 1]
tracer_name= CMB_exp
tracer_type= cmb_primary
#Has temperature?
has_t= yes
#Has polarization?
has_p= yes
#Noise level in temperature in units of uK-arcmin in the two ell regimes
sigma_t= 1.0 1.0
#Noise level in polarization in units of uK-arcmin in the two ell regimes
sigma_p= 1.4 1.4
#Beam size in units of arcmin in the two ell regimes
beam_amin= 3.0 3.0 
#Ell value that marks the transition
l_transition= 2 
#Minimum ell
lmin= 30 
#Maximum ell
lmax= 5000 
use_tracer= yes

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

[Tracer 3]
tracer_name= Galaxy_survey
tracer_type= gal_clustering
#File describing the redshift bins for this tracer
bins_file= curves/bins_c.txt
#Overall true-redshift distribution
nz_file= curves/nz_c.txt
#Clustering bias
bias_file= curves/bz_c.txt
#Magnification bias
sbias_file= curves/sz_c.txt
#Evolution bias
ebias_file= curves/ez_c.txt
use_tracer= yes

[Tracer 4]
tracer_name= Shear_survey
tracer_type= gal_shear
bins_file= curves/bins_s.txt
nz_file= curves/nz_s.txt
#Alignment bias (only needed if including IA)
abias_file= curves/az_s.txt
#Red (aligned) fraction (only needed if including IA)
rfrac_file= curves/rf_s.txt
#Intrinsic ellipticity dispertion
sigma_gamma= 0.28
use_tracer= yes

[Tracer 5]
tracer_name= Shear_survey
tracer_type= gal_clustering
bins_file= curves/bins_s.txt
nz_file= curves/nz_s.txt
bias_file= curves/bz_s.txt
sbias_file= curves/sz_s.txt
ebias_file= curves/ez_s.txt
use_tracer= yes

[CLASS parameters]
#File with bandpowers. Two columns corresponding to the l_min and l_max(inclusive) of each bandpower
bandpowers_file= curves/larr.txt
#Maximum multipole for which the CMB Cls will be computed
#(only relevant if CMB primary or CMB lensing are included)
lmax_cmb= 5000
#Maximum multipole for which the clustering and lensing Cls will be computed
#(only relevant if galaxy clustering or shear are included)
lmax_lss= 2000
#Minimum multipole from which Limber will be used
lmin_limber= 0
#Include intrinsic alignments in the shear power spectra?
include_alignment= no
#Include RSDs in the clustering power spectra?
include_rsd= no
#Include lensing magnification in the clustering power spectra?
include_magnification= no
#Include relativistic effects in the clustering power spectra?
include_gr_vel= no
include_gr_pot= no
#Command used to execute classt (it should take the CLASS param file as an argument)
exec_path= ./class_mod
#Use non-linear matter transfer function (HALOFit)?
use_nonlinear= yes
#Set to "yes" if you want to include baryonic effects in the power spectrum
use_baryons= yes
#Sky fraction (20,000 sq deg for LSST)
f_sky= 0.4

[Output parameters]
#Directory where all the data will be output
output_dir= outputs
#Prefix of all power-spectrum-related output files
output_spectra= run
#Directory where the Fisher information and plots will be stored
output_fisher= Fisher

[Behaviour parameters]
#Cosmological  model
model= LCDM
#Use [Omega_c*h^2,Omega_b*h^2,A_s]? If no, the code will use [Omega_m,Omega_b,sigma_8]
use_cmb_params= no
#Do you wanna keep the Cl files?
save_cl_files= yes
#Do you wanna keep the CLASS param files?
save_param_files= yes
#Do you want to just generate the power spectra and not the Fisher matrix?
just_run_cls= no
#File containing a "measured" power spectrum to estimate the associated parameter bias
bias_file= none
#Do you want to write the Fisher matrix in terms of power spectra?
do_pspec_fisher= yes
                                                                                                              unwrap_lft.py                                                                                       0000644 0001750 0001750 00000002031 13274606641 013424  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                import numpy as np

ncols = int(3.)

def unwrap_matrix(matrix):
    nc= len(matrix)
    vector = np.zeros(nc*(ncols+1)/2)
    ii = 0
    jj = 0
    kk = 0
    start = 0

    while True:
        if kk == len(vector):
            break
        if ii == len(matrix) - start:
            start += 1
            jj = start
            ii = 0
        vector[kk] = matrix[ii, jj]
        ii += 1
        jj += 1
        kk += 1

    return vector

def wrap_vector(vector):
    matrix = np.zeros( [ncols, ncols] )

    ii = 0
    jj = 0
    kk = 0
    start = 0

    while True:
        if kk == len(vector):
            break
        if ii == len(matrix) - start:
            start += 1
            jj = start
            ii = 0
        matrix[ii, jj] = vector[kk]
        matrix[jj, ii] = vector[kk]
        ii += 1
        jj += 1
        kk += 1

    return matrix



a = np.zeros([ncols,ncols])
for ii in range(0,ncols):
    for jj in range(0,ncols):
        a[ii,jj] = np.random.rand()

b = unwrap_matrix(a)
print a
print b
c = wrap_vector(b)
print c
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                       unwrap_lft.pyc                                                                                      0000644 0001750 0001750 00000001110 13274605462 013564  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                �
)�Zc           @   s   d  d l  Z d �  Z d S(   i����Nc         C   s�   t  |  � } t j | | d d � } d } d } d } d } x} t r� | t  | � k r^ Pn  | t  |  � k r� | d 7} | } d } n  |  | | f | | <| d 7} | d 7} | d 7} qB W| S(   Ni   i   i    (   t   lent   npt   zerost   True(   t   matrixt   ncolst   vectort   iit   jjt   kkt   start(    (    s   unwrap_lft.pyt   unwrap_matrix   s$    	
	

(   t   numpyR   R   (    (    (    s   unwrap_lft.pyt   <module>   s                                                                                                                                                                                                                                                                                                                                                                                                                                                           unwrap_matrix.py                                                                                    0000644 0001750 0001750 00000001637 13263435333 014152  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                import numpy as np

class GaussSt(object) :
    def __init__(self,nmaps) :
        self.nmaps=nmaps
        self.nvec=(self.nmaps*(self.nmaps+1))/2
        self.ind_2d=np.zeros([self.nvec,2],dtype=int)
        id1d=0
        for i in np.arange(self.nmaps) :
            for j in np.arange(self.nmaps-i) :
                self.ind_2d[id1d,:]=[j,i+j]
                id1d+=1

    def gaussian_covariance(self,mat,weight) :
        covar=np.zeros([self.nvec,self.nvec])
        for iv1 in np.arange(self.nvec) :
            i_a=self.ind_2d[iv1,0]; i_b=self.ind_2d[iv1,1];
            for iv2 in np.arange(self.nvec) :
                i_c=self.ind_2d[iv2,0]; i_d=self.ind_2d[iv2,1];
                covar[iv1,iv2]=(mat[i_a,i_c]*mat[i_b,i_d]+mat[i_a,i_d]*mat[i_b,i_c])
        covar*=weight
        return covar
                
    def unwrap_matrix(self,mat) :
        return np.array([mat[ind[0],ind[1]] for ind in self.ind_2d])
                                                                                                 unwrap_matrix.pyc                                                                                   0000644 0001750 0001750 00000002655 13263435345 014321  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                �
�:�Zc           @   s&   d  d l  Z d e f d �  �  YZ d S(   i����Nt   GaussStc           B   s#   e  Z d  �  Z d �  Z d �  Z RS(   c         C   s�   | |  _  |  j  |  j  d d |  _ t j |  j d g d t �|  _ d } xh t j |  j  � D]T } xK t j |  j  | � D]3 } | | | g |  j | d  d  � f <| d 7} q{ Wq^ Wd  S(   Ni   i   t   dtypei    (   t   nmapst   nvect   npt   zerost   intt   ind_2dt   arange(   t   selfR   t   id1dt   it   j(    (    s   unwrap_matrix.pyt   __init__   s    	!#c   
      C   s�   t  j |  j |  j g � } x� t  j |  j � D]� } |  j | d f } |  j | d f } x� t  j |  j � D]p } |  j | d f } |  j | d f }	 | | | f | | |	 f | | |	 f | | | f | | | f <qm Wq. W| | 9} | S(   Ni    i   (   R   R   R   R   R   (
   R	   t   matt   weightt   covart   iv1t   i_at   i_bt   iv2t   i_ct   i_d(    (    s   unwrap_matrix.pyt   gaussian_covariance   s      L
c         C   s5   t  j g  |  j D] } | | d | d f ^ q � S(   Ni    i   (   R   t   arrayR   (   R	   R   t   ind(    (    s   unwrap_matrix.pyt   unwrap_matrix   s    (   t   __name__t
   __module__R   R   R   (    (    (    s   unwrap_matrix.pyR       s   	
	
(   t   numpyR   t   objectR    (    (    (    s   unwrap_matrix.pyt   <module>   s                                                                                      weighting_test.py                                                                                   0000644 0001750 0001750 00000006730 13274613416 014277  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                import numpy as np
import sys

# Number of bins for each tracer
nbins_lft = int(sys.argv[1])

def unwrap_matrix(matrix):
    nc= len(matrix)
    vector = np.zeros(nc*(nc+1)/2)
    ii = 0
    jj = 0
    kk = 0
    start = 0

    while True:
        if kk == len(vector):
            break
        if ii == len(matrix) - start:
            start += 1
            jj = start
            ii = 0
        vector[kk] = matrix[ii, jj]
        ii += 1
        jj += 1
        kk += 1

    return vector

def wrap_vector(vector):
    matrix = np.zeros( [nbins_lft*2, nbins_lft*2] )

    ii = 0
    jj = 0
    kk = 0
    start = 0

    while True:
        if kk == len(vector):
            break
        if ii == len(matrix) - start:
            start += 1
            jj = start
            ii = 0
        matrix[ii, jj] = vector[kk]
        matrix[jj, ii] = vector[kk]
        ii += 1
        jj += 1
        kk += 1

    return matrix



data_matrix = np.zeros( [ nbins_lft*2, nbins_lft*2 ] )
for ii in range(0,len(data_matrix)):
    for jj in range(0, len(data_matrix)):
        data_matrix[ii, jj] = round(10*np.random.rand(), 1)
        data_matrix[jj, ii] = round(10*np.random.rand(), 1)
data_vector = unwrap_matrix(data_matrix)
cov_matrix = np.zeros([len(data_vector), len(data_vector)])
for ii in range(0,len(cov_matrix)):
    for jj in range(0, len(cov_matrix)):
        cov_matrix[ii, jj] = round(10*np.random.rand(), 1)
        cov_matrix[jj, ii] = round(10*np.random.rand(), 1)

# define a weighting vector, currently contains only 0 and 1
weight_lft = np.ones(nbins_lft*(2*nbins_lft+1))
                                                                                                                         
# find the startpoints of the individiual diagonals
startpoints_lft = np.zeros(nbins_lft*2 + 1, dtype=int)
for ii in np.arange(1, nbins_lft*2):
    startpoints_lft[ii] = startpoints_lft[ii-1] + (nbins_lft*2 - ii + 1)
startpoints_lft[-1] = startpoints_lft[-2]+1

# (1) Cut the dndn part 
dndn_scheme_lft = str(sys.argv[2])
# A ... take only leading diagonal
# B ... take all
if dndn_scheme_lft == 'A':
    for ii in np.arange(1, nbins_lft):
        weight_lft[ startpoints_lft[ii] : startpoints_lft[ii] + nbins_lft - ii ] = 0.
                                                                                                                         
# (2) Cut the e-dn part
edn_scheme_lft = str(sys.argv[3])
nrows_edn_lft = int(sys.argv[4])
# a ... take all
# b ... remove lower triangle
# c ... remove b + main diagonal
# d ... remove c + nrows_edn_lft
if edn_scheme_lft == 'b' or edn_scheme_lft == 'c' or edn_scheme_lft == 'd':
    for ii in np.arange(1, nbins_lft):
        weight_lft[ startpoints_lft[ii] + nbins_lft - ii : startpoints_lft[ii] + nbins_lft   ] = 0.
if edn_scheme_lft == 'c' or edn_scheme_lft == 'd':
    weight_lft[ startpoints_lft[nbins_lft] : startpoints_lft[nbins_lft + 1] ] = 0.
if edn_scheme_lft == 'd':
    for ii in np.arange(0, nrows_edn_lft):
        weight_lft[ startpoints_lft[nbins_lft + ii + 2] - (nrows_edn_lft - ii) : startpoints_lft[nbins_lft + ii + 2]  ] = 0.




for ii in np.arange(0, len(weight_lft)):
    if weight_lft[ii] < 0.5:
        cov_matrix[ii,:] = 0.
        cov_matrix[:,ii] = 0.
        cov_matrix[ii, ii] = 1.

for ii in range(0, len(data_vector)):
    data_vector[ii] *= weight_lft[ii]
cut_matrix = wrap_vector(data_vector)

#print 'DATA MATRIX'
#print data_matrix
print 'CUT MATRIX'
print cut_matrix
#print 'DATA VECTOR'
#print data_vector
#print 'COV MATRIX'
#print cov_matrix

                                        weighting_test.pyc                                                                                  0000644 0001750 0001750 00000002777 13274605374 014455  0                                                                                                    ustar   leander                         leander                                                                                                                                                                                                                �
�
�Zc           @   s|  d  d l  Z d Z e j e d e d g � Z d �  Z d Z e j e d e d � Z e j e d d e	 �Z
 x@ e j d e d � D]( Z e
 e d e d e d e
 e <q� Wd Z e d k rx: e j d e � D]# Z d	 e e
 e e
 e e e +q� Wn  d
 Z d Z e d k s4e d k s4e d k ryxB e j d e d � D]' Z d	 e e
 e e e e
 e e +qKWn  e d k s�e d k r�d	 e e
 e e
 e d +n  e d k rxJ e j d e � D]3 Z d	 e e
 e e d e e e
 e e d +q�Wn  xo e j d e e � � D]U Z e e d k  rd	 e e d d � f <d	 e d d � e f <d e e e f <qqWd S(   i����Ni   i   c         C   s�   t  |  � } t j | | d d � } d } d } d } d } x} t r� | t  | � k r^ Pn  | t  |  � k r� | d 7} | } d } n  |  | | f | | <| d 7} | d 7} | d 7} qB W| S(   Ni   i   i    (   t   lent   npt   zerost   True(   t   matrixt   ncolst   vectort   iit   jjt   kkt   start(    (    s   weighting_test.pyt   unwrap_matrix   s$    	
	

i    i   t   dtypet   Ag        t   at   bt   ct   dg      �?g      �?(   t   numpyR   t	   nbins_lftR   t   data_matrixR   t	   covar_lftt   onest
   weight_lftt   intt   startpoints_lftt   arangeR   t   dndn_scheme_lftt   edn_scheme_lftt   nrows_edn_lftR    t	   covar_flt(    (    (    s   weighting_test.pyt   <module>   s6   	&$$(4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 