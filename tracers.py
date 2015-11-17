import numpy as np
import matplotlib.pyplot as plt
import os as os
import sys as sys
from scipy.special import erf
from scipy.interpolate import interp1d
from scipy.integrate import quad
from scipy.interpolate import UnivariateSpline as spline
import time

NZ_IN_AMIN=True
S_GAMMA=0.3

def pdf_photo(z,z0,zf,sz) :
    denom=1./np.sqrt(2*sz*sz)
    return 0.5*(erf((zf-z)*denom)-erf((z0-z)*denom))

class NuisanceFunction :
    """ NuisanceFunction """
    name="none" #Function name
    file_name="none" #File name
    prefix="none" #Prefix path
    zbig=[] #Large z array where support is needed
    z_arr=[] #x values of the nodes
    f_arr=[] #y values of the nodes
    i_marg=[] #Whether we marginalize over this node
    df_arr=[] #Numerical derivative interval for this node
    prior_arr=[] #Prior for this node
    der_rel=0.05 #Relative interval for numerical derivatives
    der_abs=0.05 #Absolute interval for numerical derivatives

    def __init__(self,name="none",fname="none",fname_nz="none",prefix="none",typ="bias",
                 der_rel=0.05,der_abs=0.05) :
        self.prefix=prefix
        self.name=name
        self.file_name=fname
        if name!="none" :
            self.der_rel=der_rel
            self.der_abs=der_abs
            if typ=="bias" :
                self.zbig,dum=np.loadtxt(fname_nz,unpack=True)
                self.z_arr,self.f_arr,self.i_marg,self.prior_arr=np.loadtxt(self.file_name,unpack=True)
                self.df_arr=np.zeros(len(self.z_arr))
                for i in np.arange(len(self.z_arr)) :
                    f=self.f_arr[i]
                    if f>0.1: 
                        self.df_arr[i]=self.der_rel*f
                    else :
                        self.df_arr[i]=self.der_abs

                self.write_bias_file(-1,0)
                for i in np.arange(len(self.z_arr)) :
                    if self.i_marg[i]!=0 :
                        self.write_bias_file(i,1)
                        self.write_bias_file(i,-1)
            elif typ=="bins" :
                data=np.loadtxt(self.file_name,unpack=True)
                z0_arr=np.atleast_1d(data[0])
                zf_arr=np.atleast_1d(data[1])
                s_ph_arr=np.atleast_1d(data[2])
                i_marg=np.atleast_1d(data[3])
                prior_rel=np.atleast_1d(data[4])
                self.z_arr=0.5*(z0_arr+zf_arr)
                self.f_arr=s_ph_arr
                self.i_marg=i_marg
                self.prior_arr=prior_rel
                indx=np.where(prior_rel>0)[0]
                self.prior_arr[indx]=prior_rel[indx]*s_ph_arr[indx]
                self.df_arr=0.1*s_ph_arr
                self.write_bins_file(z0_arr,zf_arr,-1,0)
                for i in np.arange(len(self.z_arr)) :
                    if self.i_marg[i]!=0 :
                        self.write_bins_file(z0_arr,zf_arr,i,1)
                        self.write_bins_file(z0_arr,zf_arr,i,-1)

    def get_filename(self,i_node,sig) :
        if sig>0 :
            asig="p"
        else :
            asig="m"
            
        if i_node>=0 :
            return self.prefix+self.name+"d"+asig+"%d"%i_node
        else :
            return self.prefix+self.name+"fid"

    def write_bins_file(self,z0,zf,i_node,sig) :
        if self.name=="none" :
            return self.name

        mask=np.zeros(len(self.z_arr))
        if i_node>=0 :
            if self.i_marg[i_node]==0 :
                sys.exit("This shouldn't have happened")
            mask[i_node]=1
        
        zm=0.5*(zf+z0)
        dz=0.5*(zf-z0)
        sz=self.f_arr+mask*sig*self.df_arr
        fname=self.get_filename(i_node,sig)+"_bins.txt"
        np.savetxt(fname,np.transpose([zm,dz,sz]))
        return fname

    def write_bias_file(self,i_node,sig) :
        if self.name=="none" :
            return self.name

        mask=np.zeros(len(self.z_arr));
        if i_node>=0 :
            if self.i_marg[i_node]==0 :
                sys.exit("This shouldn't have happened")
            mask[i_node]=1
        
        xa=self.z_arr
        ya=self.f_arr+mask*sig*self.df_arr
        yi=interp1d(xa,ya)
        def yf(x) :
            if x<=xa[0] :
                return ya[0]
            elif x>=xa[-1] :
                return ya[-1]
            else :
                return yi(x)
        xb=self.zbig
        yb=np.array([yf(x) for x in xb])
        
        fname=self.get_filename(i_node,sig)
        np.savetxt(fname,np.transpose([xb,yb]))

        return fname

class Tracer :
    """ Tracer """
    name="default"
    tracer_type="gal_clustering"
    nbins=0
    bins_file="bins.txt"
    nz_file="dndz.txt"
    lmin=0
    lmax=500
    lmax_bins=[]

    #Photo-z
    nuisance_phoz =NuisanceFunction()

    #Parameters for galaxy clustering
    nuisance_bias =NuisanceFunction()
    nuisance_sbias=NuisanceFunction()
    nuisance_ebias=NuisanceFunction()

    #Parameters for galaxy shear
    sigma_gamma=S_GAMMA
    nuisance_abias=NuisanceFunction()
    nuisance_rfrac=NuisanceFunction()

    #Parameters for CMB
    sigma_t=0.
    beam_amin=0.
    has_t=True
    has_p=False
    cl_tt_u=None;
    cl_ee_u=None;
    cl_te_u=None;
    cl_bb_u=None;
    cl_tt_l=None;
    cl_ee_l=None;
    cl_te_l=None;
    cl_bb_l=None;

    #Consider this tracer?
    consider_tracer=True

    def __init__(self,par,name,type_str,bins_file,nz_file,
                 bias_file,sbias_file,ebias_file,
                 abias_file,rfrac_file,sigma_gamma,
                 has_t,has_p,sigma_t,beam_amin,
                 number,consider,lmin,lmax) :
        self.lmin=lmin
        self.lmax=lmax
        self.name=name
        self.consider_tracer=consider
        if type_str=="gal_clustering" :
            self.tracer_type="gal_clustering"
            self.bins_file=bins_file
            self.nz_file=nz_file
            self.nuisance_bias =NuisanceFunction("bias_t%d"%number,bias_file,nz_file,par.output_dir+"/","bias")
            if par.include_magnification or par.include_gr_vel or par.include_gr_pot :
                self.nuisance_sbias=NuisanceFunction("sbias_t%d"%number,sbias_file,nz_file,par.output_dir+"/","bias")
            if par.include_gr_vel or par.include_gr_pot :
                self.nuisance_ebias=NuisanceFunction("ebias_t%d"%number,ebias_file,nz_file,par.output_dir+"/","bias")
            #Get number of bins
            data=np.loadtxt(self.bins_file,unpack=True)
            z0_arr=np.atleast_1d(data[0])
            self.nbins=len(z0_arr)
            self.nuisance_phoz=NuisanceFunction("phoz_t%d"%number,self.bins_file,self.nz_file,par.output_dir+"/","bins")
        elif type_str=="gal_shear" :
            self.tracer_type="gal_shear"
            self.bins_file=bins_file
            self.nz_file=nz_file
            if par.include_alignment :
                self.nuisance_abias=NuisanceFunction("abias_t%d"%number,abias_file,nz_file,par.output_dir+"/","bias")
                self.nuisance_rfrac=NuisanceFunction("rfrac_t%d"%number,rfrac_file,nz_file,par.output_dir+"/","bias")
            self.sigma_gamma=sigma_gamma
            #Get number of bins
            data=np.loadtxt(self.bins_file,unpack=True)
            z0_arr=np.atleast_1d(data[0])
            self.nbins=len(z0_arr)
            self.nuisance_phoz=NuisanceFunction("phoz_t%d"%number,self.bins_file,self.nz_file,par.output_dir+"/","bins")
        elif type_str=="cmb_lensing" :
            self.tracer_type="cmb_lensing"
            self.beam_amin=beam_amin
            self.sigma_t=sigma_t
            self.nbins=1
        elif type_str=="cmb_primary" :
            self.tracer_type="cmb_primary"
            self.sigma_t=sigma_t
            self.beam_amin=beam_amin
            self.has_t=has_t
            self.has_p=has_p
            self.nbins=0
            if self.has_t :
                self.nbins+=1
            if self.has_p :
                self.nbins+=2
            if self.nbins==0 :
                sys.exit("CMB primary should have at least T or P")
        else :
            strout="Wrong tracer type "+type_str+"\n"
            strout+="Allowed tracer types are: gal_clustering, gal_shear, cmb_lensing, cmb_primary\n"
            sys.exit(strout)

    def get_nbins(self) :
        #Bins
        if (self.tracer_type=="gal_clustering") or (self.tracer_type=="gal_shear") :
            data=np.loadtxt(self.bins_file,unpack=True)
            if len(data)>5 :
                self.lmax_bins=np.atleast_1d(data[5])
            else :
                self.lmax_bins=self.lmax*np.ones(self.nbins)

        elif (self.tracer_type=="cmb_lensing") or (self.tracer_type=="cmb_primary") :
            self.lmax_bins=self.lmax*np.ones(self.nbins)

        return self.nbins

def get_cross_noise(tr1,tr2,lmax) :
    nbins1=tr1.nbins
    nbins2=tr2.nbins

    cl_noise=np.zeros([lmax+1,nbins1,nbins2])
    if tr1.name==tr2.name :
        if tr1.tracer_type=='gal_clustering' :
            data1=np.loadtxt(tr1.bins_file,unpack=True)
            z0_arr1=np.atleast_1d(data1[0])
            zf_arr1=np.atleast_1d(data1[1])
            sz_arr1=np.atleast_1d(data1[2])
            z_nz_arr,nz_nz_arr=np.loadtxt(tr1.nz_file,unpack=True)
            if NZ_IN_AMIN==True :
                nz_nz_arr*=(180.*60/np.pi)**2
            dz_arr=z_nz_arr[1:]-z_nz_arr[:-1]
            for i in np.arange(nbins1) :
                parr=pdf_photo(z_nz_arr,z0_arr1[i],zf_arr1[i],sz_arr1[i])
#                np.savetxt(tr1.name+"_b%d"%(i+1)+"_nz.txt",np.transpose([z_nz_arr,parr*nz_nz_arr]))
                integ=interp1d(z_nz_arr,parr*nz_nz_arr)
                ndens=quad(integ,z_nz_arr[0],z_nz_arr[-1])[0]
                cl_noise[:,i,i]=1./ndens
        elif tr1.tracer_type=='gal_shear' :
            data1=np.loadtxt(tr1.bins_file,unpack=True)
            z0_arr1=np.atleast_1d(data1[0])
            zf_arr1=np.atleast_1d(data1[1])
            sz_arr1=np.atleast_1d(data1[2])
            z_nz_arr,nz_nz_arr=np.loadtxt(tr1.nz_file,unpack=True)
            if NZ_IN_AMIN==True :
                nz_nz_arr*=(180.*60/np.pi)**2
            dz_arr=z_nz_arr[1:]-z_nz_arr[:-1]
            for i in np.arange(nbins1) :
                parr=pdf_photo(z_nz_arr,z0_arr1[i],zf_arr1[i],sz_arr1[i])
                integ=interp1d(z_nz_arr,parr*nz_nz_arr)
                ndens=quad(integ,z_nz_arr[0],z_nz_arr[-1])[0]
                cl_noise[:,i,i]=tr1.sigma_gamma**2/ndens
        elif tr1.tracer_type=='cmb_lensing' :
            #Compute temperature noise level
            sigma2_rad=(tr1.sigma_t/(2.725*1E6*180*60/np.pi))**2
            beam_rad=tr1.beam_amin*np.pi/(180*60)/(2*np.sqrt(2*np.log(2)))
            
            #Set up dictionaries
            cl_fid={}; nl_fid={}; fields=[]
            q={}; q['lMin']={}; q['lMax']={}
            
            if tr1.cl_tt_l!=None :
                lm=len(tr1.cl_tt_l)
                l=np.arange(lm-2)+2
                cl_fid['TT']=tr1.cl_tt_l[2:lm]
                nl_fid['TT']=sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                cl_fid['TT_unlen']=tr1.cl_tt_u[2:lm]
                fields.append('TT')
                q['lMin']['TT']=50
                q['lMax']['TT']=3000
            
            if tr1.cl_ee_l!=None :
                lm=len(tr1.cl_ee_l)
                l=np.arange(lm-2)+2
                cl_fid['EE']=tr1.cl_ee_l[2:lm]
                nl_fid['EE']=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                cl_fid['EE_unlen']=tr1.cl_ee_u[2:lm]
                fields.append('EE')
                q['lMin']['EE']=50
                q['lMax']['EE']=4000
            
            if tr1.cl_te_l!=None :
                lm=len(tr1.cl_te_l)
                l=np.arange(lm-2)+2
                cl_fid['TE']=tr1.cl_te_l[2:lm]
                nl_fid['TE']=0.*l
                cl_fid['TE_unlen']=tr1.cl_te_u[2:lm]
                fields.append('TE')
                q['lMin']['TE']=50
                q['lMax']['TE']=3000
            
            if tr1.cl_bb_l!=None :
                lm=len(tr1.cl_bb_l)
                l=np.arange(lm-2)+2
                cl_fid['BB']=tr1.cl_bb_l[2:lm]
                nl_fid['BB']=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                cl_fid['BB_unlen']=tr1.cl_bb_u[2:lm]
                fields.append('BB')
                q['lMin']['BB']=50
                q['lMax']['BB']=4000
            
            larr=np.arange(len(cl_fid['TT_unlen']))+2

            ell,nl=get_lensing_noise(larr,cl_fid,nl_fid,fields,q)
            cl_noise[2:,0,0]=nl[:lmax-1]
        elif tr1.tracer_type=='cmb_primary' :
            sigma2_rad=(tr1.sigma_t/(2.725*1E6*180*60/np.pi))**2
            beam_rad=tr1.beam_amin*np.pi/(180*60)/(2*np.sqrt(2*np.log(2)))
            l=np.arange(lmax+1)
            if tr1.has_t :
                cl_noise[:,0,0]=sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                if tr1.has_p :
                    cl_noise[:,1,1]=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                    cl_noise[:,2,2]=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2) 
            elif tr1.has_p :
                cl_noise[:,0,0]=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
                cl_noise[:,1,1]=2*sigma2_rad*np.exp(l*(l+1)*beam_rad**2)
        else :
            strout="Wrong tracer type "+tr1.tracer_type+"\n"
            strout+="Allowed tracer types are: gal_clustering, gal_shear, cmb_lensing, cmb_primary\n"
            sys.exit(strout)

    return cl_noise

def get_lensing_noise(ell, cl_fid, nl, fields, q):
#    print fields
    lmin = int(ell[0])
    lmax = int(ell[-1])
    
#    print 'lmin =', lmin, 'lmax =', lmax
    cl_tot = {}
    for key in ['TT','EE','BB','TE','Td','dd']:
        if key != 'dd' and key !='Td': 
            cl_tot[key] = cl_fid[key] + nl[key]

    cl_unlen = {}
    for key in ['TT','EE','BB','TE','Td','dd']:
        if key != 'dd' and key !='Td': 
            cl_unlen[key] = cl_fid[key+'_unlen']
	
## n_Ls specifies approx number of L's at which to 
## evaluate the lensing noise spectrum.
## Need unique L's for the splining used later. 
    n_Ls = 200
    LogLs = np.linspace(np.log(2),np.log(lmax+0.1), n_Ls)
    Ls = np.unique(np.floor(np.exp(LogLs)).astype(int))
		
    deltaL = 50 ## This is the Fourier space resoltuion
    deltaLx, deltaLy = deltaL, deltaL ## They don't have to be same (depends on patch geometry)
    lx = np.arange(-lmax, lmax, deltaLx)
    ly = np.arange(-lmax, lmax, deltaLy)
    lx, ly = np.meshgrid(lx, ly)
    l = np.sqrt(lx**2 + ly**2)

## Only use T,P modes with L_lens_min < l < L_lens_max 
## for the lensing reconstruction.
    l_lens_min_temp = q['lMin']['TT']
    l_lens_min_pol	= q['lMin']['EE']
    l_lens_max_temp = q['lMax']['TT']
    l_lens_max_pol  = q['lMax']['EE']
		
    clunlenFunc, clunlenArr = {}, {}
    cltotFunc, cltotArr = {}, {}
	
    for key in ['TT','EE','BB','TE']:
        clunlenFunc[key] = interp1d(ell, cl_unlen[key], kind = 'linear', bounds_error = False, fill_value = 0. )
        clunlenArr[key] = clunlenFunc[key](l)
        if key == 'TT':
            cl_tot[key][np.where(ell < l_lens_min_temp)] *= 1e6
            cl_tot[key][np.where(ell > l_lens_max_temp)] *= 1e6
        if key == 'EE' or key == 'BB':
            cl_tot[key][np.where(ell < l_lens_min_pol)] *= 1e6
            cl_tot[key][np.where(ell > l_lens_max_pol)] *= 1e6	
        if key != 'TE':
            cltotFunc[key] = interp1d(ell, cl_tot[key], kind = 'linear', bounds_error = False, fill_value = 1e50 )
        else:
            cltotFunc[key] = interp1d(ell, cl_tot[key], kind = 'linear', bounds_error = False, fill_value = 0. )	
        cltotArr[key] = cltotFunc[key](l)
	
    a = time.time()
    lensingSpec = []
    if 'TT' in fields : 
        lensingSpec += ['TT']
    if 'TE' in fields:
        lensingSpec += ['TE','TB']
    if 'EE' in fields:
        lensingSpec += ['EE', 'EB']
	
    NL = {}
    for field in lensingSpec:
        if field == 'TT':
#            print 'Field:', field
            cl1_unlen = clunlenArr['TT']
            cl1_tot = cltotArr['TT']
            NL['TT'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_unlen = clunlenFunc['TT'](l2)
                cl2_tot = cltotFunc['TT'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                kernel = ( cl1_unlen * Ldotl1 + cl2_unlen * Ldotl2 ) ** 2 / (2 * cl1_tot * cl2_tot)
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['TT'].append((integral)**(-1))
            NL['TT'] = Ls**2 * (Ls+1)**2 * np.array(NL['TT']) / 4.
		
        if field == 'TE':
#            print 'Field:', field
            cl1_unlen = clunlenArr['TE']
            cl1_tot = cltotArr['EE']
            NL['TE'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_unlen = clunlenFunc['TE'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                cos2phi = np.nan_to_num(2 * ((Ldotl1 - l**2) / (l * l2))**2 - 1.)
                
                f_l1l2 = (cl1_unlen * cos2phi * Ldotl1 + cl2_unlen * Ldotl2)
                f_l2l1 = (cl2_unlen * cos2phi * Ldotl2 + cl1_unlen * Ldotl1)
                
                F_l1l2 = (cltotArr['EE'] * cltotFunc['TT'](l2) * f_l1l2 - cltotArr['TE'] * cltotFunc['TE'](l2) * f_l2l1) / \
                    (cltotArr['TT']*cltotArr['EE']*cltotFunc['EE'](l2)*cltotFunc['TT'](l2) - (cltotArr['TE'] * cltotFunc['TE'](l2))**2)
                
                kernel = f_l1l2 * F_l1l2
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['TE'].append((integral)**(-1))
            NL['TE'] = Ls**2 * (Ls+1)**2 * np.array(NL['TE']) / 4.	
		
        if field == 'TB':
#            print 'Field:', field
            cl1_unlen = clunlenArr['TE']
            cl1_tot = cltotArr['TT']
            NL['TB'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_tot = cltotFunc['BB'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                cosphi = np.nan_to_num((Ldotl1 - l**2) / (l * l2))
                sin2phiSq = 4. * cosphi**2 * (1 - cosphi**2)
                kernel = sin2phiSq * ( cl1_unlen * Ldotl1 ) ** 2 / (cl1_tot * cl2_tot)
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['TB'].append((integral)**(-1))
            NL['TB'] = Ls**2 * (Ls+1)**2 * np.array(NL['TB']) / 4.
		
        if field == 'EE':
#            print 'Field:', field
            cl1_unlen = clunlenArr['EE']
            cl1_tot = cltotArr['EE']
            NL['EE'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_unlen = clunlenFunc['EE'](l2)
                cl2_tot = cltotFunc['EE'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                cos2phi = np.nan_to_num(2. * ((Ldotl1 - l**2) / (l * l2))**2 - 1.)
                kernel = ( (cl1_unlen * Ldotl1 + cl2_unlen * Ldotl2) * cos2phi ) ** 2 / (2 * cl1_tot * cl2_tot)
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['EE'].append((integral)**(-1))
            NL['EE'] = Ls**2 * (Ls+1)**2 * np.array(NL['EE']) / 4.
		
        if field == 'BB':
#            print 'Field:', field
            cl1_unlen = clunlenArr['BB']
            cl1_tot = cltotArr['BB']
            NL['BB'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_unlen = clunlenFunc['BB'](l2)
                cl2_tot = cltotFunc['BB'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                cos2phi = 2. * np.nan_to_num(((Ldotl1 - l**2) / (l * l2))**2) - 1.
                kernel = ( (cl1_unlen * Ldotl1 + cl2_unlen * Ldotl2) * cos2phi ) ** 2 / (2 * cl1_tot * cl2_tot)
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['BB'].append((integral)**(-1))
            NL['BB'] = Ls**2 * (Ls+1)**2 * np.array(NL['BB']) / 4.
		
        if field == 'EB':
#            print 'Field:', field
            cl1_unlen = clunlenArr['EE']
            cl1_tot = cltotArr['EE']
            NL['EB'] = []
            for L in Ls:
                l2 = np.sqrt( (L - lx)**2 + ly**2 )
                cl2_tot = cltotFunc['BB'](l2)
                cl2_unlen = clunlenFunc['BB'](l2)
                Ldotl1 = L * lx
                Ldotl2 = L * (L - lx)
                cosphi =np.nan_to_num( (Ldotl1 - l**2) / (l * l2) )
                sin2phiSq = 4. * cosphi**2 * (1 - cosphi**2)
                kernel = sin2phiSq * ( cl1_unlen * Ldotl1 - cl2_unlen * Ldotl2 ) ** 2 / (cl1_tot * cl2_tot)
                integral = np.sum(kernel) * (2 * np.pi)**(-2.) * deltaL**2
                NL['EB'].append((integral)**(-1))
            NL['EB'] = Ls**2 * (Ls+1)**2 * np.array(NL['EB']) / 4.
		
    NL_mv = np.zeros(np.shape(Ls))
## Assumes negligible correlation between noise terms
#print lensingSpec
    for field in lensingSpec:
        NL_mv += 1./NL[field]
    NL_mv = 1./NL_mv
	
    logLs = np.log(ell)
    s = spline(np.log(Ls), np.log(NL_mv), s=0)
    NL_mv = np.exp(s(logLs))
    
    return ell, NL_mv
