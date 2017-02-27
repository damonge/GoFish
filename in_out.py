import numpy as np
import struct as struct
import sys as sys
import common as com
import os as os
import matplotlib.pyplot as plt
import array as array
from scipy.interpolate import interp1d

'''
def my_parser(string) :
    lst=len(string)
    i=0
    notfound=True
    while notfound :
        if string[i]=='_' :
            notfound=False
        i+=1
        if i==lst :
            notfound=False

    if i==lst :
        prefix="none"
        i_t=-1
        i_n=-1
    else :
        prefix=string[:i-1]
        suffix=string[i:]
    
        id1=0
        if suffix[id1]!="t" :
            i_t=-1
            i_n=-1
        else :
            id2=0
            lsf=len(suffix)
            notfound=True
            while notfound :
                if suffix[id2]=="n" :
                    notfound=False
                id2+=1
                if id2==lsf :
                    notfound=False
            if id2==lsf :
                i_t=-1
                i_n=-1
            else :
                i_t=int(suffix[id1+1:id2-1])
                i_n=int(suffix[id2:])

    return prefix,i_t,i_n
'''

def my_parser(string) :
    lst=len(string)
    i=0
    notfound=True
    while notfound :
        if string[i]=='_' :
            notfound=False
        i+=1
        if i==lst :
            notfound=False

    if i==lst :
        prefix="none"
        tr_name="none"
        i_n=-1
    else :
        prefix=string[:i-1]
        suffix=string[i:]
    
        id1=0
        if ((prefix!="bias") and (prefix!="sbias") and (prefix!="ebias") and
            (prefix!="abias") and (prefix!="rfrac") and
            (prefix!="sphz") and (prefix!="bphz")) :
            tr_name="none"
            i_n=-1
        else :
            id2=0
            lsf=len(suffix)
            notfound=True
            while notfound :
                if id2==lsf-4 :
                    notfound=False
                if suffix[id2:id2+4]=="node" :
                    notfound=False
                id2+=1
            if id2==lsf-3 :
                print "shit"
                tr_name="none"
                i_n=-1
            else :
                tr_name=suffix[id1:id2-2]
                i_n=int(suffix[id2+3:])

    return prefix,tr_name,i_n

def read_cls_camb(fname) :
    data=np.loadtxt(fname,unpack=True)
    lmax=len(data[0])+1
    n_ct=len(data)

    has_tt=True; has_ee=True; has_bb=True; has_te=True;
    has_pp=False; has_tp=False; has_ep=False;
    has_td=False; has_tl=False; has_pd=False; has_pl=False;
    has_dd=False; has_dl=False; has_ll=False;

    cl_tt=None; cl_ee=None; cl_te=None; cl_bb=None;
    cl_pp=None; cl_tp=None; cl_ep=None; 
    cl_td=None; cl_tl=None; cl_pd=None; cl_pl=None;
    cl_dd=None; cl_dl=None; cl_ll=None; 

    larr=data[0]
    prefac=2*np.pi/(larr*(larr+1))
    prefac_cmbl=0.5*larr*(larr+1)
    cl_tt=np.zeros(lmax+1); cl_tt[2:]=prefac*data[1]
    cl_ee=np.zeros(lmax+1); cl_ee[2:]=prefac*data[2]
    cl_bb=np.zeros(lmax+1); cl_bb[2:]=prefac*data[3]
    cl_te=np.zeros(lmax+1); cl_te[2:]=prefac*data[4]
    if n_ct>5 :
        has_pp=True; has_tp=True; has_ep=True;
        cl_pp=np.zeros(lmax+1); cl_pp[2:]=prefac*prefac_cmbl**2*data[5]/(larr*(larr+1))
        cl_tp=np.zeros(lmax+1); cl_tp[2:]=prefac*prefac_cmbl*data[6]/(larr+0.5)
        cl_ep=np.zeros(lmax+1); cl_ep[2:]=prefac*prefac_cmbl*data[7]/(larr+0.5)
    dic={'lmax':lmax,'n_nc':0,'n_wl':0}
    dic.update({'has_tt':has_tt,'cl_tt':cl_tt})
    dic.update({'has_ee':has_ee,'cl_ee':cl_ee})
    dic.update({'has_te':has_te,'cl_te':cl_te})
    dic.update({'has_bb':has_bb,'cl_bb':cl_bb})
    dic.update({'has_pp':has_pp,'cl_pp':cl_pp})
    dic.update({'has_tp':has_tp,'cl_tp':cl_tp})
    dic.update({'has_ep':has_ep,'cl_ep':cl_ep})
    dic.update({'has_td':has_td,'cl_td':cl_td})
    dic.update({'has_tl':has_tl,'cl_tl':cl_tl})
    dic.update({'has_pd':has_pd,'cl_pd':cl_pd})
    dic.update({'has_pl':has_pl,'cl_pl':cl_pl})
    dic.update({'has_dd':has_dd,'cl_dd':cl_dd})
    dic.update({'has_dl':has_dl,'cl_dl':cl_dl})
    dic.update({'has_ll':has_ll,'cl_ll':cl_ll})

    return dic

def read_cls_class(fname) :
    f=open(fname,"rd")
    data=f.read()
    f.close()

    lmax,n_nc,n_wl=struct.unpack('3I',data[0:12])
    has_tt,has_ee,has_te,has_bb,has_pp,has_tp,has_ep=struct.unpack('7I',data[12:12+7*4])
    has_td,has_tl,has_pd,has_pl,has_dd,has_dl,has_ll=struct.unpack('7I',data[12+7*4:12+14*4])
    offset_header=17*4
    n_ct=0
    if has_tt :
        n_ct+=1
    if has_ee :
        n_ct+=1
    if has_te :
        n_ct+=1
    if has_bb :
        n_ct+=1
    if has_pp :
        n_ct+=1
    if has_tp :
        n_ct+=1
    if has_ep :
        n_ct+=1
    if has_td :
        n_ct+=n_nc
    if has_tl :
        n_ct+=n_wl
    if has_pd :
        n_ct+=n_nc
    if has_pl :
        n_ct+=n_wl
    if has_dd :
        n_ct+=n_nc*(n_nc+1)/2
    if has_dl :
        n_ct+=n_nc*n_wl
    if has_ll :
        n_ct+=n_wl*(n_wl+1)/2
    len_cls=8*n_ct*(lmax-1)

    if len(data)!=len_cls+offset_header :
        strout="Error reading CLASS file "+fname
        strout+="%.2lf"%((len(data)-offset_header)/8.)+" %.2lf"%(len_cls/8.)
        sys.exit(strout)

    dum=array.array('d')
    dum.fromstring(data[offset_header:])#.
    cldata=np.array(dum).reshape((lmax-1,n_ct))

    cl_tt=None; cl_ee=None; cl_te=None; cl_bb=None;
    cl_pp=None; cl_tp=None; cl_ep=None; 
    cl_td=None; cl_tl=None; cl_pd=None; cl_pl=None;
    cl_dd=None; cl_dl=None; cl_ll=None; 

    offset_here=0
    larr=np.arange(lmax-1)+2
    prefac=2*np.pi/(larr*(larr+1))
    prefac_cmbl=0.5*larr*(larr+1)
    if has_tt :
        cl_tt=np.zeros(lmax+1)
        cl_tt[2:]=prefac*cldata[:,offset_here]
        offset_here+=1
    if has_ee :
        cl_ee=np.zeros(lmax+1)
        cl_ee[2:]=prefac*cldata[:,offset_here]
        offset_here+=1
    if has_te :
        cl_te=np.zeros(lmax+1)
        cl_te[2:]=prefac*cldata[:,offset_here]
        offset_here+=1
    if has_bb :
        cl_bb=np.zeros(lmax+1)
        cl_bb[2:]=prefac*cldata[:,offset_here]
        offset_here+=1
    if has_pp :
        cl_pp=np.zeros(lmax+1)
        cl_pp[2:]=prefac*prefac_cmbl**2*cldata[:,offset_here]
        offset_here+=1
    if has_tp :
        cl_tp=np.zeros(lmax+1)
        cl_tp[2:]=prefac*prefac_cmbl*cldata[:,offset_here]
        offset_here+=1
    if has_ep :
        cl_ep=np.zeros(lmax+1)
        cl_ep[2:]=prefac*prefac_cmbl*cldata[:,offset_here]
        offset_here+=1
    if has_td :
        cl_td=np.zeros([lmax+1,n_nc])
        for i in np.arange(n_nc) :
            cl_td[2:,i]=prefac*cldata[:,offset_here+i]
        offset_here+=n_nc
    if has_tl :
        cl_tl=np.zeros([lmax+1,n_wl]) 
        for i in np.arange(n_wl) :
            cl_tl[2:,i]=prefac*cldata[:,offset_here+i]
        offset_here+=n_wl
    if has_pd :
        cl_pd=np.zeros([lmax+1,n_nc])
        for i in np.arange(n_nc) :
            cl_pd[2:,i]=prefac*prefac_cmbl*cldata[:,offset_here+i]
        offset_here+=n_nc
    if has_pl :
        cl_pl=np.zeros([lmax+1,n_wl])
        for i in np.arange(n_wl) :
            cl_pl[2:,i]=prefac*prefac_cmbl*cldata[:,offset_here+i]
        offset_here+=n_wl
    if has_dd :
        cl_dd=np.zeros([lmax+1,n_nc,n_nc])
        for i in np.arange(n_nc) :
            for j in np.arange(n_nc-i)+i :
                cl_dd[2:,i,j]=prefac*cldata[:,offset_here]
                if i!=j :
                    cl_dd[2:,j,i]=prefac*cldata[:,offset_here]
                offset_here+=1
    if has_dl :
        cl_dl=np.zeros([lmax+1,n_nc,n_wl])
        for i in np.arange(n_nc) :
            for j in np.arange(n_wl) :
                cl_dl[2:,i,j]=prefac*cldata[:,offset_here]
                offset_here+=1
    if has_ll :
        cl_ll=np.zeros([lmax+1,n_wl,n_wl])
        for i in np.arange(n_wl) :
            for j in np.arange(n_wl-i)+i :
                cl_ll[2:,i,j]=prefac*cldata[:,offset_here]
                if i!=j :
                    cl_ll[2:,j,i]=prefac*cldata[:,offset_here]
                offset_here+=1

    #Wrap into a dictionary
    dic={'lmax':lmax,'n_nc':n_nc,'n_wl':n_wl}
    dic.update({'has_tt':has_tt,'cl_tt':cl_tt})
    dic.update({'has_ee':has_ee,'cl_ee':cl_ee})
    dic.update({'has_te':has_te,'cl_te':cl_te})
    dic.update({'has_bb':has_bb,'cl_bb':cl_bb})
    dic.update({'has_pp':has_pp,'cl_pp':cl_pp})
    dic.update({'has_tp':has_tp,'cl_tp':cl_tp})
    dic.update({'has_ep':has_ep,'cl_ep':cl_ep})
    dic.update({'has_td':has_td,'cl_td':cl_td})
    dic.update({'has_tl':has_tl,'cl_tl':cl_tl})
    dic.update({'has_pd':has_pd,'cl_pd':cl_pd})
    dic.update({'has_pl':has_pl,'cl_pl':cl_pl})
    dic.update({'has_dd':has_dd,'cl_dd':cl_dd})
    dic.update({'has_dl':has_dl,'cl_dl':cl_dl})
    dic.update({'has_ll':has_ll,'cl_ll':cl_ll})

    return dic

def read_cls(par,fname) :
    if par.use_CAMB :
        return read_cls_camb(fname)
    else :
        return read_cls_class(fname)

def write_camb_param_file(par,param_vary,sign_vary,prefix_out) :
    """ Generates CAMB param file """
    och2,doch2,osid_och2=par.get_param_properties("och2")
    obh2,dobh2,osid_obh2=par.get_param_properties("obh2")
    hh,dhh,osid_hh=par.get_param_properties("hh")
    if par.model=='JBD' :
        obd,dobd,osid_obd=par.get_param_properties("obd")
    if ((par.model=='BADE')) :
        w0,dw0,osid_w0=par.get_param_properties("w0")
        wdeq,dwdeq,osid_wdeq=par.get_param_properties("wdeq")
        deq,ddeq,osid_deq=par.get_param_properties("deq")
    if ((par.model=='EDE')) :
        w0,dw0,osid_w0=par.get_param_properties("w0")
        wdeq,dwdeq,osid_wdeq=par.get_param_properties("wdeq")
        ede,dede,osid_ede=par.get_param_properties("ede")
    if ((par.model=='wCDM') or (par.model=='Horndeski')) :
        w0,dw0,osid_w0=par.get_param_properties("w0")
        wa,dwa,osid_wa=par.get_param_properties("wa")
    if par.model=='Horndeski' :
        bk ,dbk ,osid_bk =par.get_param_properties("bk")
        bb ,dbb ,osid_bb =par.get_param_properties("bb")
        bm ,dbm ,osid_bm =par.get_param_properties("bm")
        bt ,dbt ,osid_bt =par.get_param_properties("bt")
        ck ,dck ,osid_ck =par.get_param_properties("ck")
        cb ,dcb ,osid_cb =par.get_param_properties("cb")
        cm ,dcm ,osid_cm =par.get_param_properties("cm")
        ct ,dct ,osid_ct =par.get_param_properties("ct")
        zh ,dzh ,osid_zh =par.get_param_properties("zh")
        ezh,dezh,osid_ezh=par.get_param_properties("ezh")
        m2i,dm2i,osid_m2i=par.get_param_properties("m2i")
        kv ,dkv ,osid_rt =par.get_param_properties("kv")
    ns,dns,osid_ns=par.get_param_properties("ns")
    a_s,da_s,osid_a_s=par.get_param_properties("A_s")
    tau,dtau,osid_tau=par.get_param_properties("tau")
    mnu,dmnu,osid_mnu=par.get_param_properties("mnu")
    nnu,dnnu,osid_nnu=par.get_param_properties("nnu")
    pan,dpan,osid_pan=par.get_param_properties("pan")
    fNL,dfNL,osid_fNL=par.get_param_properties("fnl")
    rt,drt,osid_rt=par.get_param_properties("rt")
#    val,dval,osid=par.get_param_properties(param_vary)

    def add_fdiff(val,dval,sig,osd) :
        if osd==0 :
            return val+sig*dval
        else :
            return val+osd*(1.5+0.5*sig)*dval

    if param_vary=="och2" :
        och2=add_fdiff(och2,doch2,sign_vary,osid_och2)
    if param_vary=="obh2" :
        obh2=add_fdiff(obh2,dobh2,sign_vary,osid_obh2)
    if param_vary=="hh" :
        hh=add_fdiff(hh,dhh,sign_vary,osid_hh)
    if par.model=='BADE' :
        if param_vary=="w0" :
            w0=add_fdiff(w0,dw0,sign_vary,osid_w0)
        if param_vary=="wdeq" :
            wdeq=add_fdiff(wdeq,dwdeq,sign_vary,osid_wdeq)
        if param_vary=="deq" :
            deq=add_fdiff(deq,ddeq,sign_vary,osid_deq)
        cs2_lam= 0.3333333
        c2vis= 0.33333333
    else :
        w0=-1.0
        wdeq=-0.99999
        deq=0
        cs2_lam= 1.0
        c2vis= 0.0
    if par.model=='EDE' :
        if param_vary=="w0" :
            w0=add_fdiff(w0,dw0,sign_vary,osid_w0)
        if param_vary=="wdeq" :
            wdeq=add_fdiff(wdeq,dwdeq,sign_vary,osid_wdeq)
        if param_vary=="ede" :
            ede=add_fdiff(ede,dede,sign_vary,osid_ede)
        cs2_lam= 0.3333333
        c2vis= 0.33333333
    else :
        w0=-1.0
        wdeq=-0.99999
        ede=0
        cs2_lam= 1.0
        c2vis= 0.0
    if param_vary=="ns" :
        ns=add_fdiff(ns,dns,sign_vary,osid_ns)
    if param_vary=="A_s" :
        a_s=add_fdiff(a_s,da_s,sign_vary,osid_a_s)
    if param_vary=="tau" :
        tau=add_fdiff(tau,dtau,sign_vary,osid_tau)
    if param_vary=="mnu" :
        mnu=add_fdiff(mnu,dmnu,sign_vary,osid_mnu)
    if param_vary=="nnu" :
        nnu=add_fdiff(nnu,dnnu,sign_vary,osid_nnu)
    if param_vary=="pan" :
        pan=add_fdiff(pan,dpan,sign_vary,osid_pan)
    if param_vary=="fnl" :
        fNL=add_fdiff(fNL,dfNL,sign_vary,osid_fNL)
    if param_vary=="rt" :
        rt=add_fdiff(rt,drt,sign_vary,osid_rt)

    strout="#CAMB param file by GoFish\n"
    strout+="get_scalar_cls = T\n"
    strout+="get_vector_cls = F\n"
    strout+="get_tensor_cls = T\n"
    strout+="get_transfer   = F\n"
    strout+="do_lensing     = T\n"
    if par.use_nonlinear==True :
        strout+="do_nonlinear = 3\n"
    else :
        strout+="do_nonlinear = 0\n"
    strout+="l_max_scalar = %d\n"%(par.lmax_cmb+500)
    strout+="l_max_tensor = %d\n"%(par.lmax_cmb)
    strout+="k_eta_max_tensor = %d\n"%(2*par.lmax_cmb)
    strout+="use_physical = T\n"
    strout+="ombh2          = %lE\n"%obh2
    strout+="omch2          = %lE\n"%och2
    strout+="omnuh2         = %lE\n"%(mnu*0.001*0.0107524)
    strout+="omk            = 0\n"
    strout+="hubble         = %lE\n"%(100.*hh)
    strout+="w              = %lE\n"%(w0)
    strout+="w0             = %lE\n"%(wdeq)
    strout+="cs2_lam        = %lE\n"%(cs2_lam)
    strout+="c2vis          = %lE\n"%(c2vis)
    strout+="deq            = %lE\n"%(deq)
    strout+="ede            = %lE\n"%(ede)
    strout+="temp_cmb           = 2.7255\n"
    strout+="helium_fraction    = 0.24\n"
    if mnu>0 :
        strout+="massless_neutrinos = %lE\n"%(nnu)
        strout+="nu_mass_eigenstates = 1\n"
        strout+="massive_neutrinos  = 1\n"
        strout+="share_delta_neff = T\n"
        strout+="nu_mass_fractions = 1\n"
        strout+="nu_mass_degeneracies = \n"
    else :
        strout+="massless_neutrinos = %lE\n"%(nnu)
        strout+="nu_mass_eigenstates = 1\n"
        strout+="massive_neutrinos  = 0\n"
        strout+="share_delta_neff = T\n"
        strout+="nu_mass_fractions = 1\n"
        strout+="nu_mass_degeneracies = \n"
    strout+="initial_power_num         = 1\n"
    strout+="pivot_scalar              = 0.05\n"
    strout+="pivot_tensor              = 0.05\n"
    strout+="scalar_amp(1)             = %lE\n"%(a_s*1E-9)
    strout+="scalar_spectral_index(1)  = %lE\n"%(ns)
    strout+="scalar_nrun(1)            = 0\n"
    strout+="scalar_nrunrun(1)         = 0\n"
    strout+="tensor_spectral_index(1)  = 0\n"
    strout+="tensor_nrun(1)            = 0\n"
    strout+="tensor_parameterization   = 1\n"
    strout+="initial_ratio(1)          = %lE\n"%(rt)
    strout+="reionization         = T\n"
    strout+="re_use_optical_depth = T\n"
    strout+="re_optical_depth     = %lE\n"%(tau)
    strout+="RECFAST_fudge = 1.14\n"
    strout+="RECFAST_fudge_He = 0.86\n"
    strout+="RECFAST_Heswitch = 6\n"
    strout+="RECFAST_Hswitch  = T\n"
    strout+="COBE_normalize = F\n"
    strout+="CMB_outputscale = 1.0\n" #CHECK THIS
    strout+="transfer_high_precision = F\n"
    strout+="transfer_kmax           = 2\n"
    strout+="transfer_k_per_logint   = 0\n"
    strout+="transfer_num_redshifts  = 1\n"
    strout+="transfer_interp_matterpower = T\n"
    strout+="transfer_redshift(1)    = 0\n"
    strout+="transfer_filename(1)    = transfer_out.dat\n"
    strout+="transfer_matterpower(1) = matterpower.dat\n"
    strout+="transfer_power_var = 7\n"
    strout+="output_root = "+prefix_out+"\n"
    strout+="scalar_output_file = scalCls.dat\n"
    strout+="vector_output_file = vecCls.dat\n"
    strout+="tensor_output_file = tensCls.dat\n"
    strout+="total_output_file  = totCls.dat\n"
    strout+="lensed_output_file = lensedCls.dat\n"
    strout+="lensed_total_output_file  =lensedtotCls.dat\n"
    strout+="lens_potential_output_file = lenspotentialCls.dat\n"
    strout+="FITS_filename      = scalCls.fits\n"
    strout+="do_lensing_bispectrum = F\n"
    strout+="do_primordial_bispectrum = F\n"
    strout+="bispectrum_nfields = 1\n"
    strout+="bispectrum_slice_base_L = 0\n"
    strout+="bispectrum_ndelta=3\n"
    strout+="bispectrum_delta(1)=0\n"
    strout+="bispectrum_delta(2)=2\n"
    strout+="bispectrum_delta(3)=4\n"
    strout+="bispectrum_do_fisher= F\n"
    strout+="bispectrum_fisher_noise=0\n"
    strout+="bispectrum_fisher_noise_pol=0\n"
    strout+="bispectrum_fisher_fwhm_arcmin=7\n"
    strout+="bispectrum_full_output_file=\n"
    strout+="bispectrum_full_output_sparse=F\n"
    strout+="bispectrum_export_alpha_beta=F\n"
    strout+="feedback_level = 1\n"
    strout+="output_file_headers = T\n"
    strout+="derived_parameters = T\n"
    strout+="lensing_method = 1\n"
    strout+="accurate_BB = F\n"
    strout+="massive_nu_approx = 1\n"
    strout+="accurate_polarization   = T\n"
    strout+="accurate_reionization   = T\n"
    strout+="do_tensor_neutrinos     = T\n"
    strout+="do_late_rad_truncation   = T\n"
    strout+="halofit_version= \n"
    strout+="number_of_threads       = 0\n"
    strout+="high_accuracy_default=T\n"
    strout+="accuracy_boost          = 1\n"
    strout+="l_accuracy_boost        = 1\n"
    strout+="l_sample_boost          = 1\n"
    f=open(prefix_out+"_param.ini","w")
    f.write(strout)
    f.close()

def write_class_param_file(par,param_vary,sign_vary,prefix_out) :
    """ Generates CLASS param file """
    och2,doch2,osid_och2=par.get_param_properties("och2")
    obh2,dobh2,osid_obh2=par.get_param_properties("obh2")
    hh,dhh,osid_hh=par.get_param_properties("hh")
    if par.model=='JBD' :
        obd,dobd,osid_obd=par.get_param_properties("obd")
    if ((par.model=='wCDM') or (par.model=='Horndeski')) :
        w0,dw0,osid_w0=par.get_param_properties("w0")
        wa,dwa,osid_wa=par.get_param_properties("wa")
    if par.model=='Horndeski' :
        bk ,dbk ,osid_bk =par.get_param_properties("bk")
        bb ,dbb ,osid_bb =par.get_param_properties("bb")
        bm ,dbm ,osid_bm =par.get_param_properties("bm")
        bt ,dbt ,osid_bt =par.get_param_properties("bt")
        ck ,dck ,osid_ck =par.get_param_properties("ck")
        cb ,dcb ,osid_cb =par.get_param_properties("cb")
        cm ,dcm ,osid_cm =par.get_param_properties("cm")
        ct ,dct ,osid_ct =par.get_param_properties("ct")
        zh ,dzh ,osid_zh =par.get_param_properties("zh")
        ezh,dezh,osid_ezh=par.get_param_properties("ezh")
        m2i,dm2i,osid_m2i=par.get_param_properties("m2i")
        kv ,dkv ,osid_rt =par.get_param_properties("kv")
    ns,dns,osid_ns=par.get_param_properties("ns")
    a_s,da_s,osid_a_s=par.get_param_properties("A_s")
    tau,dtau,osid_tau=par.get_param_properties("tau")
    mnu,dmnu,osid_mnu=par.get_param_properties("mnu")
    pan,dpan,osid_pan=par.get_param_properties("pan")
    fNL,dfNL,osid_fNL=par.get_param_properties("fnl")
    rt,drt,osid_rt=par.get_param_properties("rt")
#    val,dval,osid=par.get_param_properties(param_vary)

    def add_fdiff(val,dval,sig,osd) :
        if osd==0 :
            return val+sig*dval
        else :
            return val+osd*(1.5+0.5*sig)*dval

    if param_vary=="och2" :
        och2=add_fdiff(och2,doch2,sign_vary,osid_och2)
    if param_vary=="obh2" :
        obh2=add_fdiff(obh2,dobh2,sign_vary,osid_obh2)
    if param_vary=="hh" :
        hh=add_fdiff(hh,dhh,sign_vary,osid_hh)
    if par.model=='JBD' :
        if param_vary=="obd" :
            obd=add_fdiff(obd,dobd,sign_vary,osid_obd)
    if ((par.model=='wCDM') or (par.model=='Horndeski')) :
        if param_vary=="w0" :
            w0=add_fdiff(w0,dw0,sign_vary,osid_w0)
        if param_vary=="wa" :
            wa=add_fdiff(wa,dwa,sign_vary,osid_wa)
    else :
        w0=-1.0
        wa= 0.0
    if par.model=='Horndeski' :
        if param_vary=="bk" :
            bk=add_fdiff(bk,dbk,sign_vary,osid_bk)
        if param_vary=="bb" :
            bb=add_fdiff(bb,dbb,sign_vary,osid_bb)
        if param_vary=="bm" :
            bm=add_fdiff(bm,dbm,sign_vary,osid_bm)
        if param_vary=="bt" :
            bt=add_fdiff(bt,dbt,sign_vary,osid_bt)
        if param_vary=="ck" :
            ck=add_fdiff(ck,dck,sign_vary,osid_ck)
        if param_vary=="cb" :
            cb=add_fdiff(cb,dcb,sign_vary,osid_cb)
        if param_vary=="cm" :
            cm=add_fdiff(cm,dcm,sign_vary,osid_cm)
        if param_vary=="ct" :
            ct=add_fdiff(ct,dct,sign_vary,osid_ct)
        if param_vary=="zh" :
            zh=add_fdiff(zh,dzh,sign_vary,osid_zh)
        if param_vary=="ezh" :
            ezh=add_fdiff(ezh,dezh,sign_vary,osid_ezh)
        if param_vary=="m2i" :
            m2i=add_fdiff(m2i,dm2i,sign_vary,osid_m2i)
        if param_vary=="kv" :
            kv=add_fdiff(kv,dkv,sign_vary,osid_kv)
    if param_vary=="ns" :
        ns=add_fdiff(ns,dns,sign_vary,osid_ns)
    if param_vary=="A_s" :
        a_s=add_fdiff(a_s,da_s,sign_vary,osid_a_s)
    if param_vary=="tau" :
        tau=add_fdiff(tau,dtau,sign_vary,osid_tau)
    if param_vary=="mnu" :
        mnu=add_fdiff(mnu,dmnu,sign_vary,osid_mnu)
    if param_vary=="pan" :
        pan=add_fdiff(pan,dpan,sign_vary,osid_pan)
    if param_vary=="fnl" :
        fNL=add_fdiff(fNL,dfNL,sign_vary,osid_fNL)
    if param_vary=="rt" :
        rt=add_fdiff(rt,drt,sign_vary,osid_rt)

    nuisance_name,tr_name,inode=my_parser(param_vary)

    n_tracers_nc=0
    n_tracers_wl=0
    photoz_nc_string=" "
    photoz_wl_string=" "
    selection_nc_string=" "
    selection_wl_string=" "
    bins_nc_string=" "
    bins_wl_string=" "
    nz_nc_string=" "
    nz_wl_string=" "
    bz_string=" "
    sz_string=" "
    ez_string=" "
    az_string=" "
    rf_string=" "
    for i in np.arange(len(par.tracers)) :
        tr=par.tracers[i]
        if tr.name==tr_name and nuisance_name=="sphz" :
                bins_fname =tr.nuisance_sphz.get_filename(inode,sign_vary)+"_bins.txt"
        elif tr.name==tr_name and nuisance_name=="bphz" :
                bins_fname =tr.nuisance_bphz.get_filename(inode,sign_vary)+"_bins.txt"
        else :
            bins_fname =tr.nuisance_sphz.get_filename(-1,sign_vary)+"_bins.txt"
        if tr.name==tr_name and nuisance_name=="bias" :
            bias_fname =tr.nuisance_bias.get_filename(inode,sign_vary)
        else :
            bias_fname =tr.nuisance_bias.get_filename(-1,sign_vary)
        if tr.name==tr_name and nuisance_name=="sbias" :
            sbias_fname =tr.nuisance_sbias.get_filename(inode,sign_vary)
        else :
            sbias_fname =tr.nuisance_sbias.get_filename(-1,sign_vary)
        if tr.name==tr_name and nuisance_name=="ebias" :
            ebias_fname =tr.nuisance_ebias.get_filename(inode,sign_vary)
        else :
            ebias_fname =tr.nuisance_ebias.get_filename(-1,sign_vary)
        if tr.name==tr_name and nuisance_name=="abias" :
            abias_fname =tr.nuisance_abias.get_filename(inode,sign_vary)
        else :
            abias_fname =tr.nuisance_abias.get_filename(-1,sign_vary)
        if tr.name==tr_name and nuisance_name=="rfrac" :
            rfrac_fname =tr.nuisance_rfrac.get_filename(inode,sign_vary)
        else :
            rfrac_fname =tr.nuisance_rfrac.get_filename(-1,sign_vary)

        if tr.tracer_type=="gal_clustering" :
            photoz_nc_string+="%d "%(tr.is_photometric)
            selection_nc_string+="tophat "
            bins_nc_string+=bins_fname+" "
            nz_nc_string+=tr.nz_file+" "
            bz_string+=bias_fname+" "
            sz_string+=sbias_fname+" "
            ez_string+=ebias_fname+" "
            n_tracers_nc+=1
        elif tr.tracer_type=="intensity_mapping" :
            photoz_nc_string+="0 "
            selection_nc_string+="tophat "
            bins_nc_string+=bins_fname+" "
            nz_nc_string+=tr.nz_file+" "
            bz_string+=bias_fname+" "
            sz_string+=sbias_fname+" "
            ez_string+=ebias_fname+" "
            n_tracers_nc+=1
        elif tr.tracer_type=="gal_shear" :
            photoz_wl_string+="1 "
            selection_wl_string+="tophat "
            bins_wl_string+=bins_fname+" "
            nz_wl_string+=tr.nz_file+" "
            az_string+=abias_fname+" "
            rf_string+=rfrac_fname+" "
            n_tracers_wl+=1

    spectra_list="mPk"
    if par.has_gal_clustering==True :
        spectra_list+=", nCl"
    if par.has_gal_shear==True :
        spectra_list+=", sCl"
    if (par.has_cmb_lensing==True) or (par.has_cmb_t==True) or (par.has_cmb_p==True) :
        spectra_list+=", lCl"
    if (par.has_cmb_t==True) or (par.has_cmb_p==True) :
        spectra_list+=", tCl, pCl"

    strout="#CLASS param file by GoFish\n"
    strout+="h = %lE\n"%hh
    strout+="T_cmb = 2.725\n"
    strout+="omega_b = %lE\n"%obh2
    strout+="omega_cdm = %lE\n"%och2
    if mnu>0 :
        strout+="N_ur = 2.0328\n"
        strout+="N_ncdm = 1\n"
        strout+="m_ncdm = %lE \n"%(mnu*0.001)
#        strout+="N_ur = 1.0196\n"
#        strout+="N_ncdm = 2\n"
#        strout+="m_ncdm = %lE, "%(mnu*0.001/2)+"%lE \n"%(mnu*0.001/2)
#        strout+="N_ur = 0.00641\n"
#        strout+="N_ncdm = 3\n"
#        strout+="m_ncdm = %lE, "%(mnu*0.001/3)+"%lE, "%(mnu*0.001/3)+"%lE \n"%(mnu*0.001/3)
    else :
        strout+="N_ur = 3.046\n"
    if ((par.model=='LCDM') or (par.model=='wCDM')) :
        if (w0!=-1.) or (wa!=0.0):
            strout+="Omega_fld = %lE\n"%(1-(och2+obh2)/hh**2)
            if wa<0 :
                strout+="w0_fld = %lE\n"%(w0-0.00001)
            elif wa>0 :
                strout+="w0_fld = %lE\n"%(w0+0.00001)
            else :
                strout+="w0_fld = %lE\n"%w0
            strout+="wa_fld = %lE\n"%wa
            strout+="cs2_fld = 1\n"
        strout+="Omega_k = 0.\n"
    if par.model=='JBD' :
        strout+="Omega_Lambda = 0.\n"
        strout+="Omega_fld = 0\n"
        strout+="Omega_smg= -1\n"
        strout+="gravity_model = brans dicke\n"
        strout+="parameters_smg = 0.7, %lE, 1., 0.\n"%(1E4/obd)
        strout+="tuning_index_smg = 0\n"
        strout+="tuning_dxdy_guess_smg = 1e-1\n"
        strout+="friedmann_branch_smg = 0\n"
        strout+="kineticity_safe_smg = 0\n" #1E-4
        strout+="cs2_safe_smg = 1e-10\n"
        strout+="skip_stability_tests_smg = no\n"
    if par.model=='Horndeski' :
        strout+="Omega_Lambda = 1E-100\n"
        strout+="Omega_fld = 0\n"
        strout+="Omega_smg= -1\n"
        if (w0!=-1.) or (wa!=0.0):
            strout+="expansion_model = cpl\n"
            strout+="expansion_smg = %lE, "%(1-(och2+obh2)/hh**2)
            if wa<0 :
                strout+="%lE, "%(w0-0.00001)
            elif wa>0 :
                strout+="%lE, "%(w0+0.00001)
            else :
                strout+="%lE, "%w0
            strout+="%lE\n"%wa
        else :
            strout+="expansion_model= lcdm\n"
            strout+="expansion_smg= %lE\n"%(1-(och2+obh2)/hh**2)
        strout+="gravity_model = threshold_alphas\n"
        strout+="parameters_smg = "
        strout+="%lE, "%bk+"%lE, "%bb+"%lE, "%bm+"%lE, "%bt
        strout+="%lE, "%ck+"%lE, "%cb+"%lE, "%cm+"%lE, "%ct
        strout+="%lE, "%zh+"%lE\n"%ezh
        strout+="kineticity_safe_smg = 1E-4\n" #NEW HICLASS
        strout+="skip_stability_tests_smg = no\n" #NEW HICLASS
        strout+="k_vainshtein = %lE\n"%kv

    strout+="f_NL = %lE\n"%fNL

    strout+="YHe = BBN\n"
    strout+="recombination = RECFAST\n"
    strout+="reio_parametrization = reio_camb\n"
    strout+="tau_reio = %lE\n"%tau
    strout+="annihilation = %lE\n"%(pan*1E-6)
    strout+="output = "+spectra_list+"\n"
    strout+="number count contributions = "+par.terms_gc+"\n"
    strout+="weak lensing contributions = "+par.terms_gs+"\n"
    if par.use_nonlinear==True :
        strout+="non linear = halofit\n"
    if (par.has_cmb_lensing==True) or (par.has_cmb_t==True) or (par.has_cmb_p==True) :
        strout+="modes = s, t\n"
        strout+="lensing = yes\n"
    else :
        strout+="modes = s\n"
    strout+="ic = ad\n"
    strout+="gauge = synchronous\n"
    strout+="P_k_ini type = analytic_Pk\n"
    strout+="k_pivot = 0.05\n"
    strout+="A_s = %lE\n"%(a_s*1E-9)
    strout+="n_s = %lE\n"%ns
    strout+="alpha_s = 0.\n"
    strout+="r = %lE\n"%(rt)
    strout+="n_t = scc\n"
    strout+="alpha_t = scc\n"
    strout+="l_max_scalars = %d\n"%(par.lmax_cmb)
    strout+="l_max_tensors = %d\n"%(par.lmax_cmb)
    strout+="l_max_lss = %d\n"%(par.lmax_lss)
    strout+="P_k_max_h/Mpc = 1.\n"
    strout+="z_pk = 0\n"
    strout+="number_of_tracers_nc = %d\n"%n_tracers_nc
    strout+="number_of_tracers_wl = %d\n"%n_tracers_wl
    strout+="selection_nc = "+selection_nc_string+"\n"
    strout+="selection_wl = "+selection_wl_string+"\n"
    strout+="selection_bins_nc = "+bins_nc_string+"\n"
    strout+="selection_bins_wl = "+bins_wl_string+"\n"
    strout+="use_photoz_nc = "+photoz_nc_string+"\n"
    strout+="use_photoz_wl = "+photoz_wl_string+"\n"
    strout+="dNdz_selection_nc ="+nz_nc_string+"\n"
    strout+="dNdz_selection_wl ="+nz_wl_string+"\n"
    strout+="bias_function = "+bz_string+"\n"
    strout+="s_bias_function = "+sz_string+"\n"
    strout+="e_bias_function = "+ez_string+"\n"
    strout+="alignment_bias_function = "+az_string+"\n"
    strout+="red_fraction = "+rf_string+"\n"
    strout+="root = "+prefix_out+"\n"
    strout+="output_binary = 1\n"
    strout+="format = class\n"
    strout+="write background = yes\n"
    strout+="write thermodynamics = no\n"
    strout+="write primordial = no\n"
    strout+="write parameters = yeap\n"
    strout+="l_switch_limber= 10.\n"   #%(par.lmin_limber)
    strout+="l_switch_limber_for_cl_density= %.1lf\n"%(par.lmin_limber)
    strout+="l_switch_limber_for_cl_lensing= %.1lf\n"%(par.lmin_limber)
    strout+="selection_sampling_bessel=6.\n"
    strout+="k_step_trans_scalars=0.4\n"
    strout+="q_linstep=0.4\n"
    strout+="k_scalar_max_tau0_over_l_max= 2.\n"
    strout+="selection_tophat_edge= 0.01\n"
    strout+="z_min_k= 0.1\n"
    strout+="perturb_sampling_stepsize=0.06\n"
    strout+="write warnings = y\n"
    strout+="input_verbose = 1\n"
    strout+="background_verbose = 1\n"
    strout+="thermodynamics_verbose = 1\n"
    strout+="perturbations_verbose = 1\n"
    strout+="transfer_verbose = 1\n"
    strout+="primordial_verbose = 1\n"
    strout+="spectra_verbose = 1\n"
    strout+="nonlinear_verbose = 1\n"
    strout+="lensing_verbose = 1\n"
    strout+="output_verbose = 1\n"
    f=open(prefix_out+"_param.ini","w")
    f.write(strout)
    f.close()

    return strout

def get_prefix(par,par_vary,sign_vary) :
    if par_vary=="none" :
        prefix_all=par.output_path+"_fid"
    else :
        if sign_vary==1 :
            prefix_all=par.output_path+"_p"+par_vary
        elif sign_vary==-1 :
            prefix_all=par.output_path+"_m"+par_vary
        else :
            sys.exit("Wrong sign!")

    return prefix_all

def get_bao_names(par,par_vary,sign_vary) :
    prefix_all=get_prefix(par,par_vary,sign_vary)
    return prefix_all+"background.dat"

def get_cl_names(par,par_vary,sign_vary) :
    prefix_all=get_prefix(par,par_vary,sign_vary)
    if par.use_CAMB :
        clfile_total =prefix_all+"_lenspotentialCls.dat"
        clfile_lensed=prefix_all+"_lensedtotCls.dat"
        clfile_scalar=prefix_all+"_scalCls.dat"
        clfile_tensor=prefix_all+"_tensCls.dat"
    else :
        clfile_total=prefix_all+"cl.dat"
        clfile_lensed=prefix_all+"cl_lensed.dat"
        clfile_scalar=prefix_all+"cls.dat"
        clfile_tensor=prefix_all+"clt.dat"

    return clfile_total,clfile_lensed,clfile_scalar,clfile_tensor

def cls_are_there(par,par_vary,sign_vary,printout) :
    """ Generate and read power spectra with CLASS """
    
    clf_total,clf_lensed,clf_scalar,clf_tensor=get_cl_names(par,par_vary,sign_vary)
    if (os.path.isfile(clf_total)==False) :
        if printout :
            print "   Couldn't find cl for "+par_vary+" %d"%sign_vary
        return False
    else :
        return True

def bao_is_there(par,par_vary,sign_vary,printout) :
    """ Generate and read power spectra with CLASS """
    
    fname_bao=get_bao_names(par,par_vary,sign_vary)
    if (os.path.isfile(fname_bao)==False) :
        if printout :
            print "   Couldn't find BAO for "+par_vary+" %d "%sign_vary+fname_bao
        return False
    else :
        return True

def start_running(par,par_vary,sign_vary) :
    """ Start running CLASS if the files aren't there """
    checkvar=False
    if par.n_tracers>0 :
        checkvar=cls_are_there(par,par_vary,sign_vary,True)
    else :
        checkvar=bao_is_there(par,par_vary,sign_vary,True)

    if checkvar :
        return True
    else :
        print "     Computing"
        prefix_all=get_prefix(par,par_vary,sign_vary)
        if par.use_CAMB :
            write_camb_param_file(par,par_vary,sign_vary,prefix_all)
        else :
            write_class_param_file(par,par_vary,sign_vary,prefix_all)
        os.system(par.exec_path+" "+prefix_all+"_param.ini ")
        if par.use_CAMB :
            os.system("mv Erminia_background.dat "+get_bao_names(par,par_vary,sign_vary))
        return False

def get_bao(par,par_vary,sign_vary) :
    fname_bao=get_bao_names(par,par_vary,sign_vary)

    data=np.loadtxt(fname_bao,unpack=True)
    h_fid,dum1,dum2=par.get_param_properties("hh")
    if par.use_CAMB :
        zarr=data[0]; daarr=data[1]; hharr=data[2]; dvarr=data[3];
    else :
        zarr=data[0]; hharr=data[3]; daarr=data[5]; dvarr=data[5];
    dafunc=interp1d(zarr,daarr*h_fid)
    hhfunc=interp1d(zarr,hharr/h_fid)
    dvfunc=interp1d(zarr,dvarr/h_fid)

    daout=np.array([dafunc(z) for z in par.z_nodes_DA])
    hhout=np.array([hhfunc(z) for z in par.z_nodes_HH])
    dvout=np.array([dvfunc(z) for z in par.z_nodes_DV])

    return daout,hhout,dvout

def get_cls(par,par_vary,sign_vary) :
    """ Generate and read power spectra with CLASS """
    prefix_all=get_prefix(par,par_vary,sign_vary)
    clf_total,clf_lensed,clf_scalar,clf_tensor=get_cl_names(par,par_vary,sign_vary)

    dict_t=read_cls(par,clf_total)
    #If we have CMB lensing, and if this is the fiducial run, we must also read Cl_CMB in order to compute the noise
    for tr in par.tracers :
        if (tr.tracer_type=='cmb_lensing') and tr.consider_tracer and (par_vary=="none"):
            dict_l=read_cls(par,clf_lensed)
            tr.cl_tt_u=dict_t['cl_tt']; tr.cl_ee_u=dict_t['cl_ee'];
            tr.cl_te_u=dict_t['cl_te']; tr.cl_bb_u=dict_t['cl_bb'];
            tr.cl_tt_l=dict_l['cl_tt']; tr.cl_ee_l=dict_l['cl_ee'];
            tr.cl_te_l=dict_l['cl_te']; tr.cl_bb_l=dict_l['cl_bb'];

    if dict_t['lmax']<par.lmax :
        sys.exit("Error reading cls")
    if par.has_gal_clustering :
        if dict_t['n_nc']!=par.nbins_gal_clustering_read :
            sys.exit("Error reading d cls %d"%(dict_t['n_nc']))
    if par.has_gal_shear :
        if dict_t['n_wl']!=par.nbins_gal_shear_read :
            sys.exit("Error reading l cls %d"%(dict_t['n_wl']))
    if par.has_cmb_t or par.has_cmb_p :
        dict_l=read_cls(par,clf_lensed)
    else :
        dict_l=dict_t
    cl_tt=dict_t['cl_tt']; cl_ee=dict_t['cl_ee']; cl_te=dict_t['cl_te']; cl_bb=dict_t['cl_bb'];
#    cl_tt=dict_l['cl_tt']; cl_ee=dict_l['cl_ee']; cl_te=dict_l['cl_te']; cl_bb=dict_l['cl_bb'];
    if par.has_cmb_p :
        cl_bb*=0
    cl_pp=dict_t['cl_pp']; cl_tp=dict_t['cl_tp']; cl_ep=dict_t['cl_ep'];
    cl_td=dict_t['cl_td']; cl_tl=dict_t['cl_tl']; cl_pd=dict_t['cl_pd']; cl_pl=dict_t['cl_pl'];
    cl_dd=dict_t['cl_dd']; cl_dl=dict_t['cl_dl']; cl_ll=dict_t['cl_ll'];

    cl_ret=np.zeros([par.lmax+1,par.nbins_total,par.nbins_total])

    nb1_sofar=0
    nbd1_sofar=0
    nbl1_sofar=0
    nbp1_sofar=0
    nbc1_sofar=0
    for tr1 in par.tracers :
        nb2_sofar=0
        nbd2_sofar=0
        nbl2_sofar=0
        nbp2_sofar=0
        nbc2_sofar=0
        nb1=tr1.nbins
        tr1clust=((tr1.tracer_type=="gal_clustering") or
                  (tr1.tracer_type=="intensity_mapping"))
        for tr2 in par.tracers :
            nb2=tr2.nbins
            lmx=np.amin(np.array([par.lmax,tr1.lmax,tr2.lmax]))
            lmn=np.amax(np.array([tr1.lmin,tr2.lmin]))
            tr2clust=((tr2.tracer_type=="gal_clustering") or
                      (tr2.tracer_type=="intensity_mapping"))
            if tr1.consider_tracer and tr2.consider_tracer :
                if (tr1.tracer_type=="cmb_primary") and (tr2.tracer_type=="cmb_primary") : #cc
                    if (tr1.has_t==True) and (tr1.has_p==False) :
                        cl_ret[lmn:lmx+1,nb1_sofar+0,nb2_sofar+0]=cl_tt[lmn:lmx+1]
                    if (tr1.has_t==False) and (tr1.has_p==True) :
                        cl_ret[lmn:lmx+1,nb1_sofar+0,nb2_sofar+0]=cl_ee[lmn:lmx+1]
                        cl_ret[lmn:lmx+1,nb1_sofar+1,nb2_sofar+1]=cl_bb[lmn:lmx+1]
                    if (tr1.has_t==True) and (tr1.has_p==True) :
                        cl_ret[lmn:lmx+1,nb1_sofar+0,nb2_sofar+0]=cl_tt[lmn:lmx+1]
                        cl_ret[lmn:lmx+1,nb1_sofar+0,nb2_sofar+1]=cl_te[lmn:lmx+1]
                        cl_ret[lmn:lmx+1,nb1_sofar+1,nb2_sofar+0]=cl_te[lmn:lmx+1]
                        cl_ret[lmn:lmx+1,nb1_sofar+1,nb2_sofar+1]=cl_ee[lmn:lmx+1]
                        cl_ret[lmn:lmx+1,nb1_sofar+2,nb2_sofar+2]=cl_bb[lmn:lmx+1]
                if (tr1.tracer_type=="cmb_lensing") and (tr2.tracer_type=="cmb_lensing") : #pp
                    cl_ret[lmn:lmx+1,nb1_sofar,nb2_sofar]=cl_pp[lmn:lmx+1]
                if (tr1.tracer_type=="cmb_lensing") and (tr2clust) : #pd
                    cl_ret[lmn:lmx+1,nb1_sofar,nb2_sofar:nb2_sofar+nb2]=cl_pd[lmn:lmx+1,nbd2_sofar:nbd2_sofar+nb2]
                if (tr1clust) and (tr2.tracer_type=="cmb_lensing") : #dp
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,nb2_sofar]=cl_pd[lmn:lmx+1,nbd1_sofar:nbd1_sofar+nb1]
                if (tr1.tracer_type=="cmb_lensing") and (tr2.tracer_type=="gal_shear") : #pl
                    cl_ret[lmn:lmx+1,nb1_sofar,nb2_sofar:nb2_sofar+nb2]=cl_pl[lmn:lmx+1,nbl2_sofar:nbl2_sofar+nb2]
                if (tr1.tracer_type=="gal_shear") and (tr2.tracer_type=="cmb_lensing") : #lp
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,nb2_sofar]=cl_pl[lmn:lmx+1,nbl1_sofar:nbl1_sofar+nb1]
                if (tr1clust) and (tr2clust) : #dd
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,
                           nb2_sofar:nb2_sofar+nb2]=cl_dd[lmn:lmx+1,nbd1_sofar:nbd1_sofar+nb1,
                                                          nbd2_sofar:nbd2_sofar+nb2]
                if (tr1clust) and (tr2.tracer_type=="gal_shear") : #dl
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,
                             nb2_sofar:nb2_sofar+nb2]=cl_dl[lmn:lmx+1,nbd1_sofar:nbd1_sofar+nb1,
                                                             nbl2_sofar:nbl2_sofar+nb2]
                if (tr1.tracer_type=="gal_shear") and (tr2clust) : #ld
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,
                             nb2_sofar:nb2_sofar+nb2]=np.transpose(cl_dl[lmn:lmx+1,nbd2_sofar:nbd2_sofar+nb2,
                                                                          nbl1_sofar:nbl1_sofar+nb1],axes=(0,2,1))
                if (tr1.tracer_type=="gal_shear") and (tr2.tracer_type=="gal_shear") : #ll
                    cl_ret[lmn:lmx+1,nb1_sofar:nb1_sofar+nb1,
                             nb2_sofar:nb2_sofar+nb2]=cl_ll[lmn:lmx+1,nbl1_sofar:nbl1_sofar+nb1,
                                                             nbl2_sofar:nbl2_sofar+nb2]
            if tr2clust :
                nbd2_sofar+=nb2
            if tr2.tracer_type=="gal_shear" :
                nbl2_sofar+=nb2
            if tr2.tracer_type=="cmb_lensing" :
                nbp2_sofar+=nb2
            if tr2.tracer_type=="cmb_primary" :
                nbc2_sofar+=nb2
            if tr2.consider_tracer :
                nb2_sofar+=nb2
        if tr1clust :
            nbd1_sofar+=nb1
        if tr1.tracer_type=="gal_shear" :
            nbl1_sofar+=nb1
        if tr1.tracer_type=="cmb_lensing" :
            nbp1_sofar+=nb1
        if tr1.tracer_type=="cmb_primary" :
            nbc1_sofar+=nb1
        if tr1.consider_tracer :
            nb1_sofar+=nb1

    if par.save_cl_files==False :
        os.system("rm -f "+clf_total)
        os.system("rm -f "+clf_lensed)
        os.system("rm -f "+clf_scalar)
        os.system("rm -f "+clf_tensor)
    if par.save_param_files==False :
        os.system("rm -f "+prefix_all+"_param.ini")
    if par.save_dbg_files==False :
        os.system("rm -f bg_test.dat Tkout.dat transfer_info.txt")
        os.system("rm -f "+prefix_all+"parameters.ini "+prefix_all+"unused_parameters")

    return cl_ret
