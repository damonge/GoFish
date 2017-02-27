import numpy as np
import matplotlib.pyplot as plt
import in_out as ino
import ConfigParser
import os as os
import sys as sys
import tracers as trc
import fisher_plot as fsh
import time

fs=16
lw=2

PARS_LCDM={'och2':[0.1197 ,0.001 , 0,'$\\omega_c$'],
           'obh2':[0.02222,0.0001, 0,'$\\omega_b$'],
           'hh'  :[0.67   ,0.01  , 0,'$h$'],
           'ns'  :[0.96   ,0.01  , 0,'$n_s$'],
           'A_s' :[2.19   ,0.01  , 0,'$A_s$'],
           'tau' :[0.06   ,0.005 , 0,'$\\tau$'],
           'mnu' :[60.0   ,5.0   , 0,'$\\sum m_\\nu$'],
           'nnu' :[2.046  ,0.1   , 0,'$\\rm{N_{eff}}$'],
           'pan' :[0.00   ,0.5   , 1,'$p_{\\rm ann}$'],
           'fnl' :[0.00   ,0.5   , 0,'$f_{\\rm NL}$'],
           'rt'  :[0.00   ,0.001 , 1,'$r$']}

PARS_BADE={'w0'  :[-0.99  ,0.01  , 0,'$w_0$'],
           'wdeq':[-0.99  ,0.01  , 0,'$w^B_0$'],
           'deq' :[0.     ,0.01  , 1,'$\\Omega_e^B$']}

PARS_EDE= {'w0'  :[-0.99  ,0.01  , 0,'$w_0$'],
           'wdeq':[-0.99  ,0.01  , 0,'$w^B_0$'],
           'ede' :[0.     ,0.01  , 1,'$\\Omega_e$']}

PARS_WCDM={'w0'  :[-1.00  ,0.01  , 1,'$w_0$'],
           'wa'  :[0.00   ,0.01  , 1,'$w_a$']}

PARS_JBD ={'obd' :[0.10   ,0.03  , 0,'$10^4/\\omega_{\\rm BD}$']}

PARS_HORN={'bk'  :[ 1E-5  ,0.05  , 1,'$b_K$'],
           'bb'  :[-1E-5  ,0.01  ,-1,'$b_B$'],
           'bm'  :[ 1E-5  ,0.007 , 1,'$b_M$'],
           'bt'  :[-1E-5  ,0.02  ,-1,'$b_T$'],
           'ck'  :[ 0.10  ,0.05  , 0,'$c_K$'],
           'cb'  :[ 0.05  ,0.01  , 0,'$c_B$'],
           'cm'  :[-0.05  ,0.007 , 0,'$c_M$'],
           'ct'  :[-0.05  ,0.01  ,-1,'$c_T$'],
           'zh'  :[ 0.5   ,0.05  , 0,'$z_H$'],
           'ezh' :[ 0.5   ,0.05  , 0,'$\\Delta z_H$'],
           'm2i' :[ 1.0   ,0.05  , 0,'$M^2_*$'],
           'kv'  :[ 0.1   ,0.02  , 0,'$k_V$']}

LMAX=10000
LMAX_CMB=10000

NLB=1


class ParamRun:
    """ Run parameters """
    #Cosmological parameters
    params_all=[] #

    #Run parameters
    lmax=LMAX #
    lmax_lss=LMAX #
    lmax_cmb=LMAX_CMB
    lmin_limber=LMAX #
    has_gal_clustering=False
    has_gal_shear=False
    has_cmb_lensing=False
    has_cmb_t=False
    has_cmb_p=False
    output_dir="test_dir" #
    output_spectra="test_prefix" #
    output_fisher="test_fisher" #
    output_path="test_dir/test_prefix" #
    exec_path="./classt" #
    fsky=1.0
    use_CAMB=False
    just_run_cls=False

    #Bin parameters
    nbins_gal_clustering=0 #
    nbins_gal_clustering_read=0 #
    nbins_gal_shear=0 #
    nbins_gal_shear_read=0 #
    nbins_cmb_primary=0 #
    nbins_cmb_primary_read=0 #
    nbins_cmb_lensing=0 #
    nbins_cmb_lensing_read=0 #
    nbins_total=1 #
    n_tracers=0 #
    tracers=[] #

    #Behavioral flags
    model='LCDM' #
    save_cl_files=True #
    save_param_files=True #
    save_dbg_files=False
    include_alignment=False #
    include_rsd=True #
    include_magnification=True #
    include_gr_vel=False #
    include_gr_pot=False #
    use_nonlinear=False
    plot_ext=".png"
    include_im_fg=False
    fit_im_fg=False

    #BAO
    n_bao=0
    include_BAO=False
    include_DA=False
    include_HH=False
    include_DV=False
    z_nodes_DA=[]
    z_nodes_HH=[]
    z_nodes_DV=[]
    e_nodes_DA=[]
    e_nodes_HH=[]
    e_nodes_DV=[]

    #Terms to include for clustering
    terms_gc="density" #
    terms_gs="lensing_shear" #

    #Power spectra
    cl_noise_arr=[]
    cl_fid_arr=[]
    dcl_arr=[]

    #Fisher param structure
    params_fshr=[] #

    #Number of parameters to vary
    npar_vary=0

    #Fisher matrix
    fshr_l=[]
    fshr_cls=[]
    fshr_bao=[]
    fshr=[]

    def __init__(self,fname) :
        #Read parameter file
        self.read_param_file(fname)

        if (self.has_cmb_lensing==False) and (self.has_cmb_t==False) and (self.has_cmb_p==False) :
            self.lmax_cmb=0
        if (self.has_gal_shear==False) and (self.has_gal_clustering==False) :
            self.lmax_lss=0
        self.lmax=max(self.lmax_lss,self.lmax_cmb)
        self.lmax=NLB*((self.lmax+1)/NLB)-1

        #List parameters to vary over
        self.npar_vary=0
        for p in self.params_all :
            if p.isfree==True :
                self.params_fshr.append(fsh.ParamFisher(p.val,p.dval,p.name,
                                                        p.label,p.isfree,p.do_plot,p.onesided))
                self.npar_vary+=1
        self.params_fshr=np.array(self.params_fshr)

        #Read bins information for each tracer
        self.nbins_total=0
        self.nbins_gal_clustering=0
        self.nbins_gal_clustering_read=0
        self.nbins_gal_shear=0
        self.nbins_gal_shear_read=0
        self.nbins_cmb_lensing=0
        self.nbins_cmb_lensing_read=0
        self.nbins_cmb_primary=0
        self.nbins_cmb_primary_read=0
        for tr in self.tracers :
            nbins=tr.get_nbins()
            if (tr.tracer_type=="gal_clustering") or (tr.tracer_type=="intensity_mapping") :
                self.nbins_gal_clustering_read+=nbins
                if tr.consider_tracer :
                    self.nbins_gal_clustering+=nbins
            elif tr.tracer_type=="gal_shear" :
                self.nbins_gal_shear_read+=nbins
                if tr.consider_tracer :
                    self.nbins_gal_shear+=nbins
            elif tr.tracer_type=="cmb_lensing" :
                if self.nbins_cmb_lensing==1 :
                    sys.exit("You can only have 1 CMB lensing!")
                else :
                    self.nbins_cmb_lensing_read=1
                    if tr.consider_tracer :
                        self.nbins_cmb_lensing=1
            elif tr.tracer_type=="cmb_primary" :
                if self.nbins_cmb_primary>0 :
                    sys.exit("You can only have 1 CMB primary!")
                else :
                    if tr.has_t :
                        self.nbins_cmb_primary_read+=1
                    if tr.has_p :
                        self.nbins_cmb_primary_read+=2
                    if tr.consider_tracer :
                        if tr.has_t :
                            self.nbins_cmb_primary+=1
                        if tr.has_p :
                            self.nbins_cmb_primary+=2
            else :
                strout="Wrong tracer type "+tr.tracer_type+"\n"
                strout+="Allowed tracer types are: gal_clustering, intensity_mapping"
                strout+=", gal_shear, cmb_lensing, cmb_primary\n"
                sys.exit(strout)
            if tr.consider_tracer :
                self.nbins_total+=nbins

        #Terms to include
        if self.include_alignment==True :
            self.terms_gs+=", intrinsic_alignment"
        if self.include_rsd==True :
            self.terms_gc+=", rsd1"
        if self.include_magnification==True :
            self.terms_gc+=", lensing"
        if self.include_gr_vel==True :
            self.terms_gc+=", rsd2, rsd3"
        if self.include_gr_pot==True :
            self.terms_gc+=", gr1, gr2, gr3, gr4, gr5"

        #Allocate power spectra
        self.cl_fid_arr=np.zeros([(self.lmax+1)/NLB,
                                  self.nbins_total,self.nbins_total])
        self.cl_noise_arr=np.zeros([(self.lmax+1)/NLB,
                                    self.nbins_total,self.nbins_total])
        self.dcl_arr=np.zeros([self.npar_vary,(self.lmax+1)/NLB,
                               self.nbins_total,self.nbins_total])

        #Allocate Fisher matrix
        self.fshr_l=np.zeros([self.npar_vary,self.npar_vary,self.lmax+1])
        self.fshr_cls=np.zeros([self.npar_vary,self.npar_vary])
        self.fshr_bao=np.zeros([self.npar_vary,self.npar_vary])
        self.fshr=np.zeros([self.npar_vary,self.npar_vary])
        self.print_params()

    def print_params(self) :
        print " <><> Run parameters"
        print "  - Overall l_max = %d"%(self.lmax)
        if self.has_gal_clustering or self.has_gal_shear :
            print "  - LSS l_max = %d"%(self.lmax_lss)
        if self.has_cmb_lensing or self.has_cmb_t or self.has_cmb_p :
            print "  - CMB l_max = %d"%(self.lmax_cmb)
        print "  - Output directory : "+self.output_dir
        print "  - Output spectra : "+self.output_path
        print "  - Output Fisher : "+self.output_fisher
        if self.has_gal_clustering :
            print "  - %d"%(self.nbins_gal_clustering)+" (%d) bins used (read) for galaxy clustering"%(self.nbins_gal_clustering_read)
        if self.has_gal_shear :
            print "  - %d"%(self.nbins_gal_shear)+" (%d) bins used (read) for galaxy shear"%(self.nbins_gal_shear_read)
        if self.has_cmb_t or self.has_cmb_p :
            print "  - %d"%(self.nbins_cmb_primary)+" (%d) bins used (read) for CMB primary"%(self.nbins_cmb_primary_read)
        if self.has_cmb_lensing :
            print "  - %d"%(self.nbins_cmb_lensing)+" (%d) bins used (read) for CMB lensing"%(self.nbins_cmb_lensing_read)
        print "  - %d bins used in total"%(self.nbins_total)
        print "  - %d tracers read : "%(self.n_tracers)
        itr=1
        for tr in self.tracers :
            print "     %d. "%itr+tr.name+". Type: "+tr.tracer_type+", %d bins."%(tr.nbins)
            itr+=1
        print "  - Overlaping sky fraction : %.3lf"%(self.fsky)
        if self.has_gal_clustering :
            print "  - Terms included for galaxy clustering : "+self.terms_gc
        if self.has_gal_shear :
            print "  - Terms included for galaxy shear : "+self.terms_gs
        print "  - %d free parameters :"%(self.npar_vary)
        for par in self.params_fshr :
            print "    * "+par.name+" : %.3lf"%(par.val)

    def read_param_file(self,fname) :
        config=ConfigParser.SafeConfigParser()
        config.read(fname)

        #Behaviour parameters
        if config.has_option('Behaviour parameters','model') :
            self.model=config.get('Behaviour parameters','model')
        if config.has_option('Behaviour parameters','use_CAMB') :
            self.use_CAMB=config.getboolean('Behaviour parameters','use_CAMB')
        if config.has_option('Behaviour parameters','save_cl_files') :
            self.save_cl_files=config.getboolean('Behaviour parameters','save_cl_files')
        if config.has_option('Behaviour parameters','save_param_files') :
            self.save_param_files=config.getboolean('Behaviour parameters',
                                                    'save_param_files')
        if config.has_option('Behaviour parameters','just_run_cls') :
            self.just_run_cls=config.getboolean('Behaviour parameters',
                                                    'just_run_cls')

        if ((self.model=='BADE') and (self.use_CAMB==False)) :
            print "Can only do BADE with CAMB"
            exit(1)

        if ((self.model=='EDE') and (self.use_CAMB==False)) :
            print "Can only do EDE with CAMB"
            exit(1)

        #Cosmological parameters
        self.params_all=[]
        def add_to_params(pars) :
            for pname in pars.keys() :
                x=pars[pname][0]
                dx=pars[pname][1]
                isfree=False
                onesided=pars[pname][2]
                if config.has_section(pname) :
                    if config.has_option(pname,'x') :
                        x=config.getfloat(pname,'x')
                    if config.has_option(pname,'dx') :
                        dx=config.getfloat(pname,'dx')
                    if config.has_option(pname,'is_free') :
                        isfree=config.getboolean(pname,'is_free')
                    if config.has_option(pname,'onesided') :
                        onesided=config.getint(pname,'onesided')
                self.params_all.append(fsh.ParamFisher(x,dx,pname,pars[pname][3],isfree,isfree*True,onesided))
        add_to_params(PARS_LCDM)
        if self.model=='LCDM' :
            pass
        if self.model=='BADE' :
            add_to_params(PARS_BADE)
        if self.model=='EDE' :
            add_to_params(PARS_EDE)
        if self.model=='wCDM' :
            add_to_params(PARS_WCDM)
        elif self.model=='JBD' :
            add_to_params(PARS_JBD)
        elif self.model=='Horndeski' :
            add_to_params(PARS_WCDM)
            add_to_params(PARS_HORN)

        #CLASS parameters
        if config.has_option('CLASS parameters','lmax_lss') :
            self.lmax_lss=config.getint('CLASS parameters','lmax_lss')
        if config.has_option('CLASS parameters','lmax_cmb') :
            self.lmax_cmb=config.getint('CLASS parameters','lmax_cmb')
        if config.has_option('CLASS parameters','lmin_limber') :
            self.lmin_limber=config.getfloat('CLASS parameters','lmin_limber')
        if config.has_option('CLASS parameters','include_alignment') :
            self.include_alignment=config.getboolean('CLASS parameters','include_alignment')
        if config.has_option('CLASS parameters','include_rsd') :
            self.include_rsd=config.getboolean('CLASS parameters','include_rsd')
        if config.has_option('CLASS parameters','include_magnification') :
            self.include_magnification=config.getboolean('CLASS parameters','include_magnification')
        if config.has_option('CLASS parameters','include_gr_vel') :
            self.include_gr_vel=config.getboolean('CLASS parameters','include_gr_vel')
        if config.has_option('CLASS parameters','include_gr_pot') :
            self.include_gr_pot=config.getboolean('CLASS parameters','include_gr_pot')
        if config.has_option('CLASS parameters','exec_path') :
            self.exec_path=config.get('CLASS parameters','exec_path')
        if config.has_option('CLASS parameters','use_nonlinear') :
            self.use_nonlinear=config.getboolean('CLASS parameters','use_nonlinear')
        if config.has_option('CLASS parameters','f_sky') :
            self.fsky=config.getfloat('CLASS parameters','f_sky')

        #Output parameters
        if config.has_option('Output parameters','output_dir') :
            self.output_dir=config.get('Output parameters','output_dir')
        if config.has_option('Output parameters','output_spectra') :
            self.output_spectra=config.get('Output parameters','output_spectra')
        if config.has_option('Output parameters','output_fisher') :
            self.output_fisher=config.get('Output parameters','output_fisher')
        #Form output prefix and create output directory
        self.output_path=self.output_dir+"/"+self.output_spectra
        os.system("mkdir -p "+self.output_dir+"/"+self.output_fisher)

        #BAO
        if config.has_section("BAO 1") :
            self.include_BAO=True
        else :
            self.include_BAO=False
        self.n_bao=0
        while config.has_section("BAO %d"%(self.n_bao+1)) :
            self.n_bao+=1
            sec_title="BAO %d"%(self.n_bao)
            if config.has_option(sec_title,'fname_da') :
                self.include_DA=True
                zn,eda=np.loadtxt(config.get(sec_title,'fname_da'),unpack=True)
                self.z_nodes_DA.extend(zn)
                self.e_nodes_DA.extend(eda)
            if config.has_option(sec_title,'fname_hh') :
                self.include_HH=True
                zn,ehh=np.loadtxt(config.get(sec_title,'fname_hh'),unpack=True)
                self.z_nodes_HH.extend(zn)
                self.e_nodes_HH.extend(ehh)
            if config.has_option(sec_title,'fname_dv') :
                self.include_DV=True
                zn,edv=np.loadtxt(config.get(sec_title,'fname_dv'),unpack=True)
                self.z_nodes_DV.extend(zn)
                self.e_nodes_DV.extend(edv)
        self.z_nodes_DA=np.array(self.z_nodes_DA)
        self.z_nodes_HH=np.array(self.z_nodes_HH)
        self.z_nodes_DV=np.array(self.z_nodes_DV)
        self.e_nodes_DA=np.array(self.e_nodes_DA)
        self.e_nodes_HH=np.array(self.e_nodes_HH)
        self.e_nodes_DV=np.array(self.e_nodes_DV)

        #Tracers
        if not config.has_section("Tracer 1") :
            print "The param file contains no Tracers"
#            sys.exit("The param file must contain at least one tracer starting with Tracer 1")
        n_tracers=0
        while config.has_section("Tracer %d"%(n_tracers+1)) :
            n_tracers+=1
            sec_title="Tracer %d"%n_tracers

            #Generic
            consider_tracer=True; lmax=self.lmax; lmin=0;
            name=config.get(sec_title,'tracer_name')
            tr_type=config.get(sec_title,'tracer_type')
            if config.has_option(sec_title,'use_tracer') :
                consider_tracer=config.getboolean(sec_title,'use_tracer')
            if config.has_option(sec_title,'lmin') :
                lmin=max(lmin,config.getint(sec_title,'lmin'))
            if config.has_option(sec_title,'lmax') :
                lmax=min(lmax,config.getint(sec_title,'lmax'))

            #CMB primary or lensing
            has_t=False; has_p=False; sigma_t=None; sigma_p=None; beam_amin=None; l_transition=None;
            if config.has_option(sec_title,'has_t') :
                has_t=config.getboolean(sec_title,'has_t')
            if config.has_option(sec_title,'has_p') :
                has_p=config.getboolean(sec_title,'has_p')
            if config.has_option(sec_title,'beam_amin') :
                beamstr=config.get(sec_title,'beam_amin')
                beam_amin=np.array(map(float,beamstr.split()))
            if config.has_option(sec_title,'sigma_t') :
                sigmastr=config.get(sec_title,'sigma_t')
                sigma_t=np.array(map(float,sigmastr.split()))
            if config.has_option(sec_title,'sigma_p') :
                sigmastr=config.get(sec_title,'sigma_p')
                sigma_p=np.array(map(float,sigmastr.split()))
            if config.has_option(sec_title,'l_transition') :
                l_transition=config.getint(sec_title,'l_transition')

            #Galaxy clustering, lensing or intensity mapping
            bins_file=None; nz_file=None;
            if config.has_option(sec_title,'bins_file') :
                bins_file=config.get(sec_title,'bins_file')
            if config.has_option(sec_title,'nz_file') :
                nz_file=config.get(sec_title,'nz_file')

            #Galaxy clustering or intensity mapping
            bias_file=None; sbias_file=None; ebias_file=None;
            if config.has_option(sec_title,'bias_file') :
                bias_file=config.get(sec_title,'bias_file')
            if config.has_option(sec_title,'sbias_file') :
                sbias_file=config.get(sec_title,'sbias_file')
            if config.has_option(sec_title,'ebias_file') :
                ebias_file=config.get(sec_title,'ebias_file')

            #Galaxy clustering
            is_photometric=True
            if config.has_option(sec_title,'is_photometric') :
                is_photometric=config.getboolean(sec_title,'is_photometric')

            #Galaxy lensing
            sigma_gamma=None; abias_file=None; rfrac_file=None;
            if config.has_option(sec_title,'sigma_gamma') :
                sigma_gamma=config.getfloat(sec_title,'sigma_gamma')
            if config.has_option(sec_title,'abias_file') :
                abias_file=config.get(sec_title,'abias_file')
            if config.has_option(sec_title,'rfrac_file') :
                rfrac_file=config.get(sec_title,'rfrac_file')

            #Intensity mapping
            tz_file=None; dish_size=None; t_inst=None; t_total=None; n_dish=None;
            area_efficiency=None; fsky_im=None; im_type=None; base_file=None;
            include_fg=False; fit_fg=False;
            a_fg=None; alp_fg=None; bet_fg=None; xi_fg=None; nux_fg=None; lx_fg=None; 
            if config.has_option(sec_title,'tz_file') :
                tz_file=config.get(sec_title,'tz_file')
            if config.has_option(sec_title,'dish_size') :
                dish_size=config.getfloat(sec_title,'dish_size')
            if config.has_option(sec_title,'t_inst') :
                t_inst=config.getfloat(sec_title,'t_inst')
            if config.has_option(sec_title,'t_total') :
                t_total=config.getfloat(sec_title,'t_total')
            if config.has_option(sec_title,'n_dish') :
                n_dish=config.getint(sec_title,'n_dish')
            if config.has_option(sec_title,'area_efficiency') :
                area_efficiency=config.getfloat(sec_title,'area_efficiency')
            if config.has_option(sec_title,'fsky_im') :
                fsky_im=config.getfloat(sec_title,'fsky_im')
            if config.has_option(sec_title,'instrument_type') :
                im_type=config.get(sec_title,'instrument_type')
            if config.has_option(sec_title,'base_file') :
                base_file=config.get(sec_title,'base_file')
            if config.has_option(sec_title,'include_foregrounds') :
                include_fg=config.getboolean(sec_title,'include_foregrounds')
                if include_fg :
                    if config.has_option(sec_title,'fit_foregrounds') :
                        fit_fg=config.getboolean(sec_title,'fit_foregrounds')
                    if config.has_option(sec_title,'A_fg') :
                        a_fg=config.getfloat(sec_title,'A_fg')
                    if config.has_option(sec_title,'alpha_fg') :
                        alp_fg=config.getfloat(sec_title,'alpha_fg')
                    if config.has_option(sec_title,'beta_fg') :
                        bet_fg=config.getfloat(sec_title,'beta_fg')
                    if config.has_option(sec_title,'xi_fg') :
                        xi_fg=config.getfloat(sec_title,'xi_fg')
                    if config.has_option(sec_title,'nux_fg') :
                        nux_fg=config.getfloat(sec_title,'nux_fg')
                    if config.has_option(sec_title,'lx_fg') :
                        lx_fg=config.getfloat(sec_title,'lx_fg')
                    if self.include_im_fg==False :
                        self.include_im_fg=include_fg
                        self.fit_im_fg=fit_fg*include_fg
                        if self.fit_im_fg :
                            self.params_all.append(fsh.ParamFisher(a_fg,0.1*a_fg,"im_fg_a_fg",
                                                                   "$A_{\\rm FG}$",True,True,0))
                            self.params_all.append(fsh.ParamFisher(alp_fg,0.1*alp_fg,"im_fg_alp_fg",
                                                                   "$\\alpha_{\\rm FG}$",True,True,0))
                            self.params_all.append(fsh.ParamFisher(bet_fg,0.1*bet_fg,"im_fg_bet_fg",
                                                                   "$\\beta_{\\rm FG}$",True,True,0))
                            self.params_all.append(fsh.ParamFisher(xi_fg,0.1*xi_fg,"im_fg_xi_fg",
                                                                   "$\\xi_{\\rm FG}$",True,True,0))

            if (tr_type=="gal_clustering") or (tr_type=="intensity_mapping") :
                self.has_gal_clustering=True
                lmax=min(lmax,self.lmax_lss)
            elif tr_type=="gal_shear" :
                self.has_gal_shear=True
                lmax=min(lmax,self.lmax_lss)
            elif tr_type=="cmb_lensing" :
                self.has_cmb_lensing=True
                lmax=min(lmax,self.lmax_cmb)
            elif tr_type=="cmb_primary" :
                if has_t :
                    self.has_cmb_t=True
                if has_p :
                    self.has_cmb_p=True
                lmax=min(lmax,self.lmax_cmb)
            else :
                strout="Wrong tracer type \""+tr_type+"\"\n"
                strout+="Allowed tracer types are: gal_clustering, intensity_mapping, "
                strout+="gal_shear, cmb_lensing, cmb_primary\n"
                sys.exit(strout)

            self.tracers.append(trc.Tracer(self,name,tr_type,bins_file,nz_file,
                                           is_photometric,bias_file,sbias_file,ebias_file,
                                           abias_file,rfrac_file,sigma_gamma,
                                           has_t,has_p,sigma_t,sigma_p,beam_amin,l_transition,
                                           tz_file,dish_size,t_inst,t_total,n_dish,
                                           area_efficiency,fsky_im,im_type,base_file,
                                           a_fg,alp_fg,bet_fg,xi_fg,nux_fg,lx_fg,
                                           n_tracers,consider_tracer,lmin,lmax))
        self.n_tracers=n_tracers

        if ((self.n_tracers==0) and (self.n_bao==0)) :
            sys.exit("You need at least one tracer or some BAO measurements")

        #Search for nuisance bias parameters
        list_nuisance=[]
        def add_nuisance(nu) :
            if len(nu.z_arr)>0 :
                n_nodes=len(nu.z_arr)
                for i in np.arange(n_nodes) :
                    if nu.i_marg[i]>0 :
                        nui_name=nu.name+"node%d"%i
                        new_nuisance=True
                        for nm in list_nuisance :
                            if nui_name==nm :
                                new_nuisance=False
                        if new_nuisance :
                            list_nuisance.append(nui_name)
                            self.params_all.append(fsh.ParamFisher(nu.f_arr[i],nu.df_arr[i],
                                                                   nui_name,nui_name,True,False,0))
        for tr in self.tracers :
            if tr.consider_tracer :
                add_nuisance(tr.nuisance_bias)
                add_nuisance(tr.nuisance_sbias)
                add_nuisance(tr.nuisance_ebias)
                add_nuisance(tr.nuisance_abias)
                add_nuisance(tr.nuisance_rfrac)
                add_nuisance(tr.nuisance_sphz)
                add_nuisance(tr.nuisance_bphz)
            
    def get_param_properties(self,parname) :
        for par in self.params_all :
            if par.name==parname :
                return par.val,par.dval,par.onesided

    def get_cls_all(self) :
        if self.n_tracers<=0 :
            return
        #Run CLASS
        allfound=True
        allfound*=ino.start_running(self,"none",0)
        for i in np.arange(self.npar_vary) :
            pname=self.params_fshr[i].name
            if pname.startswith("im_fg") :
                continue
            allfound*=ino.start_running(self,pname,1)
            allfound*=ino.start_running(self,pname,-1)
        if self.just_run_cls :
            return

        #Wait for CLASS to finish
        if not allfound :
            while(not allfound) :
                allfound=True
                allfound*=ino.cls_are_there(self,"none",0,False)
                for i in np.arange(self.npar_vary) :
                    pname=self.params_fshr[i].name
                    if pname.startswith("im_fg") :
                        continue
                    allfound*=ino.cls_are_there(self,pname,1,False)
                    allfound*=ino.cls_are_there(self,pname,-1,False)
                if not allfound :
                    time.sleep(20)
            time.sleep(5)

        print "Reading fiducial"
        self.cl_fid_arr[:,:,:]=(ino.get_cls(self,"none",0)).reshape((self.lmax+1)/NLB,NLB,self.nbins_total,self.nbins_total).mean(axis=1)

        for i in np.arange(self.npar_vary) :
            pname=self.params_fshr[i].name
            if pname.startswith("im_fg") :
                continue
            print "Deriv for "+pname
            clp=(ino.get_cls(self,pname, 1)).reshape((self.lmax+1)/NLB,NLB,self.nbins_total,self.nbins_total).mean(axis=1)
            clm=(ino.get_cls(self,pname,-1)).reshape((self.lmax+1)/NLB,NLB,self.nbins_total,self.nbins_total).mean(axis=1)
            if self.params_fshr[i].onesided==0 :
                self.dcl_arr[i,:,:,:]=(clp-clm)/(2*self.params_fshr[i].dval)
            else :
                sig=1
                if self.params_fshr[i].onesided<0 :
                    sig=-1
                cl0=self.cl_fid_arr[:,:,:]
                self.dcl_arr[i,:,:,:]=(-3*cl0+4*clm-clp)/(2*sig*self.params_fshr[i].dval)

        if self.include_im_fg :
            print "Adding IM foregrounds"
            nbt1=0
            for tr1 in self.tracers :
                if tr1.consider_tracer==False :
                    continue
                nbt2=0
                nb1=tr1.nbins
                print "   "+tr1.name
                for tr2 in self.tracers :
                    if tr2.consider_tracer==False :
                        continue
                    nb2=tr2.nbins
                    if tr1.tracer_type=='intensity_mapping' :
                        if tr2.tracer_type=='intensity_mapping' :
                            cls=(trc.get_foreground_cls(tr1,tr2,self.lmax,"none"))
                            self.cl_fid_arr[:,nbt1:nbt1+nb1,nbt2:nbt2+nb2]+=cls.reshape((self.lmax+1)/NLB,NLB,nb1,nb2).mean(axis=1)
                            if self.fit_im_fg :
                                for ip in np.arange(self.npar_vary) :
                                    pname=self.params_fshr[ip].name 
                                    if pname.startswith("im_fg") :
                                        dcls=trc.get_foreground_cls(tr1,tr2,self.lmax,pname)
                                        self.dcl_arr[ip,:,nbt1:nbt1+nb1,nbt2:nbt2+nb2]=dcls.reshape((self.lmax+1)/NLB,NLB,nb1,nb2).mean(axis=1)
                    nbt2+=nb2
                nbt1+=nb1

#        for i in np.arange(self.npar_vary) :
#            print "Deriv for "+self.params_fshr[i].name
#            toplot=[]
#            larr=np.arange(len(self.dcl_arr[i,:,0,0]))
#            toplot.append(larr)
#            cols=['r','g','b']
#            for ib in np.arange(self.nbins_total) :
#                arr_toplot=self.dcl_arr[i,:,ib,ib]/self.cl_fid_arr[:,ib,ib]*np.sqrt(larr+0.5)
#                plt.plot(larr,arr_toplot,cols[ib%3]+'-')
#                toplot.append(arr_toplot)
#            toplot=np.array(toplot)
#            np.savetxt("dcl.txt",np.transpose(toplot))
#            plt.gca().set_xscale('log')
#            plt.show()

    def get_cls_noise(self) :
        if self.n_tracers<=0 :
            return
        nbt1=0
        for i1 in np.arange(self.n_tracers) :
            if self.tracers[i1].consider_tracer==False :
                continue
            nbt2=0
            nb1=self.tracers[i1].nbins
            print "   "+self.tracers[i1].name
            for i2 in np.arange(self.n_tracers) :
                if self.tracers[i2].consider_tracer==False :
                    continue
                nb2=self.tracers[i2].nbins
                cls=trc.get_cross_noise(self.tracers[i1],self.tracers[i2],self.lmax)
                self.cl_noise_arr[:,nbt1:nbt1+nb1,nbt2:nbt2+nb2]=cls.reshape((self.lmax+1)/NLB,NLB,nb1,nb2).mean(axis=1)
                nbt2+=nb2
            nbt1+=nb1

    def join_fishers(self) :
        self.fshr=self.fshr_cls+self.fshr_bao
        np.savetxt(self.output_dir+"/"+self.output_fisher+"/fisher_raw.txt",self.fshr)

    def get_fisher_bao(self) :
        """ Compute Fisher matrix from numerical derivatives """
        if self.n_bao<=0 :
            return
        
        fname_save=self.output_dir+"/"+self.output_fisher+"/fisher_raw_bao.txt"
        if os.path.isfile(fname_save) :
            self.fshr_bao[:,:]=np.load(fname_save)
        else :
            
            allfound=True
            allfound*=ino.start_running(self,"none",0)
            for i in np.arange(self.npar_vary) :
                allfound*=ino.start_running(self,self.params_fshr[i].name,1)
                allfound*=ino.start_running(self,self.params_fshr[i].name,-1)

            if not allfound :
                while(not allfound) :
                    allfound=True
                    allfound*=ino.bao_is_there(self,"none",0,False)
                    for i in np.arange(self.npar_vary) :
                        allfound*=ino.bao_is_there(self,self.params_fshr[i].name,1,False)
                        allfound*=ino.bao_is_there(self,self.params_fshr[i].name,-1,False)

                    if not allfound :
                        time.sleep(10)
                time.sleep(5)

            da_fid_nodes,hh_fid_nodes,dv_fid_nodes=ino.get_bao(self,"none",0)
            dda_nodes=np.zeros([self.npar_vary,len(self.z_nodes_DA)])
            dhh_nodes=np.zeros([self.npar_vary,len(self.z_nodes_HH)])
            ddv_nodes=np.zeros([self.npar_vary,len(self.z_nodes_DV)])
            for i in np.arange(self.npar_vary) :
                dap,hhp,dvp=ino.get_bao(self,self.params_fshr[i].name, 1)
                dam,hhm,dvm=ino.get_bao(self,self.params_fshr[i].name,-1)
                if self.params_fshr[i].onesided==0 :
                    dda_nodes[i,:]=(dap-dam)/(2*self.params_fshr[i].dval)
                    dhh_nodes[i,:]=(hhp-hhm)/(2*self.params_fshr[i].dval)
                    ddv_nodes[i,:]=(dvp-dvm)/(2*self.params_fshr[i].dval)
                else :
                    sig=1
                    if self.params_fshr[i].onesided<0 :
                        sig=-1
                    dda_nodes[i,:]=(-3*da_fid_nodes+4*dam-dap)/(2*sig*self.params_fshr[i].dval)
                    dhh_nodes[i,:]=(-3*hh_fid_nodes+4*hhm-hhp)/(2*sig*self.params_fshr[i].dval)
                    ddv_nodes[i,:]=(-3*dv_fid_nodes+4*dvm-dvp)/(2*sig*self.params_fshr[i].dval)
                dda_nodes[i,:]/=da_fid_nodes
                dhh_nodes[i,:]/=hh_fid_nodes
                ddv_nodes[i,:]/=dv_fid_nodes
            self.fshr_bao+=np.sum(dda_nodes[:,None,:]*dda_nodes[None,:,:]/self.e_nodes_DA**2,axis=2)
            self.fshr_bao+=np.sum(dhh_nodes[:,None,:]*dhh_nodes[None,:,:]/self.e_nodes_HH**2,axis=2)
            self.fshr_bao+=np.sum(ddv_nodes[:,None,:]*ddv_nodes[None,:,:]/self.e_nodes_DV**2,axis=2)
            np.savetxt(fname_save,self.fshr_bao)

    def get_fisher_cls(self) :
        """ Compute Fisher matrix from numerical derivatives """

#        pcs=csm.PcsPar()
#        pcs.set_verbosity(1)
#        pcs.background_set(0.315,0.685,0.049,-1.,0.,0.67,2.725)
        if self.n_tracers<=0 :
            return

        fname_save=self.output_dir+"/"+self.output_fisher+"/fisher_raw_l"
        if os.path.isfile(fname_save+".npy") :
            self.fshr_l=np.load(fname_save+".npy")
        else :
            icl_arr=np.zeros_like(self.cl_fid_arr)

            lmax_arr=np.zeros(self.nbins_total)
            lmin_arr=np.zeros(self.nbins_total)
            nbt=0
            for tr in self.tracers :
                if tr.consider_tracer :
                    zarr=None
                    if ((tr.tracer_type=='gal_clustering') or
                        (tr.tracer_type=='intensity_mapping') or
                        (tr.tracer_type=='gal_shear')) :
                        data=np.loadtxt(tr.bins_file,unpack=True)
                        zarr=(data[0]+data[1])/2
                    for ib in np.arange(tr.nbins)  :
                        if zarr!=None :
#                            lmn=int(5*np.pi*pcs.hubble(1./(1+zarr[ib]))*pcs.radial_comoving_distance(1./(1+zarr[ib])))
#                            print zarr[ib], lmn, tr.lmax_bins[ib]
                            lmn=tr.lmin
                        else :
                            lmn=tr.lmin
                        lmax_arr[nbt]=min(tr.lmax_bins[ib],tr.lmax)
                        lmin_arr[nbt]=lmn
                        nbt+=1

            for lb in np.arange((self.lmax+1)/NLB) :
                for ib in np.arange(NLB) :
                    l=lb*NLB+ib
                    if l==0 :
                        continue
                    ell=float(l)
                    indices=np.where((lmin_arr<=l) & (lmax_arr>=l))[0]
                    cl_fid=self.cl_fid_arr[lb,indices,:][:,indices]
                    cl_noise=self.cl_noise_arr[lb,indices,:][:,indices]
                    icl=np.linalg.inv(cl_fid+cl_noise)
                    for i in np.arange(self.npar_vary) :
                        dcl1=self.dcl_arr[i,lb,indices,:][:,indices]
                        for j in np.arange(self.npar_vary-i)+i :
                            dcl2=self.dcl_arr[j,lb,indices,:][:,indices]
                            self.fshr_l[i,j,l]=self.fsky*(ell+0.5)*np.trace(np.dot(dcl1,
                                                                                   np.dot(icl,
                                                                                          np.dot(dcl2,
                                                                                                 icl))))
                            if i!=j :
                                self.fshr_l[j,i,l]=self.fshr_l[i,j,l]

            np.save(self.output_dir+"/"+self.output_fisher+"/fisher_raw_l",self.fshr_l)

        self.fshr_cls[:,:]=np.sum(self.fshr_l,axis=2)
        np.savetxt(self.output_dir+"/"+self.output_fisher+"/fisher_raw_cls.txt",self.fshr_cls)

    def plot_fisher(self) :
        covar=np.linalg.inv(self.fshr)
        fishfile=open(self.output_dir+"/"+self.output_fisher+"/fisher.out","w")
        print "Fisher forecast :"
        for i in np.arange(len(self.params_fshr)) :
            name=self.params_fshr[i].name
            val=self.params_fshr[i].val
            sigma_m=np.sqrt(covar[i,i])
            sigma_f=1./np.sqrt(self.fshr[i,i])
            print " - "+name+" = %.4lE"%val+" +- %.4lE(m)"%sigma_m+" +- %.4lE(f)"%sigma_f
            fishfile.write(" - "+name+" = %.4lE"%val+" +- %.4lE(m)"%sigma_m+" +- %.4lE(f)"%sigma_f+"\n")
        ##Output full Fisher matrix--------------------------------
        #for i in np.arange(len(self.params_fshr)) :
        #    for j in np.arange(len(self.params_fshr)) :
        #        fishstr="%41E"%self.fshr[i,j]
        #        if(j==0):
        #            fishmatstr=fishstr
        #        else:
        #            fishmatstr=fishmatstr+" "+fishstr
        #    fishfile.write(fishmatstr+"\n")
        ##---------------------------------------------------------
        #fishfile.close()
#        fsh.plot_fisher_all(self.params_fshr,[self.fshr],["none"],[2],
#                            ["solid"],["black"],["whatever"],3.,
#                            self.output_dir+"/"+self.output_fisher+"/fisher"+self.plot_ext)
        ##Parameter pairs, all combinations (Elisa)---------------------
        #index_plot=np.where(np.array([p.do_plot for p in self.params_fshr]))
        #n_params=len(index_plot[0])
        #param_plot=self.params_fshr[index_plot]
        #for ip in range(0,len(param_plot)) :
        #    for jp in range(ip,len(param_plot)) :
        #        plt.figure(figsize=(6,5))
        #        plt.gcf().subplots_adjust(left=0.15)
        #        ax=plt.gca()
        #        fsh.plot_fisher_two(self.params_fshr,param_plot[ip].name,param_plot[jp].name,
        #                            [self.fshr],ax,["none"],[2],["solid"],["black"],3.)
        #        plt.savefig(self.output_dir+"/"+self.output_fisher+"/fisher_"+param_plot[ip].name+"_"+
        #                    param_plot[jp].name+self.plot_ext)
        #        plt.clf()
        #        plt.close()
        ##-------------------------------------------------------------

    def plot_cls(self) :
        if self.n_tracers<=0 :
            return
        cols=['r','g','b','k','y','m','c']
        ncols=len(cols)
        ibin_tot=0
        larr=np.arange((self.lmax+1)/NLB)*NLB+0.5*(NLB-1)
        for tr in self.tracers :
            if tr.consider_tracer==False :
                continue
            plt.figure()
            plt.title(tr.name)
            for ibin in np.arange(tr.nbins) :
                print ibin,tr.lmax_bins[ibin]
                indices=np.where((larr>=tr.lmin) & (larr<=min(tr.lmax,tr.lmax_bins[ibin])))[0]
                plt.plot(larr[indices],self.cl_fid_arr[indices,ibin_tot,ibin_tot],
                         cols[ibin%ncols]+'-',label="Bin %d"%ibin)
                plt.plot(larr[indices],self.cl_noise_arr[indices,ibin_tot,ibin_tot],
                         cols[ibin%ncols]+'--')
                ibin_tot+=1
            plt.gca().set_xscale('log')
            plt.gca().set_yscale('log')
            plt.legend(loc='lower left',labelspacing=0,ncol=1+tr.nbins/5,columnspacing=0.1)
            plt.ylabel("$C_\\ell$",fontsize=fs)
            plt.xlabel("$\\ell$",fontsize=fs)
            plt.savefig(self.output_dir+"/"+self.output_fisher+"/Cls_"+tr.name+self.plot_ext)
#            plt.show()
