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

OM0=0.3
DOM=0.01
FB0=0.1156
DFB=0.001
HH0=0.67
DHH=0.01
W00=-1.
DW0=0.01
WA0=0.
DWA=0.01
NS0=0.96
DNS=0.01
AS0=2.46
DAS=0.01
FNL0=0.0
DFNL=0.1
EGR0=1.0
DEGR=0.1
LMAX=500
LMAX_CMB=4000

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
    npar_vary=5

    #Fisher matrix
    fshr_l=[]
    fshr=[]

    def __init__(self,fname) :
        #Initialize cosmological parameters
        self.params_all.append(fsh.ParamFisher(OM0,DOM,-1,"om","$\\Omega_M$",False,True))
        self.params_all.append(fsh.ParamFisher(FB0,DFB,-1,"fb","$f_b$",False,True))
        self.params_all.append(fsh.ParamFisher(HH0,DHH,-1,"hh","$h$",False,True))
        self.params_all.append(fsh.ParamFisher(W00,DW0,-1,"w0","$w_0$",False,True))
        self.params_all.append(fsh.ParamFisher(WA0,DWA,-1,"wa","$w_a$",False,True))
        self.params_all.append(fsh.ParamFisher(NS0,DNS,-1,"ns","$n_s$",False,True))
        self.params_all.append(fsh.ParamFisher(AS0,DAS,-1,"A_s","$A_s$",False,True))
        self.params_all.append(fsh.ParamFisher(FNL0,DFNL,-1,"fNL","$f_{\\rm NL}$",False,True))
        self.params_all.append(fsh.ParamFisher(EGR0,DEGR,-1,"eGR","$\\epsilon_{\\rm GR}$",False,True))
#        self.params_all=np.array(self.params_all)

        #Read parameter file
        self.read_param_file(fname)
        if (self.has_cmb_lensing==False) and (self.has_cmb_t==False) and (self.has_cmb_p==False) :
            self.lmax_cmb=0
        if (self.has_gal_shear==False) and (self.has_gal_clustering==False) :
            self.lmax_lss=0
        self.lmax=max(self.lmax_lss,self.lmax_cmb)

        #List parameters to vary over
        self.npar_vary=0
        for p in self.params_all :
            if p.isfree==True :
                self.params_fshr.append(fsh.ParamFisher(p.val,p.dval,p.prior,p.name,
                                                        p.label,p.isfree,p.do_plot))
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
            if tr.tracer_type=="gal_clustering" :
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
                strout+="Allowed tracer types are: gal_clustering, gal_shear, cmb_lensing, cmb_primary\n"
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
        self.cl_fid_arr=np.zeros([self.lmax+1,
                                  self.nbins_total,self.nbins_total])
        self.cl_noise_arr=np.zeros([self.lmax+1,
                                    self.nbins_total,self.nbins_total])
        self.dcl_arr=np.zeros([self.npar_vary,self.lmax+1,
                               self.nbins_total,self.nbins_total])

        #Allocate Fisher matrix
        self.fshr_l=np.zeros([self.npar_vary,self.npar_vary,self.lmax+1])
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
            print " - Terms included for galaxy shear : "+self.terms_gs
        print "  - %d free parameters :"%(self.npar_vary)
        for par in self.params_fshr :
            print "    * "+par.name+" : %.3lf"%(par.val)

    def read_param_file(self,fname) :
        config=ConfigParser.SafeConfigParser()
        config.read(fname)

        #Cosmological parameters
        for i in np.arange(len(self.params_all)) :
            p=self.params_all[i]
            if config.has_section(p.name) :
                x=p.val
                dx=p.dval
                isfree=p.isfree
                prior=p.prior
                if config.has_option(p.name,'x') :
                    x=config.getfloat(p.name,'x')
                if config.has_option(p.name,'dx') :
                    dx=config.getfloat(p.name,'dx')
                if config.has_option(p.name,'prior') :
                    prior=config.getfloat(p.name,'prior')
                if config.has_option(p.name,'is_free') :
                    isfree=config.getboolean(p.name,'is_free')
                self.params_all[i]=fsh.ParamFisher(x,dx,prior,p.name,p.label,isfree,isfree*True)

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
            self.include_magnification=config.getboolean('CLASS parameters',
                                                         'include_magnification')
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

        #Tracers
        if not config.has_section("Tracer 1") :
            sys.exit("The param file must contain at least one tracer starting with Tracer 1")
        n_tracers=0
        while config.has_section("Tracer %d"%(n_tracers+1)) :
            n_tracers+=1
            sec_title="Tracer %d"%n_tracers
            name=config.get(sec_title,'tracer_name')
            tr_type=config.get(sec_title,'tracer_type')

            bins_file=None; nz_file=None;
            bias_file=None; sbias_file=None; ebias_file=None;
            sigma_gamma=None; abias_file=None; rfrac_file=None;
            consider_tracer=True; lmax=self.lmax; lmin=0;
            has_t=False; has_p=False; sigma_t=None; beam_amin=None;
            if config.has_option(sec_title,'bins_file') :
                bins_file=config.get(sec_title,'bins_file')
            if config.has_option(sec_title,'nz_file') :
                nz_file=config.get(sec_title,'nz_file')
            if config.has_option(sec_title,'bias_file') :
                bias_file=config.get(sec_title,'bias_file')
            if config.has_option(sec_title,'sbias_file') :
                sbias_file=config.get(sec_title,'sbias_file')
            if config.has_option(sec_title,'ebias_file') :
                ebias_file=config.get(sec_title,'ebias_file')
            if config.has_option(sec_title,'sigma_gamma') :
                sigma_gamma=config.getfloat(sec_title,'sigma_gamma')
            if config.has_option(sec_title,'has_t') :
                has_t=config.getboolean(sec_title,'has_t')
            if config.has_option(sec_title,'has_p') :
                has_p=config.getboolean(sec_title,'has_p')
            if config.has_option(sec_title,'sigma_t') :
                sigma_t=config.getfloat(sec_title,'sigma_t')
            if config.has_option(sec_title,'beam_amin') :
                beam_amin=config.getfloat(sec_title,'beam_amin')
            if config.has_option(sec_title,'abias_file') :
                abias_file=config.get(sec_title,'abias_file')
            if config.has_option(sec_title,'rfrac_file') :
                rfrac_file=config.get(sec_title,'rfrac_file')
            if config.has_option(sec_title,'use_tracer') :
                consider_tracer=config.getboolean(sec_title,'use_tracer')
            if config.has_option(sec_title,'lmin') :
                lmin=config.getint(sec_title,'lmin')
            if config.has_option(sec_title,'lmax') :
                lmax=config.getint(sec_title,'lmax')

            if tr_type=="gal_clustering" :
                self.has_gal_clustering=True
                lmax=self.lmax_lss
            elif tr_type=="gal_shear" :
                self.has_gal_shear=True
                lmax=self.lmax_lss
            elif tr_type=="cmb_lensing" :
                self.has_cmb_lensing=True
                lmax=self.lmax_cmb
            elif tr_type=="cmb_primary" :
                if has_t :
                    self.has_cmb_t=True
                if has_p :
                    self.has_cmb_p=True
                lmax=self.lmax_cmb
            else :
                strout="Wrong tracer type \""+tr_type+"\"\n"
                strout+="Allowed tracer types are: gal_clustering, gal_shear, cmb_lensing, cmb_primary\n"
                sys.exit(strout)
            self.tracers.append(trc.Tracer(self,name,tr_type,bins_file,nz_file,
                                           bias_file,sbias_file,ebias_file,
                                           abias_file,rfrac_file,sigma_gamma,
                                           has_t,has_p,sigma_t,beam_amin,
                                           n_tracers,consider_tracer,lmin,lmax))
        self.n_tracers=n_tracers

        #Search for nuisance bias parameters
        def add_nuisance(nu) :
            if len(nu.z_arr)>0 :
                n_nodes=len(nu.z_arr)
                for i in np.arange(n_nodes) :
                    if nu.i_marg[i]>0 :
                        self.params_all.append(fsh.ParamFisher(nu.f_arr[i],nu.df_arr[i],nu.prior_arr[i],
                                                               nu.name+"n%d"%i,nu.name+"n%d"%i,True,False))
        for tr in self.tracers :
            if tr.consider_tracer :
                add_nuisance(tr.nuisance_bias)
                add_nuisance(tr.nuisance_sbias)
                add_nuisance(tr.nuisance_ebias)
                add_nuisance(tr.nuisance_abias)
                add_nuisance(tr.nuisance_rfrac)
                add_nuisance(tr.nuisance_phoz)

        #Behaviour parameters
        if config.has_option('Behaviour parameters','save_cl_files') :
            self.save_cl_files=config.getboolean('Behaviour parameters','save_cl_files')
        if config.has_option('Behaviour parameters','save_param_files') :
            self.save_param_files=config.getboolean('Behaviour parameters',
                                                    'save_param_files')
            
    def get_param_properties(self,parname) :
        for par in self.params_all :
            if par.name==parname :
                return par.val,par.dval

    def get_cls_all(self) :
        allfound=True
        allfound*=ino.start_running(self,"none",0)
        for i in np.arange(self.npar_vary) :
            allfound*=ino.start_running(self,self.params_fshr[i].name,1)
            allfound*=ino.start_running(self,self.params_fshr[i].name,-1)

        if not allfound :
            while(not allfound) :
                allfound=True
                allfound*=ino.files_are_there(self,"none",0,False)
                for i in np.arange(self.npar_vary) :
                    allfound*=ino.files_are_there(self,self.params_fshr[i].name,1,False)
                    allfound*=ino.files_are_there(self,self.params_fshr[i].name,-1,False)
                if not allfound :
                    time.sleep(20)
            time.sleep(5)

        self.cl_fid_arr[:,:,:]=ino.get_cls(self,"none",0)

        for i in np.arange(self.npar_vary) :
            clp=ino.get_cls(self,self.params_fshr[i].name,1)
            clm=ino.get_cls(self,self.params_fshr[i].name,-1)
            self.dcl_arr[i,:,:,:]=(clp-clm)/(2*self.params_fshr[i].dval)

    def get_cls_noise(self) :
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
                self.cl_noise_arr[:,nbt1:nbt1+nb1,nbt2:nbt2+nb2]=cls
                nbt2+=nb2
            nbt1+=nb1

    def get_fisher_l(self) :
        """ Compute Fisher matrix from numerical derivatives """

        icl_arr=np.zeros_like(self.cl_fid_arr)
        
        lmax_arr=np.zeros(self.nbins_total)
        lmin_arr=np.zeros(self.nbins_total)
        nbt=0
        for tr in self.tracers :
            if tr.consider_tracer :
                for lm in tr.lmax_bins :
                    lmax_arr[nbt]=min(lm,tr.lmax)
                    lmin_arr[nbt]=tr.lmin
                    nbt+=1

        for l in np.arange(self.lmax-1)+2 :
            indices=np.where((lmin_arr<=l) & (lmax_arr>=l))[0]
            cl_fid=self.cl_fid_arr[l,indices,:][:,indices]
            cl_noise=self.cl_noise_arr[l,indices,:][:,indices]
            icl=np.linalg.inv(cl_fid+cl_noise)
            for i in np.arange(self.npar_vary) :
                dcl1=self.dcl_arr[i,l,indices,:][:,indices]
                for j in np.arange(self.npar_vary-i)+i :
                    dcl2=self.dcl_arr[j,l,indices,:][:,indices]
                    self.fshr_l[i,j,l]=0.5*self.fsky*(2.*l+1)*np.trace(np.dot(dcl1,
                                                                              np.dot(icl,
                                                                                     np.dot(dcl2,
                                                                                            icl))))
                    if i!=j :
                        self.fshr_l[j,i,l]=self.fshr_l[i,j,l]
        self.fshr[:,:]=np.sum(self.fshr_l,axis=2)

        #Apply priors
        for i in np.arange(self.npar_vary) :
            dval=self.params_fshr[i].prior
            if dval>0 :
                self.fshr[i,i]+=1./dval**2

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
        #Output full Fisher matrix--------------------------------
        for i in np.arange(len(self.params_fshr)) :
            for j in np.arange(len(self.params_fshr)) :
                fishstr="%41E"%self.fshr[i,j]
                if(j==0):
                    self.fshrstr=fishstr
                else:
                    self.fshrstr=self.fshrstr+" "+fishstr
            fishfile.write(self.fshrstr+"\n")
        #---------------------------------------------------------
        fishfile.close()
        fsh.plot_fisher_all(self.params_fshr,[self.fshr],["none"],[2],
                            ["solid"],["black"],["whatever"],3.,self.output_dir+"/"+self.output_fisher+"/fisher"+self.plot_ext)
        #Parameter pairs, all combinations (Elisa)---------------------
        index_plot=np.where(np.array([p.do_plot for p in self.params_fshr]))
        n_params=len(index_plot[0])
        param_plot=self.params_fshr[index_plot]
        for ip in range(0,len(param_plot)) :
            for jp in range(ip,len(param_plot)) :
                plt.figure(figsize=(6,5))
                plt.gcf().subplots_adjust(left=0.15)
                ax=plt.gca()
                fsh.plot_fisher_two(self.params_fshr,param_plot[ip].name,param_plot[jp].name,
                                    [self.fshr],ax,["none"],[2],["solid"],["black"],3.)
                plt.savefig(self.output_dir+"/"+self.output_fisher+"/fisher_"+param_plot[ip].name+"_"+
                            param_plot[jp].name+self.plot_ext)
                plt.clf()
                plt.close()
        #-------------------------------------------------------------

    def plot_cls(self) :
        cols=['r','g','b','k','y','m','c']
        ncols=len(cols)
        ibin_tot=0
        larr=np.arange(self.lmax+1)
        for tr in self.tracers : 
            if tr.consider_tracer==False :
                continue
            plt.title(tr.name)
            for ibin in np.arange(tr.nbins) :
                plt.plot(larr[tr.lmin:min(tr.lmax,tr.lmax_bins[ibin])],
                         self.cl_fid_arr[tr.lmin:min(tr.lmax,tr.lmax_bins[ibin]),
                                         ibin_tot,ibin_tot],cols[ibin%ncols]+'-',label="Bin %d"%ibin)
                plt.plot(larr[tr.lmin:min(tr.lmax,tr.lmax_bins[ibin])],
                         self.cl_noise_arr[tr.lmin:min(tr.lmax,tr.lmax_bins[ibin]),
                                           ibin_tot,ibin_tot],cols[ibin%ncols]+'--')
                ibin_tot+=1
            plt.gca().set_xscale('log')
            plt.gca().set_yscale('log')
            plt.legend(loc='lower left',labelspacing=0,ncol=1+tr.nbins/5,columnspacing=0.1)
            plt.ylabel("$C_\\ell$",fontsize=fs)
            plt.xlabel("$\\ell$",fontsize=fs)
            plt.savefig(self.output_dir+"/"+self.output_fisher+"/Cls_"+tr.name+self.plot_ext)
            plt.show()

#        larr=np.arange(self.lmax+1)
#        for i in np.arange(self.nbins_total) :
#            for j in np.arange(self.nbins_total) :
#                plt.title("%d x "%i+"%d"%j)
#                print self.cl_fid_arr[:,i,j]
#                plt.plot(larr,np.fabs(self.cl_fid_arr[:,i,j]))
#                plt.gca().set_xscale('log')
#                plt.gca().set_yscale('log')
#                plt.show()
