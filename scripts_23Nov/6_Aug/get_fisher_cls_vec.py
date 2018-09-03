    def get_fisher_cls_vec(self) :
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
