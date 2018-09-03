import numpy as np

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
