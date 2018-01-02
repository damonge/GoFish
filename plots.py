import numpy as np
import matplotlib.pyplot as plt
import fishplot as fsh

fish_BDCM=np.load("/Users/darshkodwani/Documents/Darsh/Research/Covariance_matrix/GoFish/fisher_BDCM.npz")
fish_BICM=np.load("/Users/darshkodwani/Documents/Darsh/Research/Covariance_matrix/GoFish/fisher_BICM.npz")
fish_BDCM['fisher_tot']
fish_BDCM['names']
fish_BDCM['values']
fish_BDCM['labels']
fish_BICM['fisher_tot']
fish_BICM['names']
fish_BICM['values']
fish_BICM['labels']


pars=[]
for n,v,l in zip(fish_BDCM['names'],fish_BDCM['values'],fish_BDCM['labels']) :
    pars.append(fsh.ParamFisher(0,-1,n,l,True,True,0))
pars=np.array(pars)
#fsh.plot_fisher_all(pars,[dfy1['fisher_tot'],dfy10['fisher_tot']],['none','none'],[2,2],['solid','solid'],['r','g'],['Y1','Y10'],3.,'none')
fsh.plot_fisher_all(pars,[fish_BDCM['fisher_tot'],fish_BICM['fisher_tot']],
                    [{'col':'#F7B825','ls':'solid','lw':2,'alpha':1.0},{'col':'#0066CC','ls':'solid','lw':2,'alpha':1.0}],
                    ['BDCM','BICM'],2.,'fisher_covariance.png',do_1D=False)
plt.show()
#F7B825
#0066CC