import numpy as np
from itertools import izip, tee
def pairwise(iterable):
    "s -> (s0,s1), (s1,s2), (s2, s3), ..."
    a, b = tee(iterable)
    next(b, None)
    return izip(a, b)

for sample in ['gold','blue','red']:

	zi_arr = np.transpose(np.loadtxt("bins_" + sample + ".txt"))[0]
	zf_arr = np.transpose(np.loadtxt("bins_" + sample + ".txt"))[1]

	z_arr = list(zi_arr) + list([zf_arr[-1]])
	z_arr_mean = np.array([(x+y)/2. for x,y in pairwise(z_arr)])
	bbias_c=(1+z_arr_mean)**0.8
	marg_bias_c = np.ones(len(z_arr_mean))

	np.savetxt("bz_"+sample+"_all.txt",np.transpose([z_arr_mean,bbias_c,marg_bias_c]),fmt='%lE %lE %d')
