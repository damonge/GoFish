import numpy as np
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

