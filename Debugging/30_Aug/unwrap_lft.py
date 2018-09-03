import numpy as np

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
