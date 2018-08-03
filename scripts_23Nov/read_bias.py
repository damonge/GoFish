import sys
import numpy as np

filename = sys.argv[1]

fo = np.load(filename)

bias_vec = fo['fisher_bias'][:7]
fisher_matrix = fo['fisher_tot'][:7,:7]

inverse = np.linalg.inv(fisher_matrix)

bias = np.dot(inverse, bias_vec)

names = fo['names']

print fo.files

for ii in range(0, len(bias)):
    print names[ii] + '  ' + str(bias[ii]) + '  ' + str(np.sqrt(inverse[ii,ii])) + '\n'
