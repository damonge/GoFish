import numpy as np
import sys as sys

path = sys.argv[1];
print(path);
if (len(sys.argv)==2):
    nq = 6;
else:
    nq=int(sys.argv[2]);

fisher_object = np.load(path);
fisher_matrix = fisher_object['fisher_tot'][:nq,:nq];

inverse = np.linalg.inv(fisher_matrix);

FoM = (np.linalg.det(inverse))**(-1./nq);

print FoM;
