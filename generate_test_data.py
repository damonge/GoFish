import numpy as np
import os as os

test_path = 'ell_test/'

dndn_cases = ['A', 'B']
edn_cases = ['a', 'b', 'c', 'd']
nrows_edn = [0,1,2,3,4,5]

twonbins = 12

l = np.linspace(1, 100, num = 100)
dcl = np.sqrt(l)

for ii in range(0, len(dndn_cases)):
    for jj in range(0, len(edn_cases)):
        if edn_cases[jj] == 'd':
            for kk in range(1, len(nrows_edn)):
                folder = test_path + dndn_cases[ii] + '_' + edn_cases[jj] + str(nrows_edn[kk])
                os.system('mkdir ' + folder)
                os.system('rm ' + folder + '/*')
                for ll in range(0, twonbins):
                    for mm in range(0, twonbins):
                        fp = open(folder + '/' + str(ll) + '_' + str(mm) + '.ell', 'a+')
                        for nn in range(0, len(l)):
                            fp.write(str(l[nn]) + ',' + str(dcl[nn]) + '\n')
                        fp.close()

        else:
            for kk in range(0, 1):
                folder = test_path + dndn_cases[ii] + '_' + edn_cases[jj] + str(nrows_edn[kk])
                os.system('mkdir ' + folder)
                os.system('rm ' + folder + '/*')
                for ll in range(0, twonbins):
                    for mm in range(0, twonbins):
                        fp = open(folder + '/' + str(ll) + '_' + str(mm) + '.ell', 'a+')
                        for nn in range(0, len(l)):
                            fp.write(str(l[nn]) + ',' + str(dcl[nn]) + '\n')
                        fp.close()
