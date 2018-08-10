import numpy as np
import os

nbins = 6

param_path = 'scripts_23Nov/9_Aug/'
fisher_path = 'Fisher_results_9_Aug/'

title_yes = param_path + 'param_default_yes.ini'
title_no = param_path + 'param_default_no.ini'

outputs_yes = 'outputs_may14_magnification_yes/Fisher'
outputs_no = 'outputs_may14_magnification_no/Fisher'

dndn_cases = ['A', 'B']

edn_cases = ['a' , 'b', 'c', 'd']

nrows_edn_cases = ['0', '1', '2', '3', '4', '5']

for ii in range(0, len(dndn_cases)):
    for jj in range(0, len(edn_cases)):
        if edn_cases[jj] != 'd':
            title = dndn_cases[ii] + '_' + edn_cases[jj] + '_' + nrows_edn_cases[0]
            file_yes = param_path + 'param_yes_' + title + '.ini'
            file_no = param_path + 'param_no_' + title + '.ini'
            os.system('cp ' + title_yes + ' ' + file_yes)
            os.system('cp ' + title_no + ' ' + file_no)
            fp_yes = open(file_yes, 'a')
            fp_no = open(file_no, 'a')
            cutting_string = '\n[Cutting parameters]\n'
            cutting_string += 'nbins_lft= ' + str(nbins) + '\n'
            cutting_string += 'dndn_scheme_lft= ' + dndn_cases[ii] + '\n'
            cutting_string += 'edn_scheme_lft= ' + edn_cases[jj] + '\n'
            cutting_string += 'nrows_edn_lft= ' + nrows_edn_cases[0]
            fp_yes.write(cutting_string)
            fp_no.write(cutting_string)
            fp_yes.close()
            fp_no.close()
            os.system('mkdir ' + fisher_path + title + '_yes')
            os.system('mkdir ' + fisher_path + title + '_no')
            os.system('python main.py ' + file_yes)
            print 'Doing ' + file_yes
            os.system('python main.py ' + file_no)
            print 'Doing ' + file_no
            os.system('mv ' + outputs_yes + '/* ' + fisher_path + title + '_yes')
            os.system('mv ' + outputs_no + '/* ' + fisher_path + title + '_no')
        else:
            for kk in range(1, len(nrows_edn_cases)):
                title = dndn_cases[ii] + '_' + edn_cases[jj] + '_' + nrows_edn_cases[kk] 
                file_yes = param_path + 'param_yes_' + title + '.ini'
                file_no = param_path + 'param_no_' + title + '.ini'
                os.system('cp ' + title_yes + ' ' + file_yes)
                os.system('cp ' + title_no + ' ' + file_no)
                fp_yes = open(file_yes, 'a')
                fp_no = open(file_no, 'a')
                cutting_string = '\n[Cutting parameters]\n'
                cutting_string += 'nbins_lft= ' + str(nbins) + '\n'
                cutting_string += 'dndn_scheme_lft= ' + dndn_cases[ii] + '\n'
                cutting_string += 'edn_scheme_lft= ' + edn_cases[jj] + '\n'
                cutting_string += 'nrows_edn_lft= ' + nrows_edn_cases[kk]
                fp_yes.write(cutting_string)
                fp_no.write(cutting_string)
                fp_yes.close()
                fp_no.close()
                os.system('mkdir ' + fisher_path + title + '_yes')                    
                os.system('mkdir ' + fisher_path + title + '_no')
                os.system('python main.py ' + file_yes)
                print 'Doing ' + file_yes
                os.system('python main.py ' + file_no)
                print 'Doing ' + file_no
                os.system('mv ' + outputs_yes + '/* ' + fisher_path + title + '_yes')
                os.system('mv ' + outputs_no + '/* ' + fisher_path + title + '_no')

