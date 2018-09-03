import numpy as np
import os

# number of redshift bins
nbins = 6

# contains the default param files
param_path = '30_Aug_along_diagonal/'

# here the fisher results are stored
fisher_path = 'Fisher_results_30Aug_along_diagonal/'

# these files are the default params files, to which only the cutting scheme
# is added later
title_yes = param_path + 'param_default_yes.ini'
title_no = param_path + 'param_default_no.ini'

outputs_yes = 'outputs_may14_magnification_yes/Fisher'
outputs_no = 'outputs_may14_magnification_no/Fisher'

dndn_cases = ['G']

edn_cases = ['d']

nrows_edn_cases = ['0']#, '1', '2', '3', '4', '5']

# loop over all given cutting schemes
for ii in range(0, len(dndn_cases)):
    for jj in range(0, len(edn_cases)):
        for kk in range(0, len(nrows_edn_cases)):                                      

             # generate the params file for the current run
             title = dndn_cases[ii] + '_' + edn_cases[jj] + '_' + nrows_edn_cases[kk]    
             file_yes = param_path + 'param_yes_' + title + '.ini'                       
             file_no = param_path + 'param_no_' + title + '.ini'                         
             os.system('cp ' + title_yes + ' ' + file_yes)                               
             os.system('cp ' + title_no + ' ' + file_no)                                 

             # Now, the cutting scheme is appended to the params file
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

             # create output directories for the specific cutting scheme
             os.system('mkdir ' + fisher_path + title + '_yes')                          
             os.system('mkdir ' + fisher_path + title + '_no')                           

             # run main.py with the current param file as input
             os.system('python ../main.py ' + file_yes)                                     
             print 'Doing ' + file_yes                                                   
             os.system('python ../main.py ' + file_no)                                      
             print 'Doing ' + file_no                                                    

             # the fisher files are moved from the central fisher output directory to
             # the one specific to this cutting scheme
             os.system('mv ' + outputs_yes + '/* ' + fisher_path + title + '_yes')       
             os.system('mv ' + outputs_no + '/* ' + fisher_path + title + '_no')         
