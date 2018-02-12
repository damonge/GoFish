import sys, os
import random

batch='''#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=1
#SBATCH -t 12:00:00
#SBATCH --mem=6GB
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p hepheno
export MKL_LP64_ILP64=ilp64
source /opt/intel/compilers_and_libraries_2016.2.181/linux/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.3.048/bin64/mpivars.sh
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:/group/hepheno/heptools/MultiNest/lib/
cd /group/hepheno/smsharma/GoFish/notebooks
source activate venv_py27
'''

for red in [1,0]:
    for vary_baryons in [1,0]:
        for m_bias in [1,0]:
            for photoz in [1,0]:
                for shear_lmax5000 in [1,0]:



					batchn = batch  + "\n"
					batchn += "python run.py --red " + str(red) + " --vary_baryons " + str(vary_baryons) + " --m_bias " + str(m_bias) + " --photoz " + str(photoz) + " --shear_lmax5000 " + str(shear_lmax5000)
					fname = "batch/" + str(red) + "_" + str(vary_baryons) + "_" + str(m_bias) + "_" + str(photoz) + "_" + str(shear_lmax5000) + ".batch" 
					f=open(fname, "w")
					f.write(batchn)
					f.close()
					os.system("chmod +x " + fname);
					os.system("sbatch " + fname);
