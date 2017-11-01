#!/bin/bash

if [ "$CLUSTER" == "tiger" ]; then
cat > submit_class.batch <<EOF
#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH -t 02:05:00
#SBATCH --mem=48gb
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu

# MPI compilers
export PATH="/tigress/smsharma/anaconda2/bin:$PATH"
source activate venv
module load openmpi/gcc/1.8.8/64

cd /tigress/smsharma/GoFish/

./class_mod $1
EOF
elif [ "$CLUSTER" == "della" ]; then
cat > submit_class.batch <<EOF
#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=12
#SBATCH -t 05:05:00
#SBATCH --mem=48gb
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu

# MPI compilers
export PATH="/tigress/smsharma/anaconda2/bin:$PATH"
source activate venv
module load openmpi/gcc/1.8.8/64

cd /tigress/smsharma/GoFish/

./class_mod $1
EOF
elif [ "$CLUSTER" == "feynman" ];then
cat > submit_class.batch <<EOF
#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=16
#SBATCH -t 05:05:00
#SBATCH --mem=64gb
##SBATCH --mail-type=begin
##SBATCH --mail-type=end
##SBATCH --mail-user=smsharma@princeton.edu
#SBATCH -p dept

# MPI compilers
export MKL_LP64_ILP64=ilp64
source /opt/intel/compilers_and_libraries_2016.2.181/linux/bin/compilervars.sh intel64
source /opt/intel/impi/5.0.3.048/bin64/mpivars.sh
module load openmpi-x86_64

source activate venv_py27
cd /group/hepheno/smsharma/GoFish

./class_mod $1
EOF
fi

sbatch submit_class.batch
rm submit_class.batch
