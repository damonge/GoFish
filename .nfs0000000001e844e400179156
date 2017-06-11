#!/bin/bash
cat > submit_class.batch <<EOF
#!/bin/bash
#SBATCH -N 1   # node count
#SBATCH --ntasks-per-node=14
#SBATCH -t 01:00:00
#SBATCH --mem=50gb
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
cd /group/hepheno/smsharma/GoFish-PR

./class_mod $1

EOF

sbatch submit_class.batch
rm submit_class.batch
