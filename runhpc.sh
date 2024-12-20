#!/bin/sh

# Job name
#SBATCH --partition=cpuq
#SBATCH --account=cpuq
#SBATCH --job-name="H8n01" 

# Output and error files
#SBATCH -o out # STDOUT
#SBATCH -e err # STDERR

# Number of processor cores / tasks
#SBATCH --ntasks=100

# Wall time : maximum allowed run time
#SBATCH --time=4:00:00   

# change the working directory to current directory
echo Working directory is $SLURM_SUBMIT_DIR
cd $SLURM_SUBMIT_DIR

# Write out some information on the job
echo Running on host `hostname`
echo Time is `date`
echo $SLURM_JOB_NAME 
### Define number of processors
echo This job has allocated $SLURM_NPROCS cpus

### Make output directory
WDIR=/data/users/valamber/C19_H8_s3_n02
ODIR=$WDIR/Output
rm -rf $WDIR
if [ ! -e $WDIR ]; then
    mkdir $WDIR
fi
if [ ! -e $ODIR ]; then
    mkdir $ODIR
fi

# Tell me which nodes it is run on
echo " "
echo This jobs runs on the following processors:
echo $SLURM_JOB_NODELIST
echo " "
 
# Print out output directory
echo Output directory is $ODIR

# Run the mpi job
mpirun ./EvolvingYieldBoundary > LOG <<EOF
#output directory
$ODIR
EOF
