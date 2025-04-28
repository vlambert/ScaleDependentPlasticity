#!/bin/bash

### Number of cores
Nproc=8

### Make output directory
WDIR=/Users/valerelambert/Documents/C19_H8_s3_n02
ODIR=$WDIR/Output
rm -rf $WDIR
if [ ! -e $WDIR ]; then
    mkdir $WDIR
fi
if [ ! -e $ODIR ]; then
    mkdir $ODIR
fi
 
# Print out output directory
echo Output directory is $ODIR

# Run the mpi job
mpirun -np $Nproc ./EvolvingYieldBoundary > LOG <<EOF
#output directory
$ODIR
EOF
