
!/bin/bash -l

# Batch script to run a serial R job under SGE.

# Request ten minutes of wallclock time (format hours:minutes:seconds).
#$ -l h_rt=24:00:0

# Request 1 gigabyte of RAM
#$ -l mem=16G

# Request 15 GB of tmp space
#$ -l tmpfs=15G

# Name of the job
#$ -N eDNAplus

# Set working directory to your scratch
#$ -wd /home/ucbtbev/Scratch/eDNAPlus_stan

# Move to $TMPDIR on the compute node
cd $TMPDIR

# ----------------------------
# 1. Load modules
# ----------------------------

# Load the R module and run your R program
module -f unload compilers mpi gcc-libs
module load r/recommended


