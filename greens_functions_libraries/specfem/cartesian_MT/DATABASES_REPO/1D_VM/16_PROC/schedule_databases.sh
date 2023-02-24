#!/bin/bash
#SBATCH --job-name=specfem3D
#SBATCH --output=specfem3D.%j.output.txt
#SBATCH --time=04:00:00
#SBATCH --nodes=16
#SBATCH --ntasks-per-node=1
#SBATCH --mem-per-gpu=4096
#SBATCH --partition=snsm_itn19
#SBATCH --qos=openaccess
#SBATCH --gres=gpu:1

#export SLURM_OVERLAP=1
#ulimit unlimited

#below lines used for GCC4/OpenMPI 2.0.4 build of specfem3d_devel
module purge
module load apps/specfem3d/3.0.1i
#export CUDA_LIB=/apps/cuda/10.1/targets/x86_64-linux/lib
#export LD_LIBRARY_PATH=/usr/lib64/psm2-compat:$LD_LIBRARY_PATH
#export PSM2_MULTI_EP=1


##################################
# Explanation of above parameters
##
# job-name: What will the job be called when you use `squeue`
# output: What will the Slurm stdout file be called
# time: How much time ( hour:minutes:seconds ) will the job be allotted
# nodes: How many nodes ( max 20 when using GPUs ) will be allotted
# ntasks-per-node: How many threads will you need per node (max 2 when using GPUs)
# mem: How much memory (RAM) will you be allotted per node ( defined in MBs, default 512)
# partition: what partition of CIRCE do you want to run on? Need to specify one with GPUs if doing a GPU job
# qos=openaccess: required for GPU jobs
# gres=gpu:# : Used for GPU jobs to determine how many GPUs are assigned (max 2, should match ntasks?)
##################################

logfile=$SLURM_JOB_NAME"."$SLURM_JOB_ID.output.txt
echo "TIME BEGAN: " >> $logfile
date >> $logfile
SECONDS=0

# cleans output files
#mkdir -p OUTPUT_FILES
#rm -rf OUTPUT_FILES/*

# stores setup
cp DATA/meshfem3D_files/Mesh_Par_file OUTPUT_FILES/
cp DATA/Par_file OUTPUT_FILES/
cp DATA/CMTSOLUTION OUTPUT_FILES/
cp DATA/STATIONS OUTPUT_FILES/

echo "running example: `date`"
currentdir=`pwd`

# get the number of processors, ignoring comments in the Par_file
NPROC=`grep ^NPROC DATA/Par_file | grep -v -E '^[[:space:]]*#' | cut -d = -f 2`

BASEMPIDIR=`grep ^LOCAL_PATH DATA/Par_file | cut -d = -f 2 `


echo
echo "Running mesher"
echo

mpirun -np $NPROC xmeshfem3D

echo
echo "Generating databases"
echo

mpirun -np $NPROC xgenerate_databases

# checks exit code
if [[ $? -ne 0 ]]; then exit 1; fi

echo
echo "see results in directory: OUTPUT_FILES/"
echo
echo "done"
echo `date`


# To create a full mesh using combine_vol_data:
# cd bin
# xcombine_vol_data 0 0 vs ../OUTPUT_FILES/DATABASES_MPI/ ../OUTPUT_FILES/DATABASES_MPI/ 1
# cd ../OUTPUT_FILES/DATABASES_MPI/
# (check that mesh2vtu.pl is working)
# ../../../../../utils/Visualization/Paraview/mesh2vtu.pl -i vs.mesh -o vs.vtu