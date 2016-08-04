#!/bin/bash 
#
# bash script to call quantifyTilesDownstream.py using SLURM scheduler on Sherlock
#
# Usage: quantify_tiles.sh image_dir seq_dir roff_dir fluor_dir gv_path
#
#all commands that start with SBATCH contain commands that are just used by SLURM for scheduling  
#################
#set a job name  
#SBATCH --job-name=quantify_tiles
#################  
#a file for job output, you can check job progress
#SBATCH --output=quantify_tiles.out
#################
# a file for errors from the job
#SBATCH --error=quantify_tiles.err
#################
#time you think you need; default is one hour
#in minutes in this case, hh:mm:ss
#SBATCH --time=30:00:00
#################
#quality of service; think of it as job priority
#SBATCH --partition=owners
#SBATCH --qos=normal
#################
#number of nodes you are requesting
#SBATCH --nodes=1
#################
#tasks to run per node; a "task" is usually mapped to a MPI processes.
# for local parallelism (OpenMP or threads), use "--ntasks-per-node=1 --cpus-per-task=16" instead
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=18
#################

# Load matlab
module load matlab/R2012b


# Define paths
image_dir=$1
seq_dir=$2
roff_dir=$3
fluor_dir=$4
gv_path=$5

# Define other script parameters
num_cores="18"
data_scaling="MiSeq_to_TIRFStation1"

# Define filter subsets to use for registration
reg_subset1="MSI_4BM"
reg_subset2="MSI_0SP"
reg_subset3="MSI_1SP"
reg_subset4="MSI_2SP"
reg_subset5="MSI_3SP"
reg_subset6="MSI_4SP"

# Set outputs appropriately 
script=$(basename $0)
script_name=${script%.*}
log_dir=$image_dir/$script_name"Logs"
log_file_suffix=".log"
err_file_suffix=".err"

mkdir -p $log_dir

# Quantification using SLURM scheduler
echo "Submitting jobs via SLURM..."

for d in $image_dir/*/
do

    d_base=$(basename $d)
    log_file=$log_dir/$d_base$log_file_suffix
    err_file=$log_dir/$d_base$err_file_suffix

    start_time=$SECONDS
    echo "Starting quantification for $image_dir at timepoint $d_base..."

    srun python2.7 /share/PI/wjg/lab/array_tools/CPscripts/quantifyTilesDownstream.py \
        -id $image_dir -ftd $seq_dir -rod $roff_dir -fd $fluor_dir -n $num_cores \
        -rs $reg_subset1 -rs $reg_subset2 -rs $reg_subset3 -rs $reg_subset4 -rs $reg_subset5 -rs $reg_subset6 \
        -sf $data_scaling -gv $gv_path 1> $log_file 2> $err_file

    duration=$(( SECONDS - start_time))
    echo "Done with quantification for $image_dir at timepoint $d_base. Duration: $duration" | tee -a $log_file

done	

