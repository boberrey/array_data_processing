#!/bin/bash
#
# Quantify sets of images on sherlock
#
# Usage: quantifyAllsets.sh prefix seq_dir roff_dir fluor_dir gv_path script_dir
#
# Ben Ober-Reynolds, boberrey@stanford.edu, 20160804

# prefix is the common prefix shared by each set of images (will likely be 'set')
prefix=$1
seq_dir=$2
roff_dir=$3
fluor_dir=$4
gv_path=$5
script_dir=$6

for s in $prefix*/
do
    sbatch $script_dir/quantify_tiles.sbatch $s $seq_dir $roff_dir $fluor_dir $gv_path
done
