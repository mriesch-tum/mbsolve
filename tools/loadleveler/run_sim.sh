#!/bin/bash

# vary thread count
for threads in `seq 1 40`; do

# reproducibility
for it in `seq 1 $2`; do

out_dir=$1-$LOADL_STEP_ID/$threads/$it/

mkdir -p $out_dir

echo "Thread count: " $threads

KMP_AFFINITY=granularity=fine,proclist=[`seq -s , 1 $threads`],explicit OMP_NUM_THREADS=$threads mbsolve/build-openmp/mbsolve-tool/mbsolve-tool -m $1 -w MATLAB -o $out_dir/results.mat

done

done
