#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = fat
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/
#@ output = pade-3lvl-$(jobid).out
#@ error = pade-3lvl-$(jobid).out
#@ notification=always
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc

thread_s=40
thread_e=40
iteration=1
exp_method=pade
method=openmp-3lvl-os-red
device=song2005

# vary thread count
for threads in `seq $thread_s $thread_e`; do

# reproducibility
for it in `seq 1 $iterations`; do

out_dir=pade-3lvl-$LOADL_STEP_ID/$threads/$it/

mkdir -p $out_dir

echo "Thread count: " $threads

KMP_AFFINITY=granularity=fine,proclist=[`seq -s , 0 $(($threads - 1))`],explicit OMP_NUM_THREADS=$threads build-openmp-$exp_method/mbsolve-tool/mbsolve-tool -m $method -d $device -w matlab -o $out_dir/results.mat

done

done
