#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = fat
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/simulations/
#@ output = riesch2017b-rk-32k-$(schedd_host).$(jobid).$(stepid).out
#@ error = riesch2017b-rk-32k-$(schedd_host).$(jobid).$(stepid).out
#@ notification=always
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc

thread_s=1
thread_e=40
iterations=5

gridpoints=32768
endtime=600e-15
name=riesch2017b-rk-32k
method=openmp-2lvl-rk
device=ziolkowski1995

# vary thread count
for threads in `seq $thread_s $thread_e`; do

# reproducibility
for it in `seq 1 $iterations`; do

out_dir=$name-$LOADL_STEP_ID/$threads/$it/

mkdir -p $out_dir

echo "Thread count: " $threads

KMP_AFFINITY=granularity=fine,proclist=[`seq -s , 0 $(($threads - 1))`],explicit OMP_NUM_THREADS=$threads ../build-riesch2017b/mbsolve-tool/mbsolve-tool -m $method -d $device -w matlab -g $gridpoints -e $endtime -o $out_dir/$name.mat

done

done

# extract performance data
cat $name-$LOADL_STEP_ID.out | ../mbsolve/tools/loadleveler/generate_plot.sh > $name-$LOADL_STEP_ID/perf-$name.csv

mv $name-$LOADL_STEP_ID.out $name-$LOADL_STEP_ID/
