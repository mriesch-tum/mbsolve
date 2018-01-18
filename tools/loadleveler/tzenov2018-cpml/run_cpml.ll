#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = test
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/simulations/
#@ output = tzenov2018-cpml-$(schedd_host).$(jobid).$(stepid).out
#@ error = tzenov2018-cpml-$(schedd_host).$(jobid).$(stepid).out
#@ notification=always
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc

iterations=1
threads=32

gridpoints=32768
endtime=300e-15
name=tzenov2018-cpml
method=openmp-2lvl-os-red
device=tzenov2018-cpml

# reproducibility
for it in `seq 1 $iterations`; do

out_dir=$name-$LOADL_STEP_ID/$it/

mkdir -p $out_dir

echo "Thread count: " $threads

KMP_AFFINITY=granularity=fine,proclist=[`seq -s , 0 $(($threads - 1))`],explicit OMP_NUM_THREADS=$threads ../build-tzenov2018/mbsolve-tool/mbsolve-tool -m $method -d $device -w matlab -g $gridpoints -e $endtime -o $out_dir/$name.mat

done

# move output file to simulation folder
mv $name-$LOADL_STEP_ID.out $name-$LOADL_STEP_ID/
