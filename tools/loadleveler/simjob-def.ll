#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = fat
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/
#@ output = omp-def-$(jobid).out
#@ error = omp-def-$(jobid).out
#@ notification=always
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc

mbsolve/tools/loadleveler/run_sim.sh openmp-2lvl-pc 5 1 40
