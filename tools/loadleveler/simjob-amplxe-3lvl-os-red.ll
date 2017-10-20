#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name = mbsolve
#@ job_type = parallel
#@ class = fat
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/
#@ output = omp-os-red-$(jobid).out
#@ error = omp-os-red-$(jobid).out
#@ notification=always
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module load boost/1.61_icc
module load amplifier_xe

amplxe-cl --collect advanced-hotspots -r vtproj-fat/hotspots/0005/ -- build-openmp/mbsolve-tool/mbsolve-tool -m openmp-2lvl-pc-red -w matlab
