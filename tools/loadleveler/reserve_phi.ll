#!/bin/bash
#@ wall_clock_limit = 24:00:00
#@ job_name =TEST
#@ job_type = parallel
#@ class = phi
#@island_count = 1
#@ node = 1
#@tasks_per_node = 1
#@ node_usage = not_shared
#@ initialdir = $(home)/
#@ output = test-$(jobid).out
#@ error = test-$(jobid).out
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

sleep 36000

