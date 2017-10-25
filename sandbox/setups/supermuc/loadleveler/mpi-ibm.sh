#!/bin/bash

#@ job_type = parallel
#@ class = general
#@ node = 200 ###@ island_count= not needed for #### class general

#@ total_tasks= 3200 ## other example
##@ tasks_per_node = 16
#@ wall_clock_limit = 1:20:30
##                    1 h 20 min 30 secs
#@ job_name = mytest
#@ network.MPI  = sn_all,not_shared,us
#@ output       = job$(jobid).out
#@ error        = job$(jobid).err
#@ notification = always
#@ notify_user  = markert@ph1.uni-koeln.de
#@ queue
#@ island_count     = 1,2

. /etc/profile
. /etc/profile.d/modules.sh
poe ./myprog.exe
