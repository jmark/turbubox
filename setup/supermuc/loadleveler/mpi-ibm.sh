#!/bin/bash

#@ job_type     = parallel
#@ class        = general
#@ network.MPI  = sn_all,not_shared,us
#@ output       = out.log
#@ error        = out.log
#@ notification = always
#@ notify_user  = markert@ph1.uni-koeln.de
#@ blocking     = unlimited
#@ restart      = no
#@ initialdir   = $(pwd)

#@ total_tasks  = 3200
#@ tasks_per_node = 16
#@ island_count     = 1,2
#@ node = 200

#@ wall_clock_limit = 1:20:30

#@ job_name = mytest

#@ queue

. /etc/profile
. /etc/profile.d/modules.sh
poe ./myprog.exe
