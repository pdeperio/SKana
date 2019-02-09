#!/bin/bash

# Warning: doesn't work with date changes

#RUNS=(080148 080152 080157 080159 080161 080163 080167)
RUNS=(080249 080254 080263 080265 080269 080275 080278 080282)

START_TIMES=(`for run in ${RUNS[@]}; do summary -run $run | grep "Start time" | cut -d' ' -f11; done`)
END_TIMES=(`for run in ${RUNS[@]}; do summary -run $run | grep "End time" | cut -d' ' -f13; done`)

irun=0
for run in ${RUNS[@]}
do

    #echo $run ${START_TIMES[$irun]} ${END_TIMES[$irun]}
    SEC1=`date +%s -d ${START_TIMES[$irun]}`
    SEC2=`date +%s -d ${END_TIMES[$irun]}`

    DIFFSEC=`expr ${SEC2} - ${SEC1}`

    echo ${DIFFSEC} 

    irun=$(($irun+1))

done
