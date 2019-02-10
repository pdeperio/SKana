#!/bin/bash -f
################################################
#
# Usage:
#         ./sub_ana.sh      # For SK PMTs only
#
#         ./sub_ana.sh _hk  # For HK PMTs only
#
################################################

PMT_TYPE=$1

WORKDIR=${PWD}

RUNS=(80249 80254 80263 80265 80269 80275 80278 80282)
HV_SET=(0    +75   +50   +25    0    -25   -50   -75 )

CONNECTION_DIR=/disk02/usr6/pdeperio/precalib_2018

ihv=0
for hv in ${HV_SET[@]}; do

    if [ $hv == 0 ]; then
	CONNECTION_SET[$ihv]=${CONNECTION_DIR}/connection.super.sk-5.dat
        ln -sf ${CONNECTION_SET}

    else
	CONNECTION_SET[$ihv]=${CONNECTION_DIR}/connection.super.sk-5${hv}.dat
    fi

    ihv=$(($ihv+1))
    
done

BASEDIR=${WORKDIR}/hv_ana${PMT_TYPE}
mkdir -p ${BASEDIR}

SCRIPT_DIR=${BASEDIR}/script
mkdir -p ${SCRIPT_DIR}

LOG_DIR=${BASEDIR}/log
mkdir -p ${LOG_DIR}

irun=0
for run in ${RUNS[@]}; do

    SUBFILE=${SCRIPT_DIR}/fit_${run}.sh
    cat > ${SUBFILE}<<EOF
#!/bin/bash

source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./fit_hvscan_spe -r ${run} -c "${CONNECTION_SET[$irun]}" -h ${HV_SET[$irun]} -t ${PMT_TYPE}

EOF

    qsub -q calib -o ${LOG_DIR}/$run.log -e ${LOG_DIR}/$run.log ${SUBFILE}

    irun=$(($irun+1))
done
