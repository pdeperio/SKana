#!/bin/bash -f

WORKDIR=${PWD}

RUNS=(80249 80254 80263 80265 80269 80275 80278 80282)
HV_SET=(0    +75   +50   +25    0    -25   -50   -75 )

CONNECTION_DIR=/disk02/usr6/pdeperio/precalib_2018

ihv=0
for hv in ${HV_SET[@]}; do

    if [ $hv == 0 ]; then
	CONNECTION_SET[$ihv]=${CONNECTION_DIR}/connection.super.sk-5.dat

    else
	CONNECTION_SET[$ihv]=${CONNECTION_DIR}/connection.super.sk-5${hv}.dat
    fi

    ihv=$(($ihv+1))
    
done

mkdir -p ./hv_ana

SCRIPT_DIR=${WORKDIR}/hv_ana_script
mkdir -p ${SCRIPT_DIR}

LOG_DIR=${WORKDIR}/hv_ana_log
mkdir -p ${LOG_DIR}

irun=0
for run in ${RUNS[@]}; do

    SUBFILE=${SCRIPT_DIR}/fit_${run}.sh
    cat > ${SUBFILE}<<EOF
#!/bin/bash

source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./fit_hvscan_spe -r ${run} -c "${CONNECTION_SET[$irun]}" -h ${HV_SET[$irun]}
EOF

    qsub -q calib -o ${LOG_DIR}/$run.log -e ${LOG_DIR}/$run.log ${SUBFILE}

    irun=$(($irun+1))
done
