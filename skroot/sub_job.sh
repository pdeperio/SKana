#!/bin/bash -f

RUNS=(80249 80254 80263 80265 80269 80275 80278 80282)

datadir=/disk01/calib/sk5

BASEDIR=output

LOGDIR=${BASEDIR}/out
mkdir -p ${LOGDIR}

ERRDIR=${BASEDIR}/err
mkdir -p ${ERRDIR}

SCRIPTDIR=${BASEDIR}/script
mkdir -p ${SCRIPTDIR}

for run in ${RUNS[@]}; do

    OUTDIR=${BASEDIR}/$1
    mkdir -p ${OUTDIR}

    run=`echo $run | cut -c 1-4`

foreach rfmfile ( `ls $datadir/$1/` )

set nsub=`echo $rfmfile | cut -c 15-20`
set ofile=`echo $rfmfile | cut -c 4-25`

/bin/cp go_header.csh script/$1.$nsub.csh
echo 'source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.csh' >> script/$1.$nsub.csh
echo 'hostname' >> script/$1.$nsub.csh 
echo './sample_snld '$datadir'/'$1'/'$rfmfile output/$1/'snld'$ofile >> script/$1.$nsub.csh

qsub -q calib -o out/$1.$nsub -e err/$1.$nsub script/$1.$nsub.csh

done
