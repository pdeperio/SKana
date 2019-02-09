#!/bin/csh -f

if ( "$1" == "" ) then
    echo usage: sub_job.csh run datadir
    exit
endif

mkdir -p ./out
mkdir -p ./err
mkdir -p ./script
mkdir -p output/$1

set datadir=$2
set run=`echo $1 | cut -c 1-4`

foreach rfmfile ( `ls $datadir/$1/` )

set nsub=`echo $rfmfile | cut -c 15-20`
set ofile=`echo $rfmfile | cut -c 4-25`

/bin/cp go_header.csh script/$1.$nsub.csh
echo 'source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.csh' >> script/$1.$nsub.csh
echo 'hostname' >> script/$1.$nsub.csh 
echo './sample_snld '$datadir'/'$1'/'$rfmfile output/$1/'snld'$ofile >> script/$1.$nsub.csh

qsub -q calib -o out/$1.$nsub -e err/$1.$nsub script/$1.$nsub.csh

end
