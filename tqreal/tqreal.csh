#!/bin/csh -f

if ( "$1" == "" ) then
    echo usage: tqreal.csh run
    exit
endif

mkdir -p ./out
mkdir -p ./err
mkdir -p ./script

#set outdir=/disk02/usr6/pdeperio/calib/$1
set outdir=/disk01/calib/sk5/$1

mkdir -p $outdir

set run=`echo $1 | cut -c 1-4`

foreach rfmfile ( `ls /disk02/data7/sk5/tst/$run/$1/` )

set nsub=`echo $rfmfile | cut -c 15-20`
set ofile=`echo $rfmfile | cut -c 4-25`
/bin/cp go_header.csh script/$1.$nsub.csh

echo 'source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.csh' >> script/$1.$nsub.csh
echo 'hostname' >> script/$1.$nsub.csh 
echo './tqreal '$outdir'/ofl'$ofile '/disk02/data7/sk5/tst/'$run'/'$1'/'$rfmfile >> script/$1.$nsub.csh
 
qsub -q calib -o out/$1.$nsub -e err/$1.$nsub script/$1.$nsub.csh

end
