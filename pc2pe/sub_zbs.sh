#!/bin/bash

source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh

# pc2pe SK4+5 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
RUNS=(61889 61892 61893 61894 61895 80871 80873 80875 80877 80884 80885 80886)

WORKDIR=${PWD}

datadir="/disk0?/data?/sk?/tst"

BASEDIR=zbs

LOGDIR=${BASEDIR}/log
mkdir -p ${LOGDIR}

SCRIPTDIR=${BASEDIR}/script
mkdir -p ${SCRIPTDIR}

for run in ${RUNS[@]}; do

    run=`printf %06g $run`

    OUTDIR=${BASEDIR}/$run
    mkdir -p ${OUTDIR}

    #run=`echo $run | cut -c 1-4`

    rundir=$datadir/${run:0:4}/$run

    for rfmfile in `ls $rundir`; do

	nsub=`echo $rfmfile | cut -c 15-20`
	ofile=${rfmfile%.*}.zbs

        filepath=`ls -d $rundir/$rfmfile`

	SUBFILE=${SCRIPTDIR}/$run.$nsub.sh
	cat > ${SUBFILE}<<EOF
#!/bin/bash
source /usr/local/sklib_gcc4.8.5/skofl_19a/env.sh
hostname
ulimit -c 0

cd ${WORKDIR}

root2zbs ${filepath} ${OUTDIR}/$ofile 

EOF

	qsub -q calib -o ${LOGDIR}/$run.$nsub -e ${LOGDIR}/$run.$nsub ${SUBFILE}

    done
done
