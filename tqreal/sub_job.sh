#!/bin/bash

# pc2pe SK3+4 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
RUNS=(34202 34204 61889 61892 61893 61894 61895)

WORKDIR=${PWD}

datadir="/disk01/data?/sk?/tst"

BASEDIR=output

LOGDIR=${BASEDIR}/log
mkdir -p ${LOGDIR}

ERRDIR=${BASEDIR}/err
mkdir -p ${ERRDIR}

SCRIPTDIR=${BASEDIR}/script
mkdir -p ${SCRIPTDIR}

for run in ${RUNS[@]}; do

    OUTDIR=${BASEDIR}/$run
    mkdir -p ${OUTDIR}

    #run=`echo $run | cut -c 1-4`

    rundir=$datadir/0${run:0:3}/0$run

    for rfmfile in `ls $rundir`; do

	nsub=`echo $rfmfile | cut -c 15-20`
	ofile=tqreal`echo $rfmfile | cut -c 4-25`

        filepath=`ls -d $rundir/$rfmfile`

	SUBFILE=${SCRIPTDIR}/$run.$nsub.sh
	cat > ${SUBFILE}<<EOF
#!/bin/bash
source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./tqreal ${OUTDIR}/$ofile ${filepath} 

EOF

	qsub -q calib -o ${LOGDIR}/$run.$nsub -e ${ERRDIR}/$run.$nsub ${SUBFILE}

    done
done
