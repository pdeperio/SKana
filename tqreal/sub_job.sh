#!/bin/bash

# pc2pe SK4 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
RUNS=(61889 61892 61893 61894 61895 80871 80873 80875 80877 80884 80885 80886)

# SK5 datasets
#RUNS=(80871 80873 80875 80877 80884 80885 80886)

WORKDIR=${PWD}

datadir="/disk0?/data?/sk?/tst"

BASEDIR=output

LOGDIR=${BASEDIR}/log
mkdir -p ${LOGDIR}

ERRDIR=${BASEDIR}/err
mkdir -p ${ERRDIR}

SCRIPTDIR=${BASEDIR}/script
mkdir -p ${SCRIPTDIR}

for run in ${RUNS[@]}; do

    SK_GEOMETRY=4
    if [ $run -ge 80000 ]; then
        SK_GEOMETRY=5
    fi

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

EOF

        if [ ${SK_GEOMETRY} -eq 4 ]; then
            echo "./tqreal ${OUTDIR}/$ofile ${filepath} ${SK_GEOMETRY}" >> ${SUBFILE}
        else
            echo "./tqreal_ls_sk5 ${OUTDIR}/$ofile ${filepath}" >> ${SUBFILE}
        fi

	qsub -q calib -o ${LOGDIR}/$run.$nsub -e ${ERRDIR}/$run.$nsub ${SUBFILE}

    done
done
