#!/bin/bash

# pc2pe SK4+5 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
RUNS=(61889 61892 61893 61894 61895 80871 80873 80875 80877 80884 80885 80886)

# SK5
RUNS=(80871 80873 80875 80877 80884 80885 80886)

WORKDIR=${PWD}

datadir="/disk02/usr6/pdeperio/SKana/tqreal/output"

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

    rundir=$datadir/$run
    ConnectionFile=${SKOFL_ROOT}/const/connection.super.sk-4.dat

    # Koshio-san's SK5 files
    if [ $run -ge 80000 ]; then
        datadir=/disk01/calib/sk5/
        rundir=$datadir/0$run
        ConnectionFile=${SKOFL_ROOT}/const/connection.super.sk-5.dat
    fi

    for infile in `ls -p $rundir | grep -v '/$'`; do

	basename=`echo $infile | cut -d'_' -f2`
        basename=${basename%.*}
	ofile=pc2pe_${basename}.root

        filepath=`ls -d $rundir/$infile`

	SUBFILE=${SCRIPTDIR}/${basename}.sh
	cat > ${SUBFILE}<<EOF
#!/bin/bash
source /usr/local/sklib_gcc4.8.5/skofl_19a/env.sh
hostname

cd ${WORKDIR}

./llaser_qb_c ${filepath} ${OUTDIR}/$ofile ${run} ${ConnectionFile}

EOF

	qsub -q calib -o ${LOGDIR}/${basename}.log -e ${ERRDIR}/${basename}.err ${SUBFILE}

    done
done
