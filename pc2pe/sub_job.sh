#!/bin/bash

# pc2pe SK4 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
RUNS=(61889 61892 61893 61894 61895)

RUNS=(80871)

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

    for infile in `ls $rundir`; do

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

./llaser_qb_c ${filepath} ${OUTDIR}/$ofile

EOF

	qsub -q calib -o ${LOGDIR}/${basename}.log -e ${ERRDIR}/${basename}.err ${SUBFILE}

    done
done
