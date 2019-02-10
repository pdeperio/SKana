#!/bin/bash -f

RUNS=(080249 080254 080263 080265 080269 080275 080278 080282)

WORKDIR=${PWD}

datadir=/disk01/calib/sk5

BASEDIR=output

LOGDIR=${BASEDIR}/out
mkdir -p ${LOGDIR}

ERRDIR=${BASEDIR}/err
mkdir -p ${ERRDIR}

SCRIPTDIR=${BASEDIR}/script
mkdir -p ${SCRIPTDIR}

for run in ${RUNS[@]}; do

    OUTDIR=${BASEDIR}/$run
    mkdir -p ${OUTDIR}

    #run=`echo $run | cut -c 1-4`

    for rfmfile in `ls $datadir/$run/`; do
	
	nsub=`echo $rfmfile | cut -c 15-20`
	ofile=`echo $rfmfile | cut -c 4-25`

	SUBFILE=${SCRIPTDIR}/$run.$nsub.sh
	cat > ${SUBFILE}<<EOF
#!/bin/bash
source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./sample_snld ${datadir}/${run}/${rfmfile} ${OUTDIR}/snld${ofile}

EOF
	

	echo qsub -q calib -o ${LOGDIR}/$run.$nsub -e ${ERRDIR}/$run.$nsub ${SUBFILE}

    done
done
