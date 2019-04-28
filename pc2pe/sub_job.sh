#!/bin/bash

# pc2pe SK4+5 datasets
# https://docs.google.com/spreadsheets/d/1k3jO6WDPLnsYxZszKDaTvJiOfG688_lleVT4MxfVDJc/edit?usp=sharing
#RUNS=(61889 61892 61893 61894 61895 80871 80873 80875 80877 80884 80885 80886 81028 81030)  # Unstable
RUNS=(61889 61892 61893 61894 61895 80873 80875 80877 80885 80886 81028 81030)  # Mostly stable

# SK5
#RUNS=(80871 80873 80875 80877 80884 80885 80886)

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

    # Get start and end times
    itmp=0
    for when in Start End; do
      DATE_ARR=(`summary -run $run | grep "${when} time"`)
      TIME[itmp]=`date +%s -d"${DATE_ARR[4]} ${DATE_ARR[5]}, ${DATE_ARR[7]} ${DATE_ARR[6]}"`
      itmp=$(($itmp+1))
    done
    echo $run ${TIME[0]} ${TIME[1]}
    #continue
    
    # But overwrite with global time across runs within an SK period as follows
    TIME[0]=1224644086
    TIME[1]=1224685688

    # Koshio-san's SK5 files
    if [ $run -ge 80000 ]; then
        #datadir=/disk01/calib/sk5/
        #rundir=$datadir/0$run
        ConnectionFile=${SKOFL_ROOT}/const/connection.super.sk-5.dat
        TIME[0]=1553479750
        TIME[1]=1553601921

        if [ $run -ge 81000 ]; then
            TIME[0]=1555979550
            TIME[1]=1556031598
        fi
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

./llaser_qb_c ${filepath} ${OUTDIR}/$ofile ${run} ${ConnectionFile} ${TIME[0]} ${TIME[1]}

EOF

	qsub -q calib -o ${LOGDIR}/${basename}.log -e ${ERRDIR}/${basename}.err ${SUBFILE}

    done
done
