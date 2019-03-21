#!/bin/bash -f
################################################
#
# Usage:
#         ./sub_hvfit.sh      # For all PMTs
#
#         ./sub_hvfit.sh _sk  # For sK PMTs only
#         ./sub_hvfit.sh _hk  # For HK PMTs only
#
################################################

PMT_TYPE=$1

if [[ "${PMT_TYPE}" != "" ]]; then
    PMT_TYPE_ARG="-t ${PMT_TYPE}"
fi

WORKDIR=${PWD}

PLOTS_PER_JOB=1000
MAX_PLOTS=12000
NJOBS=$(($MAX_PLOTS / $PLOTS_PER_JOB))

BASEDIR=${WORKDIR}/hv_ana${PMT_TYPE}
mkdir -p ${BASEDIR}

SCRIPT_DIR=${BASEDIR}/script_hvfit
mkdir -p ${SCRIPT_DIR}

LOG_DIR=${BASEDIR}/log_hvfit
mkdir -p ${LOG_DIR}

for (( ijob=0; ijob<$NJOBS; ijob++ )); do

    SUBFILE=${SCRIPT_DIR}/plot_${ijob}.sh
    cat > ${SUBFILE}<<EOF
#!/bin/bash

source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./hv_curve $PMT_TYPE_ARG -l $(($ijob * ${PLOTS_PER_JOB})) -u $(($(($ijob+1)) * ${PLOTS_PER_JOB}))

EOF

    qsub -q calib -o ${LOG_DIR}/$ijob.log -e ${LOG_DIR}/$ijob.log ${SUBFILE}

done
