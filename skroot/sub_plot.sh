#!/bin/bash -f
################################################
#
# Usage:
#         ./sub_plot.sh      # For all PMTs
#
#         ./sub_plot.sh _sk  # For SK PMTs only
#         ./sub_plot.sh _hk  # For HK PMTs only
#
################################################

PMT_TYPE=$1

if [[ "${PMT_TYPE}" != "" ]]; then
    PMT_TYPE_ARG="-t ${PMT_TYPE}"
fi

WORKDIR=${PWD}

PLOTS_PER_JOB=100
MAX_PLOTS=11200
NJOBS=$(($MAX_PLOTS / $PLOTS_PER_JOB))

BASEDIR=${WORKDIR}/hv_ana${PMT_TYPE}
mkdir -p ${BASEDIR}

SCRIPT_DIR=${BASEDIR}/script_plot
mkdir -p ${SCRIPT_DIR}

LOG_DIR=${BASEDIR}/log_plot
mkdir -p ${LOG_DIR}

for (( ijob=0; ijob<$NJOBS; ijob++ )); do

    SUBFILE=${SCRIPT_DIR}/plot_${ijob}.sh
    cat > ${SUBFILE}<<EOF
#!/bin/bash

source /usr/local/sklib_gcc4.8.5/skofl-trunk/env.sh
hostname

cd ${WORKDIR}

./individual_fit $PMT_TYPE_ARG -l $(($ijob * ${PLOTS_PER_JOB})) -u $(($(($ijob+1)) * ${PLOTS_PER_JOB}))

EOF

    qsub -q calib -o ${LOG_DIR}/$ijob.log -e ${LOG_DIR}/$ijob.log ${SUBFILE}

done
