#!/bin/bash 

DIR=hv_ana
PMTTYPE="" # Set this to _sk or _hk if input files were selected

for FILEROOT in HV_Curves$PMTTYPE SPE_HV$PMTTYPE; do
    (
	echo "Merging ${FILEROOT}"
	FILES=(`ls ${DIR}/${FILEROOT}_*.pdf`)
	pdfunite ${FILES[@]} ${DIR}/${FILEROOT}.pdf
    )&
done

hadd -f ${DIR}/hvscan_parameter.root ${DIR}/hvscan_parameter_*.root &

for FILEROOT in badcables badfittings badokfitting deadchannels targethv; do
    (
        echo "Merging ${FILEROOT}"
	cat ${DIR}/${FILEROOT}_*.txt | head -1 > ${DIR}/${FILEROOT}.txt
	awk FNR!=1 ${DIR}/${FILEROOT}_*.txt >> ${DIR}/${FILEROOT}.txt
    ) &
done

wait
