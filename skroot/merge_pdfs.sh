#!/bin/bash 

DIR=hv_ana
PMTTYPE="" # Set this to _sk or _hk if input files were selected

for FILEROOT in Okay_fit$PMTTYPE Large_HV_Chi2$PMTTYPE Failed_HV_Fit$PMTTYPE; do
    (
	echo "Merging ${FILEROOT}"
	FILES=(`ls ${DIR}/${FILEROOT}_*.pdf`)
	pdfunite ${FILES[@]} ${DIR}/${FILEROOT}.pdf
    )&
done

wait
