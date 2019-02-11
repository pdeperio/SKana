#!/bin/bash 

DIR=hv_ana

for FILEROOT in Okay_fit_SK Large_HV_Chi2_SK Failed_HV_Fit_SK; do
    (
	echo "Merging ${FILEROOT}"
	FILES=(`ls ${DIR}/${FILEROOT}_*.pdf`)
	pdfunite ${FILES[@]} ${DIR}/${FILEROOT}.pdf
    )&
done

wait
