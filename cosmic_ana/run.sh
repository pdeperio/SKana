RUNANA=$1

ANADIR=/Users/pdeperio/T2K/SK/Software/SKanadaRepo/atmpd/src/analysis/cosmic_ana_timevar

# Make sure to also change variable "recoSelect" in OfficialEventParser.h
# and filenames in global.cc
RECOS=(stmu13a fqv3r0)
#RECOS=(fqv2r1 fqv3r0)

#NRGS=(3000 999999)
NRGS=(999999)

for nrg in ${NRGS[@]}
do
    
    #for (( i=0; i<3; i++ ))
    for (( i=0; i<1; i++ ))
    do
	
	#rm -f h*_cosmics_reco[0-1]_nrg${i}_Elt${nrg}.root
	
	if [[ "$RUNANA" == 1 ]]; then

	    for (( iproc=0; iproc<8; iproc++ ))
	    do
		
		( 
		    ${ANADIR}/./Official $i $nrg $iproc
		)&
		
	    done	
	    
	    wait
	    
	    ireco=0
	    for reco in ${RECOS[@]}
	    do
		
		(	    
		    hadd -f hmc_reco${ireco}_nrg${i}_Elt${nrg}.root hmc_reco${ireco}_nrg${i}_Elt${nrg}_file[0-7].root
		    rm hmc_reco${ireco}_nrg${i}_Elt${nrg}_file[0-7].root
		)&
		(
		    hadd -f h_reco${ireco}_nrg${i}_Elt${nrg}.root h_reco${ireco}_nrg${i}_Elt${nrg}_file[0-7].root
		    rm h_reco${ireco}_nrg${i}_Elt${nrg}_file[0-7].root
		)&
		
		ireco=$(( $ireco+1 ))
	    done
	    
	    wait

	fi

	for (( inorm=0; inorm<2; inorm++ ))
	do
	    
	    for (( ilog=0; ilog<2; ilog++ ))
	    do
		
		#histname="vtx_prf"
		histname="r2_nearwall"

		#${ANADIR}/./run_compare 1 $ilog $inorm $i $nrg $histname  &

	    done
	done
	
	histname=""
	${ANADIR}/./run_compare 0 0 0 $i $nrg $histname &

    done
done


