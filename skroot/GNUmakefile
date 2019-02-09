#
#

ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk

LOCAL_INC	= 

LOCAL_LIBS	= SKRootRead.o

OBJS   = SKRootRead.o sample_snld fit_hvscan_spe fit_hvscan_spe_hk hv_curve_hk individual_fit_hk

all: $(OBJS) 

clean: 
	$(RM) *.o *~ core fort.* $(OBJS)
