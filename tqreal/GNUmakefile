#
#

ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk

LOCAL_INC	=

LOCAL_LIBS      =

LOCAL_LIBS	= #sktqconv_qb.o

OBJS   = ${LOCAL_LIBS} tqreal
#OBJS = tqreal

all: $(OBJS) 

clean: 
	$(RM) *.o *~ core fort.* *.root *.log $(OBJS)
