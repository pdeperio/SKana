#
# In case of building locally, set SKOFL_ROOT variable 
#      setenv SKOFL_ROOT   ... 
#  or directly set in this makefile 
# SKOFL_ROOT = /skofl
#

ifndef SKOFL_ROOT
  SKOFL_ROOT = ../..
endif

include $(SKOFL_ROOT)/config.gmk

LOCAL_INC	=

LOCAL_LIBS	= sktqconv_qb.o

FCFLAGS += -fno-range-check

#
#  Objects
#

OBJS   = sktqconv_qb.o llaser_qb #slaser_qb llaser_led tcheck

all: $(OBJS) 

clean: 
	$(RM) *.o *~ core* fort.* *.log *.hbk $(OBJS) *.dat *.table*

install.exec: 

