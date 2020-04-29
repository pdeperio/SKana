#ifndef ATMPD_ROOT
#  ATMPD_ROOT = ../../..
#endif

#include $(ATMPD_ROOT)/config.gmk

ROOTLIBS     := $(shell root-config --glibs) -lTreePlayer
ROOTCFLAGS   := $(shell root-config --cflags) 

SITE_INCLUDES  = $(ROOTCFLAGS)
CXXINCLUDES = $(SITE_INCLUDES) -I.
CXXFLAGS =  $(CXXINCLUDES) 

SITE_LDOPTIONS += -Wl,-z -Wl,muldefs
LDOPTS	= $(SITE_LDOPTIONS)


LOCAL_INC       =
#
#  Objects
#

#$(ODIR)/global.o $(ODIR)/Official.o $(ODIR)/OfficialEventParser.o

OBJS   = Official.o global.o OfficialEventParser.o run_compare.o
OBJPLOT = run_compare.o
#OBJPLOT = global.o
#OBJCOMB = global.o 

LDLIBS_NOSK =  $(LDOPTS) $(LOCAL_LIBS) $(ROOTLIBS) $(SITE_LIBRARIES) $(SITE_LDLIBS)

all: $(OBJS) $(OBJPLOT) $(OBJCOMB) exec

clean: 
	$(RM) *.o *~ core $(OBJS) $(OBJPLOT) $(OBJCOMB) official run_compare

official: $(OBJS)
	LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o official $(OBJS) $(LDLIBS_NOSK)

run_compare: $(OBJPLOT)
	LD_RUN_PATH=$(LIBDIR) $(CXX) $(CXXFLAGS) -o run_compare $(OBJPLOT) $(LDLIBS_NOSK) 

exec: official run_compare
