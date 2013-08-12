# Makefile for the ROOT test programs.
# This Makefile shows nicely how to compile and link applications
# using the ROOT libraries on all supported platforms.
#
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Fons Rademakers, 29/2/2000

ROOTSYS = ../../..
include $(ROOTSYS)/etc/Makefile.arch

#------------------------------------------------------------------------------

ifeq ($(PLATFORM),win32)
EXTRALIBS        = "$(ROOTSYS)/lib/libLVMini.lib"
else
EXTRALIBS        = -lLVMini
CXXFLAGS += -g      
endif

# for using with MPI
ifneq ($(USE_MPI),)
CXX=mpic++
LD=mpic++
endif
ifneq ($(USE_OPENMP),)
CXXFLAGS += -D_GLIBCXX_PARALLEL -fopenmp 
LDFLAGS += -fopenmp 
endif



MINIMIZEOBJ      = testMinimize.$(ObjSuf)
MINIMIZESRC      = testMinimize.$(SrcSuf)
MINIMIZE         = testMinimize$(ExeSuf)

STRAIGHTLINEFITOBJ      = testStraightLineFit.$(ObjSuf)
STRAIGHTLINEFITSRC      = testStraightLineFit.$(SrcSuf)
STRAIGHTLINEFIT         = testStraightLineFit$(ExeSuf)

WEIBULLFITOBJ      = testWeibullFit.$(ObjSuf)
WEIBULLFITSRC      = testWeibullFit.$(SrcSuf)
WEIBULLFIT         = testWeibullFit$(ExeSuf)

OBJS          = $(MINIMIZEOBJ) $(STRAIGHTLINEFITOBJ) $(WEIBULLFITOBJ)

PROGRAMS      = $(MINIMIZE) $(STRAIGHTLINEFIT) $(WEIBULLFIT)

.SUFFIXES: .$(SrcSuf) .$(ObjSuf) $(ExeSuf)


all:            $(PROGRAMS)

$(MINIMIZE): 	$(MINIMIZEOBJ)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(EXTRALIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(STRAIGHTLINEFIT): 	$(STRAIGHTLINEFITOBJ)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(EXTRALIBS) $(OutPutOpt)$@
		@echo "$@ done"

$(WEIBULLFIT): 	$(WEIBULLFITOBJ)
		$(LD) $(LDFLAGS) $^ $(LIBS) $(EXTRALIBS) $(OutPutOpt)$@
		@echo "$@ done"

clean:
		@rm -f $(OBJS) core

distclean:      clean
		@rm -f $(PROGRAMS)  *Dict.* *.def *.exp \
		   *.root *.ps *.so *.lib *.dll *.d .def so_locations

.SUFFIXES: .$(SrcSuf)


.$(SrcSuf).$(ObjSuf):
	$(CXX) $(CXXFLAGS) -c $<
