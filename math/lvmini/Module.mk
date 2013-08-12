# Module.mk for LVMini module
# Copyright (c) 2000 Rene Brun and Fons Rademakers
#
# Author: Rene Brun, 07/05/2003

MODNAME       := lvmini
MODDIR        := $(ROOT_SRCDIR)/math/$(MODNAME)
MODDIRS       := $(MODDIR)/src
MODDIRI       := $(MODDIR)/inc

LVMINIDIR    := $(MODDIR)
LVMINIDIRS   := $(LVMINIDIR)/src
LVMINIDIRI   := $(LVMINIDIR)/inc
LVMINIDIRT   := $(call stripsrc,$(LVMINIDIR)/test)

##### libLVMini #####
LVMINIL     := $(MODDIRI)/LinkDef.h
LVMINIDS    := $(call stripsrc,$(MODDIRS)/G__LVMini.cxx)
LVMINIDO    := $(LVMINIDS:.cxx=.o)
LVMINIDH    := $(LVMINIDS:.cxx=.h)

LVMINIAH    := $(filter-out $(MODDIRI)/LinkDef%,$(wildcard $(MODDIRI)/*.h))
LVMINIBH    := $(filter-out $(MODDIRI)/LVMini/LinkDef%,$(wildcard $(MODDIRI)/LVMini/*.h))
LVMINIH     := $(LVMINIAH) $(LVMINIBH)
LVMINISCXX  := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.cxx))
LVMINISF    := $(filter-out $(MODDIRS)/G__%,$(wildcard $(MODDIRS)/*.f))
LVMINIO     := $(call stripsrc,$(LVMINISCXX:.cxx=.o)) $(call stripsrc,$(LVMINISF:.f=.o))

LVMINIDEP   := $(LVMINIO:.o=.d) $(LVMINIDO:.o=.d)

LVMINILIB   := $(LPATH)/libLVMini.$(SOEXT)
LVMINIMAP   := $(LVMINILIB:.$(SOEXT)=.rootmap)

# used in the main Makefile
ALLHDRS      += $(patsubst $(MODDIRI)/%.h,include/%.h,$(LVMINIH))
ALLLIBS      += $(LVMINILIB)
ALLMAPS      += $(LVMINIMAP)

# include all dependency files
INCLUDEFILES += $(LVMINIDEP)

##### local rules #####
.PHONY:         all-$(MODNAME) clean-$(MODNAME) distclean-$(MODNAME) \
                test-$(MODNAME)

include/LVMini/%.h: $(LVMINIDIRI)/LVMini/%.h
		@(if [ ! -d "include/LVMini" ]; then     \
		   mkdir -p include/LVMini;              \
		fi)
		cp $< $@

include/%.h:    $(LVMINIDIRI)/%.h
		cp $< $@

$(LVMINILIB):  $(LVMINIO) $(LVMINIDO) $(ORDER_) $(MAINLIBS) $(LVMINILIBDEP)
		@$(MAKELIB) $(PLATFORM) $(LD) "$(LDFLAGS)" \
		   "$(SOFLAGS)" libLVMini.$(SOEXT) $@ \
		   "$(LVMINIO) $(LVMINIDO)" \
		   "$(LVMINILIBEXTRA) $(F77LIBS)"

$(LVMINIDS):   $(LVMINIH) $(LVMINIL) $(ROOTCINTTMPDEP)
		$(MAKEDIR)
		@echo "Generating dictionary $@..."
		$(ROOTCINTTMP) -f $@ -c $(LVMINIH) $(LVMINIL)

$(LVMINIMAP):  $(RLIBMAP) $(MAKEFILEDEP) $(LVMINIL)
		$(RLIBMAP) -o $@ -l $(LVMINILIB) \
		   -d $(LVMINILIBDEPM) -c $(LVMINIL)

all-$(MODNAME):  $(LVMINILIB) $(LVMINIMAP)

test-$(MODNAME): all-$(MODNAME)
ifneq ($(ROOT_OBJDIR),$(ROOT_SRCDIR))
		@$(INSTALL) $(LVMINIDIR)/test $(LVMINIDIRT)
endif
		@cd $(LVMINIDIRT) && $(MAKE) ROOTCONFIG=../../../bin/root-config

clean-$(MODNAME):
		@rm -f $(LVMINIO) $(LVMINIDO)

clean::         clean-$(MODNAME)

distclean-$(MODNAME): clean-$(MODNAME)
		@rm -f $(LVMINIDEP) $(LVMINIDS) $(LVMINIDH) $(LVMINILIB) \
		   $(LVMINIMAP)
		@rm -rf include/LVMini
ifneq ($(ROOT_OBJDIR),$(ROOT_SRCDIR))
		@rm -rf $(LVMINIDIRT)
else
		@cd $(LVMINIDIRT) && $(MAKE) distclean ROOTCONFIG=../../../bin/root-config
endif

distclean::     distclean-$(MODNAME)

##### extra rules ######
$(LVMINIO): CXXFLAGS += -DWARNINGMSG -DUSE_ROOT_ERROR
$(LVMINIDO): CXXFLAGS += -DWARNINGMSG -DUSE_ROOT_ERROR

# Use specified fortran compiler and flags. I'm not sure why these are
# not picked up automatically. Maybe different make versions?
$(LVMINIO): FC=$(F77)
$(LVMINIO): FFLAGS += $(F77FLAGS)

# Optimize dictionary with stl containers.
$(LVMINIDO): NOOPT = $(OPT)
