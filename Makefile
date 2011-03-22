#!/bin/sh
##############################################################################
##############################################################################
#                                                                            #
#                SPHRAY Makefile - Written by Gabriel Altay                  #
#                                                                            #
##############################################################################
##############################################################################

NULL = 
OPT = $(NULL)

include Makefile.local

ifeq (OWLS,$(findstring OWLS, $(OPT)))	
   OPT += -DincHmf
   OPT += -DuseHDF5     
   OPT += -DoutGammaHI
   OPT += -DincCloudy
   OPT += -DincEOS
endif


ifeq (GIMIC,$(findstring GIMIC, $(OPT)))	
   OPT += -DincHmf
   OPT += -DuseHDF5     
   OPT += -DoutGammaHI
   OPT += -DincCloudy
   OPT += -DincEOS
endif


APPS= screen sphray 

#--------------------

READERS =  \
		readers/gadget_public_header_class.o \
		readers/ion_table_class.o \
		readers/gadget_owls_header_class.o \
		readers/gadget_public_header_hdf5_class.o \
		readers/gadget_public_input.o \
		readers/gadget_cosmoBH_input.o \
		readers/gadget_owls_input.o \
		readers/gadget_vbromm_input.o \
		readers/gadget_public_input_hdf5.o \
		$(NULL)

RATES = \
		rates/hui_gnedin_atomic_rates.o \
		rates/cen_atomic_rates.o \
		rates/hummer_atomic_rates.o  \
		rates/atomic_rates.o \
		$(NULL)

OBJS =  myf03.o mt19937.o m_mrgrnk.o \
		$(RATES) \
		gadget_constants.o physical_constants.o cosmology.o \
		b2cd.o spectra.o \
		particle_system.o octtree3.o \
		gadget_sphray_header_class.o \
		ray.o \
		global.o config.o \
		raylist.o \
		ionpar.o euler.o bdf.o ion_temperature_update.o \
		$(READERS) update_particles.o source_input.o main_input.o output.o \
		initialize.o \
		resolve.o \
		mainloop.o \
		$(NULL)

#=============================================================================
# These lines should NOT be edited if you are using the HDF5 libraries.
# Instead, simply set the variable DIRHDF, preferably in one of the 
# include files specified above.
#


ifeq (useHDF5,$(findstring useHDF5, $(OPT)))
   INCHDF= $(ISINC)$(DIRHDF)/include $(ISINC)$(DIRHDF)/lib
   MODHDF= $(ISMOD)$(DIRHDF)/include
   RLIBHDF= $(ISRLIB)$(DIRHDF)/lib
   LIBHDF= $(ISLIB)$(DIRHDF)/lib 
   OPTHDF= $(INCHDF) $(MODHDF) $(RLIBHDF) $(LIBHDF)  -lhdf5 -lhdfwrapper
endif
OPT += $(OPTHDF)


# Targets
#=============================================================================
#
# Implicit Rules
#----------------
#
# $@ = name of target 
# $< = name of first dependency
# $^ = name of all dependencies with duplicates removed
# $? = name of all dependencies newer than the target
# $+ = name of all dependencies w/o duplicates removed
#

F2OBJ= -c    # Fortran flag to compile to object without linking
FNAME= -o    # Fortran flag to name output file

all:$(APPS)

# Update the local Makefile.local
# ========================================
Makefile.local: Makefile.local.template
	@if ! [ -f Makefile.local ]; then \
		cp $< $@ && \
		touch $< && \
		echo $@ is created from Makefile.local.template; \
	else \
		echo $< updated; \
	fi;
	@echo Please edit/touch Makefile.local before make again;
	@exit 1


# Just Dummy Reporting
#=============================================================================
screen: Makefile
	@echo
	@echo "OPT=    " $(OPT)
	@echo

#
# Main SPHRAY application
#=============================================================================
sphray: $(OBJS) sphray.o
	$(FC) -I rates/ -I readers/ $(FFLAGS) $^ $(OPT) $(FNAME) $@


# HDF5 modules
#=============================================================================
gadget_input_hdf5.o: gadget_input_hdf5.F90
	$(FC) $(FFLAGS) $(OPT) $(F2OBJ) $< $(FNAME) $@


# Implicit Rules
#=============================================================================

%.o: %.F90 
	$(FC) $(FFLAGS) $(OPT) $(F2OBJ) $< $(FNAME) $@


# Standard Cleaning Targets
#=============================================================================

clean :
	rm -f $(OBJS) *.mod */*.mod

cleanall :
	rm -f $(OBJS) *.mod */*.mod $(APPS) 

tidy :
	rm -f *~ 












