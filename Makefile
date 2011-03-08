#!/bin/sh
##############################################################################
##############################################################################
#                                                                            #
#                SPHRAY Makefile - Written by Gabriel Altay                  #
#                                                                            #
##############################################################################
##############################################################################

NULL = 
# Runtime Macros 
#=============================================================================
#
# These control the particle data structure and the code operation 
# in the case of the recombination macros.  The default is for all
# these flags to be commented out. 

#OPT += -DincHe      # if you want to include Helium
#OPT += -DoutGammaHI # if you want to output HI photoionization rate

#OPT += -DincHmf     # if you want individual Hydrogen mass fractions
#OPT += -DincHemf    # if you want individual Helium mass fractions

#OPT += -DincEOS     # store equation of state variable 
#OPT += -DincSFR     # store star formation rate
#OPT += -DincCloudy  # store CLOUDY eq. ionization values for particles

#OPT += -DuseHDF5    # if you are using the hdf5 library for anything

#OPT += -DOWLS       # selects input as an OWLS snapshot
#OPT += -DGIMIC      # selects input as a GIMIC snapshot


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



# System Definitions
#=============================================================================
#
# In order for the Makefile to work you have to design an include file
# for your system that sets certain variables.  Please look in the 
# makes directory for examples.  Right now, I have examples for three 
# fortran compilers.  gfortran is in make.blue-velvet.serial.opt, ifort
# is in make.santabuntu.serial.opt, and the sun compiler is in 
# make.icc-graphics.serial.opt
#

#include makes/make.icc-graphics.serial.debug
#include makes/make.icc-graphics.serial.opt

#include makes/make.xgpc.serial.debug
#include makes/make.xgpc.serial.opt

#include makes/make.ferrari.serial.opt

#include makes/make.cosma.serial.debug
#include makes/make.cosma.serial.opt

#include makes/make.santabuntu.serial.opt
#include makes/make.santabuntu.serial.debug

#include makes/make.blue-velvet.serial.gfortran.debug

#include makes/make.blue-velvet.serial.intel.debug
#include makes/make.blue-velvet.serial.intel.opt

#include makes/make.titania.serial.debug

include makes/make.warp.opt


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
		physical_constants.o cosmology.o \
		b2cd.o spectra.o \
		particle_system.o octtree3.o \
		gadget_general_class.o gadget_sphray_header_class.o \
		global.o config.o \
		ray.o raylist.o \
		ionpar.o euler.o bdf.o ion_temperature_update.o \
		$(READERS) update_particles.o source_input.o main_input.o output.o \
		initialize.o \
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


# Just Dummy Reporting
#=============================================================================
screen: Makefile
	@echo
	@echo "OPT=    " $(OPT)
	@echo "SYSTEM= " $(SYSTEM)
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












