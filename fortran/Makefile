#-----------------------------------------------------------------------
# This is the makefile for EMC3 grid driver
#
#  - Intended to be run from the 'build' directory
#
# targets:
#  [default] -  Makes executable
#   debug    -  Makes exec with debug options
#   clean    -  Removes *.mod, *.o from build directory
#   cleaner  -  Same as clean AND removes  executables from bin directory
#   
# .mod, .o, and .a files are made in current directory
#
#	7/14/2011 - current JDL
#-----------------------------------------------------------------------

SRC_DIR = 
BIN_DIR = 

# Location of bfield library and .mod files
LIBBJDL_DIR = $(HOME)/util-library/fortran/bfield_library_jdl/build


# The following flags should work for the gfortran compiler
#FC = mpif90
FC = gfortran
FF = -O3 
FF_DEBUG = -g -frange-check -fbounds-check -Wall -Wline-truncation -Wcharacter-truncation -Wsurprising -Waliasing -Wimplicit-interface -Wunused-parameter -fcheck=all -fbacktrace

MYEXEC = $(BIN_DIR)run_shiny.exe

#-----------------------------------------------------------------------

# Libraries needed to build the executable
LIBDIR = -L$(LIBBJDL_DIR) 
LIBS = -lbjdl 

# Objects needed to build the executable that use the default rule
OBJS90 =  \
shiny_mod.o \
shiny.o

# Modules to be compiled which have identically named .f90 files
MODS90 = \
shiny_mod.mod 

#-----------------------------------------------------------------------

$(MYEXEC): $(OBJS90) $(LIBBJDL_DIR)/libbjdl.a 
	@echo "Creating "$@" in directory "$(BIN_DIR)
	${FC} -o $@ $(FF) $^ $(LIBDIR) $(LIBS)

.PHONY: debug
debug: FF=$(FF_DEBUG)
debug: $(MYEXEC)

#-----------------------------------------------------------------------
# Static patterns

$(OBJS90) : %.o:  $(SRC_DIR)%.f90
	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

$(MODS90) : %.mod:  $(SRC_DIR)%.f90
	$(FC) $(FF) -c -I$(LIBBJDL_DIR) $(LIBDIR) $(LIBS) $<

# dependencies
#shiny_mod.mod: grid_routines_3d.o

.PHONY: clean cleaner
clean:
	-rm *.o *.mod
	@echo "Cleaned object and mod files from "$(PWD)

cleaner:
	-rm *.o *.mod *.a $(MYEXEC)
	@echo "Cleaned object, mod, library and executable files from "$(PWD)
