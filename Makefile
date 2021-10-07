###############################################################
# Automatically-generated sample makefile to illustrate how to  
# link against oomph-lib from outside the automake/autoconf
# framework. Do not edit this -- make a copy first
# 
# When customising this makefile, you should only have to change
# 
# - the variable OOMPH-LIB_LIBS:
#         Add any additional oomph-lib sub-libraries that 
#         you may wish to use in your code. 
# 
# - the specific dependencies for your driver code: 
#         Include any additional local dependencies such as 
#         your own header files etc.
# 
###############################################################
 
 
# Installation-specific information -- don't change any of this! 
#-------------------------------------------------------------- 
OOMPH-LIB_INSTALL_LOCATION=/home/paul/Repositories/oomph-lib-jack
 
# Flags for C pre-processor 
AM_CPPFLAGS=-DHAVE_CONFIG_H -isystem $(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/trilinos/trilinos_default_installation/include -DOOMPH_HAS_TRILINOS -isystem $(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/boost/boost_default_installation/include -DOOMPH_HAS_BOOST -isystem $(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/gmp/gmp_default_installation/include -DOOMPH_HAS_GMP -isystem $(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/mpfr/mpfr_default_installation/include -DOOMPH_HAS_MPFR -isystem $(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/cgal/cgal_default_installation/include -DOOMPH_HAS_CGAL -DOOMPH_HAS_STACKTRACE -DOOMPH_HAS_UNISTDH -DOOMPH_HAS_FPUCONTROLH -DOOMPH_HAS_MALLOCH -DOOMPH_HAS_TRIANGLE_LIB -DOOMPH_HAS_TETGEN_LIB -DUSING_OOMPH_SUPERLU -DUSING_OOMPH_SUPERLU_DIST -I$(OOMPH-LIB_INSTALL_LOCATION)/build/include
 
# Flags for C++ compiler 
CXXFLAGS= -O3 -Wall
 
# Library include directory: This is where all the header files live
OOMPH-LIB_INCLUDE_DIR=$(OOMPH-LIB_INSTALL_LOCATION)/build/include -I$(OOMPH-LIB_INSTALL_LOCATION)
 
# Library directory: This is where all of oomph-lib's sub-libraries live
OOMPH-LIB_LIB_DIR=$(OOMPH-LIB_INSTALL_LOCATION)/build/lib
 
# These are the external (3rd party) libraries that are distributed
# with oomph-lib and that we always link against
OOMPH-LIB_EXTERNAL_LIBS=-lml -lifpack -lamesos -lanasazi -laztecoo -lepetraext -ltriutils -lepetra -lteuchosremainder -lteuchosnumerics -lteuchoscomm -lteuchosparameterlist -lteuchoscore -lpthread -lboost_thread -lboost_system -lgmp -lmpfr -lCGAL_Core -lCGAL -loomph_hsl -loomph_arpack -loomph_crbond_bessel -loomph_triangle -loomph_tetgen -loomph_superlu_4.3 -loomph_metis_from_parmetis_3.1.1 -loomph_lapack -loomph_flapack -loomph_blas
 
# This specifies where libraries built from third party 
# distributions can be found
EXTERNAL_DIST_LIBRARIES=-L$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/trilinos/trilinos_default_installation/lib -L$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/boost/boost_default_installation/lib -Wl,-rpath,$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/boost/boost_default_installation/lib -L$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/gmp/gmp_default_installation/lib -L$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/mpfr/mpfr_default_installation/lib -L$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/cgal/cgal_default_installation/lib -Wl,-rpath,$(OOMPH-LIB_INSTALL_LOCATION)/external_distributions/cgal/cgal_default_installation/lib
 
# This is additional machine-specific linking information that 
# allows mixed-language compilation/linking
FLIBS=-L/usr/lib/gcc/x86_64-linux-gnu/9 -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../x86_64-linux-gnu -L/usr/lib/gcc/x86_64-linux-gnu/9/../../../../lib -L/lib/x86_64-linux-gnu -L/lib/../lib -L/usr/lib/x86_64-linux-gnu -L/usr/lib/../lib -L/usr/lib/gcc/x86_64-linux-gnu/9/../../.. -lgfortran -lm -lquadmath
 
# Flags required for the use of shared libraries 
# Linux: Replace the above line with the following
SHARED_LIBRARY_FLAGS=-Wl,--rpath -Wl,$(OOMPH-LIB_INSTALL_LOCATION)/build/lib
 
# Mac OSX without clang(?): Replace the above line with the following
#SHARED_LIBRARY_FLAGS= --rpath=$(OOMPH-LIB_INSTALL_LOCATION)/build/lib
 
# Problem-specific information -- edit this for your driver code 
 
#---------------------------------------------------------------
# These are the specific oomph-lib sub-libraries that we have to link against
# for this driver code -- edit this according to your requirements
# but remember that the order of the libraries matters: List the
# the more specific ones before the more general ones!
OOMPH-LIB_LIBS=-lconstitutive -lsolid -lrigid_body -lfluid_interface -lgeneric 


SRC_DIR=src/
DRIVERS_DIR=drivers/
BIN_DIR=bin/
BUILD_DIR=build/
LIB_DIR=$(BUILD_DIR)lib/

INCLUDE_DIR= -I$(SRC_DIR)

SRC=$(foreach individual_dir,$(SRC_DIR),$(wildcard $(individual_dir)*.cpp))
OBJ=$(patsubst $(SRC_DIR)%.cpp,$(LIB_DIR)%.o,$(SRC))

DRIVERS=$(foreach individual_dir,$(DRIVERS_DIR),$(wildcard $(individual_dir)*.cpp))
DRIVERS_OBJ=$(patsubst $(DRIVERS_DIR)%.cpp,$(DRIVERS_DIR)%.o,$(DRIVERS))

OUT=$(patsubst $(DRIVERS_DIR)%.cpp,$(BIN_DIR)%.out,$(DRIVERS))

CC=g++
LD=g++

build: $(OUT)

test: build
	scripts/validate.sh

# Specify targets without explicit outputs as phony
.PHONY: clean clean-all doc

clean:
	rm -rf build/

clean-all: clean
	rm -rf bin/

run_bubble_unsteady:
	$(BIN_DIR)bubble_unsteady.out -n 4 -r 0.5 -c 0.01 -q 1 -o data/bubble_unsteady/ > \
		OUTPUT_bubble_unsteady

doc:
	doxygen Doxyfile

check-lib-dir: $(LIB_DIR)

$(LIB_DIR):
	mkdir -p $@

check-bin-dir: $(BIN_DIR)

$(BIN_DIR):
	mkdir -p $@

# -----------------------------------------------------------------------------

# Compile sources
$(LIB_DIR)%.o: $(SRC_DIR)%.cpp check-lib-dir
	mkdir -p $(LIB_DIR)
	$(CC) $(AM_CPPFLAGS)  $(CXXFLAGS) -c $< -o $@ \
	       -I$(OOMPH-LIB_INCLUDE_DIR) $(INCLUDE_DIR) \
		   -L$(OOMPH-LIB_LIB_DIR) $(EXTERNAL_DIST_LIBRARIES) $(OOMPH-LIB_LIBS) \
	        $(OOMPH-LIB_EXTERNAL_LIBS) $(FLIBS) 
  
# Compile drivers
$(DRIVERS_DIR)%.o: $(DRIVERS_DIR)%.cpp
	$(CC) $(AM_CPPFLAGS)  $(CXXFLAGS) -c $< -o $@ \
	       -I$(OOMPH-LIB_INCLUDE_DIR) $(INCLUDE_DIR) \
		   -L$(OOMPH-LIB_LIB_DIR) $(EXTERNAL_DIST_LIBRARIES) $(OOMPH-LIB_LIBS) \
	        $(OOMPH-LIB_EXTERNAL_LIBS) $(FLIBS) 

# Link and create executables
$(BIN_DIR)%.out: $(OBJ) $(DRIVERS_OBJ)
	mkdir -p $(BIN_DIR)
	$(LD) $(SHARED_LIBRARY_FLAGS) $^ -o $@ \
			$(INCLUDE_DIR) \
	       -L$(OOMPH-LIB_LIB_DIR) $(EXTERNAL_DIST_LIBRARIES) $(OOMPH-LIB_LIBS) \
	        $(OOMPH-LIB_EXTERNAL_LIBS) $(FLIBS) 
