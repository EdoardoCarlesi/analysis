################################################################
# Makefile for GADGET and AHF related reading routines,        #
# growth factor calculations, integrations and interpolations  #
# Edoardo Carlesi - edoardo.carlesi@uam.es 2010 - 2012	       #
################################################################
include ../Makefile.config

LIBS = libmath/libmath.a libio/libio.a libcosmo/libcosmo.a libhalo/libhalo.a
SRCS = general_variables.c general_functions.c
OBJ  = general_variables.o general_functions.o

ifeq (${WITH_MPI}, "on")
LIBS    += libparallel/libparallel.a
endif

ifeq (${ADD_LIB}, "on")
LIBS	+= -lgsl -lgslcblas -lm
endif

#EXECS = nfw_fit growth_factor halo_comparison halo_statistics halo_evolution subhalo_statistics theoretical_mass_function\
	mass_function number_density merger_tree

EXECS = halo_statistics

all:	${EXECS}
	
halo_statistics : halo_statistics.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o halo_statistics halo_statistics.c ${SRCS} ${LIBS} ${LIBS};
		mv halo_statistics ${HOME_DIR}/bin/		

test : test.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o test test.c ${SRCS} ${LIBS} ${LIBS};
		mv test ${HOME_DIR}/bin/		

clean: 
		rm -rf test *.txt *.o *~ ${HOME_DIR}/bin/*
		cd libmath ; ${MAKE} clean
		cd libio ; ${MAKE} clean
		cd libcosmo ; ${MAKE} clean
		cd libhalo ; ${MAKE} clean
		cd libparallel ; ${MAKE} clean
		
		
# Various libraries needed by the program
libmath/libmath.a: 
	cd libmath/ ; ${MAKE} libmath.a

libcosmo/libcosmo.a: 
	cd libcosmo/ ; ${MAKE} libcosmo.a

libio/libio.a: 
	cd libio/ ; ${MAKE} libio.a

libhalo/libhalo.a: 
	cd libhalo/ ; ${MAKE} libhalo.a

libparallel/libparallel.a: 
	cd libparallel/ ; ${MAKE} libparallel.a

