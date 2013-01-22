################################################################
# Makefile for GADGET and AHF related reading routines,        #
# growth factor calculations, integrations and interpolations  #
# Edoardo Carlesi - edoardo.carlesi@uam.es 2010 - 2012	       #
################################################################
include ../Makefile.config

LIBS = libmath/libmath.a libio/libio.a libcosmo/libcosmo.a libhalo/libhalo.a
SRCS = general_variables.c general_functions.c
OBJ  = general_variables.o general_functions.o
INCLUDE += general_variables.h general_functions.h

ifeq (${WITH_MPI}, "on")
LIBS    += libparallel/libparallel.a
INCLUDE += libparallel/general.h
endif

ifeq (${SYS}, "castor")
LIBS	+= -lgsl -lgslcblas -lm
endif

#EXECS = nfw_fit growth_factor halo_comparison halo_statistics halo_evolution subhalo_statistics theoretical_mass_function\
	mass_function number_density merger_tree

EXECS = growth_factor halo_statistics halo_evolution subhalo_statistics theoretical_mass_function\
	mass_function number_density

all:	${EXECS}

#nfw_fit : nfw_fit.c ${LIBS} ${SRCS}
#	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS) $(CFLAGS) ${LIBS} -o nfw_fit nfw_fit.c ${SRCS} ${LIBS} ${LIBS};
#	mv nfw_fit ${HOME_DIR}/bin/		

growth_factor : growth_factor.c ${LIBS} ${SRCS}
	$(CC) $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS) ${LIBS} -o growth_factor growth_factor.c ${SRCS} ${LIBS} ${LIBS};
	mv growth_factor ${HOME_DIR}/bin/

#halo_comparison : halo_comparison.c ${LIBS} ${SRCS}
#	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o halo_comparison halo_comparison.c  ${SRCS} ${LIBS} ${LIBS};
#		mv halo_comparison ${HOME_DIR}/bin/		
		
halo_statistics : halo_statistics.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o halo_statistics halo_statistics.c ${SRCS} ${LIBS} ${LIBS};
		mv halo_statistics ${HOME_DIR}/bin/		
		
halo_evolution : halo_evolution.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS) ${LIBS} -o halo_evolution halo_evolution.c ${SRCS} ${LIBS} ${LIBS};
		mv halo_evolution ${HOME_DIR}/bin/		

mass_function : mass_function.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o mass_function mass_function.c ${SRCS} ${LIBS} ${LIBS};
		mv mass_function ${HOME_DIR}/bin/		

#merger_tree : merger_tree.c ${LIBS} ${SRCS}
#	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o merger_tree merger_tree.c ${SRCS} ${LIBS} ${LIBS};
#		mv merger_tree ${HOME_DIR}/bin/		

number_density : number_density.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o number_density number_density.c ${SRCS} ${LIBS} ${LIBS};
		mv number_density ${HOME_DIR}/bin/		

subhalo_statistics : subhalo_statistics.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS) ${LIBS} -o subhalo_statistics subhalo_statistics.c ${SRCS} ${LIBS} ${LIBS};
		mv subhalo_statistics ${HOME_DIR}/bin/		

theoretical_mass_function : theoretical_mass_function.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS) ${LIBS} -o theoretical_mass_function theoretical_mass_function.c ${SRCS} ${LIBS} ${LIBS};
		mv theoretical_mass_function ${HOME_DIR}/bin/		

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

