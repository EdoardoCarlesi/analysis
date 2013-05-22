################################################################
# Makefile for GADGET and AHF related reading routines,        #
# growth factor calculations, integrations and interpolations  #
# Edoardo Carlesi - edoardo.carlesi@uam.es 2010 - 2012	       #
################################################################
include ../Makefile.config

LIBS = ../src/libmath/libmath.a ../src/libio/libio.a ../src/libcosmo/libcosmo.a ../src/libhalo/libhalo.a
SRCS = ../src/general_variables.c ../src/general_functions.c
OBJ  = ../src/general_variables.o ../src/general_functions.o

ifeq (${WITH_MPI}, "on")
LIBS    += ../src/libparallel/libparallel.a
endif

ifeq (${ADD_LIB}, "on")
LIBS	+= -lgsl -lgslcblas -lm
endif

EXECS = residuals

all:	${EXECS}

residuals : residuals.c ${LIBS} ${SRCS}
	$(CC)  $(OPT) $(INCLUDE) $(BASELIBS)  $(CFLAGS)  ${LIBS} -o residuals residuals.c ${SRCS} ${LIBS} ${LIBS};
		mv residuals ${HOME_DIR}/bin/		

clean: 
		rm -rf *.o *~ ${HOME_DIR}/bin/*
		rm -rf ../src/*.o 
		cd ../src/libmath ; ${MAKE} clean
		cd ../src/libio ; ${MAKE} clean
		cd ../src/libcosmo ; ${MAKE} clean
		cd ../src/libhalo ; ${MAKE} clean
		cd ../src/libparallel ; ${MAKE} clean
		
# Various libraries needed by the program
../src/libmath/libmath.a: 
	cd ../src/libmath/ ; ${MAKE} libmath.a

../src/libcosmo/libcosmo.a: 
	cd ../src/libcosmo/ ; ${MAKE} libcosmo.a

../src/libio/libio.a: 
	cd ../src/libio/ ; ${MAKE} libio.a

../src/libhalo/libhalo.a: 
	cd ../src/libhalo/ ; ${MAKE} libhalo.a

../src/libparallel/libparallel.a: 
	cd ../src/libparallel/ ; ${MAKE} libparallel.a

