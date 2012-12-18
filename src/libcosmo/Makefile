# Makefile for libcosmo.a, a library used for cosmological functions
# and calculations
include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS = cosmological_relations.o mass_function.o tinker.o number_density.o power.o

all: 	libcosmo.a

libcosmo.a: ${OBJS} 
	${AR} -r libcosmo.a ${OBJS}

clean:
	rm -rf *.o libcosmo.a
