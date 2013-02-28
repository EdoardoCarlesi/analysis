# Makefile for libhalo.a, a library used for halo analysis

include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS = halo_properties.o nfw.o subhalo_properties.o halo_subhalo_general.o

all: libhalo.a

libhalo.a: ${OBJS} 
	${AR} -r libhalo.a ${OBJS}

clean:
	rm -rf *.o libhalo.a 
