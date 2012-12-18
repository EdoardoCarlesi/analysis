# Makefile for libio.a a library to handle different file input and output
include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS = read_io.o halo_io.o power_io.o write_io.o

all: libio.a

libio.a: ${OBJS} 
	${AR} -r libio.a ${OBJS}

clean: 
	rm -rf libio.a *.o *.a
