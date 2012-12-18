# Makefile for libcosmo.a, a library used for cosmological functions
# and calculations
include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS}

AR?=ar

.PHONY: all clean;

OBJS = log_norm.o power_law.o statistics.o mathtools.o

all: libmath.a

clean:
	rm -rf libmath.a log_norm.o power_law.o mathtools.o statistics.o

libmath.a: ${OBJS} 
	${AR} -r libmath.a ${OBJS}
