include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS = general.o

all: libparallel.a

libparallel.a: ${OBJS} 
	${AR} -r libparallel.a ${OBJS}

clean: 
	rm -rf libparallel.a *.o *.a
