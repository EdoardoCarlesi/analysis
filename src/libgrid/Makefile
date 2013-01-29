include ../../Makefile.config

CFLAGS+=${OPT} ${CCFLAGS} ${DEFINEFLAGS}

MAKE?=make

AR?=ar

.PHONY: all clean;

OBJS = 

all: 	libgrid.a

libgrid.a: ${OBJS} 
	${AR} -r libgrid.a ${OBJS}

clean:
	rm -rf *.o libgrid.a
