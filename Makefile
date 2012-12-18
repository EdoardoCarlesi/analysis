######################################
# Makefile for Data Analysis Program #	
######################################
include Makefile.config

CC      = gcc
CFLAGS  = -Wall -g
LDFLAGS = 
SOURCES = main.c

OBJECTS        = $(SOURCES:.c = .o) 

all: 	
	cd src/; ${MAKE}

#EXEC        = analysis

#$(EXEC): $(OBJECTS) 
#		$(CC) $(OPT) $(LFLAGS) $(CFLAGS) $(OBJECTS) -o $@;
#	mv analysis bin
clean:
	cd src; make clean; cd ..; rm -rf analysis
