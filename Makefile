CC=g++
CFLAGS = -O2
INCLUDE = -I $(HOME)/work/EpiFire/src
LDFLAGS = $(HOME)/work/EpiFire/src/*.o
LIBS    = -lm -lgsl -lgslcblas

all: abc 

abc: AbcSmc.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) AbcSmc.cpp -o abc $(LIBS) 

