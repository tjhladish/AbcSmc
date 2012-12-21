CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2 -std=c++0x
INCLUDE = -I$(HOME)/work/lib/eigen/
#LDFLAGS = $(HOME)/work/EpiFire/src/*.o
LIBS    = -lm -lgsl -lgslcblas

all: abc 

abc: AbcSmc.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) AbcSmc.cpp utility.cpp -o abc $(LIBS) 

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)
