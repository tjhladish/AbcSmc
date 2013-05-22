CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2 -std=c++0x
HOME = /home/tjhladish/
ABCDIR = $(HOME)/work/AbcSmc/
#INCLUDE = -I$(HOME)/work/lib/eigen/
INCLUDE = -I/usr/include/eigen3/ -I$(ABCDIR)
#LDFLAGS = $(HOME)/work/EpiFire/src/*.o
LIBS    = -lm -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 

all: abc 

abc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) main.cpp -o abc $(LIBS) 

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)
