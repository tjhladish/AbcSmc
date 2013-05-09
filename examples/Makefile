CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2 -std=c++0x
ABCDIR = $(HOME)/work/AbcSmc
INCLUDE = -I/usr/include/eigen3 -I$(ABCDIR)
LIBS    = -lm -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 

all: abc dice

abc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) main.cpp -o abc $(LIBS) 

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)
