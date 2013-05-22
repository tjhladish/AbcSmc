CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2 -std=c++0x -Wall --pedantic
#CFLAGS = -O2 -std=c++0x
ABCDIR = $(HOME)/work/AbcSmc
INCLUDE = -I/usr/include/eigen3/ -I$(ABCDIR)
LIBS    = -lm -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 
SOURCES= AbcSmc.cpp utility.cpp 

LIB=libabc.a

all: $(SOURCES) $(LIB) abc

OBJECTS=$(SOURCES:.cpp=.o)
	
$(LIB): $(OBJECTS) 
	$(AR) -rv $(LIB) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ 

abc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) main.cpp -o abc $(LIBS) 

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)

clean:
	rm $(OBJECTS) $(LIB)
