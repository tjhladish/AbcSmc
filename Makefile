CC=g++
#CFLAGS = -g -std=c++0x
#CFLAGS = -O2 -std=c++0x -Wall --pedantic
CFLAGS = -O2 -std=c++0x
#INCLUDE = -I$(HOME)/work/lib/eigen/
INCLUDE = -I/usr/include/eigen3/ -I/usr/include
LIBS    = -lm -lgsl -lgslcblas -ljsoncpp
SOURCES= AbcSmc.cpp utility.cpp 

LIB=libabc.a

all: $(SOURCES) $(LIB) abc

OBJECTS=$(SOURCES:.cpp=.o)
	
$(LIB): $(OBJECTS) 
	$(AR) -rv $(LIB) $(OBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ 

abc: main.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) main.cpp AbcSmc.cpp utility.cpp -o abc $(LIBS) 

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)

clean:
	rm $(OBJECTS) $(LIB)
