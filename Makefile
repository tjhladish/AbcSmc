CC=mpicxx
#CC=g++
#CFLAGS = -g -std=c++0x
#CFLAGS = -O2 -std=c++0x -Wall --pedantic
#CFLAGS = -O2 -std=c++0x
CFLAGS =  -w0 -std=c++11 -O2 -D USING_MPI
ABCDIR = $(HOME)/AbcSmc
#INCLUDE = -I/usr/include/eigen3/ -I$(ABCDIR)
INCLUDE = -I$(HOME)/eigen -I$(ABCDIR) -I$$TACC_GSL_INC -I$(HOME)/jsoncpp/include
LIBS    = -lm -L$$TACC_GSL_LIB/ -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 
SOURCES= AbcSmc.cpp utility.cpp 

LIB=libabc.a

all: $(SOURCES) $(LIB) #abc

OBJECTS=$(SOURCES:.cpp=.o)
	
$(LIB): $(OBJECTS) 
	$(AR) -rv $(LIB) $(OBJECTS) 

.cpp.o:
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ 

#abc: main.cpp
#	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) main.cpp -o abc $(LIBS) 

#dice: dice_game.cpp
#	$(CC) $(CFLAGS) $(INCLUDE) $(LDFLAGS) dice_game.cpp -o dice_game $(LIBS)

clean:
	rm $(OBJECTS) $(LIB)
