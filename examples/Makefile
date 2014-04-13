#CC=g++
CC=mpicxx
#CFLAGS = -g -std=c++0x
#CFLAGS = -O2 -std=c++0x
CFLAGS =  -std=c++11 -O2 -w0
ABCDIR = $(HOME)/AbcSmc
ABC_INC = -I$(HOME)/eigen -I$(ABCDIR)
GSL_INC = -I$$TACC_GSL_INC
#INCLUDE = -I/usr/include/eigen3 -I$(ABCDIR)
#LIBS    = -lm -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp   
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -lgsl -lgslcblas

default: abc_mpi

all: libabc abc dice abc_mpi 

libabc:  
	$(MAKE) -C $(ABCDIR) -f Makefile

abc: main.cpp
	$(CC) $(CFLAGS) $(ABC_INC) $(GSL_INC) main.cpp -o abc $(ABC_LIB) $(GSL_LIB)

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(GSL_INC) dice_game.cpp -o dice_game $(GSL_LIB)

abc_mpi: libabc main_mpi.cpp
	$(CC) $(CFLAGS) -D USING_MPI $(ABC_INC) $(GSL_INC) main_mpi.cpp -o abc_mpi $(ABC_LIB) $(GSL_LIB)

