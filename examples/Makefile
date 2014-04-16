CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2
ABCDIR = $(HOME)/work/AbcSmc
ABC_INC = -I$(ABCDIR)
#GSL_INC = -I$$TACC_GSL_INC
GSL_INC = -I$$HPC_GSL_INC
ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp   
#GSL_LIB = -lm -L$$TACC_GSL_LIB/ -lgsl -lgslcblas
GSL_LIB = -lm -L$$HPC_GSL_LIB/ -lgsl -lgslcblas

default: all_no_mpi

all_no_mpi: CFLAGS += -std=c++0x -Wall
all_no_mpi: libabc_no_mpi dice abc 

all_mpi: CC = mpicxx
all_mpi: CFLAGS += -cxx=icc -std=c++0x -w0 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
all_mpi: libabc_mpi dice abc_mpi 

libabc_no_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

libabc_mpi:  
	$(MAKE) -C $(ABCDIR) -f Makefile all_mpi

abc: main.cpp
	$(CC) $(CFLAGS) $(ABC_INC) $(GSL_INC) main.cpp -o abc $(ABC_LIB) $(GSL_LIB)

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(GSL_INC) dice_game.cpp -o dice_game $(GSL_LIB)

abc_mpi: libabc_mpi main_mpi.cpp
	$(CC) $(CFLAGS) -D USING_MPI $(ABC_INC) $(GSL_INC) main_mpi.cpp -o abc_mpi $(ABC_LIB) $(GSL_LIB)

