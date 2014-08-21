CC=g++
#CFLAGS = -g -std=c++0x
CFLAGS = -O2
ABCDIR = $(HOME)/work/AbcSmc
ABC_INC = -I$(ABCDIR)

ifdef TACC_GSL_INC
GSL_INC = -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
GSL_INC = -I$$HPC_GSL_INC
endif

ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp
GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas

default: all_no_mpi

all_no_mpi: CFLAGS += -std=c++0x -Wall
all_no_mpi: libabc_no_mpi dice abc

all_mpi: CC = mpicxx
all_mpi: CFLAGS += -cxx=icc -std=c++0x -w0 -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
all_mpi: libabc_mpi dice abc_mpi_exec abc_mpi_simptr

all_ubuntu_mpi: CC = mpicxx
all_ubuntu_mpi: CFLAGS += -std=c++0x -Wall -DMPICH_IGNORE_CXX_SEEK -DMPICH_SKIP_MPICXX
all_ubuntu_mpi: libabc_ubuntu_mpi abc_mpi_exec abc_mpi_simptr

libabc_no_mpi:
	$(MAKE) -C $(ABCDIR) -f Makefile all_no_mpi

libabc_mpi:
	$(MAKE) -C $(ABCDIR) -f Makefile all_mpi

libabc_ubuntu_mpi:
	$(MAKE) -C $(ABCDIR) -f Makefile all_ubuntu_mpi

abc: libabc_no_mpi main_exec.cpp main_simptr.cpp
	$(CC) $(CFLAGS) $(ABC_INC) $(GSL_INC) main_exec.cpp -o abc_simulator_executable $(ABC_LIB) $(GSL_LIB)
	$(CC) $(CFLAGS) $(ABC_INC) $(GSL_INC) main_simptr.cpp -o abc_simulator_pointer $(ABC_LIB) $(GSL_LIB)

dice: dice_game.cpp
	$(CC) $(CFLAGS) $(GSL_INC) dice_game.cpp -o dice_game $(GSL_LIB)

abc_mpi_exec: libabc_mpi main_mpi_exec.cpp
	$(CC) $(CFLAGS) -D USING_MPI $(ABC_INC) $(GSL_INC) main_mpi_exec.cpp -o abc_mpi_exec $(ABC_LIB) $(GSL_LIB)

abc_mpi_simptr: libabc_mpi main_mpi_simptr.cpp
	$(CC) $(CFLAGS) -D USING_MPI $(ABC_INC) $(GSL_INC) main_mpi_simptr.cpp -o abc_mpi_simptr $(ABC_LIB) $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) clean
	rm -f dice_game abc_simulator_pointer abc_simulator_executable abc_mpi_exec abc_mpi_simptr
