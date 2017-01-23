CPP:=g++
CFLAGS = -O2 -Wall -std=c++11 --pedantic -Wno-deprecated-declarations
MKFILE_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR = $(MKFILE_PATH)/..
ABC_INC = -I$(ABCDIR) -I../sqdb -I$(ABCDIR)/gsl-2.2.1

ifdef TACC_GSL_INC
GSL_INC = -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
GSL_INC = -I$$HPC_GSL_INC
endif

ABC_LIB = -L$(ABCDIR) -labc -ljsoncpp -lsqdb ../sqlite3.o
#GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl
GSL_LIB = -lm -L../gsl-2.2.1 -lgsl_subset -lpthread -ldl

default: libabc dice abc

libabc:
	$(MAKE) -C $(ABCDIR) -f Makefile

abc_sql: libabc main_sql.cpp
	$(CPP) $(CFLAGS) $(ABC_INC) $(GSL_INC) main_sql.cpp -o abc_simulator_sql $(ABC_LIB) $(GSL_LIB)

abc: libabc main_exec.cpp
	$(CPP) $(CFLAGS) $(ABC_INC) $(GSL_INC) main_exec.cpp -o abc_simulator_executable $(ABC_LIB) $(GSL_LIB)

dice: dice_game.cpp
	$(CPP) $(CFLAGS) $(GSL_INC) dice_game.cpp -o dice_game $(GSL_LIB)

clean:
	$(MAKE) -C $(ABCDIR) clean
	rm -f dice_game abc_simulator_pointer abc_simulator_executable
