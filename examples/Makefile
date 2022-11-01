CPP:=g++
CFLAGS = -O2 -Wall -std=c++17 --pedantic -fPIC
MKFILE_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR = $(MKFILE_PATH)/..
ABC_INC = -I$(ABCDIR)

#ifdef TACC_GSL_INC
#GSL_INC = -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#GSL_INC = -I$$HPC_GSL_INC
#endif
GSL_INC = -I$(ABCDIR)/gsl_local/include

ABC_LIB = -L$(ABCDIR) -labc
#GSL_LIB = -lm -L$$TACC_GSL_LIB/ -L$$HPC_GSL_LIB/ -lgsl -lgslcblas -lpthread -ldl
GSL_LIB = -lm -lgsl -lgslcblas -lpthread -ldl

default: $(ABCDIR)/libabc.a dice abc_sql abc

$(ABCDIR)/libabc.a:
	$(MAKE) -C $(@D) -f Makefile $(@F)

abc_sql: $(ABCDIR)/libabc.a main_sql.cpp
	$(CPP) $(CFLAGS) $(ABC_INC) main_sql.cpp -o abc_simulator_sql $(ABC_LIB) $(GSL_LIB)

abc: $(ABCDIR)/libabc.a main_exec.cpp
	$(CPP) $(CFLAGS) $(ABC_INC) main_exec.cpp -o abc_simulator_executable $(ABC_LIB) $(GSL_LIB)

dice: dice_game.cpp
	$(CPP) $(CFLAGS) dice_game.cpp -o dice_game $(GSL_LIB)

clean:
	rm -f abc_simulator_sql abc_simulator_executable dice_game
