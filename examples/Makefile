# This Makefile demonstrates different options for using the AbcSmc library
# with different approaches to the simulator.
# AbcSmc supports three different ways to use the simulator:
#
# 1. The simulator is compiled into the executable (main_sql.cpp / abc_sql).
# 2. The simulator is compiled into a shared library (dice.cpp / dice.so) and loaded dynamically (main_dynamic.cpp / abc_dynamic).
# 3. The simulator is compiled into a separate executable (dice_game.cpp) and invoked externally (main_exec.cpp / abc_exec).
#


CPP := g++

CFLAGS := -O2 -Wall -std=c++17 --pedantic
SHAREDFLAGS := $(CFLAGS) -fPIC -shared -rdynamic -Wl,--no-as-needed

MKFILE_PATH = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR := $(MKFILE_PATH)/..
SQLDIR := $(ABCDIR)/sqdb
JSNDIR := $(ABCDIR)/jsoncpp/include
GSL_LIB ?= -lm -lgsl -lgslcblas -pthread

ABC_INC := -I$(ABCDIR) -I$(SQLDIR) -I$(JSNDIR) -L$(ABCDIR) -labc $(GSL_LIB) -ldl

#ifdef TACC_GSL_INC
#GSL_INC = -I$$TACC_GSL_INC
#endif
#ifdef HPC_GSL_INC
#GSL_INC = -I$$HPC_GSL_INC
#endif

default: ../libabc.a dice_game dice.so abc_sql abc_dynamic abc_exec

../libabc.a:
	$(MAKE) -C $(ABCDIR) -f Makefile

abc_%: main_%.cpp ../libabc.a
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC)

dice_game: dice_game.cpp dice.h
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC) $(GSL_LIB)

dice.so: dice.cpp dice.h ../libabc.a
	$(CPP) $(SHAREDFLAGS) $< -o $@ $(ABC_INC)
	nm -gD $@ | grep simulator # check that the simulator is exported

run_%: abc_% | %.json
	./$< $| --process --simulate -n 1000
	./$< $| --process

clean:
	$(MAKE) -C $(ABCDIR) clean
	rm -f abc_sql abc_exec abc_dynamic dice_game dice.so
