
default: all

# This Makefile demonstrates different options for using the AbcSmc library
# with different approaches to the simulator, as well as a toy simulator for use with AbcSmc.
#
# The toy simulator represents rolling a number dice, all of the same face count, then returns their
# sum and standard deviation. The simulator is defined in dice.h, and follows the function signature
# required by AbcSmc.
#
# AbcSmc supports three different ways to use a simulator:
#
# 1. The simulator is compiled into an executable along with the library (main_sql.cpp --> abc_sql).
# 2. The simulator is compiled into a shared library (dice.cpp / dice.so) and loaded dynamically (main_dynamic.cpp --> abc_dynamic).
# 3. The simulator is compiled into a separate executable (dice_game.cpp) and invoked externally (main_exec.cpp --> abc_exec).
#
# This Makefile demonstrates how to compile and run each of these options.

# make variable definitions
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

DEMOS := sql dynamic exec

all: ../libabc.a dice_game dice.so $(addprefix abc_,$(DEMOS))

abc_%.json: reference.json %_partial.json
	gojq -s add $^ > $@

.FORCE:

../libabc.a: .FORCE
	$(MAKE) -s -C $(ABCDIR) -f Makefile

# this defines how to make the various abc_(sql|exec|dynamic) executables
abc_%: main_%.cpp examples.h ../libabc.a
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC)

# this defines how to make the dice_game executable for use with abc_exec 
dice_game: dice_game.cpp dice.h
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC) $(GSL_LIB)

# this defines how to make the dice.so shared library for use with abc_dynamic
dice.so: dice.cpp dice.h ../libabc.a
	$(CPP) $(SHAREDFLAGS) $< -o $@ $(ABC_INC)
	nm -gD $@ | grep simulator # check that the simulator is exported

# for the dynamic version, need to make the library
run_dynamic: | dice.so

# for the executable version, need to make the executable
run_exec: | dice_game

abc_%.sqlite: abc_% abc_%.json
	rm -f $@
	./$^ -b -p

# this provides a convenient way to run all of the abc_(sql|exec|dynamic) executables
run_%: abc_% abc_%.json | abc_%.sqlite
	./$^ -e -v

run_all: $(addprefix run_,$(DEMOS))

.PRECIOUS: abc_%.json abc_%.sqlite abc_%

clean:
	rm -f $(addprefix abc_,$(DEMOS)) dice_game dice.so
	rm -f --interactive=once $(patsubst %,abc_%.sqlite,$(DEMOS))

cleanall: clean
	$(MAKE) -C $(ABCDIR) clean