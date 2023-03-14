
# This Makefile demonstrates different options for using the AbcSmc library
# with different approaches to the simulator, illustrated with a toy simulator.
#
# The toy simulator represents rolling a number dice, all of the same face count, then returns their
# sum and standard deviation. The simulator is defined in dice.h, and follows the function signature
# required by AbcSmc.
#
# AbcSmc supports several ways to use a simulator:
#
# 1. The simulator is compiled into an executable along with the library (main_sql.cpp --> abc_sql).
# 2. The simulator is compiled into a shared library (dice.cpp / dice.so) and loaded dynamically (main_dynamic.cpp --> abc_dynamic).
# 3. The simulator is compiled into a separate executable (dice_game.cpp) and invoked externally (main_exec.cpp --> abc_exec).
#
# This Makefile demonstrates how to compile and run each of these options.

CPP := g++

CFLAGS := -O2 -Wall -std=c++17 --pedantic
SHAREDFLAGS := $(CFLAGS) -fPIC -shared -rdynamic -Wl,--no-as-needed
LIBS := -lm -lgsl -lgslcblas -ldl

default: run_shared

demo:
	cmake -S . -B demo

demo/%: demo
	cd $^ && $(MAKE) $*

%/config.json: reference.json %/partial.json
	gojq -s add $^ > $@

# this defines how to make the dice.so shared library for use in the `shared` demo
dice.so: src/dice.cpp include/dice.h
	$(CPP) $(SHAREDFLAGS) $< -o $@ $(LIBS)
	nm -gD $@ | grep simulator # check that the simulator is exported

# this defines how to make the dice_game executable for use in the `direct` and `executable` demos 
dice_game: src/dice_game.cpp include/dice.h
	$(CPP) $(CFLAGS) $< -o $@ $(LIBS)


# for the dynamic version, need to make the shared object
run_shared: | dice.so

run_%: demo/% %/config.json
	./$^ --process --simulate -n 1000

# make variable definitions

DEMOS := sql dynamic exec

all: ../libabc.a dice_game dice.so $(addprefix abc_,$(DEMOS))

abc_%.json: reference.json %_partial.json
	gojq -s add $^ > $@

.FORCE:

../libabc.a: .FORCE
	$(MAKE) -C $(ABCDIR) -f Makefile

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

# this provides a convenient way to run all of the abc_(sql|exec|dynamic) executables
run_%: abc_% abc_%.json
	./$^ --process --simulate -n 1000

run_all: $(addprefix run_,$(DEMOS))

.PRECIOUS: abc_%.json abc_%.sqlite abc_%

clean:
	rm -f $(addprefix abc_,$(DEMOS)) dice_game dice.so
	rm -f --interactive=once $(patsubst abc_%.sqlite,$(DEMOS))

cleanall: clean
	$(MAKE) -C $(ABCDIR) clean