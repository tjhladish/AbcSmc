
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

DEMOS := direct executable shared integral

default: $(patsubst %,demo/%,$(DEMOS)) $(patsubst %,%/config.json,$(DEMOS)) dice.so

bin output:
	mkdir -p $@

demo:
	cmake -S . -B $@

demo/%: demo
	cd $^ && $(MAKE) $*

%/config.json: reference.json %/partial.json
	gojq -s add $^ > $@

# this defines how to make the dice.so shared library for use in the `shared` demo
bin/dice.so: src/dice.cpp include/dice.h | bin
	$(CPP) $(SHAREDFLAGS) -Iinclude $< -o $@ $(LIBS)
	nm -gD $@ | grep simulator # check that the simulator is exported

# this defines how to make the dice_game executable for use in the `direct` and `executable` demos 
bin/dice_game: src/dice_game.cpp include/dice.h | bin
	$(CPP) $(CFLAGS) -Iinclude $< -o $@ $(LIBS)


# for the dynamic version, need to make the shared object
run_shared: | bin/dice.so

# for the direct and executable versions, need to make the executable
run_executable: | bin/dice_game

run_direct: | bin/dice_game

# in general, run off config => sqlite approach
run_%: demo/% %/config.json | output
	./$^ --process --simulate -n 1000

# except for the direct version
run_direct: demo/direct
	./$^ 10 44



clean:
	rm -f $(addprefix abc_,$(DEMOS)) dice_game dice.so
	rm -f --interactive=once $(patsubst abc_%.sqlite,$(DEMOS))

cleanall: clean
	$(MAKE) -C $(ABCDIR) clean