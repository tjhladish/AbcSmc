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

default: dice_game abc_sql abc_executable

$(ABCDIR)/libabc.a:
	$(MAKE) -C $(@D) -f Makefile $(@F)

abc_%: main_%.cpp $(ABCDIR)/libabc.a
	$(CPP) $(CFLAGS) $(ABC_INC) $< -o $@ $(ABC_LIB) $(GSL_LIB)

dice_game: dice_game.cpp
	$(CPP) $(CFLAGS) $^ -o $@ $(GSL_LIB)

clean:
	rm -f abc_sql abc_executable dice_game
