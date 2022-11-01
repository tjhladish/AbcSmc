-include local.mk

CPP:=g++
CFLAGS = -O2 -Wall -std=c++17 --pedantic -fPIC
MKFILE_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR = $(MKFILE_PATH)
SQLDIR  = $(ABCDIR)/sqdb
JSONDIR = $(ABCDIR)/jsoncpp

INCLUDE = -I. -I$(ABCDIR) -I$(JSONDIR)/include -I$(SQLDIR)
ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

#LIBS = -lm -L$(TACC_GSL_LIB/) -L$(HPC_GSL_LIB/) -lgsl -lgslcblas
LIBS = -lm -lgsl -lgslcblas

ABCSOURCES =  AbcSmc.cpp AbcUtil.cpp CCRC32.cpp
SQLSOURCES  = $(addprefix $(SQLDIR)/,sqdb.cpp sqlite3.c)
JSONSOURCES = $(patsubst %,$(JSONDIR)/src/%.cpp,json_reader json_value json_writer)

ABCOBJECTS  = $(ABCSOURCES:.cpp=.o)
JSONOBJECTS = $(JSONSOURCES:.cpp=.o)
SQLOBJECTS  = $(SQLDIR)/sqdb.o $(SQLDIR)/sqlite3.o
ABC_HEADER = ./pls.h ./AbcUtil.h ./AbcSmc.h ./AbcSim.h

default: libabc.a

ARCHIVE ?= $(AR) -rv

libabc.a: $(filter-out CCRC32.o,$(ABCOBJECTS)) $(SQLOBJECTS) $(JSONOBJECTS)
	$(ARCHIVE) $@ $^

$(SQLDIR)/%.o:
	$(MAKE) -C $(@D) $(@F)

%.o: $(JSONDIR)/%.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

%.o: %.cpp $(ABC_HEADER)
	$(CPP) $(LIBS) $(CFLAGS) -c $(INCLUDE) $< -o $@

clean:
	rm -f *.o *.a
	$(MAKE) -C $(SQLDIR) clean
