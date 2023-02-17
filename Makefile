-include local.mk

CPP := g++
CFLAGS := -O2 -Wall -std=c++17 --pedantic -fPIC
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

ABCSOURCES =  pls.cpp AbcUtil.cpp AbcSmc.cpp CCRC32.cpp AbcMPI.cpp
JSONSOURCES = $(patsubst %,$(JSONDIR)/src/%.cpp,json_reader json_value json_writer)
SQLSOURCES  = $(addprefix $(SQLDIR)/,sqdb.cpp sqlite3.c)

ABCOBJECTS  = $(ABCSOURCES:.cpp=.o) refsql.o
JSONOBJECTS = $(JSONSOURCES:.cpp=.o)
SQLOBJECTS  = $(SQLDIR)/sqdb.o $(SQLDIR)/sqlite3.o
ABC_HEADER = ./pls.h ./AbcUtil.h ./AbcSmc.h ./AbcSim.h ./AbcMPI.h ./AbcMPIPar.h ./refsql.h

default: libabc.a

ARCHIVE ?= $(AR) -rv

libabc.a: $(filter-out CCRC32.o,$(ABCOBJECTS)) $(SQLOBJECTS) $(JSONOBJECTS)
	$(ARCHIVE) $@ $^

$(SQLDIR)/%.o:
	$(MAKE) -C $(@D) $(@F)

$(JSONDIR)/%.o: $(JSONDIR)/%.cpp
	$(CPP) $(CFLAGS) -c $< -o $@

pls.o: pls.cpp pls.h
	$(CPP) $(CFLAGS) -c -I. $< -o $@

# https://stackoverflow.com/questions/4158900/embedding-resources-in-executable-using-gcc
# bin2c available from hxtools via apt
# TODO: strip comments from sql files prior to feeding to bin2c
refsql.h: sqlviews.sql sqldynamic.sql sqlutrans.sql
	bin2c -C $@ $^

refsql.c: refsql.h

refsql.o: refsql.c
	$(CPP) $(CFLAGS) -c -I. $< -o $@

%.o: %.cpp $(ABC_HEADER)
	$(CPP) $(LIBS) $(CFLAGS) -c $(INCLUDE) $< -o $@

clean:
	rm -f *.o *.a
	$(MAKE) -C $(SQLDIR) clean

test.sqlite: sqlviews.sql test.sql sqldynamic.sql
	$(foreach sql,$^,cat $(sql) | sqlite3 $@;)