-include local.mk

CPP:=g++
CFLAGS = -O2 -Wall -std=c++11 --pedantic -Wno-deprecated-declarations -Wno-ignored-attributes -Wno-int-in-bool-context -Wno-misleading-indentation
MKFILE_PATH := $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR = $(MKFILE_PATH)
SQLDIR  = $(ABCDIR)/sqdb

INCLUDE = -I. -I$(ABCDIR) -I$(ABCDIR)/jsoncpp/include -I$(SQLDIR)
ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

#LIBS = -lm -L$(TACC_GSL_LIB/) -L$(HPC_GSL_LIB/) -lgsl -lgslcblas
LIBS = -lm -lgsl -lgslcblas

SOURCES =  AbcSmc.cpp AbcUtil.cpp CCRC32.cpp
JSONDIR = $(ABCDIR)/jsoncpp/src
JSONSOURCES = $(JSONDIR)/json_reader.cpp $(JSONDIR)/json_value.cpp $(JSONDIR)/json_writer.cpp
SQLSOURCES  = $(SQLDIR)/sqdb.cpp

LIBABC  = libabc.a
LIBJSON = libjsoncpp.a
LIBSQL  = libsqdb.a

OBJECTS     = $(SOURCES:.cpp=.o)
JSONOBJECTS = $(JSONSOURCES:.cpp=.o)
SQLOBJECTS  = $(SQLSOURCES:.cpp=.o)
ABC_HEADER = ./pls.h ./AbcUtil.h ./AbcSmc.h

default: .all
.all:  $(LIBJSON) sqlite3.o $(LIBSQL) $(SOURCES) $(LIBABC)

sqlite3.o: $(SQLDIR)/sqlite3.c $(SQLDIR)/sqlite3.h
	gcc -c $(SQLDIR)/sqlite3.c -I$(SQLDIR)

ARCHIVE ?= $(AR) -rv

$(LIBABC): $(ABC_HEADER) $(OBJECTS) $(LIBSQL)
	$(ARCHIVE) $@ $(LIBSQL) $(OBJECTS)

$(LIBJSON): $(JSONOBJECTS)
	$(ARCHIVE) $@ $(JSONOBJECTS)

$(LIBSQL): $(SQLOBJECTS)
	$(ARCHIVE) $@ $(SQLOBJECTS)

%.o: %.cpp $(ABC_HEADER)
	$(CPP) $(LIBS) $(CFLAGS) -c $(INCLUDE) $< -o $@

clean:
	rm -f $(OBJECTS) $(JSONOBJECTS) $(SQLOBJECTS) $(LIBABC) $(LIBJSON) $(LIBSQL)
