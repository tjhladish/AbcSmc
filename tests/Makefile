
default: all

# This Makefile covers test cases for various features of the AbcSmc library.

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

TESTS := storage_sqlite

all: ../libabc.a $(addprefix run_,$(TESTS))

.FORCE:

../libabc.a: .FORCE
	$(MAKE) -C $(ABCDIR) -f Makefile

try_storage_sqlite: ../AbcStorage.h

# this defines how to make the various abc_(sql|exec|dynamic) executables
try_%: test_%.cpp ../libabc.a
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC)

.PRECIOUS: test_%

run_storage_sqlite: try_storage_sqlite
	./$^ test_storage.sqlite
	sqlite3 test_storage.sqlite ".schema --indent"
	sqlite3 test_storage.sqlite "select * from par;"
#	rm test_storage.sqlite

clean:
	rm -f $(addprefix try_,$(DEMOS))

cleanall: clean
	$(MAKE) -C $(ABCDIR) clean