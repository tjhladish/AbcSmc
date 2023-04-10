
default: test_all

# make variable definitions
CPP := g++

CFLAGS := -O2 -Wall -std=c++17 --pedantic
SHAREDFLAGS := $(CFLAGS) -fPIC -shared -rdynamic -Wl,--no-as-needed

MKFILE_PATH = $(shell dirname $(realpath $(lastword $(MAKEFILE_LIST))))
ABCDIR := $(MKFILE_PATH)/..
SQLDIR := $(ABCDIR)/sqdb
JSNDIR := $(ABCDIR)/jsoncpp/include

ABC_INC := -I$(ABCDIR) -I$(SQLDIR) -I$(JSNDIR) -L$(ABCDIR) -labc $(GSL_LIB) -ldl

TESTS := ParRNG CLI

.FORCE:

%.test: %.test.cpp
	$(CPP) $(CFLAGS) -I../include $^ -o $@ -L../build -lgsl -lgslcblas -labc

test_all: $(addsuffix .test,$(TESTS))