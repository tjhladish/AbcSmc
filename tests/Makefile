
default: test_all

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

TESTS := abcutil pls

all: ../libabc.a $(TESTS)

%: %.cpp testing.h ../libabc.a
	$(CPP) $(CFLAGS) $< -o $@ $(ABC_INC)

.FORCE:

../libabc.a: .FORCE
	$(MAKE) -C $(ABCDIR) -f Makefile

test_all: $(patsubst %,test_%,$(TESTS))
	@echo tested everything
	rm $(TESTS)

test_%: %
	./$^

clean:
	rm -f $(TESTS)

cleanall: clean
	$(MAKE) -C $(ABCDIR) clean

prngtest:
	$(CPP) $(CFLAGS) -I../include ParRNG.test.cpp -o $@ -L../build -lgsl -lgslcblas -labc