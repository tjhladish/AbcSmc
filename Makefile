CC = g++
MPICC = mpicxx

CFLAGS = -O2
ABCDIR = $(HOME)/AbcSmc

#INCLUDE = -I/usr/include/eigen3/ -I$(ABCDIR)
INCLUDE = -I. -I$(ABCDIR) -I$(ABCDIR)/jsoncpp/include 
ifdef TACC_GSL_INC
INCLUDE += -I$$TACC_GSL_INC
endif
ifdef HPC_GSL_INC
INCLUDE += -I$$HPC_GSL_INC
endif

LIBS = -lm -L$(TACC_GSL_LIB/) -L$(HPC_GSL_LIB/) -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 

SOURCES =  AbcSmc.cpp AbcUtil.cpp CCRC32.cpp
JSONDIR = $(ABCDIR)/jsoncpp/src
JSONSOURCES = $(JSONDIR)/json_reader.cpp $(JSONDIR)/json_value.cpp $(JSONDIR)/json_writer.cpp

LIBABC = libabc.a
LIBJSON = libjsoncpp.a

OBJECTS = $(SOURCES:.cpp=.o)
JSONOBJECTS = $(JSONSOURCES:.cpp=.o)

default: all_no_mpi
.all:  $(LIBJSON) $(SOURCES) $(LIBABC)

all_no_mpi: CFLAGS += -Wall -std=c++11 --pedantic
all_no_mpi: .all

all_mpi: CC = $(MPICC)
all_mpi: CFLAGS += -w0 -std=c++11 -cxx=icc -D USING_MPI -D MPICH_IGNORE_CXX_SEEK -D MPICH_SKIP_MPICXX
all_mpi: .all

$(LIBABC): $(OBJECTS) 
	$(AR) -rv $(LIBABC) $(OBJECTS) 

$(LIBJSON): $(JSONOBJECTS)
	$(AR) -rv $(LIBJSON) $(JSONOBJECTS)

.cpp.o:
ifndef TACC_GSL_INC
ifndef HPC_GSL_INC
	@echo "Neither TACC_GSL_INC nor HPC_GSL_INC are defined. Do you need to run 'module load gsl'?"
endif
endif
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ 

clean:
	rm $(OBJECTS) $(JSONOBJECTS) $(LIBABC) $(LIBJSON)
