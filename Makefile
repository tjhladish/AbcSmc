CC = g++
MPICC = mpicxx

CFLAGS = -O2
ABCDIR = $(HOME)/work/AbcSmc

#INCLUDE = -I/usr/include/eigen3/ -I$(ABCDIR)
#INCLUDE = -I. -I$(ABCDIR) -I$$TACC_GSL_INC -I$(HOME)/jsoncpp/include 
INCLUDE = -I. -I$(ABCDIR) -I$(HPC_GSL_INC) -I$(ABCDIR)/jsoncpp/include 

#LIBS = -lm -L$$TACC_GSL_LIB/ -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 
LIBS = -lm -L$(HPC_GSL_LIB/) -lgsl -lgslcblas -L$(ABCDIR) -labc -ljsoncpp 

SOURCES =  AbcSmc.cpp utility.cpp 
JSONDIR = $(ABCDIR)/jsoncpp/src
JSONSOURCES = $(JSONDIR)/json_reader.cpp $(JSONDIR)/json_value.cpp $(JSONDIR)/json_writer.cpp

LIBABC = libabc.a
LIBJSON = libjsoncpp.a

OBJECTS = $(SOURCES:.cpp=.o)
JSONOBJECTS = $(JSONSOURCES:.cpp=.o)

default: all_no_mpi
.all:  $(LIBJSON) $(SOURCES) $(LIBABC)

all_no_mpi: CFLAGS += -Wall -std=c++0x --pedantic
all_no_mpi: .all

all_mpi: CC = $(MPICC)
all_mpi: CFLAGS += -w0 -std=c++0x -cxx=icc -D USING_MPI -D MPICH_IGNORE_CXX_SEEK -D MPICH_SKIP_MPICXX
all_mpi: .all

$(LIBABC): $(OBJECTS) 
	$(AR) -rv $(LIBABC) $(OBJECTS) 

$(LIBJSON): $(JSONOBJECTS)
	$(AR) -rv $(LIBJSON) $(JSONOBJECTS)

.cpp.o:
	$(CC) $(CFLAGS) -c $(INCLUDE) $< -o $@ 

clean:
	rm $(OBJECTS) $(JSONOBJECTS) $(LIBABC) $(LIBJSON)
