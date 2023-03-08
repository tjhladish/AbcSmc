
// forward definition required because AbcSmc <=> AbcMPI cyclic dependency
#ifndef ABC_MPI_PAR_H
#define ABC_MPI_PAR_H

#ifdef USING_MPI
#include <mpi.h>
#endif // USING_MPI

// add MPI management to the ABC namespace
namespace ABC {
    
    struct MPI_par {
#ifdef USING_MPI
        MPI_Comm comm;
        MPI_Info info;
        int mpi_size, mpi_rank;
#else
        const static int mpi_size = 1;
        const static int mpi_rank = 0;
#endif // USING_MPI
    };

}

#endif // ABC_MPI_PAR_H