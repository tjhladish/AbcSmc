
#ifndef ABC_MPI_H
#define ABC_MPI_H

#include "pls/pls.h" // for data types
#include "AbcSim.h" // for AbcSimFun type
#include "AbcMPIPar.h" // for the MPI_par struct

// add MPI management to the ABC namespace
namespace ABC {

    void particle_scheduler(
        Mat2D &X_orig, // metrics - should be properly sized in advance, to be filled in by this function
        const Mat2D &Y_orig, // parameters - should be filled in advance (e.g. by sample_priors)
        MPI_par *_mp
    );

    void particle_worker(
        const size_t npar, const size_t nmet,
        AbcSimFun* fptr,
        const unsigned long int seed,
        const unsigned long int serial,
        MPI_par *_mp
    );

}

#endif // ABC_MPI_H