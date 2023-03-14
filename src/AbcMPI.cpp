
#include <AbcSmc/AbcMPI.h> // pre declarations for MPI functions / types

// TODO is there a way to capture this approach in AbcSim?
// As with those functors, this is largely a wrapper around a function pointer

namespace ABC {
    void particle_scheduler(
        Mat2D &X_orig, // metrics - should be properly sized in advance, to be filled in by this function
        const Mat2D &Y_orig, // parameters - should be filled in advance (e.g. by sample_priors)
        MPI_par *_mp
    ) {
#ifdef USING_MPI
        MPI_Status status;
        assert(X_orig.rows() == Y_orig.rows());

        auto *send_data = new long double[Y_orig.size()](); // all params, flattened array
        auto *rec_data  = new long double[X_orig.size()](); // all metrics, flattened array

        const size_t _num_particles = Y_orig.rows();
        const size_t npar = Y_orig.cols();
        const size_t nmet = X_orig.cols();

        // sample parameter distributions; copy values into Y matrix and into send_data buffer
        for (size_t i = 0; i < _num_particles; i++) {
            for (size_t j = 0; j < npar; j++) { send_data[i*npar + j] = Y_orig(i, j); }
        }
        // Seed the workers with the first 'num_workers' jobs
        int particle_id = 0;
        // Don't send particles that don't exist!
        for (int rank = 1; (rank < _mp->mpi_size) and (rank <= _num_particles); ++rank) {
            particle_id = rank - 1;                       // which row in Y
            MPI_Send(&send_data[particle_id*npar],  // message buffer
                    npar,                          // number of elements
                    MPI_LONG_DOUBLE,                 // data item is a double
                    rank,                            // destination process rank
                    particle_id,                     // message tag
                    _mp->comm);                      // always use this
        }

        // Receive a result from any worker and dispatch a new work request
        long double rec_buffer[nmet];
        particle_id++; // move cursor to next particle to be sent
        while ( particle_id < _num_particles ) {
            MPI_Recv(&rec_buffer,                     // message buffer
                    nmet,                          // message size
                    MPI_LONG_DOUBLE,                 // of type double
                    MPI_ANY_SOURCE,                  // receive from any sender
                    MPI_ANY_TAG,                     // any type of message
                    _mp->comm,                       // always use this
                    &status);                        // received message info
            for (int m = 0; m<nmet; ++m) rec_data[status.MPI_TAG*nmet + m] = rec_buffer[m];

            MPI_Send(&send_data[particle_id*npar],  // message buffer
                    npar,                          // number of elements
                    MPI_LONG_DOUBLE,                 // data item is a double
                    status.MPI_SOURCE,               // send it to the rank that just finished
                    particle_id,                     // message tag
                    _mp->comm);                      // always use this
            particle_id++; // move cursor
        }

        // receive results for outstanding work requests--there are exactly 'num_workers'
        // or '_num_particles' left, whichever is smaller
        for (int rank = 1; (rank < _mp->mpi_size) and (rank <= _num_particles); ++rank) {
            MPI_Recv(&rec_buffer,
                    nmet,
                    MPI_LONG_DOUBLE,
                    MPI_ANY_SOURCE,
                    MPI_ANY_TAG,
                    _mp->comm,
                    &status);
            for (int m = 0; m<nmet; ++m) rec_data[status.MPI_TAG*nmet + m] = rec_buffer[m];
        }

        // Tell all the workers they're done for now
        for (int rank = 1; rank < _mp->mpi_size; ++rank) {
            MPI_Send(0, 0, MPI_INT, rank, STOP_TAG, _mp->comm);
        }

        vector<int> bad_particle_idx; // bandaid, in case simulator returns nonsense values
        for (size_t i = 0; i < _num_particles; i++) {
            for (size_t j = 0; j < nmet; j++) {
                const double met_val = rec_data[i*nmet + j];
                if (not isfinite(met_val)) bad_particle_idx.push_back(i);
                X_orig(i, j) = rec_data[i*nmet + j];
            }
        }

        // bandaid, in case simulator returns nonsense values
        for (size_t i = 0; i < bad_particle_idx.size(); i++) {
            X_orig.row(i).fill(numeric_limits<float_type>::min());
            Y_orig.row(i).fill(numeric_limits<float_type>::min());
        }

        delete[] send_data;
        delete[] rec_data;
#endif
    };

    void particle_worker(
        const size_t npar, const size_t nmet,
        AbcSimFun* simulator,
        const unsigned long int seed,
        const unsigned long int serial,
        MPI_par *_mp
    ) {
#ifdef USING_MPI
        MPI_Status status;
        auto *local_Y_data = new long double[npar]();
        auto *local_X_data = new long double[nmet]();

        while (1) {
            MPI_Recv(local_Y_data, npar,
                    MPI_LONG_DOUBLE,
                    mpi_root,
                    MPI_ANY_TAG,
                    _mp->comm,
                    &status);

            // Check the tag of the received message.
            if (status.MPI_TAG == STOP_TAG) {
                delete[] local_Y_data;
                delete[] local_X_data;
                return;
            }

            Row pars(npar);
            Row mets(nmet);

            for (size_t j = 0; j < npar; j++) { pars[j] = local_Y_data[j]; }

            mets = simulator(pars, seed, serial, _mp);
    //        if (not success) exit(-210); // TODO handle throws / errors?

            for (size_t j = 0; j < nmet; j++) { local_X_data[j] = mets[j]; }

            MPI_Send(local_X_data, nmet,
                    MPI_LONG_DOUBLE,
                    mpi_root,
                    status.MPI_TAG,
                    _mp->comm);
        }
#endif
    };

}