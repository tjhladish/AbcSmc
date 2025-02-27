#include <unistd.h>
#include <iostream>

#include <gsl/gsl_rng.h>
#include <AbcSmc/AbcSmc.h>

void usage() {
    cerr << "\n\tUsage: ./demo config.json\n\n";
}

int main(int argc, char* argv[]) {

    if (not (argc == 2) ) {
        usage();
        exit(100);
    }

    int buffer_size = 1;

    if ( argc == 4 ) {
        buffer_size = atoi(argv[4]);
    }

    AbcSmc* abc = new AbcSmc();
    // the config file here is setting:
    //  - the simulator shared object file (via "shared" key)
    //  - the database file (via "database_filename" key)
    abc->parse_config(string(argv[1]));

    // if the database file does not yet exist, create it
    unsigned long int rngseed = time(NULL) * getpid();
    gsl_rng * GSL_RNG = gsl_rng_alloc(gsl_rng_taus2); // RNG for AbcSmc
    gsl_rng_set(GSL_RNG, rngseed); // seed the rng using sys time and the process id
    abc->build_database(GSL_RNG);
    gsl_rng_free(GSL_RNG);

    abc->set_simulator(string(argv[2]));
    abc->simulate_next_particles(buffer_size);

    return 0;
}