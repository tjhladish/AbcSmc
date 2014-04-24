#include "AbcSmc.h"
#include "mpi.h"

using namespace std;

void setup_mpi(MPI_par &m, int &argc, char **argv) {
    /* MPI variables */
    m.comm  = MPI_COMM_WORLD;
    m.info  = MPI_INFO_NULL;

    /* Initialize MPI */
    MPI_Init(&argc, &argv);
    MPI_Comm_size(m.comm, &m.mpi_size);
    MPI_Comm_rank(m.comm, &m.mpi_rank);  
}

int main(int argc, char* argv[]) {
    MPI_par mp;
    setup_mpi(mp, argc, argv);

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc abc_config_file.json\n\n";
        return 100;
    }
    
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc();
    abc->use_MPI(mp);
    abc->parse_config(string(argv[1]));
    //abc->parse_config("abc.json");
    abc->run(RNG);

    MPI_Finalize();
    return 0;
}
