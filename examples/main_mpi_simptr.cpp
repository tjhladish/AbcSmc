#include "AbcSmc.h"
#include "mpi.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <unistd.h>

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


// wrapper for simulator
// must take vector of doubles (ABC paramters) 
// and return vector of doubles (ABC metrics)
vector<long double> simulator(vector<long double> args, const MPI_par* mp) {
    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid());
    int par1 = (int) args[0]; // number of dice
    int par2 = (int) args[1]; // number of sides on dice

//cout << "args0, par1: " << args[0] << " " << par1 << endl;
//cout << "args1, par2: " << args[1] << " " << par2 << endl;
    vector<double> results(par1,0);
    
    int sum = 0;

    for (int i = 0; i<par1; i++) {
        results[i] = gsl_rng_uniform_int(RNG, par2) + 1;
        sum += results[i];
    }
    
    vector<long double> metrics(2);
    metrics[0] = sum;
    if (par1 == 1) {
        metrics[1] = 0;
    } else {
        metrics[1] = gsl_stats_sd(&results[0], 1, par1);
    }

    return metrics;
}

int main(int argc, char* argv[]) {
    MPI_par mp;
    setup_mpi(mp, argc, argv);

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc_mpi abc_config_file.json\n\n";
        return 100;
    }

    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc(mp);
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    abc->run(RNG);

    MPI_Finalize();
    return 0;
}
