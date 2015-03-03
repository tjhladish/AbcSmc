#include "AbcSmc.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <unistd.h>

using namespace std;
const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

// wrapper for simulator
// must take vector of doubles (ABC paramters) 
// and return vector of doubles (ABC metrics)
vector<long double> simulator(vector<long double> args, const MPI_par* mp) {
    int par1 = (int) args[0]; // number of dice
    int par2 = (int) args[1]; // number of sides on dice

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

    if (argc != 3) {
        cerr << "\n\tUsage: ./abc_sql abc_config_sql.json rank\n\n";
        return 100;
    }
    
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (atoi(argv[2]) == 0) {
        // if DB exists and contains a completed set, run pls and sample for next set:
        abc->process_database(RNG); // TODO
        // else:
        //abc->build_database(RNG);
    } else {
        abc->set_simulator(simulator);
        int smc_iteration = 0;
        abc->simulate_next_particles(10);
    }
    return 0;
}
