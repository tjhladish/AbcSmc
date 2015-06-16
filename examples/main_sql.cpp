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
vector<long double> simulator(vector<long double> args, long unsigned int rng_seed, const MPI_par* mp) {
    gsl_rng_set(RNG, rng_seed);
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


void usage() {
    cerr << "\n\tUsage: ./abc_sql abc_config_sql.json --process\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --simulate -n <number of simulations per database write>\n\n";
    cerr << "\t       ./abc_sql abc_config_sql.json --process --simulate -n <number of simulations per database write>\n\n";

}


int main(int argc, char* argv[]) {

    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage();
        exit(100);
    }

    bool process_db = false;
    bool simulate_db = false;
    int buffer_size = -1;

    for (int i=2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) { 
            process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {  
            simulate_db = true;
            buffer_size = buffer_size == -1 ? 1 : buffer_size;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {  
            buffer_size = atoi(argv[++i]);
        } else {
            usage(); 
            exit(101);
        }
    }


    AbcSmc* abc = new AbcSmc();
    abc->parse_config(string(argv[1]));
    if (process_db) {
        gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id
        abc->process_database(RNG);
    } 

    if (simulate_db) {
        abc->set_simulator(simulator);
        abc->simulate_next_particles(buffer_size);
    }

    return 0;
}
