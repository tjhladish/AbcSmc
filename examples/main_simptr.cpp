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

//cout << "args0, par1: " << args[0] << " " << par1 << endl;
//cout << "args1, par2: " << args[1] << " " << par2 << endl;
    vector<double> results(par1,0);
    
    int sum = 0;

    for (int i = 0; i<par1; i++) {
        results[i] = gsl_rng_uniform_int(RNG, par2) + 1;
        sum += results[i];
    }

//cout << "Res: " << results[0] << " " << results[1] << endl;
    
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

    if (argc != 2) {
        cerr << "\n\tUsage: ./abc_simulator_pointer abc_config_simptr.json\n\n";
        return 100;
    }
    
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id

    AbcSmc* abc = new AbcSmc();
    abc->set_simulator(simulator);
    abc->parse_config(string(argv[1]));
    abc->run(RNG);

    return 0;
}
