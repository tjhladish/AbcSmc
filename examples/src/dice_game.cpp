#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <unistd.h>
#include "dice.h"

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "\n\tUsage: ./dice_game number_of_dice sides_on_die\n";
        cerr << "\tOutput: sum_of_faces stdev_of_faces\n\n";
        return 100;
    }

    const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid());
    double par1 = atof(argv[1]); // number of dice
    double par2 = atof(argv[2]); // number of sides on dice

    std::vector<double> results = simulator({ par1, par2 }, time (NULL) * getpid(), 0);
    
    cout << results[0] << " " << results[1] << endl;

    return 0;
}
