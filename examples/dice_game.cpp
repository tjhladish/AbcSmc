#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>

using namespace std;

int main(int argc, char** argv) {
    if (argc != 3) {
        cerr << "\n\tUsage: ./dice_game number_of_dice number_of_sides\n\n";
        return 100;
    }
     const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(RNG, time (NULL) * getpid());
    int par1 = atoi(argv[1]); // number of dice
    int par2 = atoi(argv[2]); // number of sides on dice

    double *results = new double[par1];
    int sum = 0;

    for (int i = 0; i<par1; i++) {
        results[i] = gsl_rng_uniform_int(RNG, par2) + 1;
        sum += results[i];
    }

    cout << sum << " " << gsl_stats_sd(results, 1, par1) << endl;

    return 0;
}
