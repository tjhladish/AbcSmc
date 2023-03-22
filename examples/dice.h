
#ifndef DICE_H
#define DICE_H

#include <gsl/gsl_statistics_double.h>
#include <vector>
#include "examples.h"
#include "AbcUtil.h"

using namespace std;

// wrapper for simulator
// must take vector of doubles (ABC paramters)
// and return vector of doubles (ABC metrics)
extern "C" std::vector<double> simulator(
    std::vector<double> args,
    const unsigned long int rng_seed,
    const unsigned long int serial,
    const ABC::MPI_par* /* mp */
) {
    
    gsl_rng_set(RNG, rng_seed); // seed the rng using the seed
    const size_t par1 = static_cast<size_t>(args[0]); // number of dice
    const size_t par2 = static_cast<size_t>(args[1]); // number of sides on dice

    vector<double> results(par1,0);

    int sum = 0;

    for (size_t i = 0; i < par1; i++) {
        results[i] = gsl_rng_uniform_int(RNG, par2) + 1;
        sum += results[i];
    }

    vector<double> metrics(2);
    metrics[0] = sum;
    if (par1 == 1) {
        metrics[1] = 0;
    } else {
        metrics[1] = gsl_stats_sd(&results[0], 1, par1);
    }

    return metrics;
};

#endif // DICE_H