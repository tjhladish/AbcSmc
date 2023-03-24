
#ifndef DICE_H
#define DICE_H

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <vector>

using namespace std;

// wrapper for simulator
// must take vector of doubles (ABC parameters)
// and return vector of doubles (ABC metrics)
extern "C" std::vector<double> simulator(
    std::vector<double> parameters,
    const unsigned long int rng_seed,
    const unsigned long int serial
) {
    gsl_rng * RNG = gsl_rng_alloc(gsl_rng_taus2);
    // seed the rng using the provided seed
    gsl_rng_set(RNG, rng_seed);
    // NB: your simulator may use whatever rng you like, seed it however you like, etc

    // a typical first step is to assign and cast parameters input to simulation-relevant variables
    const size_t num_dice  = static_cast<size_t>(parameters[0]); // number of dice
    const size_t num_faces = static_cast<size_t>(parameters[1]); // number of sides on dice

    vector<double> results(num_dice, 0);

    int sum = 0;

    for (size_t i = 0; i < num_dice; i++) {
        results[i] = gsl_rng_uniform_int(RNG, num_faces) + 1;
        sum += results[i];
    }

    vector<double> metrics(2);
    metrics[0] = sum;
    if (num_dice == 1) {
        metrics[1] = 0;
    } else {
        metrics[1] = gsl_stats_sd(&results[0], 1, num_dice);
    }

    gsl_rng_free(RNG);

    return metrics;
};

#endif // DICE_H