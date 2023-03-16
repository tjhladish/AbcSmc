
#include <AbcSmc/Priors.h>

namespace ABC {

float_type GaussianPrior::sample(
    PRNG & prng
) const override {
    return gsl_ran_gaussian(prng.rng(), sdval) + meanval;
};

float_type GaussianPrior::likelihood(
    const float_type pval
) const override {
    return gsl_ran_gaussian_pdf(pval - meanval, sdval);
};

float_type DiscreteUniformPrior::sample(
    PRNG & prng
) const override {
    return static_cast<float_type>(gsl_ran_uniform(prng.rng(), maxval - minval + 1) + minval);
};

float_type DiscreteUniformPrior::likelihood(
    const float_type pval
) const override {
    return ((pval == recast(pval)) and (minval <= pval) and (pval <= maxval)) ? 1.0 / (maxval - minval + 1) : 0.0;
};

float_type ContinuousUniformPrior::sample(
    PRNG & prng
) const override {
    return gsl_rng_uniform(prng.rng())*(maxval-minval) + minval;
};

float_type ContinuousUniformPrior::likelihood(
    const float_type pval
) const override {
    return ((minval <= pval) and (pval <= maxval)) ? 1.0 / (maxval - minval) : 0.0;
};

}