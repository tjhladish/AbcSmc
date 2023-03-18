
#ifndef ABCSMC_PRIORS_H
#define ABCSMC_PRIORS_H

#include <AbcSmc/Parameter.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

namespace ABC {

// A Prior isa Parameter that does random sampling, and has meaningful mean, sd, noise, and likelihood
// still an Abstract Base Class however!
struct Prior : public Parameter {
    Prior(
        const std::string & s, const std::string & ss,
        const float_type mv, const float_type sv
    ) : Parameter(s, ss), meanval(mv), sdval(sv) {}

    float_type noise(
        const gsl_rng * rng, const float_type mu, const float_type sigma,
        const size_t MAX_ATTEMPTS = 1000
    ) const override {
        size_t attempts = 1;
        auto dev = trynoise(rng, mu, sigma);
        while(!valid(dev) and (attempts++ < MAX_ATTEMPTS)) { dev = trynoise(rng, mu, sigma); }
        if (!valid(dev)) { 
            std::cerr << "ERROR: failed to draw valid noise from prior " << get_name() << " - returning mean value." << std::endl;
            return get_mean();
        } else {
            return dev;
        }
        return dev;
    };

    float_type get_mean() const override { return meanval; }
    float_type get_sd() const override { return sdval; }

    protected:
        const float_type meanval, sdval;
        float_type trynoise(const gsl_rng * rng, const float_type mu, const float_type sigma) const {
            return recast(gsl_ran_gaussian(rng, sigma) + mu);
        };

};

struct GaussianPrior : public Prior {
    GaussianPrior(
        const std::string & nm, const std::string & snm,
        const float_type mn, const float_type sd
    ) : Prior(nm, snm, mn, sd) {}
    
    float_type sample(PRNG & prng) const override { return gsl_ran_gaussian(prng.rng(), sdval) + meanval; };

    float_type likelihood(const float_type pval) const override {
        return gsl_ran_gaussian_pdf(pval - meanval, sdval);
    };

    float_type recast(const float_type pval) const override { return pval; }
};

struct DiscreteUniformPrior : public Prior {
    DiscreteUniformPrior(
        const std::string & nm, const std::string & snm,
        const long min, const long max
    ) : Prior(
        nm, snm,
        static_cast<float_type>(max + min) / 2.0, static_cast<float_type>(max - min) / sqrt(12.0)
    ), minval(min), maxval(max) {
        assert(min < max);
    }

    float_type sample(PRNG & prng) const override {
        return static_cast<float_type>(gsl_rng_uniform_int(prng.rng(), maxval - minval + 1) + minval);
    };

    float_type likelihood(const float_type pval) const override {
        return ((pval == recast(pval)) and (minval <= pval) and (pval <= maxval)) ? 1.0 / (maxval - minval + 1) : 0.0;
    };

    float_type recast(const float_type pval) const override { return static_cast<float_type>(std::round(pval)); }

    private:
        const long minval, maxval;
};

struct ContinuousUniformPrior : public Prior {
    ContinuousUniformPrior(
        const std::string & nm, const std::string & snm,
        const float_type min, const float_type max
    ) : Prior(
        nm, snm,
        static_cast<float_type>(max + min) / 2.0,
        static_cast<float_type>(max - min) / sqrt(12.0)
    ), minval(min), maxval(max) {
        assert(min < max);
    }

    float_type sample(PRNG & prng) const override {
        return gsl_rng_uniform(prng.rng())*(maxval-minval) + minval;
    };

    float_type likelihood(const float_type pval) const override {
        return ((minval <= pval) and (pval <= maxval)) ? 1.0 / (maxval - minval) : 0.0;
    };

    float_type recast(const float_type pval) const override { return pval; }

    private:
        const float_type minval, maxval;
};

}

#endif // ABCSMC_PRIORS_H