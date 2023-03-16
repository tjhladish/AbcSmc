
#ifndef ABCSMC_PRIORS_H
#define ABCSMC_PRIORS_H

#include <AbcSmc/Parameters.h>

namespace ABC {

// A Prior isa Parameter that does random sampling, and has meaningful mean, sd, noise, and likelihood
// still an Abstract Base Class however!
class Prior : Parameter {
    public:
        Prior(
            const std::string & s, const std::string & ss,
            const float_type mv, const float_type sv
        ) : Parameter(s, ss), meanval(mv), sdval(sv) {}

        float_type noise(
            PRNG & prng, const float_type mu, const float_type sigma,
            const size_t MAX_ATTEMPTS = 1000
        ) const override {
            size_t attempts = 1;
            auto dev = trynoise(prng, mu, sigma);
            while(!valid(dev) and (attempts++ < MAX_ATTEMPTS)) { dev = trynoise(prng, mu, sigma); }
            if (!valid(dev)) { 
                std::cerr << "ERROR: failed to draw valid noise from prior " << get_name() << " - returning mean value." std::endl;
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
        float_type trynoise(const gsl_rng * prng, const float_type mu, const float_type sigma) const {
            return recast(gsl_ran_gaussian(prng, sigma) + mu);
        };

};

class GaussianPrior : Prior, public virtual TParameter<float_type, Prior> {
    public:
        GaussianPrior(
            const std::string & nm, const std::string & snm,
            const float_type mn, const float_type sd
        ) : Prior(nm, snm, mn, sd) {}
};

class DiscreteUniformPrior : Prior, public virtual TParameter<long, Prior> {
    public:
        DiscreteUniformPrior(
            const std::string & nm, const std::string & snm,
            const long min, const long min
        ) : Prior(
            nm, snm,
            static_cast<float_type>(max + min) / 2.0, static_cast<float_type>(max - min) / sqrt(12.0)
        ), minval(min), maxval(max) {
            assert(min < max);
        }
    
    private:
        const long minval, maxval;
};

class ContinuousUniformPrior : Prior, public virtual TParameter<float_type, Prior> {
    public:
        ContinuousUniformPrior(
            const std::string & nm, const std::string & snm,
            const float_type min, const float_type min
        ) : Prior(
            nm, snm,
            static_cast<float_type>(max + min) / 2.0,
            static_cast<float_type>(max - min) / sqrt(12.0)
        ), minval(min), maxval(max) {
            assert(min < max);
        }
    private:
        const float_type minval, maxval;
};

}

#endif // ABCSMC_PRIORS_H