
#ifndef ABCSMC_INDEXEDPARS_H
#define ABCSMC_INDEXEDPARS_H

#include <AbcSmc/Parameters.h>
#include <iostream>

namespace ABC {

// A Prior isa Parameter that does random sampling, and has meaningful mean, sd, noise, and likelihood
// still an Abstract Base Class however!
class IndexedPar : Parameter, public virtual TParameter<size_t, Parameter> {
    public:
        IndexedPar(
            const std::string & s, const std::string & ss,
            const size_t maxIdx
        ) : Parameter(s, ss), max_index(maxIdx) {}

    float_type likelihood(const float_type pval) const override {
        std::cerr << "ERROR: it is an error to ask for likelihood from an IndexedPar; attempted on " << get_name() << std::endl;
        exit(-1);
    }

};

class PseudoPar : IndexedPar {
    public:
        PseudoPar(
            const std::string & s, const std::string & ss,
            const std::vector<float_type> & vals
        ) : IndexedPar(s, ss, vals.size()-1), states(vals) {}

        // TODO support other constructors?

        float_type sample(PRNG & prng) const override { return states[rng.pseudo(this)]; }

    protected:
        const std::vector<float_type> states;

};

class PosteriorPar : IndexedPar {
    PseudoPar(
        const std::string & s, const std::string & ss,
        const size_t maxIdx
    ) : IndexedPar(s, ss, maxIdx) {}

    float_type sample(PRNG & rng) const override { return static_cast<float_type>(rng.posterior()); }
};

}

#endif // ABCSMC_PRIORS_H