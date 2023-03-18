
#ifndef ABCSMC_INDEXEDPARS_H
#define ABCSMC_INDEXEDPARS_H

#include <AbcSmc/Parameter.h>
#include <iostream>

namespace ABC {

// A Prior isa Parameter that does random sampling, and has meaningful mean, sd, noise, and likelihood
// still an Abstract Base Class however!
struct IndexedPar : public Parameter {
    IndexedPar(
        const std::string & s, const std::string & ss,
        const size_t maxIdx
    ) : Parameter(s, ss, maxIdx) {}

    float_type likelihood(const float_type pval) const override {
        std::cerr << "ERROR: it is an error to ask for likelihood from an IndexedPar; attempted on " << get_name() << std::endl;
        exit(-1);
    }

    float_type recast(const float_type pval) const override {
        std::cerr << "ERROR: it is an error to attempt to recast an IndexedPar; attempted on " << get_name() << std::endl;
        exit(-1);
    }

};

struct PseudoPar : public IndexedPar {
    PseudoPar(
        const std::string & s, const std::string & ss,
        const std::vector<float_type> & vals
    ) : IndexedPar(s, ss, vals.size()-1), states(vals) {}

    // TODO support other constructors?

    float_type sample(PRNG & prng) const override { return states[prng.pseudo(this)]; }

    protected:
        const std::vector<float_type> states;

};

struct PosteriorPar : public IndexedPar {
    PosteriorPar(
        const std::string & s, const std::string & ss,
        const size_t maxIdx
    ) : IndexedPar(s, ss, maxIdx) {}

    float_type sample(PRNG & rng) const override { return static_cast<float_type>(rng.posterior()); }
    bool isPosterior() const override { return true; }
};

}

#endif // ABCSMC_INDEXEDPARS_H