#ifndef ABCSMC_PARAMETER_H
#define ABCSMC_PARAMETER_H

#include <string>
#include <concepts> // let's us declare concepts for template constraints
#include <gsl/gsl_rng.h>
#include <math.h> // round
#include <AbcSmc/ParRNG.h>
#include <AbcSmc/TypeDefs.h>

// Design goals `Parameter`s:
//  - can be integral or floating point
//  - uniform interface for Prior, Posterior, Pseudo types
//  - "easy" to implement new concrete Priors
//  - has no state (i.e. any state managed by the ABCSMC object)
//  - has no knowledge of the ABCSMC object
//  - does not need to know about other parameters
//
// Requirements
// From other elements of AbcSmc, must expose:
// - construction
// - name, short_name, (for priors: also mean and sd)
// - sample, likelihood, noise, recast + valid
// - isPosterior
// - transformation?
//
//  challenges to accomplishing this:
//  - need to transform parameters, sometimes in terms of each other
//  - have transformations managed by the ABCSMC object?

namespace ABC { 
    class Parameter;
    typedef ParRNG<const Parameter, const gsl_rng> PRNG;
}

namespace ABC {
    
    class Parameter {
        public:
            Parameter(const std::string & s, const std::string & ss, const size_t & mi = 0) :
            name(s), short_name(ss), maxIdx(mi) {
                // TODO: sanitize short_name
            }

            // get the parameter name (for printing, etc - can be whatever format)
            std::string get_name() const { return name; };
            // get the parameter short name (for storage, etc - should be short and sanitized: no spaces, symbols other than underscores)
            std::string get_short_name() const { return short_name; };

            // these are the core methods that must be implemented by a Parameter

            // if this is an integral parameter, flatten it to the appropriate integer, then recast to double
            // overriding this should be managed by mixin `TParameter`
            virtual float_type recast(const float_type pval) const = 0;


            // sample from the parameter; this side-effects the PRNG *not* the parameter
            virtual float_type sample(PRNG & /* prng */) const = 0;
            // compute the likelihood of a value
            virtual float_type likelihood(const float_type /* pval */) const = 0;

            // can define `noise`, `get_mean`, `get_sd` effectively as errors, unless overriden

            // draw noise, repeating as necessary until resulting draw is valid
            virtual float_type noise(
                const gsl_rng * /* RNG */, const float_type /* mu */, const float_type /* sigma */,
                const size_t MAX_ATTEMPTS = 1000
            ) const {
                return std::numeric_limits<float_type>::signaling_NaN();
            };
            
            virtual float_type get_mean() const { return std::numeric_limits<float_type>::signaling_NaN(); };
            virtual float_type get_sd() const { return std::numeric_limits<float_type>::signaling_NaN(); };

            // some methods, there is a meaningful default, but we might wish to override it
            virtual bool isPosterior() const { return false; }; // the typical use case of a parameter is not a posterior
                
            // some computations can be done in terms of the properly defined methods
            bool valid(const float_type pval) const { return likelihood(pval) != 0.0; };

            size_t max_index() const { return maxIdx; };     // most parameters are not indexed (pseudo or posterior) parameters

        private:
            const std::string name;
            const std::string short_name;
            const size_t maxIdx; // usually unset

    };

    // This template will be used to define the int or floatness of `Parameter`s
    // used as roughly:
    // class SomePar : Parameter, TParameter<int, Parameter> { ... }
    template <NumericType NT, class B>
    struct TParameter : public virtual B {
        float_type recast(const float_type pval) const override {
            if constexpr (std::is_integral_v<NT>) {
                return static_cast<float_type>(std::round(pval));
            } else {
                return pval;
            }
        }
    };

}

#endif // ABCSMC_PARAMETER_H