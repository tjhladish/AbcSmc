#ifndef ABCSMC_PARAMETER_H
#define ABCSMC_PARAMETER_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <concepts> // let's us declare concepts for template constraints
#include <AbcSmc/ParRNG.h>
#include <gsl/gsl_rng.h>
#include <cmath.h>

using std::cerr;
using std::endl;
using std::vector;
using std::map;
using std::is_integral_v;
using std::is_floating_point_v;

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

    // TODO relocate this
    // this defines shift/scale transformation.
    //
    // a = sum of transformed scale addends
    // b = prod of transformed scale factors
    // c = sum of untransformed scale addends
    // d = prod of untransformed scale factors
    // u = arbitrary untransform function
    // x = value on the transformed scale (i.e. the fitting scale)
    //
    // transform(x) = (u((x + a) * b) + c) * d (i.e. value on the model/natural scale)
    struct ParXform {

        ParXform(
            const float_type tp = 0.0, const float_type up = 0.0,
            const float_type tt = 1.0, const float_type ut = 1.0
        ) : tplus(tp), uplus(up), ttimes(tt), utimes(ut) {}

        template <typename T> 
        ParXform(const T & pre_shift, const T & post_shift, const T & pre_scale, const T & post_scale) :
        ParXform(
            std::accumulate(pre_shift.begin(), pre_shift.end(), 0.0),
            std::accumulate(post_shift.begin(), post_shift.end(), 0.0),
            std::accumulate(pre_scale.begin(), pre_scale.end(), 1.0, std::multiplies<float_type>()),
            std::accumulate(post_scale.begin(), post_scale.end(), 1.0, std::multiplies<float_type>())
        ) { }

        const float_type tplus, uplus, ttimes, utimes;

        float_type transform(const float_type & pval, float_type (*u)(const float_type &)) const {
            return (u((pval + tplus)*ttimes) + uplus)*utimes;
        }

    };
    
    class Parameter {
        public:
            Parameter(std::string s, std::string ss) : name(s), short_name(ss) {
                // TODO: sanitize short_name
            }

            // get the parameter name (for printing, etc - can be whatever format)
            std::string get_name() const { return name; };
            // get the parameter short name (for storage, etc - should be short and sanitized: no spaces, symbols other than underscores)
            std::string get_short_name() const { return short_name; };

            // these are the core methods that must be implemented by a Parameter

            // if this is an integral parameter, flatten it to the appropriate integer, then recast to double
            virtual float_type recast(const float_type pval) const = 0;
            // sample from the parameter; this side-effects the PRNG *not* the parameter
            virtual float_type sample(ParRNG<Parameter, const gsl_rng> & /* prng */) const = 0;
            // compute the likelihood of a value
            virtual float_type likelihood(const float_type /* pval */) const = 0;

            // can define `noise`, `get_mean`, `get_sd` effectively as errors, unless overriden

            // draw noise, repeating as necessary until resulting draw is valid
            virtual float_type noise(
                const gsl_rng* /* RNG */, const float_type /* mu */, const float_type /* sigma */,
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

        friend class ParRNG<Parameter, const gsl_rng>;

        protected:
            size_t max_index() const { return max_index; };     // most parameters are not indexed (pseudo or posterior) parameters
            const size_t max_index; // usually unset

        private:
            std::string name;
            std::string short_name;
    };

    // This template will be used to define the int or floatness of `Parameter`s
    // used as roughly:
    // class SomePar : Parameter, TParameter<int, Parameter> { ... }
    template <typename NT, typename B> requires (std::is_integral_v<NT> or std::is_floating_point_v<NT>)
    struct TParameter : virtual B {
        virtual float_type recast(const float_type pval) const override {
            if constexpr (std::is_integral_v<NT>) {
                return static_cast<float_type>(std::round(pval));
            } else {
                return pval;
            }
        }
    };

    typedef ParRNG<Parameter, const gsl_rng> PRNG;

}

#endif // ABCSMC_PARAMETER_H