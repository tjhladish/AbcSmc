#ifndef ABCSMC_PARAMETER_H
#define ABCSMC_PARAMETER_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <concepts> // let's us declare concepts for template constraints
#include <AbcSmc/ParRNG.h>
#include <gsl/gsl_rng.h>

using std::cerr;
using std::endl;
using std::vector;
using std::map;
using std::is_integral_v;
using std::is_floating_point_v;

// Design goals `Parameter`s:
//  - yields samples
//  - has no state (i.e. any state managed by the ABCSMC object)
//  - has no knowledge of the ABCSMC object
//  - does not need to know about other parameters
//
//  challenges to accomplishing this:
//  - need to transform parameters, sometimes in terms of each other
//  - sampling "posterior" or "pseudo" parameters requires state
//
//  way forward:
//  - have "sampling" a posterior / pseudo parameter increment an external state generator (analogous to the RNG for priors)?
//  - have transformations managed by the ABCSMC object?

namespace ABC {

    enum PriorType { UNIFORM, NORMAL, PSEUDO, POSTERIOR };

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
            Parameter(std::string s, std::string ss) : name(s), short_name(ss) {}

            std::string get_name() const { return name; };
            std::string get_short_name() const { return short_name; };

            virtual float_type sample(ParRNG<Parameter, const gsl_rng> & /* prng */) const = 0;
            virtual float_type likelihood(const float_type /* pval */) const = 0;

            // can define `noise`, `get_mean`, `get_sd` effectively as errors, unless overriden
            virtual float_type noise(
                const gsl_rng* /* RNG */, const float_type /* mu */, const float_type /* sigma_squared */,
                const size_t MAX_ATTEMPTS = 1000
            ) const {
                return std::numeric_limits<double>::signaling_NaN();
            };
            virtual float_type get_mean() const { return std::numeric_limits<float_type>::signaling_NaN(); };
            virtual float_type get_sd() const { return std::numeric_limits<float_type>::signaling_NaN(); };

            // some methods, there is a meaningful default, but we might wish to override it
            virtual bool isPosterior() const { return false; }; // the typical use case of a parameter is not a posterior
            virtual size_t max_index() const { return 0; };     // most parameters are not indexed (pseudo or posterior) parameters
            
            // if *not* a transforming parameter, no-op
            virtual float_type untransform(
                const float_type pval, const std::pair<float_type, float_type> & /* rescale */, const ParXform & /* xform */
            ) const {
                return pval;
            }
            // if *not* an integer type parameter, this is a no-op
            virtual float_type recast(const double pval) const { return pval; };
    
            // some computations can be done in terms of the properly defined methods
            bool valid(const float_type pval) const { return likelihood(pval) != 0.0; };

        private:
            std::string name;
            std::string short_name;
    };

    typedef ParRNG<Parameter, const gsl_rng> PRNG;

    template <typename NT> requires (std::is_integral_v<NT> or std::is_floating_point_v<NT>)
    class TParameter : public Parameter {
        public:
            TParameter(
                const std::string & s, const std::string & ss,
                const PriorType p,
                const float_type val1, const float_type val2,
                const float_type val_step,
                float_type (*u)(const float_type &)
            ) : Parameter(s, ss), ptype(p), step(val_step), untran_func(u) {
            }

            float_type sample(PRNG & prng) const override {
                if (ptype == UNIFORM) {
                    if constexpr (std::integral<NT>) {
                        // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                        return gsl_rng_uniform_int(prng.rng(), fmax-fmin + 1) + fmin;
                    } else {
                        return gsl_rng_uniform(prng.rng())*(fmax-fmin) + fmin;
                    }
                } else if (ptype == NORMAL) {
                    if constexpr (std::integral<NT>) {
                        cerr << "Integer type not supported for normal distributions.  Aborting." << endl;
                        exit(-199);
                    } else {
                        return gsl_ran_gaussian(prng.rng(), stdev) + mean;
                    }
                } else {
                    std::cerr << "Prior type " << ptype << " not supported for random sampling.  Aborting." << std::endl;
                    exit(-200);
                }
            }

            float_type likelihood(const float_type pval) const override {
                if (ptype == UNIFORM) {
                    if ((fmin <= pval) and (pval <= fmax)) {
                        return 1.0 / (fmax - fmin);
                    } else {
                        return 0.0;
                    }
                } else if (ptype == NORMAL) {
                    if constexpr (std::integral<NT>) {
                        cerr << "Integer type not supported for normal distributions.  Aborting." << endl;
                        exit(-199);
                    } else {
                        return gsl_ran_gaussian_pdf(pval - mean, stdev);
                    }
                } else {
                    std::cerr << "Prior type " << ptype << " not supported for random sampling.  Aborting." << std::endl;
                    exit(-200);
                }
            };

            // bool is_integral() const override { if constexpr (std::integral<NT>) { return true; } else { return false; } }
            
            float_type untransform(
                const float_type pval, const std::pair<float_type, float_type> & rescale, const ParXform & xform = ParXform()
            ) const override {
                return (rescale.second - rescale.first) * xform.transform(pval, untran_func) + rescale.first;
            }

        private:
            PriorType ptype;
            float_type fmin, fmax, mean, stdev, state, step;
            float_type (*untran_func) (const float_type &);
    };

    template <typename NT>
    class Prior : TParameter<NT> {
        public:
            Prior(
                const std::string & s, const std::string & ss,
                const PriorType p,
                const float_type val1, const float_type val2,
                const float_type val_step,
                float_type (*u)(const float_type &)
            ) : TParameter<NT>(s, ss, p, val1, val2, val_step, u), mean(val1), stdev(val2) {
            }

            float_type get_mean() const override { return mean; }
            float_type get_sd() const override { return stdev; }

        private:
            float_type mean, stdev;

    };

    template <typename NT>
    class Posterior : TParameter<NT> {
        public:
            Posterior(
                const std::string & s, const std::string & ss,
                const float_type val1, const float_type val2,
                const float_type val_step,
                float_type (*u)(const float_type &)
            ) : TParameter<NT>(s, ss, POSTERIOR, val1, val2, val_step, u) {}

    };

    template <typename NT>
    class Pseudo : TParameter<NT> {
        public:
            Pseudo(
                const std::string & s, const std::string & ss,
                const float_type val1, const float_type val2,
                const float_type val_step,
                float_type (*u)(const float_type &)
            ) : TParameter<NT>(s, ss, PSEUDO, val1, val2, val_step, u) {}

    };

    // this manages translating PriorTypes into particular instantiations
    template <typename NT>
    Parameter * create_parameter(
        const std::string & s, const std::string & ss,
        const PriorType p,
        const float_type val1, const float_type val2, const float_type val_step,
        float_type (*u)(const float_type &)
    ) {
        return new ABC::TParameter<NT>(s, ss, p, val1, val2, val_step, u);
    }
}

#endif // ABCSMC_PARAMETER_H