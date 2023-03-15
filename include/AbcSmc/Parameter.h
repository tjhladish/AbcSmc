#ifndef ABCSMC_PARAMETER_H
#define ABCSMC_PARAMETER_H

#include <vector>
#include <map>
#include <string>
#include <iostream>
#include <concepts> // let's us declare concepts for template constraints

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

    // this manages translating PriorTypes into particular instantiations
    Parameter * create_parameter(
        const std::string & s, const std::string & ss,
        const PriorType p,
        const double val1, const double val2, const double val_step,
        double (*u)(const double), std::pair<double, double> r,
        const map< std::string, vector<int> > & mm
    ) {
        return new ABC::TParameter<NT>(name, short_name, ptype, val1, val2, step, u, r, mm);
    };

    // this is a state-machine for use with both priors, posteriors, and pseudo parameters
    // should be passed around by reference
    // different parameter classes access this in different ways
    struct ParRNG {
        const gsl_rng* RNG;
        map<std::string, size_t> pseudo;
        bool lock_pseudo;
        size_t posterior;
    }

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
        template <typename T>
        ParXform(
            const T & pre_shift, const T & post_shift,
            const T & pre_scale, const T & post_scale
        ) : tplus(std::acumulate(pre_shift.begin(), pre_shift.end(), 0.0)),
            uplus(std::acumulate(post_shift.begin(), post_shift.end(), 0.0)),
            tplus(std::acumulate(pre_scale.begin(), pre_scale.end(), 1.0, std::multiplies<float_type>())),
            ttimes(std::acumulate(post_scale.begin(), post_scale.end(), 1.0, std::multiplies<float_type>()))
        { }

        ParXform() : tplus(0.0), uplus(0.0), ttimes(1.0), utimes(1.0) {}

        const float_type tplus, uplus, ttimes, utimes;

        float_type transform(const float_type & pval, float_type (*u)(const float_type &)) const {
            return (u((pval + tplus)*ttimes) + uplus)*utimes;
        }

    }

    class Parameter {
        public:
            Parameter(std::string s, std::string ss) : name(s), short_name(ss) {}

            std::string get_name() const { return name; };
            std::string get_short_name() const { return short_name; };

            virtual double sample(const gsl_rng* /* RNG */) const = 0;
            virtual double likelihood(const double /* pval */) const = 0;

            // can ignore defining `noise`, `get_mean`, `get_sd`, etc, if that kind of parameter never uses it
            virtual double noise(
                const gsl_rng* /* RNG */, const double /* mu */, const double /* sigma_squared */,
                const size_t MAX_ATTEMPTS = 1000
            ) const {
                return std::numeric_limits<double>::signaling_NaN();
            };
            virtual double get_mean() const { return std::numeric_limits<double>::signaling_NaN(); };
            virtual double get_sd() const { return std::numeric_limits<double>::signaling_NaN(); };

            // some methods, there is a typical, real default, but we might wish to override it
            virtual bool isPosterior() const { return false; };
            virtual bool increment_state() { return false; };
            // if *not* a transforming parameter, no-op
            virtual double untransform(
                const double pval, const std::pair<float_type, float_type> & /* rescale */, const ParXform & /* xform */
            ) const {
                return pval;
            }
            // if *not* an integer type parameter, this is a no-op
            virtual double recast(const double pval) const { return pval; };
    
            // some computations can be done in terms of the properly defined methods
            bool valid(const double pval) const { return likelihood(pval) != 0.0; };

        private:
            std::string name;
            std::string short_name;
    };

    template <typename NT> requires (std::is_integral_v<NT> or std::is_floating_point_v<NT>)
    class TParameter : public Parameter {
        public:
            TParameter(
                const std::string & s, const std::string & ss,
                const PriorType p,
                const double val1, const double val2, const double val_step,
                double (*u)(const double)
            ) : Parameter(s, ss), ptype(p), step(val_step), untran_func(u) {
                if (ptype == UNIFORM) {
                    assert(val1 < val2);
                    fmin = val1;
                    fmax = val2;
                    mean = (val2 + val1) / 2.0;
                    stdev = sqrt(pow(val2-val1,2)/12);
                    state = 0;   // dummy variable for UNIFORM
                } else if (ptype == NORMAL) {
                    fmin = std::numeric_limits<double>::lowest(); // NOT min()!!!! That's the smallest representable positive value.
                    fmax = std::numeric_limits<double>::max();
                    mean = val1;
                    stdev = val2;
                } else if (ptype == PSEUDO) {
                    fmin = val1;
                    fmax = val2;
                    mean = 0;    // dummy variable for PSEUDO
                    stdev = 0;   // dummy variable for PSEUDO
                    state = fmin;
                } else if (ptype == POSTERIOR) {
                    fmin = val1; // min index for posterior database, generally 0
                    fmax = val2; // max index for posterior database
                    mean = 0;    // dummy variable for PSEUDO
                    stdev = 0;   // dummy variable for PSEUDO
                    state = fmin;
                } else {
                    std::cerr << "Prior type " << ptype << " not supported.  Aborting." << std::endl;
                    exit(-200);
                }
            }

            double sample(const gsl_rng* RNG) override {
                if (ptype == UNIFORM) {
                    if constexpr (std::integral<NT>) {
                        // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                        return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
                    } else {
                        return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
                    }
                } else if (ptype == NORMAL) {
                    if constexpr (std::integral<NT>) {
                        cerr << "Integer type not supported for normal distributions.  Aborting." << endl;
                        exit(-199);
                    } else {
                        return gsl_ran_gaussian(RNG, stdev) + mean;
                    }
                } else {
                    std::cerr << "Prior type " << ptype << " not supported for random sampling.  Aborting." << std::endl;
                    exit(-200);
                }
            }

            double get_mean() const override { return mean; }
            double get_sd() const override { return stdev; }

            double increment_state() override { return state += step; }
            double reset_state() override { state = get_prior_min(); return state; }
            PriorType get_prior_type() const override { return ptype; }
            bool is_integral() const override { if constexpr (std::integral<NT>) { return true; } else { return false; } }
            //double untransform(const double t) const { return (rescale.second - rescale.first) * untran_func(t) + rescale.first; }
            map < std::string, vector<int> > get_par_modification_map() const override { return par_modification_map; }
            float_type untransform(
                const float_type pval, const std::pair<float_type, float_type> & rescale, const ParXform & xform = ParXform()
            ) const override {
                return (rescale.second - rescale.first) * xform.transform(pval, untran_func) + rescale.first;
            }
        }

        private:
            PriorType ptype;
            double fmin, fmax, mean, stdev, state, step;
            float_type (*untran_func) (const float_type &);
    };
}

#endif // ABCSMC_PARAMETER_H