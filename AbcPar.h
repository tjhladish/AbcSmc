
#ifndef ABCPAR_H
#define ABCPAR_H

#include <vector>
#include <limits>
#include <numeric>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <cmath>
#include <cassert>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "EnumMacros.h"

struct AbcVal {
    AbcVal(std::string nm, std::string snm = "") : name(nm), short_name(snm.empty() ? nm : snm) { };

    std::string get_name() const { return name; };
    std::string get_short_name() const { return short_name; };

    private:
        std::string name;
        std::string short_name;

};

#define NUMTYPE(M,SM) M(NumericType,SM(INT) SM(FLOAT))
CONSTRUCTENUM(NUMTYPE)

#define PARTYPE(M,SM) M(ParameterType,SM(UNIFORM) SM(GAUSSIAN) SM(NORMAL) SM(PSEUDO) SM(POSTERIOR))
CONSTRUCTENUM(PARTYPE)

#define XFORM(M,SM) M(TransformType,SM(NONE) SM(POW_10) SM(LOGISTIC))
CONSTRUCTENUM(XFORM)

struct Metric : public AbcVal {

    Metric(std::string s, std::string ss, NumericType n, double val) : AbcVal(s, ss), ntype(n), obs_val(val) {};

    NumericType get_numeric_type() const { return ntype; }
    double get_obs_val() const { return obs_val; }

    private:
        NumericType ntype;
        double obs_val;
};


// abstract base class for all Parameters
struct Parameter : public AbcVal {

    // TODO: goal of short name is provide a character limited version
    // worthwhile to trim full name when short name isn't provided?
    Parameter(std::string nm, std::string snm) : AbcVal(nm, snm) { };

    // all Parameters shall be named
    std::string get_name() const { return name; };
    std::string get_short_name() const { return short_name; };

    // must define `sample` & `likelidhood` for all concrete implementations
    virtual double sample(const gsl_rng* /* RNG */ ) const = 0;
    virtual double likelihood(const double /* pval */) const = 0;


    // can ignore defining `noise`, `get_mean`, `get_sd`, etc, if that kind of parameter never uses it
    virtual double noise(
        const double /* mu */, const double /* sigma_squared */, const gsl_rng* /* RNG */,
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
    virtual double untransform(const double pval) const { return pval; }
    // if *not* an integer type parameter, no-op
    virtual double recast(const double pval) const { return pval; };
    
    // some computations can be done in terms of the properly defined methods
    bool valid(const double pval) const { return likelihood(pval) != 0.0; };

    private:
        std::string name;
        std::string short_name;

};

// abstract base class for Priors, i.e. classes that have mean, sd, and are `noise`d in a consistent way
class Prior : public Parameter {
    public:
        Prior(std::string nm, std::string snm = "",
            double mv = std::numeric_limits<double>::signaling_NaN(),
            double sv = std::numeric_limits<double>::signaling_NaN()
        ) : Parameter(nm, snm), meanval(mv), sdval(sv) { };

        double get_mean() const override { return meanval; };
        double get_sd() const override { return sdval; };
    
        double noise(
            const double mu, const double sigma, const gsl_rng* RNG,
            const size_t MAX_ATTEMPTS = 1000
        ) const override {
            size_t attempts = 1;
            auto dev = trynoise(mu, sigma, RNG);
            while(!valid(dev) and (attempts++ < MAX_ATTEMPTS)) {
                dev = trynoise(mu, sigma, RNG);
            }
            if (!valid(dev)) { 
                exit(-300);
            }
            return dev;
        };

    protected:
        double meanval;
        double sdval;
        double trynoise(const double mu, const double sigma, const gsl_rng* RNG) const {
            return recast(gsl_ran_gaussian(RNG, sigma) + mu);
        };

};

class GaussianPrior : public Prior {

    public:
        GaussianPrior(std::string nm, std::string snm, double mn, double sd) :
            Prior(nm, snm, mn, sd) {}

        double sample(const gsl_rng* RNG ) const override {
            return gsl_ran_gaussian(RNG, sdval) + meanval;
        };

        double likelihood(const double pval) const override {
            return gsl_ran_gaussian_pdf(pval - meanval, sdval);
        };

};

template<typename NUMTYPE>
class UniformPrior : public Prior {
    public:
        UniformPrior(std::string nm, std::string snm, NUMTYPE mn, NUMTYPE mx) :
            Prior(nm, snm, (mx + mn) / 2.0, (mx - mn) / sqrt(12.0)), fmin(mn), fmax(mx) { 
                assert(mx >= mn);
            };

        double sample (const gsl_rng* RNG) const override {
            return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
        };

        double likelihood(const double pval) const override {
            return static_cast<double>((recast(pval) == pval) and (fmin <= pval) and (pval <= fmax));
        };

        double recast(const double pval) const override {
            if constexpr (std::is_integral_v<NUMTYPE>) {
                return std::round(pval);
            } else {
                return pval;
            }
        }

    private:
        NUMTYPE fmin, fmax;
};

template<typename NUMTYPE>
class PseudoParameter : public Parameter {
    public:
        PseudoParameter(
            std::string nm, std::string snm,
            std::vector<NUMTYPE> vals, bool post = false
        ) :
            Parameter(nm, snm), states(vals), posterior(post) {}
        
        double sample (const gsl_rng* /* RNG */ ) const override { return static_cast<double>(states[state]); };

        double likelihood(const double pval) const override {
            static_cast<double>(std::find(states.begin(), states.end(), static_cast<NUMTYPE>(pval)) != states.end());
        };

        bool increment_state() override {
            state++;
            if (state == states.size()) {
                state = 0;
                return false;
            } else {
                return true;
            }            
        }

        bool isPosterior() const override { return posterior; };

        double recast(const double pval) const override {
            if constexpr (std::is_integral_v<NUMTYPE>) {
                return std::round(pval);
            } else {
                return pval;
            }
        }

    private:
        size_t state = 0;
        std::vector<NUMTYPE> states;
        bool posterior;
};

typedef PseudoParameter<int> PseudoParameterInt;
typedef PseudoParameter<double> PseudoParameterDouble;

inline std::vector<int> ranks(const size_t to, const size_t from = 0) {
    assert(from < to);
    std::vector<int> res(to - from + 1);
    iota(res.begin(), res.end(), from);
    return res;
};

// TODO? posterior could be thinner, since it's always in reference to max rank => no need for vector
// but could be cleverer at the Pseudo par stage, and write in terms of begin/end iterators
// => that would allow for other container-flavored inputs, including `iota`s, which are probably lightweight
class Posterior : public PseudoParameterInt {
    public:
        Posterior(
            std::string nm, std::string snm,
            size_t maxrank, bool post = false
        ) :
            PseudoParameterInt(nm, snm, ranks(maxrank), true) {}

};

using namespace std::placeholders;

// TODO rewrite in Un-context
// wraps another Parameter in parameterized transformation
class TransformedParameter : public Parameter {
    public:
        TransformedParameter(
            Parameter * const p, double (*func) (const double, const std::vector<double> pars), const std::vector<double> ps
        ) : Parameter(p->get_name(), p->get_short_name()), pars(ps), xform(func), underlying(p) { };

        double sample(const gsl_rng* RNG ) const override { return underlying->sample(RNG); };
        double likelihood(double pval) const override { return underlying->likelihood(pval); };

        double noise(
            const double mu, const double sigma, const gsl_rng* RNG,
            const size_t MAX_ATTEMPTS = 1000
        ) const override {
            return underlying->noise(mu, sigma, RNG, MAX_ATTEMPTS);
        };

        double get_mean() const override { return underlying->get_mean(); };
        double get_sd() const override { return underlying->get_sd(); };

        bool isPosterior() const override { return underlying->isPosterior(); };
        bool increment_state() override { return underlying->increment_state(); };

        // by also calling `underlying`s untransform method, can compose transformations
        double untransform(const double pval) const { return xform(underlying->untransform(pval), pars); }
        // if *not* an integer type parameter, no-op
        double recast(const double pval) const { return underlying->recast(pval); };

    private:
        Parameter * const underlying;
        double (*xform) (const double, const std::vector<double> pars);
        std::vector<double> pars;

};

// provides common xform x => scale*(x+shift)
class LinearTransformedParameter : public TransformedParameter {
    public:
        LinearTransformedParameter(
            Parameter * const p,
            double stretch, double shift
        ) : TransformedParameter(
            p, [](double val, std::vector<double> ps) { return ps[0]*(val + ps[1]); },
            { stretch, shift }
        ) { };

};

// provides common xform x => (max - min)*val + min = 
// class RescaledParameter : public LinearTransformedParameter {
//     public:
//         RescaledParameter(
//             Parameter * const p,
//             double min, double max
//         ) : LinearTransformedParameter(
//             p, stretch, shift
//         ) { };

// };

#endif
