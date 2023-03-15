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

template<typename NT>
concept NumType = std::is_integral_v<NT> or std::is_floating_point_v<NT>;

namespace ABC {
    enum PriorType { UNIFORM, NORMAL, PSEUDO, POSTERIOR };

    class Parameter {
        public:
            virtual double sample(const gsl_rng* RNG) = 0;
            virtual void set_prior_limits(double min, double max) = 0;
            virtual double get_prior_min() const = 0;
            virtual double get_prior_max() const = 0;
            virtual double get_prior_mean() const = 0;
            virtual double get_prior_stdev() const = 0;
            virtual double get_state() const = 0;
            virtual double get_step() const = 0;
            virtual double increment_state() = 0;
            virtual double reset_state() = 0;
            virtual std::string get_name() const = 0;
            virtual std::string get_short_name() const = 0;
            virtual PriorType get_prior_type() const = 0;
            virtual bool is_integral() const = 0;
            virtual map < std::string, vector<int> > get_par_modification_map() const = 0;
            virtual double untransform(const double t, vector<double> pars) const = 0;
    };

    template <NumType NT>
    class TParameter : public Parameter {
        public:
            TParameter() {};

            TParameter(
                std::string s, std::string ss,
                PriorType p,
                double val1, double val2, double val_step,
                double (*u)(const double), std::pair<double, double> r,
                map< std::string, vector<int> > mm
            )
                : name(s), short_name(ss), ptype(p), step(val_step), untran_func(u), rescale(r), par_modification_map(mm) {
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

            void set_prior_limits(double min, double max) override { fmin = min; fmax = max; }
            double get_prior_min() const override { return fmin; }
            double get_prior_max() const override { return fmax; }
            double get_prior_mean() const override { return mean; }
            double get_prior_stdev() const override { return stdev; }
            double get_state() const override { return state; }
            double get_step() const override { return step; }
            double increment_state() override { return state += step; }
            double reset_state() override { state = get_prior_min(); return state; }
            std::string get_name() const override { return name; }
            std::string get_short_name() const override { if (short_name == "") { return name; } else { return short_name; } }
            PriorType get_prior_type() const override { return ptype; }
            bool is_integral() const override { if constexpr (std::integral<NT>) { return true; } else { return false; } }
            //double untransform(const double t) const { return (rescale.second - rescale.first) * untran_func(t) + rescale.first; }
            map < std::string, vector<int> > get_par_modification_map() const override { return par_modification_map; }
            double untransform(const double t, vector<double> pars) const override {
                double new_t = t + pars[0];
                new_t *= pars[1];
                new_t = untran_func(new_t);
                new_t += pars[2];
                new_t *= pars[3];
            return (rescale.second - rescale.first) * new_t + rescale.first;
        }

        private:
            std::string name;
            std::string short_name;
            PriorType ptype;
            double fmin, fmax, mean, stdev, state, step;
            double (*untran_func) (const double);
            std::pair<double, double> rescale;
            map < std::string, vector<int> > par_modification_map; // how this par modifies others
    };
}

#endif // ABCSMC_PARAMETER_H