#ifndef ABCSMC_H
#define ABCSMC_H

#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <float.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <math.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "pls.h"
#include "RunningStat.h"


using namespace std;

enum PriorType {UNIFORM, NORMAL};
enum NumericType {INT, FLOAT};

class Parameter {
    public:
        Parameter() {};

        Parameter( string s, PriorType p, NumericType n, double val1, double val2 ) : name(s), ptype(p), ntype(n) {
            assert(ptype == UNIFORM);
            if (ptype == UNIFORM) {
                fmin = val1;
                fmax = val2;
                mean = (val2 - val1) / 2.0;
                stdev = sqrt(pow(val2-val1,2)/12);
            } /*else if (ptype == NORMAL) {
                fmin = DBL_MIN;
                fmax = DBL_MAX;
                mean = val1;
                stdev = val2;
            }*/
        }

        double sample(const gsl_rng* RNG);/* { 
            if (ntype == INT) {
                 // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
            } else { 
                return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
            }
        }*/

        double get_doubled_variance(int t) { return doubled_variance[t]; }
        void append_doubled_variance(double v2) { doubled_variance.push_back(v2); }
        double get_min() { return fmin; }
        double get_max() { return fmax; }
        NumericType get_numeric_type() { return ntype; }

    private:
        string name;
        PriorType ptype;
        NumericType ntype;
        //int imin, imax;
        double fmin, fmax, mean, stdev;
        vector<double> doubled_variance;
};

class ModelParameters {
    public:
        ModelParameters() {};

        int size() { return _pars.size(); }

        void add_next_parameter(Parameter* p) { _pars.push_back(p); }

        void add_next_parameter(string name, PriorType ptype, NumericType ntype, double val1, double val2) {
            Parameter* p = new Parameter(name, ptype, ntype, val1, val2);
            _pars.push_back(p);
        }

        vector<Parameter*> get_parameters() { return _pars; }
        Parameter* get_parameter(int i) { return _pars[i]; }

    private:
        vector<Parameter*> _pars;
};


struct Metric {
    Metric() {};
    Metric(string n, NumericType nt, double val) : name(n), ntype(nt), obs_val(val) {};
    string name;
    NumericType ntype;
    double obs_val;
};


class ModelMetrics {
    public:
        void add_next_metric(string met, NumericType ntype, double obs_val) { _mets.push_back(new Metric(met, ntype, obs_val)); }
        vector<Metric*> get_metrics() { return _mets; }
        int size() { return _mets.size(); }

    private:
        vector<Metric*> _mets;
};


class AbcSmc {
    public:
        AbcSmc() {};
        AbcSmc(ModelParameters* pars, ModelMetrics* mets) { _model_pars = pars; _model_mets = mets; }

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_predictive_prior_size(int n) { assert(n > 0); assert(n <= _num_particles); _predictive_prior_size = n; }
        void set_predictive_prior_fraction(float f)        { assert(f > 0); assert(f <= 1); _predictive_prior_size = _num_particles * f; }
        void set_pls_validation_training_fraction(float f) { assert(f > 0); assert(f <= 1); _pls_training_set_size = _num_particles * f; }
        void set_metric_basefilename( string name ) { _metrics_filename = name; }
        //void set_executable_filename( string name ) { _executable_filename = name; }
        //void set_particle_basefilename( string name ) { _particle_filename = name; }
        //void set_predictive_prior_basefilename( string name ) { _predictive_prior_filename = name; }
        
        void run(string executable, const gsl_rng* RNG); 

        int npar() { return _model_pars->size(); }
        int nmet() { return _model_mets->size(); }

    private:
        Mat2D X_orig;
        Mat2D Y_orig;
        ModelParameters* _model_pars;
        ModelMetrics* _model_mets;
        int _num_smc_sets;
        int _num_particles;
        int _pls_training_set_size;
        int _predictive_prior_size; // number of particles that will be used to inform predictive prior
        string _executable_filename;
        string _metrics_filename;
        string _particle_filename;
        string _predictive_prior_filename;
        vector< Mat2D > _particle_metrics;
        vector< Mat2D > _particle_parameters;
        vector< vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        vector< vector<double> > _weights;

        void _populate_particles( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG ); 

        void _filter_particles ( int t, Mat2D &X_orig, Mat2D &Y_orig); 
        
        Col euclidean( Row obs_met, Mat2D sim_met ); 

        Row sample_priors( const gsl_rng* RNG );

        void calculate_doubled_variances( int t ); 

        void normalize_weights( vector<double>& weights ); 

        void calculate_predictive_prior_weights( int set_num ); 

        Row sample_predictive_priors( int set_num, const gsl_rng* RNG ); 

        Row _z_transform_observed_metrics( Row& means, Row& stdevs ); 
};

#endif
