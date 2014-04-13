#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include "utility.h"

enum PriorType {UNIFORM, NORMAL};
enum NumericType {INT, FLOAT};

class Parameter {
    public:
        Parameter() {};

        Parameter( std::string s, PriorType p, NumericType n, double val1, double val2 ) : name(s), ptype(p), ntype(n) {
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

        double sample(const gsl_rng* RNG) { 
            if (ntype == INT) {
                 // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
            } else { 
                return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
            }
        }

        double get_doubled_variance(int t) { return doubled_variance[t]; }
        void append_doubled_variance(double v2) { doubled_variance.push_back(v2); }
        double get_min() { return fmin; }
        double get_max() { return fmax; }
        std::string get_name() {return name; }
        NumericType get_numeric_type() { return ntype; }

    private:
        std::string name;
        PriorType ptype;
        NumericType ntype;
        //int imin, imax;
        double fmin, fmax, mean, stdev;
        std::vector<double> doubled_variance;
};

struct Metric {
    Metric() {};
    Metric(std::string n, NumericType nt, double val) : name(n), ntype(nt), obs_val(val) {};
    std::string name;
    NumericType ntype;
    double obs_val;
};


class AbcSmc {
    public:
        AbcSmc() { useMPI = false; _mp = NULL; };
        //AbcSmc(ModelParameters* pars, ModelMetrics* mets) { _model_pars = pars; _model_mets = mets; }

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_predictive_prior_size(int n) { assert(n > 0); assert(n <= _num_particles); _predictive_prior_size = n; }
        void set_predictive_prior_fraction(float f)        { assert(f > 0); assert(f <= 1); _predictive_prior_size = _num_particles * f; }
        void set_pls_validation_training_fraction(float f) { assert(f > 0); assert(f <= 1); _pls_training_set_size = _num_particles * f; }
        //void set_metric_basefilename( std::string name ) { _metrics_filename = name; }
        void set_executable( std::string name ) { _executable_filename = name; }
        void set_particle_basefilename( std::string name ) { _particle_filename = name; }
        void set_predictive_prior_basefilename( std::string name ) { _predictive_prior_filename = name; }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );
        void add_next_metric(std::string name, NumericType ntype, double obs_val) { 
            _model_mets.push_back(new Metric(name, ntype, obs_val)); 
        }
        void add_next_parameter(std::string name, PriorType ptype, NumericType ntype, double val1, double val2) {
            _model_pars.push_back(new Parameter(name, ptype, ntype, val1, val2));
        }
        
        bool parse_config(std::string conf_filename);
        void report_convergence_data(int);

        void run(const gsl_rng* RNG) { run(_executable_filename, RNG); }; 
        void run(std::string executable, const gsl_rng* RNG); 
           
        void use_MPI( MPI_par& mp ) { bool useMPI=true; _mp = &mp; } 
        int npar() { return _model_pars.size(); }
        int nmet() { return _model_mets.size(); }

    private:
        Mat2D X_orig;
        Mat2D Y_orig;
        std::vector<Parameter*> _model_pars;
        std::vector<Metric*> _model_mets;
        int _num_smc_sets;
        int _num_particles;
        int _pls_training_set_size;
        int _predictive_prior_size; // number of particles that will be used to inform predictive prior
        std::string _executable_filename;
        //std::string _metrics_filename;
        std::string _particle_filename;
        std::string _predictive_prior_filename;
        std::vector< Mat2D > _particle_metrics;
        std::vector< Mat2D > _particle_parameters;
        std::vector< std::vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        std::vector< std::vector<double> > _weights;

        //mpi specific variables
        bool useMPI;
        MPI_par *_mp;

        bool _populate_particles( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG ); 

        void _filter_particles ( int t, Mat2D &X_orig, Mat2D &Y_orig); 
        
        Col euclidean( Row obs_met, Mat2D sim_met ); 

        Row sample_priors( const gsl_rng* RNG );

        void calculate_doubled_variances( int t ); 

        void normalize_weights( std::vector<double>& weights ); 

        void calculate_predictive_prior_weights( int set_num ); 

        Row sample_predictive_priors( int set_num, const gsl_rng* RNG ); 

        Row _z_transform_observed_metrics( Row& means, Row& stdevs ); 
};

#endif
