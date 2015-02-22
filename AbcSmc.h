#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include "AbcUtil.h"
#include "sqdb.h"

enum PriorType {UNIFORM, NORMAL, PSEUDO};
enum NumericType {INT, FLOAT};

class Parameter {
    public:
        Parameter() {};

        Parameter( std::string s, std::string ss, PriorType p, NumericType n, double val1, double val2, double val_step ) 
            : name(s), short_name(ss), ptype(p), ntype(n), step(val_step) {
            if (ptype == UNIFORM) {
                fmin = val1;
                fmax = val2;
                mean = (val2 + val1) / 2.0;
                stdev = sqrt(pow(val2-val1,2)/12);
                state = 0; // dummy variable for UNIFORM
            } else if (ptype == PSEUDO) {
                fmin = val1;
                fmax = val2;
                mean = 0;  // dummy variable for PSEUDO 
                stdev = 0; // dummy variable for PSEUDO 
                state = fmin;
            } else {
                cerr << "Prior type " << ptype << " not supported.  Aborting." << endl;
                exit(-200);
            }

            /*else if (ptype == NORMAL) {
                fmin = DBL_MIN;
                fmax = DBL_MAX;
                mean = val1;
                stdev = val2;
            }*/
        }

        double sample(const gsl_rng* RNG) { 
            if (ptype == UNIFORM) {
                if (ntype == INT) {
                     // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                    return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
                } else { 
                    return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
                }
            } else {
                cerr << "Prior type " << ptype << " not supported for random sampling.  Aborting." << endl;
                exit(-200);
            }
        }

        // doubled variance of particles
        double get_doubled_variance(int t) const { return doubled_variance[t]; }
        void append_doubled_variance(double v2) { doubled_variance.push_back(v2); }
        double get_prior_min() const { return fmin; }
        double get_prior_max() const { return fmax; }
        double get_prior_mean() const { return mean; }
        double get_prior_stdev() const { return stdev; }
        double get_state() const { return state; }
        double increment_state() { return state += step; }
        double reset_state() { state = get_prior_min(); return state; }
        std::string get_name() const { return name; }
        std::string get_short_name() const { if (short_name == "") { return name; } else { return short_name; } }
        PriorType get_prior_type() const { return ptype; }
        NumericType get_numeric_type() const { return ntype; }

    private:
        std::string name;
        std::string short_name;
        PriorType ptype;
        NumericType ntype;
        double fmin, fmax, mean, stdev, state, step;
        std::vector<double> doubled_variance;
};

class Metric {
    public:
        Metric() {};
        Metric(std::string s, std::string ss, NumericType n, double val) : name(s), short_name(ss), ntype(n), obs_val(val) {};

        std::string get_name() const { return name; }
        std::string get_short_name() const { if (short_name == "") { return name; } else { return short_name; } }
        NumericType get_numeric_type() const { return ntype; }
        double get_obs_val() const { return obs_val; }

    private:
        std::string name;
        std::string short_name;
        NumericType ntype;
        double obs_val;
};


class AbcSmc {
    public:
        AbcSmc() { _mp = NULL; use_executable = false; use_simulator = false; resume_flag = false; resume_directory = ""; };
        AbcSmc( MPI_par &mp ) { _mp = &mp; use_executable = false; use_simulator = false; resume_flag = false; resume_directory = ""; };

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_predictive_prior_size(int n) { assert(n > 0); assert(n <= _num_particles); _predictive_prior_size = n; }
        void set_predictive_prior_fraction(float f)        { assert(f > 0); assert(f <= 1); _predictive_prior_size = _num_particles * f; }
        void set_pls_validation_training_fraction(float f) { assert(f > 0); assert(f <= 1); _pls_training_set_size = _num_particles * f; }
        //void set_metric_basefilename( std::string name ) { _metrics_filename = name; }
        void set_executable( std::string name ) { _executable_filename = name; use_executable = true; }
        void set_simulator(vector<float_type> (*simulator) (vector<float_type>, const MPI_par*)) { _simulator = simulator; use_simulator = true; }
        void set_particle_basefilename( std::string name ) { _particle_filename = name; }
        void set_predictive_prior_basefilename( std::string name ) { _predictive_prior_filename = name; }
        void set_database_filename( std::string name ) { _database_filename = name; }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );
        void add_next_metric(std::string name, std::string short_name, NumericType ntype, double obs_val) { 
            _model_mets.push_back(new Metric(name, short_name, ntype, obs_val)); 
        }
        void add_next_parameter(std::string name, std::string short_name, PriorType ptype, NumericType ntype, double val1, double val2, double step) {
            _model_pars.push_back(new Parameter(name, short_name, ptype, ntype, val1, val2, step));
        }
        
        bool parse_config(std::string conf_filename);
        void report_convergence_data(int);

        bool build_database(const gsl_rng* RNG);
        bool read_SMC_set_from_database (int t, Mat2D &X_orig, Mat2D &Y_orig);

        bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        bool do_sql(sqdb::Db &db, const char* sqlstring, const char* jobstring, Row &pars);
        bool do_sql(sqdb::Db &db, const char* sqlstring, const char* jobstring);
        bool simulate_database(const int smc_iteration, const int particle_id);

        void run(const gsl_rng* RNG) { run(_executable_filename, RNG); }; 
        void run(std::string executable, const gsl_rng* RNG); 
           
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
        vector<float_type> (*_simulator) (vector<float_type>, const MPI_par*);
        bool use_simulator;
        std::string _executable_filename;
        bool use_executable;
        bool resume_flag;
        std::string resume_directory;
        //std::string _metrics_filename;
        std::string _particle_filename;
        std::string _predictive_prior_filename;
        std::string _database_filename;
        std::vector< Mat2D > _particle_metrics;
        std::vector< Mat2D > _particle_parameters;
        std::vector< std::vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        std::vector< std::vector<double> > _weights;

        //mpi specific variables
        MPI_par *_mp;

        bool _run_simulator(Row &par, Row &met);

        bool _populate_particles( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG ); 

        bool _populate_particles_mpi( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG ); 
        void _particle_scheduler(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG);
        void _particle_worker();

        void _filter_particles ( int t, Mat2D &X_orig, Mat2D &Y_orig); 
        
        void set_resume( bool res ) { resume_flag = res; }
        bool resume() { return resume_flag; }
        void set_resume_directory( std::string res_dir ) { resume_directory = res_dir; }
        bool read_particle_set( int t, Mat2D &X_orig, Mat2D &Y_orig );
        bool read_predictive_prior( int t );

        string _build_sql_select_par_string();
        string _build_sql_select_met_string();
        string _build_sql_create_par_string(string tag);

        Col euclidean( Row obs_met, Mat2D sim_met ); 

        Row sample_priors( const gsl_rng* RNG );

        void calculate_doubled_variances( int t ); 

        void normalize_weights( std::vector<double>& weights ); 

        void calculate_predictive_prior_weights( int set_num ); 

        Row sample_predictive_priors( int set_num, const gsl_rng* RNG ); 

        Row _z_transform_observed_metrics( Row& means, Row& stdevs ); 
};

#endif
