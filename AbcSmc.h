#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include <map>
#include "AbcUtil.h"
#include "sqdb.h"
#include <json/json.h>
#include "AbcSim.h"
#include "pls.h"

enum PriorType {UNIFORM, NORMAL, PSEUDO, POSTERIOR};
enum NumericType {INT, FLOAT};
enum VerboseType {QUIET, VERBOSE};
enum FilteringType {PLS_FILTERING, SIMPLE_FILTERING};

class Parameter {
    public:
        Parameter() {};

        Parameter( std::string s, std::string ss, PriorType p, NumericType n, double val1, double val2, double val_step, double (*u)(const double), std::pair<double, double> r, std::map< std::string, std::vector<int> > mm )
            : name(s), short_name(ss), ptype(p), ntype(n), step(val_step), untran_func(u), rescale(r), par_modification_map(mm) {
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

        double sample(const gsl_rng* RNG) {
            if (ptype == UNIFORM) {
                if (ntype == INT) {
                     // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                    return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
                } else {
                    return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
                }
            } else if (ptype == NORMAL) {
                if (ntype == INT) {
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

        // doubled variance of particles
        double get_doubled_variance(int t) const { return doubled_variance[t]; }
        void append_doubled_variance(double v2) { doubled_variance.push_back(v2); }
        void set_prior_limits(double min, double max) { fmin = min; fmax = max; }
        double get_prior_min() const { return fmin; }
        double get_prior_max() const { return fmax; }
        double get_prior_mean() const { return mean; }
        double get_prior_stdev() const { return stdev; }
        double get_state() const { return state; }
        double get_step() const { return step; }
        double increment_state() { return state += step; }
        double reset_state() { state = get_prior_min(); return state; }
        std::string get_name() const { return name; }
        std::string get_short_name() const { if (short_name == "") { return name; } else { return short_name; } }
        PriorType get_prior_type() const { return ptype; }
        NumericType get_numeric_type() const { return ntype; }
        //double untransform(const double t) const { return (rescale.second - rescale.first) * untran_func(t) + rescale.first; }
        std::map < std::string, std::vector<int> > get_par_modification_map() const { return par_modification_map; }
        double untransform(const double t, vector<double> pars) const {
            double new_t = t + pars[0];
            new_t *= pars[1];
            new_t = untran_func(new_t);
            new_t += pars[2];
            new_t *= pars[3];
//std::cerr << "t | mods | new_t | scaled_t : " << t << " | "
//          << pars[0] << " " << pars[1] << " " << pars[2] << " " << pars[3] << " | "
//          << new_t << " | " << (rescale.second - rescale.first) * new_t + rescale.first << endl;
        return (rescale.second - rescale.first) * new_t + rescale.first; }

    private:
        std::string name;
        std::string short_name;
        PriorType ptype;
        NumericType ntype;
        double fmin, fmax, mean, stdev, state, step;
        std::vector<double> doubled_variance;
        double (*untran_func) (const double);
        std::pair<double, double> rescale;
        std::map < std::string, std::vector<int> > par_modification_map; // how this par modifies others
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

/*
class Particle {
//enum ParticleStatus {UNDEFINED_PARAMETERS, UNDEFINED_METRICS, PARTICLE_COMPLETE};
    public:
        Particle() {
            status = UNDEFINED_PARAMETERS; serial = -1; posterior_rank = -1; weight = -1;
        }
        Particle(int s):serial(s) {
            status = UNDEFINED_PARAMETERS; posterior_rank = -1; weight = -1;
        }
        Particle(int s, std::vector<long double> p):serial(s), pars(p) {
            status = UNDEFINED_METRICS; posterior_rank = -1; weight = -1;
        }
        Particle(int s, std::vector<long double> p, vector<long double> m):serial(s), pars(p), mets(m) {
            status = PARTICLE_COMPLETE; posterior_rank = -1; weight = -1;
        }
        Particle(int s, std::vector<long double> p, vector<long double> m, int r, double w):serial(s), pars(p), mets(m), posterior_rank(r), weight(w) {
            status = PARTICLE_COMPLETE;
        }

        std::vector<long double> get_pars() const { return pars; }
        std::vector<long double> get_mets() const { return mets; }
        void set_pars(std::vector<long double> p) { assert(status==UNDEFINED_PARAMETERS); pars = p; status = UNDEFINED_METRICS; }
        void set_mets(std::vector<long double> m) { assert(status==UNDEFINED_METRICS);    mets = m; status = PARTICLE_COMPLETE; }
        ParticleStatus get_status() const { return status; }
        void set_posterior_rank(int r) { posterior_rank = r; }
        void set_weight(int w) { weight = w; }
        bool in_posterior() const { return posterior_rank >= 0; }
        int get_posterior_rank() const { return posterior_rank; }

    private:
        int serial;
        std::vector<long double> pars;
        std::vector<long double> mets;
        int posterior_rank;
        double weight;
        ParticleStatus status;
};


class ParticleSet {
//enum AbcStatus {INCOMPLETE_SET, TOO_FEW_SETS, ABC_COMPLETE};
//enum SetStatus {UNSAMPLED_PRIOR, INCOMPLETE_PARTICLES, UNDEFINED_POSTERIOR, SET_COMPLETE};
    public:
        ParticleSet() { status = UNSAMPLED_PRIOR; }
        SetStatus get_status() const { return status; }
        void set_status(SetStatus s) { status = s; }

    private:
        std::vector<Particle*> particles;
        AbcStatus status;
};*/


class AbcSmc {
    public:
        AbcSmc() {
            _num_smc_sets = 0;
            _pls_training_fraction = 0.5;
            //int _pls_training_set_size = 0;
            //_predictive_prior_fraction = 0.05;
            //_predictive_prior_size = 0;
            //_next_predictive_prior_size = 0;
            resume_flag = false;
            resume_directory = "";
            use_transformed_pars = false;
            _retain_posterior_rank = true;
            use_mvn_noise = false;
            use_pls_filtering = true;
            _mp = NULL;
        };

        AbcSmc( ABC::MPI_par &mp ) {
            _mp = &mp;
            resume_flag = false;
            resume_directory = "";
            use_transformed_pars = false;
            use_mvn_noise = false;
            use_pls_filtering = true;
        };

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        size_t get_smc_iterations() { return _num_smc_sets; }
        //void set_smc_set_sizes(vector<int> n) { _smc_set_sizes = n; }

        size_t get_num_particles(size_t set_num, VerboseType vt = VERBOSE) {
            int set_size = 0;
            if (set_num < _smc_set_sizes.size()) {
                set_size = _smc_set_sizes[set_num];
            } else {
                set_size = _smc_set_sizes.back();
                if (vt == VERBOSE) cerr << "Assuming set size of " << set_size << endl;
            }
            return set_size;
        }

        size_t get_pred_prior_size(size_t set_num, VerboseType vt = VERBOSE) {
            int set_size = 0;
            if (set_num < _predictive_prior_sizes.size()) {
                set_size = _predictive_prior_sizes[set_num];
            } else {
                set_size = _predictive_prior_sizes.back();
                if (vt == VERBOSE) cerr << "Assuming a predictive prior size of " << set_size << endl;
            }
            return set_size;
        }

/*        void set_predictive_prior_fraction(float f) {
            assert(f > 0);
            assert(f <= 1);
            _predictive_prior_fraction = f;
        }*/

        void set_next_predictive_prior_size(int set_idx, int set_size);

        void set_pls_validation_training_fraction(float f) {
            assert(f > 0);
            assert(f <= 1);
            _pls_training_fraction = f;
        }

        void set_simulation(AbcSimFun * abcsf) { _simulator = abcsf; }
        void set_executable(std::string cmd) { set_simulation(new AbcExec(cmd)); }
        void set_simulator(AbcSimF * simulator) { set_simulation(new AbcFPtr(simulator)); }
        void set_simulator(std::string soname) { set_simulation(new AbcFPtr(soname.c_str())); }

        void set_database_filename( std::string name ) { _database_filename = name; }
        void set_posterior_database_filename( std::string name ) { _posterior_database_filename = name; }
        void set_retain_posterior_rank( std::string retain_rank ) { _retain_posterior_rank = (retain_rank == "true"); }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );
        Metric* add_next_metric(std::string name, std::string short_name, NumericType ntype, double obs_val) {
            Metric* m = new Metric(name, short_name, ntype, obs_val);
            _model_mets.push_back(new Metric(name, short_name, ntype, obs_val));
            return m;
        }
        Parameter* add_next_parameter(std::string name, std::string short_name, PriorType ptype, NumericType ntype, double val1, double val2, double step,double (*u)(const double), std::pair<double, double> r, std::map<std::string, std::vector<int> > mm) {
            Parameter* p = new Parameter(name, short_name, ptype, ntype, val1, val2, step, u, r, mm);
            _model_pars.push_back(p);
            return p;
        }

        void set_filtering_type(FilteringType ft) {
            switch(ft) {
              case PLS_FILTERING:
                use_pls_filtering = true; break;
              case SIMPLE_FILTERING:
                use_pls_filtering = false; break;
              default:
                cerr << "ERROR: Unknown FilteringType: " << ft << endl; exit(-1);
            }
        }

        FilteringType get_filtering_type() const { return use_pls_filtering ? PLS_FILTERING : SIMPLE_FILTERING; }

        void process_predictive_prior_arguments(Json::Value par);
        bool parse_config(std::string conf_filename);
        void report_convergence_data(const size_t t);


        bool build_database(const gsl_rng* RNG);
        bool process_database(const gsl_rng* RNG);
//        bool read_SMC_set_from_database (int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);
        bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);

        bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        bool fetch_particle_parameters(
            sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss,
            vector<int> &serial, vector<ABC::Row> &par_mat, vector<unsigned long int> &seeds,
            const bool verbose = false
        );
        bool update_particle_metrics(sqdb::Db &db, vector<string> &update_metrics_strings, vector<string> &update_jobs_strings);

        bool simulate_particle_by_serial(const int serial_req);
        bool simulate_particle_by_posterior_idx(const int posterior_req);
        bool simulate_next_particles(const int n = 1, const int serial_req = -1, const int posterior_req = -1); // defaults to running next particle

        Parameter* get_parameter_by_name(string name) {
            Parameter* p = nullptr;
            for (auto _p: _model_pars) { if (_p->get_name() == name) p = _p;}
            assert(p);
            return p;
        }

        size_t npar() { return _model_pars.size(); }
        size_t nmet() { return _model_mets.size(); }

        PLS_Model run_PLS(ABC::Mat2D&, ABC::Mat2D&, const int pls_training_set_size, const int ncomp);
        std::string get_database_filename()                 { return _database_filename; }
        std::vector< ABC::Mat2D > get_particle_parameters() { return _particle_parameters; }
        std::vector< ABC::Mat2D > get_particle_metrics()    { return _particle_metrics; }

    private:
        ABC::Mat2D X_orig;
        ABC::Mat2D Y_orig;
        std::vector<Parameter*> _model_pars;
        PosteriorParameters _posterior_pars;

        std::vector<Metric*> _model_mets;
        size_t _num_smc_sets;
        vector<int> _smc_set_sizes;
        //int _num_particles;
        float _pls_training_fraction;
        //int _pls_training_set_size;
        vector<int> _predictive_prior_sizes;  // TODO -- at parsing time, pred prior fractions should be converted to sizes
        //int _next_predictive_prior_size;
        //int _predictive_prior_size; // number of particles that will be used to inform predictive prior
        AbcSimFun * _simulator = new AbcSimUnset();

        bool resume_flag;
        bool use_transformed_pars;
        std::string resume_directory;
        std::string _database_filename;
        std::string _posterior_database_filename;
        bool _retain_posterior_rank;
        std::vector< ABC::Mat2D > _particle_metrics;
        std::vector< ABC::Mat2D > _particle_parameters;
        std::vector< std::vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        std::vector< std::vector<double> > _weights;
        bool use_mvn_noise;
        bool use_pls_filtering;

        //mpi specific variables
        ABC::MPI_par *_mp;

        bool _run_simulator(ABC::Row &par, ABC::Row &met, const unsigned long int rng_seed, const unsigned long int serial);

        bool _populate_particles( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG );

        bool _populate_particles_mpi( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG );
        void _particle_scheduler(int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG);
        void _particle_worker();

        void _fp_helper (const int t, const ABC::Mat2D &X_orig, const ABC::Mat2D &Y_orig, const int next_pred_prior_size, const ABC::Col& distances);
        void _filter_particles_simple ( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, int pred_prior_size);
        PLS_Model _filter_particles ( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, int pred_prior_size);
        void _print_particle_table_header();
        long double calculate_nrmse(vector<ABC::Col> posterior_mets);

        void set_resume( bool res ) { resume_flag = res; }
        bool resume() { return resume_flag; }
        void set_resume_directory( std::string res_dir ) { resume_directory = res_dir; }
        bool read_particle_set( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig );
        bool read_predictive_prior( int t );

        string _build_sql_select_par_string(string tag);
        string _build_sql_select_met_string();
        string _build_sql_create_par_string(string tag);
        string _build_sql_create_met_string(string tag);

        bool _db_execute_stringstream(sqdb::Db &db, stringstream &ss);
        bool _db_execute_strings(sqdb::Db &db, std::vector<std::string> &update_buffer);
        bool _db_tables_exist(sqdb::Db &db, std::vector<string> table_names);

        bool _update_sets_table(sqdb::Db &db, int t);
        //bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);

        ABC::Col euclidean( ABC::Row obs_met, ABC::Mat2D sim_met );

        ABC::Mat2D slurp_posterior();

        ABC::Row sample_priors( const gsl_rng* RNG, ABC::Mat2D& posterior, int &posterior_rank );

        std::vector<double> do_complicated_untransformations( std::vector<Parameter*>& _model_pars, ABC::Row& pars );

        void calculate_doubled_variances( int t );

        void normalize_weights( std::vector<double>& weights );

        void calculate_predictive_prior_weights( int set_num );

        gsl_matrix* setup_mvn_sampler(const int);
        ABC::Row sample_mvn_predictive_priors( int set_num, const gsl_rng* RNG, gsl_matrix* L );

        ABC::Row sample_predictive_priors( int set_num, const gsl_rng* RNG );

        ABC::Row _z_transform_observed_metrics( ABC::Row& means, ABC::Row& stdevs );
};

#endif
