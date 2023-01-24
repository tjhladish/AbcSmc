#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include <map>
#include "AbcUtil.h"
#include "sqdb.h"
#include <json/json.h>
#include <limits>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <cmath>

#include "AbcPar.h"

enum VerboseType { QUIET, VERBOSE };
enum FilteringType { PLS_FILTERING, SIMPLE_FILTERING };


class AbcSmc {
    public:
        AbcSmc() {
            _num_smc_sets = 0;
            _pls_training_fraction = 0.5;
            use_simulator = true;
            use_executable = false;
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
            use_executable = false;
            use_simulator = false;
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

        void set_executable( std::string name ) { _executable_filename = name; use_executable = true; }
        void set_simulator(vector<ABC::float_type> (*simulator) (vector<ABC::float_type>, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par*)) { _simulator = simulator; use_simulator = true; }
        void set_database_filename( std::string name ) { _database_filename = name; }
        void set_posterior_database_filename( std::string name ) { _posterior_database_filename = name; }
        void set_retain_posterior_rank( std::string retain_rank ) { _retain_posterior_rank = (retain_rank == "true"); }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );
        Metric* add_next_metric(Metric * m) {
            _model_mets.push_back(m);
            return m;
        }

        void add_next_parameter(Parameter * p) {
            _model_pars.push_back(p);
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
        void report_convergence_data(int);


        bool build_database(const gsl_rng* RNG);
        bool process_database(const gsl_rng* RNG);
        bool read_SMC_set_from_database (int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);

        bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        bool fetch_particle_parameters(sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss, vector<int> &serial, vector<ABC::Row> &par_mat, vector<unsigned long int> &seeds);
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

        // doubled variance of particles
        vector<double> get_doubled_variance(int t) const { return doubled_variance[t]; }

    private:
        vector<vector<double>> doubled_variance;
        ABC::Mat2D X_orig;
        ABC::Mat2D Y_orig;
        std::vector<Parameter*> _model_pars;
        std::vector<Metric*> _model_mets;
        int _num_smc_sets;
        vector<int> _smc_set_sizes;
        //int _num_particles;
        float _pls_training_fraction;
        //int _pls_training_set_size;
        vector<int> _predictive_prior_sizes;  // TODO -- at parsing time, pred prior fractions should be converted to sizes
        //int _next_predictive_prior_size;
        //int _predictive_prior_size; // number of particles that will be used to inform predictive prior
        vector<ABC::float_type> (*_simulator) (vector<ABC::float_type>, const unsigned long int rng_seed, const unsigned long int serial, const ABC::MPI_par*);
        bool use_simulator;
        std::string _executable_filename;
        bool use_executable;
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
        void _filter_particles ( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, int pred_prior_size);
        void _print_particle_table_header();
        long double calculate_nrmse(vector<ABC::Col> posterior_mets);

        void set_resume( bool res ) { resume_flag = res; }
        bool resume() { return resume_flag; }
        void set_resume_directory( std::string res_dir ) { resume_directory = res_dir; }
        bool read_particle_set( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig );
        bool read_predictive_prior( int t );

        string _build_sql_select_par_string() { return _build_sql(_model_pars, "P.", false); };
        string _build_sql_select_met_string() { return _build_sql(_model_mets, "M.", false); };
        string _build_sql_create_par_string() { return _build_sql(_model_pars, "", true); };
        string _build_sql_create_met_string() { return _build_sql(_model_mets, "", true); };

        bool _db_execute_stringstream(sqdb::Db &db, stringstream &ss);
        bool _db_execute_strings(sqdb::Db &db, std::vector<std::string> &update_buffer);
        bool _db_tables_exist(sqdb::Db &db, std::vector<string> table_names);

        bool _update_sets_table(sqdb::Db &db, int t);
        bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);

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

        template<typename ABCVAL>
        static string _build_sql(vector<ABCVAL*>& elements, const string tag, const bool creating);

};

template<typename ABCVAL>
string AbcSmc::_build_sql(vector<ABCVAL*>& elements, const string tag, const bool creating) {
    stringstream ss;
    // `tag` should be something like "P." or "par." - it's an SQL qualifier
    // so: if tag is non-empty, but doesn't end in `.`, add `.` to it
    for (size_t i = 0; i < elements.size() - 1; i++) { ss << tag << elements[i]->get_short_name() << (creating ? " real" : "") << ", "; }
    ss << tag << elements.back()->get_short_name() << (creating ? " real" : "");
    return ss.str();
}

#endif
