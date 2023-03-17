#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include <map>
#include <AbcSmc/AbcUtil.h>
#include "sqdb.h"
#include <json/json.h>
#include <AbcSmc/AbcSim.h>
#include <PLS/pls.h>
#include <AbcSmc/Parameter.h>
#include <AbcSmc/Metric.h>

class AbcLog; // forward declaration of AbcLog; see AbcLog.h

enum FilteringType {PLS_FILTERING, SIMPLE_FILTERING};

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

        size_t get_num_particles(const size_t set_num, const bool verbose = true) {
            int set_size = 0;
            if (set_num < _smc_set_sizes.size()) {
                set_size = _smc_set_sizes[set_num];
            } else {
                set_size = _smc_set_sizes.back();
                if (verbose) cerr << "Assuming set size of " << set_size << endl;
            }
            return set_size;
        }

        size_t get_pred_prior_size(size_t set_num, const bool verbose = true) {
            int set_size = 0;
            if (set_num < _predictive_prior_sizes.size()) {
                set_size = _predictive_prior_sizes[set_num];
            } else {
                set_size = _predictive_prior_sizes.back();
                if (verbose) cerr << "Assuming a predictive prior size of " << set_size << endl;
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
        void set_simulator(AbcSimMPI * simulator) { set_simulation(new AbcFPtrMPI(simulator)); }
        void set_simulator(AbcSimBase * simulator) { set_simulation(new AbcFPtrBase(simulator)); }
        void set_simulator(std::string soname) { set_simulation(new AbcFPtrBase(soname.c_str())); }

        void set_database_filename( std::string name ) { _database_filename = name; }
        void set_posterior_database_filename( std::string name ) { _posterior_database_filename = name; }
        void set_retain_posterior_rank( std::string retain_rank ) { _retain_posterior_rank = (retain_rank == "true"); }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );

        const ABC::Metric * const add_next_metric(const ABC::Metric * const m) {
            _model_mets.push_back(m);
            _met_vals.push_back(m->get_obs_val());
            return m;
        }

        // TODO model_pars => map<string, Parameter*>, check for duplicate names
        const ABC::Parameter * const add_next_parameter(const ABC::Parameter * const p) {
            _model_pars.push_back(p);
            return true;
        }

        void add_modification_map(
            const ABC::Parameter * const par,
            map<std::string, std::vector<size_t> > mod_map
        ) { _par_modification_map[par] = mod_map; };

        void add_par_rescale(
            const ABC::Parameter * const par,
            std::pair<float_type, float_type> par_rescale
        ) { _par_rescale_map[par] = par_rescale; }
    
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

        bool build_database(const gsl_rng* RNG);
        bool process_database(const gsl_rng* RNG);
//        bool read_SMC_set_from_database (int t, Mat2D &X_orig, Mat2D &Y_orig);
        bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);

        bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        bool fetch_particle_parameters(
            sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss,
            vector<int> &serial, vector<Row> &par_mat, vector<unsigned long int> &seeds,
            const bool verbose = false
        );
        bool update_particle_metrics(sqdb::Db &db, vector<string> &update_metrics_strings, vector<string> &update_jobs_strings);

        bool simulate_next_particles(const int n = 1, const int serial_req = -1, const int posterior_req = -1); // defaults to running next particle
        bool simulate_particle_by_serial(const int serial_req) { return simulate_next_particles(1, serial_req, -1); }
        bool simulate_particle_by_posterior_idx(const int posterior_req) { return simulate_next_particles(1, -1, posterior_req); }

        ABC::Parameter* get_parameter_by_name(string name) {
            ABC::Parameter* p = nullptr;
            for (auto _p: _model_pars) { if (_p->get_name() == name) p = _p;}
            assert(p);
            return p;
        }

        size_t npar() { return _model_pars.size(); }
        size_t nmet() { return _model_mets.size(); }

        std::string get_database_filename()                 { return _database_filename; }
        std::vector< Mat2D > get_particle_parameters() { return _particle_parameters; }
        std::vector< Mat2D > get_particle_metrics()    { return _particle_metrics; }

        Row get_doubled_variance(const size_t t) const { return _doubled_variance[t]; }
        void append_doubled_variance(const Row & v2) { _doubled_variance.push_back(v2); }

        map<std::string, std::vector<size_t> > get_par_modification_map(const size_t parIdx) const { return _par_modification_map.at(parIdx); }
        std::pair<float_type,float_type> get_par_rescale(const size_t parIdx) const { return _par_rescale_map.at(parIdx); }

    private:
        friend AbcLog;
        Mat2D X_orig;
        Mat2D Y_orig;
        std::vector<const ABC::Parameter*> _model_pars;
        std::vector<const ABC::Metric*> _model_mets;
        map<const ABC::Parameter*, map<std::string, std::vector<size_t> > > _par_modification_map;
        map<const ABC::Parameter*, std::pair<float_type, float_type> > _par_rescale_map; 

        Row _met_vals;
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
        std::vector< Mat2D > _particle_metrics;
        std::vector< Mat2D > _particle_parameters;
        std::vector< std::vector<size_t> > _predictive_prior; // vector of row indices for particle metrics and parameters
        std::vector<Col> _weights;
        std::vector<Row> _doubled_variance;

        bool use_mvn_noise;
        bool use_pls_filtering;

        //mpi specific variables
        ABC::MPI_par *_mp;

        bool _run_simulator(Row &par, Row &met, const unsigned long int rng_seed, const unsigned long int serial);

        bool _populate_particles( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG );

//      TODO: undefined?
//        bool _populate_particles_mpi( int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG );
        void _particle_scheduler_mpi(const size_t t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG);
        void _particle_worker_mpi(const size_t seed, const size_t serial);

        void _set_predictive_prior(const int t, const int next_pred_prior_size, const Col& distances);

        void _filter_particles_simple(int t, Mat2D &X_orig, Mat2D &Y_orig, int pred_prior_size);
        PLS_Model _filter_particles(
            const size_t t, const Mat2D & X_orig, const Mat2D & Y_orig,
            const size_t pred_prior_size, const bool verbose = true
        );

        long double calculate_nrmse(vector<Col> posterior_mets);

        void set_resume(const bool res) { resume_flag = res; }
        bool resume() { return resume_flag; }
        void set_resume_directory( std::string res_dir ) { resume_directory = res_dir; }
        bool read_particle_set( int t, Mat2D &X_orig, Mat2D &Y_orig );
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

        Mat2D slurp_posterior();

        std::vector<double> do_complicated_untransformations(const std::vector<ABC::Parameter*> & _model_pars, const Row & pars );

        void calculate_doubled_variances( const size_t previous_set );

        void calculate_predictive_prior_weights( int set_num );

};

#endif
