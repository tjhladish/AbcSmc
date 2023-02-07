#ifndef ABCDB_H
#define ABCDB_H

#include <map>
#include "sqdb.h"

const string JOB_TABLE  = "job";
const string MET_TABLE  = "met";
const string PAR_TABLE  = "par";
const string UPAR_TABLE = "upar";

const string select_pars_limited =
    "SELECT P.* FROM " + PAR_TABLE +
    " JOIN (SELECT serial FROM " + JOB_TABLE +
    " WHERE status IN ('Q','R') ORDER BY status, attempts LIMIT ?1) USING(serial);";

const string reserve_jobs =
    "UPDATE " + JOB_TABLE + " SET startTime = ?1, status = 'R', attempts = attempts + 1 WHERE serial IN (?2);";



enum PriorType {UNIFORM, NORMAL, PSEUDO, POSTERIOR};
enum NumericType {INT, FLOAT};

// TODO: make this extend generic class?
// interface would have to preclude anything DB related
// then when creating an Abc instance, hand a file / other indicator
// that would store e.g. this or an AbcCSV or somesuch and use 
// fetch / update / check commands
class AbcDB {
    public:
        AbcDB(sqdb::Db &db) : _db(db) {
            // prepare all the bound statements
            fetch_particles = db.Query(select_pars_limited.c_str());
        }
        AbcDB(string dbfile) : AbcDB(sqdb::Db(dbfile.c_str())) { }
        bool fetch_parameters(const size_t n, const vector<Row> & par_mat, const bool verbose = false);
        bool update_metrics(
            const vector<Row> & met_mat,
            const vector<int> & serials,
            const vector<int> & starts,
            const vector<int> & durs,
            const bool verbose = false
        );
        // TODO these probably become separated steps - one to write/read in
        // a standardized way (in here) and all the calculation retained in
        // in current bits
        bool build_database(const gsl_rng* RNG);
        bool process_database(const gsl_rng* RNG);

    private:
        sqdb::Db &_db;
        sqdb::Statement initialize, fetch_particles, update_metrics;
        string update_values(
            const vector<Row> & met_mat,
            const vector<int> & serials,
            const vector<int> & starts,
            const vector<int> & durs
        ) {
            stringstream uv;
            
            for (int j = 0; j<nmet(); j++) { update_metrics << met_mat[i][j] << ", "; }
            update_metrics << start_time << ", " << time(0) - start_time << ")";
            update_metrics << (i == (par_mat.size() - 1) ? "" : ", ");

            auto i = 0;
            for (; i < (met_mat.size()-1); i++) {
                uv << "(" << serial[i] << ", ";
                for (auto j = 0; j < met_mat[i].size(); j++) { uv << met_mat[i][j] << ", "; }
                uv << starts[i] << ", " << durs[i] << "), ";
            }
            uv << "(" << serial[i] << ", ";
            for (auto j = 0; j < met_mat[i].size(); j++) { uv << met_mat[i][j] << ", "; }
            uv << starts[i] << ", " << durs[i] << ")";
            return uv.str();
        }

}

class AbcSmc {
    public:
        AbcSmc() { _mp = NULL; use_executable = false; use_simulator = false; resume_flag = false; resume_directory = ""; use_transformed_pars = false; use_mvn_noise = false; };
        AbcSmc( ABC::MPI_par &mp ) { _mp = &mp; use_executable = false; use_simulator = false; resume_flag = false; resume_directory = ""; use_transformed_pars = false; use_mvn_noise = false; };

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_predictive_prior_size(int n) { assert(n > 0); assert(n <= _num_particles); _predictive_prior_size = n; }
        void set_predictive_prior_fraction(float f)        { assert(f > 0); assert(f <= 1); _predictive_prior_size = _num_particles * f; }
        void set_pls_validation_training_fraction(float f) { assert(f > 0); assert(f <= 1); _pls_training_set_size = _num_particles * f; }

        // code smell - should be a one or the other thing
        // probably manage-able by providing a AbcCLI (Script? Exec?) class
        // that the executable into the to-from interface of SimFunc
        void set_executable( std::string name ) { _executable_filename = name; use_executable = true; }
        void set_simulator(SimFunc simulator) { _simulator = simulator; use_simulator = true; }


        // these get translated to generic template class handling / assignment
        void set_database_filename( std::string name ) { _database_filename = name; }
        void set_posterior_database_filename( std::string name ) { _posterior_database_filename = name; }
        void set_retain_posterior_rank( std::string retain_rank ) { _retain_posterior_rank = (retain_rank == "true"); }
        void write_particle_file( const int t );
        void write_predictive_prior_file( const int t );
        bool build_database(const gsl_rng* RNG);
        bool process_database(const gsl_rng* RNG);
        bool read_SMC_set_from_database (int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);
        bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        bool fetch_particle_parameters(
            vector<int> &serial, vector<ABC::Row> &par_mat, vector<unsigned long int> &seeds,
            bool verbose = false
        );
        bool update_particle_metrics(
            sqdb::Db &db, string results, const bool transactional = true, const bool verbose = false
        );



        Metric* add_next_metric(std::string name, std::string short_name, NumericType ntype, double obs_val) {
            Metric* m = new Metric(name, short_name, ntype, obs_val);
            _model_mets.push_back(m);
            return m;
        }
        Parameter* add_next_parameter(std::string name, std::string short_name, PriorType ptype, NumericType ntype, double val1, double val2, double step,double (*u)(const double), std::pair<double, double> r, std::map<std::string, std::vector<int> > mm) {
            Parameter* p = new Parameter(name, short_name, ptype, ntype, val1, val2, step, u, r, mm);
            _model_pars.push_back(p);
            return p;
        }

        // TODO this seems like a generating step, not a regularly referred to item?
        // e.g. all this info could be captured in DB
        bool parse_config(std::string conf_filename);


        void report_convergence_data(int);

        bool simulate_next_particles(const int n = 1, const bool verbose = false);

        Parameter* get_parameter_by_name(string name) {
            Parameter* p = nullptr;
            for (auto _p: _model_pars) { if (_p->get_name() == name) p = _p;}
            assert(p);
            return p;
        }

        int npar() { return _model_pars.size(); }
        int nmet() { return _model_mets.size(); }

    private:
        ABC::Mat2D X_orig;
        ABC::Mat2D Y_orig;
        std::vector<Parameter*> _model_pars;
        std::vector<Metric*> _model_mets;
        int _num_smc_sets;
        int _num_particles;
        int _pls_training_set_size;
        int _predictive_prior_size; // number of particles that will be used to inform predictive prior

// TODO refactor
        SimFunc _simulator;
        bool use_simulator;
        std::string _executable_filename;
        bool use_executable;
// TODO refactor

        bool resume_flag;
        bool use_transformed_pars;
        std::string resume_directory;
// TODO extract / refactor
        std::string _database_filename;
        std::string _posterior_database_filename;
// TODO extract / refactor
        bool _retain_posterior_rank;
        std::vector< ABC::Mat2D > _particle_metrics;
        std::vector< ABC::Mat2D > _particle_parameters;
        std::vector< std::vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        std::vector< std::vector<double> > _weights;
        bool use_mvn_noise;

        //mpi specific variables
        ABC::MPI_par *_mp;

        bool _run_simulator(ABC::Row &par, ABC::Row &met, const unsigned long int rng_seed, const unsigned long int serial);

        bool _populate_particles( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG );

        bool _populate_particles_mpi( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG );
        void _particle_scheduler(int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig, const gsl_rng* RNG);
        void _particle_worker();

        void _filter_particles ( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);
        void _print_particle_table_header();
        long double calculate_nrmse(vector<ABC::Col> posterior_mets);

        void set_resume( bool res ) { resume_flag = res; }
        bool resume() { return resume_flag; }
        void set_resume_directory( std::string res_dir ) { resume_directory = res_dir; }
        bool read_particle_set( int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig );
        bool read_predictive_prior( int t );

// TODO move all this
        string _build_sql_select_par_string(string tag);
        string _build_sql_select_met_string();
        string _build_sql_create_par_string(string tag);
        string _build_sql_create_met_string(string tag);
        bool _db_execute_stringstream(sqdb::Db &db, stringstream &ss);
        bool _db_execute_strings(sqdb::Db &db, std::vector<std::string> &update_buffer);
        bool _db_tables_exist(sqdb::Db &db, std::vector<string> table_names);
        bool _update_sets_table(sqdb::Db &db, int t);
        bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);
// TODO move all this

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
