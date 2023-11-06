#ifndef ABCSMC_H
#define ABCSMC_H

#define mpi_root 0

#include <map>
#include <optional>
#include <vector>

#include <AbcSmc/AbcUtil.h>
#include <AbcSmc/AbcDB.h> // provides AbcDB object
#include <AbcSmc/AbcSim.h>
#include <PLS/types.h>
#include <AbcSmc/Parameter.h>
#include <AbcSmc/Metric.h>
#include <AbcSmc/ParXform.h>

namespace ABC {
    enum FILTER { PLS, SIMPLE };
    enum NOISE { INDEPENDENT, MULTIVARIATE };
    class AbcLog; // forward declaration of AbcLog; see AbcLog.h

// An `AbcSmc` manages the interaction of the several stateful components
// of an ABC SMC PLS analysis.
//
// Conventions:
//  - internal state fields: _field_name
//  - private methods: _method_name()
//  - public fields: field_name
//  - public methods: method_name()
class AbcSmc {
    public:
        // default constructor - should be appropriate for `projection` mode
        AbcSmc();
        // default destructor - has to be defined where Eigen is available
        ~AbcSmc();

        size_t get_smc_iterations() const { return _num_smc_sets; }

        size_t get_smc_size_at(const size_t set_num) const {
            if (set_num >= _num_smc_sets) throw std::out_of_range("set_num out of range");
            return (set_num < _smc_set_sizes.size()) ? _smc_set_sizes[set_num] : _smc_set_sizes.back();
        }

        size_t get_pred_prior_size_at(const size_t set_num) {
            if (set_num >= _num_smc_sets) throw std::out_of_range("set_num out of range");
            return (set_num < _predictive_prior_sizes.size()) ? _predictive_prior_sizes[set_num] : _predictive_prior_sizes.back();
        }

        // should be private
        void set_smc_iterations(const size_t n) { _num_smc_sets = n; }
        void set_pls_validation_training_fraction(const float_type f) {
            assert((0 < f) and (f <= 1));
            _pls_training_fraction = f;
        }

        // should all be private - set from config file OR by construction
        void set_simulation(AbcSimFun * abcsf) { _simulator = abcsf; }
        void set_executable(std::string cmd) { set_simulation(new AbcExec(cmd)); }
        void set_simulator(AbcSimMPI * simulator) { set_simulation(new AbcFPtrMPI(simulator)); }
        void set_simulator(AbcSimBase * simulator) { set_simulation(new AbcFPtrBase(simulator)); }
        void set_simulator(std::string soname) { set_simulation(new AbcFPtrBase(soname.c_str())); }

        void set_retain_posterior_rank( const bool retain_rank ) { _retain_posterior_rank = retain_rank; }

        // TODO: Metric and Parameter management should be private?
        // if the moving to the "Configuration" class approach, could still have a manual configuration
        // class / methods on that building class, and not expose this on AbcSmc
        void add_next_metric(const MetricPtr m);
        // TODO model_pars => map<string, Parameter*>, check for duplicate names
        void add_next_parameter(const ParameterPtr p);

        // should be private?
        void add_modification_map(
            const ParameterPtr par,
            const ParXform * const xform
        ) { _par_modification_map[par] = xform; };

        // should be private?
        void add_par_rescale(
            const ParameterPtr par,
            const ParRescale * const par_rescale
        ) { _par_rescale_map[par] = par_rescale; }

        // should be private?
        void set_filtering_type(const FILTER &ft) { _filtering = ft; }

        // when Config class implemented, this goes there and yields a new AbcSmc
        bool parse_config(const std::string &conf_filename);

        // @param verbose the verbosity level
        //
        // @return true if the storage is built *by this call*;
        // false if already built OR if there is an error
        bool build(
            const size_t verbose
        );

        // when Storage class implemented, these disappear
//        bool build_database(const gsl_rng* RNG);
//        void set_database_filename( std::string name ) { _database_filename = name; }
        bool process_database(const gsl_rng* RNG, const size_t verbose = 0);
        bool load_SMC_sets(std::vector<std::vector<int> > &serials);
        // bool fetch_particle_parameters(
        //     sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss,
        //     vector<int> &serial, vector<Row> &par_mat, vector<unsigned long int> &seeds,
        //     const bool verbose = false
        // );
//        bool update_particle_metrics(sqdb::Db &db, vector<string> &update_metrics_strings, vector<string> &update_jobs_strings);

        // Run available particles
        //
        // @param n: the maximum number of particles to run; if negative, run all available work
        // @param keep_going: if true, keep running if a particle fails
        // @param bulk: if true, run all particles before writing to storage
        // @param verbose: level of verbosity to use
        bool simulate(
            const int n = 1,
            const bool keep_going = false, const bool bulk = false, const size_t verbose = 0
        );
        
        bool simulate_serials(
            const std::vector<int> &serial_ids,
            const bool keep_going = false, const bool bulk = false, const size_t verbose = 0
        );

        bool simulate_posterior(
            const std::vector<size_t> &posterior_ids = {}, const int smc_set = -1,
            const bool keep_going = false, const bool bulk = false, const size_t verbose = 0
        );

        size_t npar() { return _model_pars.size(); }
        size_t nmet() { return _model_mets.size(); }

        std::vector<std::shared_ptr<Mat2D>> get_particle_parameters() { return _particle_parameters; }
        std::vector<std::shared_ptr<Mat2D>> get_particle_metrics()    { return _particle_metrics; }

    private:
        friend AbcLog;
        std::unique_ptr<AbcDB> _storage;
// CORE numerical bits:

        // model parameter containers / transformations
        ParameterVec _model_pars;
        map<const ParameterPtr, const ParXform* > _par_modification_map;
        map<const ParameterPtr, const ParRescale* > _par_rescale_map;
        
        // if there are POSTERIOR parameters, must be set during construction: they source values from this
        // otherwise, empty
        std::unique_ptr<const Mat2D> _posterior;
        
        // the model itself + internal execution; defaults to unset
        AbcSimFun * _simulator = new AbcSimUnset();
        bool _do_work(
            const Mat2D &pars,
            const Coli &serials,
            const Colsz &seeds,
            const bool keep_going,
            const bool bulk,
            const size_t verbose
        );

        // uses _par_*_maps to convert from fitting space to model space
        Row _to_model_space(const Row &pars);
        Mat2D _to_model_space(const Mat2D &pars);
        // expects _model_ space parameters
        bool _run_simulator(Row &par, Row &met, const size_t rng_seed, const size_t serial);
        // model metric containers / lookups
        MetricVec _model_mets;
        const std::unique_ptr<Row> _met_vals;

        bool _retain_posterior_rank;

// SMC / ABC / PLS management:

        // how do we decide which particles to keep (Filtering) and how to propose new particles (Noise)?
        FILTER _filtering = FILTER::PLS;
        NOISE _noise = NOISE::INDEPENDENT;

        // TODO: desire a container that works exactly like a vector,
        // but with the same values from the end of the vector until some larger size
        // smc iteration containers
        size_t _num_smc_sets;
        vector<size_t> _smc_set_sizes;
        vector<size_t> _predictive_prior_sizes;
        float_type _pls_training_fraction; // how much a set should be used training PLS model (influences how much left for test => optimal component selection => ranking)

        // everything we know the about the SMC
        // vectors indexed by SMC set; row corresponds to particle, col to parameter/metric =>
        // a Col is about particle features, a Row is about parameter or metric features
        std::vector<std::vector<size_t>> _predictive_prior; // vector of row indices for particle metrics and parameters; from storage, generally 0 to n-1; arbitrary from simulation
        std::vector<std::shared_ptr<Mat2D>> _particle_metrics; // from storage, these are generally in order of posterior rank; from simulation, they are in order of serial
        std::vector<std::shared_ptr<Mat2D>> _particle_parameters; // ibid
        std::vector<std::shared_ptr<Col>> _weights; // ordered by predictive prior
        std::vector<std::shared_ptr<Row>> _doubled_variance;

        void calculate_predictive_prior_weights( const size_t set_num );

// mpi cruft
        MPI_par *_mp;
        void _particle_scheduler_mpi(const size_t t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG);
        void _particle_worker_mpi(const size_t seed, const size_t serial);

        // if resume directory set, then we will attempt to resume / store there for resumption
        // TODO this is all private, so can just cut? not currently used, except to be set by parse_config
        std::optional<std::string> _resume_directory;
        void set_resume(const bool res) { if (not res) _resume_directory.reset(); }
        bool resume() { return _resume_directory.has_value(); }
        void set_resume_directory( const std::string &res_dir ) { _resume_directory.emplace(res_dir); }

};

} // namespace ABC
#endif
