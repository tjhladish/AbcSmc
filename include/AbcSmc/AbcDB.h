#ifndef _ABCSMC_ABCDB_H_
#define _ABCSMC_ABCDB_H_

#include <string>
#include <vector>
#include <sstream>
#include <memory>

#include <PLS/types.h>
#include <AbcSmc/Parameter.h>
#include <AbcSmc/Metric.h>

// forward declaring sqdb::Db
namespace sqdb {
    class Db;
}

// Storage classes need to implement:
//  - reading entire SMC sets: parameters, metrics, weights, seeds
//  - reading a single SMC set
//  - writing parameters, metrics, weights, seeds
//  - reading parameters, metrics, weights, seeds
//  - doing self setup
//  - indicating state for read / write / etc == completeness of a wave, initial setup complete?

namespace ABC {

class AbcDB {

    public:
        AbcDB(const std::string & path);

        // @return true if the storage container is structurally setup
        bool is_setup();

        // @param set: if -1, check last set; otherwise check only the specified set
        // @return true if the storage container is ready to provide / receive data
        bool is_set_done(const int set = -1);

        // @param pars: vector of Parameter (pointers); forms the basis for indexing, etc
        // @param mets: vector of Metric (pointers); forms the basis for indexing, etc
        // @param verbose: level of verbosity to use during setup
        // @return true if storage is ready to receive / provide data
        // n.b. does NOT schedule any jobs
        // n.b. repeated invocations should not cause overwrites
        // TODO: repeated invocations should warn
        // TODO: repeated invocations should warn especially if offered different pars / mets
        bool setup(
            const ParameterVec &pars,
            const MetricVec &mets,
            const bool has_transforms,
            const size_t verbose = 0
        );

        // @param pars: vector of parameter *values*; must be in the same order as the parameters used in setup
        // @param seeds: the associated seeds; assert length(pars) == length(seeds)
        // @param verbose: level of verbosity to use
        bool create_SMC_set(
            const std::shared_ptr<Mat2D> &pars,
            const std::vector<size_t> &seeds,
            const size_t verbose = 0
        );

        // Bulk-read an SMC set from the storage container.
        //
        // Conditions: `parameters`, `metrics`, and `weights` should be empty. If successful read, they will be filled,
        // in order matching the contents of `which_sets`
        //
        // @param parameters: vector of parameter *values*; will be in the same column order as the parameters used in setup.
        // @param metrics: vector of metric *values*; will be in the same column order as the metrics used in setup.
        // @param weights: vector of weights; will be filled with the weights for the specified set(s).
        // @param which_sets: which sets to read; if empty, read all sets.
        // "negative" indices are interpreted as "from the end" (i.e. -1 is the last set)
        // TODO: basic implementation only reads the last set with a negative index
        // @return true if the read was successful
        bool read_SMC_sets(
            std::vector<std::shared_ptr<Mat2D>> &parameters,
            std::vector<std::shared_ptr<Mat2D>> &metrics,
            std::vector<std::shared_ptr<Col>> &weights,
            const std::vector<int> &which_sets = { -1 }
        );

        // @param mets: vector of metric *values*; must be in the same order as the metric used in setup
        // @param serials: the identifying serials; assert length(mets) == length(serials)
        // @param verbose: level of verbosity to use
        bool update_SMC_set(
            const std::shared_ptr<Mat2D> &mets,
            const std::vector<size_t> &serials,
            const size_t verbose = 0
        );

        bool get_posterior_serials(
            const std::vector<size_t> &posterior, // if empty, get all > -1
            Coli &serials,
            const int smc_set = -1 // if -1, get last set
        );

        // Get parameters that need simulation results
        //
        // @param parameters: a Mat2D that will be filled with parameter values
        // @param seeds: a Col that will be filled with seeds
        // @param serials: if empty, a Col to fill; OTHERWISE, the serials to get
        // @param n: the number of parameters to get (assert: `serials` is empty)
        bool get_work(
            Mat2D &parameters,
            Colsz &seeds,
            Coli  &serials,
            const int n = 1
        );

        bool put_work(
            const Row &metrics,
            const int serial,
            const size_t start_time,
            const size_t verbose = 0
        );

    private:
        const std::string _db_name;
        const std::unique_ptr<sqdb::Db> _db;
        bool _file_exists;
        bool _db_setup;

        bool _db_execute_strings(std::vector<std::string> &update_buffer);

};

} // namespace ABC

#endif // _ABCSMC_ABCDB_H_