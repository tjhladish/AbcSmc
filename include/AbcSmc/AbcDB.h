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

namespace ABC {

class AbcDB {

    public:
        AbcDB(const std::string & path);
        //  : _db_name(path) {
        //     _file_exists = std::filesystem::exists(_db_name);
        //     _db(_db_name.c_str());
        // };

        bool is_setup() const { return _db_setup; }
        bool has_work() const;

        // @param verbose: level of verbosity to use during setup
        // @return true if storage is ready to receive / provide data
        // n.b. repeated invocations should not cause overwrites
        bool setup(
            const ParameterVec &pars,
            const MetricVec &mets,
            const bool has_transforms,
            const size_t verbose = 0
        );

        bool read_SMC_sets(
            std::vector<std::vector<int>> &serials,
            std::vector<std::shared_ptr<Mat2D>> &parameters,
            std::vector<std::shared_ptr<Mat2D>> &metrics,
            const std::vector<size_t> &which_sets = {}
        );

        bool read_SMC_set(
            std::vector<std::vector<int>> &serials,
            std::vector<std::shared_ptr<Mat2D>> &parameters,
            std::vector<std::shared_ptr<Mat2D>> &metrics,
            const size_t set_id
        );

        bool read_parameters(
            const size_t set_num,
            std::vector<Row> &par_mat,
            const vector<size_t> &posterior_ranks,
            const bool verbose = false
        );

        bool read_parameters(
            std::vector<int> &serial,
            std::vector<Row> &par_mat,
            const bool verbose = false
        );

        bool write_parameters(
            const Mat2D &pars,
            const Mat2D &upars,
            const std::vector<size_t> &seeds,
            const size_t set_num,
            const vector<size_t> &posterior_ranks,
            const size_t verbose
        );

        Mat2D read_posterior(
            const std::string &posterior_name,
            const std::vector<std::string> &_model_pars
        );

        bool read_metrics(
            const std::vector<std::string> &update_metrics_strings,
            const std::vector<std::string> &update_jobs_strings
        );

        bool write_metrics(
            const Mat2D &mets
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