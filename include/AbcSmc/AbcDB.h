#ifndef _ABCSMC_ABCDB_H_
#define _ABCSMC_ABCDB_H_

#include <string>
#include <vector>
#include <sstream>
#include <PLS/types.h>
#include <memory>

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

        // when Storage class implemented, this goes there
        bool setup();

        void read_SMC_sets(
            std::vector<std::vector<size_t>> &serials,
            std::vector<Mat2D> &parameters,
            std::vector<Mat2D> &metrics,
            const std::vector<size_t> &which_sets = {}
        );

        void read_SMC_set(
            std::vector<std::vector<size_t>> &serials,
            std::vector<std::vector<size_t>> &parameters,
            std::vector<std::vector<size_t>> &metrics,
            const size_t set_id
        ) {
            read_SMC_sets(serials, parameters, metrics, {set_id});
        };

        bool read_parameters(
            std::stringstream &select_pars_ss,
            std::stringstream &update_jobs_ss,
            std::vector<int> &serial,
            std::vector<Row> &par_mat,
            std::vector<size_t> &seeds,
            const bool verbose = false
        );

        bool write_parameters(
            const size_t verbose = 0
        );

        std::vector<std::vector<double>> read_posterior(
            const std::string &posterior_name,
            const std::vector<std::string> &_model_pars
        );

        bool write_metrics(
            const std::vector<std::string> &update_metrics_strings,
            const std::vector<std::string> &update_jobs_strings
        );

    private:
        const std::string _db_name;
        const std::unique_ptr<sqdb::Db> _db;
        bool _file_exists;
        bool _db_setup;

        bool _db_execute_stringstream(std::stringstream &ss);
        bool _db_execute_strings(std::vector<std::string> &update_buffer);

};

} // namespace ABC

#endif // _ABCSMC_ABCDB_H_