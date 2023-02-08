#ifndef ABC_STORAGE_H
#define ABC_STORAGE_H

#include <optional>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "sqdb.h"
#include "AbcUtil.h"

template<typename T>
inline std::string join(vector<T> v, const std::string & delim = ", ") {
    std::stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << delim;
        ss << *it;
    }
    return ss.str();
}

using namespace ABC;

struct AbcStorage {

    virtual bool setup(
        const std::vector<std::string> & parameters,
        const std::vector<std::string> & metrics,
        const bool transformedParameters,
        const bool overwrite,
        const bool verbose
    ) = 0;

    virtual Mat2D read_parameters(const std::vector<std::string> elements, const bool posteriorOnly) = 0;
    virtual bool write_data() = 0;
    virtual bool exists() = 0;

    const void storage_warning(
        const char *emsg, const std::string & other,
        const std::string & prompt = "WARNING: ", 
        ostream &os = std::cerr
    ) {
        os << prompt << emsg << std::endl << other << std::endl;
    };

    const void storage_error(
        const char *emsg, const std::string & other,
        const int code, ostream &os = std::cerr
    ) {
        storage_warning(emsg, other, "CAUGHT exception: ", os);
        exit(code);
    };
};

struct AbcSQLite : public AbcStorage {
    AbcSQLite(const char * dbname) : _dbname(dbname) {}
    AbcSQLite(const std::string & dbname) : AbcSQLite(dbname.c_str()) {}

    static const std::string JOB_TABLE, MET_TABLE, PAR_TABLE, UPAR_TABLE;

    bool setup(
        const std::vector<std::string> & parameters,
        const std::vector<std::string> & metrics,
        const bool transformedParameters = false,
        const bool overwrite = false,
        const bool verbose = false
    ) override {
        if (exists()) {
            return false;
        } {
            sqdb::Db db = *connection();
            stringstream ss;

            db.BeginTransaction();

            ss << "CREATE TABLE " << JOB_TABLE << " ( serial int primary key asc, smcSet int, particleIdx int, startTime int, duration real, status text, posterior int, attempts int );";
            _execute_stringstream(db, ss);

            ss << "CREATE idx1 ON " << JOB_TABLE << " (status, attempts);";
            _execute_stringstream(db, ss);

            ss << "CREATE TABLE " << PAR_TABLE << " ( serial int primary key, seed blob, " << join(parameters, " real,") << " real);";
            _execute_stringstream(db, ss);

            if (transformedParameters) {
                ss << "CREATE TABLE " << UPAR_TABLE << " ( serial int primary key, seed blob, " << join(parameters, " real,") << ");";
                _execute_stringstream(db, ss);
            }

            ss << "CREATE TABLE " << MET_TABLE << " ( serial int primary key, " << join(metrics, " real,") << ");";
            _execute_stringstream(db, ss);

            db.CommitTransaction();
            return true;
        }
    }

    bool write_data() override {
        return false;
    }

// throws SQDB exceptions
    Mat2D read_parameters(const vector<std::string> elements, const bool posteriorOnly) override {
        sqdb::Db db = *connection();
        std::string countquery = "SELECT COUNT(*) FROM " + JOB_TABLE + (posteriorOnly ? " WHERE posterior > -1;" : ";");
        auto stmt = db.Query(countquery.c_str());
        if (stmt.Next()) {
            int count = stmt.GetField(0);
            if (count) {
                Mat2D data(count, elements.size());
                auto FROM_TABLE = _check_tables(db, { UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
                std::string colquery = "SELECT " + join(elements) + " FROM " + FROM_TABLE + " JOIN ";
                return data;        
            } else {
                std::cerr << "WARNING: no qualifying parameters found for `" << countquery << "` - returning empty Mat2D."<< std::endl;
                return Mat2D(0, elements.size());
            }
        } else {
            std::cerr << "Failed to query database" << std::endl;
            exit(1);
        }
    }

    // TODO, also test if it's a valid sqlite file
    bool exists() override {
        std::ifstream infile(_dbname);
        return infile.good();
    }

    private:
        const char* _dbname;

        // TODO: should this return a pointer?
        sqdb::Db* connection() { return new sqdb::Db(_dbname); }

        bool _check_tables(sqdb::Db & db, const std::vector<std::string> table_names) {
            std::optional<const char*> errmsg;
            std::optional<const int> errcode;
            try {
                bool all_exist = true;
                for(auto table_name: table_names) {
                    if (not db.TableExists(table_name.c_str())) {
                        std::cerr << "Table " << table_name << " does not exist in database." << std::endl;
                        all_exist = false;
                    }
                }
                return all_exist;
            } catch (const sqdb::Exception& e) {
                errmsg.emplace(e.GetErrorMsg()); errcode.emplace(-212);
            } catch (const exception& e) {
                errmsg.emplace(e.what()); errcode.emplace(-213);
            }

            // only reached if exception thrown
            stringstream sserr;
            sserr << "Failed while checking whether the following tables exist:";
            for(auto table_name: table_names) sserr << " " << table_name;
            storage_error(errmsg.value(), sserr.str(), errcode.value());

            return false;
        };

        // no transaction handling here - assumed done around this in calling scope
        // this does resource management for the stringstream / error handling
        bool _execute_stringstream(sqdb::Db & db, stringstream &ss) {   
            bool db_success = false;
            std::optional<const char*> errmsg;
            try {
                // TODO - this doesn't actually mean db_success. "Next" == true means there are rows.
                db_success = db.Query(ss.str().c_str()).Next();
            } catch (const sqdb::Exception& e) {
                errmsg.emplace(e.GetErrorMsg());
            } catch (const exception& e) {
                errmsg.emplace(e.what());
            }
            if (not db_success) {
                if (!errmsg.has_value()) {
                    errmsg.emplace("Unknown, non-throwing error.");
                }
                storage_warning(errmsg.value(), std::string("Failed query:\n") + ss.str());
            }
            // clear contents & state of the stringstream
            ss.str(std::string());
            ss.clear();
            return db_success;
        }

};

const std::string AbcSQLite::JOB_TABLE = "job";
const std::string AbcSQLite::MET_TABLE = "met";
const std::string AbcSQLite::PAR_TABLE = "par";
const std::string AbcSQLite::UPAR_TABLE = "upar";

//        void set_database_filename( std::string name ) { _database_filename = name; }
//        void set_posterior_database_filename( std::string name ) { _posterior_database_filename = name; }
        // bool build_database(const gsl_rng* RNG);
        // bool process_database(const gsl_rng* RNG);
//        bool read_SMC_set_from_database (int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);
        // bool read_SMC_sets_from_database(sqdb::Db &db, std::vector<std::vector<int> > &serials);

        // bool sql_particle_already_done(sqdb::Db &db, const string sql_job_tag, string &status);
        // bool fetch_particle_parameters(
        //     sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss,
        //     vector<int> &serial, vector<ABC::Row> &par_mat, vector<unsigned long int> &seeds,
        //     const bool verbose = false
        // );
        // bool update_particle_metrics(sqdb::Db &db, vector<string> &update_metrics_strings, vector<string> &update_jobs_strings);

// inline bool update_stage(
//     sqdb::Db & db, string stage
// ) {
//     try {
//         db.Query(stage.c_str()).Next();
//     } catch (const Exception& e) {
//         cerr << "CAUGHT E: ";
//         cerr << e.GetErrorMsg() << endl;
//         cerr << "Failed on:" << endl;
//         cerr << stage << endl;
//         return false;
//     } catch (const exception& e) {
//         cerr << "CAUGHT e: ";
//         cerr << e.what() << endl;
//         cerr << "Failed on:" << endl;
//         cerr << stage << endl;
//         return false;
//     }
//     return true;
// }

#endif
