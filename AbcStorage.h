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
inline std::string join(const std::vector<T> & v, const std::string & delim = ", ") {
    std::stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << delim;
        ss << *it;
    }
    return ss.str();
}

template<typename T>
inline std::string joineach(
    const vector<vector<T>> v,
    const std::string & outerdelim = "), (",
    const std::string & innerdelim = ", "
) {
    std::stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << outerdelim;
        ss << join(*it, innerdelim);
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

    virtual bool populate(
        const std::vector<std::vector<double>> & parvalues,
        const size_t smcSet,
        const size_t cycleLength,
        const bool verbose
    ) = 0;

    // TODO: support polymorphic calls to setup; probably for AbcSmc objects (but have to beat cyclic dependency), probably for
    // AbcConfig objects
    // template <typename ABCSMC>
    // bool setup(ABCSMC & abc, const bool overwrite = false, const bool verbose = false) {
    //     return setup(
    //         abc.parameters,
    //         abc.metrics,
    //         abc.transformedParameters,
    //         overwrite,
    //         verbose
    //     );
    // }

    virtual Mat2D read_parameters(const std::vector<std::string> elements, const bool posteriorOnly) = 0;
    
    virtual bool exists() = 0;

    virtual void storage_warning (
        const char *emsg, const std::string & other,
        const bool verbose = true,
        const std::string & prompt = "WARNING: ", 
        ostream &os = std::cerr
    ) const {
        if (verbose) os << prompt << emsg << std::endl << other << std::endl;
    };

    virtual void storage_error (
        const char *emsg, const std::string & other,
        const int code, ostream &os = std::cerr
    ) const {
        storage_warning(emsg, other, true, "CAUGHT exception: ", os);
        exit(code);
    };
};

struct AbcSQLite : public AbcStorage {
    AbcSQLite(const char * dbname) : _dbname(dbname) {}
    AbcSQLite(const std::string & dbname) : AbcSQLite(dbname.c_str()) {}

    static const std::string JOB_TABLE, MET_TABLE, PAR_TABLE, UPAR_TABLE;

    // `setup` establishes an AbcSmc database, creating the tables and indices; population is *not* performed.
    bool setup(
        const std::vector<std::string> & parameters,
        const std::vector<std::string> & metrics,
        const bool transformedParameters = false,
        const bool overwrite = false, // TODO: implement overwrite
        const bool verbose = false
    ) override {
        if (exists()) {
            return false;
        } {
            _parameters = parameters; _metrics = metrics; _transformedParameters = transformedParameters;

            sqdb::Db db = connection();
            stringstream ss;

            db.BeginTransaction();
            sqdb::QueryStr qstr;

            _execute(db, qstr.Format(
              R"SQL(CREATE TABLE %s (
               serial INTEGER PRIMARY KEY ASC, smcSet INTEGER, particleIdx INTEGER,
               startTime INTEGER DEFAULT NULL, duration INTEGER DEFAULT NULL,
               status TEXT DEFAULT 'Q', posterior INTEGER,
               attempts INTEGER DEFAULT 0
              );)SQL",
              JOB_TABLE.c_str()
            ));

            _execute(db, qstr.Format(
              "CREATE INDEX idx1 ON %s (status, attempts);",
              JOB_TABLE.c_str()
            ));

            _execute(db, qstr.Format(
              R"SQL(CREATE TABLE %s (
                serial INTEGER PRIMARY KEY ASC, seed blob,
                %s real
              );)SQL",
              PAR_TABLE.c_str(),
              join(parameters, " real, ").c_str()
            ));

            if (transformedParameters) {
            _execute(db, qstr.Format(
              R"SQL(CREATE TABLE %s (
                serial INTEGER PRIMARY KEY ASC, seed blob,
                %s real
              );)SQL",
              UPAR_TABLE.c_str(),
              join(parameters, " real, ").c_str()
            ));
            }

            _execute(db, qstr.Format(
              R"SQL(CREATE TABLE %s (
                serial INTEGER PRIMARY KEY ASC,
                %s real
              );)SQL",
              MET_TABLE.c_str(),
              join(metrics, " real, ").c_str()
            ));

            db.CommitTransaction();

            return true;
        }
    }

    // `populate` populates (via INSERT) the database given parameters, and sets the status of each job to 'Q' (queued).
    bool populate(
        const std::vector<std::vector<double>> & parvalues = {{5,5}, {4, 4}},
        const size_t smcSet = 0,
        const vector<size_t> posterior_rank = {},
        // 0 => particleIdx = 0, 1, 2, ..., parvalues.size()
        // otherwise, assumes parvalues.size() == cycleLength*n,
        // and particleIdx = 0, 1, 2, ..., cycleLength-1 (repeat n times)
        const bool verbose = false
    ) override {
        if (not exists()) {
            if (verbose) {
                storage_warning("Database is not setup: ", _dbname);
            }
            return false;
        } else if (parvalues.size() != 0) {
            const size_t stride = (posterior_rank.size() == 0) ? parvalues.size() : posterior_rank.size();
            if ((parvalues.size() % stride) != 0) {
                storage_warning(
                    "Incommensurable `parvalues` size and cycle length: ", std::to_string(parvalues.size()) + " vs " + std::to_string(stride), verbose
                );
                return false;
            }

            sqdb::Db db = connection();
            sqdb::QueryStr qstr;
            db.BeginTransaction();

            // ASSERT: not required that the tables are empty.
            // if they aren't empty, however: their autoincrement is at the same place

            // insert parameter values in PAR_TABLE
            _execute(db, qstr.Format(
                SQDB_MAKE_TEXT("INSERT INTO %s (%s) VALUES (%s);"),
                PAR_TABLE.c_str(), join(_parameters).c_str(), joineach(parvalues).c_str()
            ));

            // insert parameter values in JOB_TABLE
            _execute(db, qstr.Format(
                SQDB_MAKE_TEXT("INSERT INTO %s (%s) VALUES (%s);"),
                JOB_TABLE.c_str(), join({ "smcSet", "particleIdx"}).c_str(), joineach(parvalues).c_str()
            ));

            // // capture created serials
            // auto findSerial = db.Query(qstr.Format(
            //     SQDB_MAKE_TEXT("SELECT last_insert_rowid() FROM %s;"),
            //     PAR_TABLE.c_str()
            // ));
            // findSerial.Next();
            // const int maxserial = findSerial.GetField(0);
            // const int minserial = maxserial - parvalues.size() + 1;


            // if there are transformed parameters, transform and insert into UPAR_TABLE
            // insert unfilled rows into MET_TABLE
            // insert Q'd jobs into JOB_TABLE

            db.CommitTransaction();
            return true;
        } else {
            storage_warning("Empty `parvalues` provided to populate: ", _dbname, verbose);
            return false;
        }
    }

// throws SQDB exceptions
    Mat2D read_parameters(const vector<std::string> elements, const bool posteriorOnly) override {
        sqdb::Db db = connection();
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
        std::vector<std::string> _parameters, _metrics;
        bool _transformedParameters;

        inline sqdb::Db connection() { return sqdb::Db(_dbname); }

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
        // these do error handling (and resource management stringstream)
        
        inline bool _execute(sqdb::Db & db, const char * query, const bool verbose = true) {
            std::optional<const char*> errmsg;
            try {
                db.Query(query).Next();
            } catch (const sqdb::Exception& e) {
                errmsg.emplace(e.GetErrorMsg());
            } catch (const exception& e) {
                errmsg.emplace(e.what());
            }
            if (errmsg.has_value()) {
                storage_warning(errmsg.value(), std::string("Failed query:\n") + std::string(query), verbose);
            }
            return not errmsg.has_value();
        }

        inline bool _execute(sqdb::Db & db, const std::string & query) { return _execute(db, query.c_str()); }
        inline bool _execute(sqdb::Db & db, stringstream & query) { 
            auto res = _execute(db, query.str());
            query.str(std::string());
            query.clear();
            return res;
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
