#ifndef ABC_STORAGE_H
#define ABC_STORAGE_H

#include <optional>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>

#include "sqdb.h"
#include "pls.h"

using std::string;
using std::vector;
using std::stringstream;
using std::cerr;
using std::endl;
using std::ostream;
using std::to_string;
using std::optional;
using std::ifstream;

template<typename ITERABLE>
inline string join(const ITERABLE & v, const string & delim = ", ") {
    stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << delim;
        ss << *it;
    }
    return ss.str();
}

template<typename T>
inline string joineach(
    const vector<vector<T>> & v,
    const string & outerdelim = "), (",
    const string & innerdelim = ", "
) {
    stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << outerdelim;
        ss << join(*it, innerdelim);
    }
    return ss.str();
}

// AbcStorage: abstract base class for AbcSmc storage options
// AbcStorage encapsulates different storage solutions (e.g. SQLite, ...) for AbcSmc activities
//
// In general, the steps for using AbcSmc interacts with storage as follows:
// 0. setup: create basic arrangement for analyses
//    Currently, this takes an existing AbcSmc object (having parameters, metrics, etc.). Could be doable directly from a config file?
// 1. populate: add new runs to the database
//    Writes parameters, jobs, placeholders for metrics. In general, the storage solution shouldn't care about "why" - from its perspective,
//    can ignore if this is creating a new SMC run, expanding an existing run, creating a scenario analysis plan, etc.
//    However: may make some sense to have it be somewhat aware of these distinctions, to offer convenient signatures with defaults?
//    Relative to current configuration, this might be the difficult bit to slice up - seems to happen to closely entangled with AbcSmc internals / config file?
// 2. update: add the outcome of runs to database (writes metrics + updates job status)
//    This is probably the easiest bit (at least for random-access storage solutions like SQLite)
//    Particularly, if we disentangle some of the `jobs` columns via the view approach.
// 3. read: read parameters / metrics from database
//    This might be complicated, might not - in the random access sense, it's easy (fetch these serials).
//    But, seems likely there are other standard asks (get all of a wave, get the next X open jobs)
//    And, we need to distinguish between just looking vs reading-because-we're-going-to-run-these
// 4. [introspect: determine configuration from database]
//    currently, creating an AbcSmc object has to be done from JSON config file, but for many tasks the info in the storage file should be sufficient
//    Basically, when AbcSmc isn't sampling, can just grab next-X-jobs, run them, then store the results.
struct AbcStorage {

    // given the necessary schema constraints, setup this storage solution
    // returns true if the storage solution was created. (false if fails or already exists unless overwrite is true)
    virtual bool setup(
        const vector<string> & parameters,
        const vector<string> & metrics,
        const bool transformedParameters,
        const bool overwrite,
        const bool verbose
    ) = 0;

    virtual bool populate(
        const vector<vector<double>> & parvalues,
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

    virtual Mat2D read_pars(
      const vector<size_t> serials
    ) = 0;
    
    virtual bool exists() = 0;

    virtual void storage_warning (
        const char *emsg, const string & other,
        const bool verbose = true,
        const string & prompt = "WARNING: ", 
        ostream &os = cerr
    ) const {
        if (verbose) os << prompt << emsg << endl << other << endl;
    };

    virtual void storage_error (
        const char *emsg, const string & other,
        const int code, ostream &os = cerr
    ) const {
        storage_warning(emsg, other, true, "CAUGHT exception: ", os);
        exit(code);
    };
};

struct AbcSQLite : public AbcStorage {
    AbcSQLite(const char * dbname) : _dbname(dbname) {}
    AbcSQLite(const string & dbname) : AbcSQLite(dbname.c_str()) {}

    static const string JOB_TABLE, MET_TABLE, PAR_TABLE, UPAR_TABLE;

    // `setup` establishes an AbcSmc database, creating the tables and indices; population is *not* performed.
    bool setup(
        const vector<string> & parameters,
        const vector<string> & metrics,
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
        const vector<vector<double>> & parvalues = {{5,5}, {4, 4}},
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
                    "Incommensurable `parvalues` size and cycle length: ", to_string(parvalues.size()) + " vs " + to_string(stride), verbose
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
    Mat2D read_parameters(
        const vector<size_t> & serials
    ) override {
        sqdb::Db db = connection();
        string countquery = "SELECT COUNT(*) FROM " + JOB_TABLE + (posteriorOnly ? " WHERE posterior > -1;" : ";");
        auto stmt = db.Query(countquery.c_str());
        if (stmt.Next()) {
            int count = stmt.GetField(0);
            if (count) {
                Mat2D data(count, elements.size());
                auto FROM_TABLE = _check_tables(db, { UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
                string colquery = "SELECT " + join(elements) + " FROM " + FROM_TABLE + " JOIN ";
                return data;        
            } else {
                cerr << "WARNING: no qualifying parameters found for `" << countquery << "` - returning empty Mat2D."<< endl;
                return Mat2D(0, elements.size());
            }
        } else {
            cerr << "Failed to query database" << endl;
            exit(1);
        }
    }

    // TODO, also test if it's a valid sqlite file
    bool exists() override {
        ifstream infile(_dbname);
        return infile.good();
    }

    private:
        const char* _dbname;
        vector<string> _parameters, _metrics;
        bool _transformedParameters;

        inline sqdb::Db connection() { return sqdb::Db(_dbname); }

        bool _check_tables(sqdb::Db & db, const vector<string> table_names) {
            optional<const char*> errmsg;
            optional<const int> errcode;
            try {
                bool all_exist = true;
                for(auto table_name: table_names) {
                    if (not db.TableExists(table_name.c_str())) {
                        cerr << "Table " << table_name << " does not exist in database." << endl;
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
            optional<const char*> errmsg;
            try {
                db.Query(query).Next();
            } catch (const sqdb::Exception& e) {
                errmsg.emplace(e.GetErrorMsg());
            } catch (const exception& e) {
                errmsg.emplace(e.what());
            }
            if (errmsg.has_value()) {
                storage_warning(errmsg.value(), string("Failed query:\n") + string(query), verbose);
            }
            return not errmsg.has_value();
        }

        inline bool _execute(sqdb::Db & db, const string & query) { return _execute(db, query.c_str()); }
        inline bool _execute(sqdb::Db & db, stringstream & query) { 
            auto res = _execute(db, query.str());
            query.str(string());
            query.clear();
            return res;
        }
};

const string AbcSQLite::JOB_TABLE = "job";
const string AbcSQLite::MET_TABLE = "met";
const string AbcSQLite::PAR_TABLE = "par";
const string AbcSQLite::UPAR_TABLE = "upar";

//        void set_database_filename( string name ) { _database_filename = name; }
//        void set_posterior_database_filename( string name ) { _posterior_database_filename = name; }
        // bool build_database(const gsl_rng* RNG);
        // bool process_database(const gsl_rng* RNG);
//        bool read_SMC_set_from_database (int t, ABC::Mat2D &X_orig, ABC::Mat2D &Y_orig);
        // bool read_SMC_sets_from_database(sqdb::Db &db, vector<vector<int> > &serials);

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
