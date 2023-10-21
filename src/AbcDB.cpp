
#include <filesystem> // exists
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <algorithm> // all_of
#include <assert.h>

#include <Eigen/Dense>
#include "sqdb.h"

#include <AbcSmc/AbcDB.h>

// need a positive int that is very unlikely
// to be less than the number of particles
#define STOP_TAG 10000000

using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;
using std::setw;

using namespace sqdb;

const int PREC = 5;
const string JOB_TABLE  = "job";
const string MET_TABLE  = "met";
const string PAR_TABLE  = "par";
const string UPAR_TABLE = "upar";

#define SQDB_TRYBLOCK(TRY, MESSAGE) try { \
    TRY \
} catch (const Exception &e) { \
    std::cerr << "CAUGHT E: " << e.GetErrorMsg() << std::endl; \
    std::cerr << MESSAGE << std::endl; \
    exit(-212); \
} catch (const std::exception &e) { \
    std::cerr << "CAUGHT e: " << e.what() << std::endl; \
    std::cerr << MESSAGE << std::endl; \
    exit(-213); \
}

bool _db_execute_stringstream(
    Db &db, stringstream &ss
) {
    bool db_success = false;
    SQDB_TRYBLOCK(
        {
            db_success = db.Query(ss.str().c_str()).Next();
            ss.str(string());
            ss.clear();
        },
        "Failed query:" << std::endl << ss.str()
    )
    return db_success;
}

// @param db, the target database
// @param table_names, a vector of table names to check for existence
// @param verbose, the verbosity level for output
//
// @return true if all tables exist, false otherwise
// @throws if any sql/sqlite errors occur
bool _db_tables_exist(
    Db &db,
    const std::vector<std::string> &table_names,
    const size_t verbose = 0
) {
    bool tables_exist = true;
    QueryStr qstr;
    SQDB_TRYBLOCK(
        for (auto table_name : table_names) {
            Statement s = db.Query(qstr.Format(SQDB_MAKE_TEXT(
                "SELECT COUNT(*) FROM sqlite_master WHERE type == 'table' AND name = '%s';"
            ), table_name.c_str()));           
            s.Next();
            const int count = s.GetField(0);
            if (count < 1) {
                if (verbose > 0) {
                    std::cerr << "Table " << table_name << " does not exist in database." << std::endl;
                }
                tables_exist = false;
            }
        },
        "Failed while checking whether the following tables exist:"; for (auto table_name: table_names) std::cerr << " " << table_name; std::cerr
    )
    return tables_exist;
}

string _build_sql_create_par_string(const ABC::ParameterVec &pars) {
    stringstream ss;
    for (auto par : pars) { ss << par->get_short_name() << " real, "; }
    // TODO delete the ss last two chars
    return ss.str();
}

string _build_sql_create_met_string(const ABC::MetricVec &mets) {
    stringstream ss;
    for (auto met : mets) { ss << met->get_short_name() << " real, "; }
    // TODO delete the ss last two chars
    return ss.str();
}

namespace ABC {

AbcDB::AbcDB(const std::string & path) : 
    _db_name(path), _db(new Db(_db_name.c_str())) {
    _file_exists = std::filesystem::exists(_db_name);
    _db_setup = _file_exists && _db_tables_exist(*_db, {JOB_TABLE, PAR_TABLE, MET_TABLE});
};
    
bool AbcDB::setup(
    const ParameterVec &pars,
    const MetricVec &mets,
    const bool has_transforms,
    const size_t verbose
) {
    if (not _db_tables_exist(db, {JOB_TABLE, PAR_TABLE, MET_TABLE})) {
        stringstream ss;
        _db->Query("BEGIN EXCLUSIVE;").Next();

        ss << "create table " << JOB_TABLE << " ( serial int primary key asc, smcSet int, particleIdx int, startTime int, duration real, status text, posterior int, attempts int );";
        _db_execute_stringstream(*_db, ss);

        ss << "create index idx1 on " << JOB_TABLE << " (status, attempts);";
        _db_execute_stringstream(*_db, ss);

        ss << "create table " << PAR_TABLE <<
            " ( serial int primary key, seed blob, " <<
            _build_sql_create_par_string(pars) << ");";

        _db_execute_stringstream(*_db, ss);

        if (has_transforms) {
            ss << "create table " << UPAR_TABLE <<
                " ( serial int primary key, seed blob, " <<
                _build_sql_create_par_string(pars) << ");";
            
            _db_execute_stringstream(*_db, ss);
        }

        ss << "create table " << MET_TABLE << " ( serial int primary key, " <<
            _build_sql_create_met_string(mets) << ");";
        _db_execute_stringstream(*_db, ss);

        _db->CommitTransaction();
        _db_setup = true;
    } else {
        return _db_setup;
    }
}

bool AbcDB::write_parameters(
    const Mat2D &pars,
    const Mat2D &upars,
    const std::vector<size_t> &seeds,
    const size_t set_num,
    const size_t verbose
) {

    assert(_db_setup);

    if (verbose > 0) {
        std::cerr << std::setprecision(PREC);
    }

     // TODO - get last serial from DB
    const size_t last_serial = 0;

    QueryStr qstr;

    _db->BeginTransaction();

    for (size_t i = 0; i < pars.rows(); i++) {
        const int serial = last_serial + 1 + i;
        auto par_row = pars.row(i);

        ss << "insert into " << JOB_TABLE << " values ( " << serial << ", "
                                            << set_num << ", "
                                            << i << ", "
                                            << time(NULL)
                                            << ", NULL, 'Q', -1, 0 );";

        _db_execute_stringstream(*_db, ss);

        ss << "insert into " << PAR_TABLE << " values ( " << 
            serial << ", '" << seeds[i] << "'";
            for (auto par : par_row) ss << ", " << par;
            ss << " );";

        _db_execute_stringstream(*_db, ss);

        if (upars.size() > 0) {
            auto upar_row = upars.row(i);
            ss << "insert into " << UPAR_TABLE << " values ( " << 
                serial << ", '" << seeds[i] << "'";
                for (auto upar : upar_row) ss << ", " << upar;
                ss << " );";
            _db_execute_stringstream(*_db, ss);
        }

        ss << "insert into " << MET_TABLE << " values ( " << serial << " );";
        //cerr << "attempting: " << ss.str() << endl;
        _db_execute_stringstream(*_db, ss);
    }

    _db->CommitTransaction();

    return true;
}



}

Mat2D read_posterior(
    const std::string &_posterior_database_filename,
    const std::vector<std::string> &_model_pars
) {
    // TODO - handle sql/sqlite errors
    // TODO - if "posterior" database doesn't actually have posterior values, this will fail silently

    // determine dimensions of results we're going to read in
    sqdb::Db post_db(_posterior_database_filename.c_str());
    Statement posterior_query = post_db.Query(("select count(*) from " + JOB_TABLE + " where posterior > -1;").c_str());
    posterior_query.Next();
    const int posterior_size = posterior_query.GetField(0); // num of rows

    int num_posterior_pars = 0; // num of cols
    string posterior_strings = "";
    for (auto col_name : _model_pars) {
        ++num_posterior_pars;
        if (posterior_strings != "") {
            posterior_strings += ", ";
        }
        posterior_strings += col_name;
    }

    auto posterior = std::vector<std::vector<double>>(
        posterior_size, // number of rows
        std::vector<double>(num_posterior_pars)
    );

    // Identify table to pull parameter values from
    string post_par_table = _db_tables_exist(post_db, {UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
    stringstream ss;
    ss << "select " << posterior_strings << " from " << post_par_table << " P, " << JOB_TABLE << " J where P.serial = J.serial and posterior > -1;";
    posterior_query = post_db.Query(ss.str().c_str());

    int r = 0;
    while(posterior_query.Next()) {
        for (int c = 0; c < num_posterior_pars; ++c) {
            posterior[r][c] = static_cast<double>(posterior_query.GetField(c));
        }
        ++r;
    }

    return posterior;
}

// 
// Read in existing sets and do particle filtering as appropriate
bool AbcDB::read_SMC_sets(
    std::vector<std::vector<size_t>> &serials,
    std::vector<std::vector<size_t>> &parameters,
    std::vector<std::vector<size_t>> &metrics
) {

    assert(_db_setup);

    // make sure set t is a completed set
    QueryStr qstr;
    Statement s = db.Query(("select smcSet, count(*), COUNT(case status when 'D' then 1 else null end) from " + JOB_TABLE + " group by smcSet order by smcSet;").c_str());
    
    serials.clear();
    parameters.clear();
    metrics.clear();

    bool set_size_test_success = true;
    while (s.Next()) {
        int t = s.GetField(0);
        int set_size = s.GetField(1);
        int completed_set_size = s.GetField(2);
        if (set_size != completed_set_size) {
            cerr << "ERROR: Failed to read SMC set from database because not all particles are complete in set " << t << "\n";
            return false;
        }
        const int json_set_size = get_smc_size_at(t);
        if (set_size != json_set_size) {
            cerr << "ERROR:\tSet size for one or more sets does not agree between configuration file and database:" << endl
                 << "\tSet " << t << " in configuration file has size " << json_set_size << " vs size " << set_size << " in database." << endl;
//                 << "\tNB: You may want to edit the configuration file to have an array of set sizes that reflect what is already in the database." << endl;
            set_size_test_success = false;
            break;
            //return false;
        }

        _particle_parameters.push_back( Mat2D::Zero( completed_set_size, npar() ) );
        _particle_metrics.push_back( Mat2D::Zero( completed_set_size, nmet() ) );

        // join all three tables for rows with smcSet = t, slurp and store values
        string select_str = "select J.serial, J.particleIdx, J.posterior, " + _build_sql_select_par_string("") + ", " + _build_sql_select_met_string()
                            + "from " + JOB_TABLE + " J, " + MET_TABLE + " M, " + PAR_TABLE + " P where J.serial = M.serial and J.serial = P.serial "
                            + "and J.smcSet = " + to_string((long long) t) + ";";

        serials.push_back( vector<size_t>(completed_set_size) );

        Statement s2 = db.Query( select_str.c_str() );

        int particle_counter = 0;
        vector<std::pair<int, int> > posterior_pairs;
        while (s2.Next()) {
            size_t offset = 3; // first values are J.serial, J.particleIdx, J.rank
            const int serial = s2.GetField(0);
            const int particle_idx = s2.GetField(1);
            const int posterior_rank = s2.GetField(2);

            if (particle_counter != particle_idx) cerr << "ERROR: particle_counter != particle_idx (" << particle_counter << " != " << particle_idx << ")\n";
            assert(particle_counter == particle_idx);
            serials[t][particle_counter] = serial;
            if (posterior_rank > -1) posterior_pairs.push_back(make_pair( posterior_rank, particle_idx ) );
            for(size_t i = offset; i < offset + npar(); i++) _particle_parameters[t](particle_counter,i-offset) = (double) s2.GetField(i);
            offset += npar();
            for(size_t i = offset; i < offset + nmet(); i++) _particle_metrics[t](particle_counter,i-offset) = (double) s2.GetField(i);
            particle_counter++;
        }

        //const int posterior_size = posterior_pairs.size() > 0 ? posterior_pairs.size() : _next_predictive_prior_size;
        //_predictive_prior.push_back( vector<int>(posterior_size) );
        if (posterior_pairs.size() > 0) { // This is a set that has already undergone particle filtering & ranking
            _predictive_prior.push_back( vector<size_t>(posterior_pairs.size()) );
            for (size_t i = 0; i < posterior_pairs.size(); i++) {
                const int rank = posterior_pairs[i].first;
                const int idx = posterior_pairs[i].second;
                _predictive_prior.back()[rank] = idx;
            }
        } else { // Otherwise, do the filtering now and update the DB
            // filtering: rank all particles according to fitting scheme
            switch(_filtering) {
                case FILTER::PLS: { _predictive_prior.push_back(
                    particle_ranking_PLS(_particle_metrics[t], _particle_parameters[t], _met_vals, _pls_training_fraction)
                ); break; }
                case FILTER::SIMPLE: { _predictive_prior.push_back(
                    particle_ranking_simple(_particle_metrics[t], _particle_parameters[t], _met_vals )
                ); break; }
                default: std::cerr << "ERROR: Unsupported filtering method: " << _filtering << std::endl; return false;
            }

            // trim that ranking to the size of the next predictive prior (n.b. this only drops from the indexing vector, not the actual data)
            const size_t next_pred_prior_size = get_pred_prior_size_at(t);
            _predictive_prior.back().resize(next_pred_prior_size);

            auto posterior_pars = _particle_parameters[t](_predictive_prior[t], Eigen::placeholders::all);
            auto posterior_mets = _particle_metrics[t](_predictive_prior[t], Eigen::placeholders::all);

            AbcLog::filtering_report(this, t, posterior_pars, posterior_mets);

            vector<string> update_strings(next_pred_prior_size);
            for (size_t i = 0; i < next_pred_prior_size; i++) { // best to worst performing particle in posterior?
                const int particle_idx = _predictive_prior[t][i]; // TODO how might slice work here?
                const int particle_serial = serials[t][particle_idx];
                stringstream ss;
                ss << "update " << JOB_TABLE << " set posterior = " << i << " where serial = " << particle_serial << ";";
                update_strings[i] = ss.str();
            }
            _db_execute_strings(db, update_strings);
        }

        calculate_predictive_prior_weights( t );
    }

    if (not set_size_test_success) {
        s = db.Query(("select count(*) from " + JOB_TABLE + " group by smcSet order by smcSet;").c_str());
        vector<int> db_set_sizes;
        while (s.Next()) { db_set_sizes.push_back(s.GetField(0)); }
        cerr << "\tNB: You may want to edit the configuration file to have an array of set sizes that reflect what is already in the database." << endl
             << "\t    Sizes of sets currently in database: [";
        cerr_vector(db_set_sizes, ", ");
        cerr << "]" << endl;
        exit(1);
    }

    return true;
}

string AbcSmc::_build_sql_select_par_string( string tag = "" ) {
    stringstream ss;
    for (size_t i = 0; i<npar()-1; i++) { ss << "P." << _model_pars[i]->get_short_name() << tag << ", "; }
    ss << "P." << _model_pars.back()->get_short_name() << tag << " ";
    return ss.str();
}


string AbcSmc::_build_sql_select_met_string() {
    stringstream ss;
    for (size_t i = 0; i<nmet()-1; i++) { ss << "M." << _model_mets[i]->get_short_name() << ", "; }
    ss << "M." << _model_mets.back()->get_short_name() << " ";
    return ss.str();
}

bool AbcSmc::build_database(
    const gsl_rng* RNG
) {

    // create the DB handle
    sqdb::Db db(_database_filename.c_str());

    // create the tables if they don't exist; if they already exist, this function is a no-op
    stringstream ss;
    if ( !_db_tables_exist(db, {JOB_TABLE}) and !_db_tables_exist(db, {PAR_TABLE}) and !_db_tables_exist(db, {MET_TABLE}) ) {
        db.Query("BEGIN EXCLUSIVE;").Next();
        ss << "create table " << JOB_TABLE << " ( serial int primary key asc, smcSet int, particleIdx int, startTime int, duration real, status text, posterior int, attempts int );";
        _db_execute_stringstream(db, ss);

        ss << "create index idx1 on " << JOB_TABLE << " (status, attempts);";
        _db_execute_stringstream(db, ss);

        ss << "create table " << PAR_TABLE << " ( serial int primary key, seed blob, " << _build_sql_create_par_string("") << ");";
        _db_execute_stringstream(db, ss);

        if (_par_modification_map.size()) {
            ss << "create table " << UPAR_TABLE << " ( serial int primary key, seed blob, " << _build_sql_create_par_string("") << ");";
            _db_execute_stringstream(db, ss);
        }

        ss << "create table " << MET_TABLE << " ( serial int primary key, " << _build_sql_create_met_string("") << ");";
        _db_execute_stringstream(db, ss);
        db.CommitTransaction();
    } else {
        return false;
    }

    const size_t set_num = 0;
    const size_t num_particles = get_smc_size_at(set_num);
    std::vector<size_t> posterior_ranks = {};
    Mat2D pars = sample_priors(RNG, num_particles, _posterior, _model_pars, posterior_ranks);

    db.Query("BEGIN EXCLUSIVE;").Next();

    for (size_t i = 0; i < num_particles; i++) {
        auto parrow = pars.row(i);      
        auto posterior_rank = _retain_posterior_rank ? posterior_ranks[i] : -1;  

        QueryStr qstr;

        db.Query(qstr.Format(SQDB_MAKE_TEXT("insert into %s values ( %d, %d, %d, %d, NULL, 'Q', %d, 0 );"), JOB_TABLE.c_str(), i, set_num, i, time(NULL), posterior_rank)).Next();

        Statement s = db.Query(("select last_insert_rowid() from " + JOB_TABLE + ";").c_str());
        s.Next();
        const int rowid = ((int) s.GetField(0)) - 1; // indexing should start at 0

        const unsigned long int seed = gsl_rng_get(RNG); // seed for particle
        ss << "insert into " << PAR_TABLE << " values ( " << rowid << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << parrow[j]; ss << " );";
        _db_execute_stringstream(db, ss);

        if (_par_modification_map.size()) {
            const Row upars = _to_model_space(parrow);
            ss << "insert into " << UPAR_TABLE << " values ( " << rowid << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << upars[j]; ss << " );";
            _db_execute_stringstream(db, ss);
        }

        ss << "insert into " << MET_TABLE << " values ( " << rowid; for (size_t j = 0; j < nmet(); j++) ss << ", NULL"; ss << " );";
        _db_execute_stringstream(db, ss);
    }
    db.CommitTransaction();
    return true;
}


bool AbcSmc::fetch_particle_parameters(
    sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss,
    vector<int> &serials, vector<Row> &par_mat, vector<unsigned long int> &rng_seeds,
    const bool verbose
) {
    bool db_success = false;
    try {
        if (verbose) {
            std::cerr << "Attempting: " << select_pars_ss.str() << std::endl;
        }
        db.Query("BEGIN EXCLUSIVE;").Next();
        cerr << "Lock obtained" << endl;
        Statement s = db.Query(select_pars_ss.str().c_str());
        vector<string> job_strs;

        while (s.Next()) { // grap par vals and find out which jobs need to be updated (changed to running (R), incremented attempts, updated timestamp)
            Row pars(npar());
            const int serial = (int) s.GetField(0);
            serials.push_back( serial );
            const unsigned long int seed = (unsigned long int) s.GetField(1);
            rng_seeds.push_back(seed);
            const int field_offset = 2;
            for (size_t i = 0; i < npar(); i++) pars[i] = (double) s.GetField(i+field_offset);
            par_mat.push_back(pars);

            //job_ss << serial << ";";
            string job_str = update_jobs_ss.str() + to_string((long long) serial) + ";";
            job_strs.push_back(job_str);
        }

        for (string job_str: job_strs) {
            if (verbose) {
                std::cerr << "Attempting: " << job_str << std::endl;
            }
            db.Query(job_str.c_str()).Next(); // update jobs table
        }

        db.CommitTransaction();
        db_success = true;
    } catch (const Exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while fetching particle parameters" << endl;
    } catch (const exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed while fetching particle parameters" << endl;
    }

    return db_success;
}


bool AbcSmc::update_particle_metrics(sqdb::Db &db, vector<string> &update_metrics_strings, vector<string> &update_jobs_strings) {
    bool db_success = false;

    try {
        db.Query("BEGIN EXCLUSIVE;").Next();
        for (size_t i = 0; i < update_metrics_strings.size(); ++i) {
            db.Query(update_metrics_strings[i].c_str()).Next(); // update metrics table
            db.Query(update_jobs_strings[i].c_str()).Next(); // update jobs table
        }

        db_success = true;
        db.CommitTransaction();
    } catch (const Exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while updating metrics:" << endl;
        for (size_t i = 0; i < update_metrics_strings.size(); ++i) {
            cerr << update_metrics_strings[i] << endl;
            cerr << update_jobs_strings[i] << endl;
        }
    } catch (const exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed while updating metrics:" << endl;
        for (size_t i = 0; i < update_metrics_strings.size(); ++i) {
            cerr << update_metrics_strings[i] << endl;
            cerr << update_jobs_strings[i] << endl;
        }
    }

    return db_success;
}
