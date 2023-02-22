
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <string_view>

#include "pls.h"
#include "RunningStat.h"
#include "AbcSmc.h"
#include "refsql.h"
#include "AbcMPI.h"
#include "AbcLog.h"

// need a positive int that is very unlikely
// to be less than the number of particles
#define STOP_TAG 10000000

using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;
using std::setw;
using std::cerr;

using namespace sqdb;
using namespace ABC;
using namespace chrono;

const int PREC = 5;
const string JOB_TABLE  = "job";
const string MET_TABLE  = "met";
const string MET_UNDER  = MET_TABLE + "_vals";
const string MET_REF  = MET_TABLE + "_name";

const string PAR_TABLE  = "par";
const string PAR_UNDER  = PAR_TABLE + "_vals";
const string PAR_REF  = PAR_TABLE + "_name";

const string UPAR_TABLE = "upar";
const string UPAR_UNDER = UPAR_TABLE + "_vals";
const string SEED_TABLE = "seed";

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++++++++ HELPER FUNS ++++++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

template<typename ITERABLE>
inline string join(const ITERABLE & v, const string & delim = ", ") {
    stringstream ss;
    for (auto it = v.begin(); it != v.end(); ++it) {
        if (it != v.begin()) ss << delim;
        ss << *it;
    }
    return ss.str();
}

template<typename ITERABLE>
inline string joineach(
    const vector<ITERABLE> & v,
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

vector<float> as_float_vector(Json::Value val, string key) { // not worth templating, despite appearances
    vector<float> extracted_vals;
    if ( val[key].isDouble() ) {
        extracted_vals.push_back( val[key].asFloat() );
    } else if ( val[key].isArray() ) {
        for (Json::Value jv : val[key]) extracted_vals.push_back( jv.asFloat() );
    } else {
        cerr << "Unfamiliar value type associated with " << key << " in configuration file: expecting floats or array of floats." << endl;
        exit(-216);
    }
    return extracted_vals;
}


vector<int> as_int_vector(Json::Value val, string key) {
    vector<int> extracted_vals;
    if ( val[key].isInt() ) {
        extracted_vals.push_back( val[key].asInt() );
    } else if ( val[key].isArray() ) {
        for (Json::Value jv : val[key]) extracted_vals.push_back( jv.asInt() );
    } else {
        cerr << "Unfamiliar value type associated with " << key << " in configuration file: expecting ints or array of ints." << endl;
        exit(-216);
    }
    return extracted_vals;
}


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ CONFIG PARSING +++++++++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //


void AbcSmc::process_predictive_prior_arguments(Json::Value par) {
    int arg_ct = par.isMember("predictive_prior_fraction") + par.isMember("predictive_prior_size");
    if (arg_ct == 0) {
        cerr << "Error: either predictive_prior_fraction or predictive_prior_size must be specified in configuration file." << endl;
        exit(1);
    } else if (arg_ct == 2) {
        cerr << "Error: only one of predictive_prior_fraction and predictive_prior_size may be specified in configuration file." << endl;
        exit(1);
    } else if (par.isMember("predictive_prior_fraction")) {
        vector<float> ppfs = as_float_vector(par, "predictive_prior_fraction");
        if (ppfs.size() > 1 and _smc_set_sizes.size() > 1 and ppfs.size() != _smc_set_sizes.size()) {
            cerr << "Error: If num_samples and predictive_prior_fraction both have length > 1 in configuration file, they must be equal in length." << endl;
            exit(1);
        }
        vector<int> set_sizes_copy = _smc_set_sizes;
        const int max_set = max(ppfs.size(), set_sizes_copy.size());
        ppfs.resize(max_set, ppfs.back());
        set_sizes_copy.resize(max_set, set_sizes_copy.back());
        _predictive_prior_sizes.clear();
        for (size_t i = 0; i < ppfs.size(); ++i) {
            if (ppfs[i] <= 0 or ppfs[i] > 1) {
                cerr << "Error: predictive_prior_fraction in configuration file must be > 0 and <= 1" << endl;
                exit(1);
            }
            _predictive_prior_sizes.push_back( round(ppfs[i] * set_sizes_copy[i]) );
        }
    } else if (par.isMember("predictive_prior_size")) {
        _predictive_prior_sizes = as_int_vector(par, "predictive_prior_size");
        if (_predictive_prior_sizes.size() > 1 and _smc_set_sizes.size() > 1 and _predictive_prior_sizes.size() != _smc_set_sizes.size()) {
            cerr << "Error: If num_samples and predictive_prior_size both have length > 1 in configuration file, they must be equal in length." << endl;
            exit(1);
        }
        const size_t max_set = max(_predictive_prior_sizes.size(), _smc_set_sizes.size());
        for (size_t i = 0; i < max_set; ++i) {
            if (get_pred_prior_size(i, QUIET) > get_num_particles(i, QUIET)) {
                cerr << "Error: requested predictive prior size is greater than requested SMC set size for at least one set in configuration file." << endl;
                exit(1);
            }
        }
    }
}


bool AbcSmc::parse_config(string conf_filename) {
    if (not file_exists(conf_filename.c_str())) {
        cerr << "File does not exist: " << conf_filename << endl;
        exit(1);
    }
    // TODO - Make sure any existing database actually reflects what is expected in JSON, particularly that par and met tables are legit
    Json::Value par;   // will contain the par value after parsing.
    Json::Reader reader;
    string json_data = slurp(conf_filename);

    bool parsingSuccessful = reader.parse( json_data, par );
    if ( !parsingSuccessful ) {
        // report to the user the failure and their locations in the document.
        cerr << "Failed to parse configuration\n" << reader.getFormattedErrorMessages();
        exit(1);
    }

    // TODO these should be mutually exclusive options
    std::string executable = par.get("executable", "").asString();
    if (executable != "") { set_executable( executable ); }
    std::string sharedobj = par.get("shared", "").asString();
    if (sharedobj != "") { set_simulator( sharedobj ); }

    string resume_dir = par.get("resume_directory", "").asString();
    if (resume_dir != "") {
        if (_mp->mpi_rank == mpi_root) cerr << "Resuming in directory: " << resume_dir << endl;
        set_resume_directory( resume_dir );
        set_resume( true );
    }

    set_smc_iterations( par["smc_iterations"].asInt() ); // TODO: or have it test for convergence
    _smc_set_sizes = as_int_vector(par, "num_samples");
    process_predictive_prior_arguments(par);
    // TODO--allow specification of pred prior size (single value or list of values)
    //set_predictive_prior_fraction( par["predictive_prior_fraction"].asFloat() );
    set_pls_validation_training_fraction( par["pls_training_fraction"].asFloat() ); // fraction of runs to use for training
    set_database_filename( par["database_filename"].asString() );
    // are we going to have particles that use a posterior from an earlier ABC run
    // to determine some of the parameter values?
    set_posterior_database_filename( par.get("posterior_database_filename", "").asString() );
    set_retain_posterior_rank( par.get("retain_posterior_rank", "false").asString() );
    if (_posterior_database_filename != "" and _num_smc_sets > 1) {
        cerr << "Using a posterior database as input is not currently supported with smc_iterations > 1. Aborting." << endl;
        exit(-203);
    }

    string noise = par.get("noise", "INDEPENDENT").asString();
    use_mvn_noise = (noise == "MULTIVARIATE");
    if (noise != "INDEPENDENT" and noise != "MULTIVARIATE") {
        cerr << "Unknown parameter noise type specified: " << noise << ". Aborting." << endl;
        exit(-210);
    }

    // Parse model parameters
    const Json::Value model_par = par["parameters"];
    map<string, int> par_name_idx;
    // TODO find a way to size_t this
    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Build name lookup
        string name = model_par[i]["name"].asString();
        assert(par_name_idx.count(name) == 0);
        par_name_idx[name] = i;
    }

    for ( Json::Value mpar : model_par )  {// Iterates over the sequence elements.
        string name = mpar["name"].asString();
        string short_name = mpar.get("short_name", "").asString();

        PriorType ptype = UNIFORM;
        string ptype_str = mpar["dist_type"].asString();
        if (ptype_str == "UNIFORM") {
            ptype = UNIFORM;
        } else if (ptype_str == "NORMAL" or ptype_str == "GAUSSIAN") {
            ptype = NORMAL;
        } else if (ptype_str == "PSEUDO") {
            ptype = PSEUDO;
        } else if (ptype_str == "POSTERIOR") {
            ptype = POSTERIOR;
            if (_posterior_database_filename == "") {
                cerr << "Parameter specfied as type POSTERIOR, without previously specifying a posterior_database_filename.  Aborting." << endl;
                exit(-204);
            }
        } else {
            cerr << "Unknown parameter distribution type: " << ptype_str << ".  Aborting." << endl;
            exit(-205);
        }

        NumericType ntype = INT;
        string ntype_str = mpar["num_type"].asString();
        if (ntype_str == "INT") {
            ntype = INT;
        } else if (ntype_str == "FLOAT") {
            ntype = FLOAT;
        } else {
            cerr << "Unknown parameter numeric type: " << ntype_str << ".  Aborting." << endl;
            exit(-206);
        }

        double (*_untransform_func)(const double);
        // linear rescaling [min, max] not for par as sampled, but as input to sim
        // NB: if par is not on [0, 1] after untransforming, this is still a linear rescaling, but not onto [min, max]
        pair<double, double> par_rescale = {0.0, 1.0};
        map<string, vector<int> > mod_map { {"transformed_addend", {}}, {"transformed_factor", {}}, {"untransformed_addend", {}}, {"untransformed_factor", {}} };
        //auto _untransform = [](const double t) { return t; };
        if (not mpar.isMember("untransform")) {
            _untransform_func = [](const double t) { return t; };
        } else if (mpar["untransform"].type() == Json::ValueType::stringValue) {
            string ttype_str = mpar.get("untransform", "NONE").asString();
            if (ttype_str == "NONE") { // TODO - it's possible this may not actually ever be called
                _untransform_func = [](const double t) { return t; };
                //ttype = UNTRANSFORMED;
            } else if (ttype_str == "POW_10") {
                _untransform_func = [](const double t) { return pow(10.0, t); };
                //ttype = LOG_10;
                use_transformed_pars = true;
            } else if (ttype_str == "LOGISTIC") {
                _untransform_func = [](const double t) { return ABC::logistic(t); };
                //ttype = LOGIT;
                use_transformed_pars = true;
            } else {
                cerr << "Unknown parameter transformation type: " << ttype_str << ".  Aborting." << endl;
                exit(-206);
            }

        } else if (mpar["untransform"].type() == Json::ValueType::objectValue) {
            Json::Value untransform = mpar["untransform"];
            string ttype_str = untransform["type"].asString();
            if (ttype_str != "LOGISTIC") {
                cerr << "Only type: LOGISTIC is currently supported for untransformation objects.  (NONE and POW_10 supported as untransformation strings.)\n";
                exit(-207);
            }
            par_rescale = {untransform["min"].asDouble(), untransform["max"].asDouble()};
            _untransform_func = [](const double t) { return ABC::logistic(t); };
            use_transformed_pars = true;
            //Json::ValueType mod_type = untransform["transformed_addend"].type();

            for (auto& mod_type: mod_map) {
                if (untransform.isMember(mod_type.first)) {
                    for (auto json_val: untransform[mod_type.first]) mod_type.second.push_back(par_name_idx[json_val.asString()]);
                }
            }
        } else {
            cerr << "Unsupported JSON data type associated with 'untransform' parameter key.\n";
            exit(-208);
        }

        double par1 = mpar["par1"].asDouble();
        double par2 = mpar["par2"].asDouble();
        double step = mpar.get("step", 1.0).asDouble(); // default increment is 1

        add_next_parameter(name, short_name, ptype, ntype, par1, par2, step, _untransform_func, par_rescale, mod_map);
    }

    // Parse model metrics
    const Json::Value model_met = par["metrics"];
    for ( Json::Value mmet : model_met )  {// Iterates over the sequence elements.
        string name = mmet["name"].asString();
        string short_name = mmet.get("short_name", "").asString();

        NumericType ntype = INT;
        string ntype_str = mmet["num_type"].asString();
        if (ntype_str == "INT") {
            ntype = INT;
        } else if (ntype_str == "FLOAT") {
            ntype = FLOAT;
        } else {
            cerr << "Unknown metric numeric type: " << ntype_str << ".  Aborting." << endl;
            exit(-209);
        }

        double val = mmet["value"].asDouble();

        add_next_metric(name, short_name, ntype, val);
    }

    return true;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ PARAMETER FUNCTIONS ++++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

Row AbcSmc::do_complicated_untransformations(vector<Parameter*>& _model_pars, Row& pars) {
    assert( _model_pars.size() == npar() );
    assert( static_cast<size_t>(pars.size()) == npar() );
    const vector<double> identities = {0.0, 1.0, 0.0, 1.0};
    Row upars(npar());
    for (size_t i = 0; i < npar(); ++i) {
//cerr << "Parameter " << i << ": " << _model_pars[i]->get_name() << endl;
        const Parameter* mpar = _model_pars[i];
        vector<double> modifiers(identities); // TODO -- double check that this is a legit copy constructor
        map<string, vector<int> > mod_map = mpar->get_par_modification_map();
        for (auto which : mod_map["transformed_addend"]) { modifiers[0] += pars[which]; }
        for (auto which : mod_map["transformed_factor"]) { modifiers[1] *= pars[which]; }
        for (auto which : mod_map["untransformed_addend"]) { modifiers[2] += pars[which]; }
        for (auto which : mod_map["untransformed_factor"]) { modifiers[3] *= pars[which]; }
        upars[i] = mpar->untransform(pars[i], modifiers);
    }
    return upars;
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ CORE ALGORITHM STEPS +++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

bool AbcSmc::_run_simulator(Row &par, Row &met, const unsigned long int rng_seed, const unsigned long int serial) {
    vector<float_type> met_vec = (*_simulator)( as_vector(par), rng_seed, serial, _mp );
    bool particle_success = (met_vec.size() == nmet());
    if (!particle_success) {
        cerr << "ERROR: simulator function returned the wrong number of metrics: expected " << nmet() << ", received " << met_vec.size() << endl;
    }
    met = as_row(met_vec);
    return particle_success;
}

void AbcSmc::read_SMC_complete(sqdb::Db &db, const bool all) {

    if (not all) { // ask smc_set

    } else {

    }
    
};

void AbcSmc::rank_SMC_last() {

};

void AbcSmc::summarize_SMC() {

};

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ DATABASE FUNCTIONS +++++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

bool _db_execute_strings(sqdb::Db &db, vector<string> &update_buffer) {
    bool db_success = false;
    try {
        db.Query("BEGIN EXCLUSIVE;").Next();
        for (auto buff : update_buffer) { _db_execute(db, buff); }
        db_success = true;
        db.CommitTransaction();
    } catch (const Exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        for (size_t i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    } catch (const exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        for (size_t i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    }
    return db_success;
}

bool _db_execute(sqdb::Db &db, const char * query, const bool verbose) {
    bool db_success = true;
    if (verbose) {
        cerr << query << endl;
    }
    try {
        while (*query) {
            if (verbose) { cerr << string(query).substr(0, 10) << " ... " << endl; }
            db.Query(query, &query).Next();
        }
    } catch (const Exception& e) {
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorCode() << " : " << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        cerr << query << endl;
        db_success = false;
    } catch (const exception& e) {
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        cerr << query << endl;
        db_success = false;
    }
    return db_success;
};

bool _db_execute(sqdb::Db &db, const std::string &ss, const bool verbose) {
    return _db_execute(db, ss.c_str(), verbose);
};

bool _db_execute(sqdb::Db &db, stringstream &ss, const bool verbose) {
    bool db_success = _db_execute(db, ss.str(), verbose);
    ss.str(string());
    ss.clear();
    return db_success;
};

bool _db_tables_exist(sqdb::Db &db, vector<string> table_names) {
    // Note that retval here is whether tables exist, rather than whether
    // db transaction was successful.  A failed transaction will throw
    // an exception and exit.
    bool tables_exist = true;
    try {
        for(string table_name: table_names) {
            string query_str = "SELECT COUNT(*) FROM sqlite_master WHERE type IN ('table', 'view') and name='" + table_name + "';";
            //cerr << "Attempting: " << query_str << endl;
            Statement s = db.Query( query_str.c_str() );
            s.Next();
            const int count = s.GetField(0);
            if (count < 1) {
                cerr << "Table " << table_name << " does not exist in database.\n";
                tables_exist = false;
            }
        }
    } catch (const Exception& e) {
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while checking whether the following tables exist:";
        for(string table_name: table_names) cerr << " " << table_name;
        cerr << endl;
        exit(-212);
    } catch (const exception& e) {
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed while checking whether the following tables exist:";
        for(string table_name: table_names) cerr << " " << table_name;
        cerr << endl;
        exit(-213);
    }
    return tables_exist;
}


// Build DB if it doesn't exist;
// If it does exist but more sets are needed, filter particles and sample for next set;
// If the specified number of sets already exist, exit gracefully
bool AbcSmc::process_database(const gsl_rng* RNG) {

    if (build_database(RNG)) return true; // if DB doesn't exist, create it and exit

    sqdb::Db db(_database_filename.c_str());
    _particle_parameters.clear();
    _particle_metrics.clear();
    _weights.clear();
    _predictive_prior.clear();

    vector< vector<int> > serials;
    if (not read_SMC_sets_from_database(db, serials)) return false; // slurp sets & do particle filtering (identify pred prior) if needed
    const size_t next_set = serials.size();
    assert(next_set > 0);
    const size_t last_set = next_set - 1; // this set number

    AbcLog::report_convergence_data(this, last_set);

    cerr << endl << endl;

    if (_num_smc_sets > next_set) {

        db.Query("BEGIN EXCLUSIVE;").Next();
        //ss << "insert into sets values ( 0, 'Q'"; for (int j = 0; j < npar(); j++) ss << ", NULL"; ss << ");";
        //_db_execute(db, ss);

        stringstream ss;
        Row pars;
        const int last_serial = serials.back().back();
        gsl_matrix* L = setup_mvn_sampler(next_set);

        string noise_type = use_mvn_noise ? "MULTIVARIATE" : "INDEPENDENT";
        cerr << "Populating next set using " << noise_type << " noising of parameters.\n";
        const size_t num_particles = get_num_particles(next_set);
        for (size_t i = 0; i < num_particles; i++) {
            const int serial = last_serial + 1 + i;

            if (use_mvn_noise) {
                pars = sample_mvn_predictive_priors(next_set, RNG, L);
            } else {
                pars = sample_predictive_priors(next_set, RNG);
            }
            QueryStr qstr;
// TODO here's the place to fix the NEXT SET problem
            ss << "insert into " << JOB_TABLE << " values ( " << serial << ", "
                                               << next_set << ", "
                                               << i << ", "
                                               << time(NULL)
                                               << ", NULL, 'Q', -1, 0 );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute(db, ss);

            const unsigned long int seed = gsl_rng_get(RNG); // seed for particle
            ss << "insert into " << PAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << pars[j]; ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute(db, ss);

            if (use_transformed_pars) {
                const Row upars = do_complicated_untransformations(_model_pars, pars);
                ss << "insert into " << UPAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << upars[j]; ss << " );";
                //cerr << "attempting: " << ss.str() << endl;
                _db_execute(db, ss);
            }

            ss << "insert into " << MET_TABLE << " values ( " << serial; for (size_t j = 0; j < nmet(); j++) ss << ", NULL"; ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute(db, ss);
        }
        db.CommitTransaction();
        gsl_matrix_free(L);
    } else {
        cerr << "Database already contains " << _num_smc_sets << " complete sets.\n";
    }

    return true;
}

// Read in existing sets and do particle filtering as appropriate
bool AbcSmc::read_SMC_sets_from_database (sqdb::Db &db, vector< vector<int> > &serials) {
    // make sure database looks intact
    if ( not _db_tables_exist(db, {JOB_TABLE, PAR_TABLE, MET_TABLE}) ) {
        cerr << "ERROR: Failed to read SMC set from database because one or more tables are missing.\n";
        return false;
    }

    // make sure set t is a completed set
    QueryStr qstr;
    Statement s = db.Query(qstr.Format("SELECT smcSet, total, run FROM smc_summary ORDER BY smcSet;"));

    serials.clear();
    _particle_parameters.clear();
    _particle_metrics.clear();
    bool set_size_test_success = true;

    while (s.Next()) {
        int t = s.GetField(0);
        int set_size = s.GetField(1);
        int completed_set_size = s.GetField(2);
        if (set_size != completed_set_size) {
            cerr << "ERROR: Failed to read SMC set from database because not all particles are complete in set " << t << "\n";
            return false;
        }
        const int json_set_size = get_num_particles(t, QUIET);
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

        QueryStr qstr;

        auto q = qstr.Format(SQDB_MAKE_TEXT(R"SQL(
            SELECT serial, particleIdx, posterior, P.*, M.* FROM %s
            JOIN %s USING(serial)
            JOIN %s USING(serial)
            WHERE smcSet == %d;
        )SQL"), JOB_TABLE.c_str(), PAR_TABLE.c_str(), MET_TABLE.c_str(), t);

        // join all three tables for rows with smcSet = t, slurp and store values
        string select_str = "select J.serial, J.particleIdx, J.posterior, " + _build_sql_select_par_string("") + ", " + _build_sql_select_met_string()
                            + "from " + JOB_TABLE + " J, " + MET_TABLE + " M, " + PAR_TABLE + " P where J.serial = M.serial and J.serial = P.serial "
                            + "and J.smcSet = " + to_string((long long) t) + ";";

        serials.push_back( vector<int>(completed_set_size) );

        Statement s2 = db.Query( select_str.c_str() );

        int particle_counter = 0;
        vector<pair<int, int> > posterior_pairs;
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
            //set_next_predictive_prior_size(t, _particle_parameters[t].size());
            const size_t next_pred_prior_size = get_pred_prior_size(t);
            _predictive_prior.push_back( vector<size_t>(next_pred_prior_size) );
            if (use_pls_filtering) {
                _filter_particles( t, _particle_metrics[t], _particle_parameters[t], next_pred_prior_size );
            } else {
                _filter_particles_simple( t, _particle_metrics[t], _particle_parameters[t], next_pred_prior_size );
            }
            auto posterior_pars = _particle_parameters[t](_predictive_prior[t], Eigen::all);
            auto posterior_mets = _particle_metrics[t](_predictive_prior[t], Eigen::all);
            AbcLog::filtering_report(this, t, posterior_pars, posterior_mets);

            vector<string> update_strings(next_pred_prior_size);
            for (size_t i = 0; i < next_pred_prior_size; i++) { // best to worst performing particle in posterior?
                const int particle_idx = _predictive_prior[t][i];
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


// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ MPI SUPPORT FUNCTIONS ++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

void AbcSmc::_particle_scheduler(const size_t t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {

    auto _num_particles = get_num_particles(t);

    // sample parameter distributions; copy values into Y matrix and into send_data buffer
    Row par_row;
    for (size_t i = 0; i < _num_particles; i++) {
        if (t == 0) { // sample priors
            Mat2D posterior = slurp_posterior();
            int posterior_rank = -1;
            Y_orig.row(i) = sample_priors(RNG, posterior, posterior_rank);
        } else {      // sample predictive priors
            Y_orig.row(i) = sample_predictive_priors(t, RNG);
        }
    }

    X_orig.resize(_num_particles, nmet());

    ABC::particle_scheduler(X_orig, Y_orig, _mp);

}


void AbcSmc::_particle_worker(const size_t seed, const size_t serial) {
    ABC::particle_worker(npar(), nmet(), _simulator, seed, serial, _mp);
}

// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //
// ++++++++++++++++++++ CORE ALGORITHM LOOP ++++++++++++++++++++++++++++++++++++ //
// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++ //

// 0. Setup (= construct storage from scratch - not repeated, even across invocations)
// 0.1 Populate (= Sample Priors + Generate Pseudos - not repeated, even across invocations)
// 1. Simulate (= read parameters, run sim, write metrics)
// 2. Process (= do PLS => new params)
// 3. GOTO 1.

// open the SQL database, write the static structure, write the needed simulation-specific structure,
// then write the conditional-on-specifics-but-static elements of the database
bool setup(
    const std::string dbfile,
    const std::vector<std::string> &par_names,
    const std::vector<std::string> &met_names
) {
    // NB: the cast here is to do with c-nature of sqlite3 library
    // we're doing this to make the `eval` command available in sqlite3 when we open the database
    // this is only needed when setting up the db
    sqlite3_auto_extension((void(*)())sqlite3_eval_init);

    // assert: db empty
    // TODO: error (warn? allow with overwrite flag?) if db non-empty
    sqdb::Db db(dbfile.c_str());

    stringstream ss;
    db.BeginTransaction();

    // set up the static elements of an abcsmc database
    ss << bin2c_sqlviews_sql;

    _db_execute(db, ss);

    QueryStr qs;
    ss << qs.Format(SQDB_MAKE_TEXT("INSERT INTO %s (name) VALUES ('%s');"), PAR_REF.c_str(), join(par_names, "'), ('").c_str());
    ss << qs.Format(SQDB_MAKE_TEXT("INSERT INTO %s (name) VALUES ('%s');"), MET_REF.c_str(), join(met_names, "'), ('").c_str());
    _db_execute(db, ss);

    // apply the dynamically generated elements
    ss << bin2c_sqldynamic_sql;
    _db_execute(db, ss);

    db.CommitTransaction();
        
    return true;
};

// converts a series of n serials + n value sets of k cols into long format for insertion into sql
// e.g. (serial1, col1, value1), (serial1, col2, value2), ..., (serialn, colk, valuek)
//
// TODO: most typical use case is an iota of serials - might be worth providing default / additional
// arguments / polymorphism in that direction
//
// TODO: for more generic storage, might make more sense to build an intermediate iterator
// of (serial, col, value) tuples, then use that to build the string
// that would provide a more general step that might be useable for other storage types that still
// want long data, but don't use the same serialization format as SQL(ite)
template <typename SerialIt, typename ValRowIt>
std::string serial_id_values_stream(
    const SerialIt & serials, // iterable of serials
    const ValRowIt & vals // iterable of "Rows" of values (something that has size + iterator) 
) {
    assert(serials.size() == vals.size());
    std::stringstream ss("(");
    
    for (size_t i = 0; i < serials.size(); ++i) {
        for (size_t j = 0; j < vals[i].size(); ++j) {
            if (j != 0) { ss << "), ("; };
            ss << serials[i] << ", " << j + 1 << ", " << vals[i][j];
        }
    }

    ss << ")";

    return ss.str();
};

// similar to serial_id_values_stream, but for a single value per serial - i.e. no id key
template <typename SerialIt, typename ValIt>
std::string serial_values_stream(
    const SerialIt & serials, // iterable of serials
    const ValIt & vals // iterable of values 
) {
    assert(serials.size() == vals.size());
    std::stringstream ss("(");
    
    for (size_t i = 0; i < serials.size(); ++i) {
        if (i != 0) { ss << "), ("; };
        ss << serials[i] << ", " << vals[i];
    }

    ss << ")";

    return ss.str();
};

bool populate(
    sqdb::Db &db,
    const std::vector<Row> &pars,
    const size_t smcSet = -1,
    const std::vector<size_t> &seeds = {},
    const std::vector<int> &ranks = {},
    const std::vector<Row> &upars = {}
) {
    // inserting pars.size() jobs into the db; assert matching sizes for seeds, ranks, upars
    assert((upars.size() == 0) or (upars.size() == pars.size()));
    assert((ranks.size() == 0) or (ranks.size() == pars.size()));
    assert((seeds.size() == 0) or (seeds.size() == pars.size()));

    QueryStr qstr;

// TODO: deal with alternative particleIdx? -- might solve this problem with posterior table, job becomes a view
// desire to insert a bunch of rows with some fixed values, otherwise using defaults
// https://stackoverflow.com/questions/65439177/how-to-add-specific-number-of-empty-rows-in-sqlite
// bulk insert into the job table
    db.Query(qstr.Format(R"SQL(
        WITH RECURSIVE g(n) AS ( -- this creates a temporary table, column n, values 0, 1, ..., pars.size() - 1
            VALUES (0)
            UNION ALL
            SELECT n + 1
            FROM g WHERE n < %d
        )
        INSERT INTO %s (smcSet, particleIdx) -- then inserts those temporary values into the job table
        SELECT %d, n FROM g;
    )SQL", pars.size(), JOB_TABLE.c_str(), smcSet
    )).Next();

    // now collect the newly created serials
    auto s = db.Query(qstr.Format(SQDB_MAKE_TEXT(
        "SELECT last_insert_rowid() FROM %s;"
    ), JOB_TABLE.c_str()));
    s.Next();
    const int last_serial = s.GetField(0);
    const int first_serial = last_serial - pars.size() + 1;
    std::vector<int> serials(pars.size());
    std::iota(serials.begin(), serials.end(), first_serial);

    // if we have ranks, update the job table with them
    if (ranks.size() > 0) {
        // bulk insert into the job table
        auto insert = serial_values_stream(serials, ranks);
        db.Query(qstr.Format(R"SQL(
            WITH rankdata (rserial, rank) AS (VALUES %s) -- create a temporary table of serial => rank
            UPDATE %s SET posterior = rank FROM rankdata -- update the job table with the rank
            WHERE serial = rserial; -- where the serials match
        )SQL", insert.c_str(), JOB_TABLE.c_str())).Next();
    }

    // bulk insert into the par table
    auto insert = serial_id_values_stream(serials, pars);
    db.Query(qstr.Format(SQDB_MAKE_TEXT(
        "INSERT INTO %s (serial, parIdx, value) VALUES %s ON CONFLICT DO NOTHING;"
    ), PAR_UNDER.c_str(), insert.c_str())).Next();

    // if provided seeds, bulk insert into the seed table
    if (seeds.size() > 0) {
        db.Query(qstr.Format(SQDB_MAKE_TEXT(
            "INSERT INTO %s (seed) VALUES (%s);"
        ), SEED_TABLE.c_str(), join(seeds, "), ("))).Next();
    }

    // if we have u(ntransformed )pars, bulk insert into the upar table
    if (upars.size() > 0) {
        // n.b. this does not yet allow for upars being a different size than pars
        // which might be a useful way to only xform a subset of the parameters
        insert = serial_id_values_stream(serials, upars);
        db.Query(qstr.Format(SQDB_MAKE_TEXT(
            "INSERT INTO %s (serial, parIdx, value) VALUES %s ON CONFLICT DO NOTHING;"
        ), UPAR_UNDER.c_str(), insert.c_str())).Next();
    }

    return true;
};

bool AbcSmc::_sample_priors(
    const gsl_rng* RNG, const size_t n, // number of samples
    vector<size_t> &seeds, // will be populated with seeds? will be pre-filled with seeds?
    vector<Row> &pars, // will be populated with parameters
    vector<Row> &upars, // will be populated with upars (if relevant)
    vector<int> &ranks // will be populated with ranks (if relevant)
) {

    seeds.resize(n); pars.resize(n);
    if (use_transformed_pars) { upars.resize(n); }
    if (_retain_posterior_rank) { ranks.resize(n); }

    Mat2D posterior = slurp_posterior();
    
    for (size_t i = 0; i < n; ++i) {
        seeds[i] = gsl_rng_get(RNG);
        pars[i] = Row::Zero(_model_pars.size());
        vector<int> posterior_indices;

        bool increment_nonrandom_par = true; // only one PSEUDO parameter gets incremented each time
        bool increment_posterior = true;     // posterior parameters get incremented together, when all pseudo pars reach max val
        int posterior_rank = -1;

        for (size_t j = 0; j < _model_pars.size(); j++) {
            Parameter* p = _model_pars[j];
            float_type val;
            // if it's a non-random, PSEUDO parameter
            if (p->get_prior_type() == PSEUDO) {
                // get the state now, so that the first time it will have the initialized value
                val = (float_type) p->get_state();
                // We need to imitate the way nested loops work, but in a single loop.
                // PSEUDO parameters only get incremented when any and all previous PSEUDO parameters
                // have reached their max values and are being reset
                if (increment_nonrandom_par) {
                    // This parameter has reached it's max value and gets reset to minimum
                    // Because of possible floating point errors, we check whether the state is within step/10,000 from the max
                    const float_type PSEUDO_EPSILON = 0.0001 * p->get_step();
                    if (p->get_state() + PSEUDO_EPSILON >= p->get_prior_max()) {
                        p->reset_state();
                    // otherwise, increment this one and prevent others from being incremented
                    } else {
                        p->increment_state();
                        increment_nonrandom_par = false;
                        increment_posterior     = false;
                    }
                }
            } else if (p->get_prior_type() == POSTERIOR) {
                val = 0; // will be replaced later in function
                if (posterior_indices.size() == 0) {
                    posterior_rank = (int) p->get_state();
                } else {
                    // require that posterior pars be synchronized
                    assert(posterior_rank == (int) p->get_state());
                }
                posterior_indices.push_back(j);
            } else {
                // Random parameters get sampled independently from each other, and are therefore easy
                val = p->sample(RNG);
            }

            pars[i](j) = val;
        }

        if (_posterior_database_filename != "") {

            for (size_t c = 0; c < posterior_indices.size(); ++c) {
                pars[i](posterior_indices[c]) = posterior(posterior_rank, c);
                Parameter* p = _model_pars[posterior_indices[c]];

                if (increment_posterior) {
                    if (p->get_state() >= p->get_prior_max()) {
                        p->reset_state();
                    } else {
                        p->increment_state();
                    }
                }

            }

        }

        if (_retain_posterior_rank) { ranks[i] = posterior_rank; }
        if (use_transformed_pars) { upars[i] = do_complicated_untransformations(_model_pars, pars[i]); }

    }

    return true;

};

// this is now mostly ready to have a whatever-storage backend
bool AbcSmc::build_database(const gsl_rng* RNG) {
    if (_posterior_database_filename != "" and not file_exists(_posterior_database_filename.c_str())) {
        std::cerr << "ERROR: Cannot read in posterior database: " << _posterior_database_filename << " does not exist\n";
        exit(-215);
    }

    sqdb::Db db(_database_filename.c_str());

    stringstream ss;
    if ( !_db_tables_exist(db, {JOB_TABLE}) and !_db_tables_exist(db, {PAR_TABLE}) and !_db_tables_exist(db, {MET_TABLE}) ) {
        setup(_database_filename, parameter_short_names(), metric_short_names(), use_transformed_pars);
    } else {
        return false;
    }

    // no-op if posterior database is not specified
    Mat2D posterior = slurp_posterior();

    const size_t num_particles = get_num_particles(0);

    vector<size_t> seeds; // will be populated with seeds? will be pre-filled with seeds?
    vector<Row> pars; // will be populated with parameters
    vector<Row> upars; // will be populated with upars (if relevant)
    vector<int> ranks; // will be populated with ranks (if relevant)

    _sample_priors(RNG, num_particles, seeds, pars, upars, ranks);

    db.BeginTransaction();
    populate(db, pars, 0, seeds, ranks, upars);
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
    } catch (const Exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while fetching particle parameters" << endl;
    } catch (const exception& e) {
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
        db.BeginTransaction();
        for (size_t i = 0; i < update_metrics_strings.size(); ++i) {
            _db_execute(db, update_metrics_strings[i]);
            _db_execute(db, update_jobs_strings[i]);
        }
        db_success = true;
        db.CommitTransaction();
    } catch (const Exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while updating metrics:" << endl;
        for (size_t i = 0; i < update_metrics_strings.size(); ++i) {
            cerr << update_metrics_strings[i] << endl;
            cerr << update_jobs_strings[i] << endl;
        }
    } catch (const exception& e) {
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

bool AbcSmc::simulate_particle_by_serial(const int serial_req) { return simulate_next_particles(1, serial_req, -1); }

bool AbcSmc::simulate_particle_by_posterior_idx(const int posterior_req) { return simulate_next_particles(1, -1, posterior_req); }

bool AbcSmc::simulate_next_particles(
    const int n, const int serial_req, const int posterior_req
) { // defaults are 1, -1, -1
    bool verbose = n == 1;
    assert(n == 1 or (serial_req == -1 and posterior_req == -1));
    assert(serial_req == -1 or posterior_req == -1);
//bool AbcSmc::simulate_database(const int smc_set, const int particle_id) {
    sqdb::Db db(_database_filename.c_str());
    string model_par_table = _db_tables_exist(db, {UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
    vector<Row> par_mat;  //( n, Row(npar()) ); -- expected size, if n rows are available
    vector<Row> met_mat;  //( n, Row(nmet()) );

    stringstream select_ss;
    string limit = n == -1 ? "" : "limit " + to_string(n);
    select_ss << "select J.serial, P.seed, " << _build_sql_select_par_string("");
    select_ss << "from " << model_par_table << " P, " << JOB_TABLE << " J where P.serial = J.serial ";
    // Do already running jobs as well, if there are not enough queued jobs
    // This is because we are seeing jobs fail/time out for extrinsic reasons on the stuporcomputer
    if (serial_req > -1) {
        select_ss << "and J.serial = " << serial_req << ";";
    } else if (posterior_req > -1) {
        select_ss << "and smcSet = (select max(smcSet) from job where posterior > -1) and posterior = " << posterior_req << ";";
    } else {
        select_ss << "and (J.status = 'Q' or J.status = 'R') order by J.status, J.attempts " << limit << ";";
    }
    //  line below is much faster for very large dbs, but not all particles will get run e.g. if some particles are killed by scheduler
    //  select_ss << "and J.status = 'Q' limit " << n << ";";

    const size_t overall_start_time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
    // build jobs update statement to indicate job is running
    stringstream update_ss;
    update_ss << "update " << JOB_TABLE << " set startTime = " << overall_start_time << ", status = 'R', attempts = attempts + 1 where serial = "; // we don't know the serial yet

    vector<int> serials;
    vector<unsigned long int> rng_seeds;
    bool ok_to_continue = fetch_particle_parameters(db, select_ss, update_ss, serials, par_mat, rng_seeds, verbose);
    vector<string> update_metrics_strings;
    vector<string> update_jobs_strings;
    stringstream ss;
    if (ok_to_continue) {
        for (size_t i = 0; i < par_mat.size(); ++i) {
            const high_resolution_clock::time_point start_time = high_resolution_clock::now();
            const int serial = serials[i];
            met_mat.push_back( Row(nmet()) );
            bool success = _run_simulator(par_mat[i], met_mat[i], rng_seeds[i], serial);
            if (not success) exit(-211);

            stringstream ss;
            QueryStr qstr;
            for (size_t j = 0; j < nmet(); j++) {
                ss << qstr.Format(SQDB_MAKE_TEXT(
                    "INSERT INTO %s (serial, metIdx, value) VALUES (%d, %d, %f) ON CONFLICT (serial, metIdx) DO NOTHING;"
                ), MET_UNDER.c_str(), serial, j+1, met_mat[i][j]);
            }
            update_metrics_strings.push_back(ss.str());

            ss.str(string()); ss.clear();

            const size_t time_since_unix_epoch = duration_cast<seconds>(start_time.time_since_epoch()).count();
            const duration<size_t> time_span = duration_cast<duration<size_t>>(high_resolution_clock::now() - start_time); // duration in seconds
            // build jobs update statement to indicate job is running
            ss << qstr.Format(SQDB_MAKE_TEXT(
                "UPDATE %s SET startTime = %d, duration = %d, status = 'D' WHERE serial = %d AND (status = 'R' OR status = 'Q' OR status = 'P');"
            ), JOB_TABLE.c_str(), time_since_unix_epoch, time_span.count(), serial);
            update_jobs_strings.push_back(ss.str());
            ss.str(string()); ss.clear();
        }
        update_particle_metrics(db, update_metrics_strings, update_jobs_strings);
    } else {
        cerr << "Parameter selection from database failed.\n";
    }

    return true;
}

// template <typename Orderable>
// vector<size_t> ranked_sample(
//     const Orderable & values, // the values to be ranked 
//     const size_t samplen,
//     vector<size_t> & ranking
// ) {
//     assert(values.size() >= samplen);
//     ranking = ordered(distances);
//     return vector<size_t>(ranking.begin(), ranking.begin() + samplen);
// }

void AbcSmc::_set_predictive_prior (
    const int t,
    const int next_pred_prior_size,
    const Col& distances
) {
    vector<size_t> ranking = ordered(distances);
    vector<int> sample(ranking.begin(), ranking.begin() + next_pred_prior_size); // This is the predictive prior / posterior
    _predictive_prior[t] = sample;

}

void AbcSmc::_filter_particles_simple (int t, Mat2D &X_orig, Mat2D &Y_orig, int next_pred_prior_size) {
    // x is metrics, y is parameters
    Row X_sim_means, X_sim_stdev;
    Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
    Mat2D Y = colwise_z_scores( Y_orig );
    Row obs_met_z = ABC::z_transform_vals(_met_vals, X_sim_means, X_sim_stdev);

    Col distances  = euclidean(obs_met_z, X);
    _set_predictive_prior (t, next_pred_prior_size, distances);
}

// Run PLS
// Box-Cox transform data -- TODO?
//void test_bc( Mat2D );
//test_bc(Y_orig);
// @tjh TODO: obsolete? have updated the PLS constructor to just do the calculation if presented with X/Y/components/method
// as that seems to be the most common use case
PLS_Model& AbcSmc::run_PLS(
    Mat2D &X, // predictors; rows = sample, cols = variables
    Mat2D &Y, // responses; rows = sample, cols = variables
    const size_t pls_training_set_size, // number of samples from X and Y to use
    const size_t ncomp // number of PLS components to keep
) {

    // assert same rows to X and Y? or perhaps manage that in PLS?
    assert(pls_training_set_size <= static_cast<size_t>(X.rows())); // can't train against more observations than we have
    const size_t npred = X.cols();      // number of predictor variables
    const size_t nresp = Y.cols();      // number of response variables
    return PLS_Model(npred, nresp, ncomp).plsr(
        X.topRows(pls_training_set_size),
        Y.topRows(pls_training_set_size), KERNEL_TYPE1
    );
    
}

PLS_Model AbcSmc::_filter_particles (
    int t, Mat2D &X_orig, Mat2D &Y_orig, int next_pred_prior_size,
    const bool verbose
) {
    Row X_sim_means, X_sim_stdev;
    Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
    Mat2D Y = colwise_z_scores( Y_orig );
    Row obs_met = ABC::z_transform_vals(_met_vals, X_sim_means, X_sim_stdev);

    const size_t pls_training_set_size = round(X.rows() * _pls_training_fraction);
    // @tjh TODO -- I think this may be a bug, and that ncomp should be equal to number of predictor variables (metrics in this case), not reponse variables
    size_t ncomp = npar();             // It doesn't make sense to consider more components than model parameters

    PLS_Model plsm(X.topRows(pls_training_set_size), Y.topRows(pls_training_set_size), ncomp);

    // A is number of components to use
    if (verbose) { plsm.print_explained_variance(X, Y); }

    const int test_set_size = X.rows() - pls_training_set_size; // number of observations not in training set
    auto num_components = plsm.optimal_num_components(X.bottomRows(test_set_size), Y.bottomRows(test_set_size), NEW_DATA);
    size_t num_components_used = num_components.maxCoeff();
    if (verbose) {
        cerr << "Optimal number of components for each parameter (validation method == NEW DATA):\t" << num_components << endl;
        cerr << "Using " << num_components_used << " components." << endl;
    }

    // Calculate new, orthogonal metrics (==scores) using the pls model
    // Is casting as real always safe?
    Row   obs_scores = plsm.scores(obs_met, num_components_used).row(0).real();
    Mat2D sim_scores = plsm.scores(X, num_components_used).real();
    Col   distances  = euclidean(obs_scores, sim_scores);

    _set_predictive_prior(t, next_pred_prior_size, distances);

    return plsm;
}

Col AbcSmc::euclidean( Row obs_met, Mat2D sim_met ) {
    Col distances = Col::Zero(sim_met.rows());
    for (int r = 0; r<sim_met.rows(); r++) {
        for (int c = 0; c<sim_met.cols(); c++) {
            distances(r) += pow(obs_met(c) - sim_met(r,c), 2);
        }
        distances(r) = sqrt( distances(r) );
    }
    return distances;
}

Mat2D AbcSmc::slurp_posterior(const bool verbose) {
    // TODO - handle sql/sqlite errors
    // TODO - if "posterior" database doesn't actually have posterior values, this will fail silently

    if (_posterior_database_filename == "") {
        if (verbose) { cerr << "No posterior database specified, skipping." << endl; }
        return Mat2D();
    }

    // determine dimensions of results we're going to read in
    sqdb::Db post_db(_posterior_database_filename.c_str());
    Statement posterior_query = post_db.Query(("select count(*) from " + JOB_TABLE + " where posterior > -1;").c_str());
    posterior_query.Next();
    const int posterior_size = posterior_query.GetField(0); // num of rows

    int num_posterior_pars = 0; // num of cols
    string posterior_strings = "";
    for (Parameter* p: _model_pars) {
        if (p->get_prior_type() == POSTERIOR) {
            ++num_posterior_pars;
            if (posterior_strings != "") {
                posterior_strings += ", ";
            }
            posterior_strings += p->get_short_name();
        }
    }

    Mat2D posterior = Mat2D::Zero( posterior_size, num_posterior_pars );

    // Identify table to pull parameter values from
    string post_par_table = _db_tables_exist(post_db, {UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
    stringstream ss;
    ss << "select " << posterior_strings << " from " << post_par_table << " P, " << JOB_TABLE << " J where P.serial = J.serial and posterior > -1;";
    posterior_query = post_db.Query(ss.str().c_str());

    int r = 0;
    while(posterior_query.Next()) {
        for (int c = 0; c < num_posterior_pars; ++c) {
            posterior(r,c) = (double) posterior_query.GetField(c);
        }
        ++r;
    }

    return posterior;
}

void AbcSmc::calculate_doubled_variances( int t ) {
    vector<RunningStat> stats(npar());

    for (size_t i = 0; i < _predictive_prior[t].size(); i++) {
        for (size_t j = 0; j < npar(); j++) {
            int particle_idx = _predictive_prior[t][i];
            double par_value = _particle_parameters[t](particle_idx, j);
            stats[j].Push(par_value);
        }
    }
    for (size_t j = 0; j < npar(); j++) {
        _model_pars[j]->append_doubled_variance( 2 * stats[j].Variance() );
    }
}

void AbcSmc::normalize_weights( vector<double>& weights ) {
    double total = 0;
    for (size_t i = 0; i < weights.size(); i++) {
        total += weights[i];
    }

    for (size_t i = 0; i < weights.size(); i++) {
        weights[i] /= total;
    }
}

void AbcSmc::calculate_predictive_prior_weights(int set_num) {
    // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
    calculate_doubled_variances( set_num );
    if (set_num == 0) {
        // uniform weights for set 0 predictive prior
        //_weights[set_num].resize(_predictive_prior[set_num].size(), 1.0/(double) _predictive_prior[set_num].size());
        const double uniform_wt = 1.0/(double) _predictive_prior[set_num].size();
        _weights.push_back( vector<double>(_predictive_prior[set_num].size(), uniform_wt) );
    } else if ( set_num > 0 ) {
        // weights from set - 1 are needed to calculate weights for current set
        _weights.push_back( vector<double>( _predictive_prior[set_num].size(), 0.0 ) );
        //_weights[set_num].resize( _predictive_prior[set_num].size() );

        for (size_t i = 0; i < _predictive_prior[set_num].size(); i++) {
            double numerator = 1;
            double denominator = 0.0;
            for (size_t j = 0; j < npar(); j++) {
                Parameter* par = _model_pars[j];
                const double par_value = _particle_parameters[set_num](_predictive_prior[set_num][i], j);
                if (par->get_prior_type() == NORMAL) {
                    numerator *= gsl_ran_gaussian_pdf(par_value - par->get_prior_mean(), par->get_prior_stdev());
                } else if (par->get_prior_type() == UNIFORM) {
                    // The RHS here will be 1 under normal circumstances.  If the prior has been revised during a fit,
                    // this should throw out values outside of the prior's range
                    numerator *= (int) (par_value >= par->get_prior_min() and par_value <= par->get_prior_max());
                }
            }

            for (size_t k = 0; k < _predictive_prior[set_num - 1].size(); k++) {
                double running_product = _weights[set_num - 1][k];
                for (size_t j = 0; j < npar(); j++) {
                    double par_value = _particle_parameters[set_num](_predictive_prior[set_num][i], j);
                    double old_par_value = _particle_parameters[set_num-1](_predictive_prior[set_num-1][k], j);
                    double old_doubled_variance = _model_pars[j]->get_doubled_variance(set_num-1);

                    // This conditional handles the (often improbable) case where a parameter has completely converged.
                    // It allows ABC to continue exploring other parameters, rather than causing the math
                    // to fall apart because the density at the converged value is infinite.
                    if (old_doubled_variance != 0 or par_value != old_par_value) {
                        running_product *= gsl_ran_gaussian_pdf(par_value-old_par_value, sqrt(old_doubled_variance) );
                    }
                }
                denominator += running_product;
            }
            _weights[set_num][i] = numerator / denominator;
        }
        normalize_weights( _weights[set_num] );
    }
}


gsl_matrix* AbcSmc::setup_mvn_sampler(const int set_num) {
    // ALLOCATE DATA STRUCTURES
    // variance-covariance matrix calculated from pred prior values
    // NB: always allocate this small matrix, so that we don't have to check whether it's safe to free later
    const size_t num_pars = static_cast<size_t>(_particle_parameters[set_num-1].cols());
    const size_t pred_prior_size = _predictive_prior[set_num-1].size();
    gsl_matrix* sigma_hat = gsl_matrix_alloc(num_pars, num_pars);

    if (use_mvn_noise) {
        // container for predictive prior aka posterior from last set
        gsl_matrix* posterior_par_vals = gsl_matrix_alloc(pred_prior_size, num_pars);

        // INITIALIZE DATA STRUCTURES
        for (size_t i = 0; i < pred_prior_size; ++i) {
            for (size_t j = 0; j < num_pars; ++j) {
                // copy values from pred prior into a gsl matrix
                const size_t particle_idx = _predictive_prior[set_num-1][i];
                gsl_matrix_set(posterior_par_vals, i, j, _particle_parameters[set_num-1](particle_idx, j));
            }
        }
        // calculate maximum likelihood estimate of variance-covariance matrix sigma_hat
        gsl_ran_multivariate_gaussian_vcov(posterior_par_vals, sigma_hat);

        for (size_t j = 0; j < npar(); j++) {
            // sampling is done using a kernel with a broader kernel than found in pred prior values
            const double doubled_variance = 2 * gsl_matrix_get(sigma_hat, j, j);
            gsl_matrix_set(sigma_hat, j, j, doubled_variance);
        }

        // not a nice interface, gsl.  sigma_hat is converted in place from a variance-covariance matrix
        // to the same values in the upper triangle, and the diagonal and lower triangle equal to L,
        // the Cholesky decomposition
        gsl_linalg_cholesky_decomp1(sigma_hat);

        gsl_matrix_free(posterior_par_vals);
    }

    return sigma_hat;
}


Row AbcSmc::sample_mvn_predictive_priors( int set_num, const gsl_rng* RNG, gsl_matrix* L ) {
    // container for sampled values
    Row par_values = Row::Zero(npar());
    // SELECT PARTICLE FROM PRED PRIOR TO USE AS EXPECTED VALUE OF NEW SAMPLE
    int r = gsl_rng_nonuniform_int(_weights[set_num-1], RNG);
    gsl_vector* par_val_hat = gsl_vector_alloc(npar());
    for (size_t j = 0; j < npar(); j++) {
        const int particle_idx = _predictive_prior[set_num-1][r];
        const double par_value = _particle_parameters[set_num-1](particle_idx, j);
        gsl_vector_set(par_val_hat, j, par_value);
    }
    par_values = rand_trunc_mv_normal( _model_pars, par_val_hat, L, RNG );
    gsl_vector_free(par_val_hat);

    return par_values;
}

Row AbcSmc::sample_predictive_priors( int set_num, const gsl_rng* RNG ) {
    Row par_values = Row::Zero(npar());
    // Select a particle index r to use from the predictive prior
    int r = gsl_rng_nonuniform_int(_weights[set_num-1], RNG);
    for (size_t j = 0; j < npar(); j++) {
        int particle_idx = _predictive_prior[set_num-1][r];
        double par_value = _particle_parameters[set_num-1](particle_idx, j);
        const Parameter* parameter = _model_pars[j];
        double doubled_variance = parameter->get_doubled_variance(set_num-1);
        double par_min = parameter->get_prior_min();
        double par_max = parameter->get_prior_max();
        par_values(j) = rand_trunc_normal( par_value, doubled_variance, par_min, par_max, RNG );

        if (parameter->get_numeric_type() == INT) {
            par_values(j) = (double) ((int) (par_values(j) + 0.5));
        }
    }
    return par_values;
}
