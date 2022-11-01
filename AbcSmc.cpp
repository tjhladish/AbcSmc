#include "AbcSmc.h"
#include "pls.h"
#include "RunningStat.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>

// need a positive int that is very unlikely
// to be less than the number of particles
#define STOP_TAG 10000000

using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;
using std::setw;

using namespace sqdb;
using namespace ABC;
using namespace chrono;

const string double_bar = "=========================================================================================";
const int PREC = 5;
const int WIDTH = 12;
const string JOB_TABLE  = "job";
const string MET_TABLE  = "met";
const string PAR_TABLE  = "par";
const string UPAR_TABLE = "upar";

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}


vector<float> as_float_vector(Json::Value val, string key) { // not worth templating, despite appearances
    vector<float> extracted_vals;
    if ( val[key].isDouble() ) {
        extracted_vals.push_back( val[key].asFloat() );
    } else if ( val[key].isArray() ) {
        for ( unsigned int i = 0; i < val[key].size(); ++i) extracted_vals.push_back( val[key][i].asFloat() );
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
        for ( unsigned int i = 0; i < val[key].size(); ++i) extracted_vals.push_back( val[key][i].asInt() );
    } else {
        cerr << "Unfamiliar value type associated with " << key << " in configuration file: expecting ints or array of ints." << endl;
        exit(-216);
    }
    return extracted_vals;
}


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
        for (unsigned int i = 0; i < ppfs.size(); ++i) {
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
        const int max_set = max(_predictive_prior_sizes.size(), _smc_set_sizes.size());
        for (int i = 0; i < max_set; ++i) {
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

    string executable = par.get("executable", "").asString();
    if (executable != "") { set_executable( executable ); }

    string resume_dir = par.get("resume_directory", "").asString();
    if (resume_dir != "") {
        if (_mp->mpi_rank == mpi_root) cerr << "Resuming in directory: " << resume_dir << endl;
        set_resume_directory( resume_dir );
        set_resume( true );
    }

    bool computedSamples = !par.isMember("num_samples");
    int nsamples = 1;

    if (!computedSamples) {
      set_smc_iterations( par["smc_iterations"].asInt() ); // TODO: or have it test for convergence
      _smc_set_sizes = as_int_vector(par, "num_samples");
      process_predictive_prior_arguments(par);
      set_pls_validation_training_fraction( par["pls_training_fraction"].asFloat() ); // fraction of runs to use for training
    } else {
      cerr << "Since 'num_samples' unset, assuming forecasting mode" << endl
           << "Setting smc_iterations = 1, predictive_prior_fraction = 1.0, pls_training_fraction = 1.0" << endl;
      set_smc_iterations(1);
      // set_predictive_prior_fraction( 1.0 );
      set_pls_validation_training_fraction( 1.0 ); // fraction of runs to use for training
      cerr << "Computing num_samples from product (all PSEUDO parameters) * POSTERIOR size..." << endl;
    }

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
    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Build name lookup
        string name = model_par[i]["name"].asString();
        assert(par_name_idx.count(name) == 0);
        par_name_idx[name] = i;
    }

    int posterior_size = 1;

    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Iterates over the sequence elements.
        string name = model_par[i]["name"].asString();
        string short_name = model_par[i].get("short_name", "").asString();

        PriorType ptype = UNIFORM;
        string ptype_str = model_par[i]["dist_type"].asString();
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
        string ntype_str = model_par[i]["num_type"].asString();
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
        if (not model_par[i].isMember("untransform")) {
            _untransform_func = [](const double t) { return t; };
        } else if (model_par[i]["untransform"].type() == Json::ValueType::stringValue) {
            string ttype_str = model_par[i].get("untransform", "NONE").asString();
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

        } else if (model_par[i]["untransform"].type() == Json::ValueType::objectValue) {
            Json::Value untransform = model_par[i]["untransform"];
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

        double par1 = model_par[i]["par1"].asDouble();
        double par2 = model_par[i]["par2"].asDouble();
        double step = model_par[i].get("step", 1.0).asDouble(); // default increment is 1

        if (ptype == POSTERIOR) {
          posterior_size = static_cast<int>(floor((par2-par1)/step)+1); // overwrites; TODO assert overwriting with same size
        } else if (ptype == PSEUDO and computedSamples) { // only updating pseudo samples
          int span = (par2-par1 < step) ? 1 : static_cast<int>(round((par2-par1)/step)+1);
          cerr << name << " : " << span <<  "(" << par1 << " , " << par2 << " by " << step << ")" << endl;
          nsamples *= span;
        }

        add_next_parameter(name, short_name, ptype, ntype, par1, par2, step, _untransform_func, par_rescale, mod_map);
    }

    if (computedSamples) {
      nsamples *= posterior_size;
      cerr << "The computed num_samples is " << nsamples << endl;
      _smc_set_sizes = { nsamples };
      _predictive_prior_sizes = { nsamples };
    }

    // Parse model metrics
    const Json::Value model_met = par["metrics"];
    for ( unsigned int i = 0; i < model_met.size(); ++i )  {// Iterates over the sequence elements.
        string name = model_met[i]["name"].asString();
        string short_name = model_met[i].get("short_name", "").asString();

        NumericType ntype = INT;
        string ntype_str = model_met[i]["num_type"].asString();
        if (ntype_str == "INT") {
            ntype = INT;
        } else if (ntype_str == "FLOAT") {
            ntype = FLOAT;
        } else {
            cerr << "Unknown metric numeric type: " << ntype_str << ".  Aborting." << endl;
            exit(-209);
        }

        double val = model_met[i]["value"].asDouble();

        add_next_metric(name, short_name, ntype, val);
    }

    return true;
}


vector<double> AbcSmc::do_complicated_untransformations(vector<Parameter*>& _model_pars, Row& pars) {
    assert( (signed) _model_pars.size() == npar() );
    assert( (signed) pars.size() == npar() );
    const vector<double> identities = {0.0, 1.0, 0.0, 1.0};
    vector<double> upars(npar());
    for (int i = 0; i < npar(); ++i) {
//cerr << "Parameter " << i << ": " << _model_pars[i]->get_name() << endl;
        const Parameter* mpar = _model_pars[i];
        vector<double> modifiers(identities); // TODO -- double check that this is a legit copy constructor
        map<string, vector<int> > mod_map = mpar->get_par_modification_map();
        for (unsigned int j = 0; j < mod_map["transformed_addend"].size(); ++j)   modifiers[0] += pars[mod_map["transformed_addend"][j]];
        for (unsigned int j = 0; j < mod_map["transformed_factor"].size(); ++j)   modifiers[1] *= pars[mod_map["transformed_factor"][j]];
        for (unsigned int j = 0; j < mod_map["untransformed_addend"].size(); ++j) modifiers[2] += pars[mod_map["untransformed_addend"][j]];
        for (unsigned int j = 0; j < mod_map["untransformed_factor"].size(); ++j) modifiers[3] *= pars[mod_map["untransformed_factor"][j]];
        upars[i] = mpar->untransform(pars[i], modifiers);
    }
    return upars;
}

bool AbcSmc::process_database(const unsigned long int rngseed) {
    const gsl_rng* GSL_RNG = gsl_rng_alloc (gsl_rng_taus2);
    gsl_rng_set(GSL_RNG, rngseed);
    return process_database(GSL_RNG);
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

    cerr << std::setprecision(PREC);

    vector< vector<int> > serials;
    if (not read_SMC_sets_from_database(db, serials)) return false; // slurp sets & do particle filtering (identify pred prior) if needed
    const int next_set = serials.size();
    const int last_set = next_set - 1; // this set number

//    if ( _pred_prior_size = 0) set_predictive_prior_size(last_set);

    report_convergence_data(last_set);
    cerr << endl << endl;

    if (_num_smc_sets > next_set) {

        db.Query("BEGIN EXCLUSIVE;").Next();
        //ss << "insert into sets values ( 0, 'Q'"; for (int j = 0; j < npar(); j++) ss << ", NULL"; ss << ");";
        //_db_execute_stringstream(db, ss);

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

            ss << "insert into " << JOB_TABLE << " values ( " << serial << ", "
                                               << next_set << ", "
                                               << i << ", "
                                               << time(NULL)
                                               << ", NULL, 'Q', -1, 0 );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);

            const unsigned long int seed = gsl_rng_get(RNG); // seed for particle
            ss << "insert into " << PAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (int j = 0; j<npar(); j++)  ss << ", " << pars[j]; ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);

            if (use_transformed_pars) {
                vector<double> upars = do_complicated_untransformations(_model_pars, pars);
                ss << "insert into " << UPAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (int j = 0; j<npar(); j++)  ss << ", " << upars[j]; ss << " );";
                //cerr << "attempting: " << ss.str() << endl;
                _db_execute_stringstream(db, ss);
            }

            ss << "insert into " << MET_TABLE << " values ( " << serial; for (int j = 0; j<nmet(); j++) ss << ", NULL"; ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);
        }
        db.CommitTransaction();
        gsl_matrix_free(L);
    } else {
        cerr << "Database already contains " << _num_smc_sets << " complete sets.\n";
    }

    return true;
}


void print_stats(ostream& stream, string str1, string str2, double val1, double val2, double delta, double pct_chg, string tail) {
    stream << "    " + str1 + ", " + str2 + "  ( delta, % ): "  << setw(WIDTH) << val1 << ", " << setw(WIDTH) << val2
                                                      << " ( " << setw(WIDTH) << delta << ", " << setw(WIDTH) << pct_chg  << "% )\n" + tail;
}


void AbcSmc::report_convergence_data(int t) {
    if((signed) _predictive_prior.size() <= t) {
        cerr << "ERROR: attempting to report stats for set " << t << ", but data aren't available.\n"
             << "       This can happen if --process is called on a database that is not ready to be processed.\n";
        exit(-214);
    }
    vector<double> last_means( npar(), 0 );
    vector<double> current_means( npar(), 0 );
    for (int j = 0; j < npar(); j++) {
    //cerr << "par " << j << endl;
        int N = _predictive_prior[t].size();
        for (int i = 0; i < N; i++) {
            double particle_idx = _predictive_prior[t][i];
            double par_value = _particle_parameters[t](particle_idx, j);
            current_means[j] += par_value;
        }
        current_means[j] /= N;

        if (t > 0) {
            int N2 = _predictive_prior[t-1].size();
            for (int i = 0; i < N2; i++) {
                double particle_idx = _predictive_prior[t-1][i];
                double par_value = _particle_parameters[t-1](particle_idx, j);
                last_means[j] += par_value;
            }
            last_means[j] /= N2;
        }
    }

    cerr << double_bar << endl;
    if (t == 0) {
        cerr << "Predictive prior summary statistics:\n";
    } else {
        cerr << "Convergence data for predictive priors:\n";
    }
    for (unsigned int i = 0; i < _model_pars.size(); i++) {
        const Parameter* par = _model_pars[i];
        const double current_stdev = sqrt(par->get_doubled_variance(t)/2.0);
        const double prior_mean = par->get_prior_mean();
        const double prior_mean_delta = current_means[i] - prior_mean;
        const double prior_mean_pct_chg = prior_mean != 0 ? 100 * prior_mean_delta / prior_mean : INFINITY;

        const double prior_stdev = par->get_prior_stdev();
        const double prior_stdev_delta = current_stdev - prior_stdev;
        const double prior_stdev_pct_chg = prior_stdev != 0 ? 100 * prior_stdev_delta / prior_stdev : INFINITY;
        if (t == 0) {
            cerr << "  Par " << i << ": \"" << par->get_name() << "\"\n";

            cerr << "  Means:\n";
            print_stats(cerr, "Prior", "current", prior_mean, current_means[i], prior_mean_delta, prior_mean_pct_chg, "");
            cerr << "  Standard deviations:\n";
            print_stats(cerr, "Prior", "current", prior_stdev, current_stdev, prior_stdev_delta, prior_stdev_pct_chg, "\n");
        } else {
            double last_stdev = sqrt(_model_pars[i]->get_doubled_variance(t-1)/2.0);
            double delta, pct_chg;

            cerr << "  Par " << i << ": \"" << _model_pars[i]->get_name() << "\"\n";

            delta = current_means[i] - last_means[i];
            pct_chg = last_means[i] != 0 ? 100 * delta / last_means[i] : INFINITY;
            cerr << "  Means:\n";
            print_stats(cerr, "Prior", "current", prior_mean, current_means[i], prior_mean_delta, prior_mean_pct_chg, "");
            print_stats(cerr, "Last", " current", last_means[i], current_means[i], delta, pct_chg, "\n");

            delta = current_stdev - last_stdev;
            pct_chg = last_stdev != 0 ? 100 * delta / last_stdev : INFINITY;
            cerr << "  Standard deviations:\n";
            print_stats(cerr, "Prior", "current", prior_stdev, current_stdev, prior_stdev_delta, prior_stdev_pct_chg, "");
            print_stats(cerr, "Last", " current", last_stdev, current_stdev, delta, pct_chg, "\n");
        }
    }
}

/*
bool AbcSmc::_update_sets_table(sqdb::Db &db, const int t) {
    stringstream ss;
    ss << "update sets set status = 'D', ";
    for (int j = 0; j < npar(); ++j) ss << ", " << _model_pars[j]->get_short_name() << "_dv = " << _model_pars[j]->get_doubled_variance(t);
    ss << " where smcSet = " << t << ";";
    return _db_execute_stringstream(db, ss);
}
*/

// Read in existing sets and do particle filtering as appropriate
bool AbcSmc::read_SMC_sets_from_database (sqdb::Db &db, vector< vector<int> > &serials) {
    // make sure database looks intact
    if ( not _db_tables_exist(db, {JOB_TABLE, PAR_TABLE, MET_TABLE}) ) {
        cerr << "ERROR: Failed to read SMC set from database because one or more tables are missing.\n";
        return false;
    }

    // make sure set t is a completed set
    QueryStr qstr;
    Statement s = db.Query(("select smcSet, count(*), COUNT(case status when 'D' then 1 else null end) from " + JOB_TABLE + " group by smcSet order by smcSet;").c_str());
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

        // join all three tables for rows with smcSet = t, slurp and store values
        string select_str = "select J.serial, J.particleIdx, J.posterior, " + _build_sql_select_par_string("") + ", " + _build_sql_select_met_string()
                            + "from " + JOB_TABLE + " J, " + MET_TABLE + " M, " + PAR_TABLE + " P where J.serial = M.serial and J.serial = P.serial "
                            + "and J.smcSet = " + to_string((long long) t) + ";";

        serials.push_back( vector<int>(completed_set_size) );

        Statement s2 = db.Query( select_str.c_str() );

        int particle_counter = 0;
        vector<pair<int, int> > posterior_pairs;
        while (s2.Next()) {
            int offset = 3; // first values are J.serial, J.particleIdx, J.rank
            const int serial = s2.GetField(0);
            const int particle_idx = s2.GetField(1);
            const int posterior_rank = s2.GetField(2);

            if (particle_counter != particle_idx) cerr << "ERROR: particle_counter != particle_idx (" << particle_counter << " != " << particle_idx << ")\n";
            assert(particle_counter == particle_idx);
            serials[t][particle_counter] = serial;
            if (posterior_rank > -1) posterior_pairs.push_back(make_pair( posterior_rank, particle_idx ) );
            for(int i=offset; i < offset + npar(); i++) _particle_parameters[t](particle_counter,i-offset) = (double) s2.GetField(i);
            offset += npar();
            for(int i=offset; i < offset + nmet(); i++) _particle_metrics[t](particle_counter,i-offset) = (double) s2.GetField(i);
            particle_counter++;
        }

        //const int posterior_size = posterior_pairs.size() > 0 ? posterior_pairs.size() : _next_predictive_prior_size;
        //_predictive_prior.push_back( vector<int>(posterior_size) );
        if (posterior_pairs.size() > 0) { // This is a set that has already undergone particle filtering & ranking
            _predictive_prior.push_back( vector<int>(posterior_pairs.size()) );
            for (unsigned int i = 0; i < posterior_pairs.size(); i++) {
                const int rank = posterior_pairs[i].first;
                const int idx = posterior_pairs[i].second;
                _predictive_prior.back()[rank] = idx;
            }
        } else { // Otherwise, do the filtering now and update the DB
            //set_next_predictive_prior_size(t, _particle_parameters[t].size());
            const int next_pred_prior_size = get_pred_prior_size(t);
            _predictive_prior.push_back( vector<int>(next_pred_prior_size) );
            cerr << double_bar << endl << "Set " << t << endl << double_bar << endl;
            if (use_pls_filtering) {
                _filter_particles( t, _particle_metrics[t], _particle_parameters[t], next_pred_prior_size );
            } else {
                _filter_particles_simple( t, _particle_metrics[t], _particle_parameters[t], next_pred_prior_size );
            }
            vector<string> update_strings(next_pred_prior_size);
            for (int i = 0; i < next_pred_prior_size; i++) { // best to worst performing particle in posterior?
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


/*void AbcSmc::set_next_predictive_prior_size(int set_idx, int set_size) {
    _next_predictive_prior_size = round(set_size * _predictive_prior_fraction);
}*/


void AbcSmc::_particle_scheduler(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {
#ifdef USING_MPI
    MPI_Status status;
    long double *send_data = new long double[npar() * _num_particles](); // all params, flattened array
    long double *rec_data  = new long double[nmet() * _num_particles](); // all metrics, flattened array

    // sample parameter distributions; copy values into Y matrix and into send_data buffer
    Row par_row;
    for (int i = 0; i<_num_particles; i++) {
        if (t == 0) { // sample priors
            par_row = sample_priors(RNG);
            if (i<_num_particles) Y_orig.row(i) = par_row;
        } else {      // sample predictive priors
            par_row = sample_predictive_priors(t, RNG);
            if (i<_num_particles) Y_orig.row(i) = par_row;
        }
        for (int j = 0; j<npar(); j++) { send_data[i*npar() + j] = par_row(j); }
    }
    // Seed the workers with the first 'num_workers' jobs
    int particle_id = 0;
    // Don't send particles that don't exist!
    for (int rank = 1; rank < _mp->mpi_size and rank <= _num_particles; ++rank) {
        particle_id = rank - 1;                       // which row in Y
        MPI_Send(&send_data[particle_id*npar()],  // message buffer
                 npar(),                          // number of elements
                 MPI_LONG_DOUBLE,                 // data item is a double
                 rank,                            // destination process rank
                 particle_id,                     // message tag
                 _mp->comm);                      // always use this
    }

    // Receive a result from any worker and dispatch a new work request
    long double rec_buffer[nmet()];
    particle_id++; // move cursor to next particle to be sent
    while ( particle_id < _num_particles ) {
        MPI_Recv(&rec_buffer,                     // message buffer
                 nmet(),                          // message size
                 MPI_LONG_DOUBLE,                 // of type double
                 MPI_ANY_SOURCE,                  // receive from any sender
                 MPI_ANY_TAG,                     // any type of message
                 _mp->comm,                       // always use this
                 &status);                        // received message info
        for (int m = 0; m<nmet(); ++m) rec_data[status.MPI_TAG*nmet() + m] = rec_buffer[m];

        MPI_Send(&send_data[particle_id*npar()],  // message buffer
                 npar(),                          // number of elements
                 MPI_LONG_DOUBLE,                 // data item is a double
                 status.MPI_SOURCE,               // send it to the rank that just finished
                 particle_id,                     // message tag
                 _mp->comm);                      // always use this
        particle_id++; // move cursor
    }

    // receive results for outstanding work requests--there are exactly 'num_workers'
    // or '_num_particles' left, whichever is smaller
    for (int rank = 1; rank < _mp->mpi_size and rank <= _num_particles; ++rank) {
        MPI_Recv(&rec_buffer,
                 nmet(),
                 MPI_LONG_DOUBLE,
                 MPI_ANY_SOURCE,
                 MPI_ANY_TAG,
                 _mp->comm,
                 &status);
        for (int m = 0; m<nmet(); ++m) rec_data[status.MPI_TAG*nmet() + m] = rec_buffer[m];
    }

    // Tell all the workers they're done for now
    for (int rank = 1; rank < _mp->mpi_size; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, STOP_TAG, _mp->comm);
    }

    vector<int> bad_particle_idx; // bandaid, in case simulator returns nonsense values
    for (int i = 0; i<_num_particles; i++) {
        for (int j = 0; j<nmet(); j++) {
            const double met_val = rec_data[i*nmet() + j];
            if (not isfinite(met_val)) bad_particle_idx.push_back(i);
            X_orig(i,j) = rec_data[i*nmet() + j];
        }
    }

    // bandaid, in case simulator returns nonsense values
    for (unsigned int i = 0; i < bad_particle_idx.size(); i++) {
        X_orig.row(i).fill(numeric_limits<float_type>::min());
        Y_orig.row(i).fill(numeric_limits<float_type>::min());
    }

    delete[] send_data;
    delete[] rec_data;
#endif
}


void AbcSmc::_particle_worker() {
#ifdef USING_MPI
    MPI_Status status;
    long double *local_Y_data = new long double[npar()]();
    long double *local_X_data = new long double[nmet()]();

    while (1) {
        MPI_Recv(local_Y_data,
                 npar(),
                 MPI_LONG_DOUBLE,
                 mpi_root,
                 MPI_ANY_TAG,
                 _mp->comm,
                 &status);

        // Check the tag of the received message.
        if (status.MPI_TAG == STOP_TAG) {
            delete[] local_Y_data;
            delete[] local_X_data;
            return;
        }

        Row pars(npar());
        Row mets(nmet());

        for (int j = 0; j<npar(); j++) { pars[j] = local_Y_data[j]; }

        bool success = _run_simulator(pars, mets);
        if (not success) exit(-210);

        for (int j = 0; j<nmet(); j++) { local_X_data[j] = mets[j]; }

        MPI_Send(local_X_data,
                 nmet(),
                 MPI_LONG_DOUBLE,
                 mpi_root,
                 status.MPI_TAG,
                 _mp->comm);
    }
#endif
}


bool AbcSmc::_run_simulator(Row &par, Row &met, const unsigned long int rng_seed, const unsigned long int serial) {
//cerr << par << endl;
    bool particle_success = true;
    if (use_simulator) {
        vector<float_type> met_vec = _simulator( as_vector(par), rng_seed, serial, _mp );
        if ((signed) met_vec.size() != nmet()) {
            cerr << "ERROR: simulator function returned the wrong number of metrics: expected " << nmet() << ", received " << met_vec.size() << endl;
            particle_success = false;
        }
        met = as_row(met_vec);
    } else if (use_executable) {
        string command = _executable_filename;
        for (int j = 0; j<npar(); j++) { command += " " + toString(par[j]); }

        string retval = exec(command);
        if (retval == "ERROR" or retval == "") {
            cerr << command << " does not exist or appears to be an invalid simulator on MPI rank " << _mp->mpi_rank << endl;
            particle_success = false;
        }
        stringstream ss;
        ss.str(retval);

        for (int j = 0; j<nmet(); j++) {
            if (particle_success) {
                ss >> met[j];
            } else {
                met[j] = numeric_limits<float_type>::min();
            }
        }
    } else {
        cerr << "ERROR: A pointer to a simulator function (prefered) or an external simulator executable must be defined.\n";
        particle_success = false;
    }
    return particle_success;
}


// TODO - these could likely be refactored to a single, more clever build sql string function

string AbcSmc::_build_sql_create_par_string( string tag = "" ) {
    stringstream ss;
    for (int i = 0; i<npar()-1; i++) { ss << _model_pars[i]->get_short_name() << tag << " real, "; }
    ss << _model_pars.back()->get_short_name() << tag << " real ";
    return ss.str();
}


string AbcSmc::_build_sql_create_met_string( string tag = "" ) {
    stringstream ss;
    for (int i = 0; i<nmet()-1; i++) { ss << _model_mets[i]->get_short_name() << tag << " real, "; }
    ss << _model_mets.back()->get_short_name() << tag << " real ";
    return ss.str();
}


string AbcSmc::_build_sql_select_par_string( string tag = "" ) {
    stringstream ss;
    for (int i = 0; i<npar()-1; i++) { ss << "P." << _model_pars[i]->get_short_name() << tag << ", "; }
    ss << "P." << _model_pars.back()->get_short_name() << tag << " ";
    return ss.str();
}


string AbcSmc::_build_sql_select_met_string() {
    stringstream ss;
    for (int i = 0; i<nmet()-1; i++) { ss << "M." << _model_mets[i]->get_short_name() << ", "; }
    ss << "M." << _model_mets.back()->get_short_name() << " ";
    return ss.str();
}


bool AbcSmc::_db_execute_strings(sqdb::Db &db, vector<string> &update_buffer) {
    bool db_success = false;
    try {
        db.Query("BEGIN EXCLUSIVE;").Next();
        for (unsigned int i = 0; i < update_buffer.size(); ++i) {
            db.Query(update_buffer[i].c_str()).Next();
        }
        db_success = true;
        db.CommitTransaction();
    } catch (const Exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        for (unsigned int i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    } catch (const exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        for (unsigned int i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    }
    return db_success;
}


bool AbcSmc::_db_execute_stringstream(sqdb::Db &db, stringstream &ss) {
    // We don't need BEGIN EXCLUSIVE here because the calling function has already done it
    bool db_success = false;
    try {
        db_success = db.Query(ss.str().c_str()).Next();
    } catch (const Exception& e) {
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        cerr << ss.str() << endl;
    } catch (const exception& e) {
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        cerr << ss.str() << endl;
    }
    ss.str(string());
    ss.clear();
    return db_success;
}


bool AbcSmc::_db_tables_exist(sqdb::Db &db, vector<string> table_names) {
    // Note that retval here is whether tables exist, rather than whether
    // db transaction was successful.  A failed transaction will throw
    // an exception and exit.
    bool tables_exist = true;
    try {
        for(string table_name: table_names) {
            string query_str = "select count(*) from sqlite_master where type='table' and name='" + table_name + "';";
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


bool AbcSmc::build_database(const gsl_rng* RNG) {
    if (_posterior_database_filename != "" and not file_exists(_posterior_database_filename.c_str())) {
        cerr << "ERROR: Cannot read in posterior database: " << _posterior_database_filename << " does not exist\n";
        exit(-215);
    }

    sqdb::Db db(_database_filename.c_str());

    stringstream ss;
    if ( !_db_tables_exist(db, {JOB_TABLE}) and !_db_tables_exist(db, {PAR_TABLE}) and !_db_tables_exist(db, {MET_TABLE}) ) {
        db.Query("BEGIN EXCLUSIVE;").Next();
        ss << "create table " << JOB_TABLE << " ( serial int primary key asc, smcSet int, particleIdx int, startTime int, duration real, status text, posterior int, attempts int );";
        _db_execute_stringstream(db, ss);

        ss << "create index idx1 on " << JOB_TABLE << " (status, attempts);";
        _db_execute_stringstream(db, ss);

        ss << "create table " << PAR_TABLE << " ( serial int primary key, seed blob, " << _build_sql_create_par_string("") << ");";
        _db_execute_stringstream(db, ss);

        if (use_transformed_pars) {
            ss << "create table " << UPAR_TABLE << " ( serial int primary key, seed blob, " << _build_sql_create_par_string("") << ");";
            _db_execute_stringstream(db, ss);
        }

        ss << "create table " << MET_TABLE << " ( serial int primary key, " << _build_sql_create_met_string("") << ");";
        _db_execute_stringstream(db, ss);
        db.CommitTransaction();
    } else {
        return false;
    }

    Mat2D posterior;
    if (_posterior_database_filename != "") {
        posterior = slurp_posterior();
    }

    db.Query("BEGIN EXCLUSIVE;").Next();

    Row pars;
    const int set_num = 0;
    const int num_particles = get_num_particles(set_num);
    for (int i = 0; i < num_particles; i++) {
        int posterior_rank = -1;
        pars = sample_priors(RNG, posterior, posterior_rank);
        if (not _retain_posterior_rank) posterior_rank = -1;
        QueryStr qstr;

        db.Query(qstr.Format(SQDB_MAKE_TEXT("insert into %s values ( %d, %d, %d, %d, NULL, 'Q', %d, 0 );"), JOB_TABLE.c_str(), i, set_num, i, time(NULL), posterior_rank)).Next();

        Statement s = db.Query(("select last_insert_rowid() from " + JOB_TABLE + ";").c_str());
        s.Next();
        const int rowid = ((int) s.GetField(0)) - 1; // indexing should start at 0

        const unsigned long int seed = gsl_rng_get(RNG); // seed for particle
        ss << "insert into " << PAR_TABLE << " values ( " << rowid << ", '" << seed << "'"; for (int j = 0; j<npar(); j++)  ss << ", " << pars[j]; ss << " );";
        _db_execute_stringstream(db, ss);

        if (use_transformed_pars) {
            vector<double> upars = do_complicated_untransformations(_model_pars, pars);
            ss << "insert into " << UPAR_TABLE << " values ( " << rowid << ", '" << seed << "'"; for (int j = 0; j<npar(); j++)  ss << ", " << upars[j]; ss << " );";
            _db_execute_stringstream(db, ss);
        }

        ss << "insert into " << MET_TABLE << " values ( " << rowid; for (int j = 0; j<nmet(); j++) ss << ", NULL"; ss << " );";
        _db_execute_stringstream(db, ss);
    }
    db.CommitTransaction();
    return true;
}

// alternative approach to verbose flag:
// pass in a "logstream" object
// which is cerr (or whatever log) if being verbose, and dev/null if not
bool AbcSmc::fetch_particle_parameters(sqdb::Db &db, stringstream &select_pars_ss, stringstream &update_jobs_ss, vector<int> &serials, vector<Row> &par_mat, vector<unsigned long int> &rng_seeds, bool verbose) {
    bool db_success = false;
    try {
if (verbose) cerr << "Attempting: " << select_pars_ss.str() << endl;
        db.Query("BEGIN EXCLUSIVE;").Next();
if (verbose) cerr << "Lock obtained" << endl;
        stringstream gotserials;
        Statement s = db.Query(select_pars_ss.str().c_str());
        vector<string> job_strs;

        while (s.Next()) { // grap par vals and find out which jobs need to be updated (changed to running (R), incremented attempts, updated timestamp)
            Row pars(npar());
            const int serial = (int) s.GetField(0);
            serials.push_back( serial );
            const unsigned long int seed = (unsigned long int) s.GetField(1);
            rng_seeds.push_back(seed);
            const int field_offset = 2;
            for (int i = 0; i<npar(); i++) pars[i] = (double) s.GetField(i+field_offset);
            par_mat.push_back(pars);
            gotserials << to_string((long long) serial) << ",";
        }
        auto tmp = gotserials.str();
        tmp.pop_back();
        update_jobs_ss << tmp << ");";
if (verbose) cerr << "Attempting: " << update_jobs_ss.str() << endl;
        db.Query(update_jobs_ss.str().c_str()).Next(); // update jobs table
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


bool AbcSmc::update_particle_metrics(
    sqdb::Db &db, string results
) {
    bool db_success = false;
    stringstream setup, loader, updatejob, updatemet;
    setup << "CREATE TEMPORARY TABLE metupdates (serial INT, ";
    loader << "INSERT INTO metupdates (serial, ";
    updatemet << "UPDATE met SET ";
    for (int j = 0; j<nmet(); j++) {
        setup << _model_mets[j]->get_short_name() << " " << (_model_mets[j]->get_numeric_type() == INT ? "INT" : "FLOAT") << ", ";
        loader << _model_mets[j]->get_short_name() << ", ";
        updatemet << _model_mets[j]->get_short_name() << " = tu." << _model_mets[j]->get_short_name() <<
        (j == (nmet()-1) ? "" : ", ");
    }
    setup << "startTime INT, duration INT);";
    loader << "startTime, duration) VALUES " << results << ";";
    updatemet << " FROM (SELECT tu.* FROM metupdates tu JOIN job USING(serial) WHERE status != 'D') AS tu WHERE tu.serial == met.serial;";

    updatejob << "UPDATE job SET " <<
              "startTime = tu.startTime, duration = tu.duration, status = 'D'" <<
              " FROM metupdates tu WHERE job.serial == tu.serial AND status != 'D';";

// cerr << setup.str() << endl << loader.str() << endl << updatemet.str() << endl << updatejob.str() << endl;

    db.Query("BEGIN EXCLUSIVE;").Next();

    try {
        db.Query(setup.str().c_str()).Next();
        db.Query(loader.str().c_str()).Next();
        db.Query(updatemet.str().c_str()).Next();
        db.Query(updatejob.str().c_str()).Next();
        db.CommitTransaction();
        db_success = true;
    } catch (const Exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while updating metrics:" << endl;
        cerr << results << endl;
        // for (unsigned int i = 0; i < update_metrics_strings.size(); ++i) {
        //     cerr << update_metrics_strings[i] << endl;
        //     cerr << update_jobs_strings[i] << endl;
        // }
    } catch (const exception& e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed while updating metrics:" << endl;
        cerr << results << endl;
    }

    return db_success;
}

bool AbcSmc::simulate_particle_by_serial(const int serial_req) { return simulate_next_particles(1, serial_req, -1); }

bool AbcSmc::simulate_particle_by_posterior_idx(const int posterior_req) { return simulate_next_particles(1, -1, posterior_req); }

bool AbcSmc::simulate_next_particles(
    const int n, const int serial_req, const int posterior_req, const bool verbose
) {
    // TODO refactor this n vs serial vs posterior code smell
    assert(n == 1 or (serial_req == -1 and posterior_req == -1));
    assert(serial_req == -1 or posterior_req == -1);
//bool AbcSmc::simulate_database(const int smc_set, const int particle_id) {
    sqdb::Db db(_database_filename.c_str());
    string model_par_table = _db_tables_exist(db, {UPAR_TABLE}) ? UPAR_TABLE : PAR_TABLE;
    vector<Row> par_mat;  //( n, Row(npar()) ); -- expected size, if n rows are available

    stringstream select_ss;
    // Do already running jobs as well, if there are not enough queued jobs
    // This is because we are seeing jobs fail/time out for extrinsic reasons on the stuporcomputer
    string filter = "";
    string limit = "";
    if (n > 0) {
        limit = " LIMIT " + to_string(n);  
    } else if (serial_req > -1) {
        filter = " AND serial = " + to_string(serial_req);
    } else if (posterior_req > -1) {
        filter = " AND smcSet = (SELECT max(smcSet) FROM " + JOB_TABLE + " WHERE posterior > -1) AND posterior = " + to_string(posterior_req);
    }

    select_ss << "select P.* FROM " << model_par_table << " P JOIN (SELECT * FROM " <<
                 JOB_TABLE << " WHERE status IN ('Q','R')" << filter << " ORDER BY status, attempts" << limit <<
                 ") USING(serial);"; 

    //  line below is much faster for very large dbs, but not all particles will get run e.g. if some particles are killed by scheduler
    //  select_ss << "and J.status = 'Q' limit " << n << ";";

    const size_t overall_start_time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
    // build jobs update statement to indicate job is running
    stringstream update_ss;
    update_ss << "update " << JOB_TABLE << " set startTime = " << overall_start_time << ", status = 'R', attempts = attempts + 1 where serial IN ("; // we don't know the serial yet

    vector<int> serials;
    vector<unsigned long int> rng_seeds;
    bool ok_to_continue = fetch_particle_parameters(db, select_ss, update_ss, serials, par_mat, rng_seeds, verbose);

    if (ok_to_continue) {

        stringstream update_metrics;
        vector<Row> met_mat(par_mat.size(), Row(nmet()));  //( n, Row(nmet()) );     
        
if (verbose) cerr << "Starting simulation for " << par_mat.size() << " simulations." << endl;

        for (unsigned int i = 0; i < par_mat.size(); ++i) {
            const high_resolution_clock::time_point start_time = high_resolution_clock::now();
            const int serial = serials[i];
            bool success = _run_simulator(par_mat[i], met_mat[i], rng_seeds[i], serial);
            if (not success) exit(-211);
            update_metrics << "(" << serial << ", ";
            for (int j = 0; j<nmet(); j++) { update_metrics << met_mat[i][j] << ", "; }
            update_metrics << duration_cast<seconds>(start_time.time_since_epoch()).count() << ", " <<
                duration_cast<milliseconds>(high_resolution_clock::now() - start_time).count() << ")";
            update_metrics << (i == (par_mat.size() - 1) ? "" : ", ");
        }
if (verbose) {
    cerr << "Updating with results ... " << endl;
    // for (auto ss : update_metrics_strings) {
    //     cerr << ss << endl;
    // }
    // for (auto ss : update_jobs_strings) {
    //     cerr << ss << endl;
    // }
}
        update_particle_metrics(db, update_metrics.str());
    } else {
        cerr << "Parameter selection from database failed.\n";
    }

    return true;
}


long double AbcSmc::calculate_nrmse(vector<Col> posterior_mets) {
    long double nrmse = 0.0;
    int n = 0;
    for (int i = 0; i<nmet(); i++) {
        long double obs = _model_mets[i]->get_obs_val();
        long double sim = mean(posterior_mets[i]);
        if (obs != sim) { // also means they're not both zero
            long double expected = (fabs(obs)+fabs(sim))/2;
            nrmse += pow((sim-obs)/expected, 2);
        }
        n++;
    }
    assert(n>0);
    return sqrt(nrmse/n);
}

void AbcSmc::_print_particle_table_header() {
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << _model_pars[i]->get_short_name(); } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << _model_mets[i]->get_short_name(); } cerr << endl;
}

void AbcSmc::_fp_helper (const int t, const Mat2D &X_orig, const Mat2D &Y_orig, const int next_pred_prior_size, const Col& distances) {
    vector<int> ranking = ordered(distances);

    vector<int>::iterator first = ranking.begin();
    vector<int>::iterator last  = ranking.begin() + next_pred_prior_size;
    vector<int> sample(first, last); // This is the predictive prior / posterior
    _predictive_prior[t] = sample;

    cerr << "Observed:\n";
    _print_particle_table_header();
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << "---"; } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << _model_mets[i]->get_obs_val(); } cerr << endl;

    vector<Col> posterior_pars(npar(), Col(next_pred_prior_size));
    vector<Col> posterior_mets(nmet(), Col(next_pred_prior_size));
    for (unsigned int i = 0; i < sample.size(); ++i) {
        const int idx = sample[i];
        for (int j = 0; j < npar(); j++) posterior_pars[j](i) = Y_orig(idx, j);
        for (int j = 0; j < nmet(); j++) posterior_mets[j](i) = X_orig(idx, j);
    }
    cerr << "Normalized RMSE for metric means (lower is better):  " << calculate_nrmse(posterior_mets) << "\n";
    cerr << "Posterior means:\n";
    _print_particle_table_header();
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << mean(posterior_pars[i]); } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << mean(posterior_mets[i]); } cerr << endl;

    cerr << "Posterior medians:\n";
    _print_particle_table_header();
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << median(posterior_pars[i]); } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << median(posterior_mets[i]); } cerr << endl;

    cerr << "Best five:\n";
    _print_particle_table_header();
    for (int q=0; q<5; q++) {
        const int idx = ranking[q];
        for (int i = 0; i < Y_orig.cols(); i++) { cerr << setw(WIDTH) << Y_orig(idx, i); } cerr << " | ";
        for (int i = 0; i < X_orig.cols(); i++) { cerr << setw(WIDTH) << X_orig(idx, i); } cerr << endl;
    }

    cerr << "Worst five:\n";
    _print_particle_table_header();
    for (unsigned int q=ranking.size()-1; q>=ranking.size()-5; q--) {
        const int idx = ranking[q];
        for (int i = 0; i < Y_orig.cols(); i++) { cerr << setw(WIDTH) << Y_orig(idx, i); } cerr << " | ";
        for (int i = 0; i < X_orig.cols(); i++) { cerr << setw(WIDTH) << X_orig(idx, i); } cerr << endl;
    }
}

void AbcSmc::_filter_particles_simple (int t, Mat2D &X_orig, Mat2D &Y_orig, int next_pred_prior_size) {
    // x is metrics, y is parameters
    Row X_sim_means, X_sim_stdev;
    Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
    Mat2D Y = colwise_z_scores( Y_orig );
    Row obs_met = _z_transform_observed_metrics( X_sim_means, X_sim_stdev );

    Col   distances  = euclidean(obs_met, X);
    _fp_helper (t, X_orig, Y_orig, next_pred_prior_size, distances);
}

void AbcSmc::_filter_particles (int t, Mat2D &X_orig, Mat2D &Y_orig, int next_pred_prior_size) {
    // Run PLS
    // Box-Cox transform data -- TODO?
    //void test_bc( Mat2D );
    //test_bc(Y_orig);

    Row X_sim_means, X_sim_stdev;
    Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
    Mat2D Y = colwise_z_scores( Y_orig );
    Row obs_met = _z_transform_observed_metrics( X_sim_means, X_sim_stdev );

    int nobs  = X_orig.rows();      // number of observations
    int npred = X_orig.cols();      // number of predictor variables
    int nresp = Y_orig.cols();      // number of response variables
    int ncomp = npar();             // It doesn't make sense to consider more components than model parameters
    PLS_Model plsm;
    plsm.initialize(npred, nresp, ncomp);
    const int pls_training_set_size = round(nobs * _pls_training_fraction);
    plsm.plsr(X.topRows(pls_training_set_size), Y.topRows(pls_training_set_size), KERNEL_TYPE1);

/*
//P, W, R, Q, T
cerr << "P:\n" << plsm.P << endl;
cerr << "W:\n" << plsm.W << endl;
cerr << "R:\n" << plsm.R << endl;
cerr << "Q:\n" << plsm.Q << endl;
cerr << "T:\n" << plsm.T << endl;
cerr << "coefficients:\n" << plsm.coefficients() << endl;
*/
    // A is number of components to use
    for (int A = 1; A<=ncomp; A++) {
        // How well did we do with this many components?
        cerr << setw(2) << A << " components ";
        cerr << "explained variance: " << plsm.explained_variance(X, Y, A);
        //cerr << "root mean squared error of prediction (RMSEP):" << plsm.rmsep(X, Y, A) << endl;
        cerr << " SSE: " << plsm.SSE(X,Y,A) <<  endl;
    }

    const int test_set_size = nobs - pls_training_set_size;
    Rowi num_components = plsm.optimal_num_components(X.bottomRows(test_set_size), Y.bottomRows(test_set_size), NEW_DATA);
    int num_components_used = num_components.maxCoeff();
    cerr << "Optimal number of components for each parameter (validation method == NEW DATA):\t" << num_components << endl;
    cerr << "Using " << num_components_used << " components." << endl;

    // Calculate new, orthogonal metrics (==scores) using the pls model
    // Is casting as real always safe?
    Row   obs_scores = plsm.scores(obs_met, num_components_used).row(0).real();
    Mat2D sim_scores = plsm.scores(X, num_components_used).real();
    Col   distances  = euclidean(obs_scores, sim_scores);
    _fp_helper (t, X_orig, Y_orig, next_pred_prior_size, distances);
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

Mat2D AbcSmc::slurp_posterior() {
    // TODO - handle sql/sqlite errors
    // TODO - if "posterior" database doesn't actually have posterior values, this will fail silently

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


Row AbcSmc::sample_priors(const gsl_rng* RNG, Mat2D& posterior, int &posterior_rank) {
    Row par_sample = Row::Zero(_model_pars.size());
    bool increment_nonrandom_par = true; // only one PSEUDO parameter gets incremented each time
    bool increment_posterior = true;     // posterior parameters get incremented together, when all pseudo pars reach max val
    vector<int> posterior_indices;
    // for each parameter
    for (unsigned int i = 0; i < _model_pars.size(); i++) {
        Parameter* p = _model_pars[i];
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
            posterior_indices.push_back(i);
        } else {
            // Random parameters get sampled independently from each other, and are therefore easy
            val = p->sample(RNG);
        }

        par_sample(i) = val;
    }
    if (_posterior_database_filename != "") {

        for (unsigned int c = 0; c < posterior_indices.size(); ++c) {
            par_sample(posterior_indices[c]) = posterior(posterior_rank, c);
            Parameter* p = _model_pars[posterior_indices[c]];

            if (increment_posterior) {
                if (p->get_state() >= p->get_prior_max()) {
                    p->reset_state();
                } else {
                    p->increment_state();
                }
            }

        }

//cerr << "sample: " << par_sample << endl;
    }
    return par_sample;
}

void AbcSmc::calculate_doubled_variances( int t ) {
    vector<RunningStat> stats(npar());

    for (unsigned int i = 0; i < _predictive_prior[t].size(); i++) {
        for (int j = 0; j < npar(); j++) {
            double particle_idx = _predictive_prior[t][i];
            double par_value = _particle_parameters[t](particle_idx, j);
            stats[j].Push(par_value);
        }
    }
    for (int j = 0; j < npar(); j++) {
        _model_pars[j]->append_doubled_variance( 2 * stats[j].Variance() );
    }
}

void AbcSmc::normalize_weights( vector<double>& weights ) {
    double total = 0;
    for (unsigned int i = 0; i< weights.size(); i++) {
        total += weights[i];
    }

    for (unsigned int i = 0; i< weights.size(); i++) {
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

        for (unsigned int i = 0; i < _predictive_prior[set_num].size(); i++) {
            double numerator = 1;
            double denominator = 0.0;
            for (int j = 0; j < npar(); j++) {
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

            for (unsigned int k = 0; k < _predictive_prior[set_num - 1].size(); k++) {
                double running_product = _weights[set_num - 1][k];
                for (int j = 0; j < npar(); j++) {
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
    const int num_pars = _particle_parameters[set_num-1].cols();
    const int pred_prior_size = _predictive_prior[set_num-1].size();
    gsl_matrix* sigma_hat = gsl_matrix_alloc(num_pars, num_pars);

    if (use_mvn_noise) {
        // container for predictive prior aka posterior from last set
        gsl_matrix* posterior_par_vals = gsl_matrix_alloc(pred_prior_size, num_pars);

        // INITIALIZE DATA STRUCTURES
        for (int i = 0; i < pred_prior_size; ++i) {
            for (int j = 0; j < num_pars; ++j) {
                // copy values from pred prior into a gsl matrix
                const int particle_idx = _predictive_prior[set_num-1][i];
                gsl_matrix_set(posterior_par_vals, i, j, _particle_parameters[set_num-1](particle_idx, j));
            }
        }
        // calculate maximum likelihood estimate of variance-covariance matrix sigma_hat
        gsl_ran_multivariate_gaussian_vcov(posterior_par_vals, sigma_hat);

        for (int j = 0; j<npar(); j++) {
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

Row rand_trunc_mv_normal(const vector<Parameter*> _model_pars, gsl_vector* mu, gsl_matrix* L, const gsl_rng* rng) {
    const int npar = _model_pars.size();
    Row par_values = Row::Zero(npar);
    gsl_vector* result = gsl_vector_alloc(npar);
    bool success = false;
    while (not success) {
        success = true;
        gsl_ran_multivariate_gaussian(rng, mu, L, result);
        for (int j = 0; j < npar; j++) {
            par_values[j] = gsl_vector_get(result, j);
            if (_model_pars[j]->get_numeric_type() == INT) par_values(j) = (double) ((int) (par_values(j) + 0.5));
            if (par_values[j] < _model_pars[j]->get_prior_min() or par_values[j] > _model_pars[j]->get_prior_max()) success = false;
        }
    }
    gsl_vector_free(result);
    return par_values;
}

Row AbcSmc::sample_mvn_predictive_priors( int set_num, const gsl_rng* RNG, gsl_matrix* L ) {
    // container for sampled values
    Row par_values = Row::Zero(npar());
    // SELECT PARTICLE FROM PRED PRIOR TO USE AS EXPECTED VALUE OF NEW SAMPLE
    int r = gsl_rng_nonuniform_int(_weights[set_num-1], RNG);
    gsl_vector* par_val_hat = gsl_vector_alloc(npar());
    for (int j = 0; j<npar(); j++) {
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
    for (int j = 0; j<npar(); j++) {
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

Row AbcSmc::_z_transform_observed_metrics(Row& means, Row& stdevs) {
    Row zmat = Row::Zero(nmet());
    for (int i = 0; i<nmet(); i++) { zmat(i) = (_model_mets[i]->get_obs_val() - means(i)) / stdevs(i); }
    return zmat;
}
