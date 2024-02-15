
#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>
#include <chrono>
#include <algorithm> // all_of
#include <memory> // unique_ptr / shared_ptr

#include <AbcSmc/RunningStat.h>
#include <AbcSmc/AbcSmc.h>
#include <AbcSmc/AbcMPI.h>
#include <AbcSmc/AbcLog.h>
#include <AbcSmc/Priors.h>
#include <AbcSmc/IndexedPars.h>

#include <json/json.h>
#include <Eigen/Dense>

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

const int PREC = 5;
const string JOB_TABLE  = "job";
const string MET_TABLE  = "met";
const string PAR_TABLE  = "par";
const string UPAR_TABLE = "upar";

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

template <typename T>
vector<T> as_vector(const Json::Value &val) {
    vector<T> extracted_vals;
    if (val.isArray()) { for (const Json::Value &jv : val) {
        extracted_vals.push_back( jv.as<T>() ); // NB, jsoncpp handles cast failures
    } } else {
        extracted_vals.push_back( val.as<T>() );
    }
    return extracted_vals;
}

void parse_iterations(
    const Json::Value &par, const size_t pseudosize,
    size_t * iterations,
    float_type * training_frac,
    vector<size_t> * set_sizes,
    vector<size_t> * pred_prior_sizes
) {
    if (pseudosize != 0) { // definitely in projection mode; ignore other parameters
        if (par.get("smc_iterations", 1).asInt() != 1) { // ignoring OR erroring if specified != 1
            cerr << "Cannot use smc_iterations > 1 with ONLY PSEUDO or POSTERIOR parameters.  Aborting." << endl;
            exit(-202);
        }
        if (par.isMember("num_samples")) {
            size_t checksize = (as_vector<size_t>(par["num_samples"]))[0];
            if (checksize != pseudosize) {
                std::cerr << "ERROR: `num_samples` ("<< checksize <<") does not match imputed combinations of PSEUDO and/or POSTERIOR parameters ("<< pseudosize <<")." << endl;
                exit(-201);
            } else {
                std::cerr << "WARNING: specified `num_samples` for all PSEUDO and/or POSTERIOR parameters." << std::endl;
            }
        }
        if (par.isMember("predictive_prior_fraction") or par.isMember("predictive_prior_size")) {
            std::cerr << "WARNING: ignoring `predictive_prior_*` options in projection mode." << std::endl;
        }
        *iterations = 1;
        *set_sizes = { pseudosize };
    } else { // some priors => fitting mode => predictive prior fractions
        // pred_prior_sizes is a series of sizes, but may be passed as EITHER fractions OR sizes.
        // both / neither is an error
        if (
            (par.isMember("predictive_prior_fraction") and par.isMember("predictive_prior_size")) or
            not (par.isMember("predictive_prior_fraction") or par.isMember("predictive_prior_size"))
        ) {
            std::cerr << "Error: exactly one of `predictive_prior_fraction` or `predictive_prior_size` must be specified in configuration file." << std::endl;
            exit(1);
        }

        // if doing iterations, may have a training fraction; defaults to 0.5
        *training_frac = par.get("pls_training_fraction", 0.5).as<float_type>();
        if (*training_frac <= 0 or 1 <= *training_frac) {
            std::cerr << "Error: pls_training_fraction must be in (0, 1)." << std::endl;
            exit(1);
        }

        *set_sizes = as_vector<size_t>(par["num_samples"]);

        if (par.isMember("predictive_prior_fraction")) {
            auto ppfs = as_vector<float_type>(par["predictive_prior_fraction"]);
            if (not std::all_of(ppfs.begin(), ppfs.end(), [](float_type f) { return (0 < f) and (f <= 1); })) {
                cerr << "Error: `predictive_prior_fraction`s must be in (0, 1]" << endl;
                exit(1);
            }
            auto set_sizes_copy = *set_sizes;
            const int max_set = max(ppfs.size(), set_sizes_copy.size());
            ppfs.resize(max_set, ppfs.back());
            set_sizes_copy.resize(max_set, set_sizes_copy.back());
            pred_prior_sizes->resize(max_set);
            for (size_t i = 0; i < (unsigned) max_set; ++i) { (*pred_prior_sizes)[i] = round(ppfs[i] * set_sizes_copy[i]); }
        } else { // if (par.isMember("predictive_prior_size")) {
            *pred_prior_sizes = as_vector<size_t>(par["predictive_prior_size"]);
            const size_t max_set = max(pred_prior_sizes->size(), set_sizes->size());
            size_t i = 0;
            for (; i < max_set; ++i) { if (pred_prior_sizes->at(i) > set_sizes->at(i)) {
                cerr << "Error: requested predictive prior size > SMC set size at: " << i << endl;
                exit(1);
            } }
            // now either done with both, or one has more to do, so if there *are* more...
            // EITHER more pred prior
            for (; i < pred_prior_sizes->size(); ++i) { if (pred_prior_sizes->at(i) > set_sizes->back()) {
                cerr << "Error: requested predictive prior size > SMC set size at: " << i << endl;
                exit(1);
            } }
            // OR more set sizes
            for (; i < set_sizes->size(); ++i) { if (pred_prior_sizes->back() > set_sizes->at(i)) {
                cerr << "Error: requested predictive prior size > SMC set size at: " << i << endl;
                exit(1);
            } }
        }

        *iterations = par.get("smc_iterations", max(set_sizes->size(), pred_prior_sizes->size())).asUInt64();

    }
  
}

MetricPtr parse_metric(const Json::Value &mmet) {
    const string name = mmet["name"].asString();
    const string short_name = mmet.get("short_name", name).asString();
    const float_type val = mmet["value"].asDouble();
    const string ntype_str = mmet["num_type"].asString();

    if (ntype_str == "INT") {
        return std::make_shared<TMetric<int>>(name, short_name, val);
    } else if (ntype_str == "FLOAT") {
        return std::make_shared<TMetric<float_type>>(name, short_name, val);
    } else {
        cerr << "Unknown metric numeric type: " << ntype_str << ".  Aborting." << endl;
        exit(-209);
    }

}

void parse_transform(
    const Json::Value &mparu,
    ParRescale ** pscale,
    ParXform ** pxform,
    transformer ** _untransform_func,
    const std::map<const std::string, const size_t> &par_name_idx
) {
    if (mparu.type() == Json::ValueType::stringValue) {
        string ttype_str = mparu.asString();
        if (ttype_str == "NONE") { // TODO - it's possible this may not actually ever be called
            *_untransform_func = [](const float_type &t) { return t; };
            //ttype = UNTRANSFORMED;
        } else if (ttype_str == "POW_10") {
            *_untransform_func = [](const float_type &t) { return pow(10.0, t); };
            //ttype = LOG_10;
        } else if (ttype_str == "LOGISTIC") {
            *_untransform_func = [](const float_type &t) { return ABC::logistic(t); };
            //ttype = LOGIT;
        } else {
            cerr << "Unknown parameter transformation type: " << ttype_str << ".  Aborting." << endl;
            exit(-206);
        }
        *pscale = new ParRescale();
        *pxform = new ParXform(*_untransform_func);
    } else if (mparu.type() == Json::ValueType::objectValue) {
        string ttype_str = mparu["type"].asString();
        if (ttype_str != "LOGISTIC") {
            cerr << "Only type: LOGISTIC is currently supported for untransformation objects.  (NONE and POW_10 supported as untransformation strings.)\n";
            exit(-207);
        }
        *pscale = new ParRescale(mparu["min"].asDouble(), mparu["max"].asDouble());
        *_untransform_func = [](const float_type &t) { return ABC::logistic(t); };
        //Json::ValueType mod_type = untransform["transformed_addend"].type();
        std::map<string, vector<size_t>> mod_map = {
            { "transformed_addend", {} },
            { "transformed_factor", {} },
            { "untransformed_addend", {} },
            { "untransformed_factor", {} }
        };
        for (auto &mod_type: mod_map) {
            if (mparu.isMember(mod_type.first)) {
                for (const Json::Value &json_val: mparu[mod_type.first]) {
                    const std::string par_name = json_val.asString();
                    mod_type.second.push_back(par_name_idx.at(par_name));
                }
            }
        }
        *pxform = new ParXform(*_untransform_func,
            mod_map["transformed_addend"],   mod_map["transformed_factor"],
            mod_map["untransformed_addend"], mod_map["untransformed_factor"]
        );
    } else {
        cerr << "Unsupported JSON data type associated with 'untransform' parameter key.\n";
        exit(-208);
    }
}

const ParameterPtr parse_parameter(
    const Json::Value &mpar
) {
    const string name = mpar["name"].asString();
    const string short_name = mpar.get("short_name", name).asString();

    const string ptype_str = mpar["dist_type"].asString();
    const string ntype_str = mpar["num_type"].asString();

    ParameterPtr par;

    if (not ((ntype_str == "INT") or (ntype_str == "FLOAT"))) {
        cerr << "Unknown parameter numeric type: " << ntype_str << ".  Aborting." << endl;
        exit(-206);
    }

    if (ptype_str == "UNIFORM") {
        if (ntype_str == "INT") {
            const long par1 = mpar["par1"].asInt64();
            const long par2 = mpar["par2"].asInt64();
            par = std::make_shared<DiscreteUniformPrior>(name, short_name, par1, par2);
        } else {
            const float_type par1 = mpar["par1"].asDouble();
            const float_type par2 = mpar["par2"].asDouble();
            par = std::make_shared<ContinuousUniformPrior>(name, short_name, par1, par2);
        }
    } else if (ptype_str == "NORMAL" or ptype_str == "GAUSSIAN") {
        if (ntype_str == "INT") {
            cerr << "Parameter numeric " << ntype_str << " not supported for parameter type " << ptype_str << ".  Aborting." << endl;
            exit(-206);
        }
        const float_type par1 = mpar["par1"].asDouble();
        const float_type par2 = mpar["par2"].asDouble();
        par = std::make_shared<GaussianPrior>(name, short_name, par1, par2);
    } else if (ptype_str == "PSEUDO") {
        std::vector<float_type> states;
        if (mpar.isMember("vals")) {
            states = as_vector<float_type>(mpar["vals"]);
        } else {
            const float_type smax = mpar["par2"].asDouble();
            const float_type step = mpar.get("step", 1.0).asDouble();
            if (step != 0) {
                const double EPSILON = 0.0001;
                for (float_type s = mpar["par1"].asDouble(); s <= smax + EPSILON * step; s += step) {
                    states.push_back(s);
                }
            } else {
                states.push_back(mpar["par1"].asDouble());
            }
        }
        par = std::make_shared<PseudoPar>(name, short_name, states);
    } else if (ptype_str == "POSTERIOR") {
        // TODO iota + stride?
        const size_t size = mpar["par2"].asUInt64() - mpar["par1"].asUInt64() + 1;
        par = std::make_shared<PosteriorPar>(name, short_name, size);
    } else {
        cerr << "Unknown parameter distribution type: " << ptype_str << ".  Aborting." << endl;
        exit(-205);
    }
    return par;
}

Json::Value prepare(const string &configfilename) {
    if (not file_exists(configfilename.c_str())) {
        cerr << "File does not exist: " << configfilename << endl;
        exit(1);
    }
    // TODO - Make sure any existing database actually reflects what is expected in JSON, particularly that par and met tables are legit
    Json::Value par;   // will contain the par value after parsing.
    Json::Reader reader;
    string json_data = slurp(configfilename);

    if ( !reader.parse( json_data, par ) ) {
        // report to the user the failure and their locations in the document.
        cerr << "Failed to parse configuration\n" << reader.getFormattedErrorMessages();
        exit(1);
    }
    return par;
}

bool _db_tables_exist(sqdb::Db &db, vector<string> table_names) {
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
    } catch (const Exception &e) {
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed while checking whether the following tables exist:";
        for(string table_name: table_names) cerr << " " << table_name;
        cerr << endl;
        exit(-212);
    } catch (const exception &e) {
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed while checking whether the following tables exist:";
        for(string table_name: table_names) cerr << " " << table_name;
        cerr << endl;
        exit(-213);
    }
    return tables_exist;
}

Mat2D slurp_posterior(
    const std::string &_posterior_database_filename,
    const ParameterVec &_model_pars
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
    for (auto p: _model_pars) {
        if (p->isPosterior()) {
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
            posterior(r, c) = static_cast<float_type>(posterior_query.GetField(c));
        }
        ++r;
    }

    return posterior;
}

AbcSmc::AbcSmc() : _met_vals(new Row()) {};

AbcSmc::~AbcSmc() {};

bool AbcSmc::parse_config(const string &conf_filename) {
    auto par = prepare(conf_filename);

    // if we use a posterior from an earlier ABC run to determine some of the parameter values, do we keep the rank?
    set_retain_posterior_rank( par.get("retain_posterior_rank", false).asBool() );

    // Parse model parameters
    const Json::Value model_par = par["parameters"];
    map<const string, const size_t> par_name_idx;
    // TODO find a way to size_t this?
    for (unsigned int parIdx = 0; parIdx < model_par.size(); ++parIdx)  {// Build name lookup
        const string name = model_par[parIdx]["name"].asString();
        assert(par_name_idx.count(name) == 0);
        par_name_idx.emplace(name, parIdx);
    }

    // if we're retaining the posterior rank, we definitely have a posterior; otherwise, determine from parameters
    bool anyPosterior = false;
    // if all parameters are pseudo / posterior, we can compute the size of a run
    size_t pseudosize = 1;
    size_t posterior_size = 0;

    for (const Json::Value &mpar : model_par)  { // Iterates over the sequence elements.
        const ParameterPtr par = parse_parameter(mpar);
        if (par->isPosterior()) {
            if (posterior_size == 0) {
                posterior_size = par->state_size();
                anyPosterior = true;
            } else {
                assert(par->state_size() == posterior_size);
            }
        } else {
            pseudosize *= par->state_size();
        }

        add_next_parameter(par);

        if (mpar.isMember("untransform")) {
            // declare these as pointers, pass them by reference to be created by the function
            transformer * _untransform_func;
            ABC::ParRescale * par_rescale;
            ABC::ParXform * mod_map;
            parse_transform(mpar["untransform"], &par_rescale, &mod_map, &_untransform_func, par_name_idx);
            add_modification_map(par, mod_map);
            add_par_rescale(par, par_rescale);
        }
    }

    if (anyPosterior) {
        pseudosize *= posterior_size;
        if (not par.isMember("posterior_database_filename")) {
            cerr << "Parameter specfied as type POSTERIOR, without previously specifying a posterior_database_filename.  Aborting." << endl;
            exit(-204);
        }
        if (_num_smc_sets > 1) {
            cerr << "Cannot use posterior parameters with multiple SMC sets.  Aborting." << endl;
            exit(-203);
        }
        _posterior = std::make_unique<const Mat2D>(slurp_posterior(par["posterior_database_filename"].asString(), _model_pars));
    } else {
        _posterior = std::make_unique<const Mat2D>();
    }

    for (const Json::Value &mmet : par["metrics"]) {
        auto met = parse_metric(mmet);
        add_next_metric(met);
    }

    parse_iterations(par, pseudosize, &_num_smc_sets, &_pls_training_fraction, &_smc_set_sizes, &_predictive_prior_sizes);

    // TODO these should be mutually exclusive options
    std::string executable = par.get("executable", "").asString();
    if (executable != "") { set_executable( executable ); }
    std::string sharedobj = par.get("shared", "").asString();
    if (sharedobj != "") { set_simulator( sharedobj ); }

    string resume_dir = par.get("resume_directory", "").asString();
    if (resume_dir != "") {
        if (_mp->mpi_rank == mpi_root) cerr << "Resuming in directory: " << resume_dir << endl;
        set_resume_directory( resume_dir );
    }

    // TODO--allow specification of pred prior size (single value or list of values)
    //set_predictive_prior_fraction( par["predictive_prior_fraction"].asFloat() );
    
    set_database_filename( par["database_filename"].asString() );

    const string noise = par.get("noise", "INDEPENDENT").asString();
    if (noise == "INDEPENDENT") {
        _noise = ABC::NOISE::INDEPENDENT;
    } else if (noise == "MULTIVARIATE") {
        _noise = ABC::NOISE::MULTIVARIATE;
    } else {
        cerr << "Unknown parameter noise type specified: " << noise << ". Aborting." << endl;
        exit(-210);
    }

    return true;
}

void AbcSmc::add_next_metric(const MetricPtr m) {
    _model_mets.push_back(std::shared_ptr<const Metric>(m));
    _met_vals->resize(_model_mets.size());
    _met_vals->operator[](_model_mets.size()-1) = m->get_obs_val();
}

void AbcSmc::add_next_parameter(const ParameterPtr p) {
    _model_pars.push_back(p);
}

Row AbcSmc::_to_model_space(
    const Row &fitting_space_pars
) {
    assert( _model_pars.size() == (unsigned) fitting_space_pars.size() );
    Row model_space_pars = fitting_space_pars; // copy initially - all model_space_pars == fitting_space_pars
    for (size_t parIdx = 0; parIdx < (unsigned) fitting_space_pars.size(); ++parIdx) {
        auto mpar = _model_pars[parIdx];
        if (_par_modification_map.count(mpar) == 1) { // ...if this is a modified par
            // transform => rescale => update upars
            model_space_pars[parIdx] = _par_rescale_map[mpar]->rescale(
                _par_modification_map[mpar]->transform(fitting_space_pars[parIdx], fitting_space_pars)
            );
        }
    }
    return model_space_pars;
}

// TODO figure out the smarter Eigen-y version of this
Mat2D AbcSmc::_to_model_space(
    const Mat2D &fitting_space_pars
) {
    assert( _model_pars.size() == (unsigned) fitting_space_pars.cols() );
    Mat2D model_space_pars = fitting_space_pars; // copy initially - all model_space_pars == fitting_space_pars
    for (size_t parIdx = 0; parIdx < (unsigned) fitting_space_pars.cols(); ++parIdx) {
        auto mpar = _model_pars[parIdx];
        if (_par_modification_map.count(mpar) == 1) { // ...if this is a modified par
            for (size_t rowIdx = 0; rowIdx < (unsigned) fitting_space_pars.rows(); ++rowIdx) {
                // transform => rescale => update upars
                model_space_pars(rowIdx, parIdx) = _par_rescale_map[mpar]->rescale(
                    _par_modification_map[mpar]->transform(fitting_space_pars(rowIdx, parIdx), fitting_space_pars.row(rowIdx))
                );
            }
        }
    }
    return model_space_pars;
}

// Build DB if it doesn't exist;
// If it does exist but more sets are needed, filter particles and sample for next set;
// If the specified number of sets already exist, exit gracefully
bool AbcSmc::process_database(
    const gsl_rng* RNG,
    const size_t verbose
) {

    // if this is the initial build
    if (build(verbose)) {
        std::vector<size_t> posterior_ranks = {};
        // need to insert the initial job set
        const size_t setnum = 0;
        Mat2D pars = sample_priors(
            RNG, get_smc_size_at(setnum), *_posterior, _model_pars,
            posterior_ranks
        );
        Mat2D upars;
        if (true) { upars = _to_model_space(pars); }
        vector<size_t> seeds;

        return _storage.create_SMC_set(pars, seeds, posterior_ranks, verbose);
    } else {
        _particle_parameters.clear();
        _particle_metrics.clear();
        _weights.clear();
        _predictive_prior.clear();
        vector< vector<int> > serials;

    }

    cerr << std::setprecision(PREC);

    if (not load_SMC_sets(serials)) return false; // slurp sets & do particle filtering (identify pred prior) if needed
    const size_t next_set = serials.size();
    assert(next_set > 0);
    const size_t last_set = next_set - 1; // this set number

//    if ( _pred_prior_size = 0) set_predictive_prior_size(last_set);

    AbcLog::report_convergence_data(this, last_set);

    cerr << endl << endl;

    if (_num_smc_sets > next_set) {

        stringstream ss;
        const size_t num_particles = get_smc_size_at(next_set);

        Mat2D noised_pars;
        const int last_serial = serials.back().back();

        switch (_noise) {
            case NOISE::MULTIVARIATE: {
                gsl_matrix* L = setup_mvn_sampler(
                    (*(_particle_parameters)[next_set-1])(_predictive_prior[next_set-1], Eigen::placeholders::all)
                );
                noised_pars = ABC::sample_mvn_predictive_priors(
                    RNG, num_particles, 
                    *(_weights[next_set-1]),
                    (*(_particle_parameters[next_set-1]))(_predictive_prior[next_set-1], Eigen::placeholders::all),
                    _model_pars,
                    L
                );
                gsl_matrix_free(L);
                if (verbose) std::cerr << "Populating next set using MULTIVARIATE noising of parameters." << std::endl;
                break;
            }
            case NOISE::INDEPENDENT: {
                noised_pars = ABC::sample_predictive_priors(
                    RNG, num_particles, 
                    *(_weights[next_set-1]),
                    (*(_particle_parameters[next_set-1]))(_predictive_prior[next_set-1], Eigen::placeholders::all),
                    _model_pars,
                    *(_doubled_variance[next_set-1])
                );
                if (verbose) std::cerr << "Populating next set using INDEPENDENT noising of parameters." << std::endl;
                break;
            }
            default: std::cerr << "Unknown noise type.  Aborting." << std::endl; exit(-1);
        }
        
        db.Query("BEGIN EXCLUSIVE;").Next();

        for (size_t i = 0; i < (unsigned) noised_pars.rows(); i++) {
            const int serial = last_serial + 1 + i;

            QueryStr qstr;

            ss << "insert into " << JOB_TABLE << " values ( " << serial << ", "
                                               << next_set << ", "
                                               << i << ", "
                                               << time(NULL)
                                               << ", NULL, 'Q', -1, 0 );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);

            const size_t seed = gsl_rng_get(RNG);

            ss << "insert into " << PAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << noised_pars(i, j); ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);

            if (_par_modification_map.size()) {
                const Row upars = _to_model_space(noised_pars.row(i));
                ss << "insert into " << UPAR_TABLE << " values ( " << serial << ", '" << seed << "'"; for (size_t j = 0; j < npar(); j++)  ss << ", " << upars[j]; ss << " );";
                //cerr << "attempting: " << ss.str() << endl;
                _db_execute_stringstream(db, ss);
            }

            ss << "insert into " << MET_TABLE << " values ( " << serial; for (size_t j = 0; j < nmet(); j++) ss << ", NULL"; ss << " );";
            //cerr << "attempting: " << ss.str() << endl;
            _db_execute_stringstream(db, ss);
        }
        db.CommitTransaction();
        
    } else {
        cerr << "Database already contains " << _num_smc_sets << " complete sets.\n";
    }

    return true;
}

// Read in existing sets and do particle filtering as appropriate
bool AbcSmc::load_SMC_sets(vector< vector<int> > &serials) {
    _storage.read_SMC_sets(
        serials, _particle_parameters, _particle_metrics
    );
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
        const int json_set_size = get_smc_size_at(t);
        if (set_size != json_set_size) {
            cerr << "ERROR:\tSet size for one or more sets does not agree between configuration file and database:" << endl
                 << "\tSet " << t << " in configuration file has size " << json_set_size << " vs size " << set_size << " in database." << endl;
//                 << "\tNB: You may want to edit the configuration file to have an array of set sizes that reflect what is already in the database." << endl;
            set_size_test_success = false;
            break;
            //return false;
        }
        _particle_parameters.push_back(std::make_shared<Mat2D>(completed_set_size, npar()));
        _particle_metrics.push_back(std::make_shared<Mat2D>(completed_set_size, nmet()));
        // join all three tables for rows with smcSet = t, slurp and store values
        string select_str = "select J.serial, J.particleIdx, J.posterior, " + _build_sql_select_par_string("") + ", " + _build_sql_select_met_string()
                            + "from " + JOB_TABLE + " J, " + MET_TABLE + " M, " + PAR_TABLE + " P where J.serial = M.serial and J.serial = P.serial "
                            + "and J.smcSet = " + to_string((long long) t) + ";";

        serials.push_back( vector<int>(completed_set_size, 0) );

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
            for(size_t i = offset; i < offset + npar(); i++) (_particle_parameters[t])->operator()(particle_counter, i-offset) = (double) s2.GetField(i);
            offset += npar();
            for(size_t i = offset; i < offset + nmet(); i++) (_particle_metrics[t])->operator()(particle_counter, i-offset) = (double) s2.GetField(i);
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
                    particle_ranking_PLS(*(_particle_metrics[t]), *(_particle_parameters[t]), *_met_vals, _pls_training_fraction)
                ); break; }
                case FILTER::SIMPLE: { _predictive_prior.push_back(
                    particle_ranking_simple(*(_particle_metrics[t]), *(_particle_parameters[t]), *_met_vals )
                ); break; }
                default: std::cerr << "ERROR: Unsupported filtering method: " << _filtering << std::endl; return false;
            }

            // trim that ranking to the size of the next predictive prior (n.b. this only drops from the indexing vector, not the actual data)
            const size_t next_pred_prior_size = get_pred_prior_size_at(t);
            _predictive_prior.back().resize(next_pred_prior_size);

            auto posterior_pars = (_particle_parameters[t])->operator()(_predictive_prior[t], Eigen::placeholders::all);
            auto posterior_mets = (_particle_metrics[t])->operator()(_predictive_prior[t], Eigen::placeholders::all);

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

bool AbcSmc::_run_simulator(Row &par, Row &met, const size_t rng_seed, const size_t serial) {
    vector<float_type> met_vec = (*_simulator)( as_vector(par), rng_seed, serial, _mp );
    bool particle_success = (met_vec.size() == nmet());
    if (!particle_success) {
        cerr << "ERROR: simulator function returned the wrong number of metrics: expected " << nmet() << ", received " << met_vec.size() << endl;
    }
    met = as_row(met_vec);
    return particle_success;
}


// TODO - these could likely be refactored to a single, more clever build sql string function

string AbcSmc::_build_sql_create_par_string( string tag = "" ) {
    stringstream ss;
    for (size_t i = 0; i < npar()-1; i++) { ss << _model_pars[i]->get_short_name() << tag << " real, "; }
    ss << _model_pars.back()->get_short_name() << tag << " real ";
    return ss.str();
}


string AbcSmc::_build_sql_create_met_string( string tag = "" ) {
    stringstream ss;
    for (size_t i = 0; i < nmet()-1; i++) { ss << _model_mets[i]->get_short_name() << tag << " real, "; }
    ss << _model_mets.back()->get_short_name() << tag << " real ";
    return ss.str();
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


bool AbcSmc::_db_execute_strings(sqdb::Db &db, vector<string> &update_buffer) {
    bool db_success = false;
    try {
        db.Query("BEGIN EXCLUSIVE;").Next();
        for (size_t i = 0; i < update_buffer.size(); ++i) {
            db.Query(update_buffer[i].c_str()).Next();
        }
        db_success = true;
        db.CommitTransaction();
    } catch (const Exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        for (size_t i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    } catch (const exception &e) {
        db.RollbackTransaction();
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        for (size_t i = 0; i < update_buffer.size(); ++i) cerr << update_buffer[i] << endl;
    }
    return db_success;
}


bool AbcSmc::_db_execute_stringstream(sqdb::Db &db, stringstream &ss) {
    // We don't need BEGIN EXCLUSIVE here because the calling function has already done it
    bool db_success = false;
    try {
        db_success = db.Query(ss.str().c_str()).Next();
    } catch (const Exception &e) {
        cerr << "CAUGHT E: ";
        cerr << e.GetErrorMsg() << endl;
        cerr << "Failed query:" << endl;
        cerr << ss.str() << endl;
    } catch (const exception &e) {
        cerr << "CAUGHT e: ";
        cerr << e.what() << endl;
        cerr << "Failed query:" << endl;
        cerr << ss.str() << endl;
    }
    ss.str(string());
    ss.clear();
    return db_success;
}

bool AbcSmc::build(
    const size_t verbose
) {

    // is storage *already* setup?
    if (_storage.is_setup()) return false;

    return _storage.setup(
        _model_pars, _model_mets,
        not _par_modification_map.empty(),
        verbose
    );

}

bool AbcSmc::simulate_next_particles(
    const int n, const int serial_req, const int posterior_req
) {
    return false;
}

bool AbcSmc::_do_work(
    const Mat2D &pars,
    const Coli &serials,
    const Colsz &seeds,
    const bool keep_going,
    const bool bulk,
    const size_t verbose
) {
    for (auto par_idx = 0; par_idx < pars.rows(); par_idx++) {
        const int serial = serials[par_idx];
        const size_t seed = seeds[par_idx];
        Row pars = pars.row(par_idx);
        Row mets;
        const size_t start_time = duration_cast<seconds>(high_resolution_clock::now().time_since_epoch()).count();
        bool success = _run_simulator(pars, mets, seed, serial);
        if (success) {
            _storage->put_work(mets, serial, start_time);
        } else if (not keep_going) {
            exit(-211);
        }
    }

}

bool AbcSmc::simulate(
    const int n,
    const bool keep_going,
    const bool bulk, // TODO: allow bulk writes for fast simulations?
    const size_t verbose
) {
    Mat2D pars;
    Coli serials;
    Colsz seeds;

    if (_storage->get_work(pars, seeds, serials, n)) {
        return _do_work(pars, serials, seeds, keep_going, bulk, verbose);
    }

    return false;
}

bool AbcSmc::simulate_serials(
    const vector<int> &serials,
    const bool keep_going,
    const bool bulk, // TODO: allow bulk writes for fast simulations?
    const size_t verbose
) { 
    Mat2D pars;
    Coli sercol = Coli(sercol.size());
    for (auto serial_idx = 0; serial_idx < serials.size(); serial_idx++) {
        sercol[serial_idx] = serials[serial_idx];
    }
    Colsz seeds;

    if (_storage->get_work(pars, seeds, sercol)) {
        return _do_work(pars, sercol, seeds, keep_going, bulk, verbose);
    } else { // read error fetching work
        // TODO verbosity? or leave that to storage method?
        return false;
    }

    return true;
}

bool AbcSmc::simulate_posterior(
    const vector<size_t> &posterior,
    const int smc_set,
    const bool keep_going,
    const bool bulk, // TODO: allow bulk writes for fast simulations?
    const size_t verbose
) { 
    Coli sercol;
    _storage->get_posterior_serials(posterior, sercol, smc_set);

    if (_storage->get_work(pars, seeds, sercol)) {
        return _do_work(pars, sercol, seeds, keep_going, bulk, verbose);
    } else { // read error fetching work
        // TODO verbosity? or leave that to storage method?
        return false;
    }

    return true;
}

void AbcSmc::calculate_predictive_prior_weights(
    const size_t set_num
) {
    assert(_doubled_variance.size() == set_num); // current length of doubled variance == index of previous set

    auto setview = (_particle_parameters[set_num])->operator()(_predictive_prior[set_num], Eigen::placeholders::all);
    _doubled_variance.push_back(
        ABC::calculate_doubled_variance(setview)
    );

    if (set_num == 0) {
        _weights.push_back(
            ABC::weight_predictive_prior(
                _model_pars, setview
            )
        );
    } else if ( set_num > 0 ) {
        _weights.push_back(
            ABC::weight_predictive_prior(
                _model_pars,
                setview,
                (_particle_parameters[set_num-1])->operator()(_predictive_prior[set_num-1], Eigen::placeholders::all),
                *(_weights[set_num - 1]),
                *(_doubled_variance[set_num -1])
            )
        );
    }
}

void AbcSmc::_particle_scheduler_mpi(const size_t t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {

    auto _num_particles = get_smc_size_at(t);

    // sample parameter distributions; copy values into Y matrix and into send_data buffer
    Row par_row;
    if (t == 0) {
        std::vector<size_t> dump(_num_particles);
        Y_orig = sample_priors(RNG, _num_particles, *_posterior, _model_pars, dump);
    } else {
        Y_orig = ABC::sample_predictive_priors(
            RNG, _num_particles, 
            *(_weights[t-1]),
            (_particle_parameters[t-1])->operator()(_predictive_prior[t-1], Eigen::placeholders::all),
            _model_pars,
            *(_doubled_variance[t-1])
        );
    }

    X_orig.resize(_num_particles, nmet());

    ABC::particle_scheduler(X_orig, Y_orig, _mp);

}

void AbcSmc::_particle_worker_mpi(
    const size_t seed,
    const size_t serial
) {
    ABC::particle_worker(npar(), nmet(), _simulator, seed, serial, _mp);
}
