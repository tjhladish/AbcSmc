// This header covers translating the JSON configuration file into AbcSmc object
#ifndef ABCJSON_H
#define ABCJSON_H

#include "AbcSmc.h"
#include "EnumMacros.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <string>
#include <sstream>
#include <fstream>

using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;
using std::setw;

using namespace ABC;

struct AbcJson : AbcBuilder {
    
}

bool file_exists(const char *fileName) {
    std::ifstream infile(fileName);
    return infile.good();
}

template<typename NUMTYPE>
NUMTYPE extract(const Json::Value vv) { 
    if constexpr (std::is_integral_v<NUMTYPE>) {
        return v.asInt();
    } else {
        return v.asFloat();
    }
};

template<typename NUMTYPE>
bool check(const Json::Value vv) { 
    if constexpr (std::is_integral_v<NUMTYPE>) {
        return v.isIntegral();
    } else {
        return v.isDouble();
    }
};

template<typename NUMTYPE>
vector<NUMTYPE> as_vector(Json::Value val, string key) {
    vector<NUMTYPE> extracted_vals;
    auto vv = val[key];
    string typeword;
    if constexpr (std::is_integral_v<NUMTYPE>) {
        typeword = "integers";
    } else {
        typeword = "floats";
    }
    if ( vv.isArray() ) {
        for ( unsigned int i = 0; i < vv.size(); ++i) extracted_vals.push_back( extract(vv[i]) );
    } else if ( check(vv) ) {
        extracted_vals.push_back( extract(vv) );
    } else {
        cerr << "Could not parse " << key << " as " << typeword << " or array of " << typeword << "." << endl;
        exit(-216);
    }
    return extracted_vals;
};

void AbcSmc::process_predictive_prior_arguments(Json::Value par) {
    int arg_ct = par.isMember("predictive_prior_fraction") + par.isMember("predictive_prior_size");
    if (arg_ct == 0) {
        cerr << "Error: either predictive_prior_fraction or predictive_prior_size must be specified in configuration file." << endl;
        exit(1);
    } else if (arg_ct == 2) {
        cerr << "Error: only one of predictive_prior_fraction and predictive_prior_size may be specified in configuration file." << endl;
        exit(1);
    } else if (par.isMember("predictive_prior_fraction")) {
        auto ppfs = as_vector<float>(par, "predictive_prior_fraction");
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
        _predictive_prior_sizes = as_vector<int>(par, "predictive_prior_size");
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

template<typename NT>
vector<NT> sequence(const Json::Value pardata) {
    // if the parameter provides the options, use them
    if (pardata.isMember("pars")) {
        if (pardata.isMember("par1") or pardata.isMember("par2") or pardata.isMember("step") or pardata.isMember("incs")) {
            std::cerr << "Warning: using `pars` parameter, but other sequence options provided." << std::endl;
        }
        return as_vector<NT>(pardata["pars"]);
    } else {
        NT st = extract(pardata["par1"]);
        vector<NT> res;
        NT step;
        size_t increments;

        if (pardata.isMember("step") or pardata.isMember("incs")) {
            // figure out number of steps & step size
            // then sequence = start + 0:N * del
        } else if (pardata.isMember("par1") and pardata.isMember("par2")) {
            assert(par1 < par2);
            NT del = par2 - par1;
            // assume step size of one
        }
        std::cerr << "Failed to parse sequence." << std::endl;
        return res;
    }
};

Parameter * parse_transform(Parameter * base, const Json::Value transdata) {
    double (*_untransform_func)(const double, const vector<double>);
    vector<double> ps;
    if (transdata.isString()) {
        TransformType tr = from_string(transdata.asString());
        if (tr == NONE) {
            _untransform_func = [](const double t, const vector<double> /**/) { return t; };
        } else if (tr == POW_10) {
            _untransform_func = [](const double t, const vector<double> ps) { return pow(ps[0], t); };
            ps.push_back(10.0);
        } else if (tr == LOGISTIC) {
            _untransform_func = [](const double t, const vector<double> ps) { return ABC::logistic(t); };
            ps.push_back(10.0);
        }
    } else if (transdata.isObject()) {
        TransformType tr = from_string(transdata["type"].asString());
        if (tr != LOGISTIC) {
            cerr << "Only type: LOGISTIC is currently supported for untransformation objects.  (NONE and POW_10 supported as untransformation strings.)\n";
            exit(-207);
        }
        // TODO: modifications to ps support this - feels like this is going to be breaking?
        // if inner shifts, do those first?
        _untransform_func = [](const double t, const vector<double> ps) { return ABC::logistic(t); };

    } else {
        std::cerr << "Did not understand transformation object." << std::endl;
        exit(-1);
    }

    return new TransformedParameter(base, _untransform_func, ps);
}

Parameter * parse_parameter(const Json::Value pardata) {
    string name = pardata["name"].asString();
    string short_name = pardata.get("short_name", "").asString();
    // these provide their own parsing errors
    NumericType ntype = from_string(pardata["num_type"].asString());
    ParameterType ptype = from_string(pardata["dist_type"].asString());
    Parameter * res;
    if (ptype == UNIFORM) {
        if (ntype == INT) {
            res = new UniformPrior<int>(name, short_name, pardata["par1"].asInt(), pardata["par2"].asInt());
        } else if (ntype == FLOAT) {
            res = new UniformPrior<double>(name, short_name, pardata["par1"].asDouble(), pardata["par2"].asDouble());
        }
    } else if ((ptype == NORMAL) or (ptype == GAUSSIAN)) {
        res = new GaussianPrior(name, short_name, pardata["par1"].asDouble(), pardata["par2"].asDouble());
    } else if ((ptype == PSEUDO) or (ptype == POSTERIOR)) {
        if (ntype == INT) {
            res = new PseudoParameter<int>(name, short_name, sequence<int>(pardata), (ptype == POSTERIOR));
        } else if (ntype == FLOAT) {
            res = new PseudoParameter<double>(name, short_name, sequence<double>(pardata), (ptype == POSTERIOR));
        }            
    }
    if (!pardata.isMember("untransform")) {
        return res;
    } else {
        return parse_transform(res, pardata);
    }
};

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

    set_smc_iterations( par["smc_iterations"].asInt() ); // TODO: or have it test for convergence
    _smc_set_sizes = as_vector<int>(par, "num_samples");
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
    if (noise != "INDEPENDENT" and noise != "MULTIVARIATE") {
        cerr << "Unknown parameter noise type specified: " << noise << ". Aborting." << endl;
        exit(-210);
    }
    bool use_mvn_noise = (noise == "MULTIVARIATE");

    // Parse model parameters
    const Json::Value model_par = par["parameters"];
    map<string, int> par_name_idx;
    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Build name lookup
        string name = model_par[i]["name"].asString();
        assert(par_name_idx.count(name) == 0);
        par_name_idx[name] = i;
    }

    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Iterates over the sequence elements.
        add_next_parameter(parse_parameter(model_par[i]));
    }

    // Parse model metrics
    const Json::Value model_met = par["metrics"];
    for ( unsigned int i = 0; i < model_met.size(); ++i )  {// Iterates over the sequence elements.
        string name = model_met[i]["name"].asString();
        string short_name = model_met[i].get("short_name", "").asString();

        NumericType ntype = from_string(model_met[i]["num_type"].asString());

        double val = model_met[i]["value"].asDouble();

        add_next_metric(new Metric(name, short_name, ntype, val));
    }

    return true;
}

#endif // ABCJSON_H