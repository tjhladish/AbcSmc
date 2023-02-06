
#ifndef ABCJSON_H
#define ABCJSON_H

#include <vector>
#include <type_traits>
#include <iostream>
#include <json/json.h>
#include <fstream>

Json::Value read_config(const std::string & conf_filename) {

    std::ifstream ifs(conf_filename.c_str());
    if (not ifs.good()) {
        std::cerr << "Config file is not good: " << conf_filename << std::endl;
        exit(1);
    }

    Json::Value root;
    Json::Reader reader;

    if (not reader.parse(ifs, root)) {
        std::cerr << "Failed to parse JSON file: " << conf_filename << std::endl;
        std::cerr << reader.getFormattedErrorMessages() << std::endl;
        exit(-214);
    }

    return root;
};

template<typename NUMTYPE>
NUMTYPE extract(const Json::Value & vv) { 
    if constexpr (std::is_integral_v<NUMTYPE>) {
        return vv.asInt();
    } else {
        return vv.asFloat();
    }
};

template<typename NUMTYPE>
bool check(const Json::Value & vv) { 
    if constexpr (std::is_integral_v<NUMTYPE>) {
        return vv.isIntegral();
    } else {
        return vv.isNumeric();
    }
};

template<typename NUMTYPE>
void parse_error(const std::string & key) {
    std::string typeword;
    if constexpr (std::is_integral_v<NUMTYPE>) {
        typeword = "integer(s)";
    } else {
        typeword = "float(s)";
    }
    std::cerr << "Could not parse " << key << " as " << typeword << "." << std::endl;
    exit(-216);
};

template<typename NUMTYPE>
std::vector<NUMTYPE> as_vector(const Json::Value & val, const std::string & key) {
    vector<NUMTYPE> extracted_vals;
    auto vv = val[key];
    if ( vv.isArray() ) {
        for ( unsigned int i = 0; i < vv.size(); ++i) extracted_vals.push_back( extract<NUMTYPE>(vv[i]) );
    } else if ( check<NUMTYPE>(vv) ) {
        extracted_vals.push_back( extract<NUMTYPE>(vv) );
    } else {
        parse_error<NUMTYPE>(key);
    }
    return extracted_vals;
};

template<typename NUMTYPE>
NUMTYPE parse_pseudo_step(const Json::Value & par) {
    if (par.isMember("step")) {
        if (!check<NUMTYPE>(par["step"])) {
            parse_error<NUMTYPE>("step");
        }
        auto res = extract<NUMTYPE>(par["step"]);
        if (res <= 0) {
            std::cerr << "You must specify `step` > 0 for PSEUDO parameters." << std::endl;
        }
    } else {
        return static_cast<NUMTYPE>(1);
    }
};

template<typename NUMTYPE>
size_t parse_pseudo_steps(const Json::Value & par, const NUMTYPE & step) {
    if (not (par.isMember("par1") and par.isMember("par2"))) {
        std::cerr << "You must specify `par1` and `par2` OR `values` for PSEUDO parameters." << std::endl;
        exit(-217);
    } else {
        auto par1 = extract<NUMTYPE>(par["par1"]);
        auto par2 = extract<NUMTYPE>(par["par2"]);
        if (not (par1 <= par2)) {
            std::cerr << "You must specify `par1` <= `par2` for PSEUDO parameters." << std::endl;
            exit(-218);
        } else if (par1 == par2) {
            return 0;
        } else {
            return static_cast<size_t>(std::ceil((par2 - par1) / step) + 1);
        }
    }
};

template<typename NUMTYPE>
std::vector<NUMTYPE> parse_pseudo_values(const Json::Value & par) {
    if (par.isMember("values")) {
        return as_vector<NUMTYPE>(par, "values");
    } else {
        std::vector<NUMTYPE> res;
        auto step = parse_pseudo_step<NUMTYPE>(par);
        const size_t nsteps = parse_pseudo_steps<NUMTYPE>(par, step);
        const auto par1 = extract<NUMTYPE>(par["par1"]);
        for (size_t i = 0; i < nsteps; i++) {
            res.push_back(par1 + i*step);
        }
        return res;
    };
};

// parse the predictive prior sizes by wave
// in projection mode, may be empty; 
std::vector<size_t> parse_predictive_prior(
    const Json::Value & root, const std::vector<size_t> & _smc_set_sizes,
    const bool verbose = false
) {

    std::vector<size_t> _predictive_prior_sizes;

    if (root.isMember("predictive_prior_fraction") and root.isMember("predictive_prior_size")) {
        std::cerr << "Error: Only one of predictive_prior_fraction or predictive_prior_size may be specified in configuration file." << std::endl;
        exit(1);
    } else if (root.isMember("predictive_prior_fraction")) {
        if (_smc_set_sizes.size() == 0) {
            std::cerr << "When no _smc_set_sizes, in projection mode, and predictive_prior_fraction cannot be specified." << std::endl;
            exit(1);
        }
        auto ppfs = as_vector<float>(root, "predictive_prior_fraction");

        vector<size_t> set_sizes_copy = _smc_set_sizes;
        const size_t max_set = max(ppfs.size(), set_sizes_copy.size());
        
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
    } else if (root.isMember("predictive_prior_size")) {
        _predictive_prior_sizes = as_vector<size_t>(root, "predictive_prior_size");
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
    } else {
        if (verbose) {
            std::cerr << "Neither `predictive_prior_fraction` or `predictive_prior_size` specified; assuming projection mode." << std::endl;
        }
    }
    return _predictive_prior_sizes;
}

#endif // ABCJSON_H