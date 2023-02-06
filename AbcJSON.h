
#ifndef ABCJSON_H
#define ABCJSON_H

#include <vector>
#include <type_traits>
#include <iostream>

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

#endif // ABCJSON_H