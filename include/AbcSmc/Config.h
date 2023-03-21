
#ifndef ABCSMC_CONFIG_H
#define ABCSMC_CONFIG_H

#include <vector>
#include <AbcSmc/TypeDefs.h>
#include <json/json.h>

template <NumericType T>
std::vector<T> as_vector(const Json::Value & val) {
    vector<T> extracted_vals;
    if (val.isArray()) { for (const Json::Value & jv : val) {
        extracted_vals.push_back( jv.as<T>() ); // NB, jsoncpp handles cast failures
    } } else {
        extracted_vals.push_back( val.as<T>() );
    }
    return extracted_vals;
}

namespace ABC {

struct Config {
    virtual void parse_iterations(
        const size_t pseudosize,
        size_t * iterations,
        float_type * training_frac,
        std::vector<size_t> * set_sizes,
        std::vector<size_t> * pred_prior_sizes
    ) const = 0;
};

struct JsonConfig : public Config {
    Json::Value par;
    JsonConfig(const std::string & filename);
    
    void parse_iterations(
        const size_t pseudosize,
        size_t * iterations,
        float_type * training_frac,
        std::vector<size_t> * set_sizes,
        std::vector<size_t> * pred_prior_sizes
    ) const override;

    private:
        const Json::Value _root;
};

}

#endif // ABCSMC_CONFIG_H