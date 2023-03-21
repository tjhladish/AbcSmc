
#include <AbcSmc/Config.h>
#include <iostream>
#include <math.h>

using std::cerr;
using std::endl;

namespace ABC { 

void JsonConfig::parse_iterations(
    const size_t pseudosize,
    size_t * iterations,
    float_type * training_frac,
    std::vector<size_t> * set_sizes,
    std::vector<size_t> * pred_prior_sizes
) const {
    if (pseudosize != 0) { // definitely in projection mode; ignore other parameters
        if (_root.get("smc_iterations", 1).asInt() != 1) { // ignoring OR erroring if specified != 1
            cerr << "Cannot use smc_iterations > 1 with ONLY PSEUDO or POSTERIOR parameters.  Aborting." << endl;
            exit(-202);
        }
        if (_root.isMember("num_samples")) {
            size_t checksize = (as_vector<size_t>(par["num_samples"]))[0];
            if (checksize != pseudosize) {
                std::cerr << "ERROR: `num_samples` ("<< checksize <<") does not match imputed combinations of PSEUDO and/or POSTERIOR parameters ("<< pseudosize <<")." << endl;
                exit(-201);
            } else {
                std::cerr << "WARNING: specified `num_samples` for all PSEUDO and/or POSTERIOR parameters." << std::endl;
            }
        }
        if (_root.isMember("predictive_prior_fraction") or _root.isMember("predictive_prior_size")) {
            std::cerr << "WARNING: ignoring `predictive_prior_*` options in projection mode." << std::endl;
        }
        *iterations = 1;
        *set_sizes = { pseudosize };
    } else { // some priors => fitting mode => predictive prior fractions
        // pred_prior_sizes is a series of sizes, but may be passed as EITHER fractions OR sizes.
        // both / neither is an error
        if (
            (_root.isMember("predictive_prior_fraction") and _root.isMember("predictive_prior_size")) or
            not (_root.isMember("predictive_prior_fraction") or _root.isMember("predictive_prior_size"))
        ) {
            std::cerr << "Error: exactly one of `predictive_prior_fraction` or `predictive_prior_size` must be specified in configuration file." << std::endl;
            exit(1);
        }

        // if doing iterations, may have a training fraction; defaults to 0.5
        *training_frac = _root.get("pls_training_fraction", 0.5).as<float_type>();
        if (*training_frac <= 0 or 1 <= *training_frac) {
            std::cerr << "Error: pls_training_fraction must be in (0, 1)." << std::endl;
            exit(1);
        }

        *set_sizes = as_vector<size_t>(par["num_samples"]);

        if (_root.isMember("predictive_prior_fraction")) {
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
            for (size_t i = 0; i < max_set; ++i) { (*pred_prior_sizes)[i] = round(ppfs[i] * set_sizes_copy[i]); }
        } else { // if (_root.isMember("predictive_prior_size")) {
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

        *iterations = _root.get("smc_iterations", max(set_sizes->size(), pred_prior_sizes->size())).asUInt64();

    }
  
}


}