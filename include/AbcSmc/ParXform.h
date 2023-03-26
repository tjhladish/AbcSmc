#ifndef ABCSMC_PARXFORM_H
#define ABCSMC_PARXFORM_H

#include <vector>
#include <cstddef>

namespace ABC {

typedef float_type transformer(const float_type &);

// this defines shift/scale transformation, from _fitting_ space (a/k/a transformed) to _model_ space (a/k/a untransformed).
//
// a = sum of transformed scale addends
// b = prod of transformed scale factors
// c = sum of untransformed scale addends
// d = prod of untransformed scale factors
// u = arbitrary untransform function
// x = value on the _fitting_ scale
// x' = value on the _model_ scale
//
// x' = (u((x + a) * b) + c) * d (i.e. value on the model/natural scale)
//
// When transforming, expects to be passed a vector of _fitting_ space values.
// These combined according to specification for addends and factors to produce (a,b,c,d).
struct ParXform {

    ParXform() = delete;

    ParXform(
        const transformer * u,
        const std::vector<size_t> &fitting_space_add_idx = {},
        const std::vector<size_t> &fitting_space_mul_idx = {},
        const std::vector<size_t> &model_space_add_idx = {},
        const std::vector<size_t> &model_space_mul_idx = {}
    ) : _u(u), _idx_fitting_space_addends(fitting_space_add_idx),
        _idx_fitting_space_factors(fitting_space_mul_idx),
        _idx_model_space_addends(model_space_add_idx),
        _idx_model_space_factors(model_space_mul_idx) {    
    }

    // T is some sort of random access container of float_type
    // in AbcSmc, `Row`s of parameters, generally - but could be any container
    // this enables testing with simpler containers
    template <typename T>
    float_type transform(
        const float_type &pval, const T &fitting_space_values
    ) const {
        float_type tplus = 0.0, ttimes = 1.0, uplus = 0.0, utimes = 1.0;
        for (auto i : _idx_fitting_space_addends) { tplus += fitting_space_values[i]; }
        for (auto i : _idx_fitting_space_factors) { ttimes *= fitting_space_values[i]; }
        for (auto i : _idx_model_space_addends) { uplus += fitting_space_values[i]; }
        for (auto i : _idx_model_space_factors) { utimes *= fitting_space_values[i]; }
        return (_u((pval + tplus)*ttimes) + uplus)*utimes;
    }

    private:
        const transformer * _u;
        const std::vector<size_t> _idx_fitting_space_addends;
        const std::vector<size_t> _idx_fitting_space_factors;
        const std::vector<size_t> _idx_model_space_addends;
        const std::vector<size_t> _idx_model_space_factors;
};

struct ParRescale {
    ParRescale(const float_type p1 = 0.0, const float_type p2 = 1.0) : par1(p1), par2(p2) {}
    const float_type par1, par2;
    float_type rescale(const float_type &pval) const { return (par2 - par1)*pval + par1; }
};

// TODO: options for simpler transformations / rescales
struct FunXform {
    FunXform(const transformer * u) : _u(u) {}

    template <typename T>
    float_type transform(
        const float_type &pval, const T &/* fitting_space_values */
    ) const { return _u(pval); }

    private:
        const transformer * _u;
};

struct NonRescale {
    float_type rescale(const float_type &pval) const { return pval; }
};

}

#endif // ABCSMC_PARXFORM_H