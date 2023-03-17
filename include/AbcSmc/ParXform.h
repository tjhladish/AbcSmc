#ifndef ABCSMC_PARXFORM_H
#define ABCSMC_PARXFORM_H

#include <AbcSmc/Parameter.h>
#include <vector>

namespace ABC {

// this defines shift/scale transformation.
//
// a = sum of transformed scale addends
// b = prod of transformed scale factors
// c = sum of untransformed scale addends
// d = prod of untransformed scale factors
// u = arbitrary untransform function
// x = value on the transformed scale (i.e. the fitting scale)
//
// transform(x) = (u((x + a) * b) + c) * d (i.e. value on the model/natural scale)
struct ParXform {

    ParXform(
        const float_type tp = 0.0, const float_type up = 0.0,
        const float_type tt = 1.0, const float_type ut = 1.0
    ) : tplus(tp), uplus(up), ttimes(tt), utimes(ut) {}

    template <typename T> 
    ParXform(const T & pre_shift, const T & post_shift, const T & pre_scale, const T & post_scale) :
    ParXform(
        std::accumulate(pre_shift.begin(), pre_shift.end(), 0.0),
        std::accumulate(post_shift.begin(), post_shift.end(), 0.0),
        std::accumulate(pre_scale.begin(), pre_scale.end(), 1.0, std::multiplies<float_type>()),
        std::accumulate(post_scale.begin(), post_scale.end(), 1.0, std::multiplies<float_type>())
    ) { }

    const float_type tplus, uplus, ttimes, utimes;

    float_type transform(const float_type & pval, float_type (*u)(const float_type &)) const {
        return (u((pval + tplus)*ttimes) + uplus)*utimes;
    }

    private:
        std::vector<size_t> _idx_fitting_space_addends;
        std::vector<size_t> _idx_fitting_space_factors;
        std::vector<size_t> _idx_model_space_addends;
        std::vector<size_t> _idx_model_space_factors;

};

}

#endif // ABCSMC_PARXFORM_H