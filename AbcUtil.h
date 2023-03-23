#ifndef ABCUTIL_H
#define ABCUTIL_H

#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "ranker.h"
#include <math.h>

#include "PLS/pls.h"

using namespace std;

// forward declare to avoid cyclic AbcUtil <=> AbcSmc dependency
// TODO: move Parameter into own class
class Parameter;

namespace ABC {

    struct LinearFit {
        double m;
        double b;
        double rsq;
    };

    struct LogisticFit {
        double beta0;
        double beta1;
        double simplex_size;
        int status;
        int iterations;
    };

    //using namespace std;
    using namespace Eigen;

    std::string slurp(std::string filename);
    std::string get_nth_line(const std::string& filename, int N);

    inline float_type logit(const float_type p) { assert((0.0 <= p) and (p <= 1.0)); return log( p / (1.0 - p) ); }
    inline float_type logistic(const float_type l) { return 1.0 / (1.0 + exp(-l)); }

    //int _sgn(float_type val) { return (0 < val) - (val < 0); }

    float_type median(const Col & data);

    float_type quantile(const Col & data, double q);

    float_type variance(const Col & data, const float_type _mean);

    float_type max(const Col & data);

    float_type skewness(const Col & data);

    float_type median_crossings(const Col & data, const float_type _median);

    float_type median_crossings(const Col & data);

    float_type optimize_box_cox(const Col & data, const float lambda_min, const float lambda_max, const float step);

    float_type optimize_box_cox(const Col & data);

    template <typename T>
    inline void cerr_vector(std::vector<T> & my_vector, std::string sep = " ") {
        for (size_t i = 0; i < my_vector.size() - 1; i++ ) std::cerr << my_vector[i] << sep;
        std::cerr << my_vector.back();
    }

    template <typename T>
    inline void cout_vector(std::vector<T> & my_vector, std::string sep = " ") {
        for (size_t i = 0; i < my_vector.size() - 1; i++ ) std::cout << my_vector[i] << sep;
        std::cout << my_vector.back();
    }

    inline double uniform_pdf(double a, double b) { return 1.0 / fabs(b-a); }

    std::vector<size_t> gsl_rng_nonuniform_int(const gsl_rng* RNG, const size_t num_samples, const Col & weights);

    Row gsl_ran_trunc_normal(
        const gsl_rng* RNG,
        const std::vector<Parameter*> _model_pars,
        const Row & mu, const Row & sigma_squared
    );

    Row gsl_ran_trunc_mv_normal(
        const gsl_rng* RNG,
        const vector<Parameter*> _model_pars,
        const Row & mu,
        const gsl_matrix* L
    );

    LinearFit* lin_reg(const std::vector<double> &x, const std::vector<double> &y);

    LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector< pair<int,int> > &y);

    LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector<int> &successes, const std::vector<int> &attempts);

    inline vector<float_type> as_vector(const Row data) {
        vector<float_type> vec(data.size());
        for (size_t i = 0; i < static_cast<size_t>(data.size()); i++) vec[i] = data[i];
        return vec;
    }

    inline Row as_row(const vector<float_type> & data) {
        Row row(data.size());
        for (size_t i = 0; i < data.size(); i++) row[i] = data[i];
        return row;
    }

    double calculate_nrmse(
        const Mat2D & posterior_mets,
        const Row & observed
    );

    Mat2D sample_posterior(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights,
        const Mat2D & posterior
    );

    template<typename RandomAccessible>
    gsl_vector* to_gsl_v(const RandomAccessible & from);

    gsl_matrix* to_gsl_m(const Mat2D & from);

    Col euclidean(const Mat2D & sims, const Row & ref);

    Mat2D sample_predictive_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights, const Mat2D & parameter_prior,
        const std::vector<Parameter*> & pars,
        const Row & doubled_variance
    );

    gsl_matrix* setup_mvn_sampler(
        const Mat2D & params
    );

    Mat2D sample_mvn_predictive_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights, const Mat2D & parameter_prior,
        const std::vector<Parameter*> & pars,
        const gsl_matrix* L
    );

    Mat2D sample_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Mat2D & posterior, // look up table for POSTERIOR type Parameters
        const std::vector<Parameter*> & mpars,
        std::vector<size_t> & posterior_ranks // filled in by this
    );

    std::vector<size_t> particle_ranking_simple (
        const Mat2D &X_orig, const Mat2D &Y_orig,
        const Row & target_values
    );

    std::vector<size_t> particle_ranking_PLS(
        const Mat2D & X_orig, const Mat2D & Y_orig,
        const Row & target_values,
        const float_type training_fraction
    );

    Row weight_predictive_prior(
        const std::vector<Parameter*> & mpars,
        const Mat2D & params
    );

    Row weight_predictive_prior(
        const std::vector<Parameter*> & mpars,
        const Mat2D & params,
        const Mat2D & prev_params,
        const Row & prev_weights,
        const Row & prev_doubled_variance
    );

    Row calculate_doubled_variance(
        const Mat2D & params
    );

}

#endif
