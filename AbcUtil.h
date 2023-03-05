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

#include "pls.h"

using namespace std;

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

  vector<string> split(const std::string& s, char c);

  inline int string2int(const std::string& s){ std::istringstream i(s); int x = 0; i >> x; return x; }
  inline double string2double(const std::string& s){ std::istringstream i(s); double x = 0; i >> x; return x; }
  inline float_type string2float_type(const std::string& s){ std::istringstream i(s); float_type x = 0; i >> x; return x; }

  std::string slurp(std::string filename);
  std::string get_nth_line(const std::string& filename, int N);

  inline float_type logit(const float_type p) { assert(p <= 1.0); assert (p >= 0.0); return log( p / (1.0 - p) ); }
  inline float_type logistic(const float_type l) { return 1.0 / (1.0 + exp(-l)); }

  //int _sgn(float_type val) { return (0 < val) - (val < 0); }

  Mat2D read_matrix_file(std::string filename, char sep);

  Row colwise_stdev(const Mat2D & mat, const Row & means);

  Mat2D colwise_z_scores( const Mat2D & mat, const Row & mean, const Row & stdev);
  Mat2D colwise_z_scores( const Mat2D & mat );
  Row z_scores(const Row & metobs, const Row & means, const Row & stdevs);

  float_type mean(const Col data);

  float_type median(const Col data);

  float_type quantile(const Col data, double q);

  float_type variance(const Col data, float_type _mean);

  float_type max(const Col data);

  float_type skewness(const Col data);

  float_type median_crossings(const Col data, const float_type _median);

  float_type median_crossings(const Col data);

  float optimize_box_cox (const Col data, float lambda_min, float lambda_max, float step);

  float optimize_box_cox (const Col data);

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

  template<typename Iterable>
  int gsl_rng_nonuniform_int(const Iterable & weights, const gsl_rng* rng);

  double rand_trunc_normal(double mu, double sigma_squared, double min, double max, const gsl_rng* rng);
  Row rand_trunc_mv_normal(const vector<Parameter*> _model_pars, gsl_vector* mu, gsl_matrix* L, const gsl_rng* rng);

  LinearFit* lin_reg(const std::vector<double> &x, const std::vector<double> &y);

  LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector< pair<int,int> > &y);

  LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector<int> &successes, const std::vector<int> &attempts);

  inline vector<float_type> as_vector(const Row data) {
      vector<float_type> vec(data.size());
      for (size_t i = 0; i < static_cast<size_t>(data.size()); i++) vec[i] = data[i];
      return vec;
  }

  inline Row as_row(const vector<float_type> data) {
      Row row(data.size());
      for (size_t i = 0; i < data.size(); i++) row[i] = data[i];
      return row;
  }

  double calculate_nrmse(
    const Mat2D & posterior_mets,
    const Row & observed
  );

  Row sample_posterior(
    const Col weights,
    const Mat2D posterior,
    const gsl_rng* RNG
  );

  template<typename RandomAccessible>
  gsl_vector* to_gsl_v(const RandomAccessible & from);

  gsl_matrix* to_gsl_m(const Mat2D & from);

}

#endif
