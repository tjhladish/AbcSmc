#ifndef ABCUTIL_H
#define ABCUTIL_H

#include <iostream>
#include <sstream>
#include <vector>
#include <assert.h>
#include <iomanip>
#include <fstream>
#include <Eigen/Eigenvalues>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_linalg.h>
#include "ranker.h"
#include <math.h>

using namespace std;

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

  #ifdef USING_MPI
  #include <mpi.h>
  struct MPI_par {
      MPI_Comm comm;
      MPI_Info info;
      int mpi_size, mpi_rank;
  };
  #else
  struct MPI_par {
      const static int mpi_size = 1;
      const static int mpi_rank = 0;
  };
  #endif

  //using namespace std;
  using namespace Eigen;

  #ifdef MPREAL_SUPPORT
  #include "mpreal.h"
  #include <unsupported/Eigen/MPRealSupport>
      using namespace mpfr;
      typedef mpreal float_type;
  #else
      typedef double float_type;
  #endif

  typedef Matrix<float_type,Dynamic,Dynamic> Mat2D;
  typedef Matrix<float_type, Dynamic, 1>  Col;
  typedef Matrix<float_type, 1, Dynamic>  Row;
  typedef Matrix<int, 1, Dynamic>  Rowi;
  typedef Matrix<std::complex<float_type>,Dynamic,Dynamic> Mat2Dc;
  typedef Matrix<std::complex<float_type>, Dynamic, 1>  Colc;

  vector<string> split(const std::string& s, char c);

  inline int string2int(const std::string& s){ std::istringstream i(s); int x = 0; i >> x; return x; }
  inline double string2double(const std::string& s){ std::istringstream i(s); double x = 0; i >> x; return x; }
  inline float_type string2float_type(const std::string& s){ std::istringstream i(s); float_type x = 0; i >> x; return x; }

  std::string slurp(std::string filename);
  std::string get_nth_line(const std::string& filename, int N);

  inline float_type logit(const float_type p) { assert(p <= 1.0); assert (p >= 0.0); return log( p / (1.0 - p) ); }
  inline float_type logistic(const float_type l) { return 1.0 / (1.0 + exp(-l)); }

  inline Row col_means( Mat2D mat ) { return mat.colwise().sum() / mat.rows(); }

  //int _sgn(float_type val) { return (0 < val) - (val < 0); }

  Mat2D read_matrix_file(std::string filename, char sep);

  Row col_stdev( Mat2D mat, Row means );

  float_type dominant_eigenvalue( EigenSolver<Mat2Dc> es );

  Colc dominant_eigenvector( EigenSolver<Mat2D> es );

  Mat2D colwise_z_scores( const Mat2D& mat );
  Mat2D colwise_z_scores( const Mat2D& mat, Row& means, Row& stdev );

  std::vector<int> ordered(Col const& values);

  float_type wilcoxon(const Col err_1, const Col err_2);

  float_type normalcdf(float_type z);

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

  std::string exec(std::string cmd);

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


  template <typename T>
  inline std::string toString (const T& t) {
      std::stringstream ss;
      ss << t;
      return ss.str();
  }

  inline double uniform_pdf(double a, double b) { return 1.0 / fabs(b-a); }

  int gsl_rng_nonuniform_int(std::vector<double>& weights, const gsl_rng* rng);

  double rand_trunc_normal(double mu, double sigma_squared, double min, double max, const gsl_rng* rng);

  LinearFit* lin_reg(const std::vector<double> &x, const std::vector<double> &y);

  LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector< pair<int,int> > &y);

  LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector<int> &successes, const std::vector<int> &attempts);

  inline vector<float_type> as_vector(const Row data) {
      vector<float_type> vec(data.size());
      for (int i = 0; i < data.size(); i++) vec[i] = data[i];
      return vec;
  }

  inline Row as_row(const vector<float_type> data) {
      Row row(data.size());
      for (size_t i = 0; i < data.size(); i++) row[i] = data[i];
      return row;
  }

}

#endif
