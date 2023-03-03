#include <limits>
#include "AbcUtil.h"
#include "AbcSmc.h"
#include "RunningStat.h"
#include "gsl/gsl_multimin.h"
#include "gsl/gsl_sf_gamma.h"

using std::string;
using std::vector;
using std::cout;
using std::cerr;
using std::endl;
using std::stringstream;
using std::ifstream;
using std::pow;
using std::pair;
using std::make_pair;

namespace ABC {

  vector<string> split(const string &s, char delim) {
      vector<string> tokens;
      stringstream ss(s);
      string token;
      while (getline(ss, token, delim)) {
          tokens.push_back(token);
      }
      return tokens;
  }


  string slurp(string filename) {
      ifstream ifs(filename.c_str());

      stringstream sstr;
      sstr << ifs.rdbuf();
      return sstr.str();
  }


  string get_nth_line(const std::string& filename, const size_t N) {
       ifstream in(filename.c_str());

       string s;
       //setting a max expected line length might improve performance
       //s.reserve(some_reasonable_max_line_length);

       //skip N lines
       for(size_t i = 0; i < N; ++i) std::getline(in, s);

       std::getline(in,s);
       return s;
  }


  Mat2D read_matrix_file(string filename, char sep) {
      cerr << "Loading " << filename << endl;
      ifstream myfile(filename.c_str());

      vector<vector<double> > M;
      if (myfile.is_open()) {
          string line;

          while ( getline(myfile,line) ) {
              //split string based on "," and store results into vector
              vector<string> fields = split(line, sep);

              vector<double>row(fields.size());
              for( unsigned int i=0; i < fields.size(); i++ ) {
                  row[i] = string2double(fields[i]);
              }
              M.push_back(row);
          }
      }

      Mat2D X( (int) M.size(), (int) M[0].size() );
      for( unsigned int i=0; i < M.size(); i++ ) {
          for( unsigned int j=0; j < M[i].size(); j++ ) {
              X(i,j)=M[i][j];
          }
      }
      return X;
  }

  Row col_stdev(const Mat2D& m, const Row& means) {
      return (m.rowwise() - means).array().square().colwise().mean().sqrt();
  }

  Row col_stdev(const Mat2D& m) {
      return col_stdev(m, m.colwise().mean());
  }

  Mat2D colwise_z_scores(const Mat2D& mat) {
      Row means, stdev;
      return colwise_z_scores(mat, means, stdev);
  }

  Mat2D colwise_z_scores(const Mat2D& mat, Row& means, Row& stdev) {
      // Standardize values by column, i.e. convert to Z-scores
      means = mat.colwise().mean();
      stdev = col_stdev( mat, means );

      // This is not technically correct, since z scores are undefined if the stdev is 0.
      // In our case, however, the resulting scores cannot be nan, or downstream calculations
      // are fouled up, so we basically set them to 0 as a stop-gap solution.  A better solution
      // would be to not pass this parameter in to the PLS regression.  This should be possible
      // but needs to be implemented with care.
      for (int c = 0; c<stdev.size(); c++) if (stdev[c] == 0) stdev[c] = 1;

      Mat2D zmat = Mat2D::Zero(mat.rows(), mat.cols());
      for (int r = 0; r<mat.rows(); r++) { zmat.row(r) = (mat.row(r) - means).cwiseQuotient(stdev); }
      return zmat;
  }

  Row z_transform_vals(const Row& vals, const Row& means, const Row& stdevs) {
    assert((vals.size() == means.size()) and (vals.size() == stdevs.size()));
    return (vals - means).array() / stdevs.array();
  }

  /*
  template <typename T> int sgn(T val) {
      return (T(0) < val) - (val < T(0));
  }

  // -- OR --
  template <typename T> inline constexpr
  int signum(T x, std::false_type is_signed) {
      return T(0) < x;
  }

  template <typename T> inline constexpr
  int signum(T x, std::true_type is_signed) {
      return (T(0) < x) - (x < T(0));
  }

  template <typename T> inline constexpr
  int signum(T x) {
      return signum(x, std::is_signed<T>());
  }*/

  /*
  double normal_pdf(double x, double mu, double var) {
      long double PI = 3.1415926535897932384;
      return exp(-pow(x-mu,2) / (2.0*var)) / sqrt(2*PI*var);
  }

  double normal_cdf(double x, double mu, double var) {
      x = (x-mu)/sqrt(var);
      // Abramowitz & Stegun (1964) approximation
      long double b0 = 0.2316419;
      double b1 = 0.319381530;
      double b2 = -0.356563782;
      double b3 = 1.781477937;
      double b4 = -1.821255978;
      double b5 = 1.330274429;
      if (x >= 0.0) {
          long double t = 1.0/(1.0+b0*x);
          return 1.0 - normal_pdf(x, 0, 1)*(b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
      } else {
          long double t = 1.0/(1.0-b0*x);
          return normal_pdf(x, 0, 1)*(b1*t + b2*pow(t,2) + b3*pow(t,3) + b4*pow(t,4) + b5*pow(t,5));
      }
  }*/

  float_type mean(const Col data) {
      assert( data.size() > 0 );
      return data.sum() / data.size();
  }

  float_type median(const Col data) {
      assert(data.size() > 0);
      // copy & sort data
      vector<float_type> vdata(data.data(), data.data()+data.size());
      const int n = vdata.size();
      sort(vdata.begin(), vdata.end());

      float_type median;

      if (n % 2 == 0) {
          median = (vdata[n / 2 - 1] + vdata[n / 2]) / 2;
      } else {
          median = vdata[n / 2];
      }

      return median;
  }

  float_type quantile(const Col data, double q) {
      return ::quantile(as_vector(data), q);
  }

  float_type variance(const Col data, float_type _mean) {
      if (data.size() < 2) {
          cerr << "WARNING: Variance called with " << data.size() << " data values. Returning 0." << endl;
          return 0;
      } else {
          return (data.array() - _mean).square().sum() / (data.size() - 1);
      }
  }

  float_type max(const Col data) {
      assert(data.size() > 0);
      auto colv = data.reshaped();
      return *std::max_element(colv.begin(), colv.end());
  }

  float_type skewness(const Col data) {
      float_type _x = mean(data);
      float_type _v = variance(data, _x);
      if (_v == 0) return 0; // Avoids nans.  If variance is 0, define skewness as 0
      return ((data.array() - _x).pow(3).sum() / data.size() ) / pow(_v, 1.5);
  }

  int _mc_pos(const float_type v, const float_type m) {
      enum Position {ABOVE, BELOW, AT};
      Position p;
      if (v > m) { p = ABOVE; }
      else if (v < m) { p = BELOW; }
      else { p = AT; }
      return (int) p;
  }

  float_type median_crossings(const Col data) {
      if (data.size() < 2) {
          return 0;
      } else {
          return median_crossings(data, median(data));
      }
  }

  // Calculates the number of times the data series crosses the median
  float_type median_crossings(const Col data, const float_type m) {
      int mc = 0;
      if (data.size() < 2) return mc;

      enum Position {ABOVE, BELOW, AT};
      // current and next are like cursors that trace the data series
      Position current, next;
      current = (Position) _mc_pos(data[0], m);
      if (current == AT) mc++; // increment if we're starting at the median
      for (size_t i = 1; i < static_cast<size_t>(data.size()); ++i) {
          next = (Position) _mc_pos(data[i], m);
          if (next != current and current != AT) mc++; // just crossed or touched the median
          current = next;
      }
      return ((float_type) mc)/(data.size()-1);
  }

  float optimize_box_cox (const Col data, float lambda_min, float lambda_max, float step) {
      float best_lambda = lambda_min;
      float min_skew = std::numeric_limits<float>::infinity();
      float skew;
      for (float lambda = lambda_min; lambda <= lambda_max; lambda += step) {
          if (lambda == 0) {
              skew = skewness( data.array().log() );
          } else {
              skew = skewness( (data.array().pow(lambda) - 1).matrix() / lambda );
          }
          if (abs(skew) < abs(min_skew)) {
              min_skew = skew;
              best_lambda = lambda;
          }
      }
      return best_lambda;
  }

  float optimize_box_cox (const Col data) {
      return optimize_box_cox(data, -5, 5, 0.1);
  }

  template<typename Iterable>  
  int gsl_rng_nonuniform_int(const Iterable & weights, const gsl_rng* rng) {
      // assert: weights is a CDF - i.e. sorted, first element > 0, the last element is 1.0
      double r = gsl_rng_uniform(rng);
      auto pos = std::upper_bound(weights.begin(), weights.end(), r);
      if (pos != weights.end()) {
          return std::distance(weights.begin(), pos);
      } else {
        std::cerr << "ERROR: Weights may not be a CDF." << std::endl << "\tweights.end() = " << *(weights.end()) << std::endl;
        exit(100);
      }
  }

  Row rand_trunc_mv_normal(const std::vector<Parameter*> _model_pars, const Row & parrow, const gsl_matrix* L, const gsl_rng* rng) {
      gsl_vector* mu = gsl_vector_alloc(parrow.size());
      for (size_t i = 0; i < parrow.size(); i++) { gsl_vector_set(mu, i, parrow[i]); }
      const size_t npar = _model_pars.size();
      Row par_values = Row::Zero(npar);
      gsl_vector* result = gsl_vector_alloc(npar);
      bool success = false;
      while (not success) {
          success = true;
          gsl_ran_multivariate_gaussian(rng, mu, L, result);
          for (size_t j = 0; (j < npar) and success; j++) {
              par_values[j] = gsl_vector_get(result, j);
              if (_model_pars[j]->get_numeric_type() == INT) par_values(j) = (double) ((int) (par_values(j) + 0.5));
              if (par_values[j] < _model_pars[j]->get_prior_min() or par_values[j] > _model_pars[j]->get_prior_max()) success = false;
          }
      }
      gsl_vector_free(result);
      gsl_vector_free(mu);
      return par_values;
  }

  double rand_trunc_normal(double mu, double sigma_squared, double min, double max, const gsl_rng* rng) {
      assert(min < max);
      double sigma = sqrt(sigma_squared);
      // Don't like this, but it will work
      // as long as min and max are reasonable
      // (relative to the pdf)
      while (1) {
          double dev = gsl_ran_gaussian(rng, sigma) + mu;
          if (dev >= min and dev <= max) {
              return dev;
          }
      }
  }

  LinearFit* lin_reg(const std::vector<double> &x, const std::vector<double> &y) {
      assert( x.size() == y.size() );
      LinearFit* fit = new LinearFit();
      const int n = x.size();
      double sumx = 0.0;                        /* sum of x                      */
      double sumx2 = 0.0;                       /* sum of x**2                   */
      double sumxy = 0.0;                       /* sum of x * y                  */
      double sumy = 0.0;                        /* sum of y                      */
      double sumy2 = 0.0;                       /* sum of y**2                   */

      for (int i=0; i<n; i++)   {
          sumx  += x[i];
          sumx2 += pow(x[i],2);
          sumxy += x[i] * y[i];
          sumy  += y[i];
          sumy2 += pow(y[i],2);
      }

      double denom = n * sumx2 - pow(sumx,2);
      if (denom == 0) {
          // singular matrix. can't solve the problem.
          fit->m   = 0;
          fit->b   = 0;
          fit->rsq = 0;
          return fit;
      }

      fit->m = (n * sumxy  -  sumx * sumy) / denom;
      fit->b = (sumy * sumx2  -  sumx * sumxy) / denom;
      // compute correlation coeff
      fit->rsq = pow((sumxy - sumx * sumy / n) / sqrt((sumx2 - pow(sumx,2)/n) * (sumy2 - pow(sumy,2)/n)),2);

      return fit;
  }

  struct LogisticDatum{
      LogisticDatum() : time(0.0), successes(0), attempts(0) {};
      LogisticDatum(double t, int s, int a) : time(t), successes(s), attempts(a) {};
      double time;
      unsigned int successes;
      unsigned int attempts;
  };

  double _logistic_likelihood(const gsl_vector *v, void* params) {
      long double b0, b1;
      vector<LogisticDatum*>* data = (vector<LogisticDatum*>*) params;

      b0 = gsl_vector_get(v, 0);
      b1 = gsl_vector_get(v, 1);

      double total = 0;

      // for each year, calculate probability of observing that # of severe cases; sum log probs
      for (auto p: *data) {
          const double time = p->time;
          unsigned int suc = p->successes;
          unsigned int all = p->attempts;
          const long double prob = 1.0 / (1.0 + exp(-(b0 + b1*time)));
          const long double lnchoosek = gsl_sf_lnchoose(all, suc);
          const long double ldensity = lnchoosek + suc*log(prob) + (all-suc)*log(1.0-prob);
          total += ldensity;
      }

      // total should always be <= 0.  0 would indicate a probability of 1.0
      // give it a bad score (e.g., INT_MIN) if we got garbage
      total = (std::isinf(total) or std::isnan(total)) ? INT_MIN : total;
      // gsl is expecting a minimization; we're maximizing a probability
      return -total;
  }

  LogisticFit* logistic_reg(vector<LogisticDatum*> data) {
      LogisticFit* fit = new LogisticFit();

      // Initial values for fitted parameters (beta0 and beta1)
      // 0.0, 0.0 describes a flat line with 50% probability for all x
      gsl_vector *betas;
      betas = gsl_vector_alloc (2);
      gsl_vector_set (betas, 0, 0.0); // beta 0
      gsl_vector_set (betas, 1, 0.0); // beta 1

      // Set initial step sizes to 0.01
      gsl_vector *ss;
      ss = gsl_vector_alloc (2);
      gsl_vector_set_all (ss, 0.01);

      // Initialize method and iterate
      gsl_multimin_function multimin_func;
      multimin_func.n = 2;
      multimin_func.f = _logistic_likelihood;
      multimin_func.params = &data;

      const gsl_multimin_fminimizer_type *T = gsl_multimin_fminimizer_nmsimplex2;
      gsl_multimin_fminimizer* s = gsl_multimin_fminimizer_alloc (T, 2);
      gsl_multimin_fminimizer_set (s, &multimin_func, betas, ss);

      int status = GSL_CONTINUE;
      double size = 1.0;
      for (int iter = 0; iter<10000; ++iter) {
          fit->iterations = iter;
          status = gsl_multimin_fminimizer_iterate(s);

          if (status) break;

          size = gsl_multimin_fminimizer_size (s);
          status = gsl_multimin_test_size (size, 1e-4);

          //if (status == GSL_SUCCESS) printf ("converged to minimum at\n");
          //printf ("%5d %10.3e %10.3e f() = %7.3f size = %.5f\n", iter, gsl_vector_get (s->x, 0), gsl_vector_get (s->x, 1), s->fval, size);

          if (status != GSL_CONTINUE) break;
      }

      for (auto datum: data) delete datum;

      fit->status = status; // Status should be checked, and should equal GSL_SUCCESS before using beta0 and beta1
      if (fit->status != GSL_SUCCESS) cerr << "WARNING: Logistic regression was unsuccessful (did not converge)\n";
      fit->beta0 = gsl_vector_get (s->x, 0);
      fit->beta1 = gsl_vector_get (s->x, 1);
      fit->simplex_size = size;

      gsl_vector_free(betas);
      gsl_vector_free(ss);
      gsl_multimin_fminimizer_free(s);

      return fit;
  }

  LogisticFit* logistic_reg(const vector<double> &x, const vector<int> &s, const vector<int> &a) {
      assert( x.size() == s.size() );
      assert( x.size() == a.size() );
      vector<LogisticDatum*> data;
      for (size_t i = 0; i < x.size(); ++i) {
          LogisticDatum* datum = new LogisticDatum(x[i], (unsigned) s[i], (unsigned) a[i]);
          data.push_back(datum);
      }
      return logistic_reg(data);
  }

  LogisticFit* logistic_reg(const std::vector<double> &x, const std::vector< pair<int,int> > &y) {
      assert( x.size() == y.size() );
      vector<LogisticDatum*> data;
      for (size_t i = 0; i < x.size(); ++i) {
          LogisticDatum* datum = new LogisticDatum(x[i], (unsigned) y[i].first, (unsigned) y[i].second);
          data.push_back(datum);
      }
      return logistic_reg(data);
  }

}

double ABC::calculate_nrmse(
    const Mat2D & posterior_mets, // rows are posterior samples, columns are metrics
    const Row & observed
) {
    assert(posterior_mets.cols() == observed.size());

    // get the simulated mean of each metric
    Row sim = posterior_mets.colwise().mean();
    // take the average of the observed and simulated values as the "reference"
    // tjh TODO: where does the fabs come from?
    Row expected = (observed.array().abs() + sim.array().abs()) / 2.0;
    // where sim == obs, set expected == 1 (precludes divide by zero, and otherwise doesn't matter)
    for (size_t i = 0; i < expected.size(); ++i) { if (sim[i] == observed[i]) expected[i] = 1; }

    // delta = (sim - obs) / expected => each squared => collapse to mean => sqrt
    double res = mean(((sim - observed).array() / expected.array()).square());

    return sqrt(res);

}

gsl_matrix* ABC::setup_mvn_sampler(
    const Mat2D & particle_parameters // sliced to predictive prior
) {
    // ALLOCATE DATA STRUCTURES
    // variance-covariance matrix calculated from pred prior values
    // NB: always allocate this small matrix, so that we don't have to check whether it's safe to free later
    gsl_matrix* sigma_hat = gsl_matrix_alloc(particle_parameters.cols(), particle_parameters.cols());

    // container for predictive prior aka posterior from last set
    gsl_matrix* posterior_par_vals = gsl_matrix_alloc(particle_parameters.rows(), particle_parameters.cols());

    // INITIALIZE DATA STRUCTURES
    for (size_t i = 0; i < particle_parameters.rows(); ++i) {
        auto parrow = particle_parameters.row(i);
        for (size_t j = 0; j < particle_parameters.cols(); ++j) {
            // copy values from pred prior into a gsl matrix
            gsl_matrix_set(posterior_par_vals, i, j, parrow[j]);
        }
    }
    // calculate maximum likelihood estimate of variance-covariance matrix sigma_hat
    gsl_ran_multivariate_gaussian_vcov(posterior_par_vals, sigma_hat);

    for (size_t j = 0; j < particle_parameters.cols(); j++) {
        // sampling is done using a kernel with a broader kernel than found in pred prior values
        const double doubled_variance = 2 * gsl_matrix_get(sigma_hat, j, j);
        gsl_matrix_set(sigma_hat, j, j, doubled_variance);
    }

    // not a nice interface, gsl.  sigma_hat is converted in place from a variance-covariance matrix
    // to the same values in the upper triangle, and the diagonal and lower triangle equal to L,
    // the Cholesky decomposition
    gsl_linalg_cholesky_decomp1(sigma_hat);

    gsl_matrix_free(posterior_par_vals);

    return sigma_hat;
}

void ABC::normalize_weights( std::vector<float_type> & weights ) {
    double total = std::accumulate(weights.begin(), weights.end(), 0.0);
    for (size_t i = 0; i < weights.size(); i++) {
        weights[i] /= total;
    }
}

Col ABC::CDF_weights(const Col & weights) {
    auto total = weights.sum();
    Col cdf = Col::Zero(weights.size());
    std::partial_sum(weights.begin(), weights.end(), cdf.begin());
    return cdf.array() / total;
}

Col ABC::weight_predictive_prior(
    const Mat2D & current_parameters,
    const Mat2D & prev_parameters,
    const std::vector<float_type> & prev_weights,
    const std::vector<Parameter*> & model_pars,
    const Row & prev_doubled_variances
) {

    assert(current_parameters.cols() == prev_parameters.cols());
    assert(current_parameters.cols() == model_pars.size());
    assert(prev_parameters.rows() == prev_weights.size());
    assert(prev_doubled_variances.size() == model_pars.size());

    Col weights = Col::Zero(current_parameters.rows());

    for (size_t i = 0; i < current_parameters.rows(); i++) { // for each particle in the predictive prior
    
        const Row parrow = current_parameters.row(i);
        double numerator = 1;
        double denominator = 0.0;

        // likelihood of parrow, given priors
        for (size_t j = 0; j < model_pars.size(); j++) {
            Parameter* par = model_pars[j];
            const double par_value = parrow[j];
            if (par->get_prior_type() == NORMAL) {
                numerator *= gsl_ran_gaussian_pdf(par_value - par->get_prior_mean(), par->get_prior_stdev());
            } else if (par->get_prior_type() == UNIFORM) {
                // The RHS here will be 1 under normal circumstances.  If the prior has been revised during a fit,
                // this should throw out values outside of the prior's range
                numerator *= (int) (par_value >= par->get_prior_min() and par_value <= par->get_prior_max());
            }
        }

        // likelihood of parrow, given all previous rounds particles
        for (size_t k = 0; k < prev_parameters.rows(); k++) { // for each particle in the previous predictive prior
            double running_product = prev_weights[k]; // how likely was this (previous) particle
            auto prev_parrow = prev_parameters.row(k);
            for (size_t j = 0; j < model_pars.size(); j++) {
                double par_value = parrow[j];
                double old_par_value = prev_parrow[j];
                double old_doubled_variance = prev_doubled_variances[j];

                // This conditional handles the (often improbable) case where a parameter has completely converged.
                // It allows ABC to continue exploring other parameters, rather than causing the math
                // to fall apart because the density at the converged value is infinite.
                if (old_doubled_variance != 0 or par_value != old_par_value) {
                    running_product *= gsl_ran_gaussian_pdf(par_value-old_par_value, sqrt(old_doubled_variance) );
                }
            }
            denominator += running_product;
        }

        weights[i] = numerator / denominator;
    }
    weights.normalize();
    return weights;
    
}

Row ABC::sample_predictive_prior(
    const gsl_rng* RNG,
    const Col & weights,
    const Mat2D & particle_parameters // sliced to predictive prior
) {
    assert(weights.size() == particle_parameters.rows());
    // weighted-randomly draw a particle from predictive prior to use as mean of new value
    return particle_parameters.row(gsl_rng_nonuniform_int(weights, RNG));
}

// sample_predictive_prior(RNG, weights, particle_parameters)

Row ABC::noise_mv_parameters(
    const gsl_rng* RNG,
    const Row & pred_prior_mu,
    const std::vector<Parameter*> & model_pars,
    const gsl_matrix* L
) {
    // return a new particle from a multivariate normal distribution with mean parrow and covariance L
    return rand_trunc_mv_normal(model_pars, pred_prior_mu, L, RNG);
}

Row ABC::noise_parameters(
    const gsl_rng* RNG,
    const Row & pred_prior_mu,
    const std::vector<Parameter*> & model_pars,
    const std::vector<float_type> & doubled_variances 
) {
    Row par_values = Row::Zero(model_pars.size());

    for (size_t j = 0; j < model_pars.size(); j++) {
        double par_value = pred_prior_mu[j];
        const Parameter* parameter = model_pars[j];
        double doubled_variance = doubled_variances[j];
        double par_min = parameter->get_prior_min();
        double par_max = parameter->get_prior_max();
        par_values(j) = rand_trunc_normal( par_value, doubled_variance, par_min, par_max, RNG );

        if (parameter->get_numeric_type() == INT) {
            par_values(j) = (double) ((int) (par_values(j) + 0.5));
        }
    }

    return par_values;
}

Row ABC::calculate_doubled_variances(const Mat2D & particle_parameters) {

    std::vector<RunningStat> stats(particle_parameters.cols());
    Row dvs = Row::Zero(particle_parameters.cols());

    for (size_t j = 0; j < particle_parameters.cols(); j++) { stats[j].Push( particle_parameters.col(j) ); }

    for (size_t j = 0; j < particle_parameters.cols(); j++) {
        dvs[j] =  2 * stats[j].Variance();
    }

    return dvs;
}