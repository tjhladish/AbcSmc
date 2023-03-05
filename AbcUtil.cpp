#include <limits>
#include "AbcUtil.h"
#include "AbcSmc.h"
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


  Row colwise_stdev(const Mat2D & mat, const Row & means ) {
      const float_type N = mat.rows();
      if ( N < 2 ) return Row::Zero(mat.cols());
      // N-1 for unbiased sample variance
      return ((mat.rowwise() - means).array().square().colwise().sum()/(N-1)).sqrt();
  }

  Row z_scores(
    const Row & vals, const Row & means, const Row & stdev
  ) {
    return (vals - means).cwiseQuotient(stdev);
  }

Mat2D colwise_z_scores(const Mat2D & mat, const Row & mean, const Row & stdev) {

    Row local_sd = stdev;
    // sd == 0 => implies all values the same => x_i - mean == 0
    // Technically: z scores are undefined if the stdev is 0 => this should yield nan.
    // However, scores == nan => borks downstream calculations
    // This makes the scores == 0 instead.
    // TODO: change algorithm to not pass this parameter to PLS.
    for (int c = 0; c < local_sd.size(); c++) if (local_sd[c] == 0) local_sd[c] = 1;
    Mat2D mmeans = mat.rowwise() - mean;

    Mat2D zs = mmeans.array().rowwise() / stdev.array();
    return zs;
};

Mat2D colwise_z_scores(const Mat2D & mat) {
    Row means = mat.colwise().mean();
    Row stdev = colwise_stdev(mat, means);
    return colwise_z_scores(mat, means, stdev);
};


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
  int gsl_rng_nonuniform_int(const gsl_rng* rng, const Iterable & weights) {
      // Weights must sum to 1!!
      double running_sum = 0;
      double r = gsl_rng_uniform(rng);
      for (size_t i = 0; i < weights.size(); i++) {
          running_sum += weights[i];
          if (r <= running_sum) return i;
      }
      cerr << "ERROR: Weights may not be normalized\n\t Weights summed to: " << running_sum << endl;
      exit(100);
  }

  Row rand_trunc_mv_normal(const vector<Parameter*> _model_pars, gsl_vector* mu, gsl_matrix* L, const gsl_rng* rng) {
      const size_t npar = _model_pars.size();
      Row par_values = Row::Zero(npar);
      gsl_vector* result = gsl_vector_alloc(npar);
      bool success = false;
      while (not success) {
          success = true;
          gsl_ran_multivariate_gaussian(rng, mu, L, result);
          for (size_t j = 0; j < npar; j++) {
              par_values[j] = gsl_vector_get(result, j);
              if (_model_pars[j]->get_numeric_type() == INT) par_values(j) = (double) ((int) (par_values(j) + 0.5));
              if (par_values[j] < _model_pars[j]->get_prior_min() or par_values[j] > _model_pars[j]->get_prior_max()) success = false;
          }
      }
      gsl_vector_free(result);
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

  Col euclidean(const Mat2D & sims, const Row & ref) {
    // for each row in sims, subtract the reference row
    // norm = sqrt(sum_i xi^2)
    return (sims.rowwise() - ref).rowwise().norm();
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

template<typename RandomAccessible>
gsl_vector* ABC::to_gsl_v(const RandomAccessible & from) {
    gsl_vector* res = gsl_vector_alloc(from.size());
    for (size_t i = 0; i < from.size(); i++) {
        gsl_vector_set(res, i, from[i]);
    }
    return res;
};

gsl_matrix* ABC::to_gsl_m(const Mat2D & from){
    gsl_matrix* res = gsl_matrix_alloc(from.rows(), from.cols());
    for (size_t i = 0; i < from.rows(); i++) {
        for (size_t j = 0; j < from.cols(); j++) {
            gsl_matrix_set(res, i, j, from(i,j));
        }
    }
    return res;
};

Row ABC::sample_posterior(
    const gsl_rng* RNG,
    const Col weights,
    const Mat2D posterior
) {
    return posterior(gsl_rng_nonuniform_int(RNG, weights), Eigen::placeholders::all);
};

Row ABC::sample_predictive_priors(
    const gsl_rng* RNG,
    const Col & weights, const Mat2D & parameter_prior,
    const std::vector<Parameter*> & pars,
    const Row & doubled_variance
) {
    const Row par = ABC::sample_posterior(RNG, weights, parameter_prior);
    Row new_par = Row::Zero(par.cols());
    for (size_t parIdx = 0; parIdx < par.cols(); parIdx++) {
        const Parameter* parameter = pars[parIdx];
        double par_min = parameter->get_prior_min();
        double par_max = parameter->get_prior_max();
        new_par(parIdx) = ABC::rand_trunc_normal(par[parIdx], doubled_variance[parIdx], par_min, par_max, RNG );

        if (parameter->get_numeric_type() == INT) {
            new_par(parIdx) = (double) ((int) (new_par(parIdx) + 0.5));
        }
    }
    return new_par;
}