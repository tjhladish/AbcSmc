#include <limits>
#include <compare>
#include <AbcSmc/AbcUtil.h>
#include <AbcSmc/AbcSmc.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_sf_gamma.h>
#include <AbcSmc/RunningStat.h>

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
using std::partial_ordering;

using namespace PLS;

namespace ABC {

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

    float_type median(const Col & data) {
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

    float_type quantile(const Col & data, double q) {
        return ::quantile(as_vector(data), q);
    }

    float_type variance(const Col & data, const float_type _mean) {
        if (data.size() < 2) {
            cerr << "WARNING: Variance called with " << data.size() << " data values. Returning 0." << endl;
            return 0;
        } else {
            return (data.array() - _mean).square().sum() / (data.size() - 1);
        }
    }

    float_type max(const Col & data) {
        assert(data.size() > 0);
        return data.maxCoeff();
    }

    float_type skewness(const Col & data) {
        float_type _x = data.mean();
        float_type _v = variance(data, _x);
        if (_v == 0) return 0; // Avoids nans.  If variance is 0, define skewness as 0
        return ((data.array() - _x).pow(3).sum() / data.size() ) / pow(_v, 1.5);
    }

    // Calculates the number of times the data series crosses the median
    float_type median_crossings(const Col & data, const float_type m) {
        if (data.size() < 2) return 0.0;
        
        size_t mc = 0;
 
        // current and next are like cursors that trace the data series
        // this uses the <=> operator: (a <=> b) returns equivalent if a == b, less if a < b, and greater if a > b
        partial_ordering current = data[0] <=> m;
        if (current == partial_ordering::equivalent) mc++; // increment if we're starting at the median
        for (auto pos : data) { // for each data point (n.b. the first pass here is a no-op: next == current)
            partial_ordering next = pos <=> m;
            if ((next != current) and (current != partial_ordering::equivalent)) mc++; // just crossed or touched the median
            current = next;
        }
        return (static_cast<float_type>(mc)/(data.size() - 1.0));
    }

    float_type median_crossings(const Col & data) {
        if (data.size() < 2) {
            return 0;
        } else {
            return median_crossings(data, median(data));
        }
    }

    float_type optimize_box_cox(const Col & data, const float lambda_min, const float lambda_max, const float step) {
        float_type best_lambda = lambda_min;
        float_type min_skew = std::numeric_limits<float_type>::infinity();
        float_type skew;
        for (float_type lambda = lambda_min; lambda <= lambda_max; lambda += step) {
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

    float_type optimize_box_cox(const Col & data) {
        return optimize_box_cox(data, -5, 5, 0.1);
    }

    std::vector<size_t> gsl_rng_nonuniform_int(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights
    ) {
        gsl_ran_discrete_t* gslweights = gsl_ran_discrete_preproc(weights.rows(), weights.data());
        std::vector<size_t> res(num_samples);
        for (size_t i = 0; i < res.size(); i++) { res[i] = gsl_ran_discrete(RNG, gslweights); }
        gsl_ran_discrete_free(gslweights);
        return res;
    }

    Row gsl_ran_trunc_mv_normal(
        const gsl_rng* RNG,
        const vector<const Parameter*> _model_pars,
        const Row & mu, const gsl_matrix* L
    ) {
        const size_t npar = _model_pars.size();
        Row par_values = Row::Zero(npar);
        gsl_vector* result = gsl_vector_alloc(npar);
        bool success = false;
        gsl_vector * gslmu = to_gsl_v(mu);
        while (not success) {
            success = true;
            gsl_ran_multivariate_gaussian(RNG, gslmu, L, result);
            for (size_t parIdx = 0; success and (parIdx < npar); parIdx++) {
                par_values[parIdx] = _model_pars[parIdx]->recast(gsl_vector_get(result, parIdx));
                success = _model_pars[parIdx]->valid(par_values[parIdx]);
            }
        }
        gsl_vector_free(gslmu);
        gsl_vector_free(result);
        return par_values;
    }

    Row gsl_ran_trunc_normal(
        const gsl_rng* RNG,
        const std::vector<const Parameter*> _model_pars,
        const Row & mu, const Row & sigma_squared
    ) {
        Row sigma = sigma_squared.array().sqrt();
        Row res = Row::Zero(sigma.cols());
        for (size_t parIdx = 0; parIdx < sigma.cols(); parIdx++) {
            // this manages normal noise, including finite retries; will fail to the *parameter* mean if it can't find a valid value
            res[parIdx] = _model_pars[parIdx]->noise(RNG, mu[parIdx], sigma[parIdx]);
            
        }
        return res;      
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

    float_type calculate_nrmse(
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
        float_type res = ((sim - observed).array() / expected.array()).square().mean();

        return sqrt(res);

    }

    template<typename RandomAccessible>
    gsl_vector* to_gsl_v(const RandomAccessible & from) {
        gsl_vector* res = gsl_vector_alloc(from.size());
        for (size_t i = 0; i < from.size(); i++) {
            gsl_vector_set(res, i, from[i]);
        }
        return res;
    }

    gsl_matrix* to_gsl_m(const Mat2D & from){
        gsl_matrix* res = gsl_matrix_alloc(from.rows(), from.cols());
        for (size_t i = 0; i < from.rows(); i++) {
            for (size_t j = 0; j < from.cols(); j++) {
                gsl_matrix_set(res, i, j, from(i,j));
            }
        }
        return res;
    }

    Mat2D sample_posterior(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights,
        const Mat2D & posterior
    ) {
        return posterior( // from the posterior
            gsl_rng_nonuniform_int(RNG, num_samples, weights), // get a weighted-random draw of rows
            Eigen::placeholders::all
        );
    };

    Mat2D sample_predictive_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights, const Mat2D & parameter_prior,
        const std::vector<const Parameter*> & pars,
        const Row & doubled_variance
    ) {
        const Mat2D sampled_pars = sample_posterior(RNG, num_samples, weights, parameter_prior);
        Mat2D noised_pars = Mat2D(sampled_pars.rows(), sampled_pars.cols());
        for (size_t sampIdx = 0; sampIdx < noised_pars.rows(); sampIdx++) {
            noised_pars.row(sampIdx) = gsl_ran_trunc_normal(RNG, pars, sampled_pars.row(sampIdx), doubled_variance);
        }
        return noised_pars;
    };

    Mat2D sample_mvn_predictive_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Col & weights, const Mat2D & parameter_prior,
        const std::vector<const Parameter*> & pars,
        const gsl_matrix* L
    ) {
        // SELECT PARTICLE FROM PRED PRIOR TO USE AS EXPECTED VALUE OF NEW SAMPLE
        const Mat2D sampled_pars = sample_posterior(RNG, num_samples, weights, parameter_prior);
        Mat2D noised_pars = Mat2D(sampled_pars.rows(), sampled_pars.cols());
        for (size_t sampIdx = 0; sampIdx < noised_pars.rows(); sampIdx++) {
            noised_pars.row(sampIdx) = gsl_ran_trunc_mv_normal(RNG, pars, sampled_pars.row(sampIdx), L);
        }
        return noised_pars;
    };

    // ranks predictors (X_orig) according to z-score distance
    // from the observed values (target_values)
    std::vector<size_t> particle_ranking_simple (
        const Mat2D & X_orig, const Mat2D & /* Y_orig */,
        const Row & target_values
    ) {
        Row X_sim_means = X_orig.colwise().mean();
        Row X_sim_stdev = colwise_stdev(X_orig, X_sim_means);
        Row obs_met = z_scores(target_values, X_sim_means, X_sim_stdev);

        Mat2D X = colwise_z_scores(X_orig, X_sim_means, X_sim_stdev);
    //    Mat2D Y = colwise_z_scores( Y_orig );

        Col distances  = euclidean(X, obs_met);
        return ordered(distances);
    }

    std::vector<size_t> particle_ranking_PLS(
        const Mat2D & X_orig, const Mat2D & Y_orig,
        const Row & target_values,
        const float_type training_fraction
    ) {
        assert((0 < training_fraction) and (training_fraction <= 1));

        // Box-Cox transform data -- TODO?

        Row X_sim_means = X_orig.colwise().mean();
        Row X_sim_stdev = colwise_stdev(X_orig, X_sim_means);
        Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
        Mat2D Y = colwise_z_scores( Y_orig );
        Row obs_met = z_scores(target_values, X_sim_means, X_sim_stdev);

        const size_t pls_training_set_size = round(X.rows() * training_fraction);
        // @tjh TODO -- I think this may be a bug, and that ncomp should be equal to number of predictor variables (metrics in this case), not reponse variables
        size_t ncomp = Y_orig.cols();             // It doesn't make sense to consider more components than model parameters

        // assert: X / Y have matching order; X / Y are randomly ordered
        PLS::Model plsm(X.topRows(pls_training_set_size), Y.topRows(pls_training_set_size), ncomp);

        const int test_set_size = X.rows() - pls_training_set_size; // number of observations not in training set
        auto em = plsm.error<PLS::NEW_DATA>(X.bottomRows(test_set_size), Y.bottomRows(test_set_size));
        auto num_components = PLS::optimal_num_components(em);

        size_t num_components_used = num_components.maxCoeff();

        // Calculate new, orthogonal metrics (==scores) using the pls model
        // Is casting as real always safe?
        Row   obs_scores = plsm.scores(obs_met, num_components_used).row(0).real();
        Mat2D sim_scores = plsm.scores(X, num_components_used).real();
        Col   distances  = euclidean(sim_scores, obs_scores);

        return ordered(distances);
    }

    // use with params[set_num-1] (i.e. prior set, n-1, when on set n)
    // should already be sliced to posterior - i.e. _particle_parameters[set_num-1](_predictive_prior[set_num-1], Eigen::placeholders::all)
    gsl_matrix* setup_mvn_sampler(
        const Mat2D & params
    ) {
        // ALLOCATE DATA STRUCTURES
        // variance-covariance matrix calculated from pred prior values
        // NB: always allocate this small matrix, so that we don't have to check whether it's safe to free later
        const size_t num_pars = static_cast<size_t>(params.cols());
        gsl_matrix* sigma_hat = gsl_matrix_alloc(num_pars, num_pars);
        gsl_matrix* posterior_par_vals = to_gsl_m(params);

        // calculate maximum likelihood estimate of variance-covariance matrix sigma_hat
        gsl_ran_multivariate_gaussian_vcov(posterior_par_vals, sigma_hat);

        for (size_t parIdx = 0; parIdx < params.cols(); parIdx++) {
            // sampling is done using a kernel with a broader kernel than found in pred prior values
            const double doubled_variance = 2 * gsl_matrix_get(sigma_hat, parIdx, parIdx);
            gsl_matrix_set(sigma_hat, parIdx, parIdx, doubled_variance);
        }

        // not a nice interface, gsl.  sigma_hat is converted in place from a variance-covariance matrix
        // to the same values in the upper triangle, and the diagonal and lower triangle equal to L,
        // the Cholesky decomposition
        gsl_linalg_cholesky_decomp1(sigma_hat);
        gsl_matrix_free(posterior_par_vals);

        return sigma_hat;
    }

    Mat2D sample_priors(
        const gsl_rng* RNG, const size_t num_samples,
        const Mat2D & posterior, // look up table for POSTERIOR type Parameters
        const std::vector<const Parameter*> & mpars,
        std::vector<size_t> & post_ranks // filled in by this
    ) {
        // setup sampling RNG to deal w/ mixture of prior, posterior, pseudo parameters
        ParRNG par_rng(RNG, mpars, posterior.rows());
        // return matrix of samples
        Mat2D par_samples = Mat2D::Zero(num_samples, mpars.size());
        // determine which parameters are posterior vs not
        std::vector<size_t> nonpost_indices, post_indices;
        for (size_t parIdx = 0; parIdx < mpars.size(); parIdx++) {
            if (mpars[parIdx]->isPosterior()) {
                post_indices.push_back(parIdx);
            } else {
                nonpost_indices.push_back(parIdx);
            }
        }

        if (post_indices.size()) { post_ranks.resize(num_samples); }

        // confirm posterior matrix columns == number of posterior parameters
        assert(post_indices.size() == posterior.cols());

        for (size_t sampIdx = 0; sampIdx < par_samples.rows(); sampIdx++) { // for each sample to be drawn
            par_rng.unlock(); // allow the pseudo parameter to be incremented (if present, only the first one that qualifies)
            // for each non-posterior parameter
            for (size_t parIdx : nonpost_indices) { par_samples(sampIdx, parIdx) = mpars[parIdx]->sample(par_rng); }
            // if there are any posterior parameters, get the next posterior index (increments if par_rng is still unlocked)
            if (post_indices.size()) { post_ranks[sampIdx] = mpars[post_indices[0]]->sample(par_rng); }                
        }
        // if there are any posterior parameters, fill in the posterior samples
        if (post_ranks.size() != 0) { par_samples(Eigen::placeholders::all, post_indices) = posterior(post_ranks, Eigen::placeholders::all); }
        return par_samples;

    }

    Row calculate_doubled_variance(const Mat2D & params) {
        vector<RunningStat> stats(params.cols());
        Row v2 = Row::Zero(params.cols());
        // TODO: turn this into Eigen column-wise operation?
        for (size_t parIdx = 0; parIdx < params.cols(); parIdx++) {
            stats[parIdx].Push(params.col(parIdx));
            v2[parIdx] = 2 * stats[parIdx].Variance();
        }
        return v2;
    }

    Row weight_predictive_prior(
        const std::vector<const Parameter*> & mpars,
        const Mat2D & params
    ) {
        const float_type uniform_wt = 1.0/static_cast<double>(params.rows());
        return Col::Constant(params.rows(), uniform_wt);
    }

    Row weight_predictive_prior(
        const std::vector<const Parameter*> & mpars,
        const Mat2D & params,
        const Mat2D & prev_params,
        const Row & prev_weights,
        const Row & prev_doubled_variance
    ) {
        Col weight = Col::Zero(params.rows());

        for (size_t post_rank = 0; post_rank < params.rows(); post_rank++) {
            double numerator = 1;
            double denominator = 0.0;
            for (size_t parIdx = 0; parIdx < params.cols(); parIdx++) {
                numerator *= mpars[parIdx]->likelihood( params(post_rank, parIdx) );
            }

            for (size_t prev_post_rank = 0; prev_post_rank < prev_params.rows(); prev_post_rank++) {
                double running_product = prev_weights[prev_post_rank]; // TODO: rowwise op?
                for (size_t parIdx = 0; parIdx < prev_params.cols(); parIdx++) {
                    double par_value = params(post_rank, parIdx);
                    double old_par_value = prev_params(prev_post_rank, parIdx);
                    double old_dv = prev_doubled_variance[parIdx];

                    // This conditional handles the (often improbable) case where a parameter has completely converged.
                    // It allows ABC to continue exploring other parameters, rather than causing the math
                    // to fall apart because the density at the converged value is infinite.
                    if (old_dv != 0 or par_value != old_par_value) {
                        running_product *= gsl_ran_gaussian_pdf(par_value - old_par_value, sqrt(old_dv) );
                    }
                }
                denominator += running_product;
            }

            weight[post_rank] = numerator / denominator;
        }

        weight.normalize();

        return weight;
    }

}