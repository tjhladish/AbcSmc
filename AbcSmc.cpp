#include "AbcSmc.h"
#include "pls.h"
#include "RunningStat.h"

using std::vector;
using std::string;
using std::stringstream;

void AbcSmc::run(string executable, const gsl_rng* RNG ) {
    // NEED TO ADD SANITY CHECKS HERE
    _executable_filename = executable;
    _particle_parameters.clear();
    _particle_metrics.clear();
    _weights.clear();

    _particle_parameters.resize( _num_smc_sets, Mat2D::Zero(_num_particles, npar()) );
    _particle_metrics.resize( _num_smc_sets, Mat2D::Zero(_num_particles, nmet()) );
    _weights.resize(_num_smc_sets);

    _predictive_prior.clear();
    _predictive_prior.resize(_num_smc_sets);

    for (int t = 0; t<_num_smc_sets; t++) {
        _populate_particles( t, _particle_metrics[t], _particle_parameters[t], RNG );
        // Write out particles

        // Select best scoring particles
        _filter_particles( t, _particle_metrics[t], _particle_parameters[t]);

        calculate_predictive_prior_weights( t );
        // Write out predictive prior
    }
}

void AbcSmc::_populate_particles(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {
    for (int i = 0; i<_num_particles; i++) {
        if (t == 0) { // sample priors
            Y_orig.row(i) = sample_priors(RNG);
        } else { // sample predictive priors
            Y_orig.row(i) = sample_predictive_priors(t, RNG);
        }

        //string output_filename = "tmp_sim_" + toString(t) + "_" + toString(i) + ".out"; 
        //string command = executable + " " + output_filename;
        string command = _executable_filename;
        for (int j = 0; j<npar(); j++) {
            command += " " + toString(Y_orig(i,j));
        }
        // ./executable_name summary_stats par1val par2val par3val par4val par5val ... 

        //std::cerr << command << endl;
        string retval = exec(command);
        stringstream ss;
        ss.str(retval);
    
        for (int j = 0; j<nmet(); j++) {
            ss >> X_orig(i, j);
        }
    }
}

void AbcSmc::_filter_particles (int t, Mat2D &X_orig, Mat2D &Y_orig) {
    // Run PLS
    // Box-Cox transform data?
    //void test_bc( Mat2D );
    //test_bc(Y_orig);

    Row X_sim_means, X_sim_stdev;
    Mat2D X = colwise_z_scores( X_orig, X_sim_means, X_sim_stdev );
    Mat2D Y = colwise_z_scores( Y_orig );
    Row obs_met = _z_transform_observed_metrics( X_sim_means, X_sim_stdev );

    int nobs  = X_orig.rows();      // number of observations
    int npred = X_orig.cols();      // number of predictor variables
    int nresp = Y_orig.cols();      // number of response variables
    int ncomp = npar();             // It doesn't make sense to consider more components than model parameters
    PLS_Model plsm;
    plsm.initialize(npred, nresp, ncomp);
    plsm.plsr(X.topRows(_pls_training_set_size), Y.topRows(_pls_training_set_size), plsm, KERNEL_TYPE1);

    // A is number of components to use
    for (int A = 1; A<=ncomp; A++) { 
        // How well did we do with this many components?
        cout << A << " components\t";
        cout << "explained variance: " << plsm.explained_variance(X, Y, A);
        //cout << "root mean squared error of prediction (RMSEP):" << plsm.rmsep(X, Y, A) << endl;
        cout << " SSE: " << plsm.SSE(X,Y,A) <<  endl; 
    }

    int test_set_size = nobs - _pls_training_set_size;
    Rowi num_components = plsm.optimal_num_components(X.bottomRows(test_set_size), Y.bottomRows(test_set_size), NEW_DATA);
    int max_num_components = num_components.maxCoeff();
    cout << "Optimal number of components (NEW DATA):\t" << num_components << endl;

    // Calculate new, orthogonal metrics (==scores) using the pls model
    // Is casting as real always safe?
    Row   obs_scores = plsm.scores(obs_met).row(0).real();
    Mat2D sim_scores = plsm.scores(X).real();
    Col   distances  = euclidean(obs_scores, sim_scores);
    vector<int> ranking = ordered(distances);

    vector<int>::iterator first = ranking.begin();
    vector<int>::iterator last  = ranking.begin() + _predictive_prior_size;
    vector<int> sample(first, last);
    _predictive_prior[t] = sample;

    cout << "Worst five:\n";
    for (unsigned int q=ranking.size()-1;q>=ranking.size()-5;q--) {
        int idx = ranking[q];
        cout << Y_orig.row(idx) << " | " << X_orig.row(idx) << endl;
    }

    cout << "Best five:\n";
    for (int q=0;q<5;q++) {
        int idx = ranking[q];
        cout << Y_orig.row(idx) << " | " << X_orig.row(idx) << endl;
    }
}

Col AbcSmc::euclidean( Row obs_met, Mat2D sim_met ) {
    Col distances = Col::Zero(sim_met.rows());
    for (int r = 0; r<sim_met.rows(); r++) {
        for (int c = 0; c<sim_met.cols(); c++) {
            distances(r) += pow(obs_met(c) - sim_met(r,c), 2);
        }
        distances(r) = sqrt( distances(r) );
    }
    return distances; 
}

Row AbcSmc::sample_priors(const gsl_rng* RNG) {
    Row par_sample = Row::Zero(_model_pars.size());

    for (unsigned int i = 0; i < _model_pars.size(); i++) {
        par_sample(i) =  _model_pars[i]->sample(RNG);
    }
    return par_sample;
}

void AbcSmc::calculate_doubled_variances( int t ) {
    vector<RunningStat> stats(npar());

    for (unsigned int i = 0; i < _predictive_prior[t].size(); i++) {
        for (int j = 0; j < npar(); j++) {
            double particle_idx = _predictive_prior[t][i];
            double par_value = _particle_parameters[t](particle_idx, j);
            stats[j].Push(par_value);
        }
    }

    for (int j = 0; j < npar(); j++) {
        _model_pars[j]->append_doubled_variance( 2 * stats[j].Variance() );
    }
}

void AbcSmc::normalize_weights( vector<double>& weights ) {
    double total = 0;
    for (unsigned int i = 0; i< weights.size(); i++) {
        total += weights[i];
    }

    for (unsigned int i = 0; i< weights.size(); i++) {
        weights[i] /= total;
    }
}

void AbcSmc::calculate_predictive_prior_weights(int set_num) {
    // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
    calculate_doubled_variances( set_num );
    cerr << "pred prior size in calculate_predictive_prior_weights: " << _predictive_prior[set_num].size() << endl;
    
    if (set_num == 0) {
        // uniform weights for set 0 predictive prior
        _weights[set_num].resize(_predictive_prior[set_num].size(), 1.0/(double) _predictive_prior[set_num].size());
    } else if ( set_num > 0 ) {
        // weights from set - 1 are needed to calculate weights for current set

        _weights[set_num].resize( _predictive_prior[set_num].size() );
        
        double numerator = 1;
        for (int j = 0; j<npar(); j++) {
            Parameter* par = _model_pars[j];
            numerator *= uniform_pdf(par->get_min(), par->get_max());
        }

        for (unsigned int i = 0; i < _predictive_prior[set_num].size(); i++) {
            double denominator = 0.0;
            for (int k = 0; k < _predictive_prior[set_num - 1].size(); k++) {
                double running_product = _weights[set_num - 1][k];
                for (int j = 0; j < npar(); j++) {
                    double par_value = _particle_parameters[set_num](_predictive_prior[set_num][i], j);
                    double old_par_value = _particle_parameters[set_num-1](_predictive_prior[set_num-1][k], j);
                    double old_doubled_variance = _model_pars[j]->get_doubled_variance(set_num-1);

                    running_product *= gsl_ran_gaussian_pdf(par_value-old_par_value, sqrt(old_doubled_variance) );
                }
                denominator += running_product;
            }
            _weights[set_num][i] = numerator / denominator;
        }

        normalize_weights( _weights[set_num] );
    }
}

Row AbcSmc::sample_predictive_priors(int set_num, const gsl_rng* RNG ) {
    Row par_values = Row::Zero(npar());
    // Select a particle index r to use from the predictive prior
    int r = gsl_rng_nonuniform_int(_weights[set_num-1], RNG);
    for (int j = 0; j<npar(); j++) {
        int particle_idx = _predictive_prior[set_num-1][r];
        double par_value = _particle_parameters[set_num-1](particle_idx, j);
        Parameter* parameter = _model_pars[j];
        double doubled_variance = parameter->get_doubled_variance(set_num-1);
        double par_min = parameter->get_min();
        double par_max = parameter->get_max();
        par_values(j) = rand_trunc_normal( par_value, doubled_variance, par_min, par_max, RNG );

        if (parameter->get_numeric_type() == INT) {
            par_values(j) = (double) ((int) (par_values(j) + 0.5));
        }
    }
    return par_values;
}

Row AbcSmc::_z_transform_observed_metrics(Row& means, Row& stdevs) {
    Row zmat = Row::Zero(nmet());
    for (int i = 0; i<nmet(); i++) { zmat(i) = (_model_mets[i]->obs_val - means(i)) / stdevs(i); }
    return zmat;
}

