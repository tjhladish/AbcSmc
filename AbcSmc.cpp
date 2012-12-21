#ifndef ABCSMC_H
#define ABCSMC_H

//#include "mpi.h"
#include <stdio.h>
#include <time.h>
#include <stdlib.h>
#include <sys/stat.h>
#include <vector>
#include <float.h>
#include <limits.h>
#include <string>
#include <sstream>
#include <math.h>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include "pls.h"
#include "RunningStat.h"

const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);


using namespace std;

inline double uniform_pdf(double a, double b) { return 1.0 / fabs(b-a); }

int gsl_rng_nonuniform_int(vector<double>& weights, const gsl_rng* rng) {
    // Weights must sum to 1!!
    double running_sum = 0;
    double r = gsl_rng_uniform(RNG);
    for (unsigned int i = 0; i<weights.size(); i++) {
        running_sum += weights[i];
        if (r<=running_sum) return i;
    }
    cerr << "Weights may not be normalized\n";
    exit(100);
}

double rand_trunc_normal(double mu, double sigma_squared, double min, double max, const gsl_rng* rng) {
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
 
std::string exec(std::string cmd) {
    FILE* pipe = popen(cmd.c_str(), "r");
    if (!pipe) return "ERROR";
    char buffer[512];
    std::string result = "";
    while(!feof(pipe)) {
        if(fgets(buffer, 512, pipe) != NULL)
            result += buffer;
    }
    pclose(pipe);
    return result;
}


template <typename T>
inline std::string toString (const T& t) {
    std::stringstream ss;
    ss << t;
    return ss.str();
}


enum PriorType {UNIFORM, NORMAL};
enum NumericType {INT, FLOAT};


class Parameter {
    public:
        Parameter() {};
        /*Parameter( string s, PriorType p, NumericType n, int val1, int val2 ) : name(s), ptype(p), ntype(n) {
            if (ptype == UNIFORM) {
                imin = val1;
                imax = val2;
                mean = (float) (val2 - val1) / 2.0;
                stdev = sqrt(pow(val2-val1,2)/12.0);
            } else if (ptype == NORMAL) {
                imin = INT_MIN;
                imax = INT_MAX;
                mean = val1;
                stdev = val2;
            }
        }*/

        Parameter( string s, PriorType p, NumericType n, double val1, double val2 ) : name(s), ptype(p), ntype(n) {
            assert(ptype == UNIFORM);
            if (ptype == UNIFORM) {
                fmin = val1;
                fmax = val2;
                mean = (val2 - val1) / 2.0;
                stdev = sqrt(pow(val2-val1,2)/12);
            } /*else if (ptype == NORMAL) {
                fmin = DBL_MIN;
                fmax = DBL_MAX;
                mean = val1;
                stdev = val2;
            }*/
        }

        double sample() { 
            if (ntype == INT) {
                 // + 1 makes it out of [fmin, fmax], instead of [fmin, fmax)
                return gsl_rng_uniform_int(RNG, fmax-fmin + 1) + fmin;
            } else { 
                return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin;
            }
        }

        //int sample_int() { return gsl_rng_uniform_int(RNG, imax-imin) + imin; }
        //double sample_float() { return gsl_rng_uniform(RNG)*(fmax-fmin) + fmin; }

        double get_doubled_variance(int t) { return doubled_variance[t]; }
        void append_doubled_variance(double v2) { doubled_variance.push_back(v2); }
        double get_min() { return fmin; }
        double get_max() { return fmax; }
        NumericType get_numeric_type() { return ntype; }
        

    private:
        string name;
        PriorType ptype;
        NumericType ntype;
        //int imin, imax;
        double fmin, fmax, mean, stdev;
        vector<double> doubled_variance;

//num != static_cast<int>(num)
};

class ModelParameters {
    public:
        ModelParameters() {};

        int size() { return _pars.size(); }

        //Parameter* p = new Parameter( name, PriorType, NumericType, val1, val2, );
        void add_next_parameter(Parameter* p) { _pars.push_back(p); }


        void add_next_parameter(string name, PriorType ptype, NumericType ntype, double val1, double val2) {
            Parameter* p = new Parameter(name, ptype, ntype, val1, val2);
            _pars.push_back(p);
        }

        vector<Parameter*> get_parameters() { return _pars; }
        Parameter* get_parameter(int i) { return _pars[i]; }

    private:
        vector<Parameter*> _pars;
};


struct Metric {
    Metric() {};
    Metric(string n, NumericType nt, double val) : name(n), ntype(nt), obs_val(val) {};
    string name;
    NumericType ntype;
    double obs_val;
};


class ModelMetrics {
    public:
        void add_next_metric(string met, NumericType ntype, double obs_val) { _mets.push_back(new Metric(met, ntype, obs_val)); }
        vector<Metric*> get_metrics() { return _mets; }
        int size() { return _mets.size(); }

    private:
        vector<Metric*> _mets;
};


/*class Particle {
    public:
        //add_particle(vector<double> par_vals, vector<double> met_vals) { _pars=par_vals; _mets=met_vals; }
        const inline double get_parameter_value( int i ) { return _pars[i]; } 
        vector<double> get_parameter_values() { return _pars; } 
        void set_parameter_values(vector<double> p ) { _pars = p; } 
        vector<double> get_metric_values() { return _mets; } 
        void set_metric_values(vector<double> m) { _mets = m; } 

        friend ostream& operator<<(ostream& os, Particle& par);
        //friend ostream& operator<<(ostream& os, const Particle& par);

    private:
        vector<double> _pars;
        vector<double> _mets;
};

ostream& operator<<(ostream& os, Particle& particle) {
//ostream& operator<<(ostream& os, const Particle& particle) {
    vector<double> _pars = particle.get_parameter_values();
    vector<double> _mets = particle.get_metric_values();

    for (unsigned int i = 0; i<_pars.size(); i++) os << _pars[i] << " ";
    for (unsigned int i = 0; i<_mets.size() - 1; i++) os << _mets[i] << " ";
    os << _mets.back();
    return os;
}*/



class AbcSmc {
    public:
        AbcSmc() {};
        AbcSmc(ModelParameters* pars, ModelMetrics* mets) { _model_pars = pars; _model_mets = mets; }

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_predictive_prior_size(int n) { assert(n > 0); assert(n <= _num_particles); _predictive_prior_size = n; }
        void set_predictive_prior_fraction(float f)        { assert(f > 0); assert(f <= 1); _predictive_prior_size = _num_particles * f; }
        void set_pls_validation_training_fraction(float f) { assert(f > 0); assert(f <= 1); _pls_training_set_size = _num_particles * f; }
        void set_metric_basefilename( string name ) { _metrics_filename = name; }
        //void set_executable_filename( string name ) { _executable_filename = name; }
        //void set_particle_basefilename( string name ) { _particle_filename = name; }
        //void set_predictive_prior_basefilename( string name ) { _predictive_prior_filename = name; }
        

        void run(string executable) {
        // NEED TO ADD SANITY CHECKS HERE
            _executable_filename = executable;
            _particle_parameters.clear();
            _particle_metrics.clear();
            _weights.clear();

            _particle_parameters.resize( _num_smc_sets, Mat2D::Zero(_num_particles, npar()) );
            _particle_metrics.resize( _num_smc_sets, Mat2D::Zero(_num_particles, nmet()) );
            _weights.resize(_num_smc_sets);

            //_particles.clear();
            _predictive_prior.clear();
            _predictive_prior.resize(_num_smc_sets);

            //_particles.resize( _num_smc_sets, vector<Particle*>(_num_particles) );
            //_predictive_prior.resize( _num_smc_sets, vector<Particle*>(_predictive_prior_size) );

            for (int t = 0; t<_num_smc_sets; t++) {
                //X_orig = Mat2D::Zero(_num_particles, _model_mets->size());
                //Y_orig = Mat2D::Zero(_num_particles, _model_pars->size());
cerr << "Sampling\n";
                //_populate_particles( t, X_orig, Y_orig );
                _populate_particles( t, _particle_metrics[t], _particle_parameters[t]);
cerr << "Done\n";
                // Write out particles

                // Select best scoring particles
                // "sample" contains indices (==row #) of best scoring particles
                _filter_particles( t, _particle_metrics[t], _particle_parameters[t]);
                //vector<int> sample = _filter_particles( t, X_orig, Y_orig );

cerr << "Calculating weights\n";
                calculate_predictive_prior_weights( t );
cerr << "Done\n";
                // Write out predictive prior
            }
        }

        int npar() { return _model_pars->size(); }
        int nmet() { return _model_mets->size(); }

    private:
        Mat2D X_orig;
        Mat2D Y_orig;
        ModelParameters* _model_pars;
        ModelMetrics* _model_mets;
        int _num_smc_sets;
        int _num_particles;
        int _pls_training_set_size;
        int _predictive_prior_size; // number of particles that will be used to inform predictive prior
        string _executable_filename;
        string _metrics_filename;
        string _particle_filename;
        string _predictive_prior_filename;
        vector< Mat2D > _particle_metrics;
        vector< Mat2D > _particle_parameters;
        vector< vector<int> > _predictive_prior; // vector of row indices for particle metrics and parameters
        //vector< vector<Particle*> > _particles;
        //vector< vector<Particle*> > _predictive_prior;
        vector< vector<double> > _weights;

        void _populate_particles(int t, Mat2D &X_orig, Mat2D &Y_orig) {
            for (int i = 0; i<_num_particles; i++) {
                if (t == 0) { // sample priors
                    Y_orig.row(i) = sample_priors();
                } else { // sample predictive priors
//Y_orig.row(i) = sample_priors();
                    Y_orig.row(i) = sample_predictive_priors(t);
                }

                //string output_filename = "tmp_sim_" + toString(t) + "_" + toString(i) + ".out"; 
                //string command = executable + " " + output_filename;
                string command = _executable_filename;
                for (int j = 0; j<npar(); j++) {
                    command += " " + toString(Y_orig(i,j));
                }
                // ./executable_name summary_stats par1val par2val par3val par4val par5val ... 

                //std::cerr << command << endl;
                //string retval = exec("echo 19 2.3");
                string retval = exec(command);
                stringstream ss;
                ss.str(retval);
                //std::cerr << "nmet=" << nmet() << " captured=" << retval << endl;
            
                for (int j = 0; j<nmet(); j++) {
                    ss >> X_orig(i, j);
                }
            }
        }

        void _filter_particles (int t, Mat2D &X_orig, Mat2D &Y_orig) {
        //vector<int> _filter_particles (int t, Mat2D &X_orig, Mat2D &Y_orig) {
            // Run PLS
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

            cout << "scores:\t" << plsm.scores(obs_met) << endl;
            cout << "DONE!!!\n";
            vector<int>::iterator first = ranking.begin();
            vector<int>::iterator last  = ranking.begin() + _predictive_prior_size;
            vector<int> sample(first, last);
            _predictive_prior[t] = sample;

/*cerr <<  "pred prior:\n";
            for (unsigned int i=0; i<sample.size(); i++) {
cerr << "i, sample[i], particle:" << i << " " << sample[i] << " ";
cerr << _particle_parameters[t].row(sample[i]) << " | " << _particle_metrics[t].row(sample[i]) << endl;
            }
cerr <<  "\n";
*/
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

cout << "Actual values:\n" << "100 6 | 374 1.71517\n";

            //return sample;
        }
        
        Col euclidean( Row obs_met, Mat2D sim_met ) {
            Col distances = Col::Zero(sim_met.rows());
            //vector<double> distances(sim_met.rows());
            for (int r = 0; r<sim_met.rows(); r++) {
                for (int c = 0; c<sim_met.cols(); c++) {
                    distances(r) += pow(obs_met(c) - sim_met(r,c), 2);
                }
                distances(r) = sqrt( distances(r) );
            }
            return distances; 
        }

        Row sample_priors() {
            vector<Parameter*> par_vector = _model_pars->get_parameters();
            //vector<double> par_sample;
            Row par_sample = Row::Zero(_model_pars->size());

            for (int i = 0; i < _model_pars->size(); i++) {
                par_sample(i) =  par_vector[i]->sample();
            }
            return par_sample;
        }

        void calculate_doubled_variances( int t ) {
            vector<RunningStat> stats(npar());

            for (unsigned int i = 0; i < _predictive_prior[t].size(); i++) {
                for (int j = 0; j < npar(); j++) {
                    double particle_idx = _predictive_prior[t][i];
                    //double par_value = _predictive_prior[t][i]->get_parameter_value(j);
                    double par_value = _particle_parameters[t](particle_idx, j);
                    stats[j].Push(par_value);
                }
            }

            for (int j = 0; j < npar(); j++) {
                _model_pars->get_parameter(j)->append_doubled_variance( 2 * stats[j].Variance() );
            }
        }

        void normalize_weights( vector<double>& weights ) {
            double total = 0;
            for (unsigned int i = 0; i< weights.size(); i++) {
                total += weights[i];
            }

            for (unsigned int i = 0; i< weights.size(); i++) {
                weights[i] /= total;
            }
        }

        void calculate_predictive_prior_weights(int set_num) {
            // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
            calculate_doubled_variances( set_num );
            cerr << "pred prior size in calculate_predictive_prior_weights: " << _predictive_prior[set_num].size() << endl;
            
            if (set_num == 0) {
                // uniform weights for set 0 predictive prior
                _weights[set_num].resize(_predictive_prior[set_num].size(), 1.0/(double) _predictive_prior[set_num].size());
            } else if ( set_num > 0 ) {
                // weights from set - 2 are needed to calculate weights for set - 1

                _weights[set_num].resize( _predictive_prior[set_num].size() );
                
                double numerator = 1;
                for (int j = 0; j<npar(); j++) {
                    Parameter* par = _model_pars->get_parameter(j);
                    numerator *= uniform_pdf(par->get_min(), par->get_max());
                }

                for (unsigned int i = 0; i < _predictive_prior[set_num].size(); i++) {
                    double denominator = 0.0;
                    for (int k = 0; k < _predictive_prior[set_num - 1].size(); k++) {
                        double running_product = _weights[set_num - 1][k];
                        for (int j = 0; j < npar(); j++) {
                            //double par_value = _predictive_prior[set_num][i]->get_parameter_value(j);
                            double par_value = _particle_parameters[set_num](_predictive_prior[set_num][i], j);
                            //double old_par_value = _predictive_prior[set_num-1][k]->get_parameter_value(j);
                            double old_par_value = _particle_parameters[set_num-1](_predictive_prior[set_num-1][k], j);
                            double old_doubled_variance = _model_pars->get_parameter(j)->get_doubled_variance(set_num-1);

                            //running_product *= normal_pdf(par_value, old_par_value, old_doubled_variance );
                            running_product *= gsl_ran_gaussian_pdf(par_value-old_par_value, sqrt(old_doubled_variance) );
                        }
                        denominator += running_product;
                    }
                    _weights[set_num][i] = numerator / denominator;
                }

                normalize_weights( _weights[set_num] );
            }
            //if (mpi_rank == 0) write_weights_file( weights, set_num - 1 );
        }

        Row sample_predictive_priors(int set_num) {
            // Select a particle index j to use from the predictive prior
            Row par_values = Row::Zero(npar());
            int r = gsl_rng_nonuniform_int(_weights[set_num-1], RNG);
            for (int j = 0; j<npar(); j++) {
                int particle_idx = _predictive_prior[set_num-1][r];
                //double par_value = particle->get_parameter_value(j);
                double par_value = _particle_parameters[set_num-1](particle_idx, j);
                Parameter* parameter = _model_pars->get_parameter(j);
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

        Row _z_transform_observed_metrics(Row& means, Row& stdevs) {
            Row zmat = Row::Zero(nmet());
            vector<Metric*> mets = _model_mets->get_metrics();
            for (int i = 0; i<nmet(); i++) { zmat(i) = (mets[i]->obs_val - means(i)) / stdevs(i); }
            return zmat;
        }
};

int main(int argc, char* argv[]) {
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id
    ModelParameters* par = new ModelParameters();
    par->add_next_parameter( new Parameter("ndice", UNIFORM, INT, 1, 1000) );
    par->add_next_parameter( new Parameter("sides", UNIFORM, INT, 1, 1000) );
    
    ModelMetrics* met = new ModelMetrics(); // the name and order of summary stats returned by simulator
    met->add_next_metric( "sum", INT, 374 );     // are we actually going to use the numeric type attribute?
    met->add_next_metric( "sd", FLOAT, 1.71517 );

    AbcSmc* abc = new AbcSmc(par, met);
    abc->set_smc_iterations(5); // or have it test for convergence
    abc->set_num_samples(1000);
    abc->set_predictive_prior_fraction(0.1);
    abc->set_pls_validation_training_fraction(0.5); // fraction of runs to use for training (vs. testing) pls model
    //abc->set_metric_basefilename("summary_stats");
    abc->run("/home/tjhladish/work/abc_cpp/dice_game");  // ./executable_name summary_stats par1val par2val par3val par4val par5val ...
    //abc->run("/home/tjhladish/work/abc_cpp/dice_game.py");  // ./executable_name summary_stats par1val par2val par3val par4val par5val ...





/*
    MPI_Init(&argc, &argv);
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    fprintf(stderr, " %d starting\n", rank );
    //if (rank == 0 ) { fprintf(stderr, "TEST: dynamic network seed, dynamic parameters; reusing network object\n" ); }
    
    if (argc < 2) { fprintf(stderr, "Provide configuration file name as a parameter\n"); }
    string config_file_name( argv[1] );

    // Process configuration file
    Parameters par;
    process_config_file(config_file_name, par);

    // Determine which ABC set we're supposed to run
    const int set_num = determine_set_number();

    // Determine what parameter combinations we will try
    static MTRand mtrand;
    vector<Particle> particles;
    if (set_num == 0) {
        // sample naive priors
        sample_prior(particles, par.sample_size, &mtrand);
    } else {
        // sample predictive priors based on last set
        sample_predictive_prior(set_num, particles, par.sample_size, &mtrand, rank);
    }

    map<string, vector<float> > obs_data; //string = location, float = incidence on [0,1]
    load_observed_data(par.observed_data_file, ',', obs_data);

    Network* net = new Network("EpiNet", Network::Undirected);

    for (unsigned int i = 0; i < particles.size(); i++) {

        const double R_zero       = particles[i].R0;
        const double Ih           = particles[i].Ih;
        const double h            = particles[i].h;
        const int patient_zero_ct = particles[i].P0;
        const int burnin          = par.burnin;
        const int net_size        = par.network_size;

        map<string, vector<float> > sim_data;
        map<string, vector<float> > R0_vals;

        net->clear_nodes(); // debugging
        MultiSeason_Sim* sim = new MultiSeason_Sim(net, Ih);
        generate_network(net, par, R_zero, sim);
        double new_R_zero = R_zero;

        map<string, vector<float> >::iterator it;

        for ( it = obs_data.begin(); it != obs_data.end(); it++ ) {
            double Tc_actual = sim->calc_critical_transmissibility();
            string loc = (*it).first;
            const int obs_N = obs_data[loc].size();

            sim_data[loc].resize(obs_N);
            R0_vals[loc].resize(obs_N);

            for ( int season = 0; season < burnin + obs_N; season++) {
                sim->rand_infect(patient_zero_ct);
                sim->run_simulation();
                                        
                if (season >= burnin) {
                    // epi size in percent, reduced by hospitalization factor
                    const double transmitted_size = double(sim->epidemic_size() - patient_zero_ct)/(net_size - patient_zero_ct);
                    sim_data[loc][season - burnin] =  h * transmitted_size;
                    R0_vals[loc][season - burnin]  = new_R_zero;
                }

                // now calculate what R_zero will be at the start of the next season
                vector<double> average_tk;
                double average_t = 0;
                sim->calculate_average_transmissibility(average_tk, average_t);
                new_R_zero = average_t / Tc_actual;
            }
            sim->reset();
            new_R_zero = R_zero;
        }
        // Report parameters
        report_particles(rank, i, particles[i], sim_data, R0_vals);
        delete sim;
    }

    delete net;
    fprintf( stderr, "%d done\n", rank );
    MPI_Finalize();
*/
    return 0;
}


/*void load_observed_data(string filename, char sep, map<string, vector<float> > &data);


void process_config_file (string config_filename, Parameters &p) {

    std::map<string,string> tmp_par;
    ifstream myfile(config_filename.c_str());
    std::stringstream ss;

    if (myfile.is_open()) {
        string line;
        char sep = '\t';

        while ( getline(myfile, line) ) {
            vector<string> fields;
            split(line, sep, fields);

            const char whitespace[] = " \n\t\r";
            string key =   strip( fields[0], whitespace ); 
            string value = strip( fields[1], whitespace ); 

            tmp_par[ key ] = value;
        }
    } else {
        fprintf( stderr, "Failed to open parameter file\n" );
        exit(102);
    }

    p.observed_data_file = tmp_par["observed_data_file"];
    p.network_size = to_int( tmp_par["network_size"] );
    p.burnin       = to_int( tmp_par["burnin"] );
    p.sample_size  = to_int( tmp_par["sample_size"] );  
    p.epsilon      = to_double( tmp_par["epsilon"] );
    p.threads      = to_int( tmp_par["threads"] );
    p.obs_mean     = to_double( tmp_par["obs_mean"] );
    p.obs_median   = to_double( tmp_par["obs_median"] );
    p.obs_max      = to_double( tmp_par["obs_max"] );
    p.obs_range    = to_double( tmp_par["obs_range"] );
    p.obs_sd       = to_double( tmp_par["obs_sd"] );
    p.obs_skew     = to_double( tmp_par["obs_skew"] );
    p.obs_ss       = to_double( tmp_par["obs_ss"] );
    p.obs_ll       = to_double( tmp_par["obs_ll"] );
    p.obs_sl       = to_double( tmp_par["obs_sl"] );
    p.obs_ls       = to_double( tmp_par["obs_ls"] );

    return;
}


bool fileExists(string strFilename) { 
    struct stat stFileInfo; 
    bool blnReturn; 
    int intStat; 
    // Attempt to get the file attributes 
    intStat = stat(strFilename.c_str(),&stFileInfo); 
    if(intStat == 0) { 
        blnReturn = true; 
    } else { 
        blnReturn = false;
    }
    return(blnReturn);
}

vector<Particle> read_predictive_prior_file(int set_num) {
    string filename = "predictive_prior." + to_string(set_num);
// 'R0', 'Ih', 'h', 'P0', 'mean', 'median', 'max', 'range', 'sd', 'skew', 'ss', 'll', 'sl', 'ls', 'Re'   
    ifstream myfile(filename.c_str());
    std::stringstream ss;
    vector<Particle> predictive_prior;

    if (myfile.is_open()) {
        string line;
        char sep = ' ';

        while ( getline(myfile, line) ) {
            vector<string> fields;
            split(line, sep, fields);
            Particle p;

            const char whitespace[] = " \n\t\r";
            string R0_str = strip( fields[0], whitespace ); 
            string Ih_str = strip( fields[1], whitespace ); 
            string h_str  = strip( fields[2], whitespace ); 
            string P0_str = strip( fields[3], whitespace ); 

            p.R0 = to_double( R0_str );
            p.Ih = to_double( Ih_str );
            p.h  = to_double( h_str );
            p.P0 = to_int( P0_str );

            predictive_prior.push_back(p);
        }
    } else {
        fprintf( stderr, "Failed to open parameter file\n" );
        exit(103);
    }
    return predictive_prior;
}

vector<double> read_weights_file(int set_num) {
    string filename = "weights." + to_string(set_num);

    vector<double> weights;
    ifstream myfile(filename.c_str());
    std::stringstream ss;

    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            //split string based on "," and store results into vector
            const char whitespace[] = " \n\t\r";
            string val_str = strip( line, whitespace );
        
            weights.push_back( atof( val_str.c_str() ) );
        }
    }
    return weights;
}


void write_weights_file(vector<double> weights, int set_num) {
    string filename = "weights." + to_string(set_num);

    ofstream pipe(filename.c_str(), ios::out);
    for (unsigned int i = 0; i < weights.size(); i++) {
        pipe << weights[i] << endl;
    }
    pipe.close();
    return;
}




void report_particles(int rank, int run, Particle &particle, map<string, vector<float> > sim, map<string, vector<float> > Re_values) {
    // obs_mean, obs_median, obs_max, obs_range, obs_sd, obs_skew, obs_ss, obs_ll, obs_sl, obs_ls
    vector<float> sim_flat = flatten_map( sim );
    vector<float> Re_flat  = flatten_map( Re_values );

    double mean_sim = mean(sim_flat);
    double median_sim = median(sim_flat);
    double max_sim = max_element(sim_flat); 
    double range_sim = range(sim_flat); 
    double sd_sim = stdev(sim_flat);
    double skew_sim = mean_sim - median_sim;
    vector<double> acm_sim = autocorrelation_matrix(sim);

    float epi_threshold = 0.05;
    double mean_Re = calculate_mean_effective_R (sim_flat, Re_flat, epi_threshold);
    //cout << particles[i].R0 << " " << particles[i].Ih << " " << particles[i].h << " " << particles[i].P0 << " ";
    fprintf( stderr, "%d %d   %g %g %g %d   %g %g %g %g %g %g   %g %g %g %g   %g\n", 
        rank, run,
        particle.R0, particle.Ih, particle.h, particle.P0,
        mean_sim, median_sim, max_sim, range_sim, sd_sim, skew_sim,
        acm_sim[0], acm_sim[1], acm_sim[2], acm_sim[3],
        mean_Re);
        
    return;
}
*/
/*
Determine run number: what is the last predictive_prior.X file? (or does one not exist->first run)

If first time being run, no input files exist.  Sample from uniform priors, output particles.

If second time, read in predictive prior from first run, use uniform weights. Output uniform
weights as weights for first run, output second run particles.

If third+ time, read in predictive prior from last run, weights from second to last run,
calculate and output weights for last run; output particles.

*/



/*
void load_observed_data(string filename, char sep, map<string, vector<float> > & obs_data) {
    ifstream myfile(filename.c_str());

    if (myfile.is_open()) {
        string line;

        while ( getline(myfile,line) ) {
            //split string based on "," and store results into vector
            vector<string> fields;
            split(line, sep, fields);
            const char whitespace[] = " \n\t\r";

            //format check
            if (fields.size() > 3 ) {
                fprintf( stderr, "Skipping line: too many fields: %s\n", line.c_str() );
                continue;
            } else if (fields.size() < 3 ) {
                fprintf( stderr, "Skipping line: too few fields:  %s\n", line.c_str() );
                continue;
            } else { 
                string loc   = strip(fields[0],whitespace);
                string year  = strip(fields[1],whitespace);
                float  count = to_float(strip(fields[2],whitespace));
                obs_data[loc].push_back(count); 
            }
        }
    }
    return;
}
*/



#endif
