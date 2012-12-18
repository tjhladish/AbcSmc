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

const gsl_rng* RNG = gsl_rng_alloc (gsl_rng_taus2);


using namespace std;

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
            if (ptype == UNIFORM) {
                fmin = val1;
                fmax = val2;
                mean = (val2 - val1) / 2.0;
                stdev = sqrt(pow(val2-val1,2)/12);
            } else if (ptype == NORMAL) {
                fmin = DBL_MIN;
                fmax = DBL_MAX;
                mean = val1;
                stdev = val2;
            }
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

    private:
        string name;
        PriorType ptype;
        NumericType ntype;
        //int imin, imax;
        double fmin, fmax, mean, stdev;

//num != static_cast<int>(num)
};

class ModelParameters {
    public:
        ModelParameters() {};

        int size() { return _pars.size(); }

        //Parameter* p = new Parameter( name, PriorType, NumericType, val1, val2, );
        void add_parameter(Parameter* p) { _pars.push_back(p); }


        void add_parameter(string name, PriorType ptype, NumericType ntype, double val1, double val2) {
            Parameter* p = new Parameter(name, ptype, ntype, val1, val2);
            _pars.push_back(p);
        }

        vector<Parameter*> get_parameters() { return _pars; }

    private:
        vector<Parameter*> _pars;
};


struct Metric {
    Metric() {};
    Metric(string n, NumericType nt) : name(n), ntype(nt) {};
    string name;
    NumericType ntype;
};


class ModelMetrics {
    public:
        void add_metric(string met, NumericType ntype) { _mets.push_back(new Metric(met, ntype)); }
        vector<Metric*> get_metrics() { return _mets; }
        int size() { return _mets.size(); }

    private:
        vector<Metric*> _mets;
};


class Particle {
    public:
        //add_particle(vector<double> par_vals, vector<double> met_vals) { _pars=par_vals; _mets=met_vals; }
        vector<double> get_parameter_values() { return _pars; } 
        void set_parameter_values(vector<double> p ) { _pars = p; } 
        vector<double> get_metric_values() { return _mets; } 
        void set_metric_values(vector<double> m) { _mets = m; } 

    private:
        vector<double> _pars;
        vector<double> _mets;
};


class AbcSmc {
    public:
        AbcSmc() {};
        AbcSmc(ModelParameters* pars, ModelMetrics* mets) { _model_pars = pars; _model_mets = mets; }

        void set_smc_iterations(int n) { _num_smc_sets = n; }
        void set_num_samples(int n) { _num_particles = n; }
        void set_metric_basefilename( string name ) { _metrics_filename = name; }
        //void set_executable_filename( string name ) { _executable_filename = name; }
        //void set_particle_basefilename( string name ) { _particle_filename = name; }
        //void set_predictive_prior_basefilename( string name ) { _predictive_prior_filename = name; }
        

        void run(string executable) {
            _executable_filename = executable;
            _particles.clear();
            _predictive_prior.clear();

            _particles.resize( _num_smc_sets, vector<Particle*>(_num_particles) );
            _predictive_prior.resize( _num_smc_sets, vector<Particle*>(_num_particles*_particle_threshold) );

            for (int t = 0; t<_num_smc_sets; t++) {
                X_orig = Mat2D::Zero(_num_particles, _model_mets->size());
                Y_orig = Mat2D::Zero(_num_particles, _model_pars->size());
                _populate_particles( t, X_orig, Y_orig );
                // Write out particles

                // Select best scoring particles
                vector<int> sample = _filter_particles( X_orig, Y_orig );

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
        int _particle_threshold; // fraction of particles that will be used to inform predictive prior
        string _executable_filename;
        string _metrics_filename;
        string _particle_filename;
        string _predictive_prior_filename;
        vector< vector<Particle*> > _particles;
        vector< vector<Particle*> > _predictive_prior;

        void _populate_particles(int t, Mat2D &X_orig, Mat2D &Y_orig) {
            for (int i = 0; i<_num_particles; i++) {
                if (t == 0) { // sample priors
                    Y_orig.row(i) = sample_priors();
                } else { // sample predictive priors
Y_orig.row(i) = sample_priors();
//                    Y_orig.row(i) = sample_predictive_priors();
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

        vector<int> _filter_particles (Mat2D &X_orig, Mat2D &Y_orig) {
            // Run PLS
            //void test_bc( Mat2D );
            //test_bc(Y_orig);

            Mat2D X = colwise_z_scores( X_orig );
            Mat2D Y = colwise_z_scores( Y_orig );

            PLS_Model plsm;
            int nobs  = X_orig.rows();      // number of observations
            int npred = X_orig.cols();      // number of predictor variables
            int nresp = Y_orig.cols();      // number of response variables
            int ncomp = npar();             // It doesn't make sense to consider more components than model parameters
            plsm.initialize(npred, nresp, ncomp);
            plsm.plsr(X, Y, plsm, KERNEL_TYPE1);

            // A is number of components to use
            for (int A = 1; A<=ncomp; A++) { 
                // How well did we do with this many components?
                cout << A << " components\t";
                cout << "explained variance: " << plsm.explained_variance(X, Y, A);
                //cout << "root mean squared error of prediction (RMSEP):" << plsm.rmsep(X, Y, A) << endl;
                cout << " SSE: " << plsm.SSE(X,Y,A) <<  endl; 
            }

            cout << "Validation (PRESS):\n";
            cout << plsm.loo_validation(X, Y, PRESS) << endl;

            cout << "Validation (RMSEP):\n";
            cout << plsm.loo_validation(X, Y, RMSEP) << endl;

            cout << "Optimal number of components (LOO):\t" << plsm.optimal_num_components(X, Y, LOO) << endl;

            Row obsmet_orig;
            obsmet_orig.setZero(2);
            obsmet_orig(0) = 40; // number of dice
            obsmet_orig(1) = 10; // number of faces
            Row = obsmet = colwise_z_scores // transform using mean and sd of simulated data

            cout << "scores:\t" << plsm.scores(obsmet) << endl;
            //cout << "Optimal number of components (NEW DATA):\t" << plsm.optimal_num_components(X.bottomRows(nobs - nobs*0.6), Y.bottomRows(nobs - nobs*0.6), NEW_DATA) << endl;
            cout << "DONE!!!\n";
            vector<int> DUMMY;
            return DUMMY; ////////////////FINISH IMPLEMENTING
        }
        
        vector<double> euclidean( Row obsmet, Mat2D simmet ) {
        
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

        Row sample_predictive_priors() { Row dummy; return dummy; }

/*        void sample_predictive_prior(int set_num, vector<Particle> &particles, int sample_size, MTRand* mtrand, int mpi_rank) {
            // Read in predictive prior from the last set.
            // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
            vector<Particle> predictive_prior = read_predictive_prior_file( set_num - 1 );
            Doubled_variances sampling_variance = calculate_doubled_variances( predictive_prior );
            
            vector<double> weights;
            if (set_num == 1) {
                // uniform weights for set 0 predictive prior
                weights.resize(predictive_prior.size(), 1.0/(double) predictive_prior.size());
            } else if ( set_num > 1 ) {
                // weights from set - 2 are needed to calculate weights for set - 1
                weights = calculate_predictive_prior_weights( predictive_prior, set_num );
            }
            if (mpi_rank == 0) write_weights_file( weights, set_num - 1 );

            particles.resize(sample_size);
            for (int i = 0; i < sample_size; i++) {
                // Select a particle index j to use from the predictive prior
                int j = rand_nonuniform_int( weights, mtrand );
                particles[i].R0 = rand_trunc_normal( predictive_prior[ j ].R0, sampling_variance.R0, R0_hard_min, R0_hard_max, mtrand ); 
                particles[i].Ih = rand_trunc_normal( predictive_prior[ j ].Ih, sampling_variance.Ih, Ih_hard_min, Ih_hard_max, mtrand ); 
                particles[i].h  = rand_trunc_normal( predictive_prior[ j ].h,  sampling_variance.h,  h_hard_min,  h_hard_max,  mtrand ); 
                particles[i].P0 = (int) (rand_trunc_normal( predictive_prior[ j ].P0, sampling_variance.P0, P0_hard_min, P0_hard_max, mtrand ) + 0.5);
            }
            return;
        }*/
};

int main(int argc, char* argv[]) {
    gsl_rng_set(RNG, time (NULL) * getpid()); // seed the rng using sys time and the process id
    ModelParameters* par = new ModelParameters();
    par->add_parameter( new Parameter("ndice", UNIFORM, INT, 1, 100) );
    par->add_parameter( new Parameter("sides", UNIFORM, INT, 1, 100) );
    
    ModelMetrics* met = new ModelMetrics(); // the name and order of summary stats returned by simulator
    met->add_metric( "sum", INT );
    met->add_metric( "sd", FLOAT );

    AbcSmc* abc = new AbcSmc(par, met);
    abc->set_smc_iterations(1); // or have it test for convergence
    abc->set_num_samples(1000);
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



/*
struct Parameters {
    string observed_data_file;
    int network_size, burnin, sample_size, threads;
    double epsilon, obs_mean, obs_median, obs_max, obs_range, obs_sd, obs_skew, obs_ss, obs_ll, obs_sl, obs_ls;
};

struct Particle { 
    double R0, Ih, h; 
    int P0;
};

struct Doubled_variances {
    double R0, Ih, P0, h;
};*/

/*
struct ParticleSet { 
    vector<Particle> particles;
    Particle back() { return particles.back(); }
    void resize(const int n) { particles.resize(n); }
    int size() { return (int) particles.size(); }
    Particle& operator[] ( const int nIndex) { return particles[nIndex]; }
};
*/

//vector<ParticleSet> theta; 

// Define ABC parameters
/*const double R0_min = 1.0;
const double R0_max = 8.0;
const double Ih_min = 1.0/12.0;
const double Ih_max = 100.0;
const int P0_min = 1;
const int P0_max = 256;
const double h_min  = 0.0;
const double h_max  = 1.0;

// Don't consider values outside of these bounds
// (because they would be nonsensical or would cause simulator to crash)
const double R0_hard_min = 0.0;
const double R0_hard_max = 1000.0;
const double Ih_hard_min = 0.00001;
const double Ih_hard_max = 10000.0;
const int P0_hard_min = 1;
const int P0_hard_max = 1000;
const double h_hard_min  = 0.0;
const double h_hard_max  = 1.0;
*/

/*void load_observed_data(string filename, char sep, map<string, vector<float> > &data);
Doubled_variances calculate_doubled_variances( vector<Particle> particle_set );
vector<double> calculate_predictive_prior_weights(vector<Particle> predictive_prior, int set_num); 
double uniform_pdf(double a, double b) { return 1.0 / fabs(b-a); }

vector<float> flatten_map(map<string, vector<float> > data) {
    vector<float> flat;
    map<string, vector<float> >::iterator it;

    for ( it = data.begin(); it != data.end(); it++ ) {
        string loc = it->first;
        for (unsigned int i = 0; i < data[loc].size(); i++) {
            flat.push_back(data[loc][i]);
        }
    }
    return flat;
}

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


int determine_set_number() {
    bool found_file = true;
    int set_num = 0;
    while (found_file) {
        stringstream ss;
        ss << "predictive_prior." << set_num;
        if ( fileExists(ss.str()) ) {
            set_num++;
        } else {
            found_file = false;
        }
    }
    return set_num;
}


double rand_trunc_normal(double mu, double sigma_squared, double min, double max, MTRand* mtrand) {
    double sigma = sqrt(sigma_squared);
    // Don't like this, but it should work
    while (1) {
        double dev = mtrand->randNorm(mu, sigma);
        if (dev >= min and dev <= max) {
            return dev;
        }
    }
}
                                

void sample_predictive_prior(int set_num, vector<Particle> &particles, int sample_size, MTRand* mtrand, int mpi_rank) {
    // Read in predictive prior from the last set.
    // We need to calculate the proper weights for the predictive prior so that we know how to sample from it.
    vector<Particle> predictive_prior = read_predictive_prior_file( set_num - 1 );
    Doubled_variances sampling_variance = calculate_doubled_variances( predictive_prior );
    
    vector<double> weights;
    if (set_num == 1) {
        // uniform weights for set 0 predictive prior
        weights.resize(predictive_prior.size(), 1.0/(double) predictive_prior.size());
    } else if ( set_num > 1 ) {
        // weights from set - 2 are needed to calculate weights for set - 1
        weights = calculate_predictive_prior_weights( predictive_prior, set_num );
    }
    if (mpi_rank == 0) write_weights_file( weights, set_num - 1 );

    particles.resize(sample_size);
    for (int i = 0; i < sample_size; i++) {
        // Select a particle index j to use from the predictive prior
        int j = rand_nonuniform_int( weights, mtrand );
        particles[i].R0 = rand_trunc_normal( predictive_prior[ j ].R0, sampling_variance.R0, R0_hard_min, R0_hard_max, mtrand ); 
        particles[i].Ih = rand_trunc_normal( predictive_prior[ j ].Ih, sampling_variance.Ih, Ih_hard_min, Ih_hard_max, mtrand ); 
        particles[i].h  = rand_trunc_normal( predictive_prior[ j ].h,  sampling_variance.h,  h_hard_min,  h_hard_max,  mtrand ); 
        particles[i].P0 = (int) (rand_trunc_normal( predictive_prior[ j ].P0, sampling_variance.P0, P0_hard_min, P0_hard_max, mtrand ) + 0.5);
    }
    return;
}


Doubled_variances calculate_doubled_variances( vector<Particle> particle_set ) {
    int set_size = particle_set.size();
    vector<double> R0_vec( set_size );
    vector<double> Ih_vec( set_size );
    vector<double> h_vec( set_size );
    vector<double> P0_vec( set_size );

    for (int i = 0; i < set_size; i++) {
        R0_vec[i] = particle_set[i].R0;
        Ih_vec[i] = particle_set[i].Ih;
        h_vec[i]  = particle_set[i].h;
        P0_vec[i] = particle_set[i].P0;
    }

    Doubled_variances doubled_var;
    doubled_var.R0 = 2 * variance( R0_vec );
    doubled_var.Ih = 2 * variance( Ih_vec );
    doubled_var.P0 = 2 * variance( P0_vec );
    doubled_var.h  = 2 * variance( h_vec );

    return doubled_var;
}


vector<double> calculate_predictive_prior_weights(vector<Particle> predictive_prior, int set_num) {

    const int prior_size = predictive_prior.size();
    vector<double> predictive_prior_weights( prior_size );
    const double numerator = uniform_pdf(R0_min, R0_max)
                           * uniform_pdf(Ih_min, Ih_max)
                           * uniform_pdf(P0_min, P0_max)
                           * uniform_pdf(h_min, h_max);

    vector <double>   old_predictive_prior_weights = read_weights_file( set_num - 2 ); 
    vector<Particle>  old_predictive_prior         = read_predictive_prior_file( set_num - 2 );
    const int old_prior_size                       = old_predictive_prior.size();
    Doubled_variances old_par_2var                 = calculate_doubled_variances( old_predictive_prior );


    for (int i = 0; i < prior_size; i++) {
        double denominator = 0.0;
        for (int j = 0; j < old_prior_size; j++) {
            denominator += old_predictive_prior_weights[j]
                           //double normal_pdf(double x, double mu, double var) {
                           * normal_pdf(predictive_prior[i].R0, old_predictive_prior[j].R0, old_par_2var.R0)
                           * normal_pdf(predictive_prior[i].Ih, old_predictive_prior[j].Ih, old_par_2var.Ih)
                           * normal_pdf(predictive_prior[i].P0, old_predictive_prior[j].P0, old_par_2var.P0)
                           * normal_pdf(predictive_prior[i].h,  old_predictive_prior[j].h,  old_par_2var.h);
        }
        predictive_prior_weights[i] = numerator / denominator;
    }

    return normalize_dist( predictive_prior_weights );
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
