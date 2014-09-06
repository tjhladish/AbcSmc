#include "AbcSmc.h"
#include "pls.h"
#include "RunningStat.h"

#include <iostream>
#include <cstdio>
#include <cstring>
#include <json/json.h>
#include <string>
#include <sstream>
#include <fstream>

// need a positive int that is very unlikely 
// to be less than the number of particles
#define STOP_TAG 10000000

using std::vector;
using std::string;
using std::stringstream;
using std::ofstream;
using std::setw;

const string double_bar = "=========================================================================================";
const int PREC = 5;
const int WIDTH = 12; 

bool AbcSmc::parse_config(string conf_filename) {
    Json::Value par;   // will contains the par value after parsing.
    Json::Reader reader;
    string json_data = slurp(conf_filename);

    bool parsingSuccessful = reader.parse( json_data, par );
    if ( !parsingSuccessful ) {
        // report to the user the failure and their locations in the document.
        std::cerr  << "Failed to parse configuration\n"
            << reader.getFormattedErrorMessages();
        return false;
    }

    // Parse model parameters
    const Json::Value model_par = par["parameters"];
    for ( unsigned int i = 0; i < model_par.size(); ++i )  {// Iterates over the sequence elements.
        //cerr << model_par[i] << endl;

        string name = model_par[i]["name"].asString();
        string short_name = model_par[i].get("short_name", "").asString();

        PriorType ptype = UNIFORM; 
        string ptype_str = model_par[i]["dist_type"].asString();
        if (ptype_str == "UNIFORM") {
            ptype = UNIFORM;
        } else if (ptype_str == "NORMAL" or ptype_str == "GAUSSIAN") {
            ptype = NORMAL;
        } else if (ptype_str == "PSEUDO") {
            ptype = PSEUDO;
        } else {
            cerr << "Unknown parameter distribution type: " << ptype_str << ".  Aborting." << endl;
            exit(-205);
        }

        NumericType ntype = INT; 
        string ntype_str = model_par[i]["num_type"].asString();
        if (ntype_str == "INT") {
            ntype = INT;
        } else if (ntype_str == "FLOAT") {
            ntype = FLOAT;
        } else {
            cerr << "Unknown parameter numeric type: " << ntype_str << ".  Aborting." << endl;
            exit(-206);
        }


        double par1 = model_par[i]["par1"].asDouble();
        double par2 = model_par[i]["par2"].asDouble();

        add_next_parameter( name, short_name, ptype, ntype, par1, par2);
    }

    // Parse model metrics 
    const Json::Value model_met = par["metrics"];
    for ( unsigned int i = 0; i < model_met.size(); ++i )  {// Iterates over the sequence elements.
        //cerr << model_met[i] << endl;

        string name = model_met[i]["name"].asString();
        string short_name = model_met[i].get("short_name", "").asString();

        NumericType ntype = INT; 
        string ntype_str = model_met[i]["num_type"].asString();
        if (ntype_str == "INT") {
            ntype = INT;
        } else if (ntype_str == "FLOAT") {
            ntype = FLOAT;
        } else {
            cerr << "Unknown metric numeric type: " << ntype_str << ".  Aborting." << endl;
            exit(-207);
        }


        double val = model_met[i]["value"].asDouble();

        add_next_metric( name, short_name, ntype, val);
    }

    string executable = par.get("executable", "").asString();
    if (executable != "") { set_executable( executable ); }

    string resume_dir = par.get("resume_directory", "").asString();
    if (resume_dir != "") {
        if (_mp->mpi_rank == mpi_root) cerr << "Resuming in directory: " << resume_dir << endl;
        set_resume_directory( resume_dir );
        set_resume( true );
    }

    set_smc_iterations( par["smc_iterations"].asInt() ); // or have it test for convergence
    set_num_samples( par["num_samples"].asInt() );
    set_predictive_prior_fraction( par["predictive_prior_fraction"].asFloat() );
    set_pls_validation_training_fraction( par["pls_training_fraction"].asFloat() ); // fraction of runs to use for training
    set_particle_basefilename( par["particle_basefilename"].asString() );
    set_predictive_prior_basefilename( par["predictive_prior_basefilename"].asString() );
    
    return true;

}


void AbcSmc::write_particle_file( const int t ) {
    ofstream particle_file;
    string pad = t < 10 ? "0" : "";
    particle_file.open( (_particle_filename + "." + pad + toString(t)).c_str() );
    if (particle_file.fail()) {
       cerr << "ERROR: Particle file '" << _particle_filename << "' cannot be open for writing." << endl;
       exit(-1);
    }

    // write header
    particle_file << "iteration sample ";
    for (int i = 0; i<npar(); i++) { particle_file << _model_pars[i]->get_short_name() << " "; }
    for (int i = 0; i<nmet(); i++) { particle_file << _model_mets[i]->get_short_name() << " "; } 
    particle_file << endl;

    for (int i = 0; i < _num_particles; i++) {
        particle_file << t << " " << i;
        for (int j = 0; j < npar(); j++) { particle_file << " " << _particle_parameters[t](i,j); }
        for (int j = 0; j < nmet(); j++) { particle_file << " " << _particle_metrics[t](i,j); }
        particle_file << endl;
    }
    particle_file.close();
}


void AbcSmc::write_predictive_prior_file( const int t ) {
    ofstream predictive_prior_file;
    string pad = t < 10 ? "0" : "";
    predictive_prior_file.open( (_predictive_prior_filename + "." + pad + toString(t)).c_str() );
    if (predictive_prior_file.fail()) {
       cerr << "ERROR: predictive prior file '" << _predictive_prior_filename << "' cannot be open for writing." << endl;
       exit(-1);
    }

    // write header
    predictive_prior_file << "iteration rank sample ";
    for (int i = 0; i<npar(); i++) { predictive_prior_file << _model_pars[i]->get_short_name() << " "; }
    for (int i = 0; i<nmet(); i++) { predictive_prior_file << _model_mets[i]->get_short_name() << " "; } 
    predictive_prior_file << endl;

    for (int rank = 0; rank < _predictive_prior_size; rank++) {
        const int idx = _predictive_prior[t][rank];
        predictive_prior_file << t << " " << rank << " " << idx;
        for (int j = 0; j < npar(); j++) { predictive_prior_file << " " << _particle_parameters[t](idx, j); } 
        for (int j = 0; j < nmet(); j++) { predictive_prior_file << " " << _particle_metrics[t](idx, j); } 
        predictive_prior_file << endl;
    }
    predictive_prior_file.close();
}


void AbcSmc::run(string executable, const gsl_rng* RNG) {
    // NEED TO ADD SANITY CHECKS HERE
    _executable_filename = executable;
    _particle_parameters.clear();
    _particle_metrics.clear();
    _weights.clear();

    cerr << std::setprecision(PREC);

    if (_mp->mpi_rank == 0) {
        _particle_parameters.resize( _num_smc_sets, Mat2D::Zero(_num_particles, npar()) );
        _particle_metrics.resize( _num_smc_sets, Mat2D::Zero(_num_particles, nmet()) );
        _weights.resize(_num_smc_sets);
    }

    _predictive_prior.clear();
    _predictive_prior.resize(_num_smc_sets);

    for (int t = 0; t<_num_smc_sets; t++) {

#ifdef USING_MPI
        assert(_mp->mpi_size > 1);
        bool success = _populate_particles_mpi( t, _particle_metrics[t], _particle_parameters[t], RNG );
#else
        bool success = _populate_particles( t, _particle_metrics[t], _particle_parameters[t], RNG );
#endif

        if (!success) { cerr << "MPI rank" << _mp->mpi_rank << " failed while generating particles." << endl;}
        if (_mp->mpi_rank == mpi_root) {
            write_particle_file(t);

            // Select best scoring particles
            cerr << double_bar << endl << "Set " << t << endl << double_bar << endl;
            if (resume()) {
                read_predictive_prior(t);
            } else {
                _filter_particles( t, _particle_metrics[t], _particle_parameters[t] );
            }

            calculate_predictive_prior_weights( t );

            report_convergence_data(t);
            cerr << endl << endl;

            write_predictive_prior_file(t);
        }
    }
}
void print_stats(ostream& stream, string str1, string str2, double val1, double val2, double delta, double pct_chg, string tail) {
    stream << "    " + str1 + ", " + str2 + "  ( delta, % ): "  << setw(WIDTH) << val1 << ", " << setw(WIDTH) << val2 
                                                      << " ( " << setw(WIDTH) << delta << ", " << setw(WIDTH) << pct_chg  << "% )\n" + tail;
}

void AbcSmc::report_convergence_data(int t) {
    vector<double> last_means( npar(), 0 );
    vector<double> current_means( npar(), 0 );

    for (int j = 0; j < npar(); j++) {
        int N = _predictive_prior[t].size();
        for (int i = 0; i < N; i++) {
            double particle_idx = _predictive_prior[t][i];
            double par_value = _particle_parameters[t](particle_idx, j);
            current_means[j] += par_value;
        }
        current_means[j] /= N;
        
        if (t > 0) {
            int N2 = _predictive_prior[t-1].size();
            for (int i = 0; i < N2; i++) {
                double particle_idx = _predictive_prior[t-1][i];
                double par_value = _particle_parameters[t-1](particle_idx, j);
                last_means[j] += par_value;
            }
            last_means[j] /= N2;
        }
    }

    cerr << double_bar << endl;
    if (t == 0) {
        cerr << "Predictive prior summary statistics:\n";
    } else {
        cerr << "Convergence data for predictive priors:\n";
    }
    for (unsigned int i = 0; i < _model_pars.size(); i++) {
        const Parameter* par = _model_pars[i];
        const double current_stdev = sqrt(par->get_doubled_variance(t)/2.0);
        const double prior_mean = par->get_prior_mean();
        const double prior_mean_delta = current_means[i] - prior_mean;
        const double prior_mean_pct_chg = prior_mean != 0 ? 100 * prior_mean_delta / prior_mean : INFINITY;

        const double prior_stdev = par->get_prior_stdev();
        const double prior_stdev_delta = current_stdev - prior_stdev;
        const double prior_stdev_pct_chg = prior_stdev != 0 ? 100 * prior_stdev_delta / prior_stdev : INFINITY;
        if (t == 0) {
            cerr << "  Par " << i << ": \"" << par->get_name() << "\"\n";

            cerr << "  Means:\n";
            print_stats(cerr, "Prior", "current", prior_mean, current_means[i], prior_mean_delta, prior_mean_pct_chg, "");
            cerr << "  Standard deviations:\n";
            print_stats(cerr, "Prior", "current", prior_stdev, current_stdev, prior_stdev_delta, prior_stdev_pct_chg, "\n");
        } else {
            double last_stdev = sqrt(_model_pars[i]->get_doubled_variance(t-1)/2.0);
            double delta, pct_chg;

            cerr << "  Par " << i << ": \"" << _model_pars[i]->get_name() << "\"\n";

            delta = current_means[i] - last_means[i];
            pct_chg = last_means[i] != 0 ? 100 * delta / last_means[i] : INFINITY;
            cerr << "  Means:\n";
            print_stats(cerr, "Prior", "current", prior_mean, current_means[i], prior_mean_delta, prior_mean_pct_chg, "");
            print_stats(cerr, "Last", " current", last_means[i], current_means[i], delta, pct_chg, "\n");

            delta = current_stdev - last_stdev;
            pct_chg = last_stdev != 0 ? 100 * delta / last_stdev : INFINITY;
            cerr << "  Standard deviations:\n";
            print_stats(cerr, "Prior", "current", prior_stdev, current_stdev, prior_stdev_delta, prior_stdev_pct_chg, "");
            print_stats(cerr, "Last", " current", last_stdev, current_stdev, delta, pct_chg, "\n");
        }
    }
}

bool AbcSmc::read_particle_set(int t, Mat2D &X_orig, Mat2D &Y_orig ) {
    string pad = t < 10 ? "0" : "";
    string existing_particle_filename = _particle_filename + "." + pad + toString(t);

    ifstream iss(existing_particle_filename.c_str());
    if (!iss) {
        cerr << "WARNING: " << existing_particle_filename << " not found." << endl;
        cerr << "         This is okay if all existing particle files have been read." << endl;
        set_resume( false );
        return false;
    }

    string buffer;
    vector<vector<float_type> > fields;
    int line_num=0;
    while(std::getline(iss,buffer)){
        line_num++;
        vector<string>line = split(buffer,' ');
        // not a header and #fields is correct
        // first two columns are set num (t) and iteration
        const int offset = 2;
        if( line_num > 1 && (signed) line.size() == npar() + nmet() + offset ){
            vector<float_type> row( npar() + nmet() );
            for(unsigned int i=0; i < row.size(); i++){
                row[i] = string2float_type(line[i+offset]);
            }
            fields.push_back(row);
        }
    }
    iss.close();

    if(fields.size()  == 0 ) {
        cerr << "ERROR: " << existing_particle_filename << " missing data." << endl;
        set_resume( false );
        return false;
    }

cerr << "fields size, nmet, npar: " << fields.size() << " " << nmet() << " " << npar() << endl;
    Y_orig.resize( fields.size(), npar() );
    X_orig.resize( fields.size(), nmet() );
    for ( unsigned int i=0; i < fields.size(); i++ ) {
        for ( int j=0; j < npar(); j++ ) Y_orig(i,j)=fields[i][j];
        for ( int j=npar(); j < npar()+nmet(); j++ ) X_orig(i,j-npar())=fields[i][j];
    }

    cerr << "Loaded particle set " << t << endl;
    return true;
}


bool AbcSmc::read_predictive_prior( int t ) {
    vector<int> ranking;

    string pad = t < 10 ? "0" : "";
    string existing_predictive_prior_filename = _predictive_prior_filename + "." + pad + toString(t);

    ifstream iss(existing_predictive_prior_filename.c_str());
    if (!iss) {
        cerr << "WARNING: " << existing_predictive_prior_filename << " not found." << endl;
        cerr << "         It appears a particle file exists for which there is no predictive prior file." << endl;
        set_resume( false );
        return false;
    }

    //iteration rank sample ndice sides sum sd
    //0 0 953 129 7 525 1.91697
    //0 1 503 3 10 24 2
    string buffer;
    vector<vector<float_type> > fields;
    int line_num=0;
    while(std::getline(iss,buffer)){
        line_num++;
        vector<string>line = split(buffer,' ');
        // not a header and #fields is correct
        // first two columns are set num (t) and iteration
        const int offset = 3;
        const int sample_idx = 2;
        if( line_num > 1 && (signed) line.size() == npar() + nmet() + offset ){
            ranking.push_back( string2int(line[sample_idx]) );
        }
    }
    iss.close();

    if ( ranking.size() == 0 ) {
        cerr << "ERROR: " << existing_predictive_prior_filename << " missing data." << endl;
        set_resume( false );
        return false;
    }

    _predictive_prior[t] = ranking;
    cerr << "Loaded predictive prior set " << t << endl;
    return true;
}


bool AbcSmc::_populate_particles_mpi(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {
    // this is silly -- only 'true' is ever returned from this function
    int continue_flag = (int) true;
    if (_mp->mpi_rank == mpi_root and resume() and read_particle_set( t, X_orig, Y_orig )) {
        continue_flag = (int) false;
    }
#ifdef USING_MPI
    MPI_Bcast(&continue_flag, 1, MPI_INT, mpi_root, _mp->comm);
#endif
    if (not (bool) continue_flag) return true; // resume was set, reading particles succeeded

    if (_mp->mpi_rank == mpi_root) {
        _particle_scheduler(t, X_orig, Y_orig, RNG);
    } else {
        _particle_worker();
    }

    return true;
}


void AbcSmc::_particle_scheduler(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {
#ifdef USING_MPI
    MPI_Status status;
    long double *send_data = new long double[npar() * _num_particles](); // all params, flattened array
    long double *rec_data  = new long double[nmet() * _num_particles](); // all metrics, flattened array

    // sample parameter distributions; copy values into Y matrix and into send_data buffer
    Row par_row;
    for (int i = 0; i<_num_particles; i++) {
        if (t == 0) { // sample priors
            par_row = sample_priors(RNG);
            if (i<_num_particles) Y_orig.row(i) = par_row;
        } else {      // sample predictive priors
            par_row = sample_predictive_priors(t, RNG);
            if (i<_num_particles) Y_orig.row(i) = par_row;
        }
        for (int j = 0; j<npar(); j++) { send_data[i*npar() + j] = par_row(j); }
    }
    // Seed the workers with the first 'num_workers' jobs
    int particle_id = 0;
    for (int rank = 1; rank < _mp->mpi_size; ++rank) {
        particle_id = rank - 1;                   // which row in Y
        MPI_Send(&send_data[particle_id*npar()],  // message buffer
                 npar(),                          // number of elements
                 MPI_LONG_DOUBLE,                 // data item is a double
                 rank,                            // destination process rank
                 particle_id,                     // message tag
                 _mp->comm);                      // always use this
    }
    
    // Receive a result from any worker and dispatch a new work request
    long double rec_buffer[nmet()];
    particle_id++; // move cursor to next particle to be sent
    while ( particle_id < _num_particles ) { 
        MPI_Recv(&rec_buffer,                     // message buffer
                 nmet(),                          // message size
                 MPI_LONG_DOUBLE,                 // of type double
                 MPI_ANY_SOURCE,                  // receive from any sender
                 MPI_ANY_TAG,                     // any type of message
                 _mp->comm,                       // always use this
                 &status);                        // received message info
        for (int m = 0; m<nmet(); ++m) rec_data[status.MPI_TAG*nmet() + m] = rec_buffer[m];

        MPI_Send(&send_data[particle_id*npar()],  // message buffer
                 npar(),                          // number of elements
                 MPI_LONG_DOUBLE,                 // data item is a double
                 status.MPI_SOURCE,               // send it to the rank that just finished
                 particle_id,                     // message tag
                 _mp->comm);                      // always use this
        particle_id++; // move cursor
    }
    
    // receive results for outstanding work requests--there are exactly 'num_workers' left
    for (int rank = 1; rank < _mp->mpi_size; ++rank) {
        MPI_Recv(&rec_buffer, 
                 nmet(), 
                 MPI_LONG_DOUBLE, 
                 MPI_ANY_SOURCE,
                 MPI_ANY_TAG, 
                 _mp->comm, 
                 &status);
        for (int m = 0; m<nmet(); ++m) rec_data[status.MPI_TAG*nmet() + m] = rec_buffer[m];
    }
     
    // Tell all the workers they're done for now
    for (int rank = 1; rank < _mp->mpi_size; ++rank) {
        MPI_Send(0, 0, MPI_INT, rank, STOP_TAG, _mp->comm);
    }

    vector<int> bad_particle_idx; // bandaid, in case simulator returns nonsense values
    for (int i = 0; i<_num_particles; i++) {
        for (int j = 0; j<nmet(); j++) {
            const double met_val = rec_data[i*nmet() + j];
            if (not isfinite(met_val)) bad_particle_idx.push_back(i);
            X_orig(i,j) = rec_data[i*nmet() + j];
        }
    }
    
    // bandaid, in case simulator returns nonsense values
    for (unsigned int i = 0; i < bad_particle_idx.size(); i++) {
        X_orig.row(i).fill(numeric_limits<float_type>::min());
        Y_orig.row(i).fill(numeric_limits<float_type>::min());
    }

    delete[] send_data;
    delete[] rec_data;
#endif
}


void AbcSmc::_particle_worker() {
#ifdef USING_MPI
    MPI_Status status;
    long double *local_Y_data = new long double[npar()]();
    long double *local_X_data = new long double[nmet()]();

    while (1) {
        MPI_Recv(local_Y_data, 
                 npar(), 
                 MPI_LONG_DOUBLE, 
                 mpi_root, 
                 MPI_ANY_TAG,
                 _mp->comm, 
                 &status);

        // Check the tag of the received message.
        if (status.MPI_TAG == STOP_TAG) {
            delete[] local_Y_data;
            delete[] local_X_data;
            return;
        }

        Row pars(npar());
        Row mets(nmet());

        for (int j = 0; j<npar(); j++) { pars[j] = local_Y_data[j]; }

        _run_simulator(pars, mets);

        for (int j = 0; j<nmet(); j++) { local_X_data[j] = mets[j]; }

        MPI_Send(local_X_data, 
                 nmet(), 
                 MPI_LONG_DOUBLE, 
                 mpi_root, 
                 status.MPI_TAG, 
                 _mp->comm);
    }
#endif
}


bool AbcSmc::_run_simulator(Row &par, Row &met) {
    bool particle_success = true;
    if (use_simulator) {
        vector<float_type> met_vec = _simulator( as_vector(par), _mp );
        met = as_row(met_vec);
    } else if (use_executable) {
        string command = _executable_filename;
        for (int j = 0; j<npar(); j++) { command += " " + toString(par[j]); }

        string retval = exec(command);
        if (retval == "ERROR" or retval == "") {
            cerr << command << " does not exist or appears to be an invalid simulator on MPI rank " << _mp->mpi_rank << endl;
            particle_success = false;
        }
        stringstream ss;
        ss.str(retval);

        for (int j = 0; j<nmet(); j++) {
            if (particle_success) {
                ss >> met[j];
            } else {
                met[j] = numeric_limits<float_type>::min();
            }
        }
    } else {
        cerr << "ERROR: A pointer to a simulator function (prefered) or an external simulator executable must be defined.\n";
        particle_success = false;
    }
    return particle_success;
}


bool AbcSmc::_populate_particles(int t, Mat2D &X_orig, Mat2D &Y_orig, const gsl_rng* RNG) {
    if (resume() and read_particle_set( t, X_orig, Y_orig )) {
        return true;
    }

    for (int i = 0; i<_num_particles; i++) {
        if (t == 0) { // sample priors
            Y_orig.row(i) = sample_priors(RNG);
        } else { // sample predictive priors
            Y_orig.row(i) = sample_predictive_priors(t, RNG);
        }
       
        Row pars = Y_orig.row(i);
        Row mets = X_orig.row(i);
        _run_simulator(pars, mets);
        X_orig.row(i) = mets;
    }
    return true;
}


void AbcSmc::_filter_particles (int t, Mat2D &X_orig, Mat2D &Y_orig) {
    // Run PLS
    // Box-Cox transform data -- TODO?
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
    plsm.plsr(X.topRows(_pls_training_set_size), Y.topRows(_pls_training_set_size), KERNEL_TYPE1);

    // A is number of components to use
    for (int A = 1; A<=ncomp; A++) { 
        // How well did we do with this many components?
        cerr << setw(2) << A << " components ";
        cerr << "explained variance: " << plsm.explained_variance(X, Y, A);
        //cerr << "root mean squared error of prediction (RMSEP):" << plsm.rmsep(X, Y, A) << endl;
        cerr << " SSE: " << plsm.SSE(X,Y,A) <<  endl; 
    }

    const int test_set_size = nobs - _pls_training_set_size;
    Rowi num_components = plsm.optimal_num_components(X.bottomRows(test_set_size), Y.bottomRows(test_set_size), NEW_DATA);
    int num_components_used = num_components.maxCoeff();
    cerr << "Optimal number of components for each parameter (validation method == NEW DATA):\t" << num_components << endl;
    cerr << "Using " << num_components_used << " components." << endl;

    // Calculate new, orthogonal metrics (==scores) using the pls model
    // Is casting as real always safe?
    Row   obs_scores = plsm.scores(obs_met, num_components_used).row(0).real();
    Mat2D sim_scores = plsm.scores(X, num_components_used).real();
    Col   distances  = euclidean(obs_scores, sim_scores);
    vector<int> ranking = ordered(distances);

    vector<int>::iterator first = ranking.begin();
    vector<int>::iterator last  = ranking.begin() + _predictive_prior_size;
    vector<int> sample(first, last);
    _predictive_prior[t] = sample;

    cerr << "Best five:\n";
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << _model_pars[i]->get_short_name(); } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << _model_mets[i]->get_short_name(); } cerr << endl;
    for (int q=0; q<5; q++) {
        const int idx = ranking[q];
        for (int i = 0; i < Y_orig.cols(); i++) { cerr << setw(WIDTH) << Y_orig(idx, i); } 
        cerr << " | ";
        for (int i = 0; i < X_orig.cols(); i++) { cerr << setw(WIDTH) << X_orig(idx, i); } 
        cerr << endl;
    }

    cerr << "Worst five:\n";
    for (int i = 0; i<npar(); i++) { cerr << setw(WIDTH) << _model_pars[i]->get_short_name(); } cerr << " | ";
    for (int i = 0; i<nmet(); i++) { cerr << setw(WIDTH) << _model_mets[i]->get_short_name(); } cerr << endl;
    for (unsigned int q=ranking.size()-1; q>=ranking.size()-5; q--) {
        const int idx = ranking[q];
        for (int i = 0; i < Y_orig.cols(); i++) { cerr << setw(WIDTH) << Y_orig(idx, i); } 
        cerr << " | ";
        for (int i = 0; i < X_orig.cols(); i++) { cerr << setw(WIDTH) << X_orig(idx, i); } 
        cerr << endl;
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
    bool increment_nonrandom_par = true; // only one PSEUDO parameter gets incremented each time
    // for each parameter
    for (unsigned int i = 0; i < _model_pars.size(); i++) {
        Parameter* p = _model_pars[i];
        float_type val;
        // if it's a non-random, PSEUDO parameter
        if (p->get_prior_type() == PSEUDO) {
            // get the state now, so that the first time it will have the initialized value
            val = (float_type) p->get_state();
            // We need to imitate the way nested loops work, but in a single loop.
            // PSEUDO parameters only get incremented when any and all previous PSEUDO parameters
            // have reached their max values and are being reset
            if (increment_nonrandom_par) {
                // This parameter has reached it's max value and gets reset to minimum
                if (p->get_state() == (int) p->get_prior_max()) {
                    p->reset_state();
                // otherwise, increment this one and prevent others from being incremented
                } else {
                    p->increment_state();
                    increment_nonrandom_par = false;
                }
            }
        } else {
            // Random parameters get sampled independently from each other, and are therefore easy
            val = p->sample(RNG);
        }

        par_sample(i) = val;
    }
//cerr << "sample: " << par_sample << endl;
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
    cerr << "predictive prior size (n): " << _predictive_prior[set_num].size() << endl;
    
    if (set_num == 0) {
        // uniform weights for set 0 predictive prior
        _weights[set_num].resize(_predictive_prior[set_num].size(), 1.0/(double) _predictive_prior[set_num].size());
    } else if ( set_num > 0 ) {
        // weights from set - 1 are needed to calculate weights for current set

        _weights[set_num].resize( _predictive_prior[set_num].size() );
        
        double numerator = 1;
        for (int j = 0; j<npar(); j++) {
            Parameter* par = _model_pars[j];
            numerator *= uniform_pdf(par->get_prior_min(), par->get_prior_max());
        }

        for (unsigned int i = 0; i < _predictive_prior[set_num].size(); i++) {
            double denominator = 0.0;
            for (unsigned int k = 0; k < _predictive_prior[set_num - 1].size(); k++) {
                double running_product = _weights[set_num - 1][k];
                for (int j = 0; j < npar(); j++) {
                    double par_value = _particle_parameters[set_num](_predictive_prior[set_num][i], j);
                    double old_par_value = _particle_parameters[set_num-1](_predictive_prior[set_num-1][k], j);
                    double old_doubled_variance = _model_pars[j]->get_doubled_variance(set_num-1);

                    // This conditional handles the (often improbable) case where a parameter has completely converged.
                    // It allows ABC to continue exploring other parameters, rather than causing the math
                    // to fall apart because the density at the converged value is infinite.
                    if (old_doubled_variance != 0 or par_value != old_par_value) {
                        running_product *= gsl_ran_gaussian_pdf(par_value-old_par_value, sqrt(old_doubled_variance) );
                    }
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
        const Parameter* parameter = _model_pars[j];
        double doubled_variance = parameter->get_doubled_variance(set_num-1);
        double par_min = parameter->get_prior_min();
        double par_max = parameter->get_prior_max();
        par_values(j) = rand_trunc_normal( par_value, doubled_variance, par_min, par_max, RNG );

        if (parameter->get_numeric_type() == INT) {
            par_values(j) = (double) ((int) (par_values(j) + 0.5));
        }
    }
    return par_values;
}

Row AbcSmc::_z_transform_observed_metrics(Row& means, Row& stdevs) {
    Row zmat = Row::Zero(nmet());
    for (int i = 0; i<nmet(); i++) { zmat(i) = (_model_mets[i]->get_obs_val() - means(i)) / stdevs(i); }
    return zmat;
}

