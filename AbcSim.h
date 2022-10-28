#ifndef ABCSIM_H
#define ABCSIM_H

#include "AbcUtil.h"
#include <vector>
#include <functional>

using std::vector;

using SimFunc = function<vector<ABC::float_type>(vector<ABC::float_type>, const unsigned long int, const unsigned long int, const ABC::MPI_par*)>;

struct WrapExec : SimFunc {
    string command;
    vector<ABC::float_type> operator()(
      vector<ABC::float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) {
        auto execcom = command;
        bool particle_success;
        for (auto par : pars) { execcom += " " + ABC::toString(par); }
        auto retval = ABC::exec(command);
        if (retval == "ERROR" or retval == "") {
            cerr << command << " does not exist or appears to be an invalid simulator on MPI rank " << _mp->mpi_rank << endl;
            particle_success = false;
        }
        stringstream ss;
        ss.str(retval);
        vector<ABC::float_type> mets;
// TODO deal with empty mets on !particle_success
        if (particle_success) {
            ABC::float_type met;
            while(ss >> met) mets.push_back(met);
        }

        return mets;
    }
};

#endif