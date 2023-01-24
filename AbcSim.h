#ifndef ABCSIM_H
#define ABCSIM_H

#include "AbcUtil.h"
#include <vector>
#include <functional>
#include <dlfcn.h>

using std::vector;

struct AbcSimFun {
    virtual vector<ABC::float_type> operator()(
      vector<ABC::float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const = 0;
};

struct AbcSimUnset : AbcSimFun {
    vector<ABC::float_type> operator()(
      vector<ABC::float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const {
        std::cerr << "ERROR: A pointer to a simulator function (prefered) or an external simulator executable must be defined." << std::endl;
        exit(100);
    }
};

typedef vector<ABC::float_type> AbcSimF(vector<ABC::float_type>, const unsigned long int, const unsigned long int, const ABC::MPI_par*);

inline AbcSimF * loadSO(const char * target) {
    void* handle = dlopen(target, RTLD_LAZY);
    if (!handle) {
        std::cerr << "Failed to open simulator object: " << target << " ; " << dlerror() << std::endl;
        exit(101);
    }
    auto simf = (AbcSimF*)dlsym(handle, "simulator");
    if(!simf) {
        std::cerr << "Failed to find 'simulator' function in " << target << " ; " << dlerror() << std::endl;
        dlclose(handle);
        exit(102);
    }
    return simf;
};

struct AbcFPtr : AbcSimFun {
    AbcSimF* fptr;
    AbcFPtr(AbcSimF * _fptr) : fptr(_fptr) { }
    AbcFPtr(const char * target) : fptr(loadSO(target)) { }
    vector<ABC::float_type> operator()(
      vector<ABC::float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const {
        return fptr(pars, seed, serial, _mp);
    }
};

struct AbcExec : AbcSimFun {
    const string command;
    AbcExec(string _command) : command(_command) { }

    vector<ABC::float_type> operator()(
      vector<ABC::float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const {
        auto execcom = command;
        bool particle_success;
        for (auto par : pars) { execcom += " " + ABC::toString(par); }
        auto retval = ABC::exec(command);
        if (retval == "ERROR" or retval == "") {
            std::cerr << command << " does not exist or appears to be an invalid simulator on MPI rank " << _mp->mpi_rank << std::endl;
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