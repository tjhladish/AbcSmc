#ifndef ABCSIM_H
#define ABCSIM_H

#include <vector> // container for receiving parameters / returning metrics
#include <dlfcn.h> // for dynamic version
#include <fstream> // for external executable version
#include <sstream> // for stringstream
#include <string> // for string
#include <iostream> // for cerr

#include <PLS/types.h> // for float_type
#include <AbcSmc/AbcMPIPar.h>

// TODO approach this way?
// https://stackoverflow.com/questions/58394556/c-concepts-can-i-have-a-constraint-requiring-a-function-be-present-in-a-clas

using std::vector;
using std::ostringstream;
using std::istringstream;
using std::string;

// if compiling with MPI support, have to handle slightly complex MPI objects
// this must be defined in a matching way between and AbcSmc executables *and*
// other code using this header.

// defines the core abstraction for a simulator: (1) a functor, (2) with an operator(), (3) return type vector<float> (the metrics),
// (4) arguments vector<float> (the parameters), and unsigned long int, unsigned long int, MPI_par* (the abc seed/serial/mpi information)
// implemented as a pure abstract class - i.e. must be extended by concrete implementations
// In general, the goal is that those implementations should *not* be "stateful". That is, they should not have any internal state which
// is changed when they are used.
struct AbcSimFun {
    virtual vector<float_type> operator()(
        vector<float_type> pars, const unsigned long int seed, const unsigned long int serial
    ) const = 0;

    // the default implementation is to ignore the MPI_par* argument.
    virtual vector<float_type> operator()(
        vector<float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* /* _mp */
    ) const {
        return (*this)(pars, seed, serial);
    };
};

// an AbcSimFun which throws an error if used. This is intended to be used as a default, so that if *not* replaced, the error
// will be thrown when the simulator is used.
struct AbcSimUnset : AbcSimFun {
    vector<float_type> operator()(
      vector<float_type> /*pars*/, const unsigned long int /*seed*/, const unsigned long int /*serial*/
    ) const override {
        std::cerr << "ERROR: A pointer to a simulator function (prefered) or an external simulator executable must be defined." << std::endl;
        exit(100);
    }
};

// This defines a function type, for cleaner typing when using a function pointer as a simulator
// Using a function pointer is the typical approach for both compiling the abc library + simulator together AND
// using a dynamic simulator object.
typedef vector<float_type> AbcSimMPI(vector<float_type>, const unsigned long int, const unsigned long int, const ABC::MPI_par*);
typedef vector<float_type> AbcSimBase(vector<float_type>, const unsigned long int, const unsigned long int);


// This function handles loading a shared object file -> extracting the function pointer to an AbcSimF
template <typename AbcSimType>
inline AbcSimType * loadSO(const char * target) {
    void* handle = dlopen(target, RTLD_LAZY);
    if (!handle) {
        std::cerr << "Failed to open simulator object: " << target << " ; " << dlerror() << std::endl;
        exit(101);
    }
    auto simf = (AbcSimType*)dlsym(handle, "simulator");
    if(!simf) {
        std::cerr << "Failed to find 'simulator' function in " << target << " ; " << dlerror() << std::endl;
        dlclose(handle);
        exit(102);
    }
    return simf;
}

// an AbcSimFun built around an AbcSimF pointer. That pointer can come from code compiled along with this library,
// i.e. a executable that combines a simulator and the AbcSmc code, or be loaded from a shared object file.
struct AbcFPtrMPI : AbcSimFun {
    AbcSimMPI* fptr;
    AbcFPtrMPI(AbcSimMPI * _fptr) : fptr(_fptr) { } // constructor for a function pointer directly
    AbcFPtrMPI(const char * target) : AbcFPtrMPI(loadSO<AbcSimMPI>(target)) { } // construct from a char*-style string (the file name for shared object)
    AbcFPtrMPI(const string target) : AbcFPtrMPI(target.c_str()) { } // construct from a string (the file name for shared object)
    
    // for this version, we override the MPI version of the operator()
    // rather than just having it be an ignored parameter
    vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const override {
        return fptr(pars, seed, serial, _mp);
    }

    // and we also override the non-MPI version, to throw an error, since we explicitly asked
    // for MPI support
    vector<float_type> operator()(
        vector<float_type> /*pars*/, const unsigned long int /*seed*/, const unsigned long int /*serial*/
    ) const override {
        std::cerr << "ERROR: Explicitly constructed an MPI simulator, then called it without MPI arguments." << std::endl;
        exit(100);
    };
};

struct AbcFPtrBase : AbcSimFun {
    AbcSimBase* fptr;
    AbcFPtrBase(AbcSimBase * _fptr) : fptr(_fptr) { } // constructor for a function pointer directly
    AbcFPtrBase(const char * target) : AbcFPtrBase(loadSO<AbcSimBase>(target)) { } // construct from a char*-style string (the file name for shared object)
    AbcFPtrBase(const string target) : AbcFPtrBase(target.c_str()) { } // construct from a string (the file name for shared object)

    vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int seed, const unsigned long int serial
    ) const override {
        return fptr(pars, seed, serial);
    }
};

// an AbcSimFun built around an external executable. This is constructed with a command string to be executed in a shell,
// which should receive the parameters as a sequence of command line arguments and reply on standard out with the metrics
// as a series of numbers.
struct AbcExec : AbcSimFun {
    const string command;
    AbcExec(string _command) : command(_command) { }

    vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int /*seed*/, const unsigned long int /*serial*/
    ) const override {
        ostringstream execcom(command, std::ios_base::ate);
        vector<float_type> mets;
        for (const float_type par : pars) { execcom << " " << par; }

        FILE* pipe = popen(execcom.str().c_str(), "r");
        if (!pipe) {
            std::cerr << "ERROR: Unable to create pipe to " << execcom.str() << std::endl;
            exit(103);
        }

        char buffer[512];
        string retval = "";
        while(!feof(pipe)) {
            if(fgets(buffer, 512, pipe) != NULL) { retval += buffer; }
        }
        pclose(pipe);

        if (retval == "ERROR" or retval == "") {
            std::cerr << command << " does not exist or appears to be an invalid simulator." << std::endl;
            std::cerr << "Attempted: " << execcom.str().c_str() << std::endl;
        } else {
            istringstream ss(retval);
            // TODO deal with empty mets on !particle_success
            float_type met;
            while(ss >> met) mets.push_back(met);
        }

        return mets;
    }

};

#endif
