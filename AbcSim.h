#ifndef ABCSIM_H
#define ABCSIM_H

#include <vector> // container for receiving parameters / returning metrics
#include <dlfcn.h> // for dynamic version
#include <fstream> // for external executable version
#include <sstream> // for stringstream
#include <string> // for string

#include "pls/pls.h" // for float_type
#include "AbcMPIPar.h"

using std::vector;

// if compiling with MPI support, have to handle slightly complex MPI objects
// this must be defined in a matching way between and AbcSmc executables *and*
// other code using this header.
namespace ABC {

  template <typename T>
  inline std::string toString (const T& t) {
      std::stringstream ss;
      ss << t;
      return ss.str();
  }

}

// defines the core abstraction for a simulator: (1) a functor, (2) with an operator(), (3) return type vector<float> (the metrics),
// (4) arguments vector<float> (the parameters), and unsigned long int, unsigned long int, MPI_par* (the abc seed/serial/mpi information)
// implemented as a pure abstract class - i.e. must be extended by concrete implementations
// In general, the goal is that those implementations should *not* be "stateful". That is, they should not have any internal state which
// is changed when they are used.
struct AbcSimFun {
    virtual vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const = 0;
};

// an AbcSimFun which throws an error if used. This is intended to be used as a default, so that if *not* replaced, the error
// will be thrown when the simulator is used.
struct AbcSimUnset : AbcSimFun {
    vector<float_type> operator()(
      vector<float_type> /*pars*/, const unsigned long int /*seed*/, const unsigned long int /*serial*/, const ABC::MPI_par* /*_mp*/
    ) const {
        std::cerr << "ERROR: A pointer to a simulator function (prefered) or an external simulator executable must be defined." << std::endl;
        exit(100);
    }
};

// This defines a function type, for cleaner typing when using a function pointer as a simulator
// Using a function pointer is the typical approach for both compiling the abc library + simulator together AND
// using a dynamic simulator object.
typedef vector<float_type> AbcSimF(vector<float_type>, const unsigned long int, const unsigned long int, const ABC::MPI_par*);

// This function handles loading a shared object file -> extracting the function pointer to an AbcSimF
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

// an AbcSimFun built around an AbcSimF pointer. That pointer can come from code compiled along with this library,
// i.e. a executable that combines a simulator and the AbcSmc code, or be loaded from a shared object file.
struct AbcFPtr : AbcSimFun {
    AbcSimF* fptr;
    AbcFPtr(AbcSimF * _fptr) : fptr(_fptr) { } // constructor for a function pointer directly
    AbcFPtr(const char * target) : AbcFPtr(loadSO(target)) { } // construct from a char*-style string (the file name for shared object)
    AbcFPtr(const std::string target) : AbcFPtr(target.c_str()) { } // construct from a std::string (the file name for shared object)
    vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int seed, const unsigned long int serial, const ABC::MPI_par* _mp
    ) const {
        return fptr(pars, seed, serial, _mp);
    }
};

// an AbcSimFun built around an external executable. This is constructed with a command string to be executed in a shell,
// which should receive the parameters as a sequence of command line arguments and reply on standard out with the metrics
// as a series of numbers.
struct AbcExec : AbcSimFun {
    const std::string command;
    AbcExec(std::string _command) : command(_command) { }

    vector<float_type> operator()(
      vector<float_type> pars, const unsigned long int /*seed*/, const unsigned long int /*serial*/, const ABC::MPI_par* /*_mp*/
    ) const {
        auto execcom = command;
        vector<float_type> mets;
        for (auto par : pars) { execcom += " " + ABC::toString(par); }

        FILE* pipe = popen(execcom.c_str(), "r");
        if (!pipe) {
            std::cerr << "ERROR: Unable to create pipe to " << execcom << std::endl;
            exit(103);
        }

        char buffer[512];
        std::string retval = "";
        while(!feof(pipe)) {
            if(fgets(buffer, 512, pipe) != NULL) { retval += buffer; }
        }
        pclose(pipe);

        if (retval == "ERROR" or retval == "") {
            std::cerr << command << " does not exist or appears to be an invalid simulator" << std::endl;
        } else {
            std::stringstream ss;
            ss.str(retval);
            // TODO deal with empty mets on !particle_success
            float_type met;
            while(ss >> met) mets.push_back(met);
        }

        return mets;
    }

};

#endif
