#ifndef ABCSIM_H
#define ABCSIM_H

#include <vector> // container for receiving parameters / returning metrics
#include <string> // for string

#include <PLS/types.h> // for float_type
#include <AbcSmc/AbcMPIPar.h>

// TODO approach this way?
// https://stackoverflow.com/questions/58394556/c-concepts-can-i-have-a-constraint-requiring-a-function-be-present-in-a-clas

using std::vector;
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

// This defines a function type, for cleaner typing when using a function pointer as a simulator
// Using a function pointer is the typical approach for both compiling the abc library + simulator together AND
// using a dynamic simulator object.
typedef vector<float_type> AbcSimBase(vector<float_type>, const unsigned long int, const unsigned long int);
typedef vector<float_type> AbcSimMPI(vector<float_type>, const unsigned long int, const unsigned long int, const ABC::MPI_par*);

// Lastly, here's what we expect people to actually invoke
namespace ABC {

    AbcSimFun* default_simulator();
    AbcSimFun* SO_simulator(const char * target, const bool mpi = false);
    AbcSimFun* SO_simulator(const string target, const bool mpi = false);
    AbcSimFun* ptr_simulator(AbcSimBase* fptr);
    AbcSimFun* ptr_simulator(AbcSimMPI* fptr);
    AbcSimFun* cmd_simulator(const string command);

}

#endif
