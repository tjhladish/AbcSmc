
#include <AbcSmc/AbcSim.h>
#include <iostream>
#include <dlfcn.h> // for dynamic version
#include <fstream> // for external executable version
#include <sstream> // for stringstream

using std::ostringstream;
using std::istringstream;

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


namespace ABC {

AbcSimFun* default_simulator() { return new AbcSimUnset(); }

AbcSimFun* SO_simulator(const char * target, const bool mpi = false) {
    if (mpi) {
        return new AbcFPtrMPI(target);
    } else {
        return new AbcFPtrBase(target);
    }
}

AbcSimFun* SO_simulator(const string target, const bool mpi = false) {
    if (mpi) {
        return new AbcFPtrMPI(target);
    } else {
        return new AbcFPtrBase(target);
    }
}

AbcSimFun* ptr_simulator(AbcSimBase* fptr) { return new AbcFPtrBase(fptr); };
AbcSimFun* ptr_simulator(AbcSimMPI* fptr) { return new AbcFPtrMPI(fptr); };
AbcSimFun* cmd_simulator(const string command) { return new AbcExec(command); };

}