#ifndef ABCSMC_CLI_H
#define ABCSMC_CLI_H

#include <iostream>
#include <unistd.h>
#include <vector>
#include <optional>
#include <string>
#include <gsl/gsl_rng.h>

using std::cerr;
using std::endl;
using std::string;

namespace ABC {

// A "usage" function for the built-in AbcSmc CLI
// @param cmd the name of the executable
// @param msg an optional message to print before the usage
// @param status if non-zero, will exit with this status
void usage(
    const string &cmd,
    const string &msg = "",
    const int status = 0
);

// Enum for the steps of the ABC-SMC process
// BUILD: setup a new database for this simulation configuration
// PROCESS: prepare the database for simulation; either inital sampling or performing PLS-ABC SMC step
// EVALUATE: run available simulations, potentially limited to a requested maximum
enum STEP { BUILD, PROCESS, EVALUATE };

// Stream insertion operator for ABC::STEP
std::ostream& operator<<(std::ostream &os, const STEP &step);

// Container for the parsed command line arguments
// @var config_file the path to the configuration file
// @var steps the `STEP`s to perform
// @var buffer_size the number of simulations to run in each batch (if not present, run all available)
// @var verbose the verbosity level (0 = quiet, 1 = normal, ... TBD)
struct CLIArgs {
    CLIArgs() = delete;
    CLIArgs(const string & cf) : config_file(cf) {
        // TODO assertions / error checking re config file?
    };

    string config_file;                  // based on config file ...
    std::vector<STEP> steps = {};        // ... do nothing by default
    std::optional<size_t> seed;          // ... with an unspecified seed
    std::optional<size_t> buffer_size;   // ... for all available runs
    size_t verbose = 0;                  // ... quietly
};

// parses the args passed to a typical main() function for an ABC-SMC program
// @param argc the number of arguments (per typical main() signature)
// @param argv the arguments (per typical main() signature)
CLIArgs parse_args(const size_t argc, const char * argv[]);

// Runs the ABC-SMC process, for some object that implements the ABCSMC *verbs*:
// - parse(a string [configuration file path], a size_t [verbosity level])
// - build(a verbosity level)
// - process(a random number generator, a verbosity level) n.b. RNG needed here, because we need to sample prior
// - evaluate(a buffer size, a verbosity level) n.b.: the RNG for *simulation* should be provided by the simulator
template<typename ABC>
inline void run(
    ABC* abc, const CLIArgs &args
) {

    // this should set storage, set sizes, etc 
    abc->parse(args.config_file, args.verbose);

    if (args.verbose > 0) {
        cerr << "Running ABC as: " << endl;
        for (auto it = args.steps.begin(); it != args.steps.end(); it++) {
            if (it != args.steps.begin()) { cerr << " => "; }
            cerr << *it;
        }
        cerr << endl;
    }

    for (auto step : args.steps) {
        switch(step) {
            case BUILD: abc->build(args.verbose); break;
            case PROCESS: abc->process(args.seed, args.verbose); break;
            case EVALUATE: abc->evaluate(
                args.buffer_size, // if buffer_size is an empty option, should internally do everything available
                args.verbose
            ); break;
            default:
                cerr << "Hit unimplemented STEP." << endl;
                exit(-1);
        }
    }

};


}

#endif // ABCSMC_CLI_H