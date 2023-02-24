#ifndef EXAMPLES_H
#define EXAMPLES_H

#include <iostream>
#include <gsl/gsl_rng.h>
#include <cstring>
#include <unistd.h>
#include <vector>
#include <optional>
#include <algorithm>

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

// TODO: desire to allow e.g. -bpe; also -bpepen 10 etc
// TODO: desire verbosity levels?

void usage(
    const std::string & cmd,
    const std::string & msg = "",
    const int status = 0 // if non-zero, will exit
) {
    std::cerr << msg << std::endl;
    std::cerr << "Usage: ./" << cmd << " config.json [-option|--option (each space separated)]" << std::endl << std::endl;
    std::cerr << "Core options:" << std::endl;
    std::cerr << "\t\t -(-b)uild    : setup a new database for this simulation configuration." << std::endl;
    std::cerr << "\t\t -(-p)rocess  : prepare the database for simulation; either inital sampling or performing PLS-ABC SMC step." << std::endl;
    std::cerr << "\t\t -(-e)valuate : run available simulations, potentially limited by -n flag." << std::endl;
    std::cerr << "\t\t -n 1234      : number of simulations to attempt; implies -e; when -e specified and -n unspecified, assumes attempt all available runs." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Auxilary options:" << std::endl;
    std::cerr << "\t\t -(-h)elp     : print this message; ignore all other options." << std::endl;
    std::cerr << "\t\t -(-v)erbose  : when working, be effusive." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Streamlined combination options:" << std::endl;
    std::cerr << "\t\t -(-s)etup    : implies -b -p, i.e. build the database and do initial population." << std::endl;
    std::cerr << "\t\t -(-c)ycle    : implies -e -p, i.e. evaluate (all if -n unspecified) => prepare for next round." << std::endl;
    std::cerr << "\t\t -(-a)ll      : implies build, process, evaluate; ignores -n." << std::endl;
    std::cerr << std::endl;
    std::cerr << "Example uses:" << std::endl;
    std::cerr << std::endl;
    std::cerr << "\t\t $ # For HPC based fitting:" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -s # prior to job submission" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -e -n 100 # e.g. in torque array job script, $ARRAY_SIZE * 100 ~ 1.1*total sims to do" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -p # after job completion, then resubmit torque script as desired" << std::endl;
    std::cerr << "\t\t $ # or for HPC based scenario analysis, just the first two steps." << std::endl;
    std::cerr << std::endl;
    std::cerr << "\t\t $ # Local usage:" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -a # first invocation; if just doing scenario analysis, done!" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -p # if fitting, prepare for future cycles" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -c # if fitting, do an evaluate => process cycle" << std::endl;
    std::cerr << "\t\t $ ./" << cmd << " config.json -c # if fitting, repeat as necessary" << std::endl;
    if (status != 0) { exit(status); }
}

enum STEP { BUILD, PROCESS, EVALUATE };

std::ostream& operator<<(std::ostream& os, const STEP& step) {
    switch (step) {
        case BUILD: os << "BUILD"; break;
        case PROCESS: os << "PROCESS"; break;
        case EVALUATE: os << "EVALUATE"; break;
    }
    return os;
}

struct CLIArgs {
    CLIArgs() = delete;
    CLIArgs(const std::string & cf) : config_file(cf) {
        // TODO assertions / error checking re config file?
    };

    std::string config_file;             // based on config file ...
    std::vector<STEP> steps = {};        // ... do nothing by default
    std::optional<size_t> buffer_size;   // ... for all available runs
    size_t verbose = 0;                  // ... quietly
};

bool argcheck(const char * arg, const char * short_arg, const char * long_arg) {
    return (strcmp(arg, short_arg) == 0) or (strcmp(arg, long_arg) == 0);
}

// TODO: should this be a class method?
CLIArgs parse_args(const std::string & cmd, const int argc, const char * argv[]) {

    // check for help requested
    if ((std::find(argv, argv + argc, "-h") != argv + argc) or (std::find(argv, argv + argc, "--help") != argv + argc)) {
        usage(cmd);
        exit(0);
    }

// assert argv[0] = this program
// assert argv[1] = config file

    auto args = CLIArgs(std::string(argv[1]));

    for (int i = 2; i < argc;  i++ ) {

        if ( argcheck(argv[i], "-b", "--build") ) {
            if (args.steps.size() > 0) { usage(cmd, "Error: -(-b)uild must be the first step.", 102); }
            args.steps.push_back(BUILD);
        } else if ( argcheck(argv[i], "-p", "--process") ) {
            args.steps.push_back(PROCESS);
        } else if ( argcheck(argv[i], "-e", "--evaluate") ) {
            args.steps.push_back(EVALUATE);
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            // this will occur if -n is the last argument, i.e. no number provided after
            if (i == (argc - 1)) { usage(cmd, "Error: -n must be followed by a positive integer.", 103); }
            args.buffer_size.emplace(atoi(argv[i++]));
            // this will occur if provided a negative number or a non-integer (results in 0)
            if (args.buffer_size.value() < 1) { usage(cmd, "Error: -n must be followed by a positive integer.", 103); }
        } else if ( argcheck(argv[i], "-s", "--setup") ) {
            if (args.steps.size() > 0) { usage(cmd, "Error: -(-s)etup must be the first step.", 102); }
            args.steps.push_back(BUILD);
            args.steps.push_back(PROCESS);
        } else if ( argcheck(argv[i], "-c", "--cycle") ) {
            args.steps.push_back(EVALUATE);
            args.steps.push_back(PROCESS);
        } else if ( argcheck(argv[i], "-a", "--all") ) {
            if (args.steps.size() > 0) { usage(cmd, "Error: -(-a)ll must be the first step.", 102); }
            args.steps.push_back(BUILD);
            args.steps.push_back(PROCESS);
            args.steps.push_back(EVALUATE);
        } else if ( argcheck(argv[i], "-v", "--verbose") ) {
            args.verbose += 1;
        } else {
            usage(cmd, "Error: unrecognized argument: " + string(argv[i]), 104);
        }
    }

    return args;
};

template<typename ABC>
inline void do_abc(
    ABC * abc, const CLIArgs & args, const gsl_rng * RNG
) {

    // this should set storage, set sizes, etc 
    abc->parse(args.config_file);

    if (args.verbose > 0) {
        std::cerr << "Running ABC as: " << std::endl;
        for (auto it = args.steps.begin(); it != args.steps.end(); it++) {
            if (it != args.steps.begin()) { std::cerr << " => "; }
            std::cerr << *it;
        }
        std::cerr << std::endl;
    }

    for (auto step : args.steps) {
        switch(step) {
            case BUILD: abc->build(args.verbose); break;
            case PROCESS: abc->process(RNG, args.verbose); break;
            case EVALUATE: abc->evaluate(
                RNG, args.buffer_size, // if buffer_size is an empty option, should internally do everything available
                args.verbose
            ); break;
            default:
                std::cerr << "Hit unimplemented STEP." << std::endl;
                exit(-1);
        }
    }

};

#endif // EXAMPLES_H