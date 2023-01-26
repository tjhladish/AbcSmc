
#ifndef EXAMPLES_H
#define EXAMPLES_H

#include <iostream>
#include <gsl/gsl_rng.h>
#include <cstring>
#include <unistd.h>

const gsl_rng* RNG = gsl_rng_alloc(gsl_rng_taus2);

void usage(const std::string cmd) {
    std::cerr << std::endl;
    std::cerr << "\tUsage: ./" << cmd << " config.json --process" << std::endl << std::endl;
    std::cerr << "\t       ./" << cmd << " config.json --simulate" << std::endl << std::endl;
    std::cerr << "\t       ./" << cmd << " config.json --simulate -n <number of simulations per database write>" << std::endl << std::endl;
    std::cerr << "\t       ./" << cmd << " config.json --process --simulate -n <number of simulations per database write>" << std::endl << std::endl;
}

void check_args(const std::string cmd, int argc) {
    if (not (argc == 3 or argc == 5 or argc == 6) ) {
        usage(cmd);
        exit(100);
    }
}

struct CLIArgs {
    std::string config_file;
    bool process_db = false;
    bool simulate_db = false;
    bool do_all = false;
    size_t buffer_size = 1;
};

CLIArgs parse_args(std::string cmd, int argc, char* argv[]) {
    CLIArgs args;
    args.config_file = std::string(argv[1]);

    for (int i = 2; i < argc;  i++ ) {
        if ( strcmp(argv[i], "--process") == 0  ) {
            args.process_db = true;
        } else if ( strcmp(argv[i], "--simulate") == 0  ) {
            args.simulate_db = true;
        } else if ( strcmp(argv[i], "-n" ) == 0 ) {
            args.buffer_size = atoi(argv[++i]);
        } else if ( strcmp(argv[i], "--all" ) == 0 ) {
            args.do_all = true;
        } else {
            usage(cmd);
            exit(101);
        }
    }

    return args;
};

void abc_loop(AbcSmc * abc, CLIArgs& args, const gsl_rng * RNG) {
    size_t set_count = abc->get_smc_iterations();

    for (size_t i = 0; i < set_count; ++i) {
        auto buffer_size = args.buffer_size;

        if (args.do_all) {
            buffer_size = abc->get_num_particles(i, QUIET);
        }

        if (args.process_db) {
            gsl_rng_set(RNG, time(NULL) * getpid()); // seed the rng using sys time and the process id
            abc->process_database(RNG);
        } 

        if (args.simulate_db) {
            abc->simulate_next_particles(buffer_size);
        }
    }

    if (args.do_all) {
        abc->process_database(RNG); // one last time, to get the posterior
    }
};

#endif // EXAMPLES_H