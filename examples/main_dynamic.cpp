#include <dlfcn.h>
#include "AbcSmc.h"
#include "examples.h"

int main(int argc, char* argv[]) {

    check_args("abc_dynamic", argc);

    CLIArgs args = parse_args("abc_dynamic", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(args.config_file);
    // simulator should be set from config file via shared object

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

    return 0;
}
