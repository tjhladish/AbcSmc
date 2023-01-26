#include "AbcSmc.h"

#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics_double.h>
#include <math.h>
#include <unistd.h>
#include "examples.h"
#include "dice.h"

using namespace std;

int main(int argc, char* argv[]) {

    check_args("abc_sql", argc);

    CLIArgs args = parse_args("abc_sql", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(args.config_file);
    // simulator defined in dice.h, which is compiled into this executable
    abc->set_simulator(simulator);

    // see examples.h for the core abc loop
    abc_loop(abc, args, RNG);

    return 0;
}
