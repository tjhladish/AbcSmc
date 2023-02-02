
#include "AbcSmc.h"
#include "examples.h"
#include "dice.h"

// demonstrates how to use the abc library compiled with a static simulator
int main(int argc, char* argv[]) {

    // convenience method for checking arguments / alerting usage; from examples.h
    check_args("abc_sql", argc);

    // convenience method for parsing arguments; from examples.h
    CLIArgs args = parse_args("abc_sql", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(args.config_file);
    abc->set_simulator(simulator);
    // simulator defined in dice.h and compiled into this executable
    // when using this approach, the simulator must be manually set

    // convenience method for the core abc application loop; from examples.h
    abc_loop(abc, args, RNG);

    return 0;
}
