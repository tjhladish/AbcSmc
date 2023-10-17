
#include <AbcSmc/AbcSmc.h>
#include <AbcSmc/CLI.h>
#include "dice.h" // simulator

using namespace ABC;

int main(const size_t argc, const char* argv[]) {

    // convenience method for parsing arguments; from examples.h
    CLIArgs args = ABC::parse_args(argc, argv);

    // simulator defined in dice.h and compiled into this executable
    // when using this approach, the simulator must be manually set
    AbcSmc* abc = new AbcSmc();
    abc->set_simulator(simulator);

    // convenience method for the core abc application loop; from examples.h
    ABC::run(abc, args);

    return 0;
}
