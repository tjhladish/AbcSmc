#include "AbcSmc.h"
#include <gsl/gsl_rng.h>
#include <unistd.h>
#include "examples.h"

using namespace std;

int main(int argc, char* argv[]) {

    check_args("abc_exec", argc);

    CLIArgs args = parse_args("abc_exec", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(args.config_file);
    // simulator should be set from config file - no need to invoke abc->set_simulator()
    // the simulator file is specified with the "executable" key in the exec.json file.
    // note that this version doesn't need to know anything about dice.h/cpp or dice_game.cpp

    // see examples.h for the core abc loop
    abc_loop(abc, args, RNG);

    return 0;
}
