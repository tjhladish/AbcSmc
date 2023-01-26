#include <dlfcn.h>
#include "AbcSmc.h"
#include "examples.h"

int main(int argc, char* argv[]) {

    check_args("abc_dynamic", argc);

    CLIArgs args = parse_args("abc_dynamic", argc, argv);

    AbcSmc* abc = new AbcSmc();
    abc->parse_config(args.config_file);
    // simulator should be set from config file - no need to invoke abc->set_simulator()
    // the simulator file is specified with the "shared" key in the dynamic.json file.
    // note that this version doesn't need to know anything about dice.h/cpp or dice_game.cpp

    // see examples.h for the core abc loop
    abc_loop(abc, args, RNG);

    return 0;
}
