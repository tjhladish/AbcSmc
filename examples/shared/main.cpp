#include <AbcSmc/AbcSmc.h>
#include "examples.h"

using namespace ABC;

int main(int argc, char* argv[]) {

    // convenience method for checking arguments / alerting usage; from examples.h
    check_args("shared", argc);

    // convenience method for parsing arguments; from examples.h
    CLIArgs args = parse_args("shared", argc, argv);

    auto abc = std::make_shared<AbcSmc>();
    abc->parse_config(args.config_file);
    // simulator should be set from config file - no need to invoke abc->set_simulator()
    // the simulator file is specified with the "shared" key in the dynamic.json file.
    // note that for this approach, this executable doesn't need to know anything
    // about dice.h/cpp or dice_game.cpp

    // convenience method for the core abc application loop; from examples.h
    abc_loop(abc, args, RNG);

    return 0;
}
