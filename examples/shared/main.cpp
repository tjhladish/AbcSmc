#include <AbcSmc/AbcSmc.h>
#include <AbcSmc/CLI.h>

using namespace ABC;

int main(int argc, char* argv[]) {

    auto args = ABC::parse_args(argc, argv);

    AbcSmc* abc = new AbcSmc();
    run(abc, args);

    return 0;
}
