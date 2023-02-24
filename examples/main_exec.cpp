#include "AbcSmc.h"
#include "examples.h"

int main(int argc, char* argv[]) {

    // convenience method for parsing arguments; from examples.h
    CLIArgs args = parse_args("abc_exec", argc, argv);

    AbcSmc* abc = new AbcSmc();
    do_abc(abc, args, RNG);

    return 0;
}
