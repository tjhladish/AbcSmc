#include <AbcSmc/CLI.h>

#include <cstring>
#include <unistd.h>
#include <algorithm>

using std::cerr;
using std::endl;
using std::string;

namespace ABC {

void usage(
    const string &cmd,
    const string &msg,
    const int status
) {
    const string ident = "\t";
    cerr << msg << endl;
    cerr << "Usage: " << cmd << " config.json [-option|--option (each space separated)]" << endl << endl;
    cerr << "Core options:" << endl;
    cerr << ident << "-(-b)uild    : setup a new database for this simulation configuration." << endl;
    cerr << ident << "-(-p)rocess  : prepare the database for simulation; either inital sampling or performing PLS-ABC SMC step." << endl;
    cerr << ident << "-(-e)valuate : run available simulations, potentially limited by -n flag." << endl;
    cerr << ident << "-n 1234      : number of simulations to attempt; implies -e; when -e specified and -n unspecified, assumes attempt all available runs." << endl;
    cerr << endl;
    cerr << "Auxilary options:" << endl;
    cerr << ident << "-(-h)elp     : print this message; ignore all other options." << endl;
    cerr << ident << "-(-v)erbose  : when working, be effusive." << endl;
    cerr << endl;
    cerr << "Streamlined combination options:" << endl;
    cerr << ident << "-(-s)etup    : implies -b -p, i.e. build the database and do initial population." << endl;
    cerr << ident << "-(-c)ycle    : implies -e -p, i.e. evaluate (all if -n unspecified) => prepare for next round." << endl;
    cerr << ident << "-(-a)ll      : implies build, process, evaluate; ignores -n." << endl;
    cerr << endl;
    cerr << "Example uses:" << endl;
    cerr << endl;
    cerr << "$ # For HPC based fitting:" << endl;
    cerr << "$ " << cmd << " config.json -s # prior to job submission" << endl;
    cerr << "$ " << cmd << " config.json -e -n 100 # e.g. in torque array job script, $ARRAY_SIZE * 100 ~ 1.1*total sims to do" << endl;
    cerr << "$ " << cmd << " config.json -p # after job completion, then resubmit torque script as desired" << endl;
    cerr << "$ # or for HPC based scenario analysis, just the first two steps." << endl;
    cerr << endl;
    cerr << "$ # Local usage:" << endl;
    cerr << "$ " << cmd << " config.json -a # first invocation; if just doing scenario analysis, done!" << endl;
    cerr << "$ " << cmd << " config.json -p # if fitting, prepare for future cycles" << endl;
    cerr << "$ " << cmd << " config.json -c # if fitting, do an evaluate => process cycle" << endl;
    cerr << "$ " << cmd << " config.json -c # if fitting, repeat as necessary" << endl;
    if (status != 0) { exit(status); }
}

std::ostream& operator<<(std::ostream &os, const STEP &step) {
    switch (step) {
        case BUILD: os << "BUILD"; break;
        case PROCESS: os << "PROCESS"; break;
        case EVALUATE: os << "EVALUATE"; break;
        default: os << "UNDEFINED ABC::STEP"; exit(-1);
    }
    return os;
}

// non-exported helper function for finding argument flags
bool argcheck(const char * arg, const char * short_arg, const char * long_arg) {
    return (strcmp(arg, short_arg) == 0) or (strcmp(arg, long_arg) == 0);
}

CLIArgs parse_args(const size_t argc, const char * argv[]) {

    // assert argv[0] = this program
    const string cmd = string(argv[0]);

    // check for help requested
    if ((std::find(argv, argv + argc, "-h") != argv + argc) or (std::find(argv, argv + argc, "--help") != argv + argc)) {
        usage(cmd);
        return CLIArgs("");
    }

    // assert argv[1] = config file, if not in "help" mode

    auto args = CLIArgs(string(argv[1]));

    for (size_t i = 2; i < argc; i++) {

        if (argcheck(argv[i], "-b", "--build")) {
            if (not args.steps.empty()) { 
                usage(cmd, "Error: -(-b)uild must be the first step.", 102); 
            } else {
                args.steps.push_back(BUILD);
            }
        } else if (argcheck(argv[i], "-p", "--process")) {
            if ((not args.steps.empty()) and (args.steps.back() == PROCESS)) {
                cerr << "WARNING: -(-p)rocess specified multiple times in a row; ignoring redundant invocation." << endl;
            } else {
                args.steps.push_back(PROCESS);
            }
        } else if (argcheck(argv[i], "-e", "--evaluate")) {
            if ((not args.steps.empty()) and (args.steps.back() == EVALUATE)) {
                cerr << "WARNING: -(-e)valuate specified multiple times in a row; ignoring redundant invocation." << endl;
            } else {
                args.steps.push_back(EVALUATE);
            }
        } else if (strcmp(argv[i], "-n" ) == 0) {
            // this will occur if -n is the last argument, i.e. no number provided after
            if (i == (argc - 1)) { usage(cmd, "Error: -n must be followed by a positive integer.", 103); }
            args.buffer_size.emplace(atoi(argv[i++]));
            // this will occur if provided a negative number or a non-integer (results in 0)
            if (args.buffer_size.value() < 1) { usage(cmd, "Error: -n must be followed by a positive integer.", 103); }
        } else if (argcheck(argv[i], "-s", "--setup")) {
            if (args.steps.size() > 0) { usage(cmd, "Error: -(-s)etup must be the first step.", 102); }
            args.steps.push_back(BUILD);
            args.steps.push_back(PROCESS);
        } else if (argcheck(argv[i], "-c", "--cycle")) {
            args.steps.push_back(EVALUATE);
            args.steps.push_back(PROCESS);
        } else if (argcheck(argv[i], "-a", "--all")) {
            if (args.steps.size() > 0) { usage(cmd, "Error: -(-a)ll must be the first step.", 102); }
            args.steps.push_back(BUILD);
            args.steps.push_back(PROCESS);
            args.steps.push_back(EVALUATE);
        } else if (argcheck(argv[i], "-v", "--verbose")) {
            args.verbose += 1;
        } else {
            usage(cmd, "Error: unrecognized argument: " + string(argv[i]), 104);
        }
    }

    return args;
};

} // namespace ABC