

#include <iostream>
#include <fstream>
#include "AbcStorage.h"

using namespace std;

void usage(const char * progname) {
    cerr << "\tUsage: " << progname << " $DB.sqlite" << endl;
    cerr << "\t\tn.b., $DB.sqlite should not exist." << endl;
    exit(-1);
}

int main(int argc, char* argv[]) {

    if (argc != 2) usage(argv[0]);

    AbcSQLite db(argv[1]);

    // should be able to create the storage object without requiring the file to exist
    cout << "Exists? " << db.exists() << endl;

    bool didSetup = db.setup({"ndice", "sides"}, {"sum", "sd"}, false, false, true);

    // after setup, the file should always exist
    cout << "Did setup? " << didSetup << endl << "Exists after setup? " << db.exists() << endl;

    didSetup = db.setup({}, {}, false, false, true);

    cout << "Second setup attempt should fail: " << didSetup << endl;   

    return 0;
}
