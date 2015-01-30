#include "../sqdb.h"
#include <iostream>

using namespace std;
using namespace sqdb;

int main(int argc, char** argv) {
    sqdb::Db db("db.sdb");
    // Select all tuples from table t with two columns integer and float.
    //Statement s = db.Query("select * from test where i = ?;");
    QueryStr str;
    Statement s = db.Query(str.Format(SQDB_MAKE_TEXT("select * from test where i = %s;"), argv[1])); 
    while ( s.Next() ) {
        int i = s.GetField(0);    // Get the first column.
        double d = s.GetField(1); // Get second column.
        cout << i << " " << d << endl;
    }
}
