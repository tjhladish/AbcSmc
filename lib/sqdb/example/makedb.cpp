#include "../sqdb.h"
int main() {
    sqdb::Db db("db.sdb");
    if ( !db.TableExists("test") ) {
        db.Query("create table test(i int, d double);").Next();
        db.Query("insert into test  values(1, 2.353);").Next();
        db.Query("insert into test  values(2, 3.353);").Next();
        db.Query("insert into test values(4, 5.353);").Next();
    }
}
