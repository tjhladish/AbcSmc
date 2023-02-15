
default: sqlite3.o sqdb.o pivot_vtab.so

sqlite3.o: sqlite3.c sqlite3.h
	gcc -O2 -c sqlite3.c -DSQLITE_THREADSAFE=2

REQEXTS := eval.c

# TODO: static incorporation of sqlite extensions: https://www.sqlite.org/loadext.html topic 6
# compiles sqdb library along with any required extensions into a single object
sqdb.o: sqdb.cpp sqdb.h sqlite3.o $(REQEXTS)
	g++ -O2 -c sqdb.cpp -DSQLITE_CORE


# if this worked, might prefer it over eval approach - creates virtual pivot tables, which could be inserted to perhaps more conveniently
# would also eliminate the need to drop/create met/par views when updating the definitions
# https://github.com/jakethaw/pivot_vtab
# https://stackoverflow.com/questions/65433578/assitance-in-statically-linking-a-sqlite-loadable-extension-in-a-c-program-compi
pivot_vtab.so: pivot_vtab.c
	gcc -g -fPIC -shared $^ -o $@

# from: https://www.sqlite.org/src/dir?ci=trunk&name=ext/misc
# use with https://dba.stackexchange.com/questions/266614/how-to-pivot-data-in-sqlite
# after loading parameter / metric definitions, can use this to define the views conveniently
# or update the views to use the new definitions (if pars / mets get renamed / extended / etc)
eval.so: eval.c
	gcc -g -fPIC -shared $^ -o $@

clean:
	rm -f *.o *.so
