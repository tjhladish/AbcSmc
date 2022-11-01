
sqdb.o: sqdb.cpp sqdb.h sqlite3.h
	g++ -O2 -Wall -fPIC -I. -c sqdb.cpp -o $@

sqlite3.o: sqlite3.c sqlite3.h
	gcc -O2 -Wall -fPIC -c sqlite3.c -DSQLITE_THREADSAFE=2

clean:
	rm -f *.o
	rm -f *.so
