sqlite3.o: sqlite3.c sqlite3.h
	gcc -O2 -c sqlite3.c -DSQLITE_THREADSAFE=2

sqdb.o: sqdb.cpp sqdb.h sqlite3.o
	g++ -O2 -c sqdb.cpp

clean:
	rm -f *.o
