cmake_minimum_required( VERSION ${CMAKE_VERSION} )

project( SQDB )

include_directories( include )

add_library( sqlite3 OBJECT src/sqlite3.c )
add_library( sqdb OBJECT src/sqdb.cpp )