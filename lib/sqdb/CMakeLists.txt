cmake_minimum_required( VERSION ${CMAKE_VERSION} )

project( SQDB )

add_library( sqlite3 OBJECT src/sqlite3.c )
target_include_directories( sqlite3 PUBLIC include )

add_library( sqdb OBJECT src/sqdb.cpp )
target_link_libraries(sqdb sqlite3) # sqdb depends on sqlite3
target_include_directories( sqdb PUBLIC include )