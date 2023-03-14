cmake_minimum_required( VERSION 3.24.2 )

project( AbcSmc )

include_directories( include src lib lib/CCRC32 lib/jsoncpp/include lib/pls lib/pls/eigen lib/sqdb)

add_subdirectory(lib/jsoncpp)
add_subdirectory(lib/pls)
add_subdirectory(lib/CCRC32)
add_subdirectory(lib/sqdb)

add_library(abc STATIC $<TARGET_OBJECTS:CCRC32> $<TARGET_OBJECTS:sqdb> $<TARGET_OBJECTS:sqlite3> src/AbcLog.cpp src/AbcMPI.cpp src/AbcSmc.cpp src/AbcUtil.cpp)