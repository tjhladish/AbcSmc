cmake_minimum_required( VERSION 3.24.2 )

project( AbcSmc )

include_directories( include src lib lib/CCRC32 lib/jsoncpp/include lib/pls lib/pls/eigen lib/sqdb)

add_subdirectory(lib/jsoncpp)

add_library(abc STATIC src/AbcLog.cpp src/AbcMPI.cpp src/AbcSmc.cpp src/AbcUtil.cpp)