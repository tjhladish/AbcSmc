cmake_minimum_required( VERSION 3.24.2 )

project( AbcSmc )

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_FLAGS "-Wall -Wpedantic -O2 -fPIC" )

# look in these directories for header files
include_directories(include src lib lib/jsoncpp/include)

# look in these libraries for dependencies, build their pertinent targets in the project build directory
add_subdirectory(lib/jsoncpp ${CMAKE_BINARY_DIR}/jsoncpp)
add_subdirectory(lib/PLS ${CMAKE_BINARY_DIR}/PLS)
add_subdirectory(lib/CCRC32 ${CMAKE_BINARY_DIR}/CCRC32)
add_subdirectory(lib/sqdb ${CMAKE_BINARY_DIR}/sqdb)

add_library(abc STATIC $<TARGET_OBJECTS:CCRC32> $<TARGET_OBJECTS:sqdb> $<TARGET_OBJECTS:sqlite3> src/AbcLog.cpp src/AbcMPI.cpp src/AbcSmc.cpp src/AbcUtil.cpp)
target_link_libraries(abc pls sqdb json) # the abc library objects depend on the pls, sqdb libraries