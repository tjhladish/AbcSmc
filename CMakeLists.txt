cmake_minimum_required(VERSION 3.16)

project(AbcSmc
    VERSION 1.5.0
    HOMEPAGE_URL "https://github.com/tjhladish/AbcSmc"
    LANGUAGES CXX
)

set( CMAKE_CXX_STANDARD_REQUIRED 17 )
set( CMAKE_CXX_FLAGS "-Wall -Wpedantic -O2 -fPIC" )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_MESSAGE LAZY)

add_subdirectory(lib/PLS ${CMAKE_BINARY_DIR}/PLS)
add_subdirectory(lib/CCRC32 ${CMAKE_BINARY_DIR}/CCRC32)
add_subdirectory(lib/sqdb ${CMAKE_BINARY_DIR}/sqdb)

# add_subdirectory(lib/jsoncpp ${CMAKE_BINARY_DIR}/jsoncpp)
# this is a workaround for the apparent lack of a jsoncpp library target
add_library(jsoncpp STATIC lib/jsoncpp/src/lib_json/json_reader.cpp lib/jsoncpp/src/lib_json/json_value.cpp lib/jsoncpp/src/lib_json/json_writer.cpp)
target_include_directories(jsoncpp PUBLIC lib/jsoncpp/include lib/jsoncpp/src/lib_json)

add_library(abc STATIC
    $<TARGET_OBJECTS:CCRC32> $<TARGET_OBJECTS:sqdb> $<TARGET_OBJECTS:sqlite3>
    src/AbcLog.cpp src/AbcMPI.cpp src/AbcSmc.cpp src/AbcUtil.cpp
    src/CLI.cpp
)

find_package(GSL REQUIRED) # the abc library objects depend on the GSL library
target_link_libraries(abc pls sqdb jsoncpp m GSL::gsl GSL::gslcblas ${CMAKE_DL_LIBS})
# the abc library objects depend on pls, sqdb, jsoncpp, m (aka math.h), GSL, and dl (aka dl.h) libraries

# look in these directories for header files, also make them available to other projects that link to abc
target_include_directories(abc PUBLIC include src lib)

add_subdirectory(examples)

install(TARGETS abc
  RUNTIME_DEPENDENCIES
  RUNTIME DESTINATION bin
  LIBRARY DESTINATION lib
  ARCHIVE DESTINATION lib
)