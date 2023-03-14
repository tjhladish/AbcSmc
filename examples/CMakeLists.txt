cmake_minimum_required( VERSION 3.24.2 )

project( AbcSmcExamples )

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_FLAGS "-Wall -Wpedantic -O2 -fPIC" )

# look in the parent directory for assorted dependencies
add_subdirectory(.. ${CMAKE_SOURCE_DIR}/../build)
# the example specific includes
include_directories(include)

add_executable(direct direct/main.cpp)
target_link_libraries(direct abc) # the demos depends on abc

add_executable(shared shared/main.cpp)
target_link_libraries(shared abc) # the demos depends on abc
