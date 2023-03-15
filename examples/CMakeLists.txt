cmake_minimum_required( VERSION 3.24.2 )

project( AbcSmcExamples )

set( CMAKE_CXX_STANDARD 20 )
set( CMAKE_CXX_FLAGS "-Wall -Wpedantic -O2 -fPIC" )

# look in the parent directory for assorted dependencies
add_subdirectory(.. ${CMAKE_SOURCE_DIR}/../build)
# the example specific includes
include_directories(include)

# the demos all depend on abc

# TODO: not yet working
# add_executable(direct direct/main.cpp)
# target_link_libraries(direct abc)

add_executable(shared shared/main.cpp)
target_link_libraries(shared abc)

add_executable(integral integral/main.cpp)
target_link_libraries(integral abc)

add_executable(executable executable/main.cpp)
target_link_libraries(executable abc)

# TODO this must be foreach-able?