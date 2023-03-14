cmake_minimum_required( VERSION ${CMAKE_VERSION} )

project( CCRC32 )

include_directories( include )

add_library( CCRC32 OBJECT src/CCRC32.cpp )