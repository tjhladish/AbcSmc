cmake_minimum_required( VERSION 3.16 )

project( AbcSmc VERSION 2.0.1 LANGUAGES CXX)

set( CMAKE_CXX_STANDARD 17 )
set( CMAKE_CXX_FLAGS "-Wall -Wpedantic -O2 -fPIC" )

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

set(CMAKE_INSTALL_MESSAGE LAZY)

find_package(PLS)

if (NOT PLS_FOUND)
  message(WARNING "PLS 2.0.0+ not found; installing from source...")
  include(FetchContent)
  FetchContent_Declare(
    PLS
    SOURCE_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib/pls
    GIT_REPOSITORY https://github.com/tjhladish/pls.git
  )
  FetchContent_MakeAvailable(PLS)
else()
  message(STATUS "PLS found: ${PLS_VERSION}")
endif()

add_subdirectory(lib/CCRC32 ${CMAKE_BINARY_DIR}/CCRC32)
add_subdirectory(lib/sqdb ${CMAKE_BINARY_DIR}/sqdb)

# add_subdirectory(lib/jsoncpp ${CMAKE_BINARY_DIR}/jsoncpp)
# this is a workaround for the apparent lack of a jsoncpp library target
add_library(jsoncpp STATIC lib/jsoncpp/src/lib_json/json_reader.cpp lib/jsoncpp/src/lib_json/json_value.cpp lib/jsoncpp/src/lib_json/json_writer.cpp)
target_include_directories(jsoncpp PUBLIC lib/jsoncpp/include lib/jsoncpp/src/lib_json)

add_library(abcsmc STATIC
    $<TARGET_OBJECTS:CCRC32> $<TARGET_OBJECTS:sqdb> $<TARGET_OBJECTS:sqlite3>
    src/AbcLog.cpp src/AbcMPI.cpp src/AbcSmc.cpp src/AbcUtil.cpp
)

set_target_properties(abcsmc PROPERTIES
  VERSION ${${CMAKE_PROJECT_NAME}_VERSION}
  SOVERSION ${${CMAKE_PROJECT_NAME}_VERSION_MAJOR})


find_package(GSL REQUIRED) # the abc library objects depend on the GSL library
target_link_libraries(abcsmc PLS::PLS sqdb jsoncpp m GSL::gsl GSL::gslcblas ${CMAKE_DL_LIBS})
# the abc library objects depend on pls, sqdb, jsoncpp, m (aka math.h), GSL, and dl (aka dl.h) libraries

# look in these directories for header files, also make them available to other projects that link to abc
target_include_directories(abcsmc PRIVATE include src lib)

add_executable(abc-bin src/main.cpp)
target_link_libraries(abc-bin abcsmc)
set_target_properties(abc-bin PROPERTIES OUTPUT_NAME abc)

include(CMakePackageConfigHelpers)
write_basic_package_version_file(
  "AbcSmcConfigVersion.cmake"
  VERSION ${AbcSmc_VERSION}
  COMPATIBILITY SameMajorVersion
)

message(CMAKE_CURRENT_BINARY_DIR="${CMAKE_CURRENT_BINARY_DIR}")

install(FILES "AbcSmcConfig.cmake" "${CMAKE_CURRENT_BINARY_DIR}/AbcSmcConfigVersion.cmake"
  DESTINATION lib/cmake/AbcSmc
)

install(TARGETS abc-bin
  RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
)

install(TARGETS abcsmc # EXPORT AbcSmcTargets
  LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR}
  ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
  INCLUDES DESTINATION include
)

install(DIRECTORY include/AbcSmc
  DESTINATION include
)

# install(EXPORT AbcSmcTargets
#   NAMESPACE AbcSmc::
#   DESTINATION lib/cmake/AbcSmc
# )