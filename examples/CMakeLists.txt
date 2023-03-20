
project( AbcSmcExamples )

# the example specific includes
include_directories(include)

# the demos all depend on abc

add_library(dice SHARED src/dice.cpp)
target_link_libraries(dice GSL::gsl GSL::gslcblas)

add_executable(dice_game src/dice_game.cpp)
target_link_libraries(dice_game GSL::gsl GSL::gslcblas)

# TODO: not yet working
# add_executable(direct direct/main.cpp)
# target_link_libraries(direct abc)

foreach(demo shared executable)
add_custom_command(OUTPUT ${demo}_config.json
    COMMAND gojq -s add ${CMAKE_CURRENT_SOURCE_DIR}/reference.json ${CMAKE_CURRENT_SOURCE_DIR}/${demo}/partial.json > ${demo}_config.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/reference.json ${CMAKE_CURRENT_SOURCE_DIR}/${demo}/partial.json
    VERBATIM
)
endforeach()

add_custom_command(OUTPUT pseudo_config.json
    COMMAND cp ${CMAKE_CURRENT_SOURCE_DIR}/pseudo/config.json pseudo_config.json
    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/pseudo/config.json
    VERBATIM
)

# all the demos depend on abc
foreach(demo shared integral executable pseudo)
    add_executable(${demo} ${demo}/main.cpp)
    target_link_libraries(${demo} abc)
#    add_dependencies(${demo} ${demo}/config.json)
endforeach()

# ...some depend on the shared library o
foreach(shareddemo shared pseudo)
    target_link_libraries(${shareddemo} dice)
endforeach()

# ...some depend on the executable
foreach(execdemo executable)
    add_dependencies(${execdemo} dice_game)
endforeach()

foreach(demo shared pseudo executable)
    add_custom_target(${demo}_config DEPENDS ${demo}_config.json)
    add_dependencies(${demo} ${demo}_config)
endforeach()

# TODO this must be foreach-able?