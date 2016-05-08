cmake_minimum_required (VERSION 2.8)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_LIST_DIR}/cmake")
project (libplump)

find_package(Boost 1.35.0 REQUIRED COMPONENTS program_options serialization iostreams filesystem system)

find_package(GSL REQUIRED)

include_directories ("${PROJECT_SOURCE_DIR}/src")
include_directories ("${PROJECT_SOURCE_DIR}/src/utils")
add_subdirectory(src/libplump)

add_executable(score_file src/utils/score_file.cc)
target_link_libraries(score_file plump ${Boost_LIBRARIES} ${GSL_LIBRARIES})