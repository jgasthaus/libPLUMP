include_directories("..")
file(GLOB SOURCES "*.cc")
add_library(plump SHARED ${SOURCES})