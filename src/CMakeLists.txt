project(DYNAMICVORONOI)
cmake_minimum_required(VERSION 2.6)

#set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)

SET( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../bin )
SET( CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib )
SET( CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_SOURCE_DIR}/../lib )


# COMPILER FLAGS
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O2")
set(CMAKE_CXX_FLAGS_DEBUG "-O0 -g -ggdb")

add_library(dynamicvoronoi dynamicvoronoi.cpp)

add_executable(example example.cpp)
target_link_libraries(example dynamicvoronoi)

