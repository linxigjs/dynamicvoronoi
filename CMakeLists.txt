project(dynamic_voronoi)
cmake_minimum_required(VERSION 2.6)

#set(CMAKE_BUILD_TYPE Debug)
SET(CMAKE_BUILD_TYPE Release)

# COMPILER FLAGS
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")

include_directories(include)

add_library(dynamicvoronoilib src/dynamicvoronoi.cpp
        include/bucketedqueue.h
        include/bucketedqueue.hxx
        include/dynamicvoronoi.h
        include/point.h)

add_executable(voronoi_example src/example.cpp)
target_link_libraries(voronoi_example dynamicvoronoilib)

