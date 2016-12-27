cmake_minimum_required(VERSION 3.6)
project(Material_assignment)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fpermissive -std=c++11 -lga")

set(SOURCE_FILES
        heuristic.cpp
    heuristic.h
    main.cpp
    solveAll.sh
    utils.h heuristic/CoIoTeGenome.cpp heuristic/CoIoTeGenome.h)

add_executable(Material_assignment ${SOURCE_FILES})
