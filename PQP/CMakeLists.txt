set(PROJECT_NAME PQP)
project(${PROJECT_NAME})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive")

file(GLOB DIR_SRCS
./include/*
./src/*)

add_library(${PROJECT_NAME} ${DIR_SRCS})

