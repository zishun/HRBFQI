#project ( CXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive -fPIC")

project(DataStructure)

file(GLOB DIR_SRCS
./*)

add_library(DataStructure ${DIR_SRCS})
