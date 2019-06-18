set(PROJECT_NAME polygonizer)
project(${PROJECT_NAME})
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive -fPIC")

file(GLOB DIR_SRCS
./*)

add_library(${PROJECT_NAME} ${DIR_SRCS})
target_link_libraries(${PROJECT_NAME} PQP)
