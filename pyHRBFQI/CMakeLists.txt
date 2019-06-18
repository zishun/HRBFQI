cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME pyHRBFQI)

project (${PROJECT_NAME} CXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive")
add_definitions(-DTHREADS_NUM=4)

find_package(pybind11 REQUIRED)
FIND_PACKAGE(OpenMP REQUIRED)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

add_subdirectory(../DataStructure ./DataStructure)
add_subdirectory(../PQP ./PQP)
add_subdirectory(../polygonizer ./polygonizer)

SET(SRC_FILES
./pyHRBFQI.cpp 
../HRBFQI/HRBF.cpp 
../HRBFQI/MeshCleaner.cpp 
)

pybind11_add_module(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} PRIVATE DataStructure polygonizer)

