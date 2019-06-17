cmake_minimum_required(VERSION 2.8)

set(PROJECT_NAME HRBFQI)

project (${PROJECT_NAME} CXX)
SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11 -fpermissive")
add_definitions(-DTHREADS_NUM=4)

FIND_PACKAGE(OpenMP REQUIRED)
IF(OPENMP_FOUND)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
ENDIF()

add_subdirectory(DataStructure)
add_subdirectory(PQP)
add_subdirectory(polygonizer)

SET(SRC_FILES
HRBFQI/HRBFQI.cpp 
HRBFQI/FileManager.cpp 
HRBFQI/HRBF.cpp 
HRBFQI/MeshCleaner.cpp 
)

add_executable(${PROJECT_NAME} ${SRC_FILES})
target_link_libraries(${PROJECT_NAME} DataStructure polygonizer)
