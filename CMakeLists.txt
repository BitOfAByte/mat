PROJECT(mat)
CMAKE_MINIMUM_REQUIRED(VERSION 2.8)

IF(CMAKE_COMPILER_IS_GNUCC)
 SET(CMAKE_CXX_FLAGS "-Wall -pedantic")
ENDIF()

INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include/)

ADD_EXECUTABLE(mat test.cpp)

