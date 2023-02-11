cmake_minimum_required (VERSION 3.2)

project ("LAB5")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp  ./src/SOLE.cpp ./src/Source.cpp ./src/Source1.cpp)
set(HEADERS ./include/shared.h  ./include/SOLE.h)

#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} ./src/main.cpp)

#install(TARGETS ${PROJECT_NAME} DESTINATION "LAB1/bin")