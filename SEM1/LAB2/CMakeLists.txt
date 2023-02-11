cmake_minimum_required (VERSION 3.2)

project ("LAB2")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS})
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}" )

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp  ./src/SOLE.cpp ./src/threeToNine.cpp)
set(HEADERS ./include/shared.h  ./include/SOLE.h ./include/threeToNine.h)

#find_package(GTest REQUIRED)
#include_directories(${GTEST_INCLUDE_DIRS})

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS} "include/threeToNine.h" "src/threeToNine.cpp")

#install(TARGETS ${PROJECT_NAME} DESTINATION "LAB1/bin")