cmake_minimum_required(VERSION 3.2)

project("LAB4")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp src/lab4.cpp src/SOLE.cpp src/test.cpp)
set(HEADERS ./include/SOLE.h ./include/shared.h include/lab4.h ./include/SOLE.h ./include/test.h)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
