﻿cmake_minimum_required (VERSION 3.2)

project ("LAB2")

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS ${CMAKE_CXX_FLAGS} " -O0")
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS}")

include_directories(./src ./include)

file(GLOB SOURCES *.cpp *.h)

set(SOURCES ./src/main.cpp ./src/shared.cpp src/thermalConductivity.cpp)
set(HEADERS ./include/SOLE.h ./include/shared.h include/thermalConductivity.h ./include/SOLE.h)

add_executable(${PROJECT_NAME} ${SOURCES} ${HEADERS})
