cmake_minimum_required(VERSION 3.6)

project(Serial_Mandelbrot)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z -O3 -Wall -Wextra -Werror")

set(SOURCE_FILES
	./src/main.cpp
)

set( CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR} )

add_executable(run ${SOURCE_FILES})
