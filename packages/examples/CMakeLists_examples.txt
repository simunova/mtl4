project(mtl_examples)

cmake_minimum_required(VERSION 2.8)

find_package(MTL REQUIRED)
include_directories(${MTL_INCLUDE_DIRS})
add_definitions(${MTL_CXX_DEFINITIONS})

file(GLOB SOURCES ${MTL_DIR}/examples/*.cpp)

foreach (source IN LISTS SOURCES)
get_filename_component (test ${source} NAME_WE)		# filename without extension
add_executable (${test} ${source} ${MTL_INCLUDE_DIRS})
endforeach (source IN LISTS SOURCES)
