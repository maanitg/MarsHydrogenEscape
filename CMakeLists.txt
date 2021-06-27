cmake_minimum_required(VERSION 3.4...3.18)

project(MarsHydrogenEscape)
add_subdirectory(pybind11)
pybind11_add_module(Mars src/MarsHydrogenEscape.cpp src/wrapper.cpp)

add_subdirectory(src/libsoda)
target_link_libraries(Mars PUBLIC lsoda)
target_link_libraries(Mars PUBLIC interp)

add_library(interp ${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/pwl_interp_2d_scattered.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/r8lib.cpp)
                   
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_CXX_FLAGS "-Wall")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")

install(TARGETS Mars DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})