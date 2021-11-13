cmake_minimum_required(VERSION 3.10)

project(MarsHydrogenEscape)

set(CMAKE_CXX_STANDARD 11)
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

# dependencies
add_subdirectory(pybind11)
add_subdirectory(src/libsoda)
add_library(interp ${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/pwl_interp_2d_scattered.cpp
${CMAKE_CURRENT_SOURCE_DIR}/src/interpolation/r8lib.cpp)

# Mars hydrogen escape
add_library(Mars src/MarsHydrogenEscape.cpp)
target_link_libraries(Mars PUBLIC lsoda)
target_link_libraries(Mars PUBLIC interp)

# main
add_executable(main src/main.cpp)
target_link_libraries(main Mars)

# Python module
pybind11_add_module(MarsPy src/wrapper.cpp)
set(CMAKE_POSITION_INDEPENDENT_CODE ON)
target_link_libraries(MarsPy PUBLIC Mars)
install(TARGETS MarsPy DESTINATION ${CMAKE_CURRENT_SOURCE_DIR})