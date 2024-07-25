cmake_minimum_required(VERSION 3.30)
project(Dynamic_SCC)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_CXX_COMPILER g++-14)
set(CMAKE_C_COMPIELER gcc-14)

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization mpi)

# source files
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(dynamic_scc ${SOURCES})

target_compile_options(dynamic_scc PRIVATE "-fopenmp")

target_link_libraries(dynamic_scc PRIVATE ${MPI_C_LIBRARIES} ${Boost_LIBRARIES} "/opt/homebrew/opt/libomp/lib/libomp.dylib")

#include directories
target_include_directories(dynamic_scc PRIVATE ${MPI_C_INCLUDE_PATH} ${Boost_INCLUDE_DIRS} "/opt/homebrew/opt/libomp/include" ${CMAKE_CURRENT_SOURCE_DIR}/include)
