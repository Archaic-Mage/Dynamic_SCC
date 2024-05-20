cmake_minimum_required(VERSION 3.5)
project(Dynamic_SCC)

# set(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fopenmp")

find_package(MPI REQUIRED)
find_package(Boost REQUIRED COMPONENTS serialization mpi)

# source files
file(GLOB SOURCES ${CMAKE_CURRENT_SOURCE_DIR}/src/*.cpp)
add_executable(dynamic_scc ${SOURCES})

#include boost libraries
# set(boost_serialization /lfs/usrhome/btech/cs20b073/usr/lib/libboost_serialization.a)
# set(boost_mpi /lfs/usrhome/btech/cs20b073/usr/lib/libboost_mpi.a)

# target_link_libraries(dynamic_scc PRIVATE ${MPI_C_LIBRARIES} ${boost_serialization} ${boost_mpi})
target_link_libraries(dynamic_scc PRIVATE ${MPI_C_LIBRARIES} ${Boost_LIBRARIES})
target_compile_options(dynamic_scc PRIVATE "-std=c++17" "-fopenmp")


#include directories
target_include_directories(dynamic_scc PRIVATE ${MPI_C_INCLUDE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/include /lfs/usrhome/btech/cs20b073/usr/include)
