cmake_minimum_required(VERSION 3.10)

# set the project name
project(GTasb3D VERSION 1.0 LANGUAGES CXX)

# specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)
if(NOT CMAKE_BUILD_TYPE)
	set(CMAKE_BUILD_TYPE Release)
endif()
#set(CMAKE_CXX_FLAGS "-std=c++-11 -O3")
set(CAMEK_CXX_FLAGS "-Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-g")
set(CMAKE_CXX_FLAGS_RELEASE "-O3")
#set(CMAKE_CXX_FLAGS_RELEASE "-O3 -I ~/include/armadillo/include/ -L ~/lib/lapack/ -L ~/include/armadillo/lib64/")

# add the executable
file (GLOB_RECURSE GTA_SOURCES CONFIGURE_DEPENDS "src/*.cpp")
add_executable(GTasb3D ${GTA_SOURCES})

# find package
find_package(OpenMP REQUIRED)
find_package(Armadillo REQUIRED)
if(OpenMP_CXX_FOUND)
	target_link_libraries(GTasb3D PRIVATE OpenMP::OpenMP_CXX ${ARMADILLO_LIBRARIES})
endif()
#find_package(Armadillo REQUIRED)
include_directories(${ARMADILLO_INCLUDE_DIRS})

# add directory to the list of paths to search for include files
target_include_directories(GTasb3D PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/include)


