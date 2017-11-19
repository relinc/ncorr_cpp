CMAKE_MINIMUM_REQUIRED(VERSION 3.2)
PROJECT(ncorr_library)

# Only tested for g++ on Ubuntu 12.04. This assumes all required libraries have been 
# installed, so directories to dependent libraries and their headers are not explicitly 
# included, since the install directories are searched automatically by g++.

# Set files
file(GLOB ncorr_SRC
"${CMAKE_CURRENT_SOURCE_DIR}/src/ncorr/*.cpp"
)
file(GLOB ncorr_H
"${CMAKE_CURRENT_SOURCE_DIR}/include/*.h"
)

SET(ncorr_src ${ncorr_SRC})
SET(ncorr_h ${ncorr_H})

# Set include directory
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/SuiteSparse/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/fftw/include)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/opencv/include)
SET(SuiteSparse_DIR ${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/SuiteSparse)
LIST(APPEND CMAKE_MODULE_PATH ${CMAKE_CURRENT_SOURCE_DIR})
FIND_PACKAGE(SuiteSparse)

# Set output for library
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/opencv/lib)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/fftw/lib)
LINK_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/ThirdPartyLibraries/BLAS-LAPACK)


# Add library
ADD_LIBRARY(ncorr STATIC ${ncorr_src})
find_package(SuiteSparse COMPONENTS metis suitesparse)
ADD_EXECUTABLE(ncorr_commandLine src/ncorr_commandLine)
target_link_libraries (ncorr_commandLine ncorr)
target_link_libraries(ncorr ${SuiteSparse_LIBRARIES})
target_link_libraries(ncorr opencv_world331)
target_link_libraries(ncorr libblas liblapack)
target_link_libraries(ncorr ${SuiteSparse_LAPACK_BLAS_DLL})
target_link_libraries(ncorr libfftw3-3)
# target_link_libraries(ncorr -lfftw3-3f)
# target_link_libraries(ncorr -lfftw3-3l)


# Set C++11 support
set_property(TARGET ncorr PROPERTY CXX_STANDARD 11)
set_property(TARGET ncorr PROPERTY CXX_STANDARD_REQUIRED ON)
# Set -03 optimization
INCLUDE(CheckCXXCompilerFlag)
CHECK_CXX_COMPILER_FLAG("-O3" COMPILER_SUPPORTS_O3)
if (COMPILER_SUPPORTS_O3)
	SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3")
endif()

SET(TARGET_DIR ${CMAKE_CURRENT_SOURCE_DIR}/lib)
# Disable debugging
ADD_DEFINITIONS(-DNDEBUG)

# Install library
INSTALL(TARGETS ncorr_commandLine RUNTIME  DESTINATION ${TARGET_DIR})

add_custom_command(TARGET ncorr_commandLine POST_BUILD        # Adds a post-build event to MyTest
	COMMAND ${CMAKE_COMMAND} -E copy_directory  # which executes "cmake - E copy_if_different..."
	"${PROJECT_SOURCE_DIR}/ThirdPartyLibraries/BLAS-LAPACK/"      # <--this is in-file
	"${TARGET_DIR}")                 # <--this is out-file path
add_custom_command(TARGET ncorr_commandLine POST_BUILD        # Adds a post-build event to MyTest
	COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
	"${PROJECT_SOURCE_DIR}/ThirdPartyLibraries/opencv/lib/opencv_world331.dll"      # <--this is in-file
	"${TARGET_DIR}")       
add_custom_command(TARGET ncorr_commandLine POST_BUILD        # Adds a post-build event to MyTest
	COMMAND ${CMAKE_COMMAND} -E copy_if_different  # which executes "cmake - E copy_if_different..."
	"${PROJECT_SOURCE_DIR}/ThirdPartyLibraries/fftw/lib/libfftw3-3.dll"      # <--this is in-file
	"${TARGET_DIR}") 