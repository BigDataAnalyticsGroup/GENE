cmake_minimum_required(VERSION 3.1)

project("GENE" C CXX)

# Build flags
SET(EXECUTABLE_OUTPUT_PATH      "${PROJECT_BINARY_DIR}/bin")
SET(LIBRARY_OUTPUT_PATH         "${PROJECT_BINARY_DIR}/lib")

# Compile flags
SET(CMAKE_COMPILE_FLAGS         "-W -Wall -pedantic -march=native -Wno-gnu-anonymous-struct -Wno-nested-anon-types -Wno-c99-extensions -Wno-vla-extension")
SET(CMAKE_C_FLAGS               "${CMAKE_C_FLAGS} ${CMAKE_COMPILE_FLAGS}")
SET(CMAKE_CXX_FLAGS             "-std=c++17 ${CMAKE_CXX_FLAGS} ${CMAKE_COMPILE_FLAGS} -pthread")
SET(CMAKE_CXX_FLAGS_DEBUG       "${CMAKE_CXX_FLAGS_DEBUG} -ggdb3 -fsanitize=leak")

# Include Directories
include_directories(src/)
include_directories(third_party/tlx/)
include_directories(third_party/ART/include)
include_directories(third_party/PGM-index/include)
include_directories(third_party/robin-map/include/)
include_directories(third_party/cxxopts/include/)
include_directories(third_party/blog/src)

# Subdirectories
add_subdirectory(src/)
