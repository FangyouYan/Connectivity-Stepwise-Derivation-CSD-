cmake_minimum_required(VERSION 3.28)
project(FullStepMatrixCode)

set(CMAKE_CXX_STANDARD 17)

# Pack the common ones into one
add_library(msutils
        common_utils/MsUtils.h
        common_utils/MsUtils.cpp
)
add_executable(getMSFbyCSD main/getMSFbyCSD.cpp)
add_executable(getMSFbyFloydWarshall main/getMSFbyFloydWarshall.cpp)
target_link_libraries(getMSFbyCSD msutils)
target_link_libraries(getMSFbyFloydWarshall msutils)
