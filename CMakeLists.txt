cmake_minimum_required(VERSION 3.10)
project(iopos)

find_package(OpenCV REQUIRED)
include_directories(${OpenCV_INCLUDE_DIRS})

include_directories($ENV{CMAKE_INCLUDE_PATH})

add_executable(test_imu_odo test_iopos.cpp so3.c matrix.c quaternion.c iopos.c)
target_link_libraries(test_imu_odo ${OpenCV_LIBS})