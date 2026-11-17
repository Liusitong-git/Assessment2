cmake_minimum_required(VERSION 3.8)
project(student_auto_aim)

if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 17)
  set(CMAKE_CXX_STANDARD_REQUIRED ON)
endif()
if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

find_package(ament_cmake REQUIRED)
find_package(rclcpp REQUIRED)
find_package(rm_aim_sim_bridge REQUIRED)

add_executable(aim_node src/aim_node.cpp)
ament_target_dependencies(aim_node rclcpp rm_aim_sim_bridge)

install(TARGETS aim_node DESTINATION lib/${PROJECT_NAME})

ament_package()
