cmake_minimum_required(VERSION 3.5)
project(iris_2d)

# Default to C99
if(NOT CMAKE_C_STANDARD)
  set(CMAKE_C_STANDARD 99)
endif()

# Default to C++14
if(NOT CMAKE_CXX_STANDARD)
  set(CMAKE_CXX_STANDARD 14)
endif()

if(CMAKE_COMPILER_IS_GNUCXX OR CMAKE_CXX_COMPILER_ID MATCHES "Clang")
  add_compile_options(-Wall -Wextra -Wpedantic)
endif()

# find dependencies
find_package(ament_cmake REQUIRED)
find_package(autodiff REQUIRED)
find_package(ifopt REQUIRED)
find_package(Eigen3 REQUIRED)

if(BUILD_TESTING)
  find_package(ament_lint_auto REQUIRED)
  # the following line skips the linter which checks for copyrights
  # uncomment the line when a copyright and license is not present in all source files
  #set(ament_cmake_copyright_FOUND TRUE)
  # the following line skips cpplint (only works in a git repo)
  # uncomment the line when this package is not in a git repo
  #set(ament_cmake_cpplint_FOUND TRUE)
  ament_lint_auto_find_test_dependencies()
endif()

set(dependencies
  autodiff
  ifopt
)

add_library(${PROJECT_NAME} SHARED 
  src/iris_2d.cpp
)
target_include_directories(${PROJECT_NAME}
  PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>
  ${EIGEN3_INCLUDE_DIR}
)
ament_target_dependencies(${PROJECT_NAME} 
  ${dependencies}
)
target_link_libraries(${PROJECT_NAME}
  autodiff::autodiff
)

install(
  TARGETS ${PROJECT_NAME}
  EXPORT ${PROJECT_NAME}_export
  LIBRARY DESTINATION lib
  INCLUDES DESTINATION include
  PUBLIC_HEADER DESTINATION include
)

install(
  EXPORT  ${PROJECT_NAME}_export
  FILE    ${PROJECT_NAME}-config.cmake
  DESTINATION share/cmake/${PROJECT_NAME}
  EXPORT_LINK_INTERFACE_LIBRARIES
)
install(
  DIRECTORY include/
  DESTINATION include
)

ament_export_include_directories(include ${EIGEN3_INCLUDE_DIR})
ament_export_dependencies(${dependencies})
ament_export_libraries(
  ${PROJECT_NAME}
)

function(add_iris_2d_test base_name)
	add_executable(${base_name} test/${base_name}.cpp)
	target_link_libraries(${base_name} ${PROJECT_NAME})
    target_include_directories(${base_name}
    PUBLIC
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<INSTALL_INTERFACE:include>
    )
endfunction()

add_iris_2d_test(test_ie_solver)
add_iris_2d_test(test_ccp_solver)
add_iris_2d_test(test_iris)

ament_package()
