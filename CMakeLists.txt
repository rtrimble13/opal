cmake_minimum_required(VERSION 3.16)
project(opal VERSION 0.1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

option(OPAL_BUILD_TESTS "Build unit tests" ON)
option(OPAL_BUILD_CLI "Build the opal command line tool" ON)
option(OPAL_BUILD_PYTHON "Build the Python bindings (requires pybind11)" OFF)

# Header-only library target.
add_library(opal INTERFACE)
target_include_directories(opal INTERFACE
  $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
  $<INSTALL_INTERFACE:include>)
target_compile_features(opal INTERFACE cxx_std_17)

if(OPAL_BUILD_CLI)
  add_executable(opal_cli
    cli/main.cpp)
  target_link_libraries(opal_cli PRIVATE opal)
  set_target_properties(opal_cli PROPERTIES OUTPUT_NAME opal)
  if(MSVC)
    target_compile_options(opal_cli PRIVATE /W4)
  else()
    target_compile_options(opal_cli PRIVATE -Wall -Wextra)
  endif()
  install(TARGETS opal_cli RUNTIME DESTINATION bin)
endif()

if(OPAL_BUILD_TESTS)
  enable_testing()
  add_executable(opal_tests
    tests/test_main.cpp
    tests/test_analytic.cpp
    tests/test_engines.cpp
    tests/test_rates.cpp)
  target_link_libraries(opal_tests PRIVATE opal)
  target_include_directories(opal_tests PRIVATE tests)
  if(NOT MSVC)
    target_compile_options(opal_tests PRIVATE -Wall -Wextra)
  endif()
  add_test(NAME opal_tests COMMAND opal_tests)
endif()

if(OPAL_BUILD_PYTHON)
  find_package(pybind11 CONFIG REQUIRED)
  pybind11_add_module(_opal python/bindings.cpp)
  target_link_libraries(_opal PRIVATE opal)
  # Place the extension inside the Python package so `import opal` works
  # with PYTHONPATH=python, matching the setup.py in-place build layout.
  set_target_properties(_opal PROPERTIES
    LIBRARY_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/python/opal)
endif()

install(DIRECTORY include/opal DESTINATION include)
