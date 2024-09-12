cmake_minimum_required(VERSION 3.20.0)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

project(test_numerics)

add_executable(
    ${PROJECT_NAME}
    main.cpp
)

add_library(
    numeric_backend
    STATIC IMPORTED
)

set_target_properties(
    numeric_backend
    PROPERTIES IMPORTED_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/lib/libCellSimulatorNumerics.a
)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD 17
)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
    include
    PRIVATE
    /usr/local/include
)

target_link_libraries(
    ${PROJECT_NAME}
    numeric_backend
)

install(
    TARGETS ${PROJECT_NAME}
    DESTINATION ${CMAKE_CURRENT_SOURCE_DIR}/bin
)
