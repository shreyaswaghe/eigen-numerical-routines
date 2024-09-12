cmake_minimum_required(VERSION 3.20.0)

project(CellSimulatorNumerics)
add_library(${PROJECT_NAME} STATIC)

set_target_properties(
    ${PROJECT_NAME}
    PROPERTIES
    CXX_STANDARD
    17
)

target_include_directories(
    ${PROJECT_NAME}
    PRIVATE
    /usr/local/include
)


add_subdirectory(odeintegrator)
file(GLOB odeintegrator_headers odeintegrator/*.h*)

add_subdirectory(montecarlo)
file(GLOB montecarlo_headers montecarlo/*.h*)

add_subdirectory(stochastic)
file(GLOB stochastic_headers stochastic/*.h*)

add_subdirectory(quadrature)
file(GLOB quadrature_headers quadrature/*.h*)

install(
    TARGETS ${PROJECT_NAME}
    LIBRARY DESTINATION lib
    PUBLIC_HEADER DESTINATION include
)


install(
    FILES 
    ${odeintegrator_headers}
    ${montecarlo_headers}
    ${stochastic_headers}
    ${quadrature_headers}
    DESTINATION include
)
