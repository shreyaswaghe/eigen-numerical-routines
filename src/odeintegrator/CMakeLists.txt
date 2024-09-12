file(GLOB odeintegrator_src *.cpp)
file(GLOB odeintegrator_headers *.hpp)

target_sources( 
    ${PROJECT_NAME}
    PRIVATE
    ${odeintegrator_src}
    ${odeintegrator_headers}
)
