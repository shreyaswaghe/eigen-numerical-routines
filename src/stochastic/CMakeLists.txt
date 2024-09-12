file(GLOB stochastic_src *.cpp)
file(GLOB stochastic_headers *.hpp)

target_sources(
    ${PROJECT_NAME}
    PRIVATE
    ${stochastic_src}
    ${stochastic_headers}
)

message(test2=${stochastic_headers})
