
file(GLOB quadrature_src ./*.cpp)
file(GLOB quadrature_headers ./*.hpp)

target_sources(
    ${PROJECT_NAME}
    PRIVATE
    ${quadrature_src}
    ${quadrature_headers}
)
