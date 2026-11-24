add_library(PolySolverStatic STATIC
    src/PolySolver.cpp
)

target_include_directories(PolySolverStatic PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
)

target_compile_features(PolySolverStatic PUBLIC cxx_std_20)