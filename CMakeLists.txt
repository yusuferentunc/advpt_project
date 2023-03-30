cmake_minimum_required(VERSION 3.14)

project(flow VERSION 1.0
              DESCRIPTION "Optical Flow"
              LANGUAGES CXX)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/data
        DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

#add_executable(test_tensor test_tensor.cpp)
#target_compile_features(test_tensor PRIVATE cxx_std_20)
#target_compile_options(test_tensor PRIVATE -Wall -Wextra -pedantic -Werror -pg -g)
#target_link_options(test_tensor PRIVATE -pg)

#add_executable(test_matvec test_matvec.cpp)
#target_compile_features(test_matvec PRIVATE cxx_std_20)
#target_compile_options(test_matvec PRIVATE -Wall -Wextra -pedantic -Werror -pg -g)
#target_link_options(test_matvec PRIVATE -pg)

add_executable(flow src/OpticalFlow.cpp src/mg.cpp)
target_compile_features(flow PRIVATE cxx_std_20)
target_compile_options(flow PRIVATE -O3 -fopenmp -march=native -fconcepts -pedantic -Wall -Werror -Wextra)
target_link_options(flow PRIVATE)

find_package(OpenMP)
if(OpenMP_CXX_FOUND)
        target_link_libraries(flow PUBLIC OpenMP::OpenMP_CXX)
endif()
