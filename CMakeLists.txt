cmake_minimum_required(VERSION 3.14)
project(example)
set(CMAKE_CXX_STANDARD 17)
add_subdirectory(pybind11)
add_subdirectory(matplotplusplus)

include_directories("C:/Users/pc/OneDrive/Pulpit/WYGRANE_podejcie_4_projekt_3/pybind11/include")
include_directories("C:/Program Files/Python313")
include_directories("C:/Users/pc/OneDrive/Pulpit/WYGRANE_podejcie_4_projekt_3/matplotplusplus/source/matplot")
pybind11_add_module(example example.cpp)
target_link_libraries(example PRIVATE matplot)





