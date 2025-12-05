cmake_minimum_required(VERSION 3.10)

# Назва проєкту
project(CGAL_Star_Triangulation)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release)
endif()

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)

find_package(CGAL REQUIRED COMPONENTS Qt5)
find_package(Qt5 COMPONENTS Core Gui Widgets OpenGL REQUIRED)

include(${CGAL_USE_FILE})

add_executable(main_app main.cpp)

target_compile_definitions(main_app PRIVATE CGAL_USE_BASIC_VIEWER)

target_link_libraries(main_app 
    PRIVATE 
    CGAL::CGAL 
    CGAL::CGAL_Qt5 
    Qt5::Core 
    Qt5::Gui 
    Qt5::Widgets 
    Qt5::OpenGL
)

add_executable(gen_dodeca gen_dodeca.cpp)
target_link_libraries(gen_dodeca PRIVATE CGAL::CGAL)