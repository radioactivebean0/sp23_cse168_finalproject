cmake_minimum_required(VERSION 3.10)

project(lajolla)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} "${CMAKE_SOURCE_DIR}/cmake/")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED True)

set(CMAKE_BUILD_TYPE RelWithDebInfo)

if(MSVC)
  add_compile_options(/Wall)
else()
  add_compile_options(-Wall -mavx)
endif()

include_directories(${CMAKE_SOURCE_DIR}/src)

set(SRCS src/3rdparty/miniz.h
         src/3rdparty/miniz.c
         src/3rdparty/pugiconfig.hpp
         src/3rdparty/pugixml.hpp
         src/3rdparty/pugixml.cpp
         src/3rdparty/pcg.h
         src/3rdparty/stb_image.h
         src/3rdparty/tinyexr.h
         src/3rdparty/tinyply.h
         src/sampling.h
         src/compute_normals.h
         src/ortho_basis.h
         src/camera.h
         src/camera.cpp
         src/materials.h
         src/materials.cpp
         src/lights.h
         src/shapes.h
         src/shapes.cpp
         src/flexception.h
         src/image.h
         src/matrix.h
         src/parallel.h
         src/parse_obj.h
         src/parse_ply.h
         src/parse_scene.h
         src/parse_serialized.h
         src/print_scene.h
         src/torrey.h
         src/transform.h
         src/vector.h
         src/compute_normals.cpp
         src/bvh.cpp
         src/bvh.h
         src/scene.h
         src/scene.cpp
         src/image.cpp
         src/mis.h
         src/mis.cpp
         src/next_event.h
         src/next_event.cpp
         src/path_trace.h
         src/path_trace.cpp
         src/render.h
         src/render.cpp
         src/main.cpp
         src/parallel.cpp
         src/parse_obj.cpp
         src/parse_ply.cpp
         src/parse_scene.cpp
         src/parse_serialized.cpp
         src/print_scene.cpp
         src/transform.cpp)

add_library(torrey_lib STATIC ${SRCS})
add_executable(torrey src/main.cpp)

set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)
target_link_libraries(torrey Threads::Threads)

target_link_libraries(torrey torrey_lib)
