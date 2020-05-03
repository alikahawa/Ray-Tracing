cmake_minimum_required(VERSION 3.10)
project(cg_35)

set(CMAKE_CXX_STANDARD 11)

include_directories(Code/extern)
include_directories(Code/include)

add_executable(cg_35
        Code/data/cube.mtl
        Code/data/cube.obj
        Code/data/dodgeColorTest.mtl
        Code/data/dodgeColorTest.obj
        Code/extern/stb_image_write.h
        Code/include/imageWriter.h
        Code/include/matrix.h
        Code/include/mesh.h
        Code/include/mouse.h
        Code/include/raytracing.h
        Code/include/Vec3D.h
        Code/include/Vertex.h
        Code/src/main.cpp
        Code/src/mesh.cpp
        Code/src/raytracing.cpp
        Code/CMakeLists.txt
        Code/result.bmp
        "Screenshots/basic raytracing.bmp"
        "Screenshots/diffuse lighting.bmp"
        README.md Code/include/concurrentqueue.h)
