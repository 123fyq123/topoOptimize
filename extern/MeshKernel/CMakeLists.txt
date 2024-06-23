set(LIBRARY_NAME MeshKernel)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES *.hpp)
#message(${SOURCE_FILES})

# 添加 MeshKernel 模块的源文件
set(MESH_KERNEL_SOURCES
        IO.cpp
        # 其他源文件
)

add_library(${LIBRARY_NAME} STATIC ${INCLUDE_FILES} ${SOURCE_FILES} ${MESH_KERNEL_SOURCES})

# 设置 OpenCASCADE 的路径
set(OpenCASCADE_DIR "D:/Git/LabViewer/extern/opencascade/opencascade770-install/cmake")

# 查找 OpenCASCADE 组件
find_package(OpenCASCADE COMPONENTS FoundationClasses REQUIRED)
find_package(OpenCASCADE COMPONENTS ModelingAlgorithms REQUIRED)

# 链接 OpenCASCADE 库
target_include_directories(${LIBRARY_NAME} PRIVATE ${OpenCASCADE_INCLUDE_DIR})
target_link_libraries(${LIBRARY_NAME} PRIVATE ${OpenCASCADE_LIBRARIES})
