cmake_minimum_required(VERSION 3.18 FATAL_ERROR)

#set(CUDA_LIB_DIR "$ENV{CUDA_PATH}\\bin")
#set(CMAKE_CUDA_COMPILER "${CUDA_LIB_DIR}\\nvcc.exe")

if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.18)
  cmake_policy(SET CMP0104 OLD)
endif()


project(LabViewer
        VERSION 0.2.1
        LANGUAGES C CXX CUDA)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CUDA_STANDARD 17)

#add_compile_options(-fsanitize=address)

cmake_policy(SET CMP0100 NEW)

set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

#should set for specific gpu structure
# 3090 - "86"
# 4090 - "89"
set(CMAKE_CUDA_ARCHITECTURES "70" "80" "86" "89")
set(CMAKE_CUDA_ARCHITECTURES "89")

set(CMAKE_CUDA_STANDARD_REQUIRED TRUE)
set(CMAKE_CXX_STANDARD_REQUIRED TRUE)

add_definitions(-DNO_DECREMENT_DEPRECATED_WARNINGS)
#Required For OpenMesh
add_definitions(-D_USE_MATH_DEFINES)
#Required For Ambiguous "byte" between C++17 stl and windows
add_definitions(-D_HAS_STD_BYTE=0)
if (MSVC)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /bigobj")
endif(MSVC)

#Qt
set(QT_CMAKE_PATH "D:/Dev/QT/5.14.2/msvc2017_64/lib/cmake")
set(Qt5Core_DIR ${QT_CMAKE_PATH}/Qt5Core)
set(Qt5Gui_DIR ${QT_CMAKE_PATH}/Qt5Gui)
set(Qt5OpenGL_DIR ${QT_CMAKE_PATH}/Qt5OpenGL)
set(Qt5Widgets_DIR ${QT_CMAKE_PATH}/Qt5Widgets)

find_package(OpenGL REQUIRED)
include_directories(${OPENGL_INCLUDE_DIR})
find_package(Qt5Widgets REQUIRED)
find_package(Qt5Core REQUIRED)
find_package(Qt5Gui REQUIRED)
find_package(Qt5OpenGL REQUIRED)

#should set with -DMKL_INTERFACE_FULL:STRING=intel_lp64
set(MKL_INTERFACE_FULL "intel_lp64" CACHE STRING "MKL Interface" FORCE)
set(MKL_LINK "static" CACHE STRING "MKL Linking" FORCE)
find_package(MKL REQUIRED)

add_library(MKL_header_lib INTERFACE)

if(MKL_FOUND)
  include_directories(${MKL_INCLUDE_DIRS})
  include_directories($<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)

  target_link_libraries(MKL_header_lib INTERFACE ${MKL_CORE_LIBRARY})
  target_compile_options(MKL_header_lib INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_COMPILE_OPTIONS>)
  target_include_directories(MKL_header_lib INTERFACE $<TARGET_PROPERTY:MKL::MKL,INTERFACE_INCLUDE_DIRECTORIES>)
  target_link_libraries(MKL_header_lib INTERFACE $<LINK_ONLY:MKL::MKL>)
else()
  message(WARNING "MKL libs not found")
endif()

set(BOOST_DIR "extern/boost_1_76_0")
set(Boost_INCLUDE_DIR "extern/boost_1_76_0")
set(Boost_LIBRARY_DIRS "extern/boost_1_76_0/libs")
find_package(Boost REQUIRED)
IF (Boost_FOUND)
  include_directories(${Boost_INCLUDE_DIR})
  link_directories(${Boost_LIBRARY_DIRS})
ENDIF()

find_package(OpenMP REQUIRED)

include_directories(include)
include_directories(include/rapidjson/include)
include_directories(include/spdlog/include)
include_directories(extern)
include_directories(Apps)
include_directories(extern/OpenMesh/src)
include_directories(extern/freeglut/include)
link_directories(library)
link_directories(library/opencascade)

# Auto-detect GPU architecture, sets ${CUDA_ARCHS}
include("cmake/AutoDetectCudaArch.cmake")

add_subdirectory(extern)

add_library(DIR_header_lib INTERFACE)
target_compile_definitions(DIR_header_lib
        INTERFACE INPUT_DIR=${PROJECT_SOURCE_DIR}/input/
        INTERFACE OUTPUT_DIR=${PROJECT_SOURCE_DIR}/output/
        INTERFACE RESOURCE_DIR=${PROJECT_SOURCE_DIR}/resource/
        INTERFACE PROJECT_DIR=${PROJECT_SOURCE_DIR}/
)

add_definitions(-DIGA_OUTPUT_DIR="${PROJECT_SOURCE_DIR}/output/IGADataOutput/")

#SET(CMAKE_RUNTIME_OUTPUT_DIRECTORY  ${PROJECT_BINARY_DIR}/bin)
#for windows
set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/../bin)
#for linux
#set(EXECUTABLE_OUTPUT_PATH  ${PROJECT_BINARY_DIR}/bin)

add_subdirectory(Apps)