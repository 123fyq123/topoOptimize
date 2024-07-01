set(APP_NAME iGamePoster)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB SUB_SOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/source/*.cpp)

set(QRC_SOURCE_FILES Poster.qrc)
qt5_add_resources(QRC_FILES ${QRC_SOURCE_FILES})
if (${CMAKE_GENERATOR} STREQUAL "Ninja")
    set(UI_SOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/source/*.ui)
else()
    set(UI_SOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/source/aboutDialog.ui ${CMAKE_CURRENT_LIST_DIR}/source/helpDoc.ui ${CMAKE_CURRENT_LIST_DIR}/source/iGame-Poster.ui)
endif ()
message("UI FILES ARE ${UI_SOURCE_FILES}")

qt5_wrap_ui(UI_FILES ${UI_SOURCE_FILES})

add_executable(${APP_NAME} ${SUB_SOURCE_FILES} ${SOURCE_FILES} ${QRC_FILES} ${UI_FILES}
        source/Analysis_GPU_Version.cpp
        source/IGA_CC_Helper.cpp
        source/IGA_CC_Helper.h
        source/nanoflann.h
        source/Topopt_GPU_Version.cpp
        source/Macro_Helper.h
        source/CoverageAxis.h)

if(WIN32)
    add_custom_command(TARGET ${APP_NAME}
            POST_BUILD
            COMMAND ${DEPLOYQT_EXECUTABLE} "$<TARGET_FILE:${APP_NAME}>"
            --$<LOWER_CASE:$<CONFIG>> --qmldir ${CMAKE_SOURCE_DIR}
    )
endif()

set(GUROBI_LIB_DIR "$ENV{GUROBI_HOME}\\lib")
set(GUROBI_INCLUDE_DIR "$ENV{GUROBI_HOME}\\include")

message(STATUS "GUROBI_LIB_DIR: ${GUROBI_LIB_DIR}")
message(STATUS "GUROBI_INCLUDE_DIR: ${GUROBI_INCLUDE_DIR}")

link_directories(${GUROBI_LIB_DIR})
include_directories(${GUROBI_INCLUDE_DIR})

if(CMAKE_BUILD_TYPE STREQUAL "Debug")
    set(GUROBI_CPP_LIB "${GUROBI_LIB_DIR}\\gurobi_c++mdd2017.lib")
else()
    set(GUROBI_CPP_LIB "${GUROBI_LIB_DIR}\\gurobi_c++md2017.lib")
endif()

target_link_libraries(
        ${APP_NAME}
        PUBLIC SolverGPU
        PUBLIC CommonLib
        PUBLIC MeshKernel
        PUBLIC json
        PUBLIC Tools
        PUBLIC freeglut
        PUBLIC OpenMeshCore
        PUBLIC OpenMeshTools
        PUBLIC SkeletalMesh
        PUBLIC tinyXML
        PUBLIC BSpline3D
        PUBLIC MKL_header_lib
        PUBLIC DIR_header_lib
        PUBLIC BMesh
        PUBLIC VolumeAlgorithm
        PUBLIC SurfaceAlgorithm
        PUBLIC SkeletalPipeline
        PUBLIC magic_enum::magic_enum
        PUBLIC Qt5::Widgets
        PUBLIC Qt5::Core
        PUBLIC Qt5::Gui
        PUBLIC Qt5::OpenGL
        PUBLIC "${GUROBI_LIB_DIR}\\gurobi110.lib"
        #PUBLIC "${GUROBI_LIB_DIR}\\gurobi_c++mdd2017.lib"
        PUBLIC "${GUROBI_CPP_LIB}"
)

target_include_directories(
        ${APP_NAME}
        PUBLIC source/
        PUBLIC ${CMAKE_CURRENT_LIST_DIR}/
        PUBLIC ${CMAKE_CURRENT_BINARY_DIR}/

)

file(GLOB_RECURSE EXDLL ${PROJECT_SOURCE_DIR}/ReferenceDLL/*.dll)

add_custom_command(TARGET ${APP_NAME} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${EXDLL}
        $<TARGET_FILE_DIR:${APP_NAME}>
)