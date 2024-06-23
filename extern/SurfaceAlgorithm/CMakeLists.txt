set(LIBRARY_NAME SurfaceAlgorithm)

file(GLOB SOURCE_FILES *.*)
file(GLOB SUB_SOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/OpenMeshAlgorithm/*.*)

add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES} ${SUB_SOURCE_FILES})

target_link_libraries(
        ${LIBRARY_NAME}
        PUBLIC MeshKernel
        PUBLIC Tools
        PUBLIC OpenMeshCore
        PUBLIC OpenMeshTools
)