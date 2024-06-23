set(LIBRARY_NAME Tools)

file(GLOB SOURCE_FILES *.*)
file(GLOB SUB_SOURCE_FILES ${CMAKE_CURRENT_LIST_DIR}/CrossField/*)

add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES} ${SUB_SOURCE_FILES}
        VoronoiDiagram.h)
target_link_libraries(
        ${LIBRARY_NAME}
        PUBLIC Geom
        PUBLIC MeshKernel
)