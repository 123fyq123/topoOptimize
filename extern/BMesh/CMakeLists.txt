set(LIBRARY_NAME BMesh)

file(GLOB SOURCE_FILES *.*)

add_library(${LIBRARY_NAME} STATIC ${SOURCE_FILES})

target_link_libraries(
        ${LIBRARY_NAME}
        PUBLIC CommonLib
        PUBLIC MeshKernel
        PUBLIC json
        PUBLIC Tools
        PUBLIC OpenMeshCore
        PUBLIC OpenMeshTools
        PUBLIC OpenGL::GL
        PUBLIC freeglut
        PUBLIC Qt5::Widgets
        PUBLIC Qt5::Core
        PUBLIC Qt5::Gui
        PUBLIC Qt5::OpenGL
)
