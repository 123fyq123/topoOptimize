set(LIBRARY_NAME SkeletalPipeline)

file(GLOB SOURCE_FILES *.cpp)
file(GLOB INCLUDE_FILES *.hpp)

add_library(${LIBRARY_NAME} STATIC ${INCLUDE_FILES} ${SOURCE_FILES})
include_directories(../nlopt/src/api/)
target_link_libraries(
        ${LIBRARY_NAME}
        PUBLIC MeshKernel
        PUBLIC OpenMeshCore
        PUBLIC nlopt
)