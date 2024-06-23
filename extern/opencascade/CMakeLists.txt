set(OpenCASCADE_DIR opencascade770-install/cmake)

find_package(OpenCASCADE COMPONENTS FoundationClasses REQUIRED)
find_package(OpenCASCADE COMPONENTS ModelingAlgorithms REQUIRED)
find_package(OpenCASCADE COMPONENTS ModelingData REQUIRED)

set(LIBRARY_NAME helper_opencascade)

add_library(${LIBRARY_NAME} INTERFACE)

target_include_directories(
        ${LIBRARY_NAME}
        INTERFACE ${OpenCASCADE_INCLUDE_DIR})
target_link_directories(
        ${LIBRARY_NAME}
        INTERFACE ${OpenCASCADE_LIBRARY_DIR})
target_link_libraries(
        ${LIBRARY_NAME}
        INTERFACE ${OpenCASCADE_LIBRARIES}
)