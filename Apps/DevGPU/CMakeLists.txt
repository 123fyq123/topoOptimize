set(TargetName DevGPU)

add_executable(${TargetName} test.cu)

source_group(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "${TargetName}" FILES ${SOURCE_LIST})

target_link_libraries( ${TargetName}

        #		PRIVATE RXMesh_header_lib
        #		PRIVATE RXMesh_lib
        #		PRIVATE gtest_main
        PUBLIC SolverGPU
        PRIVATE GPU_header_lib
        PRIVATE MeshKernel
        PRIVATE Hyperelasticity
        PRIVATE tinyXML
        PUBLIC TBB::tbb
        PUBLIC OpenMP::OpenMP_CXX
)

file(GLOB_RECURSE EXDLL ${PROJECT_SOURCE_DIR}/ReferenceDLL/*.dll)

add_custom_command(TARGET ${TargetName} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy_if_different
        ${EXDLL}
        $<TARGET_FILE_DIR:${TargetName}>
)
