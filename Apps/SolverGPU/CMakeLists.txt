add_library(GPU_header_lib INTERFACE)
target_compile_features(GPU_header_lib INTERFACE cxx_std_17)
target_compile_definitions(GPU_header_lib
		INTERFACE INPUT_DIR=${PROJECT_SOURCE_DIR}/input/
		INTERFACE OUTPUT_DIR=${PROJECT_SOURCE_DIR}/output/
)
file(GLOB INCLUDE_FILES *.h)
target_include_directories(GPU_header_lib
		INTERFACE ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
		${INCLUDE_FILES}
)

# CUDA and C++ compiler flags
set(cxx_flags
		$<$<CXX_COMPILER_ID:MSVC>:-D_SCL_SECURE_NO_WARNINGS /openmp /std:c++17>           #Add MSVC-specific compiler flags here
		#	$<$<CXX_COMPILER_ID:GNU>: -fopenmp -std=c++17> #Add GCC/Clang-specific compiler flags here
		$<$<CXX_COMPILER_ID:GNU>:-Wall -m64 -fopenmp -O3 -std=c++17 -Wno-unused-function> #Add GCC/Clang-specific compiler flags here
)
set(cuda_flags
		-Xcompiler=$<$<CXX_COMPILER_ID:GNU>:-Wall -O3>
		-Xcompiler=$<$<CXX_COMPILER_ID:MSVC>:${OpenMP_CXX_FLAGS}>
		#Disables warning
		#177-D "function XXX was declared but never referenced"
		-Xcudafe "--display_error_number --diag_suppress=177"
		${CUDA_ARCHS}
		-rdc=true
		-lineinfo
		--expt-extended-lambda
		-use_fast_math
		$<$<CXX_COMPILER_ID:GNU>:-O3>
		--expt-relaxed-constexpr
		-Xptxas -warn-spills -res-usage
		--ptxas-options=-v
		--threads=24
		#--split-compile=0
		#		-G
		#		-O3
)

add_library(developer_flags INTERFACE)
target_compile_options(developer_flags INTERFACE
		$<$<COMPILE_LANGUAGE:CXX>:${cxx_flags}>
		$<$<COMPILE_LANGUAGE:CUDA>:${cuda_flags}>
)

target_link_libraries(GPU_header_lib INTERFACE $<BUILD_INTERFACE:developer_flags>)
target_link_libraries(GPU_header_lib INTERFACE OpenMP::OpenMP_CXX)
target_link_libraries(GPU_header_lib INTERFACE MKL_header_lib)


set(TargetName SolverGPU)


#add_executable(${TargetName}
add_library(${TargetName}
		GP_Types.h GP_Mesh.h GP_Mesh_Static.h GP_Attribute.h GP_Handle.h
		Kernels/GP_CUDA_Attribute.cuh Util/GP_Cuda_Query.h Util/GP_Vector.h
		GP_BSplineSolid.h GP_Geometry_info.h
		GP_Spline.h GP_Spline_Static.h GP_Spline.cpp GP_IGA_kernel.cuh
		GP_Spline_Test.h Eigen_helper.h GP_Spline_Attribute.h GP_Reduce_Handle.h GP_Param.cpp GP_Spline_v0.h GP_BoundaryCPoint.h GP_BoundaryCPoint.cpp
		Cuda_Linear_Elastic_Solver.h
        Cuda_TopOpt_Pipeline.h
		Kernels/GP_Debug.cu)

set(SOURCE_LIST
		test.cu
		GP_Mesh.cpp
		Patcher/GP_Patcher.cu
		GP_Invoke.cu
		Cuda_Linear_Elastic_Solver.cu
		Cuda_Hyperelasticity_Solver.cu
		Cuda_TopOpt_Pipeline.cu
		CPU_Linear_Elastic_Solver.cpp
        CPU_Hyperelasticity_Solver.cu
		CPU_TimeIntegrator_Hyperelasticity_Solver.cpp
		GP_Spline_CPU.cpp
		GP_IGA_CPU_Kernel.cpp
		)

target_sources(${TargetName}
		PRIVATE
		${SOURCE_LIST}
		)

set_target_properties(${TargetName} PROPERTIES FOLDER "apps")

set_property(TARGET ${TargetName} PROPERTY CUDA_SEPARABLE_COMPILATION ON)

source_group(TREE ${CMAKE_CURRENT_LIST_DIR} PREFIX "${TargetName}" FILES ${SOURCE_LIST})

target_link_libraries( ${TargetName}

#		PRIVATE RXMesh_header_lib
#		PRIVATE RXMesh_lib
#		PRIVATE gtest_main
		PRIVATE GPU_header_lib
		PRIVATE MeshKernel
		PRIVATE Hyperelasticity
		PRIVATE tinyXML
		PUBLIC CommonLib
		PUBLIC TBB::tbb
		)
#gtest_discover_tests( ${TargetName} )