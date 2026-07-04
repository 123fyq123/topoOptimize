#pragma once
#define EIGEN_USE_MKL_ALL

#include <assert.h>
#include <fstream>
#include <memory>
#include <thread>

#include <cuda_profiler_api.h>
#include "IGA_helper_function.h"

#include "GP_Attribute.h"
#include "GP_Debug.h"
#include "GP_Handle.h"
#include "GP_Launch_Box.h"
#include "GP_Spline_Attribute.h"
#include "GP_Types.h"
#include "Kernels/GP_For_Each.cuh"
#include "Util/GP_Import_Obj.h"
#include "Util/GP_Log.h"
#include "Util/GP_Report.h"
#include "Util/GP_Timer.h"
#include "GP_Param.h"
#include <filesystem>
#include "Poster/source/Timer.h"
#pragma once

#include <assert.h>
#include <algorithm>
#include <fstream>
#include <memory>

#include "GP_Spline.h"

#include "Eigen_helper.h"
#include "GP_IGA_kernel.cuh"
#include "GP_IGA_kernel_optimized.cuh"
#include "GP_Reduce_Handle.h"

#include <Hyperelasticity/TimeIntegrator.h>
#include <Hyperelasticity/collision_detect.h>
#include <Hyperelasticity/utils.h>
#include <chrono>

#include <Eigen/PardisoSupport>

#include "GP_BoundaryCPoint.h"

namespace gpmesh {

	class GPSplineSolidV0 : public GPSplineSolid
	{
	private:
		SplineAttributeContainer spline_attribute_container;

	public:
		GPSplineSolidV0(const GPSplineSolidV0&) = delete;

		// debug part
		explicit GPSplineSolidV0(const bool quite = false) : GPSplineSolid() {};

		template <typename T>
		std::shared_ptr<SplineAttributeBase<T>> add_spline_cpoint_attribute(
			const std::vector<std::vector<std::vector<T>>>& local_cp_attributes,
			const std::string& name,
			layoutT                                         layout = SoA)
		{
			if (local_cp_attributes.empty()) {
				GPMESH_ERROR("Input cpoint attribute is empty");
			}

			assert(local_cp_attributes.size() == m_num_geometry);
			assert(local_cp_attributes[0].size() == 64);

			uint32_t num_attributes = local_cp_attributes[0][0].size();

			auto res_ptr =
				spline_attribute_container.template add<SplineAttributeBase<T>>(
					name.c_str(),
					m_num_geometry,
					64,
					num_attributes,
					LOCATION_ALL,
					layout);

			for (int g = 0; g < m_num_geometry; ++g) {
				for (uint16_t cp_idx = 0; cp_idx < 64; ++cp_idx) {
					for (uint32_t a = 0; a < num_attributes; ++a) {
						(*res_ptr)(g, cp_idx, a) =
							local_cp_attributes[g][cp_idx][a];
					}
				}
			}

			res_ptr->move(gpmesh::HOST, gpmesh::DEVICE);

			return res_ptr;
		}

		template <typename T>
		std::shared_ptr<SplineAttributeBase<T>> add_spline_cpoint_attribute(
			const std::string& name,
			const uint32_t     num_attributes,
			locationT          location = LOCATION_ALL,
			layoutT            layout = SoA)
		{
			return spline_attribute_container.template add<SplineAttributeBase<T>>(
				name.c_str(), m_num_geometry, 64, num_attributes, location, layout);
		}

		// one newton step performance test
		// memory optimization
		void debug_IGA_he_cg_performance_test_v1(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_condition>& m_dirichletSides,
			std::vector<boundary_condition>& m_neumannSides)
		{
			init(patch, m_dirichletSides, m_neumannSides, false);

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;


			// vel, acc, disp should be global variable
			Eigen::MatrixXd velVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
			Eigen::MatrixXd accVector;
			Eigen::MatrixXd dispVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

			Eigen::MatrixXd pre_disVector = dispVector;
			Eigen::MatrixXd pre_velVector = velVector;

			Eigen::MatrixXd pre_accVector;

			stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);
			stiff.system().matrix().setZero();
			stiff.system().rhs().setZero();

			massAssembler mass(patch, density, stiff.system());

			mass.assemble();

			Eigen::MatrixXd m_solVector = dispVector;

			Eigen::MatrixXd displacement;
			stiff.constructCoefs(m_solVector, displacement);

			stiff.stiffNonlinearMatrix(displacement, true);

			stiff.rhsVector();
			stiff.system().dataPush_Dirichlet();

			//if not using mkl, try use SimplicialLDLT
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
				mass.matrix());
			accVector = acc_solver.solve(stiff.getRhs());
			pre_accVector = accVector;

			GPMESH_INFO("Acc Solve Success!");

			uint32_t newton_step = 0;
			uint32_t newton_limit = 1;

			SparseSystem cpu_system = stiff.system();

			double absTol = 1e-12;
			double relTol = 1e-9;

			bool bQuit1 = false;
			bool bQuit2 = false;

			double initUpdateNorm;
			double initResidualNorm;

			// gpu backend
			debug_mass_matrix_assembly<double>(density);

			GPMESH_WARN("newton {} step start", newton_step);
			// stiff.assemble(m_solVector);

			cpu_system.matrix() =
				alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
			cpu_system.matrix().makeCompressed();
			cpu_system.rhs() =
				mass.matrix() *
				(alpha1(Beta, tStep) * (dispVector - m_solVector) +
					alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
				stiff.getRhs();

			m_sparse_system.matrix() = cpu_system.matrix();
			m_sparse_system.rhs() = cpu_system.rhs();

			debug_system_matrix_assembly<double>();

			// debug_test_system_matrix_gpu(cpu_system);

			Eigen::MatrixXd gpu_result;
			float           gpu_time = 0.0;

			this->debug_frame_process_helper_v0(gpu_result, gpu_time);

			Eigen::MatrixXd cpu_result;
			float           cpu_time = 0.0;
			//this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

			Eigen::MatrixXd groundtruth_result;
			float           groundtruth_time = 0.0;
			this->debug_solving_assembly_rhs_helper_groudtruth(
				groundtruth_result, groundtruth_time, false);

			double diff_cpu_res = (groundtruth_result - cpu_result).norm();
			double diff_gpu_res = (groundtruth_result - gpu_result).norm();
			double groundtruth_norm_res = groundtruth_result.norm();

			GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
			GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
			GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

			GPMESH_INFO("cpu time - {} (ms)", cpu_time);
			GPMESH_INFO("gpu time - {} (ms)", gpu_time);
			GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
		}

		void debug_IGA_pcg_performance_test_v0(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_condition>& m_dirichletSides,
			std::vector<boundary_condition>& m_neumannSides)
		{
			init(patch, m_dirichletSides, m_neumannSides, false);

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			// vel, acc, disp should be global variable
			Eigen::MatrixXd velVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
			Eigen::MatrixXd accVector;
			Eigen::MatrixXd dispVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

			Eigen::MatrixXd pre_disVector = dispVector;
			Eigen::MatrixXd pre_velVector = velVector;

			Eigen::MatrixXd pre_accVector;

			stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);
			stiff.system().matrix().setZero();
			stiff.system().rhs().setZero();

			massAssembler mass(patch, density, stiff.system());

			mass.assemble();

			Eigen::MatrixXd m_solVector = dispVector;

			Eigen::MatrixXd displacement;
			stiff.constructCoefs(m_solVector, displacement);

			stiff.stiffNonlinearMatrix(displacement, true);

			stiff.rhsVector();
			stiff.system().dataPush_Dirichlet();

			//if not using mkl, try use SimplicialLDLT
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
				mass.matrix());
			accVector = acc_solver.solve(stiff.getRhs());
			pre_accVector = accVector;

			GPMESH_INFO("Acc Solve Success!");

			uint32_t newton_step = 0;
			uint32_t newton_limit = 1;

			SparseSystem cpu_system = stiff.system();

			double absTol = 1e-12;
			double relTol = 1e-9;

			bool bQuit1 = false;
			bool bQuit2 = false;

			double initUpdateNorm;
			double initResidualNorm;

			// gpu backend
			debug_mass_matrix_assembly<double>(density);

			GPMESH_WARN("newton {} step start", newton_step);
			// stiff.assemble(m_solVector);

			cpu_system.matrix() =
				alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
			cpu_system.matrix().makeCompressed();
			cpu_system.rhs() =
				mass.matrix() *
				(alpha1(Beta, tStep) * (dispVector - m_solVector) +
					alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
				stiff.getRhs();

			m_sparse_system.matrix() = cpu_system.matrix();
			m_sparse_system.rhs() = cpu_system.rhs();

			debug_system_matrix_assembly<double>();

			// debug_test_system_matrix_gpu(cpu_system);

			Eigen::MatrixXd gpu_result;
			float           gpu_time = 0.0;

			bool bUsePCG = using_pcg();

			if (bUsePCG == false)
			{
				this->debug_frame_process_helper_v0(gpu_result, gpu_time);
			}
			else
			{
				this->debug_frame_process_PCG_jacobiPrecondition(gpu_result, gpu_time);
			}

			Eigen::MatrixXd cpu_result;
			float           cpu_time = 0.0;
			//this->debug_solving_assembly_rhs_helper_cpu_v2(cpu_result, cpu_time);

			Eigen::MatrixXd pcg_cpu_result;
			float pcg_cpu_time = 0.0;
			//this->debug_solving_assembly_rhs_helper_pcg_cpu(pcg_cpu_result, pcg_cpu_time);

			Eigen::MatrixXd groundtruth_result;
			float           groundtruth_time = 0.0;
			this->debug_solving_assembly_rhs_helper_groudtruth(
				groundtruth_result, groundtruth_time, false);

			double diff_cpu_res = (groundtruth_result - cpu_result).norm();
			double diff_pcg_res = (groundtruth_result - pcg_cpu_result).norm();
			double diff_gpu_res = (groundtruth_result - gpu_result).norm();
			double groundtruth_norm_res = groundtruth_result.norm();

			GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
			GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);

			GPMESH_INFO("pcg cpu calculation difference - {}", diff_pcg_res);

			GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

			GPMESH_INFO("cpu time - {} (ms)", cpu_time);

			if (bUsePCG)
			{
				GPMESH_INFO("gpu pcg time - {} (ms)", gpu_time);
			}
			else
			{
				GPMESH_INFO("gpu cg time - {} (ms)", gpu_time);
			}
			GPMESH_INFO("cpu pcg time - {} (ms)", pcg_cpu_time);
			GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
		}

		// local matrix assembly, local residual assembly
		void debug_IGA_hyperE_assembly(
			const vector<CBSplineSolid>& patch,
			std::vector<boundary_condition>& m_dirichletSides,
			std::vector<boundary_condition>& m_neumannSides)
		{

			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;

			init(patch, m_dirichletSides, m_neumannSides, false);

			GPUTimer timer_gpu;

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			Eigen::MatrixXd velVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
			Eigen::MatrixXd accVector;
			Eigen::MatrixXd dispVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

			Eigen::MatrixXd pre_disVector = dispVector;
			Eigen::MatrixXd pre_velVector = velVector;

			// Todo: 初始化加速度计算 - timeAssemble.initialize()
			// 加速度计算 a = F / a

			int    numIterations = 0;
			double residualNorm = 0;
			double initResidualNorm = 0;
			double updateNorm = 0;
			double initUpdateNorm = 1;


			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);

			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);

			timer_gpu.start();

			IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			timer_gpu.stop();
			float time = timer_gpu.elapsed_millis();
			GPMESH_INFO("gpu time for assembly - {} (ms)", time);

			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			timer_gpu.start();
			IGA_he_assembly_kernel_optimized_v0<double>
				<< <num_blocks, threads_per_blocks >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			timer_gpu.stop();
			time = timer_gpu.elapsed_millis();
			GPMESH_INFO("gpu time for optimized assembly - {} (ms)", time);

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

			stiffAssembler  stiff(patch, m_dirichletSides, m_neumannSides);
			Eigen::MatrixXd u;

			// seperate timeAssemble
			auto            dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
			Eigen::MatrixXd displacement;
			stiff.constructCoefs(dispSol, displacement);

			stiff.system().matrix().setZero();
			stiff.system().rhs().setZero();
			std::vector<Eigen::MatrixXd> test_localMatsVec;
			std::vector<Eigen::MatrixXd> test_localRhsVec;
			stiff.test_stiffNonlinearMatrix(
				displacement, true, test_localMatsVec, test_localRhsVec);

			std::vector<int> localIndex(64, 0);

			std::iota(localIndex.begin(), localIndex.end(), 0);

			std::vector<Eigen::MatrixXd> localMatsVec_gpu;
			localMatsVec_gpu.resize(m_num_geometry);

			// check localMats ColMajor properties
			for (int g = 0; g < m_num_geometry; ++g) {
				Eigen::MatrixXd localMats;

				localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
					m_h_assembly_info[g].local_mats, 192, 192);

				Eigen::MatrixXd localMats_v2;
				localMats_v2.setZero(192, 192);

				for (int i = 0; i < 192; ++i) {
					for (int j = 0; j < 192; ++j) {
						localMats_v2(i, j) =
							m_h_assembly_info[g].local_mats[i + 192 * j];
						auto t1 = localMats_v2(i, j);
						auto t2 = localMats(i, j);
						assert(abs(t1 - t2) < 0.001);
					}
				}

				m_sparse_system.pushToMatrix(g, localIndex, localMats);
				localMatsVec_gpu[g] = localMats_v2;
			}

			for (int g = 0; g < m_num_geometry; ++g) {
				auto mat_cpu = test_localMatsVec[g];
				auto mat_gpu = localMatsVec_gpu[g];

				for (int i = 0; i < 192; ++i) {
					for (int j = 0; j < 192; ++j) {
						assert(mat_cpu(i, j) - mat_gpu(i, j) < 0.1);
					}
				}
			}

			auto Matrix1 = m_sparse_system.test_get_bigmatrix();
			auto Matrix2 = stiff.system().test_get_bigmatrix();

			for (int i = 0; i < Matrix1.size(); ++i) {
				for (int j = 0; j < Matrix1[i].size(); ++j) {
					assert((Matrix1[i][j] - Matrix2[i][j]) < 0.0001);
				}
			}
			GPMESH_INFO("CPU and GPU Big Matrix Test Success!");
		}

		void debug_IGA_boundary_point_test_v0(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points)
		{
			init(patch, dirichlet_points, neumann_points, false);

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;



			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);
			IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);


			//Rhs Assembly
			//assemble into m_sparse_system.rhs()
			//neumann point

			m_sparse_system.rhs().setZero();

			helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
			helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);

			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			GPMESH_INFO("Use {}(ms)", gpu_pcg_time);
			string output_filename = "gpu_pcg_result2.txt";
			helper_output_matrix<double>(output_filename, gpu_pcg_result);
		}

		void debug_IGA_boundary_point_test_v0(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			double Yangs,
			double Poissons,
			Eigen::MatrixXd& U)
		{
			Timer timer;
			init(patch, dirichlet_points, neumann_points, false);

			timer.out("init");
			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			double* d_Yangs, * d_Poissons;
			cudaMalloc((void**)&d_Yangs, sizeof(double));
			cudaMalloc((void**)&d_Poissons, sizeof(double));

			cudaMemcpy(d_Yangs, &Yangs, sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(d_Poissons, &Poissons, sizeof(double), cudaMemcpyHostToDevice);


			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				Assembly_Local_Matrix_With_Yangs<double,
				gaussianPointNum>);
			Assembly_Local_Matrix_With_Yangs<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					d_Yangs,
					d_Poissons);

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			timer.out("Assembly");
			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);


			//Rhs Assembly
			//assemble into m_sparse_system.rhs()
			//neumann point

			m_sparse_system.rhs().setZero();

			helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
			helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);
			timer.out("Assembly");

			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			timer.out("PCG");
			GPMESH_INFO("Use {}(ms)", gpu_pcg_time);
			//string output_filename = "gpu_pcg_result2.txt";
			//helper_output_matrix<double>(output_filename, gpu_pcg_result);
			U = gpu_pcg_result;
		}

		void debug_IGA_linear_solver(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			double Yangs,
			double Poissons,
			Eigen::MatrixXd& U)
		{
			Timer timer;

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			double* d_Yangs, * d_Poissons;
			cudaMalloc((void**)&d_Yangs, sizeof(double));
			cudaMalloc((void**)&d_Poissons, sizeof(double));

			cudaMemcpy(d_Yangs, &Yangs, sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(d_Poissons, &Poissons, sizeof(double), cudaMemcpyHostToDevice);


			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				Assembly_Local_Matrix_With_Yangs<double,
				gaussianPointNum>);
			Assembly_Local_Matrix_With_Yangs<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					d_Yangs,
					d_Poissons);

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			timer.out("Assembly");
			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);


			//Rhs Assembly
			//assemble into m_sparse_system.rhs()
			//neumann point

			m_sparse_system.rhs().setZero();

			helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
			helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);
			timer.out("Assembly");

			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			timer.out("PCG");
			GPMESH_INFO("Use {}(ms)", gpu_pcg_time);
			//string output_filename = "gpu_pcg_result2.txt";
			//helper_output_matrix<double>(output_filename, gpu_pcg_result);
			U = gpu_pcg_result;
		}

		void TopOpt_GPU_Init(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points)
		{
			init(patch, dirichlet_points, neumann_points, false);
		}

		void Get_Local_Matrix(std::vector<Eigen::MatrixXd>& local_matrix, double Yangs, double Poissons) {
			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			double* d_Yangs, * d_Poissons;
			cudaMalloc((void**)&d_Yangs, sizeof(double));
			cudaMalloc((void**)&d_Poissons, sizeof(double));

			cudaMemcpy(d_Yangs, &Yangs, sizeof(double), cudaMemcpyHostToDevice);
			cudaMemcpy(d_Poissons, &Poissons, sizeof(double), cudaMemcpyHostToDevice);

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				Assembly_Local_Matrix_With_Yangs<double,
				gaussianPointNum>);
			Assembly_Local_Matrix_With_Yangs<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					d_Yangs,
					d_Poissons);


			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());



			double** temp_local_matrix;
			uint32_t patch_size = m_num_geometry;
			std::cout << m_num_geometry << std::endl;

			get_assembly_context().get_host_local_matrix(temp_local_matrix, m_num_geometry);


			std::cout << local_matrix.size() << std::endl;
			for (int patch_id = 0; patch_id < m_num_geometry; ++patch_id) {

				for (int i = 0; i < 192; ++i) {
					for (int j = 0; j < 192; ++j) {

						local_matrix[patch_id](i, j) = temp_local_matrix[patch_id][i + j * 192];
					}
				}

				/// todo: remove

			}
		}

		void Solver_With_Density(std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			std::vector<double>& topopt_density,
			Eigen::MatrixXd& U) {
			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			double* device_density;
			{
				int density_size = (int)topopt_density.size();
				cudaMalloc((void**)&device_density, density_size * sizeof(double));
				cudaMemcpy(device_density, topopt_density.data(), density_size * sizeof(double), cudaMemcpyHostToDevice);
			}

			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			IGA_clear_physics_displacement << <m_num_geometry, 192 >> > (get_physic_context());

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				TopOpt_Assembly_Local_Matrix<double,
				gaussianPointNum>);


			TopOpt_Assembly_Local_Matrix<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					device_density);


			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);



			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			U = gpu_pcg_result;
		}

		void Solver_With_Density_MMC(std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			std::vector<double>& topopt_density,
			Eigen::MatrixXd& U) {
			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			double* device_density;
			{
				int density_size = (int)topopt_density.size();
				cudaMalloc((void**)&device_density, density_size * sizeof(double));
				cudaMemcpy(device_density, topopt_density.data(), density_size * sizeof(double), cudaMemcpyHostToDevice);
			}
			//Stiffness Matrix Assembly
			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			IGA_clear_physics_displacement << <m_num_geometry, 192 >> > (get_physic_context());

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				TopOpt_Assembly_Local_Matrix<double,
				gaussianPointNum>);


			TopOpt_Assembly_Local_Matrix_MMC<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					device_density);


			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);



			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			U = gpu_pcg_result;
			cudaFree(device_density);
		}

		void Solver_With_MultiMaterial_Density(std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			const std::vector<double>& stiffness_coeff,
			Eigen::MatrixXd& U) {
			(void)neumann_points;
			(void)dirichlet_points;
			dim3 num_blocks_clear(m_num_geometry, 1, 1);
			dim3 num_blocks(m_num_geometry, 64, 1);
			dim3 threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;

			double* device_scale = nullptr;
			int scale_size = static_cast<int>(stiffness_coeff.size());
			cudaMalloc((void**)&device_scale, scale_size * sizeof(double));
			cudaMemcpy(device_scale, stiffness_coeff.data(), scale_size * sizeof(double), cudaMemcpyHostToDevice);

			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (m_assembly_context);
			IGA_clear_physics_displacement << <m_num_geometry, 192 >> > (get_physic_context());

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				TopOpt_Assembly_Local_Matrix_CellScale<double,
				gaussianPointNum>);

			TopOpt_Assembly_Local_Matrix_CellScale<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					device_scale);

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;
			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);

			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			U = gpu_pcg_result;
			cudaFree(device_scale);
		}

		void Solver_With_CellScale(std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			const std::vector<double>& cell_scale,
			Eigen::MatrixXd& U) {
			Solver_With_MultiMaterial_Density(neumann_points, dirichlet_points, cell_scale, U);
		}


		void debug_IGA_test_optimized_assembly_v0(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points)
		{
			init(patch, dirichlet_points, neumann_points, false);

			double timeSpan = 1;
			double timeStep = 0.1;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;



			//Stiffness Matrix Assembly
	//        IGA_clear_assembly_context<<<num_blocks_clear, threads_per_blocks>>>(
	//            m_assembly_context);
	//
	//        size_t smem_bytes_dyn_hyper = calc_shared_memory(
	//            (void*)
	//                IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
	//                                                                   gaussianPointNum>);
	//        IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
	//            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
	//                get_spline_context(),
	//                get_gaussian_integral_context(),
	//                get_physic_context(),
	//                get_assembly_context());
	//
	//        CUDA_ERROR(cudaDeviceSynchronize());
	//        CUDA_ERROR(cudaGetLastError());

			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);
			IGA_he_assembly_kernel_optimized_v0<double>
				<< <num_blocks, threads_per_blocks >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);


			//Rhs Assembly
			//assemble into m_sparse_system.rhs()
			//neumann point

			m_sparse_system.rhs().setZero();

			helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
			helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);

			Eigen::MatrixXd gpu_pcg_result;
			float gpu_pcg_time = 0.0;
			debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result, gpu_pcg_time);
			GPMESH_INFO("Use {}(ms)", gpu_pcg_time);


			//        IGA_clear_assembly_context<<<num_blocks_clear, threads_per_blocks>>>(
			//            m_assembly_context);
			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());
			//
			//        IGA_add_dirichlet_point_info_to_local_mats<<<DIVIDE_UP(Dirichlet_Point_Num, 256), 256>>>
			//            (get_spline_context(),
			//             get_assembly_context(),
			//             get_dirichelt_point_context(),
			//             Dirichlet_Point_Num);
			//
			//        Eigen::MatrixXd gpu_pcg_result2;
			//        debug_frame_process_PCG_jacobiPrecondition(gpu_pcg_result2, gpu_pcg_time);
			//        GPMESH_INFO("Use {}(ms)", gpu_pcg_time);

			GPMESH_INFO("Norm : {}", gpu_pcg_result.norm());

		}

		void helper_set_neumann_point_in_sparse_system_rhs(std::vector<boundary_point>& neumann_points)
		{
			for (auto bd_point : neumann_points)
			{
				assert(bd_point.m_condition_type == 1);

				uint32_t global_id = bd_point.m_cpoint_global_id;
				m_sparse_system.rhs()(global_id + m_num_cpoints * 0) = bd_point.m_val.getX();
				m_sparse_system.rhs()(global_id + m_num_cpoints * 1) = bd_point.m_val.getY();
				m_sparse_system.rhs()(global_id + m_num_cpoints * 2) = bd_point.m_val.getZ();
			}
		}

		void helper_set_dirichlet_point_in_sparse_system_rhs(std::vector<boundary_point>& dirichlet_points)
		{
			const double       MAXNUM_DIRICHLET = 1e10;
			for (auto bd_point : dirichlet_points)
			{
				assert(bd_point.m_condition_type == 0);

				uint32_t global_id = bd_point.m_cpoint_global_id;
				m_sparse_system.rhs()(global_id + m_num_cpoints * 0) = bd_point.m_val.getX() * MAXNUM_DIRICHLET;
				m_sparse_system.rhs()(global_id + m_num_cpoints * 1) = bd_point.m_val.getY() * MAXNUM_DIRICHLET;
				m_sparse_system.rhs()(global_id + m_num_cpoints * 2) = bd_point.m_val.getZ() * MAXNUM_DIRICHLET;
			}
		}

		//CG solver: solving one newton step
		void debug_frame_process_helper_v0(Eigen::MatrixXd& res,
			float& time,
			bool             bSilence = true)
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_num_geometry;

			Eigen::MatrixXd updateVector_cps;
			updateVector_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> updateVector_tripleVec;
			helper_convert_cps_to_tripleVec<double>(updateVector_cps,
				updateVector_tripleVec);

			auto updateVector_gpu_global_ptr =
				this->add_spline_cpoint_attribute_global<double>(
					updateVector_tripleVec, "updateVector_global");

			// X_global is the initial solution, influence the efficiency of solver
			auto X_global = updateVector_gpu_global_ptr;
			X_global->reset_global_part(0.0, LOCATION_ALL);
			// GPMESH_INFO("X Build success");

			auto S_global = this->add_spline_cpoint_attribute_global<double>(
				"S_global", 3, LOCATION_ALL, gpmesh::SoA);
			S_global->reset_local_part(0.0, LOCATION_ALL);
			// GPMESH_INFO("S local part reset success");
			S_global->reset_global_part(0.0, LOCATION_ALL);

			// GPMESH_INFO("S Build success");

			auto R_global = this->add_spline_cpoint_attribute_global<double>(
				"R_global", 3, LOCATION_ALL, gpmesh::SoA);
			R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
			R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
			auto P_global = this->add_spline_cpoint_attribute_global<double>(
				"P_global", 3, LOCATION_ALL, gpmesh::SoA);
			P_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
			P_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);

			// GPMESH_INFO("R, P Build success");

			//            ::getchar();
			//            ::getchar();
			// init_B
			Eigen::MatrixXd B_cps;
			B_cps.noalias() = m_sparse_system.rhs();
			//
			//            ::getchar();
			//            //transfer m_sparse_system.rhs() to gpu
			std::vector<std::vector<std::vector<double>>> B_tripleVec;
			helper_convert_cps_to_tripleVec<double>(B_cps, B_tripleVec);

			// CUDA_ERROR(cudaDeviceSynchronize());

			// load information to local part
			auto B_global = this->add_spline_cpoint_attribute_global<double>(
				B_tripleVec, "B_global");

			B_global->load_global_part(B_cps);

			// GPMESH_INFO("init B success");

			GPUTimer timer_gpu;

			timer_gpu.start();


			// S_global = A * X_global
			IGA_matvec_A_multiply_v1<double, 192> << <geometryNum, 192 >> > (
				get_spline_context(), get_assembly_context(), *X_global, *S_global);

			S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
			S_global->move_local_to_global_gpu(m_attr_mapping_context);

			// GPMESH_INFO("S calculate success");

			// init PR
			size_t cps_length = m_num_cpoints * 3;
			size_t init_PR_threads_num = 256;
			size_t init_PR_blocks_num = DIVIDE_UP(cps_length, 256);

			// R_global = B_global - S_global
			// P_global = R_global
			init_PR<double> << <init_PR_blocks_num, init_PR_threads_num >> > (
				B_global->m_d_global_attr,
				S_global->m_d_global_attr,
				R_global->m_d_global_attr,
				P_global->m_d_global_attr,
				cps_length);

			R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
			P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
			// GPMESH_INFO("init P R success");
			// CUDA_ERROR(cudaDeviceSynchronize());
			// CUDA_ERROR(cudaGetLastError());
			// 新方法计算Delta_new

			ReduceHandle<double> reduce_handle(*R_global);

			double Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

			double Delta_0_gpu = Delta_new_gpu;

			GPMESH_INFO("Delta_0 = {}", Delta_0_gpu);

			// move P_global global part to local part
			P_global->move_global_to_local_gpu(m_attr_mapping_context);
			S_global->reset_local_part(0.0, gpmesh::DEVICE);

			// CUDA_ERROR(cudaDeviceSynchronize());
			// CUDA_ERROR(cudaGetLastError());

			double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

			uint32_t iterNum = 0;
			uint32_t maxIter = 100000;

			while (iterNum <= maxIter) {


				S_global->reset_local_part(0.0, gpmesh::DEVICE);
				S_global->reset_global_part(0.0, gpmesh::DEVICE);

				IGA_matvec_A_multiply_v1<double, 192>
					<< <geometryNum, 192 >> > (get_spline_context(),
						get_assembly_context(),
						*P_global,
						*S_global);

				S_global->move_local_to_global_gpu(m_attr_mapping_context);
				S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

				Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

				Alpha_gpu = Delta_new_gpu / Alpha_gpu;

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*X_global,
						*P_global,
						Alpha_gpu,
						1.0,
						X_global->m_num_cpoints * 3);

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*R_global,
						*S_global,
						-Alpha_gpu,
						1.0,
						X_global->m_num_cpoints * 3);

				Delta_old_gpu = Delta_new_gpu;

				Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

				CUDA_ERROR(cudaStreamSynchronize(0));

				// 1e-16
				// 1e-13
				if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
					break;
				}

				Beta_gpu = Delta_new_gpu / Delta_old_gpu;

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*P_global,
						*R_global,
						1.0,
						Beta_gpu,
						X_global->m_num_cpoints * 3);

				P_global->move_global_to_local_gpu(m_attr_mapping_context);

				CUDA_ERROR(cudaStreamSynchronize(0));
				iterNum++;
			}

			timer_gpu.stop();
			time = timer_gpu.elapsed_millis();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			GPMESH_INFO("GPU takes {} iterations", iterNum);

			X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
			X_global->move_global_to_local_gpu(m_attr_mapping_context);

			for (int d = 0; d < m_num_cpoints * 3; ++d) {
				updateVector_cps(d, 0) = X_global->m_h_global_attr[d];
			}

			auto gpu_result = updateVector_cps;

			res = gpu_result;

			// in cpu, updateVector = updateVector_cps;
			//  m_solVector += updateVector;

			// in gpu, use X_global to update
			// m_physic_context.m_d_cpoints_displacement input X_global,
			// m_physic_context, m_attr_mapping_context

			debug_update_cpoint_displacement_in_physic_context(*X_global);

			// stiff matrix assembly verify should do in outside
		}

		void debug_hyperelasticity_init(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			double _timeStep
		)
		{
			assert(bInitPatch && bInitDirichlet);

			if (!bInitPatch)
			{
				throw std::runtime_error("Not Init Geometry");
			}

			if (!bInitDirichlet)
			{
				throw std::runtime_error("Not init Dirichlet");
			}


			double timeStep = _timeStep;
			double density = 3320;

			// gpu backend
			debug_mass_matrix_assembly_dirichlet_point<double>(density);

			//debug_system_matrix_assembly<double>();
			//debug_system_matrix_assembly_dirichlet_point<double>(timeStep);
		}

		void debug_hyperelasticity_process(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			double _timeStep
		)
		{
			assert(bInitPatch && bInitDirichlet);

			if (!bInitPatch)
			{
				throw std::runtime_error("Not Init Geometry");
			}

			if (!bInitDirichlet)
			{
				throw std::runtime_error("Not init Dirichlet");
			}

			double timeSpan = 1;
			double timeStep = _timeStep;
			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = _timeStep;


			// vel, acc, disp should be global variable
			Eigen::MatrixXd velVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
			Eigen::MatrixXd accVector;
			Eigen::MatrixXd dispVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

			Eigen::MatrixXd pre_disVector = dispVector;
			Eigen::MatrixXd pre_velVector = velVector;

			Eigen::MatrixXd pre_accVector;

			stiffAssembler stiff(patch, dirichlet_points, neumann_points);
			stiff.system().matrix().setZero();
			stiff.system().rhs().setZero();

			massAssembler mass(patch, density, stiff.system());

			mass.assemble();

			Eigen::MatrixXd m_solVector = dispVector;

			Eigen::MatrixXd displacement;
			stiff.constructCoefs(m_solVector, displacement);

			stiff.stiffNonlinearMatrix(displacement, true);

			stiff.rhsVector();
			stiff.system().dataPush_Dirichlet();

			//if not using mkl, try use SimplicialLDLT
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
				mass.matrix());
			accVector = acc_solver.solve(stiff.getRhs());
			pre_accVector = accVector;

			GPMESH_INFO("Acc Solve Success!");

			uint32_t newton_step = 0;
			uint32_t newton_limit = 1;

			SparseSystem cpu_system = stiff.system();

			double absTol = 1e-12;
			double relTol = 1e-9;

			bool bQuit1 = false;
			bool bQuit2 = false;

			double initUpdateNorm;
			double initResidualNorm;

			// gpu backend
			debug_mass_matrix_assembly_dirichlet_point<double>(density);

			GPMESH_WARN("newton {} step start", newton_step);
			// stiff.assemble(m_solVector);

			cpu_system.matrix() =
				alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
			cpu_system.matrix().makeCompressed();
			cpu_system.rhs() =
				mass.matrix() *
				(alpha1(Beta, tStep) * (dispVector - m_solVector) +
					alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
				stiff.getRhs();

			m_sparse_system.matrix() = cpu_system.matrix();
			m_sparse_system.rhs() = cpu_system.rhs();

			//debug_system_matrix_assembly<double>();
			debug_system_matrix_assembly_dirichlet_point<double>(timeStep);
			debug_test_system_matrix_gpu(cpu_system);

			Eigen::MatrixXd gpu_result;
			float           gpu_time = 0.0;

			this->debug_frame_process_PCG_jacobiPrecondition(gpu_result, gpu_time);

			Eigen::MatrixXd cpu_result;
			float           cpu_time = 0.0;
			//this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

			Eigen::MatrixXd groundtruth_result;
			float           groundtruth_time = 0.0;
			this->debug_solving_assembly_rhs_helper_groudtruth(
				groundtruth_result, groundtruth_time, false);

			double diff_cpu_res = (groundtruth_result - cpu_result).norm();
			double diff_gpu_res = (groundtruth_result - gpu_result).norm();
			double groundtruth_norm_res = groundtruth_result.norm();

			GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
			GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
			GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

			GPMESH_INFO("cpu time - {} (ms)", cpu_time);
			GPMESH_INFO("gpu time - {} (ms)", gpu_time);
			GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
		}

		void debug_IGA_hyperelasticity_frame_process_v0(
			vector<CBSplineSolid>& patch,
			std::vector<boundary_point>& neumann_points,
			std::vector<boundary_point>& dirichlet_points,
			double _timeStep
		)
		{
			init(patch, dirichlet_points, neumann_points, false);


			double density = 3320;

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = _timeStep;

			// vel, acc, disp should be global variable
			Eigen::MatrixXd velVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
			Eigen::MatrixXd accVector;
			Eigen::MatrixXd dispVector =
				Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

			Eigen::MatrixXd pre_disVector = dispVector;
			Eigen::MatrixXd pre_velVector = velVector;

			Eigen::MatrixXd pre_accVector;

			stiffAssembler stiff(patch, dirichlet_points, neumann_points);
			stiff.system().matrix().setZero();
			stiff.system().rhs().setZero();

			massAssembler mass(patch, density, stiff.system());

			mass.assemble();

			Eigen::MatrixXd m_solVector = dispVector;

			Eigen::MatrixXd displacement;
			stiff.constructCoefs(m_solVector, displacement);

			stiff.stiffNonlinearMatrix(displacement, true);

			stiff.rhsVector_NeumannPoint();
			stiff.system().dataPush_Dirichlet();

			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
				mass.matrix());
			accVector = acc_solver.solve(stiff.getRhs());

			pre_accVector = accVector;

			GPMESH_INFO("acc norm = {}", accVector.norm());

			uint32_t newton_step = 0;
			uint32_t newton_limit = 10;

			double absTol = 1e-12;
			double relTol = 1e-9;

			// gpu backend
			debug_mass_matrix_assembly_dirichlet_point(density);

			SparseSystem cpu_system = stiff.system();

			while (newton_step < newton_limit) {
				GPMESH_INFO("newton {} step start", newton_step);

				stiff.system().matrix().setZero();
				stiff.system().rhs().setZero();
				displacement.setZero();
				stiff.constructCoefs(m_solVector, displacement);
				stiff.stiffNonlinearMatrix(displacement, true);
				stiff.rhsVector_NeumannPoint();
				stiff.system().dataPush_Dirichlet();

				cpu_system.matrix() =
					alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();

				cpu_system.rhs() =
					mass.matrix() *
					(alpha1(Beta, tStep) * (dispVector - m_solVector) +
						alpha2(Beta, tStep) * velVector +
						alpha3(Beta) * accVector) +
					stiff.getRhs();

				debug_system_matrix_assembly_dirichlet_point<double>(tStep);

				m_sparse_system.matrix() = cpu_system.matrix();
				m_sparse_system.rhs() = cpu_system.rhs();

				Eigen::MatrixXd gpu_result;

				float gpu_time;

				debug_frame_process_PCG_jacobiPrecondition(gpu_result, gpu_time, false);
				CUDA_ERROR(cudaDeviceSynchronize());

				Eigen::MatrixXd groundtruth_result;

				float cpu_time;
				debug_frame_process_Pardiso(groundtruth_result, cpu_time);
				m_solVector += groundtruth_result;

				GPMESH_INFO(
					"updateVector Difference between groundtruth and gpu = {}",
					(groundtruth_result - gpu_result).norm());

				GPMESH_INFO("gpu update vector norm is {}", gpu_result.norm());

				Eigen::MatrixXd gpu_disp = debug_get_gpu_displacement();

				GPMESH_INFO("solVector difference between cpu and gpu - {}",
					(m_solVector - gpu_disp).norm());
				GPMESH_INFO("m_solVector norm is {}", m_solVector.norm());
				GPMESH_INFO("gpu_disp norm is {}", gpu_disp.norm());

				newton_step++;
			}
		};

		void debug_frame_process_Pardiso(Eigen::MatrixXd& cpu_result, float& cpu_time)
		{
			CPUTimer cpu_timer;


			auto& matrix = m_sparse_system.matrix();
			auto& rhs = m_sparse_system.rhs();

			cpu_timer.start();
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> pardiso_solver;
			pardiso_solver.compute(matrix);

			cpu_result = pardiso_solver.solve(rhs);

			cpu_timer.stop();

			cpu_time = cpu_timer.elapsed_millis();
		}

		void debug_solving_assembly_rhs_helper_cpu(Eigen::MatrixXd& res,
			float& time,
			bool             bSilence = true)
		{

			Eigen::SparseMatrix<double> A = m_sparse_system.matrix();

			Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);

			Eigen::VectorXd b = m_sparse_system.rhs();

			GPUTimer timer_cpu;
			timer_cpu.start();

			Eigen::VectorXd r = b - A * x;

			Eigen::VectorXd p = r;

			double delta_new = r.squaredNorm();
			double delta_0(delta_new);

			if (!bSilence)
				GPMESH_TRACE("delta_0 in CPU is {}", delta_0);

			Eigen::VectorXd Ap = A * p;

			if (!bSilence)
				GPMESH_TRACE("delta_new in CPU iteration - {} is {}", 1, delta_new);
			if (!bSilence)
				GPMESH_TRACE("p.dot(Ap) in CPU iteration-{} is {}", 1, p.dot(Ap));

			double alpha = delta_new / p.dot(Ap);
			if (!bSilence)
				GPMESH_TRACE("alpha in CPU iteration-{} is {}", 1, alpha);

			x += alpha * p;
			r -= alpha * Ap;

			double delta_old = delta_new;
			delta_new = r.squaredNorm();

			if (!bSilence)
				GPMESH_TRACE("delta_new in CPU iteration{} is {}", 1, delta_new);

			double beta = delta_new / delta_old;
			if (!bSilence)
				GPMESH_TRACE("beta in CPU iteration-{} is {}", 1, beta);

			p = r + beta * p;

			uint32_t iterNum = 1;

			uint32_t iterMax = 3;
			while (iterNum < iterMax) {
				// CPU Iter
				Ap = A * p;

				alpha = delta_new / p.dot(Ap);
				if (!bSilence)
					GPMESH_TRACE("CPU Iteration {} : ", iterNum);
				if (!bSilence)
					GPMESH_TRACE(
						"delta_old={}, p.dot(Ap)={}", delta_new, p.dot(Ap));

				x += alpha * p;
				r -= alpha * Ap;

				delta_old = delta_new;
				delta_new = r.squaredNorm();

				beta = delta_new / delta_old;

				p = r + beta * p;

				if (!bSilence)
					GPMESH_TRACE("delta_new={}", delta_new);

				if (delta_new < 1e-15 * delta_0) {
					break;
				}

				iterNum++;
			}
			timer_cpu.stop();

			time = timer_cpu.elapsed_millis();
			GPMESH_INFO("cpu solving take {} iterations", iterNum);
			auto cpu_result = x;

			res = cpu_result;
		}

		template <typename T>
		std::shared_ptr<SplineAttributeGlobal<T>>
			add_spline_cpoint_attribute_global(
				const std::vector<std::vector<std::vector<T>>>& local_cp_attributes,
				const std::string& name,
				layoutT                                         layout = SoA)
		{
			if (local_cp_attributes.empty()) {
				GPMESH_ERROR("Input cpoint attribute is empty");
				return nullptr;
			}

			assert(local_cp_attributes.size() == m_num_geometry);
			assert(local_cp_attributes[0].size() == 64);

			uint32_t num_attributes = local_cp_attributes[0][0].size();

			auto attr_global_ptr =
				spline_attribute_container.add_global_attr<T>(name.c_str(),
					m_num_geometry,
					64,
					num_attributes,
					m_num_cpoints,
					LOCATION_ALL,
					layout);

			for (int g = 0; g < m_num_geometry; ++g) {
				for (uint16_t cp_idx = 0; cp_idx < 64; ++cp_idx) {
					for (uint32_t a = 0; a < num_attributes; ++a) {
						(*attr_global_ptr)(g, cp_idx, a) =
							local_cp_attributes[g][cp_idx][a];
					}
				}
			}

			attr_global_ptr->move_local_part(gpmesh::HOST, gpmesh::DEVICE);

			return attr_global_ptr;
		}

		template <typename T>
		std::shared_ptr<SplineAttributeGlobal<T>>
			add_spline_cpoint_attribute_global(const std::string& name,
				const uint32_t     num_attributes,
				locationT location = LOCATION_ALL,
				layoutT   layout = SoA)
		{
			auto attr_global_ptr =
				spline_attribute_container.add_global_attr<T>(name.c_str(),
					m_num_geometry,
					64,
					num_attributes,
					m_num_cpoints,
					LOCATION_ALL,
					layout);

			return attr_global_ptr;
		}

		void helper_performance_test_for_matvec_A_multiply()
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_geometry_size;

			Eigen::MatrixXd Disp_cps;
			Disp_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
			helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

			auto Disp_gpu_ptr =
				this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp0");

			auto X = Disp_gpu_ptr;

			auto S = this->add_spline_cpoint_attribute<double>(
				"S0", 3, LOCATION_ALL, gpmesh::SoA);
			S->reset(0.0, LOCATION_ALL);

			GPUTimer timer_gpu;

			timer_gpu.start();

			IGA_debug_matvec_A_Multiply<double, localCPointNum>
				<< <geometryNum, localCPointNum >> > (
					get_spline_context(), get_assembly_context(), *X, *S);

			timer_gpu.stop();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			GPMESH_INFO("IGA_debug_matvec_A_Multiply Run Success!");
			GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
		}

		void helper_performance_test_for_matvec_A_multiply_v0()
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_geometry_size;

			Eigen::MatrixXd Disp_cps;
			Disp_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
			helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

			auto Disp_gpu_ptr =
				this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp1");

			auto X = Disp_gpu_ptr;

			auto S = this->add_spline_cpoint_attribute<double>(
				"S1", 3, LOCATION_ALL, gpmesh::SoA);
			S->reset(0.0, LOCATION_ALL);

			GPUTimer timer_gpu;

			timer_gpu.start();

			IGA_matvec_A_multiply_v1<double, 192> << <geometryNum, 192 >> > (
				get_spline_context(), get_assembly_context(), *X, *S);

			timer_gpu.stop();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			GPMESH_INFO("IGA_matvec_A_multiply_v1 Run Success!");
			GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
		}

		void helper_performance_test_for_matvec_A_multiply_v1()
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_geometry_size;

			Eigen::MatrixXd Disp_cps;
			Disp_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
			helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

			auto Disp_gpu_ptr =
				this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp2");

			auto X = Disp_gpu_ptr;

			auto S = this->add_spline_cpoint_attribute<double>(
				"S2", 3, LOCATION_ALL, gpmesh::SoA);
			S->reset(0.0, LOCATION_ALL);

			GPUTimer timer_gpu;

			timer_gpu.start();

			IGA_matvec_A_multiply_v2<double, 192> << <geometryNum, 192 >> > (
				get_spline_context(), get_assembly_context(), *X, *S);

			timer_gpu.stop();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			GPMESH_INFO("IGA_matvec_A_multiply_v2 Run Success!");
			GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
		}

		void debug_solving_assembly_rhs_helper_gpu(Eigen::MatrixXd& res,
			float& time,
			bool             bSilence = true)
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_geometry_size;

			Eigen::MatrixXd Disp_cps;
			Disp_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
			helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

			auto Disp_gpu_ptr =
				this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp");

			auto X = Disp_gpu_ptr;

			auto S = this->add_spline_cpoint_attribute<double>(
				"S", 3, LOCATION_ALL, gpmesh::SoA);
			S->reset(0.0, LOCATION_ALL);

			GPUTimer timer_gpu;

			timer_gpu.start();

			//            IGA_debug_matvec_A_Multiply<double, localCPointNum><<<
			//            geometryNum,
			//            localCPointNum>>>(
			//                    get_spline_context(),
			//                    get_assembly_context(),
			//                    *X,
			//                    *S
			//            );

			IGA_matvec_A_multiply_v1<double, 192> << <geometryNum, 192 >> > (
				get_spline_context(), get_assembly_context(), *X, *S);

			S->move(gpmesh::DEVICE, gpmesh::HOST);

			Eigen::MatrixXd B_cps = m_sparse_system.rhs();

			Eigen::MatrixXd S_cps;
			S_cps.setZero(m_dims * m_num_cpoints, 1);

			// convert S(gpu) to S_cps(cpu);
			for (int g = 0; g < m_num_geometry; ++g) {
				for (uint16_t local_idx = 0; local_idx < localCPointNum;
					++local_idx) {
					uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
					for (uint32_t d = 0; d < 3; ++d) {
						S_cps(global_idx + d * m_sparse_system.getControlNums()) +=
							S->operator()(g, local_idx, d);
					}
				}
			}

			Eigen::MatrixXd R_cps = B_cps - S_cps;
			Eigen::MatrixXd P_cps = R_cps;

			std::vector<std::vector<std::vector<double>>> B_tripleVec;
			helper_convert_cps_to_tripleVec(B_cps, B_tripleVec);

			std::vector<std::vector<std::vector<double>>> P_tripleVec;
			helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

			double Delta_new = 0;

			Delta_new = R_cps.squaredNorm();

			double Delta_0(Delta_new);

			if (!bSilence)
				GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

			auto B = this->add_spline_cpoint_attribute(B_tripleVec, "B");
			auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

			S->reset(0.0, gpmesh::DEVICE);

			// GPU CG Step 2: First Iteration
			IGA_debug_matvec_A_Multiply<double, localCPointNum>
				<< <geometryNum, localCPointNum >> > (
					get_spline_context(), get_assembly_context(), *P, *S);

			S->move(gpmesh::DEVICE, gpmesh::HOST);

			convert_spline_attribute_to_eigen_matrix(*S, S_cps);

			double Alpha = 0.0;

			for (int i = 0; i < S_cps.rows(); ++i) {
				Alpha += S_cps(i, 0) * P_cps(i, 0);
			}

			Alpha = Delta_new / Alpha;

			if (!bSilence)
				GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

			Disp_cps += Alpha * P_cps;
			R_cps -= Alpha * S_cps;

			double Delta_old = Delta_new;

			Delta_new = R_cps.squaredNorm();

			if (!bSilence)
				GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);

			double Beta = Delta_new / Delta_old;

			if (!bSilence)
				GPMESH_TRACE("Beta in GPU iteration-{} is {}", 1, Beta);

			P_cps = R_cps + Beta * P_cps;
			hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
			P->move(gpmesh::HOST, gpmesh::DEVICE);

			uint32_t iterNum = 1;

			while (iterNum < 100000) {
				iterNum++;

				S->reset(0.0, gpmesh::DEVICE);
				//                IGA_debug_matvec_A_Multiply<double,
				//                localCPointNum><<<geometryNum, localCPointNum>>>
				//                        (
				//                                get_spline_context(),
				//                                get_assembly_context(),
				//                                *P,
				//                                *S
				//                        );

				IGA_matvec_A_multiply_v1<double, 192> << <geometryNum, 192 >> > (
					get_spline_context(), get_assembly_context(), *P, *S);

				S->move(gpmesh::DEVICE, gpmesh::HOST);

				convert_spline_attribute_to_eigen_matrix(*S, S_cps);

				Alpha = 0.0;

				for (int i = 0; i < S_cps.rows(); ++i) {
					Alpha += S_cps(i, 0) * P_cps(i, 0);
				}
				if (!bSilence)
					GPMESH_TRACE("GPU Iteration {} : ", iterNum);
				if (!bSilence)
					GPMESH_TRACE("Delta_old={}, P.dot(AP)={}", Delta_new, Alpha);
				Alpha = Delta_new / Alpha;

				Disp_cps += Alpha * P_cps;

				R_cps -= Alpha * S_cps;

				Delta_old = Delta_new;

				Delta_new = R_cps.squaredNorm();

				Beta = Delta_new / Delta_old;

				P_cps = R_cps + Beta * P_cps;


				if (!bSilence)
					GPMESH_TRACE("Delta_new={}", Delta_new);

				hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
				P->move(gpmesh::HOST, gpmesh::DEVICE);

				if (Delta_new < 1e-12 * Delta_0) {
					break;
				}

				CUDA_ERROR(cudaStreamSynchronize(0));
			}
			timer_gpu.stop();

			time = timer_gpu.elapsed_millis();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			auto gpu_result = Disp_cps;

			res = gpu_result;
		}

		//PCG solver: solving one newton step
		void debug_frame_process_PCG_jacobiPrecondition(Eigen::MatrixXd& res,
			float& time,
			bool             bSilence = true)
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;
			size_t             geometryNum = m_num_geometry;
			size_t m_num_threads_for_global = m_num_cpoints * 3;

			Eigen::MatrixXd updateVector_cps;
			updateVector_cps.setZero(m_num_cpoints * m_dims, 1);

			//[i][j][d]
			std::vector<std::vector<std::vector<double>>> updateVector_tripleVec;
			helper_convert_cps_to_tripleVec<double>(updateVector_cps,
				updateVector_tripleVec);

			auto updateVector_gpu_global_ptr =
				this->add_spline_cpoint_attribute_global<double>(
					updateVector_tripleVec, "updateVector_global");

			// X_global is the initial solution, influence the efficiency of solver
			auto X_global = updateVector_gpu_global_ptr;
			X_global->reset_global_part(0.0, LOCATION_ALL);
			// GPMESH_INFO("X Build success");

			auto S_global = this->add_spline_cpoint_attribute_global<double>(
				"S_global", 3, LOCATION_ALL, gpmesh::SoA);
			S_global->reset_local_part(0.0, LOCATION_ALL);
			// GPMESH_INFO("S local part reset success");
			S_global->reset_global_part(0.0, LOCATION_ALL);

			// GPMESH_INFO("S Build success");

			auto R_global = this->add_spline_cpoint_attribute_global<double>(
				"R_global", 3, LOCATION_ALL, gpmesh::SoA);
			R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
			R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
			auto P_global = this->add_spline_cpoint_attribute_global<double>(
				"P_global", 3, LOCATION_ALL, gpmesh::SoA);
			P_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
			P_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);

			// GPMESH_INFO("R, P Build success");

			//            ::getchar();
			//            ::getchar();
			// init_B
			Eigen::MatrixXd B_cps;
			B_cps.noalias() = m_sparse_system.rhs();
			//
			//            ::getchar();
			//            //transfer m_sparse_system.rhs() to gpu
			std::vector<std::vector<std::vector<double>>> B_tripleVec;
			helper_convert_cps_to_tripleVec<double>(B_cps, B_tripleVec);

			// CUDA_ERROR(cudaDeviceSynchronize());

			// load information to local part
			auto B_global = this->add_spline_cpoint_attribute_global<double>(
				B_tripleVec, "B_global");

			B_global->load_global_part(B_cps);

			// initialize Jacobi Precondition
			auto Preconditioner = this->add_spline_cpoint_attribute_global<double>(
				"Preconditioner", 3, LOCATION_ALL, gpmesh::SoA);

			Preconditioner->reset_global_part(0, gpmesh::LOCATION_ALL);
			Preconditioner->reset_local_part(0, gpmesh::LOCATION_ALL);

			build_preconditioner_from_system_matrix(Preconditioner);

			{
				//test for preconditioner
				Preconditioner->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

				//和 m_sparse_system.matrix()进行比较
				debug_test_jacobi_precondition(Preconditioner);
			}

			GPUTimer timer_gpu;

			timer_gpu.start();


			// S_global = A * X_global
			IGA_matvec_A_multiply_v1<double, 192> << <geometryNum, 192 >> > (
				get_spline_context(), get_assembly_context(), *X_global, *S_global);

			S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
			S_global->move_local_to_global_gpu(m_attr_mapping_context);

			// GPMESH_INFO("S calculate success");

			// init PR
			size_t cps_length = m_num_cpoints * 3;
			size_t init_PR_threads_num = 256;
			size_t init_PR_blocks_num = DIVIDE_UP(cps_length, 256);

			// R_global = B_global - S_global
			// P_global = R_global
			init_PR<double> << <init_PR_blocks_num, init_PR_threads_num >> > (
				B_global->m_d_global_attr,
				S_global->m_d_global_attr,
				R_global->m_d_global_attr,
				P_global->m_d_global_attr,
				cps_length);

			IGA_kernel_CWise_operation_dot<double>
				<< <DIVIDE_UP(m_num_threads_for_global, 256), 256 >> >
				(
					*P_global,
					*Preconditioner,
					m_num_threads_for_global
					);


			//Todo List:
			//1. line to line compare PCG Process,
			// referring to GPSplineSolidTest.debug_IGA_solving_process_test_v3

			//2. Folloing pcg in cpu, build pcg in gpu

			// 新方法计算Delta_new

			ReduceHandle<double> reduce_handle(*R_global);

			double Delta_new_gpu = reduce_handle.dot(*R_global, *P_global);

			double Delta_0_gpu = Delta_new_gpu;

			GPMESH_INFO("Delta_0 = {}", Delta_0_gpu);

			// move P_global global part to local part
			P_global->move_global_to_local_gpu(m_attr_mapping_context);
			S_global->reset_local_part(0.0, gpmesh::DEVICE);

			// CUDA_ERROR(cudaDeviceSynchronize());
			// CUDA_ERROR(cudaGetLastError());

			double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

			uint32_t iterNum = 0;
			uint32_t maxIter = 5000;

			while (iterNum <= maxIter) {


				S_global->reset_local_part(0.0, gpmesh::DEVICE);
				S_global->reset_global_part(0.0, gpmesh::DEVICE);

				IGA_matvec_A_multiply_v1<double, 192>
					<< <geometryNum, 192 >> > (get_spline_context(),
						get_assembly_context(),
						*P_global,
						*S_global);

				S_global->move_local_to_global_gpu(m_attr_mapping_context);
				//S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

				Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

				Alpha_gpu = Delta_new_gpu / Alpha_gpu;

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*X_global,
						*P_global,
						Alpha_gpu,
						1.0,
						X_global->m_num_cpoints * 3);

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*R_global,
						*S_global,
						-Alpha_gpu,
						1.0,
						X_global->m_num_cpoints * 3);

				Delta_old_gpu = Delta_new_gpu;

				IGA_kernel_CWise_operation_dot<double>
					<< <DIVIDE_UP(m_num_threads_for_global, 256), 256 >> >
					(*S_global,
						*Preconditioner,
						*R_global,
						m_num_threads_for_global
						);

				Delta_new_gpu = reduce_handle.dot(*S_global, *R_global);

				CUDA_ERROR(cudaStreamSynchronize(0));

				// 1e-16
				// 1e-13
				if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
					break;
				}

				Beta_gpu = Delta_new_gpu / Delta_old_gpu;

				axpy<double> << <reduce_handle.get_num_block(),
					reduce_handle.get_num_thread() >> > (
						*P_global,
						*S_global,
						1.0,
						Beta_gpu,
						X_global->m_num_cpoints * 3);

				P_global->move_global_to_local_gpu(m_attr_mapping_context);

				CUDA_ERROR(cudaStreamSynchronize(0));
				iterNum++;
			}

			timer_gpu.stop();
			time = timer_gpu.elapsed_millis();
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
			CUDA_ERROR(cudaProfilerStop());

			GPMESH_INFO("GPU takes {} iterations", iterNum);

			X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
			X_global->move_global_to_local_gpu(m_attr_mapping_context);

			for (int d = 0; d < m_num_cpoints * 3; ++d) {
				updateVector_cps(d, 0) = X_global->m_h_global_attr[d];
			}

			auto gpu_result = updateVector_cps;

			res = gpu_result;

			// in cpu, updateVector = updateVector_cps;
			//  m_solVector += updateVector;

			// in gpu, use X_global to update
			// m_physic_context.m_d_cpoints_displacement input X_global,
			// m_physic_context, m_attr_mapping_context

			debug_update_cpoint_displacement_in_physic_context(*X_global);

			// stiff matrix assembly verify should do in outside
		}

		template <typename T>
		void debug_push_gpu_local_mats_to_sparse_system()
		{
			m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);
			std::vector<int> localIndex(64, 0);
			std::iota(localIndex.begin(), localIndex.end(), 0);

			for (int g = 0; g < m_num_geometry; ++g) {
				Eigen::MatrixXd localMats;

				localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
					m_h_assembly_info[g].local_mats, 192, 192);
				m_sparse_system.pushToMatrix(g, localIndex, localMats);
			}

			m_sparse_system.dataPush();
		}

		template <typename T>
		void debug_mass_matrix_assembly(T density)
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;

			dim3 num_blocks(m_num_geometry, 64, 1);
			dim3 threads_per_blocks(16, 16, 1);

			const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
			const int GeometrySurfNum = 16;

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>);

			IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_mass_context(),
					density);

			IGA_add_dirichlet_to_mass_local_mats << <DirichNum, GeometrySurfNum >> > (
				get_spline_context(), get_mass_context(), get_dirichlet_context());
		}

		template <typename T>
		void debug_mass_matrix_assembly_dirichlet_point(T density)
		{
			constexpr uint32_t gaussianPointNum = 64;
			const int          localCPointNum = 64;
			const double       MAXNUM_DIRICHLET = 1e10;

			dim3 num_blocks(m_num_geometry, 64, 1);
			dim3 threads_per_blocks(16, 16, 1);

			const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
			const int GeometrySurfNum = 16;

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>);

			IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_mass_context(),
					density);

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_to_mass_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_mass_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);

		}


		template <typename T>
		void debugging_system_matrix_assembly()
		{

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;


			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);

			IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
			const int GeometrySurfNum = 16;
			IGA_add_dirichlet_to_local_mats << <DirichNum, GeometrySurfNum >> > (
				get_spline_context(),
				get_assembly_context(),
				get_dirichlet_context());

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			IGA_add_mass_matrix_to_stiff_matrix<double, 256>
				<< <m_num_geometry, 256 >> > (get_assembly_context(),
					get_mass_context(),
					alpha1(Beta, tStep));
		}

		template <typename T>
		void debug_system_matrix_assembly()
		{

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;


			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);

			IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
			const int GeometrySurfNum = 16;
			IGA_add_dirichlet_to_local_mats << <DirichNum, GeometrySurfNum >> > (
				get_spline_context(),
				get_assembly_context(),
				get_dirichlet_context());

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = 0.1;

			IGA_add_mass_matrix_to_stiff_matrix<double, 256>
				<< <m_num_geometry, 256 >> > (get_assembly_context(),
					get_mass_context(),
					alpha1(Beta, tStep));
		}

		template <typename T>
		void debug_system_matrix_assembly_dirichlet_point(double timeStep)
		{

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;


			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);

			IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = timeStep;

			IGA_add_mass_matrix_to_stiff_matrix<double, 256>
				<< <m_num_geometry, 256 >> > (get_assembly_context(),
					get_mass_context(),
					alpha1(Beta, tStep));
		}

		template <typename T>
		void debug_system_matrix_assembly_dirichlet_point_v2(double timeStep)
		{

			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;


			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
				gaussianPointNum>);

			IGA_hyperelasticity_assembly_debug_kernel_v3<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context(),
					get_lambda(),
					get_mu());

			//        IGA_he_assembly_kernel_optimized_v0<double>
			//            <<<num_blocks, threads_per_blocks>>>(
			//                get_spline_context(),
			//                get_gaussian_integral_context(),
			//                get_physic_context(),
			//                get_assembly_context());

			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());

			const uint32_t Dirichlet_Point_Num = m_dirichlet_point_context.m_num_dirichlet_point;

			IGA_add_dirichlet_point_info_to_local_mats << <DIVIDE_UP(Dirichlet_Point_Num, 256), 256 >> >
				(get_spline_context(),
					get_assembly_context(),
					get_dirichelt_point_context(),
					Dirichlet_Point_Num);

			double Beta = 0.25;
			double Gamma = 0.5;
			double tStep = timeStep;

			IGA_add_mass_matrix_to_stiff_matrix<double, 256>
				<< <m_num_geometry, 256 >> > (get_assembly_context(),
					get_mass_context(),
					alpha1(Beta, tStep));
		}

		template <typename T>
		void build_preconditioner_from_system_matrix(std::shared_ptr<SplineAttributeGlobal<T>> ptr_attribute_global)
		{
			IGA_kernel_build_preconditioner_from_system_matrix<double> << <m_num_geometry, 192 >> > (*ptr_attribute_global,
				get_assembly_context());
			ptr_attribute_global->move_local_to_global_gpu(m_attr_mapping_context);

			//get reciporal value of ptr_attribute_global;
			uint32_t num_threads = m_num_cpoints * 3;
			IGA_kernel_convert_preconditioner_reciprocal << <DIVIDE_UP(num_threads, 256), 256 >> > (*ptr_attribute_global, num_threads);

		}

		template <typename T>
		void debug_stiff_matrix_assembly()
		{
			dim3               num_blocks_clear(m_num_geometry, 1, 1);
			dim3               num_blocks(m_num_geometry, 64, 1);
			dim3               threads_per_blocks(16, 16, 1);
			constexpr uint32_t gaussianPointNum = 64;


			IGA_clear_assembly_context << <num_blocks_clear, threads_per_blocks >> > (
				m_assembly_context);

			this->m_quite = true;
			size_t smem_bytes_dyn_hyper = calc_shared_memory(
				(void*)
				IGA_hyperelasticity_assembly_debug_kernel_v2<double,
				gaussianPointNum>);

			IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
				<< <num_blocks, threads_per_blocks, smem_bytes_dyn_hyper >> > (
					get_spline_context(),
					get_gaussian_integral_context(),
					get_physic_context(),
					get_assembly_context());

			const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
			const int GeometrySurfNum = 16;
			IGA_add_dirichlet_to_local_mats << <DirichNum, GeometrySurfNum >> > (
				get_spline_context(),
				get_assembly_context(),
				get_dirichlet_context());
		}

		template <typename T>
		void debug_update_cpoint_displacement_in_physic_context(
			SplineAttributeGlobal<T>& X)
		{
			IGA_update_physics_context_by_spline_attr_global << <m_num_geometry,
				192 >> > (
					X, m_physic_context, m_attr_mapping_context);
			CUDA_ERROR(cudaDeviceSynchronize());
			CUDA_ERROR(cudaGetLastError());
		}

		void debug_solving_assembly_rhs_helper_groudtruth(Eigen::MatrixXd& res,
			float& time,
			bool bSilence = true)
		{
			Eigen::VectorXd             u_conjgrad;
			Eigen::SparseMatrix<double> A;
			A = m_sparse_system.matrix();
			Eigen::VectorXd B;
			B.noalias() = m_sparse_system.rhs();
			u_conjgrad.resize(m_dims * m_num_cpoints);
			u_conjgrad.setZero();

			CPUTimer timer_cpu;

			timer_cpu.start();
			//conjgrad<double>(A, B, u_conjgrad);
			timer_cpu.stop();

			GPMESH_INFO("conjgrad takes {}(ms)", timer_cpu.elapsed_millis());

			timer_cpu.start();

			Eigen::VectorXd u_cg;

			Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
			Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> pardisoSolve;

			//        // cg.setMaxIterations(50);
			//        cg.compute(m_sparse_system.matrix());
			//        u_cg = cg.solve(m_sparse_system.rhs());
			//        if (cg.info() == Eigen::NoConvergence ||
			//            cg.info() == Eigen::NumericalIssue) {
			//            if (!bSilence)
			//                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
			//            if (!bSilence)
			//                GPMESH_TRACE(
			//                    "Iteration {}, Error {}", cg.iterations(), cg.error());
			//        } else {
			//            if (!bSilence)
			//                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
			//            if (!bSilence)
			//                GPMESH_TRACE(
			//                    "Iteration {}, Error {}", cg.iterations(), cg.error());
			//        }

			// cg.setMaxIterations(50);
			pardisoSolve.compute(m_sparse_system.matrix());
			u_cg = pardisoSolve.solve(m_sparse_system.rhs());
			if (pardisoSolve.info() == Eigen::NoConvergence ||
				pardisoSolve.info() == Eigen::NumericalIssue) {
			}
			else {
				if (!bSilence)
					GPMESH_TRACE("Pardiso LDLT Solve Success!");
			}

			if (!bSilence)
				GPMESH_INFO("Difference between u_cg and u_conjgrad = {}",
					(u_cg - u_conjgrad).norm());
			// assert((u_cg - u_conjgrad).norm() < 1e-6);
			if (!bSilence)
				GPMESH_INFO("conjgrad Solver Test Success!");

			timer_cpu.stop();

			time = timer_cpu.elapsed_millis();

			res = u_cg;
		}

		void debug_solving_assembly_rhs_helper_cpu_v2(Eigen::MatrixXd& res,
			float& time,
			bool bSilence = true)
		{
			Eigen::VectorXd             u_conjgrad;
			Eigen::SparseMatrix<double> A;
			A = m_sparse_system.matrix();
			Eigen::VectorXd B;
			B.noalias() = m_sparse_system.rhs();
			u_conjgrad.resize(m_dims * m_num_cpoints);
			u_conjgrad.setZero();

			CPUTimer timer_cpu;

			timer_cpu.start();
			conjgrad<double>(A, B, u_conjgrad);
			timer_cpu.stop();

			GPMESH_INFO("conjgrad takes {}(ms)", timer_cpu.elapsed_millis());

			res = u_conjgrad;
			time = timer_cpu.elapsed_millis();
		}

		void debug_solving_assembly_rhs_helper_pcg_cpu(Eigen::MatrixXd& res,
			float& time,
			bool bSilence = true)
		{
			Eigen::VectorXd             u_conjgrad;
			Eigen::SparseMatrix<double> A;
			A = m_sparse_system.matrix();
			Eigen::VectorXd B;
			B.noalias() = m_sparse_system.rhs();
			u_conjgrad.resize(m_dims * m_num_cpoints);
			u_conjgrad.setZero();

			CPUTimer timer_cpu;

			timer_cpu.start();
			preconditioned_conjgrad<double>(A, B, u_conjgrad);
			timer_cpu.stop();

			GPMESH_INFO("pcg takes {}(ms)", timer_cpu.elapsed_millis());

			res = u_conjgrad;
			time = timer_cpu.elapsed_millis();
		}

		Eigen::MatrixXd debug_get_gpu_displacement()
		{
			m_physic_context.update_host_info(m_h_cpoints_displacement_array,
				m_num_geometry);

			Eigen::MatrixXd res_disp;
			res_disp.setZero(m_num_cpoints * m_dims, 1);

			for (int g = 0; g < m_num_geometry; ++g) {
				for (int i = 0; i < 64; ++i) {
					int global_id = m_h_geometry_ltog_cp[g][i];

					for (int d = 0; d < 3; ++d) {
						res_disp(global_id + d * m_num_cpoints, 0) =
							m_h_cpoints_displacement_array[g](i, d);
					}
				}
			}

			return res_disp;
		}

		// attention this function change m_sparse_system.matrix
		void debug_test_system_matrix_gpu(SparseSystem& sparse_system)
		{
			m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

			vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

			std::vector<int> localIndex(64, 0);
			std::iota(localIndex.begin(), localIndex.end(), 0);

			for (int g = 0; g < m_num_geometry; ++g) {
				Eigen::MatrixXd localMats;

				localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
					m_h_assembly_info[g].local_mats, 192, 192);
				m_sparse_system.pushToMatrix(g, localIndex, localMats);
			}
			m_sparse_system.matrix().setZero();
			// here don`t add dirichlet information
			m_sparse_system.dataPush();

			GPMESH_INFO("Input sparse system matrix with gpu matrix");
			GPMESH_INFO("difference = {}",
				(m_sparse_system.matrix() - sparse_system.matrix()).norm());
		}

		void debug_test_jacobi_precondition(shared_ptr<SplineAttributeGlobal<double>> preconditioner)
		{
			Eigen::VectorXd test;
			test.resize(m_num_cpoints * 3);

			for (int i = 0; i < test.size(); ++i)
			{
				double val = m_sparse_system.matrix().coeffRef(i, i);

				if (abs(val) < 1e-10)
				{
					val = 1;
				}
				else
				{
					val = 1.0 / val;
				}
				test(i) = val;
			}

			debug_test_vector_gpu(test, preconditioner);
		}

		void debug_test_vector_gpu(Eigen::VectorXd& x, shared_ptr<SplineAttributeGlobal<double>> X)
		{
			//X->move_local_to_global_gpu(m_attr_mapping_context);
			//X->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

			double test = 0;

			for (int i = 0; i < x.rows(); ++i)
			{
				double Xi = X->m_h_global_attr[i];
				double diff = abs(Xi - x(i));

				if (diff > 1e-10)
				{
					//GPMESH_WARN("data[{}] show difference = {}", i, diff);
				}

				test += diff * diff;
			}

			GPMESH_INFO("data all errors = {}", test);
		}

		// compare local_mats between cpu and gpu, not adding dirichlet boundary
		// information
		void debug_test_localMatsVecCPU(
			std::vector<Eigen::MatrixXd>& local_mats_cpu)
		{
			assert(local_mats_cpu.size() == m_num_geometry);

			m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

			std::vector<Eigen::MatrixXd> local_mats_gpu;
			local_mats_gpu.resize(m_num_geometry);

			m_sparse_system.matrix().setZero();

			// check localMats ColMajor properties
			for (int g = 0; g < m_num_geometry; ++g) {
				Eigen::MatrixXd localMats;

				localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
					m_h_assembly_info[g].local_mats, 192, 192);

				Eigen::MatrixXd localMats_v2;
				localMats_v2.setZero(192, 192);

				std::vector<int> localIndex(64, 0);

				std::iota(localIndex.begin(), localIndex.end(), 0);

				for (int i = 0; i < 192; ++i) {
					for (int j = 0; j < 192; ++j) {
						localMats_v2(i, j) =
							m_h_assembly_info[g].local_mats[i + 192 * j];
						auto t1 = localMats_v2(i, j);
						auto t2 = localMats(i, j);
						assert(abs(t1 - t2) < 0.001);
					}
				}

				m_sparse_system.pushToMatrix(g, localIndex, localMats);
				local_mats_gpu[g] = localMats_v2;
			}

			for (int g = 0; g < m_num_geometry; ++g) {
				auto mat_cpu = local_mats_cpu[g];
				auto mat_gpu = local_mats_gpu[g];

				for (int i = 0; i < 192; ++i) {
					for (int j = 0; j < 192; ++j) {
						double diff = mat_cpu(i, j) - mat_gpu(i, j);
						if (abs(diff) > 0.001) {
							GPMESH_ERROR(
								"local mats check error for G-{}({},{}) = {}",
								g,
								i,
								j,
								diff);
						}
					}
				}
			}
		}



		size_t calc_shared_memory(const void* kernel)
		{
			// check if total shared memory (static + dynamic) consumed by
			// k_base_query are less than the max shared per block

			// here especially set shared memory per block to 96kb
			int maxbytes = 98304;
			cudaFuncSetAttribute(
				kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);

			cudaFuncAttributes func_attr = cudaFuncAttributes();
			CUDA_ERROR(cudaFuncGetAttributes(&func_attr, kernel));


			size_t smem_bytes_static = func_attr.sharedSizeBytes;
			size_t num_reg_per_thread = static_cast<uint32_t>(func_attr.numRegs);
			int    device_id;
			CUDA_ERROR(cudaGetDevice(&device_id));
			cudaDeviceProp devProp;
			CUDA_ERROR(cudaGetDeviceProperties(&devProp, device_id));

			uint16_t max_num_cp = 0;

			for (int i = 0; i < m_num_geometry; ++i) {
				max_num_cp = std::max(max_num_cp, m_h_geometry_info[i].num_CPoints);
			}

			m_max_cp_per_geometry = static_cast<uint32_t>(max_num_cp);

			size_t local_mats_size = 192 * 192;

			size_t smem_bytes_dyn = max_num_cp * m_dims * sizeof(double);

			if (!this->m_quite) {
				GPMESH_TRACE(
					"GPMeshStatic::check_shared_memory() user function requires "
					"shared memory = {} (dynamic) + {} (static) = {} (bytes) and "
					"{} registers per thread",
					smem_bytes_dyn,
					smem_bytes_static,
					smem_bytes_dyn + smem_bytes_static,
					num_reg_per_thread);

				GPMESH_TRACE(
					"GPMeshStatic::check_shared_memory() available total shared "
					"memory per block = {} (bytes) = {} (Kb)",
					devProp.sharedMemPerBlock,
					float(devProp.sharedMemPerBlock) / 1024.0f);
			}

			if (smem_bytes_static + smem_bytes_dyn > devProp.sharedMemPerBlock) {
				GPMESH_ERROR(
					" GPMeshStatic::check_shared_memory() shared memory needed for"
					" input function ({} bytes) exceeds the max shared memory "
					"per block on the current device ({} bytes)",
					smem_bytes_static + smem_bytes_dyn,
					devProp.sharedMemPerBlock);
				exit(EXIT_FAILURE);
			}

			return smem_bytes_dyn;
		}

		void output_cpu_localMats(std::vector<std::vector<double>>& localMats,
			std::string                       filename)
		{
			std::string output_path = STRINGIFY(OUTPUT_DIR) "record/";
			output_path += filename;

			FILE* pFile;

			const char* name = output_path.c_str();

			fopen_s(&pFile, name, "w");

			// std::ofstream fout(output_path, std::ios::out);

			vector<int> iterGroup = { 0, 64, 128 };

			for (auto i : iterGroup) {
				for (int j = 0; j < 64 * 3; ++j) {
					// fout << "localMats(" << i << "," << j << ")= " <<
					// localMats[i][j] << std::endl;
					fprintf(
						pFile, "localMats(%d,%d)=%.5f\n", i, j, localMats[i][j]);
				}
			}

			fclose(pFile);
			// fout.close();
		}

		template <typename T>
		void output_cpu_localMats(Eigen::MatrixX<T>& localMats,
			std::string        filename)
		{
			std::string output_path = STRINGIFY(OUTPUT_DIR) "record/";
			output_path += filename;

			FILE* pFile;

			const char* name = output_path.c_str();

			fopen_s(&pFile, name, "w");

			// std::ofstream fout(output_path, std::ios::out);

			vector<int> iterGroup = { 0, 64, 128 };

			for (auto i : iterGroup) {
				for (int j = 0; j < 64 * 3; ++j) {
					// fout << "localMats(" << i << "," << j << ")= " <<
					// localMats[i][j] << std::endl;
					fprintf(
						pFile, "localMats(%d,%d)=%.5f\n", i, j, localMats(i, j));
				}
			}

			fclose(pFile);
			// fout.close();
		}

		template <typename T>
		void convert_spline_attribute_to_eigen_matrix(
			gpmesh::SplineAttributeBase<T>& spline_attr,
			Eigen::MatrixXd& out)
		{
			out.resize(m_num_cpoints * m_dims, 1);
			out.setZero();

			for (uint32_t g = 0; g < m_num_geometry; ++g) {
				for (uint16_t local_idx = 0; local_idx < 64; ++local_idx) {
					uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
					for (uint32_t d = 0; d < m_dims; ++d) {
						out(global_idx + d * m_num_cpoints, 0) +=
							spline_attr.operator()(g, local_idx, d);
					}
				}
			}
		}

		template <typename T>
		void pass_eigen_matrix_data_to_spline_attribute(
			gpmesh::SplineAttributeBase<T>& out,
			Eigen::MatrixXd& in_coefs,
			Eigen::MatrixXd& in_cps)
		{
			for (uint32_t g = 0; g < m_num_geometry; ++g) {
				for (uint16_t local_id = 0; local_id < 64; ++local_id) {
					uint32_t global_id = m_h_geometry_ltog_cp[g][local_id];
					auto     local_coefs_ids = m_h_geometry_gtol_cp[global_id];

					bool isFind = false;

					for (auto local_coefs_id : local_coefs_ids) {
						if (local_id == local_coefs_id % 64 &&
							g == local_coefs_id / 64) {
							isFind = true;
						}
					}

					assert(isFind);

					for (uint32_t d = 0; d < m_dims; ++d) {

						assert(in_cps(global_id + d * m_num_cpoints, 0) ==
							in_coefs(g * 64 + local_id, d));
						out.operator()(g, local_id, d) =
							in_cps(global_id + d * m_num_cpoints);
					}
				}
			}
		}

		template <typename T>
		void hepler_pass_eigen_cps_to_spline_attribute_host(
			gpmesh::SplineAttributeBase<T>& out,
			Eigen::MatrixX<T>& in_cps)
		{
			for (uint32_t g = 0; g < m_num_geometry; ++g) {
				for (uint16_t local_id = 0; local_id < 64; ++local_id) {
					uint32_t global_id = m_h_geometry_ltog_cp[g][local_id];
					auto     local_coefs_ids = m_h_geometry_gtol_cp[global_id];

					bool isFind = false;

					for (auto local_coefs_id : local_coefs_ids) {
						if (local_id == local_coefs_id % 64 &&
							g == local_coefs_id / 64) {
							isFind = true;
						}
					}

					assert(isFind);

					for (uint32_t d = 0; d < m_dims; ++d) {
						out.operator()(g, local_id, d) =
							in_cps(global_id + d * m_num_cpoints, 0);
					}
				}
			}
		}

		template <typename T>
		void helper_output_matrix(string& filename, Eigen::MatrixXd data)
		{
			std::string output_path = STRINGIFY(OUTPUT_DIR) "record/";
			output_path += filename;

			FILE* pFile;

			const char* name = output_path.c_str();

			fopen_s(&pFile, name, "w");

			// std::ofstream fout(output_path, std::ios::out);

			vector<int> iterGroup = { 0, 64, 128 };

			for (size_t i = 0; i < data.rows(); ++i)
			{
				fprintf(pFile, "Row-%d: ", i);
				for (size_t j = 0; j < data.cols(); ++j)
				{
					fprintf(pFile, "%.5f= ", data(i, j));
				}
				fprintf(pFile, "\n");
			}

			fclose(pFile);
			// fout.close();
		}

		template <typename T>
		void helper_convert_tripVec_to_cps(
			std::vector<std::vector<std::vector<T>>>& tripleVec,
			Eigen::MatrixX<T>& cps)
		{
			cps.setZero(m_dims * m_num_cpoints, 1);

			for (int i = 0; i < tripleVec.size(); ++i) {
				for (int j = 0; j < tripleVec[i].size(); ++j) {
					int global_index = m_h_geometry_ltog_cp[i][j];

					for (int k = 0; k < tripleVec[i][j].size(); ++k) {
						cps(global_index + k * m_num_cpoints, 0) =
							tripleVec[i][j][k];
					}
				}
			}
		}

		template <typename T>
		void helper_convert_cps_to_tripleVec(
			Eigen::MatrixX<T>& cps,
			std::vector<std::vector<std::vector<T>>>& tripleVec)
		{
			tripleVec.resize(m_num_geometry);

			for (int i = 0; i < m_num_geometry; ++i) {
				tripleVec[i].resize(64);

				for (int j = 0; j < 64; ++j) {
					int global_index = m_h_geometry_ltog_cp[i][j];

					tripleVec[i][j].resize(m_dims);

					for (int d = 0; d < m_dims; ++d) {
						tripleVec[i][j][d] =
							cps(global_index + d * m_num_cpoints, 0);
					}
				}
			}
		}

		template <typename T>
		void conjgrad(const Eigen::SparseMatrix<double>& A,
			const Eigen::VectorXd& b,
			Eigen::VectorXd& x)
		{
			Eigen::VectorXd r = b - A * x;
			Eigen::VectorXd p = r;
			double          rsold = r.squaredNorm();
			double          tol = 1e-6;

			double rs_0 = rsold;

			for (int i = 0; i < b.size(); ++i) {
				Eigen::VectorXd Ap = A * p;

				double alpha = rsold / (p.template dot(Ap));

				x.noalias() = x + alpha * p;
				r.noalias() = r - alpha * Ap;
				double rsnew = r.squaredNorm();

				if (rsnew < tol * tol * rs_0) {
					GPMESH_TRACE("conjgrad solve success for iter {}", i);
					break;
				}

				p = r + (rsnew / rsold) * p;
				rsold = rsnew;
			}

			GPMESH_TRACE("conjgrad solve result rs={}", rsold);
		}

		template <typename T>
		void preconditioned_conjgrad(const Eigen::SparseMatrix<double>& A,
			const Eigen::VectorXd& b,
			Eigen::VectorXd& x)
		{
			Eigen::VectorXd r = b - A * x;

			Eigen::VectorXd Preconditioner;
			Preconditioner = helper_get_preconditioner_from_matrix(A);
			Eigen::VectorXd p = Preconditioner.array() * r.array();

			{
				for (int i = 0; i < p.size(); ++i)
				{
					double val = p(i) - Preconditioner(i) * r(i);
					if (abs(val) > 1e-15)
					{
						GPMESH_ERROR("CWise Operation failed!");
					}
				}
			}

			double          rsold = p.dot(r);
			double          tol = 1e-6;

			double rs_0 = rsold;

			for (int i = 0; i < b.size(); ++i) {
				Eigen::VectorXd Ap = A * p;

				double alpha = rsold / (p.template dot(Ap));

				x.noalias() = x + alpha * p;
				r.noalias() = r - alpha * Ap;

				Eigen::VectorXd s = Preconditioner.array() * r.array();

				double rsnew = r.dot(s);

				if (rsnew < tol * tol * 1e-1 * rs_0) {
					GPMESH_TRACE("PCG solve success for iter {}", i);
					break;
				}

				p = s + (rsnew / rsold) * p;

				rsold = rsnew;
			}

			GPMESH_TRACE("conjgrad solve result rs={}", rsold);
		}

		Eigen::VectorXd helper_get_preconditioner_from_matrix(const Eigen::SparseMatrix<double>& A)
		{
			Eigen::VectorXd M;
			assert(A.rows() == A.cols());

			M.resize(A.rows());

			for (int i = 0; i < A.rows(); ++i)
			{
				double val = A.coeff(i, i);
				if (abs(val) < 1e-10)
				{
					val = 1;
				}
				else
				{
					val = 1.0 / val;
				}

				M(i) = val;
			}

			return M;
		}
	};

}  // namespace gpmesh
