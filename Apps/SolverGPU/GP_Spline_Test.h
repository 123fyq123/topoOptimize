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

#include <filesystem>

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

namespace gpmesh {

class GPSplineSolidTest : public GPSplineSolid
{
   private:
    SplineAttributeContainer spline_attribute_container;

   public:
    GPSplineSolidTest(const GPSplineSolidTest&) = delete;

    // debug part
    explicit GPSplineSolidTest(const bool quite = false) : GPSplineSolid(){};

    template <typename T>
    std::shared_ptr<SplineAttributeBase<T>> add_spline_cpoint_attribute(
        const std::vector<std::vector<std::vector<T>>>& local_cp_attributes,
        const std::string&                              name,
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
        layoutT            layout   = SoA)
    {
        return spline_attribute_container.template add<SplineAttributeBase<T>>(
            name.c_str(), m_num_geometry, 64, num_attributes, location, layout);
    }

    template <typename T>
    std::shared_ptr<SplineAttributeGlobal<T>>
    add_spline_cpoint_attribute_global(
        const std::vector<std::vector<std::vector<T>>>& local_cp_attributes,
        const std::string&                              name,
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
                                       layoutT   layout   = SoA)
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

    void debug_IGA_test(const vector<CBSplineSolid>&     patch,
                        std::vector<boundary_condition>& m_dirichletSides,
                        std::vector<boundary_condition>& m_neumannSides)
    {
        init(patch, m_dirichletSides, m_neumannSides, false);

        constexpr uint32_t gaussianPointNum = 64;

        size_t smem_bytes_dyn = calc_shared_memory(
            (void*)IGA_debug_kernel<double, gaussianPointNum>);
        size_t geomtryNum = m_geometry_size;

        GPUTimer timer;
        timer.start();

        IGA_debug_kernel<double, gaussianPointNum>
            <<<geomtryNum, gaussianPointNum, smem_bytes_dyn>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());
        timer.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_kernel Run Success!");
        GPMESH_INFO("Kernel took {} (ms)", timer.elapsed_millis());

        // move data from gpu to cpu

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        std::vector<int> localIndex(64, 0);

        std::iota(localIndex.begin(), localIndex.end(), 0);

        std::vector<Eigen::MatrixXd> localMatsVec_gpu;
        localMatsVec_gpu.resize(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);

            m_sparse_system.pushToMatrix(g, localIndex, localMats);

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
            localMatsVec_gpu[g] = localMats_v2;
        }

        auto Matrix1 = m_sparse_system.test_get_bigmatrix();
        // test matrix
        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // stiff.assemble();
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();

        std::vector<Eigen::MatrixXd> test_localMatsVec;
        stiff.test_stiffMatrix(test_localMatsVec);

        // 计算rhs
        stiff.rhsVector();

        auto Matrix2 = stiff.system().test_get_bigmatrix();

        for (int i = 0; i < Matrix1.size(); ++i) {
            for (int j = 0; j < Matrix1[i].size(); ++j) {
                assert((Matrix1[i][j] - Matrix2[i][j]) < 0.0001);
            }
        }
        GPMESH_INFO("CPU and GPU Big Matrix Test Success!");

        // cpu version Solver for Linear Elastic Equation
        {
            m_sparse_system.dataPush();

            // Try "ConjugateGradient.h" Version

            uint32_t                    i = 0;
            Eigen::SparseMatrix<double> A;
            A = m_sparse_system.matrix();
            Eigen::VectorXd B;
            B.noalias() = stiff.getRhs();


            Eigen::VectorXd little_x;
            little_x.resize(m_dims * m_num_cpoints);
            little_x.setZero();

            conjgrad<double>(A, B, little_x);

            Eigen::SparseMatrix<double> AT = A.transpose();

            for (int i = 0; i < AT.rows(); ++i) {
                for (int j = 0; j < AT.cols(); ++j) {
                    if (A.coeff(i, j) - A.coeff(j, i) < 0.0001) {

                    } else {
                        GPMESH_ERROR(
                            "noequal for i-%d, j-%d with A(i,j) = {}, A(j, i) "
                            "= {}",
                            i,
                            j,
                            A.coeff(i, j),
                            A.coeff(j, i));
                    }
                }
            }

            Eigen::VectorXd X;
            X.resize(m_dims * m_num_cpoints);
            X.setZero();

            Eigen::VectorXd r, d;

            r.noalias() = B - A * X;
            d.noalias() = r;

            double Delta_new, Delta_0, Delta_old;

            Delta_new = r.squaredNorm();

            Delta_0 = Delta_new;

            double tol = 1e-6;

            while (i < 10) {
                Eigen::VectorXd q = A * d;

                double Alpha = Delta_new / (d.dot(q));

                X += Alpha * d;
                r -= Alpha * q;

                Delta_old = Delta_new;

                Delta_new = r.squaredNorm();

                double Beta = Delta_new / Delta_old;

                d.noalias() = r + Beta * d;

                i++;

                if (Delta_new < tol * Delta_0) {
                    GPMESH_INFO(
                        "ConjugateGradient.h copy work find for {} iters", i);
                    break;
                }
            }

            // Simple LDLT Solver

            CPUTimer timer_cpu;

            timer_cpu.start();

            Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
                m_sparse_system.matrix());
            Eigen::MatrixXd u_ldlt = solver.solve(stiff.getRhs());

            timer_cpu.stop();
            GPMESH_INFO(
                "Eigen SimplicaialLDLT solver take {} (ms) to solve Ku=F",
                timer_cpu.elapsed_millis());

            std::vector<double> u_vec(m_dims * m_num_cpoints, 0);

            for (int i = 0; i < m_dims * m_num_cpoints; ++i) {
                u_vec[i] = u_ldlt(i, 0);
            }

            // std::cout << u_vec[m_dims * m_num_cpoints - 1];

            std::string output_path = STRINGIFY(OUTPUT_DIR) "IGADataOutput/H";

            stiff.WriteParaview(const_cast<vector<CBSplineSolid>&>(patch),
                                output_path + "_CPU_LinearSolver",
                                u_ldlt,
                                1000,
                                false);


            // Eigen Conjugate Gradient Solver
            timer_cpu.start();

            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(m_sparse_system.matrix());
            Eigen::MatrixXd u_cg = cg.solve(stiff.getRhs());

            timer_cpu.stop();
            GPMESH_INFO("Eigen CG solver take {} (ms) to solve Ku=F",
                        timer_cpu.elapsed_millis());

            assert((u_ldlt - u_cg).norm() < 0.001);

            // CG Solver (copy from MCG.cu)
            Eigen::VectorXd X_cps;
            X_cps.resize(m_dims * m_num_cpoints);
            X_cps.setZero();

            Eigen::VectorXd B_cps, R_cps(X_cps), P_cps(X_cps), S_cps(X_cps);


            assert(stiff.getRhs().size() == m_dims * m_num_cpoints);

            B_cps.resize(stiff.getRhs().size());

            for (int i = 0; i < stiff.getRhs().size(); ++i) {
                B_cps(i) = stiff.getRhs()(i, 0);
            }

            // S_cps = A * X_cps;
            m_sparse_system.matrix();
            S_cps = m_sparse_system.matrix() * X_cps;

            bool bVerify = true;

            R_cps           = B_cps - S_cps;
            P_cps.noalias() = R_cps;

            double       delta_new = R_cps.dot(R_cps);
            double       delta_old = 0.0;
            const double delta_0(delta_new);

            float start_residual = delta_0;

            uint32_t iter = 0;

            while (iter < 10) {
                S_cps = m_sparse_system.matrix() * P_cps;

                double alpha = S_cps.dot(P_cps);

                alpha = delta_new / alpha;

                X_cps = X_cps + alpha * P_cps;

                if (iter % 50 == 0) {
                    R_cps = B_cps - m_sparse_system.matrix() * X_cps;
                } else {
                    R_cps = R_cps - alpha * S_cps;
                }


                GPMESH_INFO("Delta_new {} in iter {}", delta_new, iter);

                delta_old = delta_new;

                delta_new = R_cps.dot(R_cps);

                double beta = delta_new / delta_old;

                if (delta_new < 1e-6 * 1e-6 * delta_0) {
                    break;
                }

                P_cps = R_cps + beta * P_cps;

                iter++;
            }

            GPMESH_INFO("CG solver take {} iteration to solve Ku=F", iter);
        }

        spline_attribute_container.add<SplineAttributeBase<double>>(
            "X", m_num_geometry, 64, 3, LOCATION_ALL, AoS);

        Eigen::MatrixXd displacement_cps;
        Eigen::MatrixXd displacement_coefs;

        displacement_cps.resize(m_num_cpoints * 3, 1);
        displacement_cps.setZero();
        m_sparse_system.constructCoefs(displacement_cps, displacement_coefs);

        std::vector<std::vector<std::vector<double>>> displacement_tripleVec;
        convert_displacement_coefs_to_tripleVec<double>(displacement_coefs,
                                                        displacement_tripleVec);

        auto displacement_attrib_ptr = this->add_spline_cpoint_attribute(
            displacement_tripleVec, "displacement");

        auto X = displacement_attrib_ptr;

        auto S = this->add_spline_cpoint_attribute<double>(
            "S", 3, LOCATION_ALL, gpmesh::SoA);

        S->reset(0.0, LOCATION_ALL);

        const int localCPointsNum = 64;

        timer.start();
        IGA_debug_matvec_A_Multiply<double, localCPointsNum>
            <<<geomtryNum, localCPointsNum>>>(
                get_spline_context(), get_assembly_context(), *X, *S);
        timer.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_matVec run success!");
        GPMESH_INFO("Kernel took {} (ms)", timer.elapsed_millis());

        //        for (int g = 0; g < m_num_geometry; ++g) {
        //            for (int local_idx = 0; local_idx < localCPointsNum;
        //            ++local_idx) {
        //                for (int d = 0; d < 3; ++d) {
        //                    assert(abs(S->operator()(g, local_idx, d) - 0) <
        //                    0.0001);
        //                }
        //            }
        //        }

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        for (int g = 0; g < m_num_geometry; ++g) {
            for (int local_idx = 0; local_idx < localCPointsNum; ++local_idx) {
                for (int d = 0; d < 3; ++d) {
                    if (abs(S->operator()(g, local_idx, d) - 0) > 0.0001) {
                        std::cout << "Error here : " << g << " " << local_idx
                                  << " " << d << std::endl;
                    }
                    assert(abs(S->operator()(g, local_idx, d) - 0) < 0.0001);
                }
            }
        }

        // Todo:
        // epilogue(local to global) and then prologue(global to local)
        //  r = b - s = b - Ax
        //  p = rk
        //  delta_new = <r,r>
        //  delta_0 = delta_new
        //  while(iterNum < setLimit)

        // 1. epilogue
        // cpu version

        Eigen::MatrixXd b_cps = stiff.system().rhs();

        Eigen::MatrixXd S_cps;
        S_cps.resize(
            m_sparse_system.getDim() * m_sparse_system.getControlNums(), 1);
        S_cps.setZero();

        for (int g = 0; g < m_num_geometry; ++g) {
            for (uint16_t local_idx = 0; local_idx < localCPointsNum;
                 ++local_idx) {
                uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
                for (uint32_t d = 0; d < 3; ++d) {
                    S_cps(global_idx + d * m_sparse_system.getControlNums()) +=
                        S->operator()(g, local_idx, d);
                }
            }
        }

        // calculate the differences

        Eigen::MatrixXd r_cps = b_cps - S_cps;
        Eigen::MatrixXd p_cps = r_cps;

        Eigen::MatrixXd p_coefs;
        m_sparse_system.constructCoefs(p_cps, p_coefs);
        std::vector<std::vector<std::vector<double>>> p_tripleVec;
        convert_displacement_coefs_to_tripleVec(p_coefs, p_tripleVec);

        auto P = this->add_spline_cpoint_attribute(p_tripleVec, "P");

        double delta_new = 0;
        for (int i = 0; i < r_cps.rows(); ++i) {
            delta_new += r_cps(i, 0) * r_cps(i, 0);
        }

        double delta_0 = delta_new;

        uint32_t num_cg_iter_taken = 0;
        uint32_t max_num_cg_iter   = 50;

        double alpha(0.0);
        double delta_old(0.0);
        double beta(0.0);

        while (num_cg_iter_taken < max_num_cg_iter) {

            // S = AP : here 'P'=='d'
            S->reset(0.0, gpmesh::DEVICE);
            IGA_debug_matvec_A_Multiply<double, localCPointsNum>
                <<<geomtryNum, localCPointsNum>>>(
                    get_spline_context(), get_assembly_context(), *P, *S);

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            alpha = 0.0;

            for (int i = 0; i < S_cps.rows(); ++i) {
                alpha += S_cps(i, 0) * p_cps(i, 0);
            }

            alpha = delta_new / alpha;

            displacement_cps = displacement_cps + alpha * p_cps;

            r_cps = r_cps - alpha * S_cps;

            delta_old = delta_new;

            delta_new = 0;

            for (int i = 0; i < r_cps.rows(); ++i) {
                delta_new += r_cps(i, 0) * r_cps(i, 0);
            }

            GPMESH_TRACE(
                "Delta Value {} in iterNum {}", delta_new, num_cg_iter_taken);

            if (delta_new < 1e-6 * delta_0) {
                GPMESH_INFO("CG Solver Quit With iterNum={}",
                            num_cg_iter_taken);
                break;
            }

            beta = delta_new / delta_old;


            p_cps = beta * p_cps + r_cps;

            m_sparse_system.constructCoefs(p_cps, p_coefs);
            pass_eigen_matrix_data_to_spline_attribute(*P, p_coefs, p_cps);
            // push p_cps back to P

            num_cg_iter_taken++;
        }

        GPMESH_INFO("CG Solver Quit With iterNum={}", num_cg_iter_taken);
    }

    void debug_IGA_cg_test(const vector<CBSplineSolid>&     patch,
                           std::vector<boundary_condition>& m_dirichletSides,
                           std::vector<boundary_condition>& m_neumannSides)
    {

        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        size_t smem_bytes_dyn = calc_shared_memory(
            (void*)IGA_debug_kernel<double, gaussianPointNum>);
        size_t geomtryNum = m_geometry_size;

        GPUTimer timer_gpu;
        timer_gpu.start();
        IGA_debug_kernel<double, gaussianPointNum>
            <<<geomtryNum, gaussianPointNum, smem_bytes_dyn>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_kernel run success!");
        GPMESH_INFO("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        // move data from gpu to cpu
        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        std::vector<int> localIndex(64, 0);

        std::iota(localIndex.begin(), localIndex.end(), 0);

        std::vector<Eigen::MatrixXd> localMatsVec_gpu;
        localMatsVec_gpu.resize(m_num_geometry);

        // verify  Eigen::ColMajor
        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);

            m_sparse_system.pushToMatrix(g, localIndex, localMats);

            Eigen::MatrixXd localMats_v2;
            localMats_v2.setZero(192, 192);

            // internal test
            for (int i = 0; i < 192; ++i) {
                for (int j = 0; j < 192; ++j) {
                    localMats_v2(i, j) =
                        m_h_assembly_info[g].local_mats[i + 192 * j];
                    auto t1 = localMats_v2(i, j);
                    auto t2 = localMats(i, j);
                    assert(abs(t1 - t2) < 0.001);
                }
            }
            localMatsVec_gpu[g] = localMats_v2;
        }

        auto Matrix1 = m_sparse_system.test_get_bigmatrix();

        // reference assembler
        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();

        std::vector<Eigen::MatrixXd> test_localMatsVec;
        stiff.test_stiffMatrix(test_localMatsVec);
        // 计算rhs
        stiff.rhsVector();

        auto Matrix2 = stiff.system().test_get_bigmatrix();


        // output_cpu_localMats(localMatsVec_gpu[0], "localMats_cpu1.txt");
        // output_cpu_localMats(Matrix1, "localMats_gpu1.txt");

        for (int i = 0; i < Matrix1.size(); ++i) {
            for (int j = 0; j < Matrix1[i].size(); ++j) {
                assert((Matrix1[i][j] - Matrix2[i][j]) < 0.0001);

                if (i == 0 || i == 64 || i == 128) {
                    // printf("localMats(%d, %d)= %.5f\n", i, j, Matrix1[i][j]);
                }
            }
        }
        GPMESH_INFO("CPU and GPU Big Matrix Test Success!");

        {
            stiff.system().dataPush();

            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();

            stiff.stiffMatrix();
            // stiff.test_stiffMatrix(test_localMatsVec);

            stiff.rhsVector();

            stiff.system().dataPush();

            // notice m_sparse_system.rhs() is not calculated here
            m_sparse_system.dataPush();

            for (int i = 0; i < m_sparse_system.matrix().rows(); ++i) {
                for (int j = 0; j < m_sparse_system.matrix().cols(); ++j) {
                    assert(abs(m_sparse_system.matrix().coeffRef(i, j) -
                               stiff.system().matrix().coeff(i, j)) < 1e-6);
                }
            }

            GPMESH_INFO("CPU and GPU Sparse Matrix Test Success!");
        }


        // CG Solver Test
        {

            bool bCGTest = false;

            if (bCGTest) {
                // “ConjugateGradient.h" Version
                Eigen::VectorXd u_conjgrad;
                {
                    uint32_t                    i = 0;
                    Eigen::SparseMatrix<double> A;
                    A = m_sparse_system.matrix();
                    Eigen::VectorXd B;
                    B.noalias() = stiff.getRhs();
                    u_conjgrad.resize(m_dims * m_num_cpoints);
                    u_conjgrad.setZero();

                    conjgrad<double>(A, B, u_conjgrad);
                }


                // Eigen Preconditioned Conjugate Gradient
                Eigen::VectorXd u_cg;
                {
                    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
                    cg.compute(m_sparse_system.matrix());
                    u_cg = cg.solve(stiff.getRhs());
                    GPMESH_TRACE("Eigen CG Solver success for iter {}",
                                 cg.iterations());
                    GPMESH_TRACE("Eigen CG Solver result error={}", cg.error());
                }

                GPMESH_TRACE("Difference between u_cg and u_conjgrad = {}",
                             (u_cg - u_conjgrad).norm());
                assert((u_cg - u_conjgrad).norm() < 1e-6);
                GPMESH_INFO("conjgrad Solver Test Success!");
            }

            // GPU version


            Eigen::MatrixXd Disp_cps;

            Disp_cps.setZero(m_num_cpoints * m_dims, 1);

            std::vector<std::vector<std::vector<double>>> Disp_tripleVec;

            helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

            auto Disp_gpu_ptr = this->add_spline_cpoint_attribute<double>(
                Disp_tripleVec, "disp");

            auto X = Disp_gpu_ptr;

            auto S = this->add_spline_cpoint_attribute<double>(
                "S", 3, LOCATION_ALL, gpmesh::SoA);
            S->reset(0.0, LOCATION_ALL);

            timer_gpu.start();

            // CG First Step

            IGA_debug_matvec_A_Multiply<double, localCPointNum>
                <<<geomtryNum, localCPointNum>>>(
                    get_spline_context(), get_assembly_context(), *X, *S);
            timer_gpu.stop();
            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());

            GPMESH_INFO("IGA_debug_matVec run success!");
            GPMESH_INFO("GPU Kernel took {} (ms)", timer_gpu.elapsed_millis());

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            // Test for first matvec
            {
                for (int g = 0; g < m_num_geometry; ++g) {
                    for (int local_idx = 0; local_idx < localCPointNum;
                         ++local_idx) {
                        for (int d = 0; d < 3; ++d) {
                            if (abs(S->operator()(g, local_idx, d) - 0) >
                                0.0001) {
                                std::cout << "Error here : " << g << " "
                                          << local_idx << " " << d << std::endl;
                            }
                            assert(abs(S->operator()(g, local_idx, d) - 0) <
                                   0.0001);
                        }
                    }
                }
            }

            // First Iteration in CG
            Eigen::MatrixXd B_cps = stiff.system().rhs();

            Eigen::MatrixXd S_cps;
            S_cps.setZero(m_dims * m_num_cpoints, 1);

            // convert S(gpu) to S_cps(cpu);
            for (int g = 0; g < m_num_geometry; ++g) {
                for (uint16_t local_idx = 0; local_idx < localCPointNum;
                     ++local_idx) {
                    uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
                    for (uint32_t d = 0; d < 3; ++d) {
                        S_cps(global_idx +
                              d * m_sparse_system.getControlNums()) +=
                            S->operator()(g, local_idx, d);
                    }
                }
            }

            Eigen::MatrixXd R_cps = B_cps - S_cps;
            Eigen::MatrixXd P_cps = R_cps;

            std::vector<std::vector<std::vector<double>>> P_tripleVec;
            helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

            Eigen::MatrixXd P_verify_cps;
            helper_convert_tripVec_to_cps<double>(P_tripleVec, P_verify_cps);

            assert((P_verify_cps - P_cps).norm() < 1e-8);

            double Delta_new = 0;

            Delta_new = R_cps.squaredNorm();

            double Delta_0(Delta_new);

            GPMESH_TRACE("Delta_0 from GPU is {}", Delta_0);

            auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

            S->reset(0.0, gpmesh::DEVICE);

            IGA_debug_matvec_A_Multiply<double, localCPointNum>
                <<<geomtryNum, localCPointNum>>>(
                    get_spline_context(), get_assembly_context(), *P, *S);

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            double Alpha = 0.0;

            for (int i = 0; i < S_cps.rows(); ++i) {
                Alpha += S_cps(i, 0) * P_cps(i, 0);
            }

            Alpha = Delta_new / Alpha;

            GPMESH_TRACE("Alpha from GPU First Iteration is {}", Alpha);

            Disp_cps += Alpha * P_cps;
            R_cps -= Alpha * S_cps;

            double Delta_old = Delta_new;

            Delta_new = R_cps.squaredNorm();

            GPMESH_TRACE("Delta_new in GPU 1st Iter is {}", Delta_new);

            double Beta = Delta_new / Delta_old;

            // CPU test

            Eigen::SparseMatrix<double> A = m_sparse_system.matrix();
            Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);
            Eigen::VectorXd b = stiff.getRhs();

            CPUTimer timer_cpu;

            timer_cpu.start();
            Eigen::VectorXd r = b - A * x;
            timer_cpu.stop();

            GPMESH_INFO("CPU matvec took {} (ms)", timer_cpu.elapsed_millis());

            Eigen::VectorXd p = r;

            double delta_new = r.squaredNorm();
            double delta_0(delta_new);

            GPMESH_TRACE("delta_0 from cpu is {}", delta_new);

            Eigen::VectorXd diff_tmp = p - P_cps;

            for (int i = 0; i < m_num_cpoints * m_dims; ++i) {
                if (abs(diff_tmp(i)) > 1e-6) {
                    // std::cout << "Error here : " << i << std::endl;
                    printf("Error here : %d for p-%f, P_cps-%f\n",
                           i,
                           p(i),
                           P_cps(i, 0));

                    // obvious this test is past
                    // because the test is equaly to "stiff.getRhs() ==
                    // stiff.getRhs()"
                } else {
                    assert(abs(diff_tmp(i)) < 1e-6);
                }
            }

            // move gpu code here to keep P_cps in diff calculation
            P_cps = R_cps + Beta * P_cps;
            hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
            P->move(gpmesh::HOST, gpmesh::DEVICE);

            Eigen::VectorXd Ap = A * p;

            diff_tmp = Ap - S_cps;

            for (int i = 0; i < m_num_cpoints * m_dims; ++i) {
                if (abs(diff_tmp(i)) > 1e-6) {
                    printf("Error here : %d for Ap-%f, S_cps-%f\n",
                           i,
                           Ap(i),
                           S_cps(i, 0));
                }
                assert(diff_tmp(i) < 1e-6);
            }

            GPMESH_INFO("CPU and GPU \"A * p\" test Success!");

            double alpha = delta_new / p.dot(Ap);

            GPMESH_TRACE("alpha from cpu First Iteration is {}", alpha);

            x += alpha * p;
            r -= alpha * Ap;

            double delta_old = delta_new;
            delta_new        = r.squaredNorm();

            GPMESH_TRACE("delta_new in cpu 1st Iter is {}", delta_new);

            double beta = delta_new / delta_old;

            p = r + beta * p;

            assert((x - Disp_cps).norm() < 1e-6);

            // CPU and GPU Solver Iteration

            uint32_t iterNum = 1;

            while (iterNum < 20) {
                iterNum++;
                // GPU Iter

                S->reset(0.0, gpmesh::DEVICE);
                IGA_debug_matvec_A_Multiply<double, localCPointNum>
                    <<<geomtryNum, localCPointNum>>>(
                        get_spline_context(), get_assembly_context(), *P, *S);

                S->move(gpmesh::DEVICE, gpmesh::HOST);

                convert_spline_attribute_to_eigen_matrix(*S, S_cps);

                Alpha = 0.0;

                for (int i = 0; i < S_cps.rows(); ++i) {
                    Alpha += S_cps(i, 0) * P_cps(i, 0);
                }

                Alpha = Delta_new / Alpha;

                Disp_cps += Alpha * P_cps;

                R_cps -= Alpha * S_cps;

                Delta_old = Delta_new;

                Delta_new = R_cps.squaredNorm();

                Beta = Delta_new / Delta_old;

                P_cps = R_cps + Beta * P_cps;

                GPMESH_TRACE(
                    "GPU Iteration {} : Alpha={}, Delta_new={}, Beta={}",
                    iterNum,
                    Alpha,
                    Delta_new,
                    Beta);

                hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
                P->move(gpmesh::HOST, gpmesh::DEVICE);

                // CPU Iter
                Ap = A * p;

                alpha = delta_new / p.dot(Ap);

                x += alpha * p;
                r -= alpha * Ap;

                delta_old = delta_new;
                delta_new = r.squaredNorm();

                beta = delta_new / delta_old;

                p = r + beta * p;

                GPMESH_TRACE(
                    "CPU Iteration {} : alpha={}, delta_new={}, beta={}",
                    iterNum,
                    alpha,
                    delta_new,
                    beta);

                // assert((x - Disp_cps).norm() < 1e-6);
                // assert((p - P_cps).norm() < 1e-6);
            }
        }
    }

    void debug_IGA_dirichlet_test(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // stiff.assemble();

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        // stiff.stiffMatrix();
        stiff.rhsVector();

        size_t smem_bytes_dyn = calc_shared_memory(
            (void*)IGA_debug_kernel<double, gaussianPointNum>);
        size_t geomtryNum = m_geometry_size;

        GPUTimer timer_gpu;
        timer_gpu.start();
        IGA_debug_kernel<double, gaussianPointNum>
            <<<geomtryNum, gaussianPointNum, smem_bytes_dyn>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_kernel run success!");
        GPMESH_INFO("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        // move data from gpu to cpu
        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        std::vector<Eigen::MatrixXd> localMatsVec_gpu;
        localMatsVec_gpu.resize(m_num_geometry);

        std::vector<int> localIndex(64, 0);

        std::iota(localIndex.begin(), localIndex.end(), 0);

        // verify  Eigen::ColMajor
        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);

            m_sparse_system.pushToMatrix(g, localIndex, localMats);

            Eigen::MatrixXd localMats_v2;
            localMats_v2.setZero(192, 192);

            // internal test
            for (int i = 0; i < 192; ++i) {
                for (int j = 0; j < 192; ++j) {
                    localMats_v2(i, j) =
                        m_h_assembly_info[g].local_mats[i + 192 * j];
                    auto t1 = localMats_v2(i, j);
                    auto t2 = localMats(i, j);
                    assert(abs(t1 - t2) < 0.001);
                }
            }
            localMatsVec_gpu[g] = localMats_v2;
        }

        auto Matrix1 = m_sparse_system.test_get_bigmatrix();

        std::vector<Eigen::MatrixXd> test_localMatsVec;
        stiff.test_stiffMatrix(test_localMatsVec);

        auto Matrix2 = stiff.system().test_get_bigmatrix();

        for (int i = 0; i < Matrix1.size(); ++i) {
            for (int j = 0; j < Matrix1[i].size(); ++j) {
                assert((Matrix1[i][j] - Matrix2[i][j]) < 0.0001);

                if (i == 0 || i == 64 || i == 128) {
                    // printf("localMats(%d, %d)= %.5f\n", i, j, Matrix1[i][j]);
                }
            }
        }
        GPMESH_INFO("CPU and GPU Big Matrix Test Success!");

        // test dirichlet boundary info
        {
            assert(m_dirichletSides.size() ==
                   m_dirichlet_context.m_num_dirichlet_boundary);

            DirichletInfo* m_dirichlet_info = (DirichletInfo*)malloc(
                m_dirichlet_context.m_num_dirichlet_boundary *
                sizeof(DirichletInfo));

            CUDA_ERROR(cudaMemcpy(m_dirichlet_info,
                                  m_dirichlet_context.m_d_dirichlet_info,
                                  m_dirichlet_context.m_num_dirichlet_boundary *
                                      sizeof(DirichletInfo),
                                  cudaMemcpyDeviceToHost));

            for (int i = 0; i < m_dirichlet_context.m_num_dirichlet_boundary;
                 ++i) {
                assert(m_dirichletSides[i].getDir() == m_dirichlet_info[i].dim);
                assert(m_dirichletSides[i].getValue() ==
                       m_dirichlet_info[i].value);
                assert(m_dirichletSides[i].side().patch() ==
                       m_dirichlet_info[i].geometry_id);
                assert(static_cast<uint32_t>(
                    m_dirichletSides[i].side().index() ==
                    static_cast<uint32_t>(m_dirichlet_info[i].side)));
                assert(m_dirichlet_info[i].dirichlet_boundary_id == i);
            }

            GPMESH_INFO("Dirichlet Context Test Success!");

            const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
            const int GeometrySurfNum = 16;

            IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
                get_spline_context(),
                get_assembly_context(),
                get_dirichlet_context());

            GPMESH_INFO("Add Dirichlet Boundary Info to Mats");


            stiff.system().dataPush_Dirichlet();
            // stiff.system().dataPush();

            Eigen::VectorXd u_stiff;

            Eigen::SparseMatrix<double> A_stiff;
            A_stiff = stiff.getMatrix();
            Eigen::VectorXd B_stiff;
            B_stiff = stiff.getRhs();

            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;

            cg.compute(A_stiff);

            u_stiff = cg.solve(B_stiff);

            if (cg.info() != Eigen::Success) {
                GPMESH_INFO("Eigen CG Solver Failed for iter {} and error {}",
                            cg.iterations(),
                            cg.error());
            } else {
                GPMESH_TRACE("Eigen CG Solver success for iter {}",
                             cg.iterations());
                GPMESH_TRACE("Eigen CG Solver result error={}", cg.error());
            }

            // Test for stiff.system().matrix with local_mat assembly;
            Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;

            m_assembly_context.update_host_info(m_h_assembly_info,
                                                m_num_geometry);

            auto m_p_hash = [&](const pair<uint32_t, uint32_t>& p) -> size_t {
                return p.first + p.second * m_dims * m_num_cpoints;
            };


            // assembly local mats to global mats
            unordered_map<pair<uint32_t, uint32_t>, double, decltype(m_p_hash)>
                global_mat_map(m_num_geometry * 192 * 192, m_p_hash);

            for (uint32_t geom_id = 0; geom_id < m_num_geometry; ++geom_id) {
                for (uint16_t local_i = 0; local_i < localCPointNum;
                     ++local_i) {
                    uint32_t global_i = m_h_geometry_ltog_cp[geom_id][local_i];
                    for (uint16_t i_d = 0; i_d < 3; ++i_d) {
                        for (uint16_t local_j = 0; local_j < localCPointNum;
                             ++local_j) {
                            uint32_t global_j =
                                m_h_geometry_ltog_cp[geom_id][local_j];
                            for (uint16_t j_d = 0; j_d < 3; ++j_d) {
                                pair<int, int> key =
                                    make_pair(global_i + i_d * m_num_cpoints,
                                              global_j + j_d * m_num_cpoints);

                                global_mat_map[key] +=
                                    m_h_assembly_info[geom_id]
                                        .local_mats[(local_i + i_d * 64) +
                                                    192 * (local_j + j_d * 64)];
                            }
                        }
                    }
                }
            }

            std::vector<Eigen::Triplet<double>> triplets_vec;

            for (auto& pr : global_mat_map) {
                triplets_vec.emplace_back(
                    pr.first.first, pr.first.second, pr.second);
            }

            Eigen::SparseMatrix<double> SparseMatrixFromGPU;
            SparseMatrixFromGPU.resize(m_dims * m_num_cpoints,
                                       m_dims * m_num_cpoints);
            SparseMatrixFromGPU.setFromTriplets(triplets_vec.begin(),
                                                triplets_vec.end());
            SparseMatrixFromGPU.makeCompressed();


            double diff_s_m =
                (GroundTruthSparseMatrix - SparseMatrixFromGPU).norm();

            // find Accuracy Problem Here

            //            for (int i = 0; i < GroundTruthSparseMatrix.rows();
            //            ++i)
            //            {
            //                for (int j = 0; j <
            //                GroundTruthSparseMatrix.cols(); ++j)
            //                {
            //                    double local_diff =
            //                    GroundTruthSparseMatrix.coeff(i, j) -
            //                    SparseMatrixFromGPU.coeff(i, j);
            //
            //                    if(local_diff > 1e-16)
            //                    {
            //                        GPMESH_INFO("Error Here : diff = {} for
            //                        ({},{}) = {} and {}", local_diff, i, j,
            //                        GroundTruthSparseMatrix.coeff(i, j),
            //                        SparseMatrixFromGPU.coeff(i, j));
            //                    }
            //                }
            //            }

            GPMESH_INFO("Sparse Matrix diff = {}", diff_s_m);
        }
    }

    // try to use CG to solve Ax=b
    void debug_IGA_cg_solve_test_v0(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        size_t smem_bytes_dyn = calc_shared_memory(
            (void*)IGA_debug_kernel<double, gaussianPointNum>);
        size_t geomtryNum = m_geometry_size;

        GPUTimer timer_gpu;
        timer_gpu.start();
        IGA_debug_kernel<double, gaussianPointNum>
            <<<geomtryNum, gaussianPointNum, smem_bytes_dyn>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_kernel Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;

        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        stiff.stiffMatrix();
        // 在cpu端计算rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        bool bCPU_CG_Test = true;

        Eigen::VectorXd disp_pre_set;
        Eigen::VectorXd GroundTruthResult;
        if (bCPU_CG_Test) {
            Eigen::VectorXd u_conjgrad;
            {
                uint32_t                    i = 0;
                Eigen::SparseMatrix<double> A;
                A = stiff.system().matrix();
                Eigen::VectorXd B;
                B.noalias() = stiff.getRhs();
                u_conjgrad.resize(m_dims * m_num_cpoints);
                u_conjgrad.setZero();

                conjgrad<double>(A, B, u_conjgrad);
            }

            // Eigen Preconditioned Conjugate Gradient
            Eigen::VectorXd u_cg;
            {
                Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
                // cg.setMaxIterations(50);
                cg.compute(stiff.system().matrix());
                u_cg = cg.solve(stiff.getRhs());
                if (cg.info() == Eigen::NoConvergence ||
                    cg.info() == Eigen::NumericalIssue) {
                    GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                    GPMESH_TRACE(
                        "Iteration {}, Error {}", cg.iterations(), cg.error());
                } else {
                    GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                    GPMESH_TRACE(
                        "Iteration {}, Error {}", cg.iterations(), cg.error());
                }
            }

            GPMESH_INFO("Difference between u_cg and u_conjgrad = {}",
                        (u_cg - u_conjgrad).norm());
            // assert((u_cg - u_conjgrad).norm() < 1e-6);
            GPMESH_INFO("conjgrad Solver Test Success!");

            disp_pre_set      = u_cg;
            GroundTruthResult = u_cg;
        }

        Eigen::MatrixXd Disp_cps;

        if (bCPU_CG_Test) {
            // Disp_cps = disp_pre_set;
            Disp_cps.setZero(m_num_cpoints * m_dims, 1);
        } else {
            Disp_cps.setZero(m_num_cpoints * m_dims, 1);
        }

        std::vector<std::vector<std::vector<double>>> Disp_tripleVec;

        helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

        auto Disp_gpu_ptr =
            this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp");

        auto X = Disp_gpu_ptr;

        auto S = this->add_spline_cpoint_attribute<double>(
            "S", 3, LOCATION_ALL, gpmesh::SoA);
        S->reset(0.0, LOCATION_ALL);

        timer_gpu.start();

        // GPU CG Step 1: S = A * X

        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geomtryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *X, *S);

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_matvec_A_Multiply Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        Eigen::MatrixXd B_cps = stiff.system().rhs();

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

        std::vector<std::vector<std::vector<double>>> P_tripleVec;
        helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

        double Delta_new = 0;

        Delta_new = R_cps.squaredNorm();

        double Delta_0(Delta_new);

        GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

        auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

        S->reset(0.0, gpmesh::DEVICE);


        // GPU CG Step 2: First Iteration
        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geomtryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *P, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        convert_spline_attribute_to_eigen_matrix(*S, S_cps);

        double Alpha = 0.0;

        for (int i = 0; i < S_cps.rows(); ++i) {
            Alpha += S_cps(i, 0) * P_cps(i, 0);
        }

        GPMESH_TRACE("Delta_new in CPU iteration - {} is {}", 1, Delta_new);
        GPMESH_TRACE("P.dot(AP) in CPU iteration-{} is {}", 1, Alpha);

        Alpha = Delta_new / Alpha;

        GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

        Disp_cps += Alpha * P_cps;
        R_cps -= Alpha * S_cps;

        double Delta_old = Delta_new;

        Delta_new = R_cps.squaredNorm();

        GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);

        double Beta = Delta_new / Delta_old;

        GPMESH_TRACE("Beta in GPU iteration-{} is {}", 1, Beta);

        P_cps = R_cps + Beta * P_cps;
        hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
        P->move(gpmesh::HOST, gpmesh::DEVICE);

        // CPU CG Step 1:

        Eigen::SparseMatrix<double> A = stiff.system().matrix();
        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);
        if (bCPU_CG_Test) {
            // x = disp_pre_set;
        } else {
        }
        Eigen::VectorXd b = stiff.getRhs();

        CPUTimer timer_cpu;

        timer_cpu.start();
        Eigen::VectorXd r = b - A * x;
        timer_cpu.stop();

        GPMESH_INFO("Eigen r = b - A * x Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_cpu.elapsed_millis());

        Eigen::VectorXd p = r;

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);
        GPMESH_TRACE("delta_0 in CPU is {}", delta_0);

        Eigen::VectorXd Ap = A * p;

        GPMESH_TRACE("delta_new in CPU iteration - {} is {}", 1, delta_new);
        GPMESH_TRACE("p.dot(Ap) in CPU iteration-{} is {}", 1, p.dot(Ap));

        double alpha = delta_new / p.dot(Ap);
        GPMESH_TRACE("alpha in CPU iteration-{} is {}", 1, alpha);

        x += alpha * p;
        r -= alpha * Ap;

        double delta_old = delta_new;
        delta_new        = r.squaredNorm();

        GPMESH_TRACE("delta_new in CPU iteration{} is {}", 1, delta_new);

        double beta = delta_new / delta_old;
        GPMESH_TRACE("beta in CPU iteration-{} is {}", 1, beta);

        p = r + beta * p;

        uint32_t iterNum = 1;

        while (iterNum < 200) {
            iterNum++;

            S->reset(0.0, gpmesh::DEVICE);
            IGA_debug_matvec_A_Multiply<double, localCPointNum>
                <<<geomtryNum, localCPointNum>>>(
                    get_spline_context(), get_assembly_context(), *P, *S);

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            Alpha = 0.0;

            for (int i = 0; i < S_cps.rows(); ++i) {
                Alpha += S_cps(i, 0) * P_cps(i, 0);
            }
            GPMESH_TRACE("GPU Iteration {} : ", iterNum);
            GPMESH_TRACE("Delta_old={}, P.dot(AP)={}", Delta_new, Alpha);
            Alpha = Delta_new / Alpha;

            Disp_cps += Alpha * P_cps;

            R_cps -= Alpha * S_cps;

            Delta_old = Delta_new;

            Delta_new = R_cps.squaredNorm();

            Beta = Delta_new / Delta_old;

            P_cps = R_cps + Beta * P_cps;


            GPMESH_TRACE("Delta_new={}", Delta_new);

            hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
            P->move(gpmesh::HOST, gpmesh::DEVICE);

            // CPU Iter
            Ap = A * p;

            alpha = delta_new / p.dot(Ap);
            GPMESH_TRACE("CPU Iteration {} : ", iterNum);
            GPMESH_TRACE("delta_old={}, p.dot(Ap)={}", delta_new, p.dot(Ap));

            x += alpha * p;
            r -= alpha * Ap;

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            beta = delta_new / delta_old;

            p = r + beta * p;


            GPMESH_TRACE("delta_new={}", delta_new);

            if (delta_new < 1e-12 * delta_0) {
                GPMESH_INFO("CPU calculation should break in iteration {}!",
                            iterNum);
            }

            if (Delta_new < 1e-12 * Delta_0) {
                break;
            }
        }

        auto cpu_result = x;
        auto gpu_result = Disp_cps;

        GPMESH_INFO("CPU calculation error = {}",
                    (GroundTruthResult - cpu_result).norm());
        GPMESH_INFO("GPU calculation error = {}",
                    (GroundTruthResult - gpu_result).norm());
    }

    void debug_IGA_hyperelasticity_v0(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan         = 1;
        double timeStep         = 0.1;
        double density          = 3320;
        double gravitationalAcc = 9.8;

        double Beta = 0.25, Gamma = 0.5;

        double alpha1 = 1. / Beta / pow(timeStep, 2);
        double alpha2 = 1. / Beta / timeStep;
        double alpha3 = (1 - 2 * Beta) / 2 / Beta;
        double alpha4 = Gamma / Beta / timeStep;
        double alpha5 = 1 - Gamma / Beta;
        double alpha6 = (1 - Gamma / Beta / 2) * timeStep;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        int    numIterations    = 0;
        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 1, 1);
        dim3 threads_per_blocks(64, 1, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)IGA_hyperelasticity_assemble_debug_kernel<double,
                                                             gaussianPointNum>);

        for (int t = 0; t < timeSpan; ++t) {
            std::cout << t << "\n";

            // make one time step
            // 上一步的dispVector作为初始解
            // 牛顿法解
            int maxIter = 1;

            while (numIterations < maxIter) {
                // 进行计算
                timer_gpu.start();
                IGA_hyperelasticity_assemble_debug_kernel<double,
                                                          gaussianPointNum>
                    <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                        get_spline_context(),
                        get_gaussian_integral_context(),
                        get_physic_context(),
                        get_assembly_context());
                timer_gpu.stop();
                CUDA_ERROR(cudaDeviceSynchronize());
                CUDA_ERROR(cudaGetLastError());
                CUDA_ERROR(cudaProfilerStop());
                GPMESH_INFO(
                    "IGA_hyperelasticity_assemble_debug_kernel Run Success!");
                GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
                ++numIterations;
            }

            Eigen::MatrixXd newDispVector;

            Eigen::MatrixXd oldVelVector = velVector;

            velVector = alpha4 * (newDispVector - dispVector) +
                        alpha5 * oldVelVector + alpha6 * accVector;
            accVector = alpha1 * (newDispVector - dispVector) -
                        alpha2 * oldVelVector - alpha3 * accVector;
            dispVector = newDispVector;
        }
    }

    void debug_IGA_hyperelasticity_v1(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan         = 1;
        double timeStep         = 0.1;
        double density          = 3320;
        double gravitationalAcc = 9.8;

        double Beta = 0.25, Gamma = 0.5;

        double alpha1 = 1. / Beta / pow(timeStep, 2);
        double alpha2 = 1. / Beta / timeStep;
        double alpha3 = (1 - 2 * Beta) / 2 / Beta;
        double alpha4 = Gamma / Beta / timeStep;
        double alpha5 = 1 - Gamma / Beta;
        double alpha6 = (1 - Gamma / Beta / 2) * timeStep;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        int    numIterations    = 0;
        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(64, 1, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)IGA_hyperelasticity_assemble_debug_kernel<double,
                                                             gaussianPointNum>);

        for (int t = 0; t < timeSpan; ++t) {
            std::cout << t << "\n";

            // make one time step
            // 上一步的dispVector作为初始解
            // 牛顿法解
            int maxIter = 1;

            while (numIterations < maxIter) {
                // 进行计算
                timer_gpu.start();
                IGA_hyperelasticity_assembly_debug_kernel_v1<double,
                                                             gaussianPointNum>
                    <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                        get_spline_context(),
                        get_gaussian_integral_context(),
                        get_physic_context(),
                        get_assembly_context());
                timer_gpu.stop();
                CUDA_ERROR(cudaDeviceSynchronize());
                CUDA_ERROR(cudaGetLastError());
                CUDA_ERROR(cudaProfilerStop());
                GPMESH_INFO(
                    "IGA_hyperelasticity_assemble_debug_kernel Run Success!");
                GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
                ++numIterations;
            }

            Eigen::MatrixXd newDispVector;

            Eigen::MatrixXd oldVelVector = velVector;

            velVector = alpha4 * (newDispVector - dispVector) +
                        alpha5 * oldVelVector + alpha6 * accVector;
            accVector = alpha1 * (newDispVector - dispVector) -
                        alpha2 * oldVelVector - alpha3 * accVector;
            dispVector = newDispVector;
        }
    }


    void debug_IGA_hyperelasticity_v2(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan         = 1;
        double timeStep         = 0.1;
        double density          = 3320;
        double gravitationalAcc = 9.8;

        double Beta = 0.25, Gamma = 0.5;

        double alpha1 = 1. / Beta / pow(timeStep, 2);
        double alpha2 = 1. / Beta / timeStep;
        double alpha3 = (1 - 2 * Beta) / 2 / Beta;
        double alpha4 = Gamma / Beta / timeStep;
        double alpha5 = 1 - Gamma / Beta;
        double alpha6 = (1 - Gamma / Beta / 2) * timeStep;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        int    numIterations    = 0;
        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        stiffAssembler  stiff(patch, m_dirichletSides, m_neumannSides);
        Eigen::MatrixXd u;

        massAssembler mass(patch, 3320, stiff.system());

        timeIntegrator timeAssemble(stiff, mass);

        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();

        // seperate timeAssemble
        auto            dispSol = timeAssemble.displacementVector();
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        stiff.test_stiffNonlinearMatrix(displacement, true, test_localMatsVec);

        for (int t = 0; t < timeSpan; ++t) {
            std::cout << t << "\n";

            // make one time step
            // 上一步的dispVector作为初始解
            // 牛顿法解
            int maxIter = 1;

            while (numIterations < maxIter) {
                // 进行计算
                timer_gpu.start();
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>
                    <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                        get_spline_context(),
                        get_gaussian_integral_context(),
                        get_physic_context(),
                        get_assembly_context());
                timer_gpu.stop();
                CUDA_ERROR(cudaDeviceSynchronize());
                CUDA_ERROR(cudaGetLastError());
                CUDA_ERROR(cudaProfilerStop());
                GPMESH_INFO(
                    "IGA_hyperelasticity_assembly_debug_kernel_v2 Run "
                    "Success!");
                GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
                ++numIterations;
            }

            m_assembly_context.update_host_info(m_h_assembly_info,
                                                m_num_geometry);

            std::vector<int> localIndex(64, 0);

            std::iota(localIndex.begin(), localIndex.end(), 0);

            std::vector<Eigen::MatrixXd> localMatsVec_gpu;
            localMatsVec_gpu.resize(m_num_geometry);

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
                localMatsVec_gpu[g] = localMats_v2;
            }


            Eigen::MatrixXd newDispVector;

            Eigen::MatrixXd oldVelVector = velVector;

            velVector = alpha4 * (newDispVector - dispVector) +
                        alpha5 * oldVelVector + alpha6 * accVector;
            accVector = alpha1 * (newDispVector - dispVector) -
                        alpha2 * oldVelVector - alpha3 * accVector;
            dispVector = newDispVector;
        }
    }

    void helper_performance_test_for_matvec_A_multiply()
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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
            <<<geometryNum, localCPointNum>>>(
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
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
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
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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

        IGA_matvec_A_multiply_v2<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X, *S);

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_matvec_A_multiply_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
    }

    void debug_solving_assembly_rhs_helper_gpu()
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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

        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geometryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *X, *S);

        //            IGA_matvec_A_multiply_v1<double, 192><<<
        //            geometryNum,
        //            192>>>(
        //                    get_spline_context(),
        //                    get_assembly_context(),
        //                    *X,
        //                    *S
        //            );

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_matvec_A_Multiply Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

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

        std::vector<std::vector<std::vector<double>>> P_tripleVec;
        helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

        double Delta_new = 0;

        Delta_new = R_cps.squaredNorm();

        double Delta_0(Delta_new);

        GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

        auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

        S->reset(0.0, gpmesh::DEVICE);

        // GPU CG Step 2: First Iteration
        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geometryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *P, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        convert_spline_attribute_to_eigen_matrix(*S, S_cps);

        double Alpha = 0.0;

        for (int i = 0; i < S_cps.rows(); ++i) {
            Alpha += S_cps(i, 0) * P_cps(i, 0);
        }

        Alpha = Delta_new / Alpha;

        GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

        Disp_cps += Alpha * P_cps;
        R_cps -= Alpha * S_cps;

        double Delta_old = Delta_new;

        Delta_new = R_cps.squaredNorm();

        GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);

        double Beta = Delta_new / Delta_old;

        GPMESH_TRACE("Beta in GPU iteration-{} is {}", 1, Beta);

        P_cps = R_cps + Beta * P_cps;
        hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
        P->move(gpmesh::HOST, gpmesh::DEVICE);

        uint32_t iterNum = 1;

        while (iterNum < 100000) {
            iterNum++;

            S->reset(0.0, gpmesh::DEVICE);
            IGA_debug_matvec_A_Multiply<double, localCPointNum>
                <<<geometryNum, localCPointNum>>>(
                    get_spline_context(), get_assembly_context(), *P, *S);

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            Alpha = 0.0;

            for (int i = 0; i < S_cps.rows(); ++i) {
                Alpha += S_cps(i, 0) * P_cps(i, 0);
            }
            GPMESH_TRACE("GPU Iteration {} : ", iterNum);
            GPMESH_TRACE("Delta_old={}, P.dot(AP)={}", Delta_new, Alpha);
            Alpha = Delta_new / Alpha;

            Disp_cps += Alpha * P_cps;

            R_cps -= Alpha * S_cps;

            Delta_old = Delta_new;

            Delta_new = R_cps.squaredNorm();

            Beta = Delta_new / Delta_old;

            P_cps = R_cps + Beta * P_cps;


            GPMESH_TRACE("Delta_new={}", Delta_new);

            hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
            P->move(gpmesh::HOST, gpmesh::DEVICE);


            if (Delta_new < 1e-12 * Delta_0) {
                break;
            }
        }

        auto gpu_result = Disp_cps;
    }

    void debug_solving_assembly_rhs_helper_gpu(Eigen::MatrixXd& res)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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

        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geometryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *X, *S);

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_debug_matvec_A_Multiply Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

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

        std::vector<std::vector<std::vector<double>>> P_tripleVec;
        helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

        double Delta_new = 0;

        Delta_new = R_cps.squaredNorm();

        double Delta_0(Delta_new);

        GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

        auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

        S->reset(0.0, gpmesh::DEVICE);

        // GPU CG Step 2: First Iteration
        IGA_debug_matvec_A_Multiply<double, localCPointNum>
            <<<geometryNum, localCPointNum>>>(
                get_spline_context(), get_assembly_context(), *P, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        convert_spline_attribute_to_eigen_matrix(*S, S_cps);

        double Alpha = 0.0;

        for (int i = 0; i < S_cps.rows(); ++i) {
            Alpha += S_cps(i, 0) * P_cps(i, 0);
        }

        Alpha = Delta_new / Alpha;

        GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

        Disp_cps += Alpha * P_cps;
        R_cps -= Alpha * S_cps;

        double Delta_old = Delta_new;

        Delta_new = R_cps.squaredNorm();

        GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);

        double Beta = Delta_new / Delta_old;

        GPMESH_TRACE("Beta in GPU iteration-{} is {}", 1, Beta);

        P_cps = R_cps + Beta * P_cps;
        hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
        P->move(gpmesh::HOST, gpmesh::DEVICE);

        uint32_t iterNum = 1;

        while (iterNum < 100000) {
            iterNum++;

            S->reset(0.0, gpmesh::DEVICE);
            IGA_debug_matvec_A_Multiply<double, localCPointNum>
                <<<geometryNum, localCPointNum>>>(
                    get_spline_context(), get_assembly_context(), *P, *S);

            S->move(gpmesh::DEVICE, gpmesh::HOST);

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            Alpha = 0.0;

            for (int i = 0; i < S_cps.rows(); ++i) {
                Alpha += S_cps(i, 0) * P_cps(i, 0);
            }
            GPMESH_TRACE("GPU Iteration {} : ", iterNum);
            GPMESH_TRACE("Delta_old={}, P.dot(AP)={}", Delta_new, Alpha);
            Alpha = Delta_new / Alpha;

            Disp_cps += Alpha * P_cps;

            R_cps -= Alpha * S_cps;

            Delta_old = Delta_new;

            Delta_new = R_cps.squaredNorm();

            Beta = Delta_new / Delta_old;

            P_cps = R_cps + Beta * P_cps;


            GPMESH_TRACE("Delta_new={}", Delta_new);

            hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
            P->move(gpmesh::HOST, gpmesh::DEVICE);


            if (Delta_new < 1e-12 * Delta_0) {
                break;
            }
        }

        auto gpu_result = Disp_cps;

        res = gpu_result;
    }

    void debug_solving_assembly_rhs_helper_gpu(Eigen::MatrixXd& res,
                                               float&           time,
                                               bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
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
            <<<geometryNum, localCPointNum>>>(
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

            IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
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

    void debug_solving_assembly_rhs_helper_gpu_optimization_v0(
        Eigen::MatrixXd& res,
        float&           time,
        bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

        Eigen::MatrixXd Disp_cps;
        Disp_cps.setZero(m_num_cpoints * m_dims, 1);

        //[i][j][d]
        std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
        helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

        auto Disp_gpu_ptr =
            this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp");

        auto Disp_gpu_global_ptr =
            this->add_spline_cpoint_attribute_global<double>(Disp_tripleVec,
                                                             "disp_global");

        auto X = Disp_gpu_ptr;

        auto X_global = Disp_gpu_global_ptr;

        auto S = this->add_spline_cpoint_attribute<double>(
            "S", 3, LOCATION_ALL, gpmesh::SoA);
        S->reset(0.0, LOCATION_ALL);

        auto S_global = this->add_spline_cpoint_attribute_global<double>(
            "S_global", 3, LOCATION_ALL, gpmesh::SoA);
        S_global->reset_local_part(0.0, LOCATION_ALL);
        S_global->reset_global_part(0.0, LOCATION_ALL);

        GPUTimer timer_gpu;
        GPUTimer timer_gpu_rt;
        CPUTimer timer_cpu_rt;

        timer_gpu.start();

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);


        // S_global = A * X_global
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        // Vector
        Eigen::MatrixXd B_cps = m_sparse_system.rhs();

        // transfer m_sparse_system.rhs() to gpu
        std::vector<std::vector<std::vector<double>>> B_tripleVec;
        helper_convert_cps_to_tripleVec<double>(B_cps, B_tripleVec);
        // load information to local part
        auto B_global = this->add_spline_cpoint_attribute_global<double>(
            B_tripleVec, "B_global");

        B_global->load_global_part(B_cps);

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

        // local part test
        Eigen::MatrixXd S_global_cps_test0;
        S_global_cps_test0.setZero(m_dims * m_num_cpoints, 1);
        for (int g = 0; g < m_num_geometry; ++g) {
            for (uint16_t local_idx = 0; local_idx < localCPointNum;
                 ++local_idx) {
                uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
                for (uint32_t d = 0; d < 3; ++d) {
                    S_global_cps_test0(global_idx +
                                       d * m_sparse_system.getControlNums()) +=
                        S_global->operator()(g, local_idx, d);
                }
            }
        }

        S_global->move_local_to_global_gpu(m_attr_mapping_context);
        S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        // global part test
        Eigen::MatrixXd S_global_cps_test1;
        S_global_cps_test1.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                S_global_cps_test1(global_d_id) =
                    S_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO("difference between S_cps and S_global_cps_test0 = {}",
                    (S_cps - S_global_cps_test0).norm());
        GPMESH_INFO("difference between S_cps and S_global_cps_test1 = {}",
                    (S_cps - S_global_cps_test1).norm());

        Eigen::MatrixXd R_cps = B_cps - S_cps;
        Eigen::MatrixXd P_cps = R_cps;

        auto R_global = this->add_spline_cpoint_attribute_global<double>(
            "R_global", 3, LOCATION_ALL, gpmesh::SoA);
        R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
        auto P_global = this->add_spline_cpoint_attribute_global<double>(
            "P_global", 3, LOCATION_ALL, gpmesh::SoA);
        P_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        P_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
            B_global->m_d_global_attr,
            S_global->m_d_global_attr,
            R_global->m_d_global_attr,
            P_global->m_d_global_attr,
            cps_length);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        vector<double> b_vec(m_num_cpoints * 3);

        for (int i = 0; i < m_num_cpoints * 3; ++i) {
            b_vec[i] = B_cps.coeff(i, 0);
        }

        // global part test
        Eigen::MatrixXd R_global_cps_test0;
        R_global_cps_test0.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                R_global_cps_test0(global_d_id) =
                    R_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO("difference between R_cps and R_global_cps_test0 = {}",
                    (R_cps - R_global_cps_test0).norm());

        std::vector<std::vector<std::vector<double>>> P_tripleVec;
        helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

        double Delta_new = 0;

        Delta_new = R_cps.squaredNorm();

        // 新方法计算Delta_new
        double               Delta_new_calculated_in_gpu = 0.0;
        ReduceHandle<double> reduce_handle(*R_global);

        Delta_new_calculated_in_gpu = reduce_handle.squaredNorm(*R_global);

        double test_p_global = reduce_handle.squaredNorm(*P_global);

        assert(Delta_new - test_p_global < 1e-5);

        GPMESH_INFO("Delta_new in CPU = {}, Delta_new in GPU is {}",
                    Delta_new,
                    Delta_new_calculated_in_gpu);
        GPMESH_INFO("Delta_new diff = {}",
                    Delta_new_calculated_in_gpu - Delta_new);

        double Delta_0(Delta_new);

        if (!bSilence)
            GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

        auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

        S->reset(0.0, gpmesh::DEVICE);

        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        // GPU CG Step 2: First Iteration
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *P, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *P_global, *S_global);

        S_global->move_local_to_global_gpu(m_attr_mapping_context);
        S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        convert_spline_attribute_to_eigen_matrix(*S, S_cps);

        S_global_cps_test1.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                S_global_cps_test1(global_d_id) =
                    S_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO(
            "difference between S_cps and S_global_cps_test1 in CG Step 2 = {}",
            (S_cps - S_global_cps_test1).norm());

        double Alpha = 0.0;

        timer_cpu_rt.start();
        for (int i = 0; i < S_cps.rows(); ++i) {
            Alpha += S_cps(i, 0) * P_cps(i, 0);
        }
        timer_cpu_rt.stop();
        GPMESH_INFO("CPU Dot time={}", timer_cpu_rt.elapsed_millis());


        double sp_dot_product_in_cpu = Alpha;

        Alpha = Delta_new / Alpha;

        double sp_dot_product_in_gpu = 0.0;

        timer_gpu_rt.start();
        sp_dot_product_in_gpu = reduce_handle.dot(*S_global, *P_global);
        timer_gpu_rt.stop();

        GPMESH_INFO("GPU Reduction Dot time={}", timer_gpu_rt.elapsed_millis());


        GPMESH_INFO("sp dot in CPU={}, in GPU={}",
                    sp_dot_product_in_cpu,
                    sp_dot_product_in_gpu);
        GPMESH_INFO("sp dot Diff = {}",
                    sp_dot_product_in_cpu - sp_dot_product_in_gpu);

        double Alpha_in_gpu =
            Delta_new_calculated_in_gpu / sp_dot_product_in_gpu;

        GPMESH_INFO("Alpha in CPU={}, in GPU={}", Alpha, Alpha_in_gpu);
        GPMESH_INFO("Alpha diff = {}", Alpha - Alpha_in_gpu);

        if (!bSilence)
            GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);
        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        Disp_cps += Alpha * P_cps;
        R_cps -= Alpha * S_cps;

        double Delta_old = Delta_new;

        double Delta_old_in_gpu = Delta_new_calculated_in_gpu;

        Delta_new = R_cps.squaredNorm();

        Delta_new_calculated_in_gpu = reduce_handle.squaredNorm(*R_global);


        GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);
        GPMESH_TRACE("Delta_new_calculated_in_gpu in GPU iteration-{} is {}",
                     1,
                     Delta_new_calculated_in_gpu);

        GPMESH_TRACE("Delta_new diff={}",
                     Delta_new - Delta_new_calculated_in_gpu);

        double Beta = Delta_new / Delta_old;

        double Beta_gpu = Delta_new_calculated_in_gpu / Delta_old_in_gpu;

        GPMESH_TRACE(
            "Beta in GPU iteration-{}: cpu={}, gpu={}", 1, Beta, Beta_gpu);
        GPMESH_TRACE("Beta difference={}", Beta - Beta_gpu);

        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

        P_cps = R_cps + Beta * P_cps;
        hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
        P->move(gpmesh::HOST, gpmesh::DEVICE);

        uint32_t iterNum = 1;

        while (iterNum < 100000) {
            iterNum++;

            S->reset(0.0, gpmesh::DEVICE);
            IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
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

    void debug_solving_assembly_rhs_helper_gpu_optimization_v1(
        Eigen::MatrixXd& res,
        float&           time,
        bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

        Eigen::MatrixXd Disp_cps;
        Disp_cps.setZero(m_num_cpoints * m_dims, 1);

        //[i][j][d]
        std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
        helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

        auto Disp_gpu_ptr =
            this->add_spline_cpoint_attribute<double>(Disp_tripleVec, "disp");

        auto Disp_gpu_global_ptr =
            this->add_spline_cpoint_attribute_global<double>(Disp_tripleVec,
                                                             "disp_global");

        auto X = Disp_gpu_ptr;

        auto X_global = Disp_gpu_global_ptr;

        auto S = this->add_spline_cpoint_attribute<double>(
            "S", 3, LOCATION_ALL, gpmesh::SoA);
        S->reset(0.0, LOCATION_ALL);

        auto S_global = this->add_spline_cpoint_attribute_global<double>(
            "S_global", 3, LOCATION_ALL, gpmesh::SoA);
        S_global->reset_local_part(0.0, LOCATION_ALL);
        S_global->reset_global_part(0.0, LOCATION_ALL);

        GPUTimer timer_gpu;
        GPUTimer timer_gpu_rt;
        CPUTimer timer_cpu_rt;

        timer_gpu.start();

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);


        // S_global = A * X_global
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        // Vector
        Eigen::MatrixXd B_cps = m_sparse_system.rhs();

        // transfer m_sparse_system.rhs() to gpu
        std::vector<std::vector<std::vector<double>>> B_tripleVec;
        helper_convert_cps_to_tripleVec<double>(B_cps, B_tripleVec);
        // load information to local part
        auto B_global = this->add_spline_cpoint_attribute_global<double>(
            B_tripleVec, "B_global");

        B_global->load_global_part(B_cps);

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

        // local part test
        Eigen::MatrixXd S_global_cps_test0;
        S_global_cps_test0.setZero(m_dims * m_num_cpoints, 1);
        for (int g = 0; g < m_num_geometry; ++g) {
            for (uint16_t local_idx = 0; local_idx < localCPointNum;
                 ++local_idx) {
                uint32_t global_idx = m_h_geometry_ltog_cp[g][local_idx];
                for (uint32_t d = 0; d < 3; ++d) {
                    S_global_cps_test0(global_idx +
                                       d * m_sparse_system.getControlNums()) +=
                        S_global->operator()(g, local_idx, d);
                }
            }
        }

        S_global->move_local_to_global_gpu(m_attr_mapping_context);
        S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        // global part test
        Eigen::MatrixXd S_global_cps_test1;
        S_global_cps_test1.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                S_global_cps_test1(global_d_id) =
                    S_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO("difference between S_cps and S_global_cps_test0 = {}",
                    (S_cps - S_global_cps_test0).norm());
        GPMESH_INFO("difference between S_cps and S_global_cps_test1 = {}",
                    (S_cps - S_global_cps_test1).norm());

        Eigen::MatrixXd R_cps = B_cps - S_cps;
        Eigen::MatrixXd P_cps = R_cps;

        auto R_global = this->add_spline_cpoint_attribute_global<double>(
            "R_global", 3, LOCATION_ALL, gpmesh::SoA);
        R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
        auto P_global = this->add_spline_cpoint_attribute_global<double>(
            "P_global", 3, LOCATION_ALL, gpmesh::SoA);
        P_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        P_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
            B_global->m_d_global_attr,
            S_global->m_d_global_attr,
            R_global->m_d_global_attr,
            P_global->m_d_global_attr,
            cps_length);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        vector<double> b_vec(m_num_cpoints * 3);

        for (int i = 0; i < m_num_cpoints * 3; ++i) {
            b_vec[i] = B_cps.coeff(i, 0);
        }

        // global part test
        Eigen::MatrixXd R_global_cps_test0;
        R_global_cps_test0.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                R_global_cps_test0(global_d_id) =
                    R_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO("difference between R_cps and R_global_cps_test0 = {}",
                    (R_cps - R_global_cps_test0).norm());

        std::vector<std::vector<std::vector<double>>> P_tripleVec;
        helper_convert_cps_to_tripleVec(P_cps, P_tripleVec);

        double Delta_new = 0;

        Delta_new = R_cps.squaredNorm();

        // 新方法计算Delta_new
        double               Delta_new_calculated_in_gpu = 0.0;
        ReduceHandle<double> reduce_handle(*R_global);

        Delta_new_calculated_in_gpu = reduce_handle.squaredNorm(*R_global);

        double test_p_global = reduce_handle.squaredNorm(*P_global);

        assert(Delta_new - test_p_global < 1e-5);

        GPMESH_INFO("Delta_new in CPU = {}, Delta_new in GPU is {}",
                    Delta_new,
                    Delta_new_calculated_in_gpu);
        GPMESH_INFO("Delta_new diff = {}",
                    Delta_new_calculated_in_gpu - Delta_new);

        double Delta_0(Delta_new);

        double Delta_0_gpu = Delta_new_calculated_in_gpu;

        if (!bSilence)
            GPMESH_TRACE("Delta_0 in GPU is {}", Delta_0);

        auto P = this->add_spline_cpoint_attribute(P_tripleVec, "P");

        S->reset(0.0, gpmesh::DEVICE);

        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        // GPU CG Step 2: First Iteration
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *P, *S);

        S->move(gpmesh::DEVICE, gpmesh::HOST);

        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *P_global, *S_global);

        S_global->move_local_to_global_gpu(m_attr_mapping_context);
        S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        convert_spline_attribute_to_eigen_matrix(*S, S_cps);

        S_global_cps_test1.setZero(m_dims * m_num_cpoints, 1);
        for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
            for (int d = 0; d < 3; ++d) {
                int global_d_id = cpIdx + d * m_sparse_system.getControlNums();
                S_global_cps_test1(global_d_id) =
                    S_global->m_h_global_attr[global_d_id];
            }
        }

        GPMESH_INFO(
            "difference between S_cps and S_global_cps_test1 in CG Step 2 = {}",
            (S_cps - S_global_cps_test1).norm());

        double Alpha = 0.0;

        timer_cpu_rt.start();
        for (int i = 0; i < S_cps.rows(); ++i) {
            Alpha += S_cps(i, 0) * P_cps(i, 0);
        }
        timer_cpu_rt.stop();
        GPMESH_INFO("CPU Dot time={}", timer_cpu_rt.elapsed_millis());


        double sp_dot_product_in_cpu = Alpha;

        Alpha = Delta_new / Alpha;

        double sp_dot_product_in_gpu = 0.0;

        timer_gpu_rt.start();
        sp_dot_product_in_gpu = reduce_handle.dot(*S_global, *P_global);
        timer_gpu_rt.stop();

        GPMESH_INFO("GPU Reduction Dot time={}", timer_gpu_rt.elapsed_millis());


        GPMESH_INFO("sp dot in CPU={}, in GPU={}",
                    sp_dot_product_in_cpu,
                    sp_dot_product_in_gpu);
        GPMESH_INFO("sp dot Diff = {}",
                    sp_dot_product_in_cpu - sp_dot_product_in_gpu);

        double Alpha_in_gpu =
            Delta_new_calculated_in_gpu / sp_dot_product_in_gpu;

        GPMESH_INFO("Alpha in CPU={}, in GPU={}", Alpha, Alpha_in_gpu);
        GPMESH_INFO("Alpha diff = {}", Alpha - Alpha_in_gpu);

        if (!bSilence)
            GPMESH_TRACE("Alpha in GPU iteration-{} is {}", 1, Alpha);

        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        Disp_cps += Alpha * P_cps;
        R_cps -= Alpha * S_cps;

        double Delta_old = Delta_new;

        double Delta_old_in_gpu = Delta_new_calculated_in_gpu;

        Delta_new = R_cps.squaredNorm();

        Delta_new_calculated_in_gpu = reduce_handle.squaredNorm(*R_global);


        GPMESH_TRACE("Delta_new in GPU iteration-{} is {}", 1, Delta_new);
        GPMESH_TRACE("Delta_new_calculated_in_gpu in GPU iteration-{} is {}",
                     1,
                     Delta_new_calculated_in_gpu);

        GPMESH_TRACE("Delta_new diff={}",
                     Delta_new - Delta_new_calculated_in_gpu);

        double Beta = Delta_new / Delta_old;

        double Beta_gpu = Delta_new_calculated_in_gpu / Delta_old_in_gpu;

        GPMESH_TRACE(
            "Beta in GPU iteration-{}: cpu={}, gpu={}", 1, Beta, Beta_gpu);
        GPMESH_TRACE("Beta difference={}", Beta - Beta_gpu);

        axpy<double>
            <<<reduce_handle.get_num_block(), reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

        P_global->move_global_to_local_gpu(m_attr_mapping_context);

        P_cps = R_cps + Beta * P_cps;
        hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
        P->move(gpmesh::HOST, gpmesh::DEVICE);

        uint32_t iterNum = 1;
        uint32_t maxIter = 10000;
        while (iterNum < maxIter) {
            iterNum++;

            S->reset(0.0, gpmesh::DEVICE);
            IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
                get_spline_context(), get_assembly_context(), *P, *S);
            S->move(gpmesh::DEVICE, gpmesh::HOST);

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            for (int cpIdx = 0; cpIdx < m_num_cpoints; ++cpIdx) {
                for (int d = 0; d < 3; ++d) {
                    int global_d_id =
                        cpIdx + d * m_sparse_system.getControlNums();
                    S_global_cps_test1(global_d_id) =
                        S_global->m_h_global_attr[global_d_id];
                }
            }

            convert_spline_attribute_to_eigen_matrix(*S, S_cps);

            GPMESH_INFO(
                "difference between S_cps and S_global_cps_test1 in CG Step 2 "
                "= {}",
                (S_cps - S_global_cps_test1).norm());

            Alpha = 0.0;

            Alpha_in_gpu = 0.0;
            for (int i = 0; i < S_cps.rows(); ++i) {
                Alpha += S_cps(i, 0) * P_cps(i, 0);
            }


            Alpha_in_gpu = reduce_handle.dot(*S_global, *P_global);

            GPMESH_TRACE("GPU Iteration {} : ", iterNum);
            GPMESH_TRACE("Delta_old={}, dot(P,AP)={}", Delta_new, Alpha);
            GPMESH_TRACE("Delta_old_gpu={}, dot_gpu(P,AP)={}",
                         Delta_new_calculated_in_gpu,
                         Alpha_in_gpu);
            Alpha        = Delta_new / Alpha;
            Alpha_in_gpu = Delta_new_calculated_in_gpu / Alpha_in_gpu;
            GPMESH_TRACE("Alpha in cpu={}", Alpha);
            GPMESH_TRACE("Alpha in gpu={}", Alpha_in_gpu);

            Disp_cps += Alpha * P_cps;

            R_cps -= Alpha * S_cps;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_in_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            GPMESH_INFO("X norm cpu={}, gpu={}",
                        Disp_cps.squaredNorm(),
                        reduce_handle.squaredNorm(*X_global));

            Delta_old = Delta_new;

            Delta_old_in_gpu = Delta_new_calculated_in_gpu;

            Delta_new = R_cps.squaredNorm();
            GPMESH_TRACE("Delta_new={}", Delta_new);

            Delta_new_calculated_in_gpu = reduce_handle.squaredNorm(*R_global);
            GPMESH_TRACE("Delta_new_gpu={}", Delta_new_calculated_in_gpu);

            Beta = Delta_new / Delta_old;

            Beta_gpu = Delta_new_calculated_in_gpu / Delta_old_in_gpu;

            GPMESH_TRACE("Beta={}", Beta);
            GPMESH_TRACE("Beta_gpu={}", Beta_gpu);

            P_cps = R_cps + Beta * P_cps;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            hepler_pass_eigen_cps_to_spline_attribute_host(*P, P_cps);
            P->move(gpmesh::HOST, gpmesh::DEVICE);

            if (Delta_new < 1e-12 * Delta_0) {
                // break;
            }

            if (Delta_new_calculated_in_gpu < 1e-12 * Delta_0_gpu) {
                break;
            }

            CUDA_ERROR(cudaStreamSynchronize(0));

            GPMESH_INFO("--------{} step finish---------", iterNum);
        }
        timer_gpu.stop();

        time = timer_gpu.elapsed_millis();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            Disp_cps(d, 0) = X_global->m_h_global_attr[d];
        }

        auto gpu_result = Disp_cps;

        res = gpu_result;
    }

    void debug_solving_assembly_rhs_helper_gpu_optimization_v2(
        Eigen::MatrixXd& res,
        float&           time,
        bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

        Eigen::MatrixXd Disp_cps;
        Disp_cps.setZero(m_num_cpoints * m_dims, 1);

        //[i][j][d]
        std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
        helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

        auto Disp_gpu_global_ptr =
            this->add_spline_cpoint_attribute_global<double>(Disp_tripleVec,
                                                             "disp_global");

        auto X_global = Disp_gpu_global_ptr;

        auto S_global = this->add_spline_cpoint_attribute_global<double>(
            "S_global", 3, LOCATION_ALL, gpmesh::SoA);
        S_global->reset_local_part(0.0, LOCATION_ALL);
        S_global->reset_global_part(0.0, LOCATION_ALL);

        auto R_global = this->add_spline_cpoint_attribute_global<double>(
            "R_global", 3, LOCATION_ALL, gpmesh::SoA);
        R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
        auto P_global = this->add_spline_cpoint_attribute_global<double>(
            "P_global", 3, LOCATION_ALL, gpmesh::SoA);
        P_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        P_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);

        GPMESH_INFO("R, P Build success");

        // init_B
        Eigen::MatrixXd B_cps = m_sparse_system.rhs();

        // transfer m_sparse_system.rhs() to gpu
        std::vector<std::vector<std::vector<double>>> B_tripleVec;
        helper_convert_cps_to_tripleVec<double>(B_cps, B_tripleVec);
        // load information to local part
        auto B_global = this->add_spline_cpoint_attribute_global<double>(
            B_tripleVec, "B_global");

        B_global->load_global_part(B_cps);

        GPMESH_INFO("init B success");

        GPUTimer timer_gpu;

        timer_gpu.start();


        // S_global = A * X_global
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
            B_global->m_d_global_attr,
            S_global->m_d_global_attr,
            R_global->m_d_global_attr,
            P_global->m_d_global_attr,
            cps_length);

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        GPMESH_INFO("init P R success");
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        // 新方法计算Delta_new

        ReduceHandle<double> reduce_handle(*R_global);

        double Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

        double Delta_0_gpu = Delta_new_gpu;


        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());

        double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

        uint32_t iterNum = 0;
        uint32_t maxIter = 10000;

        while (iterNum < maxIter) {
            iterNum++;

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            if (Delta_new_gpu < 1e-12 * Delta_0_gpu) {
                break;
            }

            CUDA_ERROR(cudaStreamSynchronize(0));
        }

        timer_gpu.stop();
        time = timer_gpu.elapsed_millis();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("GPU takes {} iterations", iterNum);

        X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            Disp_cps(d, 0) = X_global->m_h_global_attr[d];
        }

        auto gpu_result = Disp_cps;

        res = gpu_result;
    }

    // same to computeUpdate, get updateVector as res
    void debug_solving_assembly_rhs_helper_gpu_optimization_v3(
        Eigen::MatrixXd& res,
        float&           time,
        bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

        Eigen::MatrixXd Disp_cps;
        Disp_cps.setZero(m_num_cpoints * m_dims, 1);

        //[i][j][d]
        std::vector<std::vector<std::vector<double>>> Disp_tripleVec;
        helper_convert_cps_to_tripleVec<double>(Disp_cps, Disp_tripleVec);

        auto Disp_gpu_global_ptr =
            this->add_spline_cpoint_attribute_global<double>(Disp_tripleVec,
                                                             "disp_global");

        auto X_global = Disp_gpu_global_ptr;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
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


        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        // CUDA_ERROR(cudaDeviceSynchronize());
        // CUDA_ERROR(cudaGetLastError());

        double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

        uint32_t iterNum = 0;
        uint32_t maxIter = 10000;

        while (iterNum < maxIter) {
            iterNum++;

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            if (Delta_new_gpu < 1e-11 * Delta_0_gpu) {
                break;
            }

            // CUDA_ERROR(cudaStreamSynchronize(0));
        }

        timer_gpu.stop();
        time = timer_gpu.elapsed_millis();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("GPU takes {} iterations", iterNum);

        X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            Disp_cps(d, 0) = X_global->m_h_global_attr[d];
        }

        auto gpu_result = Disp_cps;

        res = gpu_result;
    }

    //CG solver: solving one newton step
    void debug_frame_process_helper_v0(Eigen::MatrixXd& res,
                                       float&           time,
                                       bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_num_geometry;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
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
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
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

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
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

    void debug_frame_process_helper_v1(Eigen::MatrixXd& res,
                                       float&           time,
                                       bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

        Eigen::MatrixXd updateVector_cps;
        updateVector_cps.setZero(m_num_cpoints * m_dims, 1);

        //[i][j][d]
        std::vector<std::vector<std::vector<double>>> updateVector_tripleVec;
        helper_convert_cps_to_tripleVec<double>(updateVector_cps,
                                                updateVector_tripleVec);

        auto updateVector_gpu_global_ptr =
            this->add_spline_cpoint_attribute_global<double>(
                updateVector_tripleVec, "updateVector_global_1");

        // X_global is the initial solution, influence the efficiency of solver
        auto X_global = updateVector_gpu_global_ptr;
        X_global->reset_global_part(0.0, LOCATION_ALL);
        // GPMESH_INFO("X Build success");

        auto S_global = this->add_spline_cpoint_attribute_global<double>(
            "S_global_1", 3, LOCATION_ALL, gpmesh::SoA);
        S_global->reset_local_part(0.0, LOCATION_ALL);
        // GPMESH_INFO("S local part reset success");
        S_global->reset_global_part(0.0, LOCATION_ALL);

        // GPMESH_INFO("S Build success");

        auto R_global = this->add_spline_cpoint_attribute_global<double>(
            "R_global_1", 3, LOCATION_ALL, gpmesh::SoA);
        R_global->reset_local_part(0.0, gpmesh::LOCATION_ALL);
        R_global->reset_global_part(0.0, gpmesh::LOCATION_ALL);
        auto P_global = this->add_spline_cpoint_attribute_global<double>(
            "P_global_1", 3, LOCATION_ALL, gpmesh::SoA);
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
            B_tripleVec, "B_global_1");

        B_global->load_global_part(B_cps);

        // GPMESH_INFO("init B success");

        GPUTimer timer_gpu;

        timer_gpu.start();


        // S_global = A * X_global
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
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


        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        // CUDA_ERROR(cudaDeviceSynchronize());
        // CUDA_ERROR(cudaGetLastError());

        double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

        uint32_t iterNum = 0;
        uint32_t maxIter = 10000;

        while (iterNum < maxIter) {
            iterNum++;

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
                break;
            }

            // CUDA_ERROR(cudaStreamSynchronize(0));
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

    void debug_frame_process_helper_v2(Eigen::MatrixXd& res,
                                       float&           time,
                                       bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
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
        uint32_t maxIter = 10000;

        while (iterNum < maxIter) {

            GPMESH_INFO("GPU Iter {}", iterNum);

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            GPMESH_INFO("Ap * p = {}", Alpha_gpu);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            GPMESH_INFO("Alpha = {}", Alpha_gpu);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            GPMESH_INFO("Delta_new = {}", Delta_new_gpu);

            CUDA_ERROR(cudaStreamSynchronize(0));

            // 1e-16
            // 1e-13
            if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
                break;
            }

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;

            GPMESH_INFO("Beta = {}", Beta_gpu);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
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

    //PCG solver: solving one newton step
    void debug_frame_process_PCG_jacobiPrecondition(Eigen::MatrixXd& res,
                                       float&           time,
                                       bool             bSilence = true)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_num_geometry;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
                get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
                B_global->m_d_global_attr,
                S_global->m_d_global_attr,
                R_global->m_d_global_attr,
                P_global->m_d_global_attr,
                cps_length);

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);


         GPMESH_INFO("init P R success");
         CUDA_ERROR(cudaDeviceSynchronize());
         CUDA_ERROR(cudaGetLastError());

         //add Jacobi Precondition
         auto Preconditioner = this->add_spline_cpoint_attribute_global<double>(
                 "Preconditioner", 3, LOCATION_ALL, gpmesh::SoA);

        build_preconditioner_from_system_matrix(Preconditioner);

        {
            //test for preconditioner

        }

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
            <<<geometryNum, 192>>>(get_spline_context(),
                                   get_assembly_context(),
                                   *P_global,
                                   *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            axpy<double><<<reduce_handle.get_num_block(),
            reduce_handle.get_num_thread()>>>(
                    *X_global,
                    *P_global,
                    Alpha_gpu,
                    1.0,
                    X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
            reduce_handle.get_num_thread()>>>(
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

            axpy<double><<<reduce_handle.get_num_block(),
            reduce_handle.get_num_thread()>>>(
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
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>);

        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);

        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
    }

    template <typename T>
    void debug_system_matrix_assembly()
    {

        dim3               num_blocks_clear(m_num_geometry, 1, 1);
        dim3               num_blocks(m_num_geometry, 64, 1);
        dim3               threads_per_blocks(16, 16, 1);
        constexpr uint32_t gaussianPointNum = 64;


        IGA_clear_assembly_context<<<num_blocks_clear, threads_per_blocks>>>(
            m_assembly_context);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double,
                                                             gaussianPointNum>);

        IGA_hyperelasticity_assembly_debug_kernel_v2_debug<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
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
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
    }

    template <typename T>
    void build_preconditioner_from_system_matrix(std::shared_ptr<SplineAttributeGlobal<T>> ptr_attribute_global)
    {
        IGA_kernel_build_preconditioner_from_system_matrix<double><<<m_num_geometry, 192>>>(*ptr_attribute_global,
                                                                                 get_assembly_context());
        ptr_attribute_global->move_local_to_global_gpu(m_attr_mapping_context);

        //get reciporal value of ptr_attribute_global;

        ptr_attribute_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
    }

    template <typename T>
    void debug_stiff_matrix_assembly()
    {
        dim3               num_blocks_clear(m_num_geometry, 1, 1);
        dim3               num_blocks(m_num_geometry, 64, 1);
        dim3               threads_per_blocks(16, 16, 1);
        constexpr uint32_t gaussianPointNum = 64;


        IGA_clear_assembly_context<<<num_blocks_clear, threads_per_blocks>>>(
            m_assembly_context);

        this->m_quite               = true;
        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
    }

    template <typename T>
    void debug_update_cpoint_displacement_in_physic_context(
        SplineAttributeGlobal<T>& X)
    {
        IGA_update_physics_context_by_spline_attr_global<<<m_num_geometry,
                                                           192>>>(
            X, m_physic_context, m_attr_mapping_context);
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
    }

    void debug_solving_assembly_rhs_helper_groudtruth()
    {
        Eigen::VectorXd             u_conjgrad;
        Eigen::SparseMatrix<double> A;
        A = m_sparse_system.matrix();
        Eigen::VectorXd B;
        B.noalias() = m_sparse_system.rhs();
        u_conjgrad.resize(m_dims * m_num_cpoints);
        u_conjgrad.setZero();

        conjgrad<double>(A, B, u_conjgrad);

        Eigen::VectorXd u_cg;

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
        // cg.setMaxIterations(50);
        cg.compute(m_sparse_system.matrix());
        u_cg = cg.solve(m_sparse_system.rhs());
        if (cg.info() == Eigen::NoConvergence ||
            cg.info() == Eigen::NumericalIssue) {
            GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
            GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
        } else {
            GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
            GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
        }

        GPMESH_INFO("Difference between u_cg and u_conjgrad = {}",
                    (u_cg - u_conjgrad).norm());
        // assert((u_cg - u_conjgrad).norm() < 1e-6);
        GPMESH_INFO("conjgrad Solver Test Success!");
    }

    void debug_solving_assembly_rhs_helper_groudtruth(Eigen::MatrixXd& res)
    {
        Eigen::VectorXd             u_conjgrad;
        Eigen::SparseMatrix<double> A;
        A = m_sparse_system.matrix();
        Eigen::VectorXd B;
        B.noalias() = m_sparse_system.rhs();
        u_conjgrad.resize(m_dims * m_num_cpoints);
        u_conjgrad.setZero();

        conjgrad<double>(A, B, u_conjgrad);

        Eigen::VectorXd u_cg;

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
        // cg.setMaxIterations(50);
        cg.compute(m_sparse_system.matrix());
        u_cg = cg.solve(m_sparse_system.rhs());
        if (cg.info() == Eigen::NoConvergence ||
            cg.info() == Eigen::NumericalIssue) {
            GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
            GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
        } else {
            GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
            GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
        }

        GPMESH_INFO("Difference between u_cg and u_conjgrad = {}",
                    (u_cg - u_conjgrad).norm());
        // assert((u_cg - u_conjgrad).norm() < 1e-6);
        GPMESH_INFO("conjgrad Solver Test Success!");

        res = u_cg;
    }

    void debug_solving_assembly_rhs_helper_groudtruth(Eigen::MatrixXd& res,
                                                      float&           time,
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
        } else {
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

    void debug_solving_assembly_rhs_helper_cpu()
    {
        Eigen::SparseMatrix<double> A = m_sparse_system.matrix();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);

        Eigen::VectorXd b = m_sparse_system.rhs();

        CPUTimer timer_cpu;

        timer_cpu.start();
        Eigen::VectorXd r = b - A * x;
        timer_cpu.stop();

        GPMESH_INFO("Eigen r = b - A * x Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_cpu.elapsed_millis());

        Eigen::VectorXd p = r;

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);
        GPMESH_TRACE("delta_0 in CPU is {}", delta_0);

        Eigen::VectorXd Ap = A * p;

        GPMESH_TRACE("delta_new in CPU iteration - {} is {}", 1, delta_new);
        GPMESH_TRACE("p.dot(Ap) in CPU iteration-{} is {}", 1, p.dot(Ap));

        double alpha = delta_new / p.dot(Ap);
        GPMESH_TRACE("alpha in CPU iteration-{} is {}", 1, alpha);

        x += alpha * p;
        r -= alpha * Ap;

        double delta_old = delta_new;
        delta_new        = r.squaredNorm();

        GPMESH_TRACE("delta_new in CPU iteration{} is {}", 1, delta_new);

        double beta = delta_new / delta_old;
        GPMESH_TRACE("beta in CPU iteration-{} is {}", 1, beta);

        p = r + beta * p;

        uint32_t iterNum = 0;

        while (iterNum < 10000) {
            // CPU Iter
            Ap = A * p;

            alpha = delta_new / p.dot(Ap);
            GPMESH_TRACE("CPU Iteration {} : ", iterNum);
            GPMESH_TRACE("delta_old={}, p.dot(Ap)={}", delta_new, p.dot(Ap));

            x += alpha * p;
            r -= alpha * Ap;

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            beta = delta_new / delta_old;

            p = r + beta * p;

            GPMESH_TRACE("delta_new={}", delta_new);

            if (delta_new < 1e-12 * delta_0) {
                break;
            }
            iterNum++;
        }

        auto cpu_result = x;
    }

    void debug_solving_assembly_rhs_helper_cpu(Eigen::MatrixXd& res)
    {
        Eigen::SparseMatrix<double> A = m_sparse_system.matrix();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);

        Eigen::VectorXd b = m_sparse_system.rhs();

        CPUTimer timer_cpu;

        timer_cpu.start();
        Eigen::VectorXd r = b - A * x;
        timer_cpu.stop();

        GPMESH_INFO("Eigen r = b - A * x Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_cpu.elapsed_millis());

        Eigen::VectorXd p = r;

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);
        GPMESH_TRACE("delta_0 in CPU is {}", delta_0);

        Eigen::VectorXd Ap = A * p;

        GPMESH_TRACE("delta_new in CPU iteration - {} is {}", 1, delta_new);
        GPMESH_TRACE("p.dot(Ap) in CPU iteration-{} is {}", 1, p.dot(Ap));

        double alpha = delta_new / p.dot(Ap);
        GPMESH_TRACE("alpha in CPU iteration-{} is {}", 1, alpha);

        x += alpha * p;
        r -= alpha * Ap;

        double delta_old = delta_new;
        delta_new        = r.squaredNorm();

        GPMESH_TRACE("delta_new in CPU iteration{} is {}", 1, delta_new);

        double beta = delta_new / delta_old;
        GPMESH_TRACE("beta in CPU iteration-{} is {}", 1, beta);

        p = r + beta * p;

        uint32_t iterNum = 1;

        while (iterNum < 100000) {
            // CPU Iter
            Ap = A * p;

            alpha = delta_new / p.dot(Ap);
            GPMESH_TRACE("CPU Iteration {} : ", iterNum);
            GPMESH_TRACE("delta_old={}, p.dot(Ap)={}", delta_new, p.dot(Ap));

            x += alpha * p;
            r -= alpha * Ap;

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            beta = delta_new / delta_old;

            p = r + beta * p;

            GPMESH_TRACE("delta_new={}", delta_new);

            if (delta_new < 1e-12 * delta_0) {
                break;
            }
        }

        auto cpu_result = x;

        res = cpu_result;
    }

    void debug_solving_assembly_rhs_helper_cpu(Eigen::MatrixXd& res,
                                               float&           time,
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
        delta_new        = r.squaredNorm();

        if (!bSilence)
            GPMESH_TRACE("delta_new in CPU iteration{} is {}", 1, delta_new);

        double beta = delta_new / delta_old;
        if (!bSilence)
            GPMESH_TRACE("beta in CPU iteration-{} is {}", 1, beta);

        p = r + beta * p;

        uint32_t iterNum = 1;

        while (iterNum < 100000) {
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

    void debug_solving_assembly_rhs_helper_cpu_v1(Eigen::MatrixXd& res,
                                                  float&           time,
                                                  bool bSilence = true)
    {

        Eigen::SparseMatrix<double> A = m_sparse_system.matrix();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);

        Eigen::VectorXd b = m_sparse_system.rhs();

        Eigen::VectorXd s;

        GPUTimer timer_cpu;
        timer_cpu.start();

        s = A * x;

        Eigen::VectorXd r = b - s;

        Eigen::VectorXd p = r;

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);

        GPMESH_INFO("delta_0 in CPU is {}", delta_0);

        double alpha(0), beta(0), delta_old(0);

        uint32_t iterNum = 0;

        while (iterNum < 100000) {

            GPMESH_INFO("CPU Iter {}", iterNum);

            // CPU Iter
            s = A * p;

            alpha = p.dot(s);

            GPMESH_INFO("ap * p = {}", alpha);

            alpha = delta_new / alpha;

            GPMESH_INFO("alpha = {}", alpha);

            x += alpha * p;
            r -= alpha * s;

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            GPMESH_INFO("delta_new = {}", delta_new);

            beta = delta_new / delta_old;

            p = r + beta * p;

            GPMESH_INFO("beta = {}", beta);

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

    void debug_solving_assembly_rhs_helper_cpu_v2(Eigen::MatrixXd& res,
                                                  float&           time,
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

        res = u_conjgrad;
        time = timer_cpu.elapsed_millis();
    }

    // hyperelasticity stiffmatrix assembly verify;
    void debug_IGA_hyperE_test(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        int    numIterations    = 0;
        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
        ++numIterations;

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
        CPUTimer                     timer_cpu;
        timer_cpu.start();
        stiff.test_stiffNonlinearMatrix(displacement, true, test_localMatsVec);
        timer_cpu.stop();

        GPMESH_TRACE("stiff.test_stiffNonlinearMatrix took {} (ms)",
                     timer_cpu.elapsed_millis());

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

    // local matrix assembly, local residual assembly
    void debug_IGA_hyperE_assembly(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        int    numIterations    = 0;
        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());
        ++numIterations;

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

    // try to get full local mats and local rhs
    void debug_IGA_hyperelasticity_assembly_dirichlet_test(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);


        // stiff matrix operation
        stiffAssembler  stiff(patch, m_dirichletSides, m_neumannSides);
        Eigen::MatrixXd u;

        auto            dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);

        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        Eigen::SparseMatrix<double>  A_stiff = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;

        std::vector<int> localIndex(64, 0);

        std::iota(localIndex.begin(), localIndex.end(), 0);

        // check localMats ColMajor properties
        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);

            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        m_sparse_system.dataPush();

        auto m_p_hash = [&](const pair<uint32_t, uint32_t>& p) -> size_t {
            return p.first + p.second * m_dims * m_num_cpoints;
        };

        // assembly local mats to global mats
        unordered_map<pair<uint32_t, uint32_t>, double, decltype(m_p_hash)>
            global_mat_map(m_num_geometry * 192 * 192, m_p_hash);

        for (uint32_t geom_id = 0; geom_id < m_num_geometry; ++geom_id) {
            for (uint16_t local_i = 0; local_i < localCPointNum; ++local_i) {
                uint32_t global_i = m_h_geometry_ltog_cp[geom_id][local_i];
                for (uint16_t i_d = 0; i_d < 3; ++i_d) {
                    for (uint16_t local_j = 0; local_j < localCPointNum;
                         ++local_j) {
                        uint32_t global_j =
                            m_h_geometry_ltog_cp[geom_id][local_j];
                        for (uint16_t j_d = 0; j_d < 3; ++j_d) {
                            pair<int, int> key =
                                make_pair(global_i + i_d * m_num_cpoints,
                                          global_j + j_d * m_num_cpoints);
                            global_mat_map[key] +=
                                m_h_assembly_info[geom_id]
                                    .local_mats[(local_i + i_d * 64) +
                                                192 * (local_j + j_d * 64)];
                        }
                    }
                }
            }
        }

        std::vector<Eigen::Triplet<double>> triplets_vec;

        for (auto& pr : global_mat_map) {
            triplets_vec.emplace_back(
                pr.first.first, pr.first.second, pr.second);
        }

        Eigen::SparseMatrix<double> SparseMatrixFromGPU;
        SparseMatrixFromGPU.resize(m_dims * m_num_cpoints,
                                   m_dims * m_num_cpoints);
        SparseMatrixFromGPU.setFromTriplets(triplets_vec.begin(),
                                            triplets_vec.end());
        SparseMatrixFromGPU.makeCompressed();

        Eigen::SparseMatrix<double> SparseMatrixFromGPU2 =
            m_sparse_system.matrix();

        double diff_sparse_matrix1 =
            (GroundTruthSparseMatrix - SparseMatrixFromGPU).norm();

        double diff_sparse_matrix2 =
            (GroundTruthSparseMatrix - SparseMatrixFromGPU2).norm();

        GPMESH_INFO("Sparse Matrix Diff1 = {}", diff_sparse_matrix1);
        GPMESH_INFO("Sparse Matrix Diff2 = {}", diff_sparse_matrix2);

        // find Accuracy Problem Here
        bool bMatDirichletVerify = true;
        for (int i = 0; i < GroundTruthSparseMatrix.rows(); ++i) {
            for (int j = 0; j < GroundTruthSparseMatrix.cols(); ++j) {
                double local_diff = GroundTruthSparseMatrix.coeff(i, j) -
                                    SparseMatrixFromGPU.coeff(i, j);

                if (local_diff > 1e-12) {
                    bMatDirichletVerify = false;
                    GPMESH_INFO(
                        "Error Here : diff = {} for ({},{}) = {} and {}",
                        local_diff,
                        i,
                        j,
                        GroundTruthSparseMatrix.coeff(i, j),
                        SparseMatrixFromGPU.coeff(i, j));
                }
            }
        }
        if (bMatDirichletVerify) {
            GPMESH_INFO("Sparse Matrix with Dirichlet Boundary Verify Success");
        } else {
            GPMESH_INFO("Sparse Matrix with Dirichlet Boundary Verify Fail");
        }

        // GPMESH_INFO("Sparse Matrix diff = {}", diff_s_m);
    }


    void debug_IGA_hyperelasticity_solving_v0(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        timer_gpu.start();
        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);
        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_mass_assembly_debug_kernel_v0 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        timer_gpu.start();
        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
        timer_gpu.stop();

        GPMESH_INFO("IGA_add_dirichlet_to_mass_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_mass_context.update_host_info(m_h_mass_info, m_num_geometry);

        timer_gpu.start();
        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
        timer_gpu.stop();
        GPMESH_INFO("IGA_add_mass_matrix_to_stiff_matrix Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            mass_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_mass_info[g].mass_local_mats, 192, 192);
        }

        vector<Eigen::MatrixXd> system_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            system_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
        }

        std::vector<int> localIndex(64, 0);
        std::iota(localIndex.begin(), localIndex.end(), 0);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        auto SystemBigMatrxiFromGPU = m_sparse_system.test_get_bigmatrix();
        // here do not adding dirichlet boundary information again
        m_sparse_system.dataPush();
        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        auto velVec  = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        // displacment = solution = m_solVec
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        // mass matrix test
        bool bMassVerify = true;
        for (int g = 0; g < m_num_geometry; ++g) {
            double diff_mass_matrix =
                (mass_matrix_gpu[g] - mass_matrix_cpu[g]).norm();
            if (diff_mass_matrix > 1e-10) {
                bMassVerify = false;
                break;
            }
        }

        if (!bMassVerify)
            GPMESH_ERROR("mass matrix verify fail");


        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd              accVec      = solver.solve(stiff.getRhs());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * displacement +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);

        // this->helper_performance_test_for_matvec_A_multiply();
        // this->helper_performance_test_for_matvec_A_multiply_v0();

        Eigen::MatrixXd gpu_result;
        float           gpu_time;
        this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        Eigen::MatrixXd cpu_result;
        float           cpu_time;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);
        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time;
        this->debug_solving_assembly_rhs_helper_groudtruth(groundtruth_result,
                                                           groundtruth_time);

        double diff_cpu_res = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res = (groundtruth_result - gpu_result).norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
    }

    // try to solving Assembly * X = Rhs
    void debug_IGA_hyperelasticity_solving_v1(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        timer_gpu.start();
        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_mass_assembly_debug_kernel_v0 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        timer_gpu.start();
        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
        timer_gpu.stop();

        GPMESH_INFO("IGA_add_dirichlet_to_mass_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_mass_context.update_host_info(m_h_mass_info, m_num_geometry);

        timer_gpu.start();
        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
        timer_gpu.stop();
        GPMESH_INFO("IGA_add_mass_matrix_to_stiff_matrix Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            mass_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_mass_info[g].mass_local_mats, 192, 192);
        }

        vector<Eigen::MatrixXd> system_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            system_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
        }

        std::vector<int> localIndex(64, 0);
        std::iota(localIndex.begin(), localIndex.end(), 0);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        auto SystemBigMatrxiFromGPU = m_sparse_system.test_get_bigmatrix();
        m_sparse_system.dataPush();
        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        auto velVec  = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        // displacment = solution = m_solVec
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        // mass matrix test
        bool bMassVerify = true;
        for (int g = 0; g < m_num_geometry; ++g) {
            double diff_mass_matrix =
                (mass_matrix_gpu[g] - mass_matrix_cpu[g]).norm();
            if (diff_mass_matrix > 1e-10) {
                bMassVerify = false;
                break;
            }
        }

        if (!bMassVerify)
            GPMESH_ERROR("mass matrix verify fail");


        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd              accVec      = solver.solve(stiff.getRhs());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * displacement +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);
    }

    void debug_IGA_hyperelasticity_solving_v2(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        timer_gpu.start();
        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);
        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_mass_assembly_debug_kernel_v0 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        timer_gpu.start();
        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
        timer_gpu.stop();

        GPMESH_INFO("IGA_add_dirichlet_to_mass_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_mass_context.update_host_info(m_h_mass_info, m_num_geometry);

        timer_gpu.start();
        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
        timer_gpu.stop();
        GPMESH_INFO("IGA_add_mass_matrix_to_stiff_matrix Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            mass_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_mass_info[g].mass_local_mats, 192, 192);
        }

        vector<Eigen::MatrixXd> system_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            system_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
        }

        std::vector<int> localIndex(64, 0);
        std::iota(localIndex.begin(), localIndex.end(), 0);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        auto SystemBigMatrxiFromGPU = m_sparse_system.test_get_bigmatrix();
        // here do not adding dirichlet boundary information again
        m_sparse_system.dataPush();
        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        auto velVec  = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        // displacment = solution = m_solVec
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        // mass matrix test
        bool bMassVerify = true;
        for (int g = 0; g < m_num_geometry; ++g) {
            double diff_mass_matrix =
                (mass_matrix_gpu[g] - mass_matrix_cpu[g]).norm();
            if (diff_mass_matrix > 1e-10) {
                bMassVerify = false;
                break;
            }
        }

        if (!bMassVerify)
            GPMESH_ERROR("mass matrix verify fail");


        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd              accVec      = solver.solve(stiff.getRhs());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * displacement +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);

        // this->helper_performance_test_for_matvec_A_multiply();
        // this->helper_performance_test_for_matvec_A_multiply_v0();

        Eigen::MatrixXd gpu_result;
        float           gpu_time = 0.0;
        // this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        this->debug_solving_assembly_rhs_helper_gpu_optimization_v3(gpu_result,
                                                                    gpu_time);

        Eigen::MatrixXd cpu_result;
        float           cpu_time = 0.0;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time = 0.0;
        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        double diff_cpu_res = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res = (groundtruth_result - gpu_result).norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
    }

    void debug_IGA_hyperelasticity_frame_process_v0(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        timer_gpu.start();
        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);
        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_mass_assembly_debug_kernel_v0 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        timer_gpu.start();
        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
        timer_gpu.stop();

        GPMESH_INFO("IGA_add_dirichlet_to_mass_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_mass_context.update_host_info(m_h_mass_info, m_num_geometry);

        timer_gpu.start();
        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
        timer_gpu.stop();
        GPMESH_INFO("IGA_add_mass_matrix_to_stiff_matrix Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            mass_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_mass_info[g].mass_local_mats, 192, 192);
        }

        vector<Eigen::MatrixXd> system_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            system_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
        }

        std::vector<int> localIndex(64, 0);
        std::iota(localIndex.begin(), localIndex.end(), 0);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        auto SystemBigMatrxiFromGPU = m_sparse_system.test_get_bigmatrix();
        // here do not adding dirichlet boundary information again
        m_sparse_system.dataPush();
        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble
        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        // m_solVector == solution
        Eigen::MatrixXd m_solVector = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd velVec      = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        Eigen::MatrixXd displacement;
        stiff.constructCoefs(m_solVector, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        // mass matrix test
        bool bMassVerify = true;
        for (int g = 0; g < m_num_geometry; ++g) {
            double diff_mass_matrix =
                (mass_matrix_gpu[g] - mass_matrix_cpu[g]).norm();
            if (diff_mass_matrix > 1e-10) {
                bMassVerify = false;
                break;
            }
        }

        // cpu assembly process don`t add dirichlet to local_mat_matrix
        // so show difference
        if (!bMassVerify)
            GPMESH_ERROR("mass matrix verify fail");


        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd accVec = solver.solve(stiff.getRhs());
        GPMESH_INFO("accVec norm = {}", accVec.norm());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * dispVector +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);


        Eigen::MatrixXd gpu_result;
        float           gpu_time = 0.0;
        // this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        this->debug_frame_process_helper_v0(gpu_result, gpu_time);

        Eigen::MatrixXd cpu_result;
        float           cpu_time = 0.0;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time = 0.0;
        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);

        // second iteration

        // in cpu, use m_solVector to do assembly again
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        assert(m_solVector.rows() == groundtruth_result.rows());
        assert(m_solVector.cols() == groundtruth_result.cols());
        m_solVector += gpu_result;
        test_localMatsVec.clear();
        test_localRhsVec.clear();

        displacement.setZero();
        stiff.constructCoefs(m_solVector, displacement);

        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();
        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
            stiff.getRhs();

        timeAssemble.assemble(m_solVector, true);

        diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
        GPMESH_INFO(
            "In second iteration\n"
            "diff_mat = {}, diff_rhs = {}",
            diff_mat,
            diff_rhs);


        // in gpu, use physic_context to assembly stiff matrix
        // physic_context have been updated in "debug_frame_process_helper_v0"
        timer_gpu.start();
        debug_system_matrix_assembly<double>();
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO("debug_test_localMatsVecCPU Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());


        // test for localMats
        //        debug_test_localMatsVecCPU(test_localMatsVec);

        // test for system matrix
        //        debug_test_system_matrix_gpu(cpu_system);
        //        debug_test_system_matrix_gpu(timeAssemble.system());

        // second solving iteration

        // m_sparse_system.matrix() = timeAssemble.getMatrix();
        debug_test_system_matrix_gpu(timeAssemble.system());

        m_sparse_system.rhs() = timeAssemble.getRhs();

        Eigen::MatrixXd gpu_result_v0;
        this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
        GPMESH_INFO("second iter gpu v0 time - {} (ms)", gpu_time);

        this->debug_frame_process_helper_v1(gpu_result, gpu_time);
        GPMESH_INFO("second iter gpu v1 time - {} (ms)", gpu_time);

        this->debug_frame_process_helper_v1(gpu_result, gpu_time);
        GPMESH_INFO("second iter gpu v1 time - {} (ms)", gpu_time);

        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        GPMESH_INFO("second iter cpu groundtruth time - {} (ms)",
                    groundtruth_time);

        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

        GPMESH_INFO("second iter cpu time - {}(ms)", cpu_time);

        diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu v1 calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("gpu v0 calculation difference - {}",
                    (groundtruth_result - gpu_result_v0).norm());

        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);
    }

    // encapsulate gpu code
    void debug_IGA_hyperelasticity_frame_process_v1(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        debug_mass_matrix_assembly<double>(density);

        debug_system_matrix_assembly<double>();

        debug_push_gpu_local_mats_to_sparse_system<double>();

        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble
        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        // m_solVector == solution
        Eigen::MatrixXd m_solVector = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd velVec      = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        Eigen::MatrixXd displacement;
        stiff.constructCoefs(m_solVector, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd accVec = solver.solve(stiff.getRhs());
        GPMESH_INFO("accVec norm = {}", accVec.norm());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * dispVector +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);


        Eigen::MatrixXd gpu_result;
        float           gpu_time = 0.0;
        // this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        this->debug_frame_process_helper_v0(gpu_result, gpu_time);

        Eigen::MatrixXd cpu_result;
        float           cpu_time = 0.0;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time = 0.0;
        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);

        // second iteration

        // in cpu, use m_solVector to do assembly again
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        assert(m_solVector.rows() == groundtruth_result.rows());
        assert(m_solVector.cols() == groundtruth_result.cols());
        m_solVector += gpu_result;
        test_localMatsVec.clear();
        test_localRhsVec.clear();

        displacement.setZero();
        stiff.constructCoefs(m_solVector, displacement);

        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();
        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
            stiff.getRhs();

        timeAssemble.assemble(m_solVector, true);

        diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
        GPMESH_INFO(
            "In second iteration\n"
            "diff_mat = {}, diff_rhs = {}",
            diff_mat,
            diff_rhs);


        // in gpu, use physic_context to assembly stiff matrix
        // physic_context have been updated in "debug_frame_process_helper_v0"
        timer_gpu.start();
        debug_system_matrix_assembly<double>();
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO("debug_test_localMatsVecCPU Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());


        // test for localMats
        //        debug_test_localMatsVecCPU(test_localMatsVec);

        // test for system matrix
        //        debug_test_system_matrix_gpu(cpu_system);
        //        debug_test_system_matrix_gpu(timeAssemble.system());

        // second solving iteration

        // m_sparse_system.matrix() = timeAssemble.getMatrix();
        debug_test_system_matrix_gpu(timeAssemble.system());

        m_sparse_system.rhs() = timeAssemble.getRhs();

        Eigen::MatrixXd gpu_result_v0;
        this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
        GPMESH_INFO("second iter gpu v0 time - {} (ms)", gpu_time);

        this->debug_frame_process_helper_v1(gpu_result, gpu_time);
        GPMESH_INFO("second iter gpu v1 time - {} (ms)", gpu_time);

        this->debug_frame_process_helper_v1(gpu_result, gpu_time);
        GPMESH_INFO("second iter gpu v1 time - {} (ms)", gpu_time);

        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        GPMESH_INFO("second iter cpu groundtruth time - {} (ms)",
                    groundtruth_time);

        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

        GPMESH_INFO("second iter cpu time - {}(ms)", cpu_time);

        diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu v1 calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("gpu v0 calculation difference - {}",
                    (groundtruth_result - gpu_result_v0).norm());

        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);
    }

    // use newtonStep to do iteration
    void debug_IGA_hyperelasticity_frame_process_v2(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // gpu mass matrix assembly
        debug_mass_matrix_assembly<double>(density);

        // gpu system matrix assembly
        debug_system_matrix_assembly<double>();

        debug_push_gpu_local_mats_to_sparse_system<double>();

        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble
        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);
        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());
        massAssembler mass2(patch, density, stiff.system());

        // m_solVector == solution
        Eigen::MatrixXd m_solVector = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd velVec      = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        Eigen::MatrixXd displacement;
        stiff.constructCoefs(m_solVector, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd accVec = solver.solve(stiff.getRhs());
        GPMESH_INFO("accVec norm = {}", accVec.norm());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * dispVector +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(dispVector);
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        debug_test_system_matrix_gpu(timeAssemble.system());
        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);


        Eigen::MatrixXd gpu_result;
        float           gpu_time = 0.0;
        // this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        this->debug_frame_process_helper_v0(gpu_result, gpu_time);

        Eigen::MatrixXd cpu_result;
        float           cpu_time = 0.0;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time = 0.0;
        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);

        // second iteration

        // in cpu, use m_solVector to do assembly again
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        assert(m_solVector.rows() == groundtruth_result.rows());
        assert(m_solVector.cols() == groundtruth_result.cols());
        m_solVector += groundtruth_result;
        test_localMatsVec.clear();
        test_localRhsVec.clear();

        displacement.setZero();
        stiff.constructCoefs(m_solVector, displacement);

        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();
        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
            stiff.getRhs();

        timeAssemble.assemble(m_solVector, true);

        diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
        GPMESH_INFO(
            "In second iteration\n"
            "diff_mat = {}, diff_rhs = {}",
            diff_mat,
            diff_rhs);


        // in gpu, use physic_context to assembly stiff matrix
        // physic_context have been updated in "debug_frame_process_helper_v0"
        timer_gpu.start();
        debug_system_matrix_assembly<double>();
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO("debug_test_localMatsVecCPU Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());


        // test for localMats
        // debug_test_localMatsVecCPU(test_localMatsVec);

        // test for system matrix
        //        debug_test_system_matrix_gpu(cpu_system);
        //        debug_test_system_matrix_gpu(timeAssemble.system());

        // second solving iteration

        // m_sparse_system.matrix() = timeAssemble.getMatrix();
        // CUDA_ERROR(cudaDeviceSynchronize());
        debug_test_system_matrix_gpu(timeAssemble.system());

        m_sparse_system.rhs() = timeAssemble.getRhs();

        Eigen::MatrixXd gpu_result_v0;
        this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
        GPMESH_INFO("second iter gpu v0 time - {} (ms)", gpu_time);

        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        GPMESH_INFO("second iter cpu groundtruth time - {} (ms)",
                    groundtruth_time);

        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

        GPMESH_INFO("second iter cpu time - {}(ms)", cpu_time);

        diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu v0 calculation difference - {}",
                    (groundtruth_result - gpu_result_v0).norm());

        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);


        // third cpu solving iteration

        assert(m_solVector.rows() == groundtruth_result.rows());
        assert(m_solVector.cols() == groundtruth_result.cols());
        m_solVector += groundtruth_result;
        test_localMatsVec.clear();
        test_localRhsVec.clear();

        displacement.setZero();

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();

        stiff.constructCoefs(m_solVector, displacement);
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();
        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
            stiff.getRhs();

        timeAssemble.assemble(m_solVector, true);

        diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
        GPMESH_INFO(
            "In third iteration\n"
            "diff_mat = {}, diff_rhs = {}",
            diff_mat,
            diff_rhs);


        // in gpu, use physic_context to assembly stiff matrix
        // physic_context have been updated in "debug_frame_process_helper_v0"
        timer_gpu.start();
        debug_system_matrix_assembly<double>();
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());

        // test for localMats
        // debug_test_localMatsVecCPU(test_localMatsVec);

        // test for system matrix
        //        debug_test_system_matrix_gpu(cpu_system);
        //        debug_test_system_matrix_gpu(timeAssemble.system());

        // second solving iteration

        // m_sparse_system.matrix() = timeAssemble.getMatrix();
        debug_test_system_matrix_gpu(timeAssemble.system());

        // comment this
        // sovling process is ok
        // error in assembly process
        m_sparse_system.matrix() = timeAssemble.getMatrix();
        m_sparse_system.rhs()    = timeAssemble.getRhs();

        this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
        GPMESH_INFO("third iter gpu v0 time - {} (ms)", gpu_time);

        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        GPMESH_INFO("third iter cpu groundtruth time - {} (ms)",
                    groundtruth_time);

        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

        GPMESH_INFO("third iter cpu time - {}(ms)", cpu_time);

        diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        // GPMESH_INFO("gpu v1 calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("gpu v0 calculation difference - {}",
                    (groundtruth_result - gpu_result_v0).norm());

        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);


        // fourth cpu solving iteration
        assert(m_solVector.rows() == groundtruth_result.rows());
        assert(m_solVector.cols() == groundtruth_result.cols());
        m_solVector += groundtruth_result;
        test_localMatsVec.clear();
        test_localRhsVec.clear();

        displacement.setZero();

        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();

        stiff.constructCoefs(m_solVector, displacement);
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();
        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
            stiff.getRhs();

        timeAssemble.assemble(m_solVector, true);

        diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
        GPMESH_INFO(
            "In third iteration\n"
            "diff_mat = {}, diff_rhs = {}",
            diff_mat,
            diff_rhs);


        // in gpu, use physic_context to assembly stiff matrix
        // physic_context have been updated in "debug_frame_process_helper_v0"
        timer_gpu.start();
        debug_system_matrix_assembly<double>();
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());

        // test for localMats
        // debug_test_localMatsVecCPU(test_localMatsVec);

        // test for system matrix
        //        debug_test_system_matrix_gpu(cpu_system);
        //        debug_test_system_matrix_gpu(timeAssemble.system());

        // second solving iteration

        // m_sparse_system.matrix() = timeAssemble.getMatrix();
        debug_test_system_matrix_gpu(timeAssemble.system());

        // comment this
        // sovling process is ok
        // error in assembly process
        m_sparse_system.matrix() = timeAssemble.getMatrix();
        m_sparse_system.rhs()    = timeAssemble.getRhs();

        this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
        GPMESH_INFO("fourth iter gpu v0 time - {} (ms)", gpu_time);

        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        GPMESH_INFO("fourth iter cpu groundtruth time - {} (ms)",
                    groundtruth_time);

        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

        GPMESH_INFO("fourth iter cpu time - {}(ms)", cpu_time);

        diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        // GPMESH_INFO("gpu v1 calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("gpu v0 calculation difference - {}",
                    (groundtruth_result - gpu_result_v0).norm());

        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        int iter = 0;

        while (iter < 10) {
            iter++;
            m_solVector += groundtruth_result;
            test_localMatsVec.clear();
            test_localRhsVec.clear();

            displacement.setZero();

            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();

            stiff.constructCoefs(m_solVector, displacement);
            stiff.test_stiffNonlinearMatrix(
                displacement, true, test_localMatsVec, test_localRhsVec);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            timeAssemble.assemble(m_solVector, true);

            diff_mat = (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
            diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();
            GPMESH_INFO(
                "In num-{} iteration\n"
                "diff_mat = {}, diff_rhs = {}",
                iter,
                diff_mat,
                diff_rhs);


            // in gpu, use physic_context to assembly stiff matrix
            // physic_context have been updated in
            // "debug_frame_process_helper_v0"
            timer_gpu.start();
            debug_system_matrix_assembly<double>();
            timer_gpu.stop();
            CUDA_ERROR(cudaDeviceSynchronize());

            // test for localMats
            // debug_test_localMatsVecCPU(test_localMatsVec);

            // test for system matrix
            //        debug_test_system_matrix_gpu(cpu_system);
            //        debug_test_system_matrix_gpu(timeAssemble.system());

            // second solving iteration

            // m_sparse_system.matrix() = timeAssemble.getMatrix();
            debug_test_system_matrix_gpu(timeAssemble.system());

            // comment this
            // sovling process is ok
            // error in assembly process
            m_sparse_system.matrix() = timeAssemble.getMatrix();
            m_sparse_system.rhs()    = timeAssemble.getRhs();

            this->debug_frame_process_helper_v0(gpu_result_v0, gpu_time);
            GPMESH_INFO("{} iter gpu v0 time - {} (ms)", iter, gpu_time);

            this->debug_solving_assembly_rhs_helper_groudtruth(
                groundtruth_result, groundtruth_time, false);

            GPMESH_INFO("{} iter cpu groundtruth time - {} (ms)",
                        iter,
                        groundtruth_time);

            this->debug_solving_assembly_rhs_helper_cpu(
                cpu_result, cpu_time, true);

            GPMESH_INFO("{} iter cpu time - {}(ms)", iter, cpu_time);

            diff_cpu_res         = (groundtruth_result - cpu_result).norm();
            diff_gpu_res         = (groundtruth_result - gpu_result).norm();
            groundtruth_norm_res = groundtruth_result.norm();

            GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
            // GPMESH_INFO("gpu v1 calculation difference - {}", diff_gpu_res);
            GPMESH_INFO("gpu v0 calculation difference - {}",
                        (groundtruth_result - gpu_result_v0).norm());

            GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);
        }
    }

    // refactor makeTimeStep into frame process
    void debug_IGA_hyperelasticity_frame_process_v3(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 100;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();


        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        while (newton_step < newton_limit) {

            // stiff.assemble(m_solVector)
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            displacement.setZero();
            stiff.constructCoefs(m_solVector, displacement);
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();

            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVector +
                     alpha3(Beta) * accVector) +
                stiff.getRhs();

            Eigen::VectorXd                                       cpu_result;
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            cpu_result = cg.solve(cpu_system.rhs());

            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            }

            double updateNorm   = cpu_result.norm();
            double residualNorm = cpu_system.rhs().norm();

            if (newton_step == 0) {
                initUpdateNorm   = updateNorm;
                initResidualNorm = residualNorm;
            }

            m_solVector += cpu_result;

            bool bQuit1Cur = false;

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm / initResidualNorm < relTol ||
                updateNorm / initUpdateNorm < relTol) {
                bQuit1Cur = true;
            }

            if (!bQuit1 && bQuit1Cur) {
                GPMESH_INFO("refactor solver should quit in newton step {}",
                            newton_step);
            }


            he_solver.computeUpdate();

            bool bQuit2Cur = false;

            if (he_solver.m_residualNorm < absTol ||
                he_solver.m_updateNorm < absTol ||
                he_solver.m_residualNorm / he_solver.m_initResidualNorm <
                    relTol ||
                he_solver.m_updateNorm / he_solver.m_initUpdateNorm < relTol) {
                bQuit2Cur = true;
            }

            if (!bQuit2 && bQuit2Cur) {
                GPMESH_INFO("he_solver should quit in newton step {}",
                            newton_step);
            }

            bQuit1 = bQuit1 || bQuit1Cur;
            bQuit2 = bQuit2 || bQuit2Cur;
            GPMESH_INFO("In newton step {}", newton_step);
            GPMESH_INFO("Difference between solVector = {}",
                        (m_solVector - he_solver.solution()).norm());

            if (bQuit1 && bQuit2) {
                break;
            }
            newton_step++;
        }

        GPMESH_INFO("Newton Method finish in step {}", newton_step);
    }

    // compare one frame cpu, gpu calculation
    void debug_IGA_hyperelasticity_frame_process_v4(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 10;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();


        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        CUDA_ERROR(cudaDeviceSynchronize());
        while (newton_step < newton_limit) {
            GPMESH_INFO("newton {} step start", newton_step);
            // stiff.assemble(m_solVector)
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            displacement.setZero();
            stiff.constructCoefs(m_solVector, displacement);
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();

            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVector +
                     alpha3(Beta) * accVector) +
                stiff.getRhs();

            Eigen::VectorXd groundtruth_result;
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            groundtruth_result = cg.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            }

            double updateNorm   = groundtruth_result.norm();
            double residualNorm = cpu_system.rhs().norm();

            if (newton_step == 0) {
                initUpdateNorm   = updateNorm;
                initResidualNorm = residualNorm;
            }

            m_solVector += groundtruth_result;

            bool bQuit1Cur = false;

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm / initResidualNorm < relTol ||
                updateNorm / initUpdateNorm < relTol) {
                bQuit1Cur = true;
            }

            if (!bQuit1 && bQuit1Cur) {
                GPMESH_INFO("refactor solver should quit in newton step {}",
                            newton_step);
            }

            he_solver.computeUpdate();

            bool bQuit2Cur = false;

            if (he_solver.m_residualNorm < absTol ||
                he_solver.m_updateNorm < absTol ||
                he_solver.m_residualNorm / he_solver.m_initResidualNorm <
                    relTol ||
                he_solver.m_updateNorm / he_solver.m_initUpdateNorm < relTol) {
                bQuit2Cur = true;
            }

            if (!bQuit2 && bQuit2Cur) {
                GPMESH_INFO("he_solver should quit in newton step {}",
                            newton_step);
            }

            bQuit1 = bQuit1 || bQuit1Cur;
            bQuit2 = bQuit2 || bQuit2Cur;
            // GPMESH_INFO("In newton step {}", newton_step);
            // GPMESH_INFO("solVector Difference between refactor and he_solver
            // = {}", (m_solVector - he_solver.solution()).norm());

            debug_system_matrix_assembly<double>();

            m_sparse_system.matrix() = cpu_system.matrix();
            m_sparse_system.rhs()    = cpu_system.rhs();

            // Todo: fix these problems
            // 1.
            // comment out debug_test_system_matrix_gpu cause gpu iteration
            // failure CUDA_ERROR(cudaDeviceSynchronize());

            // 2.
            // after 2-3 iteration,
            // system sparse matrix show big difference
            // debug_test_system_matrix_gpu(cpu_system);
            m_sparse_system.matrix() = cpu_system.matrix();

            Eigen::MatrixXd gpu_result;

            float gpu_time;

            debug_frame_process_helper_v0(gpu_result, gpu_time);
            CUDA_ERROR(cudaDeviceSynchronize());

            GPMESH_INFO(
                "updateVector Difference between groundtruth and gpu = {}",
                (groundtruth_result - gpu_result).norm());
            GPMESH_INFO("updateVector gpu norm = {}", gpu_result.norm());
            Eigen::MatrixXd cpu_result;
            float           cpu_time;
            debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

            GPMESH_INFO(
                "updateVector Difference between groundtruth and cpu = {}",
                (groundtruth_result - cpu_result).norm());
            GPMESH_INFO("updateVector cpu norm = {}", cpu_result.norm());

            // check physic context m_d_cpoints_displacement
            Eigen::MatrixXd gpu_disp = debug_get_gpu_displacement();

            GPMESH_INFO("solVector difference between cpu and gpu - {}",
                        (m_solVector - gpu_disp).norm());
            GPMESH_INFO("m_solVector norm is {}", m_solVector.norm());
            GPMESH_INFO("gpu_disp norm is {}", gpu_disp.norm());
            newton_step++;

            if (bQuit1 && bQuit2) {
                break;
            }
        }

        GPMESH_INFO("Newton Method finish in step {}", --newton_step);

        // first frame output
        string output = STRINGIFY(OUTPUT_DIR) "IGADataOutput/Konglong/";

        // output m_solVector (GroundTruth), and gpu_disp
        Eigen::MatrixXd u = m_solVector;
        stiff.constructCoefs(u, displacement);

        int begin;

        begin = 0;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement.block(begin, 0, num, 3));
            begin += num;
        }
        stiff.WriteParaview(patch, output + "kl_CPU_", u, 1000, false);
        Eigen::MatrixXd gpu_disp = debug_get_gpu_displacement();
        stiff.WriteParaview(patch, output + "kl_GPU_", gpu_disp, 1000, false);
        Eigen::MatrixXd diff_disp = u - gpu_disp;
        stiff.WriteParaview(patch, output + "kl_Diff_", diff_disp, 1000, false);

        begin = 0;

        displacement *= -1;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement.block(begin, 0, num, 3));
            begin += num;
        }
    }

    // compare one frame cpu, gpu calculation
    // check quit point
    void debug_IGA_hyperelasticity_frame_process_v5(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 10;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();


        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuitCPU1 = false;
        bool bQuitCPU2 = false;
        bool bQuitGPU  = false;

        double initUpdateNorm;
        double initResidualNorm;

        double initUpdateNorm_gpu, initResidualNorm_gpu;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        CUDA_ERROR(cudaDeviceSynchronize());
        while (newton_step < newton_limit) {
            GPMESH_INFO("newton {} step start", newton_step);
            // stiff.assemble(m_solVector)
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            displacement.setZero();
            stiff.constructCoefs(m_solVector, displacement);
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();

            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVector +
                     alpha3(Beta) * accVector) +
                stiff.getRhs();

            Eigen::VectorXd groundtruth_result;
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            groundtruth_result = cg.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            }

            double updateNorm   = groundtruth_result.norm();
            double residualNorm = cpu_system.rhs().norm();

            if (newton_step == 0) {
                initUpdateNorm   = updateNorm;
                initResidualNorm = residualNorm;
            }

            m_solVector += groundtruth_result;

            bool bQuitCPU1_iter = false;

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm / initResidualNorm < relTol ||
                updateNorm / initUpdateNorm < relTol) {
                bQuitCPU1_iter = true;
            }

            if (!bQuitCPU1 && bQuitCPU1_iter) {
                GPMESH_INFO("refactor solver should quit in newton step {}",
                            newton_step);
            }

            he_solver.computeUpdate();

            bool bQuitCPU2_iter = false;

            if (he_solver.m_residualNorm < absTol ||
                he_solver.m_updateNorm < absTol ||
                he_solver.m_residualNorm / he_solver.m_initResidualNorm <
                    relTol ||
                he_solver.m_updateNorm / he_solver.m_initUpdateNorm < relTol) {
                bQuitCPU2_iter = true;
            }

            if (!bQuitCPU2 && bQuitCPU2_iter) {
                GPMESH_INFO("he_solver should quit in newton step {}",
                            newton_step);
            }

            bQuitCPU1 = bQuitCPU1 || bQuitCPU1_iter;
            bQuitCPU2 = bQuitCPU2 || bQuitCPU2_iter;
            // GPMESH_INFO("In newton step {}", newton_step);
            // GPMESH_INFO("solVector Difference between refactor and he_solver
            // = {}", (m_solVector - he_solver.solution()).norm());

            debug_system_matrix_assembly<double>();

            m_sparse_system.matrix() = cpu_system.matrix();
            m_sparse_system.rhs()    = cpu_system.rhs();

            // Todo: fix these problems
            // 1.
            // comment out debug_test_system_matrix_gpu cause gpu iteration
            // failure CUDA_ERROR(cudaDeviceSynchronize());

            // 2.
            // after 2-3 iteration,
            // system sparse matrix show big difference
            // debug_test_system_matrix_gpu(cpu_system);
            m_sparse_system.matrix() = cpu_system.matrix();

            Eigen::MatrixXd gpu_result;

            float gpu_time;

            debug_frame_process_helper_v0(gpu_result, gpu_time);
            CUDA_ERROR(cudaDeviceSynchronize());

            GPMESH_INFO(
                "updateVector Difference between groundtruth and gpu = {}",
                (groundtruth_result - gpu_result).norm());
            GPMESH_INFO("updateVector gpu norm = {}", gpu_result.norm());
            Eigen::MatrixXd cpu_result;
            float           cpu_time;
            debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time, true);

            GPMESH_INFO(
                "updateVector Difference between groundtruth and cpu = {}",
                (groundtruth_result - cpu_result).norm());
            GPMESH_INFO("updateVector cpu norm = {}", cpu_result.norm());

            double updateNorm_gpu   = gpu_result.norm();
            double residualNorm_gpu = m_sparse_system.rhs().norm();

            if (newton_step == 0) {
                initUpdateNorm_gpu   = updateNorm_gpu;
                initResidualNorm_gpu = residualNorm_gpu;
            }

            bool bQuitGPU_iter = false;

            if (residualNorm_gpu < absTol || updateNorm_gpu < absTol ||
                residualNorm_gpu / initResidualNorm_gpu < relTol ||
                updateNorm_gpu / initUpdateNorm_gpu < relTol) {
                bQuitGPU_iter = true;
            }

            if (!bQuitGPU && bQuitGPU_iter) {
                GPMESH_INFO("GPU Solver should quit in newton step {}",
                            newton_step);
            }

            bQuitGPU = bQuitGPU || bQuitGPU_iter;

            // check physic context m_d_cpoints_displacement
            Eigen::MatrixXd gpu_disp = debug_get_gpu_displacement();

            GPMESH_INFO("solVector difference between cpu and gpu - {}",
                        (m_solVector - gpu_disp).norm());
            GPMESH_INFO("m_solVector norm is {}", m_solVector.norm());
            GPMESH_INFO("gpu_disp norm is {}", gpu_disp.norm());
            newton_step++;

            if (bQuitCPU1 && bQuitCPU2 && bQuitGPU) {
                break;
            }
        }

        GPMESH_INFO("Newton Method finish in step {}", --newton_step);

        // first frame output

        string product_name = "konglong1";
        string output = STRINGIFY(OUTPUT_DIR) "IGADataOutput/";

        output += product_name + "/";

        {
            auto        end      = std::chrono::system_clock::now();
            std::time_t end_time = std::chrono::system_clock::to_time_t(end);
            std::ostringstream time_oss;
            time_oss << std::put_time(std::localtime(&end_time),
                                      "%d-%m-%Y-%H-%M-%S");
            output += time_oss.str() + "/";
        }

        filesystem::path p(output);
        filesystem::create_directories(p);

        // output m_solVector (GroundTruth), and gpu_disp
        Eigen::MatrixXd u = m_solVector;
        stiff.constructCoefs(u, displacement);

        int begin;

        begin = 0;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement.block(begin, 0, num, 3));
            begin += num;
        }
        stiff.WriteParaview(patch, output + product_name + "_CPU_", u, 1000, false);
        Eigen::MatrixXd gpu_disp = debug_get_gpu_displacement();
        stiff.WriteParaview(patch, output + product_name + "_GPU_", gpu_disp, 1000, false);
        Eigen::MatrixXd diff_disp = u - gpu_disp;
        stiff.WriteParaview(patch, output + product_name + "_Diff_", diff_disp, 1000, false);

        begin = 0;

        displacement *= -1;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement.block(begin, 0, num, 3));
            begin += num;
        }
    }

    // one newton step performance test
    void debug_IGA_he_cg_performance_test_v0(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        GPMESH_INFO("Acc Solve Success!");

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();

        timeAssemble.assemble(timeAssemble.displacementVector());

        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

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

        m_sparse_system.matrix() = timeAssemble.getMatrix();
        m_sparse_system.rhs()    = timeAssemble.getRhs();

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

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
    }

    // one newton step performance test
    // memory optimization
    void debug_IGA_he_cg_performance_test_v1(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        GPMESH_INFO("Acc Solve Success!");

        uint32_t newton_step  = 0;
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
        m_sparse_system.rhs()    = cpu_system.rhs();

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

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
    }

    // one newton step performance test
    void debug_IGA_he_system_matrix_assembly_test_v0(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();

        timeAssemble.assemble(timeAssemble.displacementVector());

        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        // stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
            stiff.getRhs();

        GPMESH_INFO("System Matrix Difference is {}",
                    (cpu_system.matrix() - timeAssemble.getMatrix()).norm());

        m_sparse_system.matrix() = timeAssemble.getMatrix();
        m_sparse_system.rhs()    = timeAssemble.getRhs();

        debug_system_matrix_assembly<double>();

        // debug_test_system_matrix_gpu(cpu_system);

        Eigen::MatrixXd gpu_result;
        float           gpu_time = 0.0;
        // this->debug_solving_assembly_rhs_helper_gpu(gpu_result, gpu_time);
        cudaDeviceSynchronize();
        cudaStreamSynchronize(0);
        this->debug_frame_process_helper_v0(gpu_result, gpu_time);

        Eigen::MatrixXd cpu_result;
        float           cpu_time = 0.0;
        this->debug_solving_assembly_rhs_helper_cpu(cpu_result, cpu_time);

        Eigen::MatrixXd groundtruth_result;
        float           groundtruth_time = 0.0;
        this->debug_solving_assembly_rhs_helper_groudtruth(
            groundtruth_result, groundtruth_time, false);

        double diff_cpu_res         = (groundtruth_result - cpu_result).norm();
        double diff_gpu_res         = (groundtruth_result - gpu_result).norm();
        double groundtruth_norm_res = groundtruth_result.norm();

        GPMESH_INFO("cpu calculation difference - {}", diff_cpu_res);
        GPMESH_INFO("gpu calculation difference - {}", diff_gpu_res);
        GPMESH_INFO("ground truth norm - {}", groundtruth_norm_res);

        GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        GPMESH_INFO("groundtruth time - {} (ms)", groundtruth_time);
    }

    void debug_IGA_solving_process_test_v0(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();

        timeAssemble.assemble(timeAssemble.displacementVector());

        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        std::this_thread::sleep_for(std::chrono::milliseconds(1000));

        while (newton_step < newton_limit) {
            GPMESH_INFO("newton {} step start", newton_step);

            stiff.assemble(m_solVector);

            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVector +
                     alpha3(Beta) * accVector) +
                stiff.getRhs();

            debug_system_matrix_assembly<double>();

            GPMESH_INFO("System Matrix Assembly in GPU success");

            debug_test_system_matrix_gpu(cpu_system);

            // m_assembly_context.update_host_info(m_h_assembly_info,
            // m_num_geometry); debug_test_system_matrix_gpu(cpu_system);

            m_sparse_system.matrix() = cpu_system.matrix();
            m_sparse_system.rhs()    = cpu_system.rhs();

            // cg steps

            constexpr uint32_t gaussianPointNum = 64;
            const int          localCPointNum   = 64;
            const double       MAXNUM_DIRICHLET = 1e10;
            size_t             geometryNum      = m_geometry_size;

            Eigen::MatrixXd updateVector_cps;
            updateVector_cps.setZero(m_num_cpoints * m_dims, 1);

            std::vector<std::vector<std::vector<double>>>
                updateVector_tripleVec;
            helper_convert_cps_to_tripleVec<double>(updateVector_cps,
                                                    updateVector_tripleVec);

            auto updateVector_gpu_global_ptr =
                this->add_spline_cpoint_attribute_global<double>(
                    updateVector_tripleVec, "updateVector_global");

            // X_global is the initial solution, influence the efficiency of
            // solver
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
            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *X_global,
                                       *S_global);

            S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
            S_global->move_local_to_global_gpu(m_attr_mapping_context);

            // GPMESH_INFO("S calculate success");

            // init PR
            size_t cps_length          = m_num_cpoints * 3;
            size_t init_PR_threads_num = 256;
            size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

            // R_global = B_global - S_global
            // P_global = R_global
            init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
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


            // move P_global global part to local part
            P_global->move_global_to_local_gpu(m_attr_mapping_context);
            S_global->reset_local_part(0.0, gpmesh::DEVICE);

            // CUDA_ERROR(cudaDeviceSynchronize());
            // CUDA_ERROR(cudaGetLastError());

            double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

            uint32_t iterNum = 0;
            uint32_t maxIter = 10000;


            while (iterNum <= maxIter) {

                GPMESH_INFO("cg {} step start", iterNum);

                S_global->reset_local_part(0.0, gpmesh::DEVICE);
                S_global->reset_global_part(0.0, gpmesh::DEVICE);

                IGA_matvec_A_multiply_v1<double, 192>
                    <<<geometryNum, 192>>>(get_spline_context(),
                                           get_assembly_context(),
                                           *P_global,
                                           *S_global);

                S_global->move_local_to_global_gpu(m_attr_mapping_context);
                S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

                Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

                GPMESH_TRACE("AP * P = {}", Alpha_gpu);

                Alpha_gpu = Delta_new_gpu / Alpha_gpu;

                GPMESH_INFO("Alpha = {}", Alpha_gpu);


                axpy<double><<<reduce_handle.get_num_block(),
                               reduce_handle.get_num_thread()>>>(
                    *X_global,
                    *P_global,
                    Alpha_gpu,
                    1.0,
                    X_global->m_num_cpoints * 3);

                axpy<double><<<reduce_handle.get_num_block(),
                               reduce_handle.get_num_thread()>>>(
                    *R_global,
                    *S_global,
                    -Alpha_gpu,
                    1.0,
                    X_global->m_num_cpoints * 3);

                Delta_old_gpu = Delta_new_gpu;

                Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

                GPMESH_INFO("Del ta = {}", Delta_new_gpu);

                CUDA_ERROR(cudaStreamSynchronize(0));

                // 1e-16
                // 1e-13
                if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
                    break;
                }

                Beta_gpu = Delta_new_gpu / Delta_old_gpu;

                GPMESH_INFO("Beta = {}", Beta_gpu);

                axpy<double><<<reduce_handle.get_num_block(),
                               reduce_handle.get_num_thread()>>>(
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


            debug_update_cpoint_displacement_in_physic_context(*X_global);

            debug_solving_assembly_rhs_helper_cpu();

            ++newton_step;
        }
    }

    void debug_IGA_solving_process_test_v1(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();


        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        CUDA_ERROR(cudaDeviceSynchronize());
        while (newton_step < newton_limit) {
            GPMESH_INFO("newton {} step start", newton_step);
            // stiff.assemble(m_solVector)
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            displacement.setZero();
            stiff.constructCoefs(m_solVector, displacement);
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();

            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                     alpha2(Beta, tStep) * velVector +
                     alpha3(Beta) * accVector) +
                stiff.getRhs();

            Eigen::VectorXd groundtruth_result;
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            groundtruth_result = cg.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE(
                    "Iteration {}, Error {}", cg.iterations(), cg.error());
            }

            double updateNorm   = groundtruth_result.norm();
            double residualNorm = cpu_system.rhs().norm();

            if (newton_step == 0) {
                initUpdateNorm   = updateNorm;
                initResidualNorm = residualNorm;
            }

            m_solVector += groundtruth_result;

            bool bQuit1Cur = false;

            if (residualNorm < absTol || updateNorm < absTol ||
                residualNorm / initResidualNorm < relTol ||
                updateNorm / initUpdateNorm < relTol) {
                bQuit1Cur = true;
            }

            if (!bQuit1 && bQuit1Cur) {
                GPMESH_INFO("refactor solver should quit in newton step {}",
                            newton_step);
            }

            he_solver.computeUpdate();

            bool bQuit2Cur = false;

            if (he_solver.m_residualNorm < absTol ||
                he_solver.m_updateNorm < absTol ||
                he_solver.m_residualNorm / he_solver.m_initResidualNorm <
                    relTol ||
                he_solver.m_updateNorm / he_solver.m_initUpdateNorm < relTol) {
                bQuit2Cur = true;
            }

            if (!bQuit2 && bQuit2Cur) {
                GPMESH_INFO("he_solver should quit in newton step {}",
                            newton_step);
            }

            bQuit1 = bQuit1 || bQuit1Cur;
            bQuit2 = bQuit2 || bQuit2Cur;
            // GPMESH_INFO("In newton step {}", newton_step);
            // GPMESH_INFO("solVector Difference between refactor and he_solver
            // = {}", (m_solVector - he_solver.solution()).norm());

            m_sparse_system.matrix() = cpu_system.matrix();
            m_sparse_system.rhs()    = cpu_system.rhs();

            debug_system_matrix_assembly<double>();

            // debug_test_system_matrix_gpu(cpu_system);
            // m_assembly_context.update_host_info(m_h_assembly_info,
            // m_num_geometry); cudaDeviceSynchronize();


            Eigen::MatrixXd gpu_result;

            float gpu_time;

            debug_frame_process_helper_v2(gpu_result, gpu_time);
            //            CUDA_ERROR(cudaDeviceSynchronize());
            //
            //            GPMESH_INFO("updateVector Difference between
            //            groundtruth and gpu = {}", (groundtruth_result -
            //            gpu_result).norm()); GPMESH_INFO("updateVector gpu
            //            norm = {}", gpu_result.norm()); Eigen::MatrixXd
            //            cpu_result; float cpu_time;
            //            debug_solving_assembly_rhs_helper_cpu_v1(cpu_result,
            //            cpu_time, true);
            //
            //            GPMESH_INFO("updateVector Difference between
            //            groundtruth and cpu = {}", (groundtruth_result -
            //            cpu_result).norm()); GPMESH_INFO("updateVector cpu
            //            norm = {}", cpu_result.norm());
            //
            //            //check physic context m_d_cpoints_displacement
            //            Eigen::MatrixXd gpu_disp =
            //            debug_get_gpu_displacement();
            //
            //            GPMESH_INFO("solVector difference between cpu and gpu
            //            - {}", (m_solVector - gpu_disp).norm());
            //            GPMESH_INFO("m_solVector norm is {}",
            //            m_solVector.norm()); GPMESH_INFO("gpu_disp norm is
            //            {}", gpu_disp.norm());
            newton_step++;

            if (bQuit1 && bQuit2) {
                break;
            }
        }

        GPMESH_INFO("Newton Method finish in step {}", --newton_step);
        //
        //        //first frame output
        //        string output = STRINGIFY(OUTPUT_DIR)
        //        "IGADataOutput/Konglong/";
        //
        //        //output m_solVector (GroundTruth), and gpu_disp
        //        Eigen::MatrixXd u = m_solVector;
        //        stiff.constructCoefs(u, displacement);
        //
        //        int begin;
        //
        //        begin = 0;
        //        for (size_t i = 0; i < patch.size(); ++i) {
        //            int num = patch[i].getSizeControlPoints();
        //            patch[i].updateCoefs(displacement.block(begin, 0, num,
        //            3)); begin += num;
        //        }
        //        stiff.WriteParaview(patch, output + "kl_CPU_", u, 1000,
        //        false); Eigen::MatrixXd gpu_disp =
        //        debug_get_gpu_displacement(); stiff.WriteParaview(patch,
        //        output + "kl_GPU_", gpu_disp, 1000, false); Eigen::MatrixXd
        //        diff_disp = u - gpu_disp; stiff.WriteParaview(patch, output +
        //        "kl_Diff_", diff_disp, 1000, false);
        //
        //        begin = 0;
        //
        //        displacement *= -1;
        //        for (size_t i = 0; i < patch.size(); ++i) {
        //            int num = patch[i].getSizeControlPoints();
        //            patch[i].updateCoefs(displacement.block(begin, 0, num,
        //            3)); begin += num;
        //        }
    }

    // split cpu and gpu solving, try to find bugs
    // three test
    // might exceed memory limit
    void debug_IGA_solving_process_test_v2(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
        accVector     = acc_solver.solve(stiff.getRhs());
        pre_accVector = accVector;

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        SparseSystem cpu_system = stiff.system();

        stiffAssembler stiff2(patch, m_dirichletSides, m_neumannSides);
        massAssembler  mass2(patch, density, stiff.system());

        timeIntegrator timeAssemble(stiff2, mass2);
        timeAssemble.setInitialDisplacement(m_solVector);
        timeAssemble.setInitialVelocity(velVector);
        timeAssemble.initialize();


        double absTol = 1e-12;
        double relTol = 1e-9;

        Hyperelastic he_solver(&timeAssemble,
                               timeAssemble.displacementVector());

        bool bQuit1 = false;
        bool bQuit2 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        // stiff.assemble(m_solVector)
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        displacement.setZero();
        stiff.constructCoefs(m_solVector, displacement);
        stiff.stiffNonlinearMatrix(displacement, true);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
            stiff.getRhs();

        Eigen::VectorXd groundtruth_result;

        bool bUsePardiso = true;

        if(bUsePardiso)
        {
            Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> pardisoSolver;
            pardisoSolver.compute(cpu_system.matrix());
            groundtruth_result = pardisoSolver.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (pardisoSolver.info() == Eigen::NoConvergence ||
                pardisoSolver.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Pardiso LDLT Solver Failed!");
            } else {
                GPMESH_TRACE("Eigen Pardiso LDLT Solve Success!");
            }
        }
        else
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            groundtruth_result = cg.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
            }
        }

        double updateNorm   = groundtruth_result.norm();
        double residualNorm = cpu_system.rhs().norm();

        if (newton_step == 0) {
            initUpdateNorm   = updateNorm;
            initResidualNorm = residualNorm;
        }

        m_solVector += groundtruth_result;

        bool bQuit1Cur = false;

        if (residualNorm < absTol || updateNorm < absTol ||
            residualNorm / initResidualNorm < relTol ||
            updateNorm / initUpdateNorm < relTol) {
            bQuit1Cur = true;
        }

        if (!bQuit1 && bQuit1Cur) {
            GPMESH_INFO("refactor solver should quit in newton step {}",
                        newton_step);
        }

        he_solver.computeUpdate();

        bool bQuit2Cur = false;

        if (he_solver.m_residualNorm < absTol ||
            he_solver.m_updateNorm < absTol ||
            he_solver.m_residualNorm / he_solver.m_initResidualNorm < relTol ||
            he_solver.m_updateNorm / he_solver.m_initUpdateNorm < relTol) {
            bQuit2Cur = true;
        }

        if (!bQuit2 && bQuit2Cur) {
            GPMESH_INFO("he_solver should quit in newton step {}", newton_step);
        }

        bQuit1 = bQuit1 || bQuit1Cur;
        bQuit2 = bQuit2 || bQuit2Cur;
        // GPMESH_INFO("In newton step {}", newton_step);
        // GPMESH_INFO("solVector Difference between refactor and he_solver =
        // {}", (m_solVector - he_solver.solution()).norm());

        m_sparse_system.matrix() = cpu_system.matrix();
        m_sparse_system.rhs()    = cpu_system.rhs();

        Eigen::SparseMatrix<double> A = m_sparse_system.matrix();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);

        Eigen::VectorXd b = m_sparse_system.rhs();

        Eigen::VectorXd s;

        Eigen::MatrixXd gpu_result;

        float gpu_time;

        debug_system_matrix_assembly<double>();

        // debug_test_system_matrix_gpu(cpu_system);
        //  m_assembly_context.update_host_info(m_h_assembly_info,
        //  m_num_geometry); m_num_geometry); cudaDeviceSynchronize();


        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_geometry_size;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        s = A * x;

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
            B_global->m_d_global_attr,
            S_global->m_d_global_attr,
            R_global->m_d_global_attr,
            P_global->m_d_global_attr,
            cps_length);

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        Eigen::VectorXd r = b - s;

        //        //test for R
        //        double test = 0;
        //        for (int i = 0; i < r.rows(); ++i)
        //        {
        //            double Ri = R_global->m_h_global_attr[i];
        //            double diff = Ri - r(i,0);
        //            diff = abs(diff);
        //            if(diff > 1e-16)
        //            {
        //                GPMESH_ERROR("R and r show diff={}", diff);
        //            }
        //
        //            test += Ri * Ri;
        //        }
        //
        //        GPMESH_INFO("R,r test success");
        //        GPMESH_INFO("test sn = {}", test);

        Eigen::VectorXd p = r;


        // GPMESH_INFO("init P R success");
        // CUDA_ERROR(cudaDeviceSynchronize());
        // CUDA_ERROR(cudaGetLastError());
        // 新方法计算Delta_new

        ReduceHandle<double> reduce_handle(*R_global);

        double Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

        double Delta_0_gpu = Delta_new_gpu;

        GPMESH_INFO("Delta_0 = {}", Delta_0_gpu);

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);

        GPMESH_INFO("delta_0 = {}", delta_0);

        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);
        S_global->reset_local_part(0.0, gpmesh::DEVICE);

        // CUDA_ERROR(cudaDeviceSynchronize());
        // CUDA_ERROR(cudaGetLastError());

        double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

        double alpha(0), beta(0), delta_old(0);

        uint32_t iterNum = 0;
        uint32_t maxIter = 1;

        while (iterNum < maxIter) {

            GPMESH_INFO("GPU Iter {}", iterNum);

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

            s = A * p;

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            GPMESH_INFO("Ap * p = {}", Alpha_gpu);

            alpha = p.dot(s);

            GPMESH_INFO("ap * p = {}", alpha);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            GPMESH_INFO("Alpha = {}", Alpha_gpu);

            alpha = delta_new / alpha;

            GPMESH_INFO("alpha = {}", alpha);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            x += alpha * p;
            r -= alpha * s;

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            GPMESH_INFO("Delta_new = {}", Delta_new_gpu);

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            GPMESH_INFO("delta_new = {}", delta_new);

            if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
                break;
            }

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;
            GPMESH_INFO("Beta = {}", Beta_gpu);

            beta = delta_new / delta_old;
            GPMESH_INFO("beta = {}", beta);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            p = r + beta * p;

            iterNum++;
        }

        timer_gpu.stop();
        gpu_time = timer_gpu.elapsed_millis();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("GPU takes {} iterations", iterNum);

        X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        X_global->move_global_to_local_gpu(m_attr_mapping_context);

        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            updateVector_cps(d, 0) = X_global->m_h_global_attr[d];
        }

        gpu_result = updateVector_cps;
        debug_update_cpoint_displacement_in_physic_context(*X_global);
    }

    // split cpu and gpu solving, try to find bugs
    // two test
    void debug_IGA_solving_process_test_v3(
        vector<CBSplineSolid>&           patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {

        init(patch, m_dirichletSides, m_neumannSides, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
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

        {
            Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
                mass.matrix());
            accVector     = acc_solver.solve(stiff.getRhs());
            pre_accVector = accVector;
        }

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 1;

        SparseSystem& cpu_system = stiff.system();

        double absTol = 1e-12;
        double relTol = 1e-9;

        bool bQuit1 = false;

        double initUpdateNorm;
        double initResidualNorm;

        // gpu backend
        debug_mass_matrix_assembly<double>(density);

        // stiff.assemble(m_solVector)
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        displacement.setZero();
        stiff.constructCoefs(m_solVector, displacement);
        stiff.stiffNonlinearMatrix(displacement, true);
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        cpu_system.matrix() =
            alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
        cpu_system.matrix().makeCompressed();
        cpu_system.rhs() =
            mass.matrix() *
                (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                 alpha2(Beta, tStep) * velVector + alpha3(Beta) * accVector) +
            stiff.getRhs();

        Eigen::VectorXd groundtruth_result;

        bool bUsePardiso = true;

        if(bUsePardiso)
        {
            Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> pardisoSolver;
            pardisoSolver.compute(cpu_system.matrix());
            groundtruth_result = pardisoSolver.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (pardisoSolver.info() == Eigen::NoConvergence ||
                pardisoSolver.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Pardiso LDLT Solver Failed!");
            } else {
                GPMESH_TRACE("Eigen Pardiso LDLT Solve Success!");
            }
        }
        else
        {
            Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> cg;
            cg.compute(cpu_system.matrix());
            groundtruth_result = cg.solve(cpu_system.rhs());
            GPMESH_INFO("updateVector groundtruth norm = {}",
                        groundtruth_result.norm());
            if (cg.info() == Eigen::NoConvergence ||
                cg.info() == Eigen::NumericalIssue) {
                GPMESH_ERROR("Eigen Conjugate Gradient Solver Failed!");
                GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
            } else {
                GPMESH_TRACE("Eigen Conjugate Gradient Solve Success!");
                GPMESH_TRACE("Iteration {}, Error {}", cg.iterations(), cg.error());
            }
        }

        double updateNorm   = groundtruth_result.norm();
        double residualNorm = cpu_system.rhs().norm();

        if (newton_step == 0) {
            initUpdateNorm   = updateNorm;
            initResidualNorm = residualNorm;
        }

        m_solVector += groundtruth_result;

        bool bQuit1Cur = false;

        if (residualNorm < absTol || updateNorm < absTol ||
            residualNorm / initResidualNorm < relTol ||
            updateNorm / initUpdateNorm < relTol) {
            bQuit1Cur = true;
        }

        if (!bQuit1 && bQuit1Cur) {
            GPMESH_INFO("refactor solver should quit in newton step {}",
                        newton_step);
        }

        bQuit1 = bQuit1 || bQuit1Cur;

        // GPMESH_INFO("In newton step {}", newton_step);
        // GPMESH_INFO("solVector Difference between refactor and he_solver =
        // {}", (m_solVector - he_solver.solution()).norm());

        m_sparse_system.matrix() = cpu_system.matrix();
        m_sparse_system.rhs()    = cpu_system.rhs();

        Eigen::SparseMatrix<double>& A = m_sparse_system.matrix();

        Eigen::VectorXd x = Eigen::VectorXd::Zero(m_dims * m_num_cpoints);
        //x = 0.000000001 * Eigen::VectorXd::Ones(m_num_cpoints * m_dims);

        Eigen::VectorXd b = m_sparse_system.rhs();

        Eigen::VectorXd s;

        Eigen::MatrixXd gpu_result;

        float gpu_time;

        debug_system_matrix_assembly<double>();

        //debug_test_system_matrix_gpu(cpu_system);
        //  m_assembly_context.update_host_info(m_h_assembly_info,
        //  m_num_geometry); m_num_geometry); cudaDeviceSynchronize();


        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;
        size_t             geometryNum      = m_num_geometry;

        Eigen::MatrixXd updateVector_cps;
        updateVector_cps.setZero(m_num_cpoints * m_dims, 1);
        //updateVector_cps.setOnes(m_num_cpoints * m_dims, 1);
        //updateVector_cps *= 0.000000001;

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
        IGA_matvec_A_multiply_v1<double, 192><<<geometryNum, 192>>>(
            get_spline_context(), get_assembly_context(), *X_global, *S_global);

        S_global->move_local_part(gpmesh::DEVICE, gpmesh::HOST);
        S_global->move_local_to_global_gpu(m_attr_mapping_context);

        s = A * x;

        // GPMESH_INFO("S calculate success");

        // init PR
        size_t cps_length          = m_num_cpoints * 3;
        size_t init_PR_threads_num = 256;
        size_t init_PR_blocks_num  = DIVIDE_UP(cps_length, 256);

        // R_global = B_global - S_global
        // P_global = R_global
        init_PR<double><<<init_PR_blocks_num, init_PR_threads_num>>>(
            B_global->m_d_global_attr,
            S_global->m_d_global_attr,
            R_global->m_d_global_attr,
            P_global->m_d_global_attr,
            cps_length);

        R_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        P_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        Eigen::VectorXd r = b - s;

        Eigen::VectorXd p = r;

        ReduceHandle<double> reduce_handle(*R_global);

        double Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

        double Delta_0_gpu = Delta_new_gpu;

        GPMESH_INFO("Delta_0 = {}", Delta_0_gpu);

        double delta_new = r.squaredNorm();
        double delta_0(delta_new);

        GPMESH_INFO("delta_0 = {}", delta_0);

        // move P_global global part to local part
        P_global->move_global_to_local_gpu(m_attr_mapping_context);

        // CUDA_ERROR(cudaDeviceSynchronize());
        // CUDA_ERROR(cudaGetLastError());

        double Alpha_gpu(0), Beta_gpu(0), Delta_old_gpu(0);

        double alpha(0), beta(0), delta_old(0);

        uint32_t iterNum = 0;
        uint32_t maxIter = 1;

        s.setZero();

        while (iterNum < maxIter) {

            GPMESH_INFO("GPU Iter {}", iterNum);

            S_global->reset_local_part(0.0, gpmesh::DEVICE);
            S_global->reset_global_part(0.0, gpmesh::DEVICE);

            debug_test_system_matrix_gpu(cpu_system);
            GPMESH_INFO("test for s");
            //debug_test_vector_gpu(s, S_global);
            GPMESH_INFO("test for p");
            //debug_test_vector_gpu(p, P_global);

            IGA_matvec_A_multiply_v1<double, 192>
                <<<geometryNum, 192>>>(get_spline_context(),
                                       get_assembly_context(),
                                       *P_global,
                                       *S_global);

            S_global->move_local_to_global_gpu(m_attr_mapping_context);
            S_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);


            s = A * p;
            //debug_test_vector_gpu(s, S_global);

            Alpha_gpu = reduce_handle.dot(*S_global, *P_global);

            GPMESH_INFO("Ap * p = {}", Alpha_gpu);

            alpha = p.dot(s);

            GPMESH_INFO("ap * p = {}", alpha);

            Alpha_gpu = Delta_new_gpu / Alpha_gpu;

            GPMESH_INFO("Alpha = {}", Alpha_gpu);

            alpha = delta_new / alpha;

            GPMESH_INFO("alpha = {}", alpha);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *X_global,
                *P_global,
                Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *R_global,
                *S_global,
                -Alpha_gpu,
                1.0,
                X_global->m_num_cpoints * 3);

            x += alpha * p;
            r -= alpha * s;

            Delta_old_gpu = Delta_new_gpu;

            Delta_new_gpu = reduce_handle.squaredNorm(*R_global);

            GPMESH_INFO("Delta_new = {}", Delta_new_gpu);

            delta_old = delta_new;
            delta_new = r.squaredNorm();

            GPMESH_INFO("delta_new = {}", delta_new);

            if (Delta_new_gpu < 1e-13 * Delta_0_gpu) {
                break;
            }

            Beta_gpu = Delta_new_gpu / Delta_old_gpu;
            GPMESH_INFO("Beta = {}", Beta_gpu);

            beta = delta_new / delta_old;
            GPMESH_INFO("beta = {}", beta);

            axpy<double><<<reduce_handle.get_num_block(),
                           reduce_handle.get_num_thread()>>>(
                *P_global,
                *R_global,
                1.0,
                Beta_gpu,
                X_global->m_num_cpoints * 3);

            P_global->move_global_to_local_gpu(m_attr_mapping_context);

            p = r + beta * p;

            iterNum++;
        }

        timer_gpu.stop();
        gpu_time = timer_gpu.elapsed_millis();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("GPU takes {} iterations", iterNum);

        X_global->move_global_part(gpmesh::DEVICE, gpmesh::HOST);
        X_global->move_global_to_local_gpu(m_attr_mapping_context);

        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            updateVector_cps(d, 0) = X_global->m_h_global_attr[d];
        }

        gpu_result = updateVector_cps;
        debug_update_cpoint_displacement_in_physic_context(*X_global);
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

    void debug_test_vector_gpu(Eigen::VectorXd& x, shared_ptr<SplineAttributeGlobal<double>> X)
    {
        //X->move_local_to_global_gpu(m_attr_mapping_context);
        //X->move_global_part(gpmesh::DEVICE, gpmesh::HOST);

        double test = 0;

        for (int i = 0; i < x.rows(); ++i)
        {
            double Xi = X->m_h_global_attr[i];
            double diff = Xi - x(i);

            if(diff > 1e-10)
            {
                GPMESH_WARN("data[{}] show difference = {}", i, diff);
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

    void debug_matvec_peformance_test(
        const vector<CBSplineSolid>&     patch,
        std::vector<boundary_condition>& m_dirichletSides,
        std::vector<boundary_condition>& m_neumannSides)
    {
        constexpr uint32_t gaussianPointNum = 64;
        const int          localCPointNum   = 64;
        const double       MAXNUM_DIRICHLET = 1e10;

        init(patch, m_dirichletSides, m_neumannSides, false);

        GPUTimer timer_gpu;

        double timeSpan = 1;
        double timeStep = 0.1;
        double density  = 3320;

        double Beta  = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);
        Eigen::MatrixXd accVector;
        Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_num_cpoints * m_dims, 1);

        Eigen::MatrixXd pre_disVector = dispVector;
        Eigen::MatrixXd pre_velVector = velVector;

        // Todo: 初始化加速度计算 - timeAssemble.initialize()
        // 加速度计算 a = F / a

        double residualNorm     = 0;
        double initResidualNorm = 0;
        double updateNorm       = 0;
        double initUpdateNorm   = 1;

        // using saint_venant_kirchhoff material

        dim3 num_blocks(m_num_geometry, 64, 1);
        dim3 threads_per_blocks(16, 16, 1);

        size_t smem_bytes_dyn_hyper = calc_shared_memory(
            (void*)
                IGA_hyperelasticity_assembly_debug_kernel_v2<double,
                                                             gaussianPointNum>);

        timer_gpu.start();
        IGA_hyperelasticity_assembly_debug_kernel_v2<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_assembly_context());

        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());
        GPMESH_INFO(
            "IGA_hyperelasticity_assembly_debug_kernel_v2 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        const int DirichNum = m_dirichlet_context.m_num_dirichlet_boundary;
        const int GeometrySurfNum = 16;
        timer_gpu.start();
        IGA_add_dirichlet_to_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(),
            get_assembly_context(),
            get_dirichlet_context());
        timer_gpu.stop();
        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_add_dirichlet_to_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        timer_gpu.start();
        IGA_mass_assembly_debug_kernel_v0<double, gaussianPointNum>
            <<<num_blocks, threads_per_blocks, smem_bytes_dyn_hyper>>>(
                get_spline_context(),
                get_gaussian_integral_context(),
                get_physic_context(),
                get_mass_context(),
                density);
        timer_gpu.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        GPMESH_INFO("IGA_mass_assembly_debug_kernel_v0 Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        timer_gpu.start();
        IGA_add_dirichlet_to_mass_local_mats<<<DirichNum, GeometrySurfNum>>>(
            get_spline_context(), get_mass_context(), get_dirichlet_context());
        timer_gpu.stop();

        GPMESH_INFO("IGA_add_dirichlet_to_mass_local_mats Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_mass_context.update_host_info(m_h_mass_info, m_num_geometry);

        timer_gpu.start();
        IGA_add_mass_matrix_to_stiff_matrix<double, 256>
            <<<m_num_geometry, 256>>>(get_assembly_context(),
                                      get_mass_context(),
                                      alpha1(Beta, tStep));
        timer_gpu.stop();
        GPMESH_INFO("IGA_add_mass_matrix_to_stiff_matrix Run Success!");
        GPMESH_TRACE("Kernel took {} (ms)", timer_gpu.elapsed_millis());

        m_assembly_context.update_host_info(m_h_assembly_info, m_num_geometry);

        vector<Eigen::MatrixXd> mass_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            mass_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_mass_info[g].mass_local_mats, 192, 192);
        }

        vector<Eigen::MatrixXd> system_matrix_gpu(m_num_geometry);

        for (int g = 0; g < m_num_geometry; ++g) {
            system_matrix_gpu[g] = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
        }

        std::vector<int> localIndex(64, 0);
        std::iota(localIndex.begin(), localIndex.end(), 0);

        for (int g = 0; g < m_num_geometry; ++g) {
            Eigen::MatrixXd localMats;

            localMats = Eigen::Map<Eigen::MatrixXd, Eigen::ColMajor>(
                m_h_assembly_info[g].local_mats, 192, 192);
            m_sparse_system.pushToMatrix(g, localIndex, localMats);
        }

        auto SystemBigMatrxiFromGPU = m_sparse_system.test_get_bigmatrix();
        // here do not adding dirichlet boundary information again
        m_sparse_system.dataPush();
        auto SystemSparseMatrixFromGPU = m_sparse_system.matrix();

        // cpu assemble

        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);

        // 质量矩阵
        massAssembler mass(patch, density, stiff.system());

        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        auto velVec  = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);

        // displacment = solution = m_solVec
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);

        // initialize

        vector<Eigen::MatrixXd> mass_matrix_cpu;

        mass.assemble(mass_matrix_cpu);

        // mass matrix test
        bool bMassVerify = true;
        for (int g = 0; g < m_num_geometry; ++g) {
            double diff_mass_matrix =
                (mass_matrix_gpu[g] - mass_matrix_cpu[g]).norm();
            if (diff_mass_matrix > 1e-10) {
                bMassVerify = false;
                break;
            }
        }

        if (!bMassVerify)
            GPMESH_ERROR("mass matrix verify fail");


        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
            displacement, true, test_localMatsVec, test_localRhsVec);


        // 外力rhs
        stiff.rhsVector();
        stiff.system().dataPush_Dirichlet();

        // a = F / m
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> solver(
            mass.matrix());
        Eigen::MatrixXd              accVec      = solver.solve(stiff.getRhs());
        Eigen::MatrixXd              pre_accVect = accVec;
        Eigen::SparseMatrix<double>  A_stiff     = stiff.getMatrix();
        Eigen::SparseMatrix<double>& GroundTruthSparseMatrix = A_stiff;


        // 实现一个时间步的弹性体模拟

        cout << "Step 0!" << endl;

        // Hyperelastic.solve()

        int local_num_iters = 0;
        int max_iters       = 1;

        SparseSystem cpu_system;

        while (local_num_iters < max_iters) {
            // 进行一次牛顿迭代求解计算

            // 1.重新进行assemble
            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector();
            stiff.system().dataPush_Dirichlet();
            cpu_system.matrix() =
                alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();
            cpu_system.matrix().makeCompressed();
            cpu_system.rhs() =
                mass.matrix() *
                    (alpha1(Beta, tStep) * displacement +
                     alpha2(Beta, tStep) * velVec + alpha3(Beta) * accVec) +
                stiff.getRhs();

            local_num_iters++;
        }

        // 和 timeIntegrator 进行对比
        timeIntegrator timeAssemble(stiff, mass);
        timeAssemble.setInitialDisplacement(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.setInitialVelocity(
            Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
        timeAssemble.initialize();
        timeAssemble.assemble(timeAssemble.displacementVector(), true);

        // 对比 cpu_system 和 timeAssemble.system()

        double diff_mat =
            (timeAssemble.getMatrix() - cpu_system.matrix()).norm();
        double diff1_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - timeAssemble.getMatrix()).norm();
        double diff2_system_mat_between_cpu_and_gpu =
            (SystemSparseMatrixFromGPU - cpu_system.matrix()).norm();
        double diff_rhs = (timeAssemble.getRhs() - cpu_system.rhs()).norm();

        m_sparse_system.rhs() = timeAssemble.getRhs();

        GPMESH_INFO("diff_mat = {}, diff_rhs = {}", diff_mat, diff_rhs);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff1_system_mat_between_cpu_and_gpu);
        GPMESH_INFO("system mat difference between cpu and gpu is {}",
                    diff2_system_mat_between_cpu_and_gpu);

        this->helper_performance_test_for_matvec_A_multiply();
        this->helper_performance_test_for_matvec_A_multiply_v0();
        this->helper_performance_test_for_matvec_A_multiply_v1();
        //            Eigen::MatrixXd gpu_result;
        //            float gpu_time;
        //            this->debug_solving_assembly_rhs_helper_gpu(gpu_result,
        //            gpu_time); Eigen::MatrixXd cpu_result; float cpu_time;
        //            this->debug_solving_assembly_rhs_helper_cpu(cpu_result,
        //            cpu_time); Eigen::MatrixXd groundtruth_result; float
        //            groundtruth_time;
        //            this->debug_solving_assembly_rhs_helper_groudtruth(groundtruth_result,
        //            groundtruth_time);
        //
        //            double diff_cpu_res = (groundtruth_result -
        //            cpu_result).norm(); double diff_gpu_res =
        //            (groundtruth_result - gpu_result).norm();
        //
        //            GPMESH_INFO("cpu calculation difference - {}",
        //            diff_cpu_res); GPMESH_INFO("gpu calculation difference -
        //            {}", diff_gpu_res);
        //
        //            GPMESH_INFO("cpu time - {} (ms)", cpu_time);
        //            GPMESH_INFO("gpu time - {} (ms)", gpu_time);
        //            GPMESH_INFO("groundtruth time - {} (ms)",
        //            groundtruth_time);
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


        size_t smem_bytes_static  = func_attr.sharedSizeBytes;
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

        vector<int> iterGroup = {0, 64, 128};

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

        vector<int> iterGroup = {0, 64, 128};

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
        Eigen::MatrixXd&                out)
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
        Eigen::MatrixXd&                in_coefs,
        Eigen::MatrixXd&                in_cps)
    {
        for (uint32_t g = 0; g < m_num_geometry; ++g) {
            for (uint16_t local_id = 0; local_id < 64; ++local_id) {
                uint32_t global_id       = m_h_geometry_ltog_cp[g][local_id];
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
        Eigen::MatrixX<T>&              in_cps)
    {
        for (uint32_t g = 0; g < m_num_geometry; ++g) {
            for (uint16_t local_id = 0; local_id < 64; ++local_id) {
                uint32_t global_id       = m_h_geometry_ltog_cp[g][local_id];
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
    void helper_convert_tripVec_to_cps(
        std::vector<std::vector<std::vector<T>>>& tripleVec,
        Eigen::MatrixX<T>&                        cps)
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
        Eigen::MatrixX<T>&                        cps,
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
                  const Eigen::VectorXd&             b,
                  Eigen::VectorXd&                   x)
    {
        Eigen::VectorXd r     = b - A * x;
        Eigen::VectorXd p     = r;
        double          rsold = r.squaredNorm();
        double          tol   = 1e-6;

        double rs_0 = rsold;

        for (int i = 0; i < b.size(); ++i) {
            Eigen::VectorXd Ap = A * p;

            double alpha = rsold / (p.template dot(Ap));

            x.noalias()  = x + alpha * p;
            r.noalias()  = r - alpha * Ap;
            double rsnew = r.squaredNorm();

            if (rsnew < tol * tol * rs_0) {
                GPMESH_TRACE("conjgrad solve success for iter {}", i);
                break;
            }

            p     = r + (rsnew / rsold) * p;
            rsold = rsnew;
        }

        GPMESH_TRACE("conjgrad solve result rs={}", rsold);
    }

    template <typename T>
    void preconditioned_conjgrad(const Eigen::SparseMatrix<double>& A,
                  const Eigen::VectorXd&             b,
                  Eigen::VectorXd&                   x)
    {
        Eigen::VectorXd r     = b - A * x;
        Eigen::VectorXd p     = r;
        double          rsold = r.squaredNorm();
        double          tol   = 1e-6;

        Eigen::VectorXd M;
        M.resize(p.size());

        double rs_0 = rsold;

        for (int i = 0; i < b.size(); ++i) {
            Eigen::VectorXd Ap = A * p;

            double alpha = rsold / (p.template dot(Ap));

            x.noalias()  = x + alpha * p;
            r.noalias()  = r - alpha * Ap;
            double rsnew = r.squaredNorm();

            if (rsnew < tol * tol * rs_0) {
                GPMESH_TRACE("conjgrad solve success for iter {}", i);
                break;
            }

            p     = r + (rsnew / rsold) * p;
            rsold = rsnew;
        }

        GPMESH_TRACE("conjgrad solve result rs={}", rsold);
    }
};

}  // namespace gpmesh