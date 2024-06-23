//#include<Hyperelasticity/utils.h>
//#include<Hyperelasticity/TimeIntegrator.h>
//#include<chrono>
//#include<Hyperelasticity/collision_detect.h>
//#include <iostream>
//#include "Util/GP_Log.h"
//#include "Util/GP_Macros.h"
//#include "Util/GP_Util.h"
//#include "GP_Local.h"
//#include "GP_Patch_Info.h"
//#include "GP_Context.h"
//#include "Patcher/GP_Patcher.h"
//#include "GP_Mesh_Static.h"
//#include "GP_vertex_normal_kernel.cuh"
//#include "GP_Spline.h"
//#include "GP_Spline_Static.h"
//#include "GP_Spline_Test.h"
//#include "GP_Spline_v0.h"
//
//#include "GP_Invoke.h"
//#include "GP_Spline_t.h"
//#include "GP_Param.h"
//#include "GP_BoundaryCPoint.h"
//
//int dim = 3;
//
//int main() {
//
//    Eigen::setNbThreads(24);
//    Eigen::initParallel();
//
//    double timeSpan = 1;
//    double timeStep = 0.1;
//    double density = 3320;
//    double gravitationalAcc = 9.8;
//    std::vector<CBSplineSolid> patch;
//    bool bH = true;
//    if(bH)
//    {
//        ToolFunc::readFileXml(get_model_name(0), patch);
//    }
//    else
//    {
//        ToolFunc::readFileXml(get_model_name(0), patch);
//    }
//
//    std::vector<boundary_condition> m_dirichletSides, m_neumannSides;
//
//    bool bHyperElasticity = false;
//
//    if (bHyperElasticity) {
//        CVector traction(0, 1500, -1500);
//        CVector traction1(10000, 0, 0);
//        CVector traction2(3000, -3000, 0);
//
//        if (bH) {
//            //H model
//            for (int d = 0; d < 3; ++d) {
//                m_dirichletSides.push_back(boundary_condition(0, boundary::back, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(1, boundary::front, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(1, boundary::north, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(1, boundary::south, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(1, boundary::back, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(2, boundary::back, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(6, boundary::back, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(7, boundary::back, condition_type::dirichlet, 0, d));
//            }
//            m_neumannSides.push_back(boundary_condition(3, boundary::front, condition_type::neumann, traction1));
//            //m_neumannSides.push_back(boundary_condition(9, boundary::north, condition_type::neumann, traction1));
//            //m_neumannSides.push_back(boundary_condition(4, boundary::front, condition_type::neumann, traction1));
//            //m_neumannSides.push_back(boundary_condition(10, boundary::north, condition_type::neumann, traction1));
//            m_neumannSides.push_back(boundary_condition(5, boundary::front, condition_type::neumann, traction1));
//        }
//        else
//        {
//
//            //konglong model
//            for (int d = 0; d < 3; ++d) {
//                m_dirichletSides.push_back(boundary_condition(31, boundary::south, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(31, boundary::north, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(31, boundary::front, condition_type::dirichlet, 0, d));
//                m_dirichletSides.push_back(boundary_condition(31, boundary::back, condition_type::dirichlet, 0, d));
//            }
//
//            m_neumannSides.push_back(boundary_condition(36, boundary::front, condition_type::neumann, traction2));
//            m_neumannSides.push_back(boundary_condition(37, boundary::front, condition_type::neumann, traction2));
//        }
//
//    } else {
//        CVector traction(0, 15, -15);
//        CVector traction1(0, -15, -15);
//
//        for (int d = 0; d < 3; ++d) {
//
//
//            m_dirichletSides.push_back(boundary_condition(0, boundary::west, condition_type::dirichlet, 0, d));
//            //m_dirichletSides.push_back(boundary_condition(5, boundary::west, condition_type::dirichlet, 0, d));
//
//        }
//
//        //for  IGA/test.xml
//        m_neumannSides.push_back(boundary_condition(0, boundary::east, condition_type::neumann, traction));
//    }
//
//    bool bUseBoundaryPoint =  true;
//
//    bool bOriginTest = false;
//
//    if (bOriginTest) {
//        stiffAssembler stiff(patch, m_dirichletSides, m_neumannSides);
//
//        //Eigen::MatrixXd displacement;
//        Eigen::MatrixXd u;
//
//        massAssembler mass(patch, density, stiff.system());
//
//
//        timeIntegrator timeAssemble(stiff, mass);
//        timeAssemble.setInitialDisplacement(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
//        timeAssemble.setInitialVelocity(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
//        timeAssemble.initialize();
//
//        timeAssemble.assemble(timeAssemble.displacementVector(), true);
//
//        //stiff.WriteParaview("D:/Resource/FileDataOutput/disp/konglong",u,1000,false );
//        Eigen::MatrixXd displacement1;
//        //cout << "start to output" << endl;
//        //stiff.WriteParaview("D:/Resource/FileDataOutput/xie/Lmesh", u, 1000, false);
//        //////Timer<std::chrono::milliseconds> timer1;
//
//        for (int t = 0; t < (int) timeSpan / timeStep; ++t) {
//            std::cout << t << "\n";
//            timeAssemble.makeTimeStep(timeStep);
//            //std::cout << .topRows(10);
//            stiff.constructCoefs(timeAssemble.displacementVector(), displacement1);
//            //stiff.constructCoefs(u, displacement1);
//            u = timeAssemble.displacementVector();
//            int begin = 0;
//            for (size_t i = 0; i < patch.size(); ++i) {
//                int num = patch[i].getSizeControlPoints();
//                patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
//                begin += num;
//            }
//            //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/xiao/threeblock2_disp");
//            //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/simulate/C_deal_fail");
//            u = timeAssemble.displacementVector();
//            stiff.WriteParaview(patch, "B:/Win11Projects/OpenCascadeUse/FileDataOutput/disp/H" + to_string(t) + "_", u,
//                                1000, false);
//            begin = 0;
//            displacement1 *= -1;
//            for (size_t i = 0; i < patch.size(); ++i) {
//                int num = patch[i].getSizeControlPoints();
//                patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
//                begin += num;
//            }
//        }
//
//        bool runSolver = false;
//
//        if (runSolver) {
//            for (int t = 0; t < (int) timeSpan / timeStep; ++t) {
//                std::cout << t << "\n";
//                timeAssemble.makeTimeStep(timeStep);
//                //std::cout << .topRows(10);
//                stiff.constructCoefs(timeAssemble.displacementVector(), displacement1);
//                //stiff.constructCoefs(u, displacement1);
//                u = timeAssemble.displacementVector();
//                int begin = 0;
//                for (size_t i = 0; i < patch.size(); ++i) {
//                    int num = patch[i].getSizeControlPoints();
//                    patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
//                    begin += num;
//                }
//                //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/xiao/threeblock2_disp");
//                //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/simulate/C_deal_fail");
//                u = timeAssemble.displacementVector();
//                std::string output_path = STRINGIFY(OUTPUT_DIR) "IGADataOutput/H";
//                stiff.WriteParaview(patch, output_path + to_string(t) + "_", u, 1000, false);
//                begin = 0;
//                displacement1 *= -1;
//                for (size_t i = 0; i < patch.size(); ++i) {
//                    int num = patch[i].getSizeControlPoints();
//                    patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
//                    begin += num;
//                }
//            }
//        } else {
//
//        }
//    }
//
//    gpmesh::GPLog::init();
//    GPMESH_TRACE("test");
//
//    CUDA_ERROR(cudaDeviceReset());
//
//    //gpmesh::GPSplineSolidTest m_test_spline;
//
//    //m_test_spline.debug_IGA_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_cg_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_dirichlet_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_cg_solve_test_v0(patch, m_dirichletSides, m_neumannSides);
//
//    //m_test_spline.debug_IGA_hyperE_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_hyperelasticity_assembly_dirichlet_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_hyperelasticity_solving_v0(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_matvec_peformance_test(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_he_cg_performance_test_v0(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_hyperelasticity_frame_process_v5(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_hyperelasticity_frame_process_v4(patch, m_dirichletSides, m_neumannSides);
//
//    //m_test_spline.debug_IGA_solving_process_test_v3(patch, m_dirichletSides, m_neumannSides);
//    //u = timeAssemble.displacementVector();
//    //Eigen::MatrixXd displacewment;
//    //int begin = 0;
//    //stiff.constructCoefs(u, displacement);
//    //for (size_t i = 0; i < patch.size(); ++i) {
//    //	int num = patch[i].getSizeControlPoints();
//    //	patch[i].updateCoefs(displacement.block(begin, 0, num, 3));
//    //	begin += num;
//    //}
//    //stiff.WriteParaview(patch, "D:/Resource/FileDataOutput/disp/H", u, 1000, false);
//
//    std::vector<boundary_point> neumann_points;
//    std::vector<boundary_point> dirichlet_points;
//
//    if(bUseBoundaryPoint)
//    {
//        helper_provide_neumann_dirichlet_points(neumann_points, dirichlet_points);
//    }
//
//    gpmesh::GPSplineSolidV0 m_test_spline;
//    //m_test_spline.debug_IGA_hyperE_assembly(patch, m_dirichletSides, m_neumannSides);
//    //m_test_spline.debug_IGA_test_optimized_assembly_v0(patch, neumann_points, dirichlet_points);
//    m_test_spline.debug_IGA_boundary_point_test_v0(patch, neumann_points, dirichlet_points);
//}