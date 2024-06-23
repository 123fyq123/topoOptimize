#include<Hyperelasticity/utils.h>
#include<Hyperelasticity/TimeIntegrator.h>
#include<chrono>
#include<Hyperelasticity/collision_detect.h>
#include <iostream>
#include "Util/GP_Log.h"
#include "Util/GP_Macros.h"

#include "Util/GP_Util.h"
#include "GP_Local.h"
#include "GP_Patch_Info.h"
#include "GP_Context.h"
#include "Patcher/GP_Patcher.h"
#include "GP_Mesh_Static.h"
#include "GP_vertex_normal_kernel.cuh"
#include "GP_Spline.h"
#include "GP_Spline_Static.h"

#include "GP_Spline_Test.h"
#include "GP_Spline_v0.h"
#include "GP_Invoke.h"
#include "GP_Spline_t.h"
#include "GP_Param.h"
#include "GP_BoundaryCPoint.h"
#include "GP_nsp.h"
#include "Poster/source/Timer.h"
#include "omp.h"
#include "Poster/source/nanoflann.h"

#include <algorithm>
#include <vector>


#include "Cuda_Hyperelasticity_Solver.h"

#include <ios>

void Cuda_Hyperelasticity_Solver(const std::string& input_model_file,
                                 const std::string& boundary_file,
                                 const std::string& output_file,
                                 double Yangs, double Poissons,
                                 int frameCount, double timeStep,
                                 HexV3fMesh* mesh,
                                 std::vector<double> analysis_data[3],
                                 std::atomic_bool &is_need_update_model)
{
    Timer timer;
    Eigen::setNbThreads(24);
    Eigen::initParallel();

    std::vector<CBSplineSolid> patch;
    bool bH = true;


    ToolFunc::readFileXml(input_model_file, patch);
    timer.out("read_xml");

    bool bUseBoundaryPoint = true;

    std::vector<CPoint> bdy_points;
    std::vector<int> bdy_type;
    std::vector<CVector> bdy_value;

    // add bdy info
    {
        int control_size = 0;
        std::ifstream in(boundary_file);
        in >> control_size;
        std::cout <<"control_size = " << control_size <<"\n";

        for(int i = 0; i < control_size; ++i) {
            double x, y, z;
            int cur_type;
            double vx, vy, vz;
            in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
            CPoint cp(x, y, z);
            CVector cv(vx, vy, vz);
            bdy_points.push_back(cp);
            bdy_value.push_back(cv);
            bdy_type.push_back(cur_type);
        }
    }
    timer.out("add bdy info");
    gpmesh::GPLog::init(spdlog::level::trace);

    CUDA_ERROR(cudaDeviceReset());
    std::vector<boundary_point> neumann_points;
    std::vector<boundary_point> dirichlet_points;
    gpmesh::GPSplineSolidV0 m_spline;
//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
//    timer.out("TopOpt_GPU_Init");
    m_spline.init_patch(patch);
    timer.out("Linear Elasticity GPU Geometry Init");

    timer.out("before sys");
    const SparseSystem& temp_sys = m_spline.m_sparse_system;
    timer.out("after sys");
    {


        std::vector<gpmesh::Point> A, B;
        std::vector<int> global_ids;

        for(auto& p : bdy_points){
            A.emplace_back(gpmesh::Point{p.getX(), p.getY(), p.getZ()});
        }
        for(int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                auto& pp = patch[p_id].coefs()[local_id];
                B.emplace_back(gpmesh::Point{pp.getX(), pp.getY(), pp.getZ()});
            }
        }

        GetNearestIds(A, B, global_ids);
        timer.out("after GetNearestIds");
        std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
        auto& global_local_map = temp_sys.get_l2g_map();

        for(int i = 0; i < bdy_points.size(); ++i){
            auto& p = bdy_points[i];
            int patch_id = global_ids[i] / 64;
            int local_id = global_ids[i] % 64;
            int global_id = global_local_map[patch_id][local_id];

            if(global_id == -1){
                std::cerr <<"cant find control points\n";
                //std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
                std::cout << p.getX() <<" " <<p.getY() <<" " << p.getZ() <<"\n";
            }else{
                if(bdy_type[i] == 2){
                    neumann_points.emplace_back(boundary_point(condition_type::neumann, bdy_value[i],global_id));
                }else if(bdy_type[i] == 1){
                    dirichlet_points.emplace_back(boundary_point(condition_type::dirichlet, CVector(0, 0, 0) ,global_id));
                }
            }
        }
        timer.out("after bdy_points");
    }

    m_spline.init_dirichlet_points(dirichlet_points);

    m_spline.setYoungPossion(Yangs, Poissons);

    Eigen::MatrixXd UU;

    double density  = 3320;


    timer.out("Solver Mass Init");


        // vel, acc, disp should be global variable
    Eigen::MatrixXd velVector =
            Eigen::MatrixXd::Zero(m_spline.get_num_cpoints() * 3, 1);
    Eigen::MatrixXd accVector;
    Eigen::MatrixXd dispVector =
            Eigen::MatrixXd::Zero(m_spline.get_num_cpoints()  * 3, 1);

    stiffAssembler stiff(patch, dirichlet_points, neumann_points);
    stiff.setYangsPossion(Yangs, Poissons);

    stiff.system().matrix().setZero();
    stiff.system().rhs().setZero();

    massAssembler mass(patch, density, stiff.system());

    mass.assemble();

    Eigen::MatrixXd displacement;
    stiff.constructCoefs(dispVector, displacement);

    stiff.stiffNonlinearMatrix(displacement, true);
    stiff.rhsVector_NeumannPoint();
    stiff.system().dataPush_Dirichlet();

    Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> acc_solver(
            mass.matrix());
    accVector = acc_solver.solve(stiff.getRhs());

    Eigen::MatrixXd m_solVector = dispVector;
    Eigen::MatrixXd cpuDisp = dispVector;
    double Beta  = 0.25;
    double Gamma = 0.5;
    double tStep = timeStep;

    m_spline.debug_mass_matrix_assembly_dirichlet_point(density);

    for(int frameId = 0; frameId < frameCount; ++frameId)
    {

        uint32_t newton_step  = 0;
        uint32_t newton_limit = 3;                                                                                                                                                                                                                                                                             ;

        double absTol = 1e-12;
        double relTol = 1e-9;

        while(newton_step < newton_limit)
        {
            newton_step++;

            displacement.setZero();
            stiff.constructCoefs(m_solVector, displacement);

            stiff.system().matrix().setZero();
            stiff.system().rhs().setZero();
            //You should notice that nonlinear matrix assembly affect rhs
            stiff.stiffNonlinearMatrix(displacement, true);
            stiff.rhsVector_NeumannPoint();
            stiff.system().dataPush_Dirichlet();

            m_spline.m_sparse_system.matrix() = alpha1(Beta, tStep) * mass.matrix() + stiff.getMatrix();

            m_spline.m_sparse_system.rhs() = mass.matrix() *
                                             (alpha1(Beta, tStep) * (dispVector - m_solVector) +
                                              alpha2(Beta, tStep) * velVector +
                                              alpha3(Beta) * accVector) +
                                             stiff.getRhs();

            m_spline.debug_system_matrix_assembly_dirichlet_point_v2<double>(tStep);
            Eigen::MatrixXd gpu_result;
            float gpu_time;
            m_spline.debug_frame_process_PCG_jacobiPrecondition(gpu_result, gpu_time, false);

            Eigen::MatrixXd cpu_result;
            float cpu_time;

            m_spline.debug_frame_process_Pardiso(cpu_result, cpu_time);
            GPMESH_INFO("update vector norm is {}", gpu_result.norm());
            GPMESH_INFO("update vector difference between cpu and gpu - {}",
                        (gpu_result - cpu_result).norm());

            m_solVector += gpu_result;
            cpuDisp += cpu_result;

            GPMESH_INFO("Disp difference between cpu and gpu - {}",
                        (m_solVector - cpuDisp).norm());
        }

        GPMESH_INFO("Newton Step-{} Finish", newton_step);

        GPMESH_INFO("Frame-{} Calculation Done", frameId);
        Eigen::MatrixXd newDispVector = m_solVector;

        Eigen::MatrixXd tempVelVector = velVector;
        velVector = alpha4(Beta, Gamma, tStep) * (newDispVector - dispVector) + alpha5(Beta, Gamma) * tempVelVector + alpha6(Beta, Gamma, tStep) * accVector;
        accVector = alpha1(Beta, tStep) * (newDispVector - dispVector) - alpha2(Beta, tStep) * tempVelVector - alpha3(Beta) * accVector;
        dispVector = newDispVector;


        Eigen::MatrixXd UU = dispVector;
        {
            std::ofstream out(output_file,std::ofstream::out | std::ofstream::trunc);
            for(int p_id = 0; p_id < patch.size(); ++p_id){
                for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
                    int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

                    out << UU(global_id, 0) <<" " << UU(global_id + m_spline.m_sparse_system.getControlNums(), 0)
                        <<" " << UU(global_id+ 2 * m_spline.m_sparse_system.getControlNums(), 0) <<"\n";
                }
            }
            out.close();
        }

        Eigen::MatrixXd displacement1;
        stiff.constructCoefs(UU, displacement1);
        int begin = 0;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
            begin += num;
        }

        std::string output_path = std::string(IGA_OUTPUT_DIR) + "ResultsGPU/GPU_" + to_string(frameId) + "_";

        stiff.WriteParaview(patch, output_path, UU,
                            1000, false);
        begin = 0;
        displacement1 *= -1;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
            begin += num;
        }

        uint32_t patch_size = patch.size();

        std::vector<std::vector<CPoint>> uu(patch_size, std::vector<CPoint>(64));

        for(uint32_t p_id = 0; p_id < patch_size; ++p_id)
        {
            for(int local_id = 0; local_id < 64; ++local_id)
            {
                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];
                uu[p_id][local_id].setX(UU.coeff(global_id + 0 * m_spline.m_sparse_system.getControlNums()));
                uu[p_id][local_id].setY(UU.coeff(global_id + 1 * m_spline.m_sparse_system.getControlNums()));
                uu[p_id][local_id].setZ(UU.coeff(global_id + 2 * m_spline.m_sparse_system.getControlNums()));
            }
        }

        std::vector<double>& m_bdycoef_x = analysis_data[0];
        std::vector<double>& m_bdycoef_y = analysis_data[1];
        std::vector<double>& m_bdycoef_z = analysis_data[2];

        m_bdycoef_x.resize(mesh->vertices_.size());
        m_bdycoef_y.resize(mesh->vertices_.size());
        m_bdycoef_z.resize(mesh->vertices_.size());


        std::vector<int> mm = { 51, 3, 15, 63, 48, 0, 12, 60 };
        for (int i = 0; i < mesh->cells_.size(); ++i)
        {
            for (int j = 0; j < 8; ++j)
            {
                auto vh = mesh->cells_[i].vertex_[j];

                /*m_bdycoef_x[vh] = m_bdycoef[i][3 * j];
                m_bdycoef_y[vh] = m_bdycoef[i][3 * j + 1];
                m_bdycoef_z[vh] = m_bdycoef[i][3 * j + 2];*/
                int jj = mm[j];
                m_bdycoef_x[vh] = fabs(uu[i][jj].getX());
                m_bdycoef_y[vh] = fabs(uu[i][jj].getY());
                m_bdycoef_z[vh] = fabs(uu[i][jj].getZ());
                if (i == 0) {
                    std::cout << "vh = " << vh << std::endl;
                    std::cout << "solution_z[0] = " << m_bdycoef_z[vh] << std::endl;
                }
            }
        }

        is_need_update_model = true;
    }

    UU = dispVector;

    // 输出
    {

        std::ofstream out(output_file,std::ofstream::out | std::ofstream::trunc);
        for(int p_id = 0; p_id < patch.size(); ++p_id){
            for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

                out << UU.coeff(global_id, 0) <<" " << UU.coeff(global_id + m_spline.m_sparse_system.getControlNums(), 0)
                    <<" " << UU.coeff(global_id+ 2 * m_spline.m_sparse_system.getControlNums(), 0) <<"\n";
            }
        }
        out.close();
    }
}
