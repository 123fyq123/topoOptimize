#include "CPU_Hyperelasticity_Solver.h"

#include "IGA_helper_function.h"
#include "Util/GP_Timer.h"

#include "Poster/source/Timer.h"
#include <Hyperelasticity/TimeIntegrator.h>
#include <Hyperelasticity/collision_detect.h>
#include "GP_nsp.h"
#include "Eigen/PardisoSupport"

void CPU_Hyperelasticity_TimeIntegrator_Solver(const std::string& input_model_file,
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

    uint32_t patch_size = patch.size();

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
    SparseSystem sparse_system;
    sparse_system.init_geometry(patch, 3);

    timer.out("Hyperelasticity GPU Geometry Init");

    timer.out("before sys");
    const SparseSystem& temp_sys = sparse_system;
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
    timer.out("deal with bdy info");

    sparse_system.init_dirichlet(dirichlet_points);

    stiffAssembler stiff(patch, dirichlet_points, neumann_points);

    stiff.setYangsPossion(Yangs, Poissons);

    double density = 3320;
    massAssembler mass(patch, density, stiff.system());

    int MaxFrame = frameCount;

    double tStep = timeStep;
    double Beta = 0.25;
    double Gamma = 0.5;

    timeIntegrator solver(stiff, mass);

    solver.setInitialDisplacement(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
    solver.setInitialVelocity(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
    solver.initialize();

    if(solver.bUsePoint)
    {
        solver.assemble_points(solver.displacementVector(), true);
    }
    else
    {
        solver.assemble(solver.displacementVector(), true);
    }

    Eigen::MatrixXd displacement1;

    Eigen::MatrixXd UU = Eigen::MatrixXd::Zero(stiff.system().getControlNums() * stiff.system().getDim(), 1);

    std::string origin_output_path = std::string(IGA_OUTPUT_DIR) + "Origin/Origin_";

    stiff.WriteParaview(patch, origin_output_path, UU,
                        1000, false);

    for (int frameId = 0; frameId < MaxFrame; ++frameId) {
        std::cout << frameId << "\n";
        solver.makeTimeStep(tStep);
        //std::cout << .topRows(10);
        stiff.constructCoefs(solver.displacementVector(), displacement1);
        //stiff.constructCoefs(u, displacement1);
        UU = solver.displacementVector();
        int begin = 0;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
            begin += num;
        }
        //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/xiao/threeblock2_disp");
        //ToolFunc::outputXml(patch, "D:/Resource/FileDataInput/simulate/C_deal_fail");
        UU = solver.displacementVector();

        std::string output_path = std::string(IGA_OUTPUT_DIR) + "Test/test_" + to_string(frameId) + "_";

        stiff.WriteParaview(patch, output_path, UU,
                            1000, false);
        begin = 0;
        displacement1 *= -1;
        for (size_t i = 0; i < patch.size(); ++i) {
            int num = patch[i].getSizeControlPoints();
            patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
            begin += num;
        }

        std::vector<std::vector<CPoint>> uu(patch_size, std::vector<CPoint>(64));

        for(uint32_t p_id = 0; p_id < patch_size; ++p_id)
        {
            for(int local_id = 0; local_id < 64; ++local_id)
            {
                int global_id = stiff.system().get_l2g_map()[p_id][local_id];
                uu[p_id][local_id].setX(UU.coeff(global_id + 0 * stiff.system().getControlNums()));
                uu[p_id][local_id].setY(UU.coeff(global_id + 1 * stiff.system().getControlNums()));
                uu[p_id][local_id].setZ(UU.coeff(global_id + 2 * stiff.system().getControlNums()));
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

}