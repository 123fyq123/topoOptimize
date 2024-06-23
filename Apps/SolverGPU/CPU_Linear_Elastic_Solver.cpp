#include "CPU_Linear_Elastic_Solver.h"
#include "GP_Spline_CPU.h"
#include "Poster/source/Timer.h"
#include <Hyperelasticity/TimeIntegrator.h>
#include <Hyperelasticity/collision_detect.h>
#include "GP_nsp.h"
#include <Eigen/PardisoSupport>

void CPU_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
                                const std::string& output_file, double Yangs, double Poissons)
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
    timer.out("before sys");
    SparseSystem temp_sys(patch, dirichlet_points, 3);
    timer.out("after sys");
    // 处理边界条件
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

    gpmesh::GPSplineSolidCPU<double> m_spline;
    Eigen::MatrixXd UU;
    m_spline.linear_elasticity_pcg_solving(patch,  dirichlet_points, neumann_points, UU, Yangs, Poissons);
    timer.out("solver");
    // 输出
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
}

void CPU_Linear_Elastic_Solver_Paridiso(const std::string& input_model_file, const std::string& boundary_file,
                                        const std::string& output_file, double Yangs, double Poissons)
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
    timer.out("before sys");
    SparseSystem temp_sys(patch, dirichlet_points, 3);
    timer.out("after sys");
    // 处理边界条件
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

    gpmesh::GPSplineSolidCPU<double> m_spline;
    Eigen::MatrixXd UU;
    m_spline.linear_elasticity_pardiso(patch,  dirichlet_points, neumann_points, UU, Yangs, Poissons);
    timer.out("solver");
    // 输出
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
}

void CPU_Linear_Elastic_Solver_Eigen(const std::string& input_model_file, const std::string& boundary_file,
                                     const std::string& output_file, double Yangs, double Poissons, bool bUseCG)
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
    timer.out("before sys");
    SparseSystem temp_sys(patch, dirichlet_points, 3);
    timer.out("after sys");
    // 处理边界条件
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

    gpmesh::GPSplineSolidCPU<double> m_spline;
    Eigen::MatrixXd UU;
    m_spline.linear_elasticity_Eigen(patch,  dirichlet_points, neumann_points, UU, Yangs, Poissons, bUseCG);
    timer.out("solver");
    // 输出
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
}

void CPU_Linear_Elastic_Solver_WN(const std::string& input_model_file, const std::string& boundary_file,
                                     const std::string& output_file, double Yangs, double Poissons, bool bUseCG)
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
    timer.out("before sys");
    SparseSystem temp_sys(patch, dirichlet_points, 3);
    timer.out("after sys");
    // 处理边界条件
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

    stiffAssembler stiff(patch, dirichlet_points, neumann_points);

    stiff.setYangsPossion(Yangs, Poissons);
    auto initDisp = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
    stiff.assemble_points();

    Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> Solver (stiff.getMatrix());

    auto res = Solver.solve(stiff.getRhs());

    timer.out("solver");
    // 输出
    Eigen::MatrixXd UU;
    UU.setZero(stiff.system().getControlNums() * 3, 1);

    {
        std::ofstream out(output_file,std::ofstream::out | std::ofstream::trunc);
        for(int p_id = 0; p_id < patch.size(); ++p_id){
            for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
                int global_id = stiff.system().get_l2g_map()[p_id][local_id];

                out << UU(global_id, 0) <<" " << UU(global_id + stiff.system().getControlNums(), 0)
                    <<" " << UU(global_id+ 2 * stiff.system().getControlNums(), 0) <<"\n";
            }
        }
        out.close();
    }
}