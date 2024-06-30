//
// Created by XieJin on 10/27/2023.
//
#include "Cuda_TopOpt_Pipeline.h"
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
#include "Poster/source/Timer.h"
#include "omp.h"
#include "Poster/source/nanoflann.h"

#include <algorithm>
#include <vector>

namespace Cuda_TopOpt_Pipeline_nsp{
    using namespace nanoflann;
    struct Point{
        double x, y, z;
    };
    template <typename T>
    struct PointCloud
    {

        using coord_t = T;  //!< The type of each coordinate

        std::vector<Point> pts;

        // Must return the number of data points
        inline size_t kdtree_get_point_count() const { return pts.size(); }

        // Returns the dim'th component of the idx'th point in the class:
        // Since this is inlined and the "dim" argument is typically an immediate
        // value, the
        //  "if/else's" are actually solved at compile time.
        inline T kdtree_get_pt(const size_t idx, const size_t dim) const
        {
            if (dim == 0)
                return pts[idx].x;
            else if (dim == 1)
                return pts[idx].y;
            else
                return pts[idx].z;
        }

        // Optional bounding-box computation: return false to default to a standard
        // bbox computation loop.
        //   Return true if the BBOX was already computed by the class and returned
        //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
        //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
        template <class BBOX>
        bool kdtree_get_bbox(BBOX& /* bb */) const
        {
            return false;
        }
    };

    typedef KDTreeSingleIndexAdaptor<
            L2_Simple_Adaptor<double, PointCloud<double>>,
            PointCloud<double>,
            3 /* dim */
    > KDTree;

    void GetNearestIds(std::vector<Point>& p, std::vector<Point>& q, std::vector<int>& ans){
        ans.resize(p.size());
        PointCloud<double> A,B;
        A.pts = p;
        B.pts = q;
        KDTree index(3, B, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
        index.buildIndex();
        int a_index = 0;
        for (const auto& a : A.pts) {
            size_t ret_index;
            double out_dist_sqr;
            nanoflann::KNNResultSet<double> resultSet(1);
            resultSet.init(&ret_index, &out_dist_sqr);
            index.findNeighbors(resultSet, &a.x);
            ans[a_index++] = ret_index;
        }
    }

}


void GPU_OC(std::vector<double>& Density, const double TargetRate, const std::vector<double>& dc, const std::vector<int>& isPicked) {
    double a = 0.0;
    double b = 10000000000;
    double MoveStep = 0.2;

    std::vector<double> TempDensity = Density;

    int max_iter = 50;
    while (b - a > 1e-10 && max_iter-- > 0) {
        std::cout << "[OC]: " << "(a,b) = " << a << "," << b << "\n";
        double Lmid = 0.5 * (a + b);
        double SumDensity = 0.0;
        for (int i = 0; i < Density.size(); ++i) {
            if (isPicked[i] == 1) TempDensity[i] = 1, std::cout <<"111\n" << std::endl;
            else TempDensity[i] = fmax(0.0001, fmax(Density[i] - MoveStep, fmin(1.0, fmin(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
            //TempDensity[i] = fmax(0.001, fmax(Density[i] - MoveStep, fmin(1.0, fmin(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
            SumDensity += TempDensity[i];
        }
        if (SumDensity - TargetRate * Density.size() > 0) {
            a = Lmid;
        }
        else {
            b = Lmid;
        }
    }
    double SumDensityOutput = 0.0;
    for (int i = 0; i < Density.size(); ++i) {
        SumDensityOutput += TempDensity[i];
    }
    std::cout << "SumDensityOutput = " << SumDensityOutput / Density.size() << std::endl;
    std::cout <<"TargetRate = " << TargetRate << std::endl;
    Density = TempDensity;
}

void Cuda_TopOpt_Pipeline(const std::string& input_model_file,
                          const std::string& bdy_file,
                          const std::vector<int>& protect_cell_ids,
                          const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
                          double target_rate,
                          double max_iter,
                          double Yangs,
                          double Poissons,
                          std::vector<double>& density,
                          std::atomic_bool& is_need_update_model){
    Timer timer;
    Eigen::setNbThreads(24);
    Eigen::initParallel();

    std::vector<CBSplineSolid> patch;

    ToolFunc::readFileXml(input_model_file, patch);
    timer.out("read_xml");

    std::vector<CPoint> bdy_points;
    std::vector<int> bdy_type;
    std::vector<CVector> bdy_value;

    // add bdy info
    {
        int control_size = 0;
        std::ifstream in(bdy_file);
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

    timer.out("deal with bdy info");
    std::vector<double> dc_pre(patch.size(), 0);
    auto exist_num = (double)patch.size();
    auto tot_num = (double)patch.size();
    //while(1){
    // data init
    gpmesh::GPSplineSolidV0 m_spline;
//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
//    timer.out("TopOpt_GPU_Init");
    m_spline.init_patch(patch);
    timer.out("TopOpt_GPU_Geometry_Init");

    timer.out("before sys");
    const SparseSystem& temp_sys = m_spline.m_sparse_system;
    timer.out("after sys");
    {
        std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
        std::vector<int> global_ids;

        for(auto& p : bdy_points){
            A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{p.getX(), p.getY(), p.getZ()});
        }
        for(int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                auto& pp = patch[p_id].coefs()[local_id];
                B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{pp.getX(), pp.getY(), pp.getZ()});
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
                    neumann_points.emplace_back(condition_type::neumann, bdy_value[i],global_id);
                }else if(bdy_type[i] == 1){
                    dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0) ,global_id);
                }
            }
        }
    }

    m_spline.init_dirichlet_points(dirichlet_points);

    std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
    // first assembly
    m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
    timer.out("Get_Local_Matrix");
        m_spline.m_sparse_system.rhs().setZero();

        m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
        m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


    int topopt_iter = 0;

    std::vector<int> isPicked(patch.size(), 0);
    for(auto& xx : protect_cell_ids){
        isPicked[xx] = 1;
    }

    while(topopt_iter < max_iter){

        std::cout << topopt_iter++ << std::endl;
        int reduce_num = (int) (0.02 * (exist_num - tot_num * target_rate));
        if(reduce_num < 1) break;
        Eigen::MatrixXd U;
        /// todo: need to check


        m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

        Eigen::MatrixXd UU(patch.size(), 64 * 3);

        for(int p_id = 0; p_id < patch.size(); ++p_id){
            for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

                UU(p_id, local_id) = U(global_id, 0);
                UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
                UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
            }
        }

        std::vector<double> Ce(patch.size());
        std::vector<double> dc(patch.size());
        for(int p_id = 0; p_id < patch.size(); ++p_id){
            Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
            if(Ce[p_id] < 0) std::cout <<"Ce is wrong!\n";
            dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
        }


        // fliter
        {
            auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
            int n = patch.size();                        //�ܵ����ɶ���Ŀ
            for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
                double weight = 0;
                double temp_dx = 0;
                int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
                for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
                    weight += neighbor_info[i][j].second;
                    temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
                }
                temp_dx /= weight;
                Tempdc[i] = temp_dx;          //temp_dc��ֵ
            }
            dc = Tempdc;       //������Ľ������
        }

//        {
//            for(int i = 0; i < tot_num; ++i){
//                dc[i] = (dc[i] + dc_pre[i]) / 2;
//                dc_pre[i] = dc[i];
//            }
//        }

        GPU_OC(density, target_rate, dc, isPicked);

//        std::vector<std::pair<double, int>> dc_ids;
//        for(int i = 0; i < patch.size(); ++i){
//            dc_ids.emplace_back(-dc[i], i);
//            //std::cout <<density[i] << " " <<dc[i] <<" ";
//        }
//
//        std::sort(dc_ids.begin(), dc_ids.end());
//
//        int remove_num = 0;
//
//        int kk = -1;
//        while(remove_num < reduce_num && kk < (int)dc_ids.size()){
//            kk++;
//            int id = dc_ids[kk].second;
//           //std::cout <<"dc_ids_id = " << id << std::endl;
//            if(density[id] < 0.8) continue;
//            //std::cout <<"NOOOOOOOOOOOOOOO\n";
//            if(std::find(protect_cell_ids.begin(), protect_cell_ids.end(), id) != protect_cell_ids.end()) continue;
//           //std::cout <<"YESSSSSSSSSSSSSSSSSSS\n";
//            density[id] = 0.0001;
//            remove_num++;
//            exist_num--;
//        }
//        std::cout << "reduce_num = " << reduce_num << std::endl;
//        std::cout << exist_num << std::endl;
        is_need_update_model = true;
    }
}


void Cuda_TopOpt_Discrete(const std::string& input_model_file,
    const std::string& bdy_file,
    const std::vector<int>& protect_cell_ids,
    const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
    double target_rate,
    double max_iter,
    double Yangs,
    double Poissons,
    std::vector<double>& density,
    std::atomic_bool& is_need_update_model) {
    Timer timer;
    Eigen::setNbThreads(24);
    Eigen::initParallel();

    std::vector<CBSplineSolid> patch;

    ToolFunc::readFileXml(input_model_file, patch);
    timer.out("read_xml");

    std::vector<CPoint> bdy_points;
    std::vector<int> bdy_type;
    std::vector<CVector> bdy_value;

    // add bdy info
    {
        int control_size = 0;
        std::ifstream in(bdy_file);
        in >> control_size;
        std::cout << "control_size = " << control_size << "\n";

        for (int i = 0; i < control_size; ++i) {
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

    timer.out("deal with bdy info");
    std::vector<double> dc_pre(patch.size(), 0);
    auto exist_num = (double)patch.size();
    auto tot_num = (double)patch.size();
    //while(1){
    // data init
    gpmesh::GPSplineSolidV0 m_spline;
    //    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
    //    timer.out("TopOpt_GPU_Init");
    m_spline.init_patch(patch);
    timer.out("TopOpt_GPU_Geometry_Init");

    timer.out("before sys");
    const SparseSystem& temp_sys = m_spline.m_sparse_system;
    timer.out("after sys");
    {
        std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
        std::vector<int> global_ids;

        for (auto& p : bdy_points) {
            A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
        }
        for (int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                auto& pp = patch[p_id].coefs()[local_id];
                B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
            }
        }

        GetNearestIds(A, B, global_ids);
        timer.out("after GetNearestIds");
        std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
        auto& global_local_map = temp_sys.get_l2g_map();

        for (int i = 0; i < bdy_points.size(); ++i) {
            auto& p = bdy_points[i];
            int patch_id = global_ids[i] / 64;
            int local_id = global_ids[i] % 64;
            int global_id = global_local_map[patch_id][local_id];

            if (global_id == -1) {
                std::cerr << "cant find control points\n";
                //std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
                std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
            }
            else {
                if (bdy_type[i] == 2) {
                    neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
                }
                else if (bdy_type[i] == 1) {
                    dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
                }
            }
        }
    }

    m_spline.init_dirichlet_points(dirichlet_points);

    std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
    // first assembly
    m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
    timer.out("Get_Local_Matrix");
    m_spline.m_sparse_system.rhs().setZero();

    m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
    m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


    int topopt_iter = 0;

    std::vector<int> isPicked(patch.size(), 0);
    for (auto& xx : protect_cell_ids) {
        isPicked[xx] = 1;
    }

    Eigen::MatrixXd U;
    m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

    Eigen::MatrixXd UU(patch.size(), 64 * 3);

    for (int p_id = 0; p_id < patch.size(); ++p_id) {
        for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
            int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

            UU(p_id, local_id) = U(global_id, 0);
            UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
            UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
        }
    }

    std::vector<double> Ce(patch.size());
    std::vector<double> dc(patch.size());
    for (int p_id = 0; p_id < patch.size(); ++p_id) {
        Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
        if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
        dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
    }


    // fliter
    {
        auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
        int n = patch.size();                        //�ܵ����ɶ���Ŀ
        for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
            double weight = 0;
            double temp_dx = 0;
            int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
            for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
                weight += neighbor_info[i][j].second;
                temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
            }
            temp_dx /= weight;
            Tempdc[i] = temp_dx;          //temp_dc��ֵ
        }
        dc = Tempdc;       //������Ľ������
    }

    double maxdc = fabs(dc[0]);
    for (int i = 1; i < dc.size(); i++) {
        if (fabs(dc[i]) > maxdc) {
            maxdc = fabs(dc[i]);
        }
    }
    for (int i = 0; i < dc.size(); i++) {
        dc[i] /= maxdc;
    }

    double c_d;
    double c_old = 0, c_new = 0, sumnew = 0;
    c_d = density.size();

    double beta, v;
    std::vector<double> sigma;
    std::vector<double> olddensity = density;
    std::vector<double> newdensity = density;

    sigma.resize(density.size(), 0);
    beta = 2000;
    v = 1;
    double lambda = 1e-3;

    double w1 = 1e-8, Vr;
    Vr = density.size() * 0.98;

    bool flag = false;
    int pareto = 0;

    while (topopt_iter < max_iter) {

        std::cout << topopt_iter++ << std::endl;
        int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
        if (reduce_num < 1) break;

        /// todo: need to check

        olddensity = density;



        while (true) {
            double suma = 0, sumb = 0;
            std::vector<double> s(density.size());
            for (int i = 0; i < density.size(); ++i) {
                s[i] = lambda + dc[i];
            }

            for (int i = 0; i < s.size(); ++i) {
                if (fabs(s[i]) <= 1E-8) {
                    s[i] = 1E-8 * (s[i] >= 0 ? 1 : -1);
                }
            }
            // 计算 delta
            std::complex<double> delta = pow(beta, 2) / 27;
            // 计算 r
            std::vector<std::complex<double>> r(s.size());
            for (int i = 0; i < s.size(); ++i) {
                std::complex<double> sqrt_term = sqrt(std::complex<double>(pow(s[i], 2) * (pow(s[i], 2) - delta)));
                r[i] = pow(delta, -1.0 / 3) * pow(std::complex<double>(2 * pow(s[i], 2) - delta + sqrt_term + sqrt_term), 1.0 / 3);
            }
            // 计算 sigma
            for (int i = 0; i < sigma.size(); ++i) {
                sigma[i] = (1.0 / 6) * beta * (-1 + r[i].real() * 2);
                suma += v - dc[i] / sigma[i];
                sumb += 1 / sigma[i];
            }

            lambda = (suma - 2 * Vr) / sumb;

            c_old = 0;
            c_new = 0;
            sumnew = 0;
            for (int i = 0; i < density.size(); i++) {
                if (isPicked[i] == 1) {
                    newdensity[i] = 1;
                }
                else {
                    newdensity[i] = 0.5 * (1.0 - (lambda + dc[i]) / sigma[i]);
                }
                
                sumnew += newdensity[i];
                c_old += olddensity[i] * dc[i];
                c_new += newdensity[i] * dc[i];
            }
            if (fabs(c_old - c_new) / fabs(c_old) <= w1) {
                break;
            }
            else {
                olddensity = newdensity;
            }
        }

        for (int i = 0; i < newdensity.size(); i++) {
            if (newdensity[i] <= 0) {
                newdensity[i] = 0.001;
            }
            else if (newdensity[i] >= 1) {
                newdensity[i] = 1;
            }

            newdensity[i] = round(newdensity[i]);
        }

        //density = newdensity;
        //Vr *= 0.9;
        m_spline.Solver_With_Density(neumann_points, dirichlet_points, newdensity, U);

        for (int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

                UU(p_id, local_id) = U(global_id, 0);
                UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
                UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
            }
        }
        for (int p_id = 0; p_id < patch.size(); ++p_id) {
            Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
            if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
            dc[p_id] = (-3 * newdensity[p_id] * newdensity[p_id] * (1 - 0.0001)) * Ce[p_id];
        }

        // fliter
        {
            auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
            int n = patch.size();                        //�ܵ����ɶ���Ŀ
            for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
                double weight = 0;
                double temp_dx = 0;
                int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
                for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
                    weight += neighbor_info[i][j].second;
                    temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
                }
                temp_dx /= weight;
                Tempdc[i] = temp_dx;          //temp_dc��ֵ
            }
            dc = Tempdc;       //������Ľ������
        }

        double maxdc = fabs(dc[0]);
        for (int i = 1; i < dc.size(); i++) {
            if (fabs(dc[i]) > maxdc) {
                maxdc = fabs(dc[i]);
            }
        }
        for (int i = 0; i < dc.size(); i++) {
            dc[i] /= maxdc;
        }

        if (fabs(c_new - c_d) / fabs(c_d) > 5e-3) {
            density = newdensity;
            c_d = c_new;
        }
        else {
            if (flag) {
                break;
            }
            density = newdensity;
            if (Vr / density.size() <= target_rate) {
                Vr = density.size() * max(0.3, Vr / density.size());
                flag = true;
                continue;
            }
            else {
                Vr *= 0.98;
            }

        }

        is_need_update_model = true;
    }
    std::cout << "Finsh" << std::endl;
}


void Cuda_TopOpt_Multi_Pipeline(const std::string& input_model_file,
                          const std::string& bdy_file,
                          const std::vector<int>& protect_cell_ids,
                          const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
                          int sub_time,
                          double target_rate,
                          double max_iter,
                          double Yangs,
                          double Poissons,
                          std::vector<double>& density,
                          std::atomic_bool& is_need_update_model){
    Timer timer;
    Eigen::setNbThreads(24);
    Eigen::initParallel();

    std::vector<CBSplineSolid> patch;

    ToolFunc::readFileXml(input_model_file, patch);
    timer.out("read_xml");

    std::vector<CPoint> bdy_points;
    std::vector<int> bdy_type;
    std::vector<CVector> bdy_value;

    // add bdy info
    {
        int control_size = 0;
        std::ifstream in(bdy_file);
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

    timer.out("deal with bdy info");
    std::vector<double> dc_pre(patch.size(), 0);
    auto exist_num = (double)patch.size();
    auto tot_num = (double)patch.size();

    {
        assert(sub_time < 5);
        int sub_size = 1 << sub_time;
        std::vector<std::vector<Eigen::MatrixXd>> sub_local_matrix(sub_size, std::vector<Eigen::MatrixXd>(patch.size(), Eigen::MatrixXd::Zero(192, 192)));

    }





    //while(1){
    // data init
    gpmesh::GPSplineSolidV0 m_spline;
//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
//    timer.out("TopOpt_GPU_Init");
    m_spline.init_patch(patch);
    timer.out("TopOpt_GPU_Geometry_Init");

    timer.out("before sys");
    const SparseSystem& temp_sys = m_spline.m_sparse_system;
    timer.out("after sys");
    {
        std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
        std::vector<int> global_ids;

        for(auto& p : bdy_points){
            A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{p.getX(), p.getY(), p.getZ()});
        }
        for(int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                auto& pp = patch[p_id].coefs()[local_id];
                B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{pp.getX(), pp.getY(), pp.getZ()});
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
                    neumann_points.emplace_back(condition_type::neumann, bdy_value[i],global_id);
                }else if(bdy_type[i] == 1){
                    dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0) ,global_id);
                }
            }
        }
    }

    m_spline.init_dirichlet_points(dirichlet_points);

    std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
    // first assembly
    m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
    timer.out("Get_Local_Matrix");
    m_spline.m_sparse_system.rhs().setZero();

    m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
    m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


    int topopt_iter = 0;

    std::vector<int> isPicked(patch.size(), 0);
    for(auto& xx : protect_cell_ids){
        isPicked[xx] = 1;
    }

    while(topopt_iter < max_iter){

        std::cout << topopt_iter++ << std::endl;
        int reduce_num = (int) (0.02 * (exist_num - tot_num * target_rate));
        if(reduce_num < 1) break;
        Eigen::MatrixXd U;
        /// todo: need to check


        m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

        Eigen::MatrixXd UU(patch.size(), 64 * 3);

        for(int p_id = 0; p_id < patch.size(); ++p_id){
            for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

                UU(p_id, local_id) = U(global_id, 0);
                UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
                UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
            }
        }

        std::vector<double> Ce(patch.size());
        std::vector<double> dc(patch.size());
        for(int p_id = 0; p_id < patch.size(); ++p_id){
            Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
            if(Ce[p_id] < 0) std::cout <<"Ce is wrong!\n";
            dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
        }


        // fliter
        {
            auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
            int n = patch.size();                        //�ܵ����ɶ���Ŀ
            for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
                double weight = 0;
                double temp_dx = 0;
                int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
                for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
                    weight += neighbor_info[i][j].second;
                    temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
                }
                temp_dx /= weight;
                Tempdc[i] = temp_dx;          //temp_dc��ֵ
            }
            dc = Tempdc;       //������Ľ������
        }

//        {
//            for(int i = 0; i < tot_num; ++i){
//                dc[i] = (dc[i] + dc_pre[i]) / 2;
//                dc_pre[i] = dc[i];
//            }
//        }

        GPU_OC(density, target_rate, dc, isPicked);

//        std::vector<std::pair<double, int>> dc_ids;
//        for(int i = 0; i < patch.size(); ++i){
//            dc_ids.emplace_back(-dc[i], i);
//            //std::cout <<density[i] << " " <<dc[i] <<" ";
//        }
//
//        std::sort(dc_ids.begin(), dc_ids.end());
//
//        int remove_num = 0;
//
//        int kk = -1;
//        while(remove_num < reduce_num && kk < (int)dc_ids.size()){
//            kk++;
//            int id = dc_ids[kk].second;
//           //std::cout <<"dc_ids_id = " << id << std::endl;
//            if(density[id] < 0.8) continue;
//            //std::cout <<"NOOOOOOOOOOOOOOO\n";
//            if(std::find(protect_cell_ids.begin(), protect_cell_ids.end(), id) != protect_cell_ids.end()) continue;
//           //std::cout <<"YESSSSSSSSSSSSSSSSSSS\n";
//            density[id] = 0.0001;
//            remove_num++;
//            exist_num--;
//        }
//        std::cout << "reduce_num = " << reduce_num << std::endl;
//        std::cout << exist_num << std::endl;
        is_need_update_model = true;
    }
}


void Cuda_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
                                const std::string& output_file, double Yangs, double Poissons){
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
        std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
        std::vector<int> global_ids;

        for(auto& p : bdy_points){
            A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{p.getX(), p.getY(), p.getZ()});
        }
        for(int p_id = 0; p_id < patch.size(); ++p_id) {
            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
                auto& pp = patch[p_id].coefs()[local_id];
                B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{pp.getX(), pp.getY(), pp.getZ()});
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
                    neumann_points.emplace_back(condition_type::neumann, bdy_value[i],global_id);
                }else if(bdy_type[i] == 1){
                    dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0) ,global_id);
                }
            }
        }
    }

    m_spline.init_dirichlet_points(dirichlet_points);

    Eigen::MatrixXd UU;
    m_spline.debug_IGA_linear_solver(patch, neumann_points, dirichlet_points, Yangs, Poissons, UU);
    timer.out("solver");

    stiffAssembler stiff(patch, dirichlet_points, neumann_points);

    Eigen::MatrixXd originUU = Eigen::MatrixXd::Zero(stiff.system().getControlNums() * stiff.system().getDim(), 1);

    std::string origin_output_path = std::string(IGA_OUTPUT_DIR) + "Origin/Origin_";

    stiff.WriteParaview(patch, origin_output_path, UU,
                        1000, false);

    std::string output_path = std::string(IGA_OUTPUT_DIR) + "Test/test_linear" + to_string(0) + "_";

    Eigen::MatrixXd displacement1;

    int begin = 0;
    stiff.constructCoefs(UU, displacement1);

    for (size_t i = 0; i < patch.size(); ++i) {
        int num = patch[i].getSizeControlPoints();
        patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
        begin += num;
    }

    stiff.WriteParaview(patch, output_path, UU,
                        1000, false);

    begin = 0;
    displacement1 *= -1;
    for (size_t i = 0; i < patch.size(); ++i) {
        int num = patch[i].getSizeControlPoints();
        patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
        begin += num;
    }

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
