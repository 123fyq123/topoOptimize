////
//// Created by XieJin on 10/19/2023.
////
//
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
//#include "Poster/source/Timer.h"
//#include "omp.h"
//#include "Poster/source/nanoflann.h"
//
//int dim = 3;
//
//namespace{
//    using namespace nanoflann;
//    struct Point{
//        double x, y, z;
//    };
//    template <typename T>
//    struct PointCloud
//    {
//
//        using coord_t = T;  //!< The type of each coordinate
//
//        std::vector<Point> pts;
//
//        // Must return the number of data points
//        inline size_t kdtree_get_point_count() const { return pts.size(); }
//
//        // Returns the dim'th component of the idx'th point in the class:
//        // Since this is inlined and the "dim" argument is typically an immediate
//        // value, the
//        //  "if/else's" are actually solved at compile time.
//        inline T kdtree_get_pt(const size_t idx, const size_t dim) const
//        {
//            if (dim == 0)
//                return pts[idx].x;
//            else if (dim == 1)
//                return pts[idx].y;
//            else
//                return pts[idx].z;
//        }
//
//        // Optional bounding-box computation: return false to default to a standard
//        // bbox computation loop.
//        //   Return true if the BBOX was already computed by the class and returned
//        //   in "bb" so it can be avoided to redo it again. Look at bb.size() to
//        //   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
//        template <class BBOX>
//        bool kdtree_get_bbox(BBOX& /* bb */) const
//        {
//            return false;
//        }
//    };
//
//    typedef KDTreeSingleIndexAdaptor<
//            L2_Simple_Adaptor<double, PointCloud<double>>,
//            PointCloud<double>,
//            3 /* dim */
//    > KDTree;
//
//    void GetNearestIds(std::vector<Point>& p, std::vector<Point>& q, std::vector<int>& ans){
//        ans.resize(p.size());
//        PointCloud<double> A,B;
//        A.pts = p;
//        B.pts = q;
//        KDTree index(3, B, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
//        index.buildIndex();
//        int a_index = 0;
//        for (const auto& a : A.pts) {
//            size_t ret_index;
//            double out_dist_sqr;
//            nanoflann::KNNResultSet<double> resultSet(1);
//            resultSet.init(&ret_index, &out_dist_sqr);
//            index.findNeighbors(resultSet, &a.x);
//            ans[a_index++] = ret_index;
//        }
//    }
//
//}
//
//
//
//void Cuda_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
//                                const std::string& output_file){
//    Timer timer;
//    Eigen::setNbThreads(24);
//    Eigen::initParallel();
//
//    std::vector<CBSplineSolid> patch;
//    bool bH = true;
//
//
//    ToolFunc::readFileXml(input_model_file, patch);
//    timer.out("read_xml");
//
//    bool bUseBoundaryPoint = true;
//
//    std::vector<CPoint> bdy_points;
//    std::vector<int> bdy_type;
//    std::vector<CVector> bdy_value;
//
//    // add bdy info
//    {
//        int control_size = 0;
//        std::ifstream in(boundary_file);
//        in >> control_size;
//        std::cout <<"control_size = " << control_size <<"\n";
//
//        for(int i = 0; i < control_size; ++i) {
//            double x, y, z;
//            int cur_type;
//            double vx, vy, vz;
//            in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
//            CPoint cp(x, y, z);
//            CVector cv(vx, vy, vz);
//            bdy_points.push_back(cp);
//            bdy_value.push_back(cv);
//            bdy_type.push_back(cur_type);
//        }
//    }
//    timer.out("add bdy info");
//    gpmesh::GPLog::init();
//
//    CUDA_ERROR(cudaDeviceReset());
//    std::vector<boundary_point> neumann_points;
//    std::vector<boundary_point> dirichlet_points;
//    timer.out("before sys");
//    SparseSystem temp_sys(patch, dirichlet_points, 3);
//    timer.out("after sys");
//    // 处理边界条件
//    {
//
//
//        std::vector<Point> A, B;
//        std::vector<int> global_ids;
//
//        for(auto& p : bdy_points){
//            A.emplace_back(Point{p.getX(), p.getY(), p.getZ()});
//        }
//        for(int p_id = 0; p_id < patch.size(); ++p_id) {
//            for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
//                auto& pp = patch[p_id].coefs()[local_id];
//                B.emplace_back(Point{pp.getX(), pp.getY(), pp.getZ()});
//            }
//        }
//
//        GetNearestIds(A, B, global_ids);
//        timer.out("after GetNearestIds");
//        std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
//        auto& global_local_map = temp_sys.get_l2g_map();
//
//        for(int i = 0; i < bdy_points.size(); ++i){
//            auto& p = bdy_points[i];
//            int patch_id = global_ids[i] / 64;
//            int local_id = global_ids[i] % 64;
//            int global_id = global_local_map[patch_id][local_id];
//
//            if(global_id == -1){
//                std::cerr <<"cant find control points\n";
//                //std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
//                std::cout << p.getX() <<" " <<p.getY() <<" " << p.getZ() <<"\n";
//            }else{
//                if(bdy_type[i] == 2){
//                    neumann_points.emplace_back(boundary_point(condition_type::neumann, bdy_value[i],global_id));
//                }else if(bdy_type[i] == 1){
//                    dirichlet_points.emplace_back(boundary_point(condition_type::dirichlet, CVector(0, 0, 0) ,global_id));
//                }
//            }
//        }
//        timer.out("after bdy_points");
//    }
//
//    timer.out("deal with bdy info");
//
//    gpmesh::GPSplineSolidV0 m_spline;
//    Eigen::MatrixXd UU;
//    m_spline.debug_IGA_boundary_point_test_v0(patch, neumann_points, dirichlet_points, UU);
//    timer.out("solver");
//    // 输出
//    {
//        std::ofstream out(output_file,std::ofstream::out | std::ofstream::trunc);
//        for(int p_id = 0; p_id < patch.size(); ++p_id){
//            for(int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id){
//                int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];
//
//                out << UU(global_id, 0) <<" " << UU(global_id + m_spline.m_sparse_system.getControlNums(), 0)
//                <<" " << UU(global_id+ 2 * m_spline.m_sparse_system.getControlNums(), 0) <<"\n";
//            }
//        }
//        out.close();
//    }
//}
