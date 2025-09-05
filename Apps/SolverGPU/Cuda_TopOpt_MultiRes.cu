#include<chrono>
#include<Hyperelasticity/collision_detect.h>
#include <iostream>
#include <algorithm>
#include <vector>
#include <execution>
#include "Util/GP_Log.h"
#include "GP_Spline_Test.h"
#include "GP_Spline_v0.h"
#include "GP_BoundaryCPoint.h"
#include "Poster/source/Timer.h"
#include "Poster/source/nanoflann.h"
#include "Cuda_TopOpt_MultiRes.h"

typedef nanoflann::KDTreeSingleIndexAdaptor<nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
        PointCloud<double>,
        3 /* dim */
> KDTree;

void GetNearestIds(std::vector<MPoint>& p, std::vector<MPoint>& q, std::vector<int>& ans){
    ans.resize(p.size());
    PointCloud<double> A,B;
    A.pts = p;
    B.pts = q;
    KDTree index(3, B, nanoflann::KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
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

void GPU_OCC(std::vector<double>& Density, const double TargetRate,
             const std::vector<double>& dc,
             const std::vector<int>& isPicked,
             const std::vector<double>& volumes)
{
    double a = 0.0, b = 1e10, MoveStep = 0.1;
    std::vector<double> TempDensity = Density;

    int max_iter = 50;
    while (b - a > 1e-10 && max_iter-- > 0) {
        double Lmid = 0.5 * (a + b);
        double SumVolRho = 0.0, TotalVol = 0.0;
        for (int i = 0; i < (int)Density.size(); ++i) {
            if (isPicked[i] == 1) TempDensity[i] = 1.0;
            else {
                double rho = Density[i];
                double val = rho * std::sqrt(std::max(1e-20, -dc[i]) / Lmid);
                double low = std::max(0.0001, rho - MoveStep);
                double upp = std::min(1.0,     rho + MoveStep);
                TempDensity[i] = std::min(upp, std::max(low, val));
            }
            SumVolRho += TempDensity[i] * volumes[i];
            TotalVol  += volumes[i];
        }
        double curRate = SumVolRho / std::max(1e-20, TotalVol);
        if (curRate - TargetRate > 0) a = Lmid; else b = Lmid;
    }
    Density.swap(TempDensity);
}

void ComputeFineElementEnergies(
        //const std::vector<Eigen::MatrixXd>& coarse_local_matrix,
        const Eigen::MatrixXd& coarse_UU,
        const std::vector<CBSplineSolid>& coarse_patch,
        const double Youngs, const double Poissons,
        MultiResMapping& mapping,
        std::vector<FineElement>& fine_elements)
{
    if (!mapping.isValid() || coarse_UU.rows() != mapping.coarse_size) {
        std::cerr << "Invalid input for energy computation (mapping/UU).\n";
        return;
    }
    fine_elements.resize(mapping.fine_size);

    const double gp[4] = { -0.8611363115940526, -0.3399810435848563,
                           0.3399810435848563,  0.8611363115940526 };
    const double gw[4] = {  0.3478548451374538,  0.6521451548625461,
                            0.6521451548625461,  0.3478548451374538 };

    const int linear_subdiv = (int)std::round(std::pow(mapping.subdiv_factor, 1.0/3.0));

    const double nu = Poissons, E = Youngs;
    const double lambda = (E*nu)/((1+nu)*(1-2*nu));
    const double mu     = E/(2*(1+nu));
    Eigen::Matrix<double,6,6> D = Eigen::Matrix<double,6,6>::Zero();
    D(0,0)=D(1,1)=D(2,2)=lambda+2.0*mu;
    D(0,1)=D(0,2)=D(1,0)=D(1,2)=D(2,0)=D(2,1)=lambda;
    D(3,3)=D(4,4)=D(5,5)=mu;

    auto bern3 = [](double t, double B[4]){
        double s = 1.0 - t;
        B[0]=s*s*s; B[1]=3.0*t*s*s; B[2]=3.0*t*t*s; B[3]=t*t*t;
    };
    auto dbern3_dt = [](double t, double dB[4]){
        double s = 1.0 - t;
        dB[0] = -3.0*s*s;
        dB[1] =  3.0*s*s - 6.0*t*s;
        dB[2] =  6.0*t*s - 3.0*t*t;
        dB[3] =  3.0*t*t;
    };
    auto lid = [](int i,int j,int k)->int { return k*16 + j*4 + i; }; // i,j,k ∈ [0,3]
    std::for_each(std::execution::par_unseq, fine_elements.begin(), fine_elements.end(),
                  [](FineElement& e){ e.energy=0.0; e.volume=0.0; e.sensitivity=0.0; });

    std::vector<int> coarse_ids(mapping.coarse_size);
    std::iota(coarse_ids.begin(), coarse_ids.end(), 0);

    std::for_each(std::execution::par_unseq, coarse_ids.begin(), coarse_ids.end(),
                  [&](int c){
                      if (c >= (int)coarse_patch.size()) return;
                      if (c >= coarse_UU.rows()) return;

                      const auto& coefs = coarse_patch[c].coefs();
                      if ((int)coefs.size() != 64) return;

                      Eigen::VectorXd uc = coarse_UU.row(c).transpose();
                      if (uc.size() != 192) return;

                      std::vector<double> e_f(mapping.subdiv_factor, 0.0);
                      std::vector<double> v_f(mapping.subdiv_factor, 0.0);

                      for (int a=0;a<4;++a) for (int b=0;b<4;++b) for (int d=0;d<4;++d){
                                  const double xi   = gp[a],  wi = gw[a];
                                  const double eta  = gp[b],  wj = gw[b];
                                  const double zeta = gp[d],  wk = gw[d];
                                  const double wgt  = wi*wj*wk;
                                  const double u = 0.5*(xi   + 1.0);
                                  const double v = 0.5*(eta  + 1.0);
                                  const double w = 0.5*(zeta + 1.0);

                                  double Bu[4], Bv[4], Bw[4];
                                  double dBu_dt[4], dBv_dt[4], dBw_dt[4];
                                  bern3(u, Bu);      dbern3_dt(u, dBu_dt);
                                  bern3(v, Bv);      dbern3_dt(v, dBv_dt);
                                  bern3(w, Bw);      dbern3_dt(w, dBw_dt);

                                  double dBu_dxi[4], dBv_deta[4], dBw_dzeta[4];
                                  for (int m=0;m<4;++m) {
                                      dBu_dxi[m]   = 0.5 * dBu_dt[m];
                                      dBv_deta[m]  = 0.5 * dBv_dt[m];
                                      dBw_dzeta[m] = 0.5 * dBw_dt[m];
                                  }

                                  Eigen::Vector3d X_xi(0,0,0), X_eta(0,0,0), X_zet(0,0,0);
                                  double dN_dxi[64], dN_deta[64], dN_dzet[64];

                                  int idx=0;
                                  for (int kk=0; kk<4; ++kk)
                                      for (int jj=0; jj<4; ++jj)
                                          for (int ii=0; ii<4; ++ii, ++idx){
                                              const double Nxi  = dBu_dxi[ii]*Bv[jj]*Bw[kk];
                                              const double Neta = Bu[ii]*dBv_deta[jj]*Bw[kk];
                                              const double Nzet = Bu[ii]*Bv[jj]*dBw_dzeta[kk];

                                              dN_dxi[idx]  = Nxi;
                                              dN_deta[idx] = Neta;
                                              dN_dzet[idx] = Nzet;

                                              const int L = lid(ii,jj,kk);
                                              const double X = coefs[L].getX();
                                              const double Y = coefs[L].getY();
                                              const double Z = coefs[L].getZ();
                                              X_xi  += Eigen::Vector3d(X,Y,Z) * Nxi;
                                              X_eta += Eigen::Vector3d(X,Y,Z) * Neta;
                                              X_zet += Eigen::Vector3d(X,Y,Z) * Nzet;
                                          }
                                  Eigen::Matrix3d J;
                                  J.col(0)=X_xi; J.col(1)=X_eta; J.col(2)=X_zet;

                                  const double detJ = J.determinant();
                                  if (!(detJ > 1e-14)) continue;
                                  const Eigen::Matrix3d Jinv = J.inverse();

                                  // 小应变 eps（直接累计，避免显式 B 矩阵）
                                  Eigen::Matrix<double,6,1> eps = Eigen::Matrix<double,6,1>::Zero();

                                  idx=0;
                                  for (int kk=0; kk<4; ++kk)
                                      for (int jj=0; jj<4; ++jj)
                                          for (int ii=0; ii<4; ++ii, ++idx){
                                              const int L = lid(ii,jj,kk);

                                              Eigen::Vector3d gparam(dN_dxi[idx], dN_deta[idx], dN_dzet[idx]);
                                              Eigen::Vector3d gxyz   = Jinv * gparam;
                                              const double dNx=gxyz(0), dNy=gxyz(1), dNz=gxyz(2);

                                              const double ux = uc(L);
                                              const double uy = uc(L + 64);
                                              const double uz = uc(L + 128);

                                              eps(0) += dNx * ux;
                                              eps(1) += dNy * uy;
                                              eps(2) += dNz * uz;
                                              eps(3) += dNx * uy + dNy * ux;
                                              eps(4) += dNy * uz + dNz * uy;
                                              eps(5) += dNz * ux + dNx * uz;
                                          }

                                  const double e_point = std::max(0.0, eps.dot(D*eps));

                                  int fi = std::min(linear_subdiv-1, (int)std::floor( (0.5*(xi+1.0))*linear_subdiv ));
                                  int fj = std::min(linear_subdiv-1, (int)std::floor( (0.5*(eta+1.0))*linear_subdiv ));
                                  int fk = std::min(linear_subdiv-1, (int)std::floor( (0.5*(zeta+1.0))*linear_subdiv ));
                                  int fine_local = fk*linear_subdiv*linear_subdiv + fj*linear_subdiv + fi;
                                  if (fine_local<0 || fine_local>=mapping.subdiv_factor) continue;

                                  e_f[fine_local] += e_point * detJ * wgt;
                                  v_f[fine_local] += detJ * wgt;
                              }

                      const auto& c2f = mapping.coarse_to_fine[c];
                      const int nloc = std::min<int>((int)c2f.size(), mapping.subdiv_factor);
                      for (int t=0; t<nloc; ++t) {
                          const int fgid = c2f[t];
                          if (fgid>=0 && fgid<(int)fine_elements.size()){
                              fine_elements[fgid].energy = e_f[t];
                              fine_elements[fgid].volume = std::max(v_f[t], 1e-15);
                          }
                      }
                  });

    mapping.update_fine_volumes_from_elements(fine_elements);
}

void ComputeSensitivitiesParallel(
        std::vector<FineElement>& fine_elements,
        const std::vector<double>& density)
{
    std::vector<int> ids(fine_elements.size());
    std::iota(ids.begin(), ids.end(), 0);

    std::for_each(std::execution::par_unseq, ids.begin(), ids.end(),
                  [&](int f){
                      if (f >= (int)density.size()) return;
                      const double rho = std::clamp(density[f], 1e-4, 1.0);
                      fine_elements[f].sensitivity = (-3.0 * rho * rho) * (fine_elements[f].energy);
                  });
}

void MultiRes_TopOpt_Pipeline(
        const std::string& input_model_file,
        const std::string& bdy_file,
        const std::vector<int>& protect_cell_ids,
        const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
        double target_rate,
        double max_iter,
        double Yangs,
        double Poissons,
        int subdivision_level,
        std::vector<double>& density,
        std::atomic_bool& is_need_update_model)
{
    Timer timer;
    Eigen::setNbThreads(std::min(24, (int)std::thread::hardware_concurrency()));
    Eigen::initParallel();

    // 1) 读取粗网格
    std::vector<CBSplineSolid> coarse_patch;
    ToolFunc::readFileXml(input_model_file, coarse_patch);
    timer.out("read_coarse_xml");

    // 2) 建立映射
    MultiResMapping mapping;
    if (!mapping.init((int)coarse_patch.size(), subdivision_level)) {
        std::cerr << "Failed to initialize mapping\n"; return;
    }
    timer.out("create_mapping");

    // 3) 初始化细密度
    if ((int)density.size() != mapping.fine_size) {
        density.assign(mapping.fine_size, 1.0);
    }
    for (auto& d : density) d = std::clamp(d, 1e-4, 1.0);

    const double tot_num = (double)mapping.fine_size;
    double exist_num = (double)mapping.fine_size;

    // 4) 读取边界条件
    std::vector<CPoint>  bdy_points;
    std::vector<int>     bdy_type;
    std::vector<CVector> bdy_value;
    {
        int control_size = 0;
        std::ifstream in(bdy_file);
        in >> control_size;
        for (int i=0;i<control_size;++i){
            double x,y,z; int cur_type; double vx,vy,vz;
            in>>x>>y>>z>>cur_type>>vx>>vy>>vz;
            bdy_points.emplace_back(CPoint(x,y,z));
            bdy_value.emplace_back(CVector(vx,vy,vz));
            bdy_type.push_back(cur_type);
        }
    }
    timer.out("add bdy info");
    gpmesh::GPLog::init(spdlog::level::trace);

    // 5) 边界点映射
    gpmesh::GPSplineSolidV0 coarse_spline;
    coarse_spline.init_patch(coarse_patch);

    std::vector<boundary_point> neumann_points, dirichlet_points;
    std::vector<MPoint> A, B;
    std::vector<int> global_ids;

    A.reserve(bdy_points.size());
    for (const auto& p : bdy_points) A.push_back({p.getX(),p.getY(),p.getZ()});

    size_t total_cp = 0;
    for (const auto& patch : coarse_patch) total_cp += patch.coefs().size();
    B.reserve(total_cp);
    for (const auto & patch : coarse_patch)
        for (const auto & coef : patch.coefs())
            B.push_back({coef.getX(), coef.getY(), coef.getZ()});

    GetNearestIds(A, B, global_ids);
    if ((int)global_ids.size() != (int)bdy_points.size()) { std::cerr<<"GetNearestIds failed\n"; return; }

    const auto& global_local_map = coarse_spline.m_sparse_system.get_l2g_map();
    for (int i=0;i<(int)bdy_points.size();++i){
        int gid = global_ids[i];
        int patch_id = gid / 64;
        int local_id = gid % 64;
        if (patch_id >= (int)global_local_map.size()) continue;
        if (local_id >= (int)global_local_map[patch_id].size()) continue;
        int global_id = global_local_map[patch_id][local_id];
        if (global_id != -1 && i<(int)bdy_type.size() && i<(int)bdy_value.size()){
            if (bdy_type[i]==2)      neumann_points.emplace_back(condition_type::neumann,   bdy_value[i], global_id);
            else if (bdy_type[i]==1) dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0,0,0), global_id);
        }
    }
    coarse_spline.init_dirichlet_points(dirichlet_points);
    timer.out("setup_boundary");

    // 6) 预计算粗单元局部刚度
    std::vector<Eigen::MatrixXd> coarse_local_matrix(coarse_patch.size(), Eigen::MatrixXd::Zero(192,192));
    coarse_spline.Get_Local_Matrix(coarse_local_matrix, Yangs, Poissons);
    timer.out("compute_coarse_K0");

    // 7) 设置 RHS
    coarse_spline.m_sparse_system.rhs().setZero();
    coarse_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
    coarse_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);

    // 8) 保护单元映射到细桶
    std::vector<int> fine_protect(mapping.fine_size, 0);
    for (auto c_id : protect_cell_ids){
        if (c_id>=0 && c_id<(int)mapping.coarse_to_fine.size())
            for (auto fid : mapping.coarse_to_fine[c_id])
                if (fid>=0 && fid<(int)fine_protect.size()) fine_protect[fid] = 1;
    }

    std::vector<FineElement> fine_elements(mapping.fine_size);

    // 9) 主迭代
    int topopt_iter = 0;
    while (topopt_iter < max_iter) {
        std::cout << "Iteration: " << topopt_iter++ << "\n";

        int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
        if (reduce_num < 1) break;

        // 9.1 等效粗密度
        std::vector<double> coarse_density(coarse_patch.size(), 0.0);
        for (int c=0;c<(int)coarse_patch.size();++c){
            double wsum=0.0, vsum=0.0;
            for (auto f : mapping.coarse_to_fine[c]) {
                if (f>=0 && f<(int)density.size() && f<(int)mapping.fine_volumes.size()){
                    wsum += density[f] * mapping.fine_volumes[f];
                    vsum += mapping.fine_volumes[f];
                }
            }
            double rhobar = (vsum>1e-12) ? (wsum/vsum) : 0.0;
            coarse_density[c] = std::clamp(rhobar, 1e-4, 1.0);
        }

        // 9.2 粗求解
        Eigen::MatrixXd coarse_U;
        coarse_spline.Solver_With_Density(neumann_points, dirichlet_points, coarse_density, coarse_U);

        // 9.3 粗局部位移展开
        Eigen::MatrixXd coarse_UU(coarse_patch.size(), 192);
        coarse_UU.setZero();
        const auto& l2g_map = coarse_spline.m_sparse_system.get_l2g_map();
        int control_nums = coarse_spline.m_sparse_system.getControlNums();
        for (int p_id=0; p_id<(int)coarse_patch.size(); ++p_id){
            if (p_id >= (int)l2g_map.size()) continue;
            for (int local_id=0; local_id<64 && local_id<(int)l2g_map[p_id].size(); ++local_id){
                int gid = l2g_map[p_id][local_id];
                if (gid>=0 && gid<coarse_U.rows()){
                    coarse_UU(p_id, local_id       ) = coarse_U(gid, 0);
                    if (gid+control_nums          < coarse_U.rows()) coarse_UU(p_id, local_id+ 64) = coarse_U(gid+control_nums, 0);
                    if (gid+2*control_nums        < coarse_U.rows()) coarse_UU(p_id, local_id+128) = coarse_U(gid+2*control_nums, 0);
                }
            }
        }

        // 9.4 高斯积分
        std::for_each(std::execution::par_unseq, fine_elements.begin(), fine_elements.end(),
                      [](FineElement& e){ e.energy=0.0; e.volume=0.0; e.sensitivity=0.0; });

        ComputeFineElementEnergies(coarse_UU, coarse_patch, Yangs, Poissons, mapping, fine_elements);
        timer.out("compute_fine_energies_with_quadrature");

        // 9.5 细尺度灵敏度
        ComputeSensitivitiesParallel(fine_elements, density);
        std::vector<double> dc(mapping.fine_size);
        std::transform(std::execution::par_unseq, fine_elements.begin(), fine_elements.end(), dc.begin(),
                       [](const FineElement& e){ return e.sensitivity; });
        timer.out("compute_sensitivities");

        // 9.6 灵敏度滤波
        auto Tempdc = dc;
        std::vector<int> fids(dc.size());
        std::iota(fids.begin(), fids.end(), 0);
        std::for_each(std::execution::par_unseq, fids.begin(), fids.end(), [&](int i){
            if (i<(int)neighbor_info.size() && !neighbor_info[i].empty()){
                double wsum=0.0, acc=0.0;
                double rho_i = std::clamp(density[i], 1e-4, 1.0);
                for (const auto& nb : neighbor_info[i]){
                    int j = nb.first;
                    if (j>=0 && j<(int)dc.size() && j<(int)density.size()){
                        double rho_j = std::clamp(density[j], 1e-4, 1.0);
                        double w = nb.second * rho_j;
                        wsum += w; acc += dc[j]*w;
                    }
                }
                if (wsum>1e-12) Tempdc[i] = acc/(wsum * rho_i);
            }
        });
        dc.swap(Tempdc);
        timer.out("filter_sensitivity");

        // 9.7 OC 更新
        GPU_OCC(density, target_rate, dc, fine_protect, mapping.fine_volumes);

        double min_rho = *std::min_element(density.begin(), density.end());
        double max_rho = *std::max_element(density.begin(), density.end());
        std::cout << "[Density Check] Iter " << (topopt_iter-1)
                  << " | min=" << min_rho << " | max=" << max_rho << "\n";

        double totalV=0.0, rhoV=0.0;
        for (int i=0;i<(int)density.size();++i){
            totalV += mapping.fine_volumes[i];
            rhoV   += mapping.fine_volumes[i]*density[i];
        }
        double curFrac = rhoV/std::max(1e-20,totalV);
        std::cout << "[Volume Check] Iter " << (topopt_iter-1)
                  << " | Current=" << curFrac
                  << " | Target="  << target_rate
                  << " | Diff="    << std::abs(curFrac - target_rate) << "\n";

        is_need_update_model = true;
        timer.out("iteration_"+std::to_string(topopt_iter-1));
    }

    std::cout << "Optimization completed successfully\n";
}