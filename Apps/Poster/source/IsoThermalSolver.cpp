 #include "IsoThermalSolver.h"
 #include <cmath>
 #include <fstream>

 IsoThermalSolver::IsoThermalSolver(double alpha, double dt, double total_time)
     : alpha_(alpha), dt_(dt), total_time_(total_time) {}

 // 添加Bezier体
 void IsoThermalSolver::addBezierSolid(const CubicBezierSolid& solid) {
     solids_.push_back(solid);
 }

 // 生成全局自由度映射
 void IsoThermalSolver::buildGlobalDOFMap() {
     global_ctrl_pts_.clear();
     global_dof_indices_.clear();
     size_t global_id = 0;

     for (const auto& solid : solids_) {
         for (const auto& pt : solid.control_points) {
             std::string key = hashControlPoint(pt);
             if (!global_ctrl_pts_.count(key)) {
                 global_ctrl_pts_[key] = { pt, global_id++ };
             }
             global_dof_indices_.push_back(global_ctrl_pts_[key].global_id);
         }
     }

     // 初始化线性系统
     size_t num_dof = global_ctrl_pts_.size();
     global_A_.resize(num_dof, num_dof);
     global_A_.setZero();
     global_b_.resize(num_dof);
     U_old_.resize(num_dof);
     U_old_.setZero();

     for (const auto& entry : global_ctrl_pts_) {
         const auto& pt = entry.second.coord;
         //U_old_(entry.second.global_id) = 20.0;
         if (std::abs(pt.x()) < 0.1) { // 固定x=0边界
             size_t dof = entry.second.global_id;
             U_old_(dof) = 100.0;
         }
     }
 }

 // 控制点哈希键生成
 std::string IsoThermalSolver::hashControlPoint(const Eigen::Vector4d& pt) const {
     const double precision = 1e-7;
     long x = std::lround(pt.x() / precision);
     long y = std::lround(pt.y() / precision);
     long z = std::lround(pt.z() / precision);
     return std::to_string(x) + "_" + std::to_string(y) + "_" + std::to_string(z);
 }

 // 雅可比矩阵计算
 Eigen::Matrix3d IsoThermalSolver::computeJacobian(const CubicBezierSolid& solid, const Eigen::Vector3d& param) const {
     Eigen::Matrix3d J = Eigen::Matrix3d::Zero();
     Eigen::VectorXd N;
     std::vector<Eigen::Vector3d> grad_N;

     computeBasisFunctions(param, N, grad_N);

     for (size_t i = 0; i < solid.control_points.size(); ++i) {
         const auto& pt = solid.control_points[i];
         const auto& grad = grad_N[i];
         J(0, 0) += grad.x() * pt.x();
         J(0, 1) += grad.y() * pt.x();
         J(0, 2) += grad.z() * pt.x();
         J(1, 0) += grad.x() * pt.y();
         J(1, 1) += grad.y() * pt.y();
         J(1, 2) += grad.z() * pt.y();
         J(2, 0) += grad.x() * pt.z();
         J(2, 1) += grad.y() * pt.z();
         J(2, 2) += grad.z() * pt.z();
     }

     return J;
 }

 // 三次Bernstein基函数及其梯度计算
 void IsoThermalSolver::computeBasisFunctions(const Eigen::Vector3d& param, Eigen::VectorXd& N, std::vector<Eigen::Vector3d>& grad_N) const {
     auto bernstein = [](int i, double t) {
         return (i == 0) ? (1 - t) * (1 - t) * (1 - t) :
             (i == 1) ? 3 * t * (1 - t) * (1 - t) :
             (i == 2) ? 3 * t * t * (1 - t) :
             t * t * t;
         };

     auto d_bernstein = [](int i, double t) {
         return (i == 0) ? -3 * (1 - t) * (1 - t) :
             (i == 1) ? 3 * (1 - t) * (1 - 3 * t) :
             (i == 2) ? 3 * t * (2 - 3 * t) :
             3 * t * t;
         };

     double u = param.x(), v = param.y(), w = param.z();
     Eigen::Vector4d Bu, Bv, Bw, dBu, dBv, dBw;

     for (int i = 0; i < 4; ++i) {
         Bu[i] = bernstein(i, u);
         Bv[i] = bernstein(i, v);
         Bw[i] = bernstein(i, w);
         dBu[i] = d_bernstein(i, u);
         dBv[i] = d_bernstein(i, v);
         dBw[i] = d_bernstein(i, w);
     }

     N.resize(64);
     grad_N.resize(64);
     int idx = 0;
     for (int i = 0; i < 4; ++i) {
         for (int j = 0; j < 4; ++j) {
             for (int k = 0; k < 4; ++k) {
                 N[idx] = Bu[i] * Bv[j] * Bw[k];
                 grad_N[idx] = Eigen::Vector3d(
                     dBu[i] * Bv[j] * Bw[k],
                     Bu[i] * dBv[j] * Bw[k],
                     Bu[i] * Bv[j] * dBw[k]
                 );
                 idx++;
             }
         }
     }
 }

 // 系统矩阵组装
 void IsoThermalSolver::assembleSystem() {
     global_A_.setZero();
     global_b_.setZero();

     // 4x4x4高斯积分规则
     const int num_gauss = 4;
     const double gauss_points[num_gauss] = {
         0.0694318442029, 0.3300094782075, 0.6699905217925, 0.9305681557971
     };
     const double gauss_weights[num_gauss] = {
         0.1739274225687, 0.3260725774313, 0.3260725774313, 0.1739274225687
     };



     for (size_t solid_idx = 0; solid_idx < solids_.size(); ++solid_idx) {
         const auto& solid = solids_[solid_idx];
         size_t local_start = solid_idx * solid.control_points.size();
         std::vector<size_t> local_dof(solid.control_points.size());
         for (size_t i = 0; i < solid.control_points.size(); ++i) {
             local_dof[i] = global_dof_indices_[local_start + i];
         }

         for (int u = 0; u < num_gauss; ++u) {
             for (int v = 0; v < num_gauss; ++v) {
                 for (int w = 0; w < num_gauss; ++w) {
                     Eigen::Vector3d param(gauss_points[u], gauss_points[v], gauss_points[w]);
                     double weight = gauss_weights[u] * gauss_weights[v] * gauss_weights[w];

                     Eigen::Matrix3d J = computeJacobian(solid, param);
                     double detJ = J.determinant();

                     Eigen::VectorXd N;
                     std::vector<Eigen::Vector3d> grad_N;
                     computeBasisFunctions(param, N, grad_N);

                     Eigen::MatrixXd grad_phys(grad_N.size(), 3);
                     for (size_t i = 0; i < grad_N.size(); ++i) {
                         grad_phys.row(i) = (J.inverse().transpose() * grad_N[i]).transpose();
                     }

                     for (size_t i = 0; i < N.size(); ++i) {
                         size_t global_i = local_dof[i];
                         for (size_t j = 0; j < N.size(); ++j) {
                             size_t global_j = local_dof[j];
                             double K_ij = alpha_ * dt_ * grad_phys.row(i).dot(grad_phys.row(j)) * detJ * weight;
                             double M_ij = N(i) * N(j) * detJ * weight;
                             global_A_.coeffRef(global_i, global_j) += K_ij + M_ij;
                         }
                         global_b_(global_i) += U_old_(global_i) * N(i) * detJ * weight;
                     }
                 }
             }
         }
     }
 }

 // 应用Dirichlet边界条件
 void IsoThermalSolver::applyBoundaryConditions() {
     for (const auto& entry : global_ctrl_pts_) {
         const auto& pt = entry.second.coord;
         if (std::abs(pt.x()) < 0.1) { // 固定x=0边界
             size_t dof = entry.second.global_id;
             // 清除行 dof 的所有非零元素
             for (Eigen::SparseMatrix<double>::InnerIterator it(global_A_, dof); it; ++it) {
                 global_A_.coeffRef(it.row(), it.col()) = 0.0;
             }
             // 设置对角线为 1.0
             global_A_.coeffRef(dof, dof) = global_b_(dof) / 100.0;
             //global_b_(dof) = 100.0;
         }
     }
 }

 // 求解时间步
 void IsoThermalSolver::solveTimeStep() {
     Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
     solver.compute(global_A_);
     U_old_ = solver.solve(global_b_);
 }

 void IsoThermalSolver::solveTimeStep_V2() {
     //U_old_ = global_A_.fullPivHouseholderQr().solve(global_b_);

 }

 // 输出VTK文件
 void IsoThermalSolver::outputVTK(int step) const {
     std::ofstream file("ThermalSolver/output_" + std::to_string(step) + ".vtk");
     file << "# vtk DataFile Version 3.0\n";
     file << "Isogeometric Thermal Solution\n";
     file << "ASCII\n";
     file << "DATASET UNSTRUCTURED_GRID\n";

     // 输出控制点坐标
     file << "POINTS " << global_ctrl_pts_.size() << " double\n";
     for (const auto& entry : global_ctrl_pts_) {
         const auto& pt = entry.second.coord;
         file << pt.x() << " " << pt.y() << " " << pt.z() << "\n";
     }

     // 输出温度场
     file << "POINT_DATA " << global_ctrl_pts_.size() << "\n";
     file << "SCALARS Temperature double 1\n";
     file << "LOOKUP_TABLE default\n";
     for (size_t i = 0; i < U_old_.size(); ++i) {
         file << U_old_(i) << "\n";
     }
 }

 // 主运行流程
 void IsoThermalSolver::run() {
     buildGlobalDOFMap();
     int num_steps = total_time_ / dt_;

     std::cout << "start..." << std::endl;
     for (int step = 0; step < num_steps; ++step) {
         assembleSystem();
         applyBoundaryConditions();
         solveTimeStep();
         outputVTK(step);
         std::cout << "Runing...  " << step + 1 << "/" << num_steps << std::endl;
     }
     std::cout << "done" << std::endl;
 }