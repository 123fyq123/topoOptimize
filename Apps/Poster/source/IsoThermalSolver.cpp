#include "IsoThermalSolver.h"
#include <cmath>
#include <fstream>

IsoThermalSolver::IsoThermalSolver(double alpha, double dt, double total_time, double dx)
    : alpha_(alpha), dt_(dt), total_time_(total_time), dx_(dx) {
}

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
    U_new_.resize(num_dof);
    U_new_.setZero();

    const double PI = 3.14159265358979323846;
    for (const auto& entry : global_ctrl_pts_) {
        const auto& pt = entry.second.coord;
        size_t dof = entry.second.global_id;
        U_old_(dof) = 100 * sin(PI * pt.x()) * sin(PI * pt.y()) * sin(PI * pt.z());
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
                        //grad_phys.row(i) = (J.inverse().transpose() * grad_N[i]).transpose();
                        grad_phys.row(i) = (J.transpose().lu().solve(grad_N[i])).transpose();
                    }

                    for (size_t i = 0; i < N.size(); ++i) {
                        size_t global_i = local_dof[i];
                        for (size_t j = 0; j < N.size(); ++j) {
                            size_t global_j = local_dof[j];
                            double K_ij = alpha_ * dt_ * grad_phys.row(i).dot(grad_phys.row(j)) * detJ * weight;
                            double M_ij = N(i) * N(j) * detJ * weight;
                            global_A_.coeffRef(global_i, global_j) += K_ij + M_ij;
                            //global_A_.coeffRef(global_i, global_j) += K_ij;
                            global_b_(global_i) += N(i) * N(j) * detJ * weight * U_old_(global_i);
                        }
                        //global_b_(global_i) += U_old_(global_i) * N(i) * detJ * weight;
                    }
                }
            }
        }
    }
}


void IsoThermalSolver::assembleSystem_V2() {
    global_A_.setZero();
    global_b_.setZero();

    for (const auto& entry1 : global_ctrl_pts_) {
        const auto& pt1 = entry1.second.coord;
        size_t dof1 = entry1.second.global_id;

        //if (std::abs(pt1.x()) < 0.01) {
        //    continue;
        //}
        double diag = 0.0;

        for (const auto& entry2 : global_ctrl_pts_) {
            const auto& pt2 = entry2.second.coord;
            size_t dof2 = entry2.second.global_id;

            if (pt1 == pt2) continue;
            double x1 = entry1.second.coord.x();
            double y1 = entry1.second.coord.y();
            double z1 = entry1.second.coord.z();
            double x2 = entry2.second.coord.x();
            double y2 = entry2.second.coord.y();
            double z2 = entry2.second.coord.z();
            double dist = sqrt((x1 - x2) * (x1 - x2) + (y1 - y2) * (y1 - y2) + (z1 - z2) * (z1 - z2));
            if (dist >= 0.09) continue;
            double val = alpha_ * dt_ / (dist * dist);
            global_A_.coeffRef(dof1, dof2) -= val;
            //global_A_.coeffRef(dof1, dof1) += val;
            diag += val;
        }
        global_A_.coeffRef(dof1, dof1) = 1.0 + diag;

    }

    for (int i = 0; i < U_old_.size(); ++i) {
        global_b_(i) = U_old_(i);
    }

    // 示例：假设每个Bezier体有27个控制点
    //for (const auto& solid : solids_) {
    //    const int num_ctrl = solid.control_points.size();
    //    Eigen::MatrixXd A_local = Eigen::MatrixXd::Zero(num_ctrl, num_ctrl);
    //    Eigen::VectorXd b_local = Eigen::VectorXd::Zero(num_ctrl);

    //    // 伪代码：在每个单元中进行积分（高斯点积分）
    //    for (int i = 0; i < num_ctrl; ++i) {
    //        b_local(i) += U_old_(global_dof_indices_[i]);  // 显式项（右端）
    //    }

    //    // 组装进全局
    //    for (int i = 0; i < num_ctrl; ++i) {
    //        size_t row = global_dof_indices_[i];
    //        global_b_(row) += b_local(i);
    //        for (int j = 0; j < num_ctrl; ++j) {
    //            size_t col = global_dof_indices_[j];
    //            global_A_.coeffRef(row, col) += A_local(i, j);
    //        }
    //    }
    //}
}


void IsoThermalSolver::assembleSystem_V3() {
    U_new_.setZero();

    //const double dx = 0.0833333;

    for (const auto& entry1 : global_ctrl_pts_) {
        const auto& pt1 = entry1.second.coord;
        size_t dof1 = entry1.second.global_id;

        if (std::abs(pt1.x()) < 0.0001 || std::abs(pt1.x() - 1) < 0.0001 || std::abs(pt1.y()) < 0.0001 || std::abs(pt1.y() - 1) < 0.0001 || std::abs(pt1.z()) < 0.0001 || std::abs(pt1.z() - 1) < 0.0001) {
            continue;
        }
        double diag = 0.0;
        for (const auto& entry2 : global_ctrl_pts_) {
            const auto& pt2 = entry2.second.coord;
            size_t dof2 = entry2.second.global_id;

            if (dof1 == dof2) continue;
            double x1 = entry1.second.coord.x();
            double y1 = entry1.second.coord.y();
            double z1 = entry1.second.coord.z();
            double x2 = entry2.second.coord.x();
            double y2 = entry2.second.coord.y();
            double z2 = entry2.second.coord.z();
            double dist = std::abs(x1 - x2) + std::abs(y1 - y2) + std::abs(z1 - z2);
            if (dist > dx_ + 1e-4) continue;

            diag += U_old_(dof2);
        }

        U_new_(dof1) = U_old_(dof1) + alpha_ * dt_ * (diag - 6.0 * U_old_(dof1)) / (dx_ * dx_);

    }
    U_old_ = U_new_;
}

// 应用Dirichlet边界条件
void IsoThermalSolver::applyBoundaryConditions() {
    for (const auto& entry : global_ctrl_pts_) {
        const auto& pt = entry.second.coord;
        if (std::abs(pt.x()) < 0.0001 || std::abs(pt.x() - 1) < 0.0001 || std::abs(pt.y()) < 0.0001 || std::abs(pt.y() - 1) < 0.0001 || std::abs(pt.z()) < 0.0001 || std::abs(pt.z() - 1) < 0.0001) { // 固定x=0边界
            size_t dof = entry.second.global_id;
            // 清除行 dof 的所有非零元素
            for (Eigen::SparseMatrix<double>::InnerIterator it(global_A_, dof); it; ++it) {
                global_A_.coeffRef(it.row(), it.col()) = 0.0;
            }
            // 设置对角线为 1.0
            global_A_.coeffRef(dof, dof) = 1.0;
            global_b_(dof) = 0.0;
        }
    }
}

// 求解时间步
void IsoThermalSolver::solveTimeStep() {
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> solver;
    solver.compute(global_A_);
    U_old_ = solver.solve(global_b_);
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
    //for (size_t i = 0; i < U_old_.size(); ++i) {
    //    file << U_old_(i) << "\n";
    //}

    for (const auto& entry : global_ctrl_pts_) {
        const auto& id = entry.second.global_id;
        file << U_old_(id) << "\n";
    }
    file.close();
}

// 主运行流程
void IsoThermalSolver::run(std::vector<std::pair<Eigen::Vector4d, double>>& TemperatureField) {
    buildGlobalDOFMap();
    int num_steps = total_time_ / dt_;

    std::cout << "start..." << std::endl;
    for (int step = 0; step < num_steps; ++step) {
        //assembleSystem();
        //assembleSystem_V2();
        assembleSystem_V3();
        //applyBoundaryConditions();
        //solveTimeStep();
        //solveTimeStep_V2();
        outputVTK(step);
        std::cout << "Running...  " << step + 1 << "/" << num_steps << std::endl;
    }

    size_t max_idx = 0;
    for (const auto& entry : global_ctrl_pts_) {
        max_idx = std::max(max_idx, entry.second.global_id);
    }
    TemperatureField = std::vector<std::pair<Eigen::Vector4d, double>>(max_idx + 1, {{0,0,0,0}, 0.0});
    for (const auto& entry : global_ctrl_pts_) {
        const auto& pt_info = entry.second;
        size_t id = pt_info.global_id;
        TemperatureField[id] = std::make_pair(pt_info.coord, U_old_(id));
    }

    std::cout << "done" << std::endl;
}