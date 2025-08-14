#pragma once
#include <iostream>
#include <vector>
#include <unordered_map>
#include <string>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Common/GeometryKernel.h>


// 三次Bezier体定义
struct CubicBezierSolid {
    std::vector<Eigen::Vector4d> control_points; // 控制点 (x,y,z,w)
};

class IsoThermalSolver {
public:
    IsoThermalSolver(double alpha, double dt, double total_time, double dx);
    void addBezierSolid(const CubicBezierSolid& solid);
    void run(std::vector<std::pair<Eigen::Vector4d, double>>& TemperatureField);

private:
    // 物理参数
    double alpha_;        // 热扩散系数
    double dt_;           // 时间步长
    double total_time_;   // 总时间
    double dx_;

    // Bezier体集合
    std::vector<CubicBezierSolid> solids_;

    // 全局自由度管理
    struct GlobalControlPoint {
        Eigen::Vector4d coord;  // 控制点坐标 (x,y,z,w)
        size_t global_id;        // 全局唯一ID
    };
    std::unordered_map<std::string, GlobalControlPoint> global_ctrl_pts_;
    std::vector<size_t> global_dof_indices_; // 局部到全局ID映射
    std::vector<std::pair<size_t, size_t>> edges_;
    // 线性系统
    Eigen::SparseMatrix<double> global_A_;  // 系统矩阵
    Eigen::VectorXd global_b_;              // 右端向量
    Eigen::VectorXd U_old_;                 // 上一时刻解向量
    Eigen::VectorXd U_new_;

    // 核心函数
    void buildGlobalDOFMap();
    void assembleSystem();
    void assembleSystem_V2();
    void assembleSystem_V3();
    void applyBoundaryConditions();

    void solveTimeStep();
    void outputVTK(int step) const;

    // 等几何辅助函数
    std::string hashControlPoint(const Eigen::Vector4d& pt) const;
    Eigen::Matrix3d computeJacobian(const CubicBezierSolid& solid, const Eigen::Vector3d& param) const;
    void computeBasisFunctions(const Eigen::Vector3d& param, Eigen::VectorXd& N, std::vector<Eigen::Vector3d>& grad_N) const;
};