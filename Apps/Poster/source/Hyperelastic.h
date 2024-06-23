#pragma once
#include"stiffAssembler.h"
namespace Poster_Needed {
    class Hyperelastic {

    public:
        Hyperelastic(stiffAssembler *assembleInstance);

        Hyperelastic(stiffAssembler *assembleInstance, const Eigen::MatrixXd initSolution);

        std::string status();

        void reset();  // 状态重置

        bool computeUpdate(); // 单步更新

        void solve();         // 开始求解

        void setSolverType(int solveType);

        solveInfo defaultInfo();

        void setMaterial(const pair<int, vector<double>> material);

        Eigen::MatrixXd solution();

        void setInfoNone();

    protected:
        stiffAssembler *m_assembleInstance = nullptr;           // 装配实例
        Eigen::MatrixXd m_solVector;                 // 位移解

        size_t m_numIterations;

        double m_residualNorm;

        double m_initResidualNorm;

        double m_updateNorm;

        double m_initUpdateNorm;

        newton_status m_status;                       // 收敛状态

        std::pair<int, vector<double>> m_material; // 材料参数

        const double tol = 1e-10;

        solveInfo info;

        std::ofstream myfile;
    };
}