#pragma once
#include <Eigen/Sparse>
#include <Eigen/PardisoSupport>
#include <Eigen/Dense>
#include <vector>
#include <map>
#include "Bdy.h"

namespace CCBezier
{
    struct CatmullClarkBezier {
    public:
        // 控制顶点  M
        std::vector<Eigen::Vector3d> control_vertices;
        // 全局映射  patch_size * 64
        std::vector<std::vector<int>> global_id;
        // 局部映射
        int get_local_id(int u_id, int v_id, int w_id) {
            return 16 * w_id + 4 * v_id + u_id;
        }
    };
    class CatmullClarkBezierLinerElasticlSolver
    {
    public:
        CatmullClarkBezierLinerElasticlSolver(CatmullClarkBezier _model, IGA::CCBezierBdy _bdy);
        void Solving();
        void SolvingWithDensity(const std::vector<double>& density, std::vector<std::vector<std::vector<double>>>& local_matrixs);
        Eigen::Vector3d getU(int patch_id, double u, double v, double w);
        std::vector<double> m_u;
    private:
        void getC();
        void initAll();
        void pushToMatrix(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localMats);
        void getJacobi(const Eigen::MatrixXd& basisValue, const std::vector<int>& localIndex, Eigen::MatrixXd& jac);
        void PreComputeBasisValue(std::vector<double>& p_u, std::vector<double>& p_v, std::vector<double>& p_w, std::vector<Eigen::MatrixXd>& basisValue);
        void EvalCatmullClarkBezier2Ders(double& u,
                                         double& v,
                                         double& w,
                                         std::vector<Eigen::VectorXd>& singleValue);
        const double E = 1e9; const double pr = 0.3;
        const double lambda = E * pr / ((1. + pr) * (1. - 2. * pr));
        const double mu = E / (2. * (1. + pr));
        CatmullClarkBezier m_model;
        IGA::CCBezierBdy m_bdy;

        Eigen::MatrixXd C;  // 材料矩阵
        std::vector<std::vector<double>> m_bigmatrix;
        std::vector<std::map<int, int>> m_colIDs;
        std::vector<double>  m_bigrhs;

    };
};



