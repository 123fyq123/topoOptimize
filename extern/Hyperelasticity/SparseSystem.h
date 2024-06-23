#pragma once
#include"head.h"
class SparseSystem {
    int controlpoints = 0;
    int dim = 0;
    int coefsNum = 0;

protected:
    Eigen::SparseMatrix<double> m_matrix;  // 刚度矩阵
    Eigen::MatrixXd m_rhs;                 // 荷载向量
    std::vector<std::unordered_set<int>> m_ConnectPatch;

    std::vector<std::vector<int>> m_globalIndex;  // 局部到全局索引 片号 局部索引号 -> 全局索引号

    std::vector<std::vector<int>> m_gl_coefs_index; //全局到局部索引

    std::vector<std::vector<double>> m_bigmatrix;     // 存储稀疏值

    std::vector<std::map<int,int>> m_colids;           // 存储稀疏值对应的列号

    std::vector<bool> m_globalDirichlet;       // 标记所有索引号是否为迪利克雷点
    std::vector<int> m_globalDirichlet_count;    // 索引号为迪利克雷点的标记次数

    std::vector<double> m_dirichletValue;        // 记录迪利克雷点索引值

    std::vector<std::vector<pair<int, int>>> m_inter;  // 对接边

    std::vector<uint32_t> m_cpoints_inclusive_sum;
public:
    SparseSystem(){}
    SparseSystem(bool b) {}

    //empty SparseSystem分成两步初始化 for dirichlet point
    void init_geometry(const std::vector<CBSplineSolid>& geometry, int dims);
    void init_dirichlet(const std::vector<boundary_point>& dirichletPoints);

    SparseSystem(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_condition>& dirichletSides, int dims);
    SparseSystem(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_point>& dirichletPoints, int dims);

    void computerInterface(const std::vector<CBSplineSolid>& geometry);
    void computerIndex(const std::vector<CBSplineSolid>& geometry);

    void setDirichlet(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_condition>& dirichletSides);
    void setDirichlet(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_point>& dirichletPoints);

    int numControlPoints() { return controlpoints; }

    void createStore(const std::vector<CBSplineSolid>& geometry);  // 稀疏矩阵存储空间预计算


    void constructCoefs(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement);

    std::vector<std::vector<double>>& test_get_bigmatrix()
    {
        return m_bigmatrix;
    }

    Eigen::SparseMatrix<double>& matrix()
    {
        return m_matrix;
    }

    const Eigen::SparseMatrix<double>& matrix() const
    {
        return m_matrix;
    }

    Eigen::MatrixXd& rhs() {
        return m_rhs;
    }

    const Eigen::MatrixXd& rhs() const {
        return m_rhs;
    }

    inline std::vector<std::unordered_set<int>> getConnectPatch() const { return m_ConnectPatch; }

    inline const std::vector<std::vector<int>>& get_l2g_map() const { return m_globalIndex; };
    inline const std::vector<std::vector<int>>& get_g2l_map() const { return m_gl_coefs_index; };
    inline std::vector<std::vector<pair<int, int>>> getInterface() const { return m_inter; }
    inline std::vector<uint32_t>& get_cpoints_inclusive_sum() { return m_cpoints_inclusive_sum; };

    int getDim() {
        return dim;
    }

    int getControlNums() const {
        return controlpoints;
    }

    void pushToMatrix(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localMats);
    void pushToMatrix(uint32_t patch_id, void* data);
    void pushToRhs(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localRhs);
    void pushToRhs(const std::vector<pair<int,double>>& f);

    void pushToSystem(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localRhs, const Eigen::MatrixXd& localMats, bool isAssembleMatrix);
    void dataPush(bool isOne = true);
    void dataPush_Dirichlet(bool isOne = true);
    void dataPush_Dirichlet_V2();

    void dataPush(Eigen::MatrixXd& extraRhs);

    int getIndexById(int patch, int index);

    bool IsDirichlet(int patch, int index);

    bool IsDirichlet(uint32_t global_id);

    void constructNewMatrix(Eigen::VectorXd& extraRhs, Eigen::SparseMatrix<double>& newMatrix, Eigen::MatrixXd& newRhs);

};