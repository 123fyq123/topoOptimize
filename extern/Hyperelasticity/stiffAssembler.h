#pragma once

#include"assembler.h"
#include"ToolFunc.h"
#include<unordered_map>
#include<omp.h>

class stiffAssembler : public assembler {
    int dim = 3;
    double lambda, mu;        // 拉梅常数
    double E = 1e3, pr = 0.3;            // 杨氏模量、泊松率
    std::vector<double> mat;        // 材料模型系数
    double forceScaling = 1.;
    const double tol = 1e-15;
    std::vector<int> controlPointsSumNum;
    bool isCollision = false;
    Eigen::MatrixXd extraRhs;

    SparseSystem m_stiff_system;

protected:
    std::vector<CBSplineSolid> geometry;
    std::vector<boundary_condition> dirichletSides, neumanSides;
    std::vector<boundary_point> dirichletPoints, neumannPoints;
    //Eigen::SparseMatrix<double> matrix;  // 刚度矩阵
    //Eigen::VectorXd rhs;                 // 荷载向量

    //boundary_condition bodyforce = { 0,CVector(0.,0.,-1*9.8*3320) };
    boundary_condition bodyforce = {0, CVector(0., 0., 0.)};

    std::pair<int, vector<double>> m_materLaw;      // 材料

    std::vector<std::vector<Eigen::MatrixXd>> basisValueDer;
    //std::vector<Eigen::MatrixXd> basisValueDer;

    bool isPreCompute = false;

public:
    bool openCollide = false;
    typedef assembler Base;

    // 构造样条体
    stiffAssembler(const std::vector<CBSplineSolid> &m_patch, const std::vector<boundary_condition> &m_dirichletSides,
                   const std::vector<boundary_condition> &m_neumannSides);

    stiffAssembler(const std::vector<CBSplineSolid> &m_patch, const std::vector<boundary_point> &m_dirichletPoints,
                   const std::vector<boundary_point> &m_neumannPoints);

    stiffAssembler() = default;

    void refresh_sides();

    void refresh_points();

    int numDofs();

    // 装配 线性
    void stiffMatrix();

    void test_stiffMatrix(std::vector<Eigen::MatrixXd> &localMats);

    //装配 非线性
    void stiffNonlinearMatrix(const Eigen::MatrixXd &displacement, bool isAssembleMatrix);

    void test_stiffNonlinearMatrix(const Eigen::MatrixXd &displacement,
                                   bool isAssembleMatrix,
                                   std::vector<Eigen::MatrixXd> &localMats);

    void test_stiffNonlinearMatrix(const Eigen::MatrixXd &displacement,
                                   bool isAssembleMatrix,
                                   std::vector<Eigen::MatrixXd> &localMats,
                                   std::vector<Eigen::MatrixXd> &localRhs);

    // 装配
    void rhsVector();

    void rhsVector_NeumannPoint();

    void assemble();

    void assemble_points();

    virtual void assemble_points(const Eigen::VectorXd &solution, bool isAssemblerMatrix = true);

    virtual void assemble(const Eigen::VectorXd &solution, bool isAssemblerMatrix = true);

    void
    assembleNoRhs(const Eigen::VectorXd &solution, std::vector<pair<int, double>> f, bool isAssemblerMatrix = true);

    void setMaterialLaw(const pair<int, vector<double>> material);

    void setYangsPossion(double Yangs, double Possion);

    // 系数更新
    void updateCoefs(Eigen::MatrixXd &coefs);

    // 构造位移矩阵
    void constructCoefs(const Eigen::VectorXd &solution, Eigen::MatrixXd &displacement);

    virtual const Eigen::SparseMatrix<double> &getMatrix() const;   // 获取刚度矩阵
    virtual const Eigen::MatrixXd &getRhs() const;             // 获取右端项


    inline std::vector<std::unordered_set<int>> getConnectPatch() const { return m_system.getConnectPatch(); }

    inline std::vector<std::vector<pair<int, int>>> getInterface() const { return m_system.getInterface(); }

    void preComputeBasis();

    int getIndexById(int patch, int index) { return m_system.getIndexById(patch, index); }

    bool
    isCollide(const Eigen::MatrixXd &dispVector, Eigen::SparseMatrix<double> &newMatrix, Eigen::MatrixXd &extraRhs);

    void clearCollid() {
        extraRhs.setZero(m_system.rhs().rows(), 1);
        isCollision = false;
    }

    void
    WriteParaview(std::vector<CBSplineSolid> &patch, const std::string &fn, Eigen::MatrixXd &solution, unsigned npts,
                  bool isStress = false);

    void WriteSingleParaview(const unsigned patchNum, const std::string &fn, unsigned npts, CBSplineSolid &displacement,
                             CBSplineSolid &dispVector, Eigen::MatrixXd &solution, bool isStress = false);

    Eigen::MatrixXd eval_stress(const Eigen::MatrixXd &u, const Eigen::MatrixXd &displacement, int patchId);

public:
    std::vector<pair<int, double>> Force;
};