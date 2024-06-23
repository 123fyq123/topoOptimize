#pragma once
#include"head.h"
//#include"tinyXMl/tinyxml.h"
#include<fstream>
#include<regex>

class ToolFunc {
public:
	static void transformGradients(Eigen::MatrixXd jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, int dim);

	//基函数梯度+原值
	static void transformG(std::vector<CVector> jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, Eigen::MatrixXd& G_hat, int dim);

	//高斯点获取
	static void gaussIntegalPoint(int N, double lowpoint, double uppoint, std::vector<double>& point, std::vector<double>& weight);

	//材料矩阵
	static void setC(Eigen::MatrixXd& C, const Eigen::MatrixXd& R, double lambda, double mu);

	//B
	static void setB(Eigen::MatrixXd& B, const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad);

	//删除矩阵第i,j后计算determine
	static void firstMinor(Eigen::MatrixXd& input, int i, int j, Eigen::MatrixXd& result);

	//
	static void outerNormal(const Eigen::MatrixXd& jac, boxSide s, Eigen::VectorXd& result);

	static void outCsv(const Eigen::MatrixXd& info, string filename);

	static void readFileXml(const string filename, std::vector<CBSplineSolid>& geometry);

    static void voigtStressToK(const Eigen::MatrixXd& S, Eigen::VectorXd& Svec);
};


struct time_integration
{
    enum scheme
    {
        explicit_ = 0,  // 显式积分
        implicit_linear = 1,  /// 隐式线性积分
        implicit_nonlinear = 2 /// 隐式非线性积分
    };
};


// 稀疏系统求解方案
struct linear_solver
{
    enum solver
    {
        LU = 0, 
        LDLT = 1, 
        CGDiagonal = 2, 
        BiCGSTABDiagonal = 3 
    };
};


struct newton_verbosity
{
    enum verbosity
    {
        none = 0,  
        some = 1,  
        all = 2  /// output everything
    };
};



enum class newton_status {
    converged, /// 收敛
    interrupted, /// 超出迭代次数
    working, /// 正在求解
    bad_solution     ///失败的解
}; 


struct stress_type
{
    enum type
    {
        von_mises = 0,  /// 冯米斯应力 
        all_2D = 1,  /// 二维方向的应力 两个正应力 一个切应力 （s11,s22,s12)
        normal_3D = 2,  /// 三维正应力 (s11 s22 s33)
        shear_3D = 3   /// 三维切应力 (s12 s13 s23)
    };
};

struct material_law
{
    enum type
    {
        saint_venant_kirchhoff = 0,         /// S = 2*mu*E + lambda*tr(E)*I
        neo_hooke_ln = 1,                   /// S = lambda*ln(J)*C^-1 + mu*(I-C^-1)
        mooney_rivlin = 2,                  
        yeoh = 3,
        bidtanerman = 4
    };
};

struct solveInfo {
    int SolveType;
    int MaxIters;
    double AbsTol;
    double RelTol;
    int Vervosity;
};