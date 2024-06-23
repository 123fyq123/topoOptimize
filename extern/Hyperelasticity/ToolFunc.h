#pragma once
#include"head.h"
#include"tinyXMl/tinyxml.h"
#include<fstream>
#include<regex>
#include <cstdlib>
#include <immintrin.h>


struct CollideInfo {
    CPoint src;
    CPoint move;
};



class ToolFunc {
public:

    static void outputXml(const std::vector<CBSplineSolid>& geometry, string filename);

	static void transformGradients(const Eigen::MatrixXd& jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, int dim);

	//基函数梯度+原值
	static void transformG(std::vector<CVector> jac, int q, const Eigen::MatrixXd& allGrads, Eigen::MatrixXd& trfGradsK, Eigen::MatrixXd& G_hat, int dim);

	//高斯点获取
	static void gaussIntegalPoint(int N, double lowpoint, double uppoint, std::vector<double>& point, std::vector<double>& weight);

	//材料矩阵
	static void setC(Eigen::MatrixXd& C, const Eigen::MatrixXd& R, double lambda, double mu);

    template <int m>
    static void setC(double C[m][m], const Eigen::MatrixXd& R, double lambda, double mu) {
        double a = lambda, b = mu;
        short a_dim = R.cols();
        short dimTensor = m;
        // voigt indices
        Eigen::MatrixX2i v(dimTensor, 2);
        if (a_dim == 2)
            v << 0, 0,
            1, 1,
            0, 1;
        if (a_dim == 3)
            v << 0, 0,
            1, 1,
            2, 2,
            0, 1,
            1, 2,
            0, 2;
        // componentwise definition
        for (short i = 0; i < dimTensor; ++i)
            for (short j = 0; j < dimTensor; ++j)
                C[i][j] = a * R(v(i, 0), v(i, 1)) * R(v(j, 0), v(j, 1)) +
                b * (R(v(i, 0), v(j, 0)) * R(v(i, 1), v(j, 1)) +
                    R(v(i, 0), v(j, 1)) * R(v(i, 1), v(j, 0)));
    }

	//B
	static void setB(Eigen::MatrixXd& B, const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad);


    static void setBT(double B[][6], const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad);
   

    static void setB(double B[][3], const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad);
   

    //
    //static void setBBlock(Eigen::MatrixXd& B_left, Eigen::MatrixXd& B_right, const Eigen::MatrixXd& F, const Eigen::VectorXd& bGrad);

	//删除矩阵第i,j后计算determine
	static void firstMinor(Eigen::MatrixXd& input, int i, int j, Eigen::MatrixXd& result);

	//
	static void outerNormal(const Eigen::MatrixXd& jac, boxSide s, Eigen::VectorXd& result);

	static void outCsv(const Eigen::MatrixXd& info, string filename);

	static void readFileXml(const string filename, std::vector<CBSplineSolid>& geometry);

    static void voigtStressToK(const Eigen::MatrixXd& S, Eigen::VectorXd& Svec);

    template<int m, int n, int k>
    static void matMutly(double A[m][n], double B[n][k], double res[m][k]) {

        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < k; ++j)
                res[i][j] = 0.0;
        }

        for (int aa = 0; aa < n; aa++) {
            for (int mm = 0; mm < m; mm++) {
                for (int nn = 0; nn < k; nn++) {
                    res[mm][nn] += A[mm][aa] * B[aa][nn];
                }
            }
        } // 切线刚度矩阵第一项 位移项

    }


    static std::vector<int> uniformSampleCount(const Eigen::VectorXd& lower,
        const Eigen::VectorXd& upper,
        int numPoints)
    {
        const unsigned d = lower.rows();
        assert(d == upper.rows());

        Eigen::VectorXd span = upper - lower;
        const double volume = span.prod();
        const double h = pow(volume / numPoints, 1.0 / d);

        std::vector<int> np(d);

        for (unsigned i = 0; i < d; ++i)
        {
            np[i] = static_cast<int>(ceil(span[i] / h));
            //assert(np[i] > 0, "Something went wrong, number of points is zero..");
            assert(np[i] > 0);
        }

        return np;
    }
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

class Matrix {
private:
    double** _Matrix;
    size_t _Row, _Column;
public:
    Matrix() :_Matrix(nullptr), _Row(0), _Column(0) {}//默认构造
    Matrix(size_t r, size_t c);
    Matrix(size_t r, size_t c, const double init);
    Matrix(const Matrix& B);
    ~Matrix();
    double& operator()(size_t i, size_t j) { return _Matrix[j][i]; }//访问第i行、第j列的元素
    const double operator()(size_t i, size_t j)const { return _Matrix[j][i]; }//访问第i行、第j列的元素
    Matrix& operator=(Matrix&& B);
    Matrix operator*(const Matrix& B);
    void multi4kernel(double** c, double** a, double** b, int row, int col);
};
