#pragma once
#include"head.h"
#include"ToolFunc.h"
#include<unordered_map>
#include<omp.h>

class stiffAssembler {
	int dim = 3;
	double lambda, mu;		// 拉梅常数
	double E = 74e9, pr = 0.33;			// 杨氏模量、泊松率
	std::vector<double> mat;		// 材料模型系数
	int coefsNum = 0;
	double forceScaling;
	const double tol = 1e-15;
protected:
	std::vector<CBSplineSolid> geometry;
	std::vector<boundary_condition> dirichletSides, neumanSides;
	Eigen::SparseMatrix<double> matrix;  // 刚度矩阵
	Eigen::VectorXd rhs;                 // 荷载向量
	std::vector<std::vector<int>> globalIndex;  // 局部到全局索引 片号 局部索引号 -> 全局索引号
	int controlpoints = 0;          // 控制点数量
	std::vector<std::vector<int>> gl_coefs_index; // 全局到局部 作用是方便更改控制点坐标

	std::vector<int> globalDirichlet;       // 标记所有索引号是否为迪利克雷点，是 1 否 0
	map<int, double> dir_index;             // 记录迪利克雷点的索引号，为右端项赋值  

	std::vector<unordered_map<int, double>> sparseValue;		// 行优先 记录稀疏矩阵所有非零值
	boundary_condition bodyforce = { 0,CVector(0.,0.,0.) };

	std::pair<int, vector<double>> m_materLaw;      // 材料

public:
	// 构造样条体
	stiffAssembler(std::vector<CBSplineSolid> m_patch, std::vector<boundary_condition> m_dirichletSides,
		std::vector<boundary_condition> m_neumannSides, std::vector<interFace> m_inf);
	// 根据控制点记录局部到全局的索引
	void computerIndex();

	int numDofs();

	//bool computerTopogy(std::vector<interFace>& inf);
	//void findBoundaryLoc(int blockIndex, std::vector<CFace>& blockFace);

	// 装配 线性
	void stiffMatrix();


	//装配 非线性
	void stiffNonlinearMatrix(const Eigen::MatrixXd& displacement, bool isAssembleMatrix);

	//


	// 装配
	void rhsVector();

	// 预计算机迪利克雷点
	void setDirichlet();

	void assemble();

	void assemble(const Eigen::VectorXd& solution, bool isAssemblerMatrix = true);

	void setMaterialLaw(const pair<int, vector<double>> material);

	// 系数更新
	void updateCoefs(Eigen::MatrixXd& coefs);

	// 构造位移矩阵
	void constructCoefs(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement);

	const Eigen::SparseMatrix<double>& getMatrix() const;   // 获取刚度矩阵
	const Eigen::VectorXd& getRhs() const;             // 获取右端项
};