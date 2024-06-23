#pragma once
#include "Common/GeometryKernel.h"
#include "mkl.h"

typedef std::vector< std::vector<double> > Matrix;

class BSplineSolid
{
public:
	std::vector<HexV3fMesh> Local_BSplineSolid;//局部的小的每个单元
	std::vector<std::vector<std::vector<std::vector<VertexHandle>>>> Control_Point_id;//每一个块按顺序排序控制点的下标
	std::vector<std::vector<VertexHandle>> CP_Local_to_Global;//局部的小单元的每个点对应于全局模型中的点
	std::vector<std::vector<int>> Local_Order;//局部的每个小单元的u,v,w三个方向的次数
	std::vector<std::vector<std::vector<int>>> Local_Node;//局部的每个小单元的u,v,w三个方向的节点向量
public:
	HexV3fMesh Global_BSplineSolid;//全局的大的整个模型
public://新的结构
	std::vector<std::vector<std::vector<std::vector<int>>>> Control_Point;//每一个小块的控制点
	std::vector<std::vector<V3f>> Local_V3f;
	std::vector<std::vector<int>> Local_to_Global;
	std::vector<std::vector<bool>> local_bdy;
	int Global_size;
	int Global_size_without_bdy;
public:
	BSplineSolid() {};
	~BSplineSolid()
	{
		Local_BSplineSolid.clear();
		Control_Point_id.clear();
		CP_Local_to_Global.clear();
		Local_Order.clear();
		Local_Node.clear();

		Control_Point.clear();
		Local_V3f.clear();
		Local_to_Global.clear();
		local_bdy.clear();
	}
	void init(int n)//初始化BSplineSolid
	{
		Local_BSplineSolid.resize(n);
		Control_Point_id.resize(n);
		CP_Local_to_Global.resize(n);
		Local_Order.resize(n);
		Local_Node.resize(n);

		Control_Point.resize(n);//每一个小块的控制点
		Local_V3f.resize(n);
		Local_to_Global.resize(n);
		local_bdy.resize(n);
		Global_size = 0;
		Global_size_without_bdy = 0;
	}
	BSplineSolid& operator=(const BSplineSolid& bsplinesolid)
	{
		this->Local_BSplineSolid = bsplinesolid.Local_BSplineSolid;
		this->Control_Point_id = bsplinesolid.Control_Point_id;
		this->CP_Local_to_Global = bsplinesolid.CP_Local_to_Global;
		this->Local_Order = bsplinesolid.Local_Order;
		this->Local_Node = bsplinesolid.Local_Node;
		this->Global_BSplineSolid = bsplinesolid.Global_BSplineSolid;

		this->Control_Point = bsplinesolid.Control_Point;
		this->Local_V3f = bsplinesolid.Local_V3f;
		this->Local_to_Global = bsplinesolid.Local_to_Global;
		this->local_bdy = bsplinesolid.local_bdy;
		this->Global_size = bsplinesolid.Global_size;
		this->Global_size_without_bdy = bsplinesolid.Global_size_without_bdy;
		return *this;
	}
public:
	//求点x在T中的位置
	int X_LOCATION(int solid_id, int uvw, double x);
	//计算几何变换.X=G(u,v),也就是曲面在u,v的值。
	V3f Geometry(int solid_id, double u, double v, double w);
	//m表示节点向量的最大下标，即节点向量共有m+1个节点，i表示基Ni（i从0开始算），
	//n表示求导到n次，结果在ders中（ders共n+1个元素）。
	//void DersOneBasisFun(int degree,int m,vector<double> U, int i,double u,int n,double ders[])
	void DersOneBasisFun(int solid_id, int uvw, int i, double u, int n, double ders[]);
	//对边界函数近似求基的系数，m_bdy[0]=1表示u方向，m_bdy[1]=1表示v方向，m_bdy[2]=1表示w方向, m_para表示另一个方向的参数值
	void SolveCoef(int solid_id, std::vector<double> m_scoef, std::vector<int> m_para, std::vector<int> m_bdy, std::vector<double>& m_coef);
	//寻找结点向量中所含的区域个数
	void IntegralArea(int solid_id, std::vector<double>& u_interval, std::vector<double>& v_interval, std::vector<double>& w_interval);
	//计算一系列点上的雅可比矩阵的转置的逆及行列式
	void JacobiMatrix(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
		std::vector<double>& jacobidet, std::vector<Matrix>& jacobi);
	double ExactSol(double x, double y, double z);
	double RHS(double x, double y, double z);
	double RHS1(double u, double v, double w);
public:
	void JacobiMatrix_new(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
		std::vector<double>& jacobidet, std::vector<Matrix>& jacobi);
	void JacobiMatrix_new(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
		std::vector<double>& jacobidet, std::vector<Matrix>& jacobi, std::vector <double>& PreCalDerBasisFunction, bool isPreCal);
	V3f Geometry_new(int solid_id, double u, double v, double w);
	V3f Geometry_new(int solid_id, int uid, int vid, int wid, const std::vector<double>& PreCalDerBasisFunction);
public:
	/*-----------------求解面PDE-----------------*/
	const int U = 0, V = 1, W = 2;
	void JacobiMatrixFace(int solid_id, int uvw, int uvw_id, std::vector<double> u, std::vector<double> v,
		std::vector<double>& jacobidet, std::vector<Matrix>& jacobi);
};


