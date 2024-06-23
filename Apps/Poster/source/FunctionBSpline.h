#pragma once
#include "BSpline3D/BSplineSolid.h"
#include "Common/GeometryKernel.h"
#include<unordered_map>
#include<unordered_set>
class BdyBase
{
public:
	BdyBase(int _local_id, int _global_id, int _solid_id, int _ty = 1) :local_id(_local_id), global_id(_global_id), solid_id(_solid_id), bdy_type(_ty) {};
public:
	int local_id;
	int global_id;
	int solid_id;
	int bdy_type;
};

class simpsort
{
public:
	double senitiv;
	int cell_id;
};

class cube4sort
{
public:
	cube4sort(int _local_id, int _global_id, int _solid_id, double _key) :local_id(_local_id), global_id(_global_id), solid_id(_solid_id)
		, key(_key) {};
public:
	int local_id;
	int global_id;
	int solid_id;
	double key;
};

class FunctionBSpline
{
public:
	FunctionBSpline(BSplineSolid _m_solids) : m_solids(_m_solids) {};
	FunctionBSpline(BSplineSolid _m_solids, std::vector<double> _material_distribution) : m_solids(_m_solids),
		material_distribution(_material_distribution) {
	};
	~FunctionBSpline() {};
	void update_material_distribution(std::vector<double>& _material_distribution)
	{
		material_distribution.clear();
		material_distribution.resize(_material_distribution.size(), 0);
		for (int i = 0; i < _material_distribution.size(); i++)
		{
			//std::cout << "Schedule of update material distribution " << i + 1 << ":" << _material_distribution.size() << std::endl;
			material_distribution[i] = _material_distribution[i];
		}
	}
public:
	BSplineSolid m_solids;
	std::vector<std::vector<double>> PreCalDerBasisFunction;
	std::vector<double> material_distribution;//材料分布
	const double umin = 0.0000000001;
	const double p = 3;//惩罚因子
	std::vector<std::unordered_map<int, int>> colsID;
	std::vector<int> colsize;
public:
	void BdyFunLine(std::vector< std::vector<double> >& m_coef);
	void Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs);
	void Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector< std::vector<double>>& local_bigrhs);
	void Global_MatrixAndRhs1(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector< std::vector<double>>& local_bigrhs, std::map<int, int>& base_check);
	//求高斯积分点, n表示积分点的个数，w是求出的权因子，gausspoint是求出的高斯点
	void GaussIntegalPoint(int n, std::vector<double>& w, std::vector<double>& gausspoint);
	//求解线性方程组
	void SolveLinearEquation(int mtype, const std::vector<std::vector<double>>& coeff_matrix,
		const std::vector<double>& load_vect, std::vector<double>& m_solute);
	//计算数值解的Ｌ２误差
	void Error(std::vector<double> m_coef, std::vector<double>& m_errors);
	double ExactSol(double x, double y, double z);
	double RHS(double x, double y, double z);
	double RHS1(double u, double v, double w);
	double uTKu(std::vector<double> u, Matrix K0);
	double PTu(std::vector<double> P, std::vector<double> u);
public:
	/*-----------------求解面PDE-----------------*/
	const int U = 0, V = 1, W = 2;
	void BdyFunFace(std::vector< std::vector<double> >& m_coef);
	//uvw判断是哪一个面
	void Global_MatrixAndRhsFace(int solid_id, int uvw, int uvw_id, Matrix& m_bigmatrix, std::vector<double>& m_bigrhs);
	void facePDE(int solid_id, int uvw, int uvw_id, std::vector< std::vector<double> >& m_coef,
		std::vector< std::vector<double> > m_coeftest);
public:
	void Global_MatrixAndRhs1_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::map<int, int>& base_check);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<Matrix>& local_matrix1);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs1);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<double>& m_bigrhs1);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector<std::vector<double>>& local_bigrhs);
	void Global_MatrixAndRhs_new(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs);
	void Global_MatrixAndRhs(Matrix& m_bigmatrix, std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs
		, std::vector<std::vector<double>>& local_bigrhs);
public:
	void Global_MatrixAndRhs_Sparse(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<double>& m_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals);
	void Global_MatrixAndRhs_SparseOld(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<double>& m_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals);
	void Global_MatrixAndRhs_Sparse(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);
	void Global_MatrixAndRhs_Sparse_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);
	/*--------------线弹性------------------*/
	const double E = 2 * 1.0e11;
	const double v = 0.3;
	const double lbda = v * E / ((1 + v) * (1 - 2 * v));
	const double mu = E / (2 * (1 + v));
	void Global_MatrixAndRhs_liner(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector<std::vector<double>>& local_bigrhs);
	void Global_MatrixAndRhs_liner_quick(Matrix& m_bigmatrix, std::vector<double>& m_bigrhs, std::vector<Matrix>& local_matrix,
		std::vector<std::vector<double>>& local_bigrhs);
	bool InVecID(int m_id, std::vector<int>& m_vecid)
	{
		if (std::find(m_vecid.begin(), m_vecid.end(), m_id) == m_vecid.end())
			return false;
		else
			return true;
		// for (int i = 0; i < m_vecid.size(); i++)
		// {
			// if (m_id == m_vecid[i])
			// {
				// return true;
			// }
		// }
		// return false;
	}
	int InVecIndex(int m_id, const std::vector<int>& m_vecid)
	{
		

		int ans = 0;
		for (int i = 0; i < m_vecid.size(); ++i) {
			if (m_vecid[i] == m_id) return ans;
			ans++;
		}
		return ans;
		//std::vector<int>::iterator it = std::find(m_vecid.begin(), m_vecid.end(), m_id);
		//return it - m_vecid.begin();
		// for (int i = 0; i < m_vecid.size(); i++)
		// {
			// if (m_id == m_vecid[i])
			// {
				// return i;
			// }
		// }
		// return -1;
	}
	void SolveLinearEquation(int mtype, int GlobalSize, std::vector<int>& rows, std::vector<int>& cols,
		std::vector<double>& elems, std::vector<double>& load_vect, std::vector<double>& m_solute);
	int fount = 0;
	double f1(double x, double y, double z)
	{
		// if(fabs(x-9)<0.01)
			// return 0;
		// else
		return 1;
	}
	double f2(double x, double y, double z)
	{
		// if(fabs(x-9)<0.01)
			// return 0;
		// else
		return 1;
	}
	double f3(double x, double y, double z)
	{
		// if(fabs(x-9)<0.01)
			// return 0;
		// else
		return 1;
	}
	/*double fx(int solid_id, int local_id)
	{
		if(forceorder == 0)
		{
			for(int i = 0;i<bdybaseid.size();++i)
			{
				if(bdybaseid[i].local_id == local_id && bdybaseid[i].solid_id == solid_id)
				{
					fount++;
					return force;
				}
			}
			return 0;
		}
		else
			return 0;
	}
	double fy(int solid_id, int local_id)
	{
		if(forceorder == 1)
		{
			for(int i = 0;i<bdybaseid.size();++i)
			{
				if(bdybaseid[i].local_id == local_id && bdybaseid[i].solid_id == solid_id)
				{
					fount++;
					return force;
				}
			}
			return 0;
		}
		else
			return 0;
	}
	double fz(int solid_id, int local_id)
	{
		if(forceorder == 2)
		{
			for(int i = 0;i<bdybaseid.size();++i)
			{
				if(bdybaseid[i].local_id == local_id && bdybaseid[i].solid_id == solid_id)
				{
					fount++;
					return force;
				};
			}
			return 0;
		}
		else
			return 0;
	}
	std::vector<BdyBase> bdybaseid;
	double force;*/
	int subnum;
	int forceorder = 1;
	double martix_size = 8;

	void Single_Resolution_TopOpt_Init_Step(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
		std::vector<Matrix>& local_matrix, std::vector<double>& m_bigrhs, std::vector<std::vector<int>>& m_colids, std::vector<double>& uu);

	void Single_Resolution_TopOpt_Loop_Step(std::vector<Matrix>& local_matrixs, std::vector<double>& m_bigrhs,
		std::vector<std::vector<int>>& m_colids, std::vector<BdyBase>& bdybaseid, std::vector<double>& u);

	void Global_MatrixAndRhs_Sparse_liner(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);

	void Global_MatrixAndRhs_Sparse_liner_V1(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix, 
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids, 
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick, std::vector<double>& u);


	void Global_MatrixAndRhs_Sparse_liner_V2(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick, std::vector<double>& u);

	void Global_MatrixAndRhs_Sparse_liner_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
		std::vector<double>& u);

	void Global_MatrixAndRhs_Sparse_liner(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
		std::string name);
	void Global_MatrixAndRhs_Sparse_liner1(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid, std::vector<Matrix>& local_matrix,
		std::vector<double>& m_bigrhs, std::vector<std::vector<double>>& local_bigrhs, std::vector<int>& r_rowids,
		std::vector<int>& r_colids, std::vector<double>& r_vals, std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);
public:
	std::vector<std::vector<double>> material_distribution_mutires;///材料分布
	std::vector<double> material_distribution_mutires1;
	/*void Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
		std::vector<std::vector<Matrix>> &local_matrixs, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
		std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
		std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick,
		std::vector<std::vector<V3f>> &centersm,std::vector<std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>> &matrixs_quick,
		std::vector<std::vector<std::vector<int>>> &matrixs_quick_check);

	void Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>> &m_bdycoef, std::vector<BdyBase> &bdybaseid,
		std::vector<std::vector<Matrix>> &local_matrixs, std::vector<double> &m_bigrhs, std::vector<std::vector<std::vector<double>>> &local_bigrhs,
		std::vector<int> &r_rowids, std::vector<int> &r_colids, std::vector<double> &r_vals,
		std::vector<std::vector<int>> &m_colids, std::vector<std::vector<std::pair<int, int>>> &m_bdyquick,
		std::vector<std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>> &matrixs_quick,std::vector<std::vector<std::vector<int>>> &matrixs_quick_check);*/
public:
	void Global_MatrixAndRhs_Sparse_liner_mutires(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
		std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
		std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
		std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
		std::vector<std::vector<V3f>>& centers);
	void Global_MatrixAndRhs_Sparse_liner_mutires_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
		std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
		std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
		std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);

public:
	void Global_MatrixAndRhs_Sparse_hot_mutires(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
		std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
		std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
		std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick,
		std::vector<std::vector<V3f>>& centers);
	void Global_MatrixAndRhs_Sparse_hot_mutires_quick(std::vector<std::vector<double>>& m_bdycoef, std::vector<BdyBase>& bdybaseid,
		std::vector<std::vector<Matrix>>& local_matrixs, std::vector<double>& m_bigrhs, std::vector<std::vector<std::vector<double>>>& local_bigrhs,
		std::vector<int>& r_rowids, std::vector<int>& r_colids, std::vector<double>& r_vals,
		std::vector<std::vector<int>>& m_colids, std::vector<std::vector<std::pair<int, int>>>& m_bdyquick);
};


