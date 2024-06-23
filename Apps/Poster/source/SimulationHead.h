#pragma once
/* std */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include<unordered_set>
#include<set>
#include <cstdio>
/* Eigen */
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

#include "CBSplineSolid.h"




struct Solver
{
	/// 求解器
	enum type {
		LDLT = 1, LU = 2
	};
};
struct boundary
{
	/// Identifiers for topological sides.
	enum side {
		west = 1, east = 2, south = 3, north = 4, front = 5, back = 6,
		stime = 7, etime = 8
	};
};
struct condition_type
{
	enum type
	{
		unknownType = -1,
		dirichlet = 0, ///< Dirichlet type
		neumann = 1, ///< Neumann type
		robin = 2  ///< Robin type
	};
};
struct boxSide {
public:
	int m_index;
	int m_patch;
public:
	boxSide(int index, int patch) :m_index(index), m_patch(patch) {}
	int direction() const { return (m_index - 1) / 2; }
	int index() const { return m_index; }
	int sideOrientation()
	{
		return ((m_index + (m_index + 1) / 2) % 2 ? 1 : -1);
	}
};
struct boundary_condition {

	int m_condition_type;
	CVector m_value;
	boxSide m_side;
	int m_dir = 0;
	double m_vec;
public:
	//诺依曼条件
	boundary_condition(int patch, int index, int type, CVector value) :m_side(index, patch), m_condition_type(type),
		m_value(value), m_vec(0.0) {}
	//狄利克雷条件
	boundary_condition(int patch, int index, int type, double vec, int dir) :m_side(index, patch), m_condition_type(type),
		m_vec(vec), m_dir(dir) {}

	boundary_condition(int patch, CVector value) :m_side(0, patch), m_value(value), m_condition_type(condition_type::neumann),
		m_vec(0.0) {}

	int patchSide() const { return m_side.m_patch; }

	double getValue() const { return m_vec; }

	int getDir() const { return m_dir; }

	boxSide side() const { return m_side; }

	void eval(std::vector<double> u, std::vector<double> v, Eigen::MatrixXd& res) {
		int len = u.size() * v.size();
		res.resize(2, len);
		for (int i = 0; i < len; ++i) {
			res(0, i) = m_value.getX();
			res(1, i) = m_value.getY();
		}
	}
	void eval(std::vector<double> u, std::vector<double> v, std::vector<double> w, Eigen::MatrixXd& res) {
		int len = u.size() * v.size() * w.size();
		res.resize(3, len);
		for (int i = 0; i < len; ++i) {
			res(0, i) = m_value.getX();
			res(1, i) = m_value.getY();
			res(2, i) = m_value.getZ();
		}
	}
};

struct interFace {  // 供二维使用
public:
	int p1, p2;     // 相交的片号
	int dir1, dir2; //相交的面号
	//std::vector<int> block1; // p1相交的局部索引
	//std::vector<int> block2; // p2相交的局部索引
	//std::vector<CPoint> loc;
};
struct LocInd {
	double x, y, z;
	int index;     // 索引号
	int p;         // 片号
	bool operator < (const LocInd& _rhs) const {
		if (fabs(x - _rhs.x) > 1e-15) return x < _rhs.x;
		if (fabs(y - _rhs.y) > 1e-15) return y < _rhs.y;
		return z < _rhs.z;
	}
	bool operator == (const LocInd& _rhs) const {
		return fabs(x - _rhs.x) < 1e-15 && fabs(y - _rhs.y) < 1e-15 && fabs(z - _rhs.z) < 1e-15;
	}
};
struct CFace {
	int p;			// 块号
	int dir;        // 面号
	std::vector<LocInd> locInfo;
public:
	bool operator==(const CFace& it) {
		if (locInfo.size() != it.locInfo.size()) return false;
		for (int i = 0; i < locInfo.size(); ++i) {
			if (locInfo[i] == it.locInfo[i])
				continue;
			return false;
		}
		return true;
	}
	bool operator<(const CFace& it) {
		int len = std::min(locInfo.size(), it.locInfo.size());
		for (int i = 0; i < len; ++i)
		{
			if (locInfo[i] == it.locInfo[i]) {
				continue;
			}
			else
				return locInfo[i] < it.locInfo[i];
		}
		return false;
	}
};