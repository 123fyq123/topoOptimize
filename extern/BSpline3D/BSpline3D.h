#pragma once
#include"BSplineSolid.h"

typedef std::vector<std::vector<double>> KnotVec;
typedef std::vector<std::vector<std::vector<V3f>>> ControlPoints;
typedef std::vector<int> Orders;

class BSpline3D {
private:
	// n 个节点
	// p次
	// m个控制顶点
	// m + p + 1 = n
	KnotVec knot_vec;                   // 节点向量  n1, n2, n3
	ControlPoints control_points;       // 控制顶点  m1, m2, m3;
	Orders orders;                      // 每一维度的阶数  p1, p2, p3;
	const int dim = 3;
public:
	BSpline3D(const KnotVec& _kv, const ControlPoints& _cp, const Orders& _od) {
		knot_vec = _kv;
		control_points = _cp;
		orders = _od;
	}


};