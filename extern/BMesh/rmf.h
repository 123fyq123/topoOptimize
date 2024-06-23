#pragma once
#include "Common/GeometryKernel.h"
#include "BMesh/quat.hpp"
#include "BMesh/frame.hpp"

class Rmf
{
public:
	Rmf();
	~Rmf();


	//class Frame {
	//public:
	//	V3f r, s, t;
	//	V3f center; // optional

	//	Frame() { r = V3f(1, 0, 0); s = V3f(0, 1, 0); t = V3f(0, 0, 1); }
	//	Frame(const V3f& R, const V3f& S, const V3f& T) { r = R; s = S; t = T; normalize(); }

	//	static Frame fromTR(const V3f& T, const V3f& R) { return Frame(R, T % R, T); }
	//	static Frame fromRS(const V3f& R, const V3f& S) { return Frame(R, S, R % S); }
	//	static Frame fromST(const V3f& S, const V3f& T) { return Frame(S % T, S, T); }
	//	static Frame fromT(const V3f& T) {
	//		V3f Temp = T;
	//		Temp = Temp.normalize();
	//		V3f R = orthogonalVector(Temp).normalize();
	//		return fromTR(T, R);
	//	}

	//	static V3f orthogonalVector(const V3f& n) {
	//		if ((abs(n.y) >= 0.9 * abs(n.x)) && abs(n.z) >= 0.9 * abs(n.x)) return V3f(0.0, -n.z, n.y);
	//		else if (abs(n.x) >= 0.9 * abs(n.y) && abs(n.z) >= 0.9 * abs(n.y)) return V3f(-n.z, 0.0, n.x);
	//		else return V3f(-n.y, n.x, 0.0);
	//	}

	//	void normalize() { r = r.normalize(); s = s.normalize(); t = t.normalize(); };
	//};
	bool Run();
	BMesh::Lvector findR0(BMesh::Frame& c);
	void compute(BMesh::Lvector firstR);
	void getResult();
private:
	std::string algorithmName;
public:

	std::vector<BMesh::Frame> U;
	BMesh::Frame c;
	std::vector<BMesh::Lvector> in_points;
};