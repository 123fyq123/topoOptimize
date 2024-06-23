#include "rmf.h"
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#define ZERO_NORM 1e-7

Rmf::Rmf() {}
Rmf::~Rmf() {}

/**
 * @description: find R0
 * @param {type} c endFrame of a skel. V3f ts of the line
 * @return {type}
 */

//BMesh::Lvector Rmf::findR0(BMesh::Frame& c) {
//	BMesh::Lvector firstT = (in_points[1] - in_points[0]).normalize();
//	BMesh::Lvector T(1, 1, 1);
//	T = T.normalize();
//	BMesh::Lvector R(1, -1, 0);
//	R = R.normalize();
//	BMesh::Quat q;
//	q.setAxisAngle(BMesh::Lvector(T.x, T.y, T.z), M_PI / 180.00);
//	std::vector<BMesh::Lvector> Rs;
//	//Rs.push_back(BMesh::Lvector(R.x, R.y, R.z));
//	BMesh::Lvector rs = BMesh::Lvector(R.x, R.y, R.z);
//	for (int i = 0; i < 360; i++) {
//		rs = q * rs;
//		Rs.push_back(rs);
//	}
//	// test point
//	/*std::ofstream out("testpoints.txt");
//	for (int i = 0; i < Rs.size(); i++) {
//		out << Rs[i].x << " " << Rs[i].y << " " << Rs[i].z << "\n";
//	}
//	out.close();*/
//	BMesh::Quat qs;
//	BMesh::Lvector ts(T.x, T.y, T.z);
//	qs.QuatWithtwoVec(BMesh::Lvector(T.x, T.y, T.y), BMesh::Lvector(firstT.x, firstT.y, firstT.z));
//	double minValue = INT_MIN;
//	BMesh::Lvector rlast;
//	for (int i = 0; i < 360; i++) {
//		BMesh::Lvector r = BMesh::Lvector(Rs[i].x, Rs[i].y, Rs[i].z);
//		BMesh::Lvector s = ts % r;
//		double countValue = 0;
//		for (int j = 0; j < 3; j++) {
//			if (j == 0) {
//				countValue =  fabs(r * c.r) ;
//			}
//			else if (j == 1) {
//				countValue =  fabs(r * c.s);
//			}
//			else if (j == 2) {
//				countValue = fabs(r * c.t) ;
//			}
//			if (minValue < countValue) {
//				rlast = r;
//				minValue = countValue;
//			}
//		}
//	}
//	std::cout << "[DDBUG] " << minValue << std::endl;
//	return rlast;
//}

BMesh::Lvector Rmf::findR0(BMesh::Frame& c) {
	BMesh::Lvector firstT = (in_points[1] - in_points[0]).normalize();
	BMesh::Lvector T(1, 1, 1);
	T = T.normalize();
	BMesh::Lvector R(1, -1, 0);
	R = R.normalize();
	BMesh::Quat q;
	q.setAxisAngle(BMesh::Lvector(T.x, T.y, T.z), M_PI / 180.00);
	std::vector<BMesh::Lvector> Rs;
	//Rs.push_back(BMesh::Lvector(R.x, R.y, R.z));
	BMesh::Lvector rs = BMesh::Lvector(R.x, R.y, R.z);
	for (int i = 0; i < 360; i++) {
		rs = q * rs;
		Rs.push_back(rs);
	}
	// test point
	std::ofstream out("testpoints.txt");
	for (int i = 0; i < Rs.size(); i++) {
		out << Rs[i].x << " " << Rs[i].y << " " << Rs[i].z << "\n";
	}
	out.close();
	BMesh::Quat qs;
	BMesh::Lvector ts(T.x, T.y, T.z);
	qs.QuatWithtwoVec(BMesh::Lvector(T.x, T.y, T.y), BMesh::Lvector(firstT.x, firstT.y, firstT.z));
	for (int i = 0; i < 360; i++) {
		Rs[i] = qs * Rs[i];
	}
	/*std::ofstream out1("testpoints1.txt");
	for (int i = 0; i < Rs.size(); i++) {
	 out1 << Rs[i].x << " " << Rs[i].y << " " << Rs[i].z << "\n";
	}
	out1.close();*/

	std::ofstream out2("testpoints2.txt");
	out2 << "0 0 0" << "\n";
	out2 << firstT.x << " " << firstT.y << " " << firstT.z << "\n";
	out2.close();
	double minValue = INT_MIN;
	BMesh::Lvector rlast;
	for (int i = 0; i < 360; i++) {
		BMesh::Lvector r = BMesh::Lvector(Rs[i].x, Rs[i].y, Rs[i].z);
		BMesh::Lvector s = ts % r; // ts 放在这里好像不太对
		double countValue = 0;
		for (int j = 0; j < 3; j++) {
			if (j == 0) {
				countValue = fabs(r * c.r);
			}
			else if (j == 1) {
				countValue = fabs(r * c.s);
			}
			else if (j == 2) {
				countValue = fabs(r * c.t);
			}
			if (minValue < countValue) {
				rlast = r;
				minValue = countValue;
			}
		}
	}
	std::cout << "[DDBUG] " << minValue << std::endl;
	return rlast;
}

bool Rmf::Run() {
	if (in_points.size() <= 1) return false;
	if (in_points.size() > 1 && (in_points[0] - in_points[1]).norm() > ZERO_NORM) {
		compute(findR0(c));
	}
	
	return true;
}


static inline BMesh::Lvector pointOnPlane(BMesh::Lvector p, BMesh::Lvector plane_normal, double plane_d = 0)
{
	double t = plane_normal * p - plane_d;
	return p - BMesh::Lvector(t * plane_normal.x, t * plane_normal.y, t * plane_normal.z);
}

void Rmf::compute(BMesh::Lvector firstR)
{
	// Reset computation
	std::vector<BMesh::Lvector> tangent;

	// Estimate tangents
	for (int i = 0; i < in_points.size() - 1; i++) {
		BMesh::Lvector t = in_points[i + 1] - in_points[i];
		if (i > 0 && t.norm() < ZERO_NORM)
			tangent.push_back(tangent.back());
		else
			tangent.push_back((t).normalize());
	}
	tangent.push_back(tangent.back());

	BMesh::Lvector firstT = tangent.front().normalize();

	// Make sure firstR is perpendicular to plane of first tangent
	firstR = pointOnPlane(firstR, firstT);
	if (firstR.norm() > ZERO_NORM)
		firstR.normalize();
	else {
		firstR = BMesh::Frame::orthogonalVector(firstT);
	}


	// First frame
	BMesh::Frame firstFrame = BMesh::Frame::fromTR(firstT, firstR.normalize());

	U.clear();
	U.push_back(firstFrame);

	// Double reflection method: compute rotation minimizing frames
	for (int i = 0; i < in_points.size() - 1; i++)
	{
		BMesh::Lvector ri = U.back().r, ti = U.back().t, tj = tangent[i + 1];

		/*1 */ BMesh::Lvector v1 = in_points[i + 1] - in_points[i];
		if (v1.norm() < ZERO_NORM) { U.push_back(U.back()); continue; }
		/*2 */ double c1 = v1 * v1;
		/*3 */ BMesh::Lvector rLi = ri - (v1 * (v1 * ri)) * (2.0 / c1);
		/*4 */ BMesh::Lvector tLi = ti - (v1 * (2.0 / c1)) * (v1 * ti);
		/*5 */ BMesh::Lvector v2 = tj - tLi;
		/*6 */ double c2 = v2 * v2;
		/*7 */ BMesh::Lvector rj = rLi - (v2 * (2.0 / c2)) * (v2 * rLi);
		/*8 */ BMesh::Lvector sj = tj % rj;
		U.push_back(BMesh::Frame::fromST(sj.normalize(), tj.normalize()));
	}

	// RMF Visualization
	for (int i = 0; i < (int)in_points.size(); i++)
		U[i].center = BMesh::Point3D(in_points[i].x, in_points[i].y, in_points[i].z);


	//getResult();

}

void Rmf::getResult() {
	//if (getData()->edges.size() == 0) {
	//	getData()->edges.push_back(std::vector<V3f>());
	//}
	//getData()->edges.push_back(std::vector<V3f>());
	//for (int i = 0; i < U.size(); i++) {
	//	getData()->edges[1].push_back(U[i].center);
	//	getData()->edges[1].push_back(U[i].center + U[i].t);
	//	getData()->edges[1].push_back(U[i].center);
	//	getData()->edges[1].push_back(U[i].center + U[i].r);
	//	getData()->edges[1].push_back(U[i].center);
	//	getData()->edges[1].push_back(U[i].center + U[i].s);
	//}
}