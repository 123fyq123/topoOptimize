#pragma once
#include "BMesh/local_vector.hpp"
namespace BMesh {
	class Frame {
	public:
		BMesh::Lvector r, s, t;
		BMesh::Point3D center; // optional

		Frame() { r = BMesh::Lvector(1, 0, 0); s = BMesh::Lvector(0, 1, 0); t = BMesh::Lvector(0, 0, 1); }
		Frame(const BMesh::Lvector& R, const BMesh::Lvector& S, const BMesh::Lvector& T) { r = R; s = S; t = T; normalize(); }

		static Frame fromTR( BMesh::Lvector& T, BMesh::Lvector& R) { return Frame(R, T % R, T); }
		static Frame fromTR(BMesh::Lvector& T, const BMesh::Lvector& R) { return Frame(R, T % R, T); }
		static Frame fromRS(BMesh::Lvector& R, BMesh::Lvector& S) { return Frame(R, S, R % S); }
		static Frame fromST(BMesh::Lvector& S, BMesh::Lvector& T) { return Frame(S % T, S, T); }
		static Frame fromST(const BMesh::Lvector& S, const BMesh::Lvector& T) { return Frame(S % T, S, T); }
		static Frame fromT(BMesh::Lvector& T) {
			BMesh::Lvector Temp = T;
			Temp = Temp.normalize();
			BMesh::Lvector R = orthogonalVector(Temp).normalize();
			return fromTR(T, R);
		}

		static BMesh::Lvector orthogonalVector(const BMesh::Lvector& n) {
			if ((abs(n.y) >= 0.9 * abs(n.x)) && abs(n.z) >= 0.9 * abs(n.x)) return BMesh::Lvector(0.0, -n.z, n.y);
			else if (abs(n.x) >= 0.9 * abs(n.y) && abs(n.z) >= 0.9 * abs(n.y)) return BMesh::Lvector(-n.z, 0.0, n.x);
			else return BMesh::Lvector(-n.y, n.x, 0.0);
		}

		void normalize() { r = r.normalize(); s = s.normalize(); t = t.normalize(); };
	};
}
