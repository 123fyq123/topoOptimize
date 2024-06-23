#include "Analysis_HexMesh.h"

void Analysis_HexMesh::calculateCellQuality_Jacobian(std::vector<double>& scaled_jacobian) {

	int ccnt = mesh.csize();
	scaled_jacobian.reserve(ccnt);
	auto& cells = mesh.allcells();

	for (auto& cp : cells) {
		CH ch = cp.first;
		auto& c = mesh.cells(ch);
		double min_det = 2.f;// the value of Scaled Jacobian is the	determinant of the matrix Jacobian
		std::unordered_map<int, std::vector<int>> neighbor;
		const auto& _ehs = c.getEdgeHandle();
		for (auto& _eh : _ehs) {
			auto& _e = mesh.edges(_eh);
			auto& vh1 = _e.vh1();
			auto& vh2 = _e.vh2();
			neighbor[vh1].push_back(vh2);
			neighbor[vh2].push_back(vh1);
		}
		const auto& _vhs = c.getVertexHandle();
		for (auto& _vh : _vhs) {
			auto& v = mesh.vertices(_vh);
			auto adjvh = neighbor[_vh];
			assert(adjvh.size() == 3);
			auto v1 = mesh.vertices(VH(adjvh[0]));
			auto v2 = mesh.vertices(VH(adjvh[1]));
			auto v3 = mesh.vertices(VH(adjvh[2]));
			auto center = (v1 + v2 + v3) / 3;
			auto vec = (center - v).normalized();

			// ajust oridering
			auto vec12 = v2 - v1, vec13 = v3 - v1;
			auto normal = (vec12 % vec13).normalized();
			double cosine = vec * normal;
			if (cosine < 0) std::swap(v2, v3);

			std::vector<Eigen::Vector3d> ev(3);
			auto tmp = v1 - v;
			ev[0] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v2 - v;
			ev[1] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();
			tmp = v3 - v;
			ev[2] = (Eigen::Vector3d(tmp.x(), tmp.y(), tmp.z())).normalized();

			Eigen::Matrix3d J;
			J << ev[0][0], ev[1][0], ev[2][0],
				ev[0][1], ev[1][1], ev[2][1],
				ev[0][2], ev[1][2], ev[2][2];
			//std::cout << J << endl;
			/*
			For an element, the full range of the Scaled Jacobian value is from −1 to + 1.
			And fscaled_jacobian = 1 if the element is an ideal element,
			fscaled_jacobian = −1 if the element is a worst distorted element
			*/
			min_det = std::min(min_det, J.determinant());
		}

		scaled_jacobian.push_back(min_det);

	}

}