#include "Smoothing.h"

HexMesh_Smoothing::HexMesh_Smoothing(MeshKernel::VolumeMesh& _mesh) : mesh(_mesh) {

}

void HexMesh_Smoothing::Laplacian_Smoothing_HC(std::vector<int> feature_vhs, int iter_max) {

	double alpha = 0, beta = 0.5;
	double volume_origin = MeshMath::get_volume_hexahedral_mesh(mesh);
	int n_vertices = mesh.vsize();
	exterior_vhs.clear();
	interior_vhs.clear();
	vector<iGameVertex> vertices_origin;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) exterior_vhs.push_back(vp.first);
		else interior_vhs.push_back(vp.first);
		vertices_origin.push_back(vp.second);
	}

	int exterior_cnt = exterior_vhs.size(), interior_cnt = interior_vhs.size();

	for (int it = 0; it < iter_max; ++it) {

		vector<iGameVertex> differ(n_vertices);
		vector<iGameVertex> vertices_smoothed(n_vertices);

#pragma omp parallel for
		for (int i = 0; i < exterior_cnt; ++i) {
			auto vh = exterior_vhs[i];
			int ex_cnt = 0;
			iGameVertex adj_v(0, 0, 0);
			for (auto& adjvh : mesh.NeighborVh(vh)) {
				if (mesh.isOnBoundary(adjvh)) {
					adj_v += mesh.vertices(adjvh);
					ex_cnt++;
				}
			}
			vertices_smoothed[vh] = adj_v / ex_cnt;
			differ[vh] = vertices_smoothed[vh] - (vertices_origin[vh] * alpha + mesh.vertices(vh) * (1 - alpha));
		}

#pragma omp parallel for
		for (int i = 0; i < exterior_cnt; ++i) {
			auto vh = exterior_vhs[i];
			int ex_cnt = 0;
			iGameVertex adj_diff(0, 0, 0);
			for (auto& adjvh : mesh.NeighborVh(vh)) {
				if (mesh.isOnBoundary(adjvh)) {
					adj_diff += differ[adjvh];
					ex_cnt++;
				}
			}
			adj_diff /= ex_cnt;
			vertices_smoothed[vh] = vertices_smoothed[vh] - (differ[vh] * beta + adj_diff * (1 - beta));
			auto& v = mesh.vertices(vh);
			v = vertices_smoothed[vh];
		}

#pragma omp parallel for
		for (int i = 0; i < interior_cnt; ++i) {
			auto vh = interior_vhs[i];
			const auto& adjvhs = mesh.NeighborVh(vh);
			iGameVertex adj_v(0, 0, 0);
			for (auto& adjvh : adjvhs) {
				adj_v += mesh.vertices(adjvh);
			}
			vertices_smoothed[vh] = adj_v / adjvhs.size();
		}

#pragma omp parallel for
		for (int i = 0; i < interior_cnt; ++i) {
			auto vh = interior_vhs[i];
			auto& v = mesh.vertices(vh);
			v = vertices_smoothed[vh];
		}

	}

	double volume_current = MeshMath::get_volume_hexahedral_mesh(mesh);

	std::cout << "Smooth(HC): Volume_origin = " << volume_origin << ", Volume_current = " << volume_current << std::endl;

}

void HexMesh_Smoothing::Laplacian_Smoothing_Diffusion(int max_iter, std::vector<int> feature_cells, std::vector<int> feature_faces,
	std::vector<int> feature_edges, std::vector<int> feature_vertices) {

	double original_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
	int n_vertices = mesh.vsize();
	vector<bool> is_feature(n_vertices, false);
	for (auto& ch : feature_cells) {
		iGameCellHandle _ch(ch);
		auto& cell = mesh.cells(_ch);
		for (auto& vh : cell.getVertexHandle()) {
			is_feature[vh] = true;
		}
	}
	for (auto& fh : feature_faces) {
		iGameFaceHandle _fh(fh);
		auto& face = mesh.faces(_fh);
		for (auto& vh : face.getVertexHandle()) {
			is_feature[vh] = true;
		}
	}
	for (auto& eh : feature_edges) {
		iGameEdgeHandle _eh(eh);
		auto& edge = mesh.edges(_eh);
		is_feature[edge.vh1()] = true;
		is_feature[edge.vh2()] = true;
	}
	for (auto& vh : feature_vertices) {
		is_feature[iGameVertexHandle(vh)] = true;
	}

	vector<double> weight_s(n_vertices, 1.f);// 采用渐进扩散，越靠近特征顶点该权重越小
	vector<bool> visited(n_vertices, false);
	queue<int> que;
	int feature_cnt = 0;
	for (int i = 0; i < n_vertices; ++i) {
		if (is_feature[i]) {
			feature_cnt++;
			visited[i] = true;
			weight_s[i] = 0.f;
			que.push(i);
		}
	}
	while (!que.empty()) {
		auto vid = que.front();
		que.pop();
		MeshKernel::iGameVertexHandle vh(vid);
		double weight_cur = weight_s[vid];
		for (auto& adjvh : mesh.NeighborVh(vh)) {
			if (visited[adjvh]) continue;
			visited[adjvh] = true;
			weight_s[adjvh] = std::min((double)1.f, weight_cur + 0.01f);
			if (weight_s[adjvh] < 1.f - (1E-4F)) {
				que.push(adjvh);
			}
		}
	}
	/*for (int i = 0; i < n_vertices; ++i) {
		if (weight_s[i] < 1.f)
			std::cout << "#V " << i << ": " << weight_s[i] << std::endl;
	}*/


	unordered_map<int, bool> is_on_boundary;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			exterior_vhs.push_back(vp.first);
			is_on_boundary[vp.first] = true;
		} else {
			interior_vhs.push_back(vp.first);
		}
	}
	int exterior_vcnt = exterior_vhs.size(), interior_vcnt = interior_vhs.size();

	for (int it = 0; it < max_iter; ++it) {

		vector<iGameVertex> total_vertices;
		vector<iGameVertex> positions(exterior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			auto& vh = exterior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			int vcnt = 0;
			for (auto& adjvh : adjvhs) {
				if (is_on_boundary.count(adjvh)) {
					pos = pos + mesh.vertices(adjvh);
					vcnt++;
				}
			}
			positions[i] = pos / vcnt;
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			if (is_feature[exterior_vhs[i]]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(exterior_vhs[i]);
			v = v * (1 - weight_s[exterior_vhs[i]]) + positions[i] * weight_s[exterior_vhs[i]];
		}

		positions.resize(interior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			auto& vh = interior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			for (auto& adjvh : adjvhs) {
				pos = pos + mesh.vertices(adjvh);
			}
			positions[i] = pos / adjvhs.size();
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			if (is_feature[interior_vhs[i]]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(interior_vhs[i]);
			v = v * (1 - weight_s[interior_vhs[i]]) + positions[i] * weight_s[interior_vhs[i]];
		}

		// rescale
		double current_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
		double scale = std::pow(original_volume / current_volume, 1.f / 3);
#pragma omp parallel for
		for (int i = 0; i < mesh.vsize(); ++i) {
			if (is_feature[i]) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(MeshKernel::iGameVertexHandle(i));
			v *= scale;
		}

	}

	std::cout << "Smoothing success. Feature vertices size = " << feature_cnt << std::endl;
}

void HexMesh_Smoothing::Laplacian_Smoothing(int max_iter, std::vector<int> feature_cells, std::vector<int> feature_faces, 
	std::vector<int> feature_edges, std::vector<int> feature_vertices) {

	double volume_origin = MeshMath::get_volume_hexahedral_mesh(mesh);

	unordered_map<int, bool> feature_vhs;
	for (auto& ch : feature_cells) {
		iGameCellHandle _ch(ch);
		auto& cell = mesh.cells(_ch);
		for (auto& vh : cell.getVertexHandle()) {
			feature_vhs[vh] = true;
		}
	}
	for (auto& fh : feature_faces) {
		iGameFaceHandle _fh(fh);
		auto& face = mesh.faces(_fh);
		for (auto& vh : face.getVertexHandle()) {
			feature_vhs[vh] = true;
		}
	}
	for (auto& eh : feature_edges) {
		iGameEdgeHandle _eh(eh);
		auto& edge = mesh.edges(_eh);
		feature_vhs[edge.vh1()] = true;
		feature_vhs[edge.vh2()] = true;
	}
	for (auto& vh : feature_vertices) {
		feature_vhs[iGameVertexHandle(vh)] = true;
	}

	unordered_map<int, bool> is_on_boundary;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.isOnBoundary(vp.first)) {
			exterior_vhs.push_back(vp.first);
			is_on_boundary[vp.first] = true;
		} else {
			interior_vhs.push_back(vp.first);
		}
	}
	int exterior_vcnt = exterior_vhs.size(), interior_vcnt = interior_vhs.size();

	for (int it = 0; it < max_iter; ++it) {

		vector<iGameVertex> total_vertices;
		vector<iGameVertex> positions(exterior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			auto& vh = exterior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			int vcnt = 0;
			for (auto& adjvh : adjvhs) {
				if (is_on_boundary.count(adjvh)) {
					pos = pos + mesh.vertices(adjvh);
					vcnt++;
				}
			}
			positions[i] = pos / vcnt;
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < exterior_vcnt; ++i) {
			if (feature_vhs.count(exterior_vhs[i])) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(exterior_vhs[i]);
			v = positions[i];
		}

		positions.resize(interior_vcnt);

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			auto& vh = interior_vhs[i];
			iGameVertex pos(0, 0, 0);
			auto adjvhs = mesh.NeighborVh(vh);
			for (auto& adjvh : adjvhs) {
				pos = pos + mesh.vertices(adjvh);
			}
			positions[i] = pos / adjvhs.size();
		}
		total_vertices.insert(total_vertices.end(), positions.begin(), positions.end());

#pragma omp parallel for
		for (int i = 0; i < interior_vcnt; ++i) {
			if (feature_vhs.count(interior_vhs[i])) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(interior_vhs[i]);
			v = positions[i];
		}

		// rescale
		double current_volume = MeshMath::get_volume_hexahedral_mesh(mesh);
		double scale = std::pow(volume_origin / current_volume, 1.f / 3);
#pragma omp parallel for
		for (int i = 0; i < mesh.vsize(); ++i) {
			if (feature_vhs.count(i)) continue;// 特征顶点不修改位置
			auto& v = mesh.vertices(MeshKernel::iGameVertexHandle(i));
			v *= scale;
		}

	}
	double volume_current = MeshMath::get_volume_hexahedral_mesh(mesh);

	std::cout << "Smooth: Volume_origin = " << volume_origin << ", Volume_current = " << volume_current << std::endl;

	//std::cout << "Smoothing success. Feature vertices size = " << feature_vhs.size() << std::endl;
}

double HexMesh_Smoothing::CalculateCellQuality(MeshKernel::iGameCellHandle ch) {
	auto& c = mesh.cells(ch);
	double min_det = 2.f;// the value of Scaled Jacobian is the	determinant of the matrix Jacobian
	unordered_map<int, vector<int>> neighbor;
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
		auto v1 = mesh.vertices(iGameVertexHandle(adjvh[0]));
		auto v2 = mesh.vertices(iGameVertexHandle(adjvh[1]));
		auto v3 = mesh.vertices(iGameVertexHandle(adjvh[2]));
		auto center = (v1 + v2 + v3) / 3;
		auto vec = (center - v).normalized();

		// ajust oridering
		auto vec12 = v2 - v1, vec13 = v3 - v1;
		auto normal = (vec12 % vec13).normalized();
		double cosine = vec * normal;
		if (cosine < 0) swap(v2, v3);

		vector<Eigen::Vector3d> ev(3);
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
		min_det = min(min_det, J.determinant());
	}
	return min_det;
}