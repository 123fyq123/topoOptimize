#include "Denoising_BNF.h"

void Denoising_BNF::Execute() {

	
	initPositions();
	initTrianglesData();
	initSigmaCenter();
	for (int i = 0; i < it_normal; ++i) {
		BilateralNormalFiltering();
		//printf("BNF %d success\n", i);
	}
	for (int i = 0; i < it_vertex; ++i) {
		VertexUpdate();
		//printf("VexUpdate %d success\n", i);
	}
	for (auto& vp : mesh.allvertices()) {
		MeshKernel::iGameVertex& vex = mesh.vertices(vp.first);
		vex.setPosition(positions[vp.first][0], positions[vp.first][1], positions[vp.first][2]);
	}
	
	printf("denoising success\n");
}

void Denoising_BNF::BilateralNormalFiltering() {
	auto fcnt = mesh.FaceSize();
	float double_SigmaCenter2 = 2.f * SigmaCenter * SigmaCenter;
	float double_SigmaNormal2 = 2.f * SigmaNormal * SigmaNormal;
	std::vector<Eigen::Vector3f> new_normals;
	new_normals.resize(fcnt, Eigen::Vector3f::Zero());
	for (auto& fp : mesh.allfaces()) {
		int fi = fp.first;
		auto adjF = mesh.Neighbor2Fh(fp.first);// share common vertex
		float sum_weight = 0.f;
		for (auto& fj : adjF) {
			float delta_center = (centers[fi] - centers[fj]).norm();
			float delta_normal = (normals[fi] - normals[fj]).norm();
			float weight_center = std::exp(-(delta_center * delta_center) / double_SigmaCenter2);
			float weight_normal = std::exp(-(delta_normal * delta_normal) / double_SigmaNormal2);
			float wj = areas[fj] * weight_center * weight_normal;
			new_normals[fi] += (wj * normals[fj]);
			sum_weight += wj;
		}
		assert(sum_weight > static_cast<float>(0));
		new_normals[fi] /= sum_weight;
		new_normals[fi].normalize();
	}
	normals = new_normals;
}

void Denoising_BNF::VertexUpdate() {
	for (auto& vp : mesh.allvertices()) {
		auto& vh = vp.first;
		Eigen::Vector3f move = Eigen::Vector3f::Zero();
		auto adj_faces = mesh.NeighborFh(vh);
		for (auto& fh : adj_faces) {
			float wj = normals[fh].transpose() * (centers[fh] - positions[vh]);
			move += (wj * normals[fh]);
		}
		move /= adj_faces.size();
		positions[vh] += move;
	}
}

void Denoising_BNF::initSigmaCenter() {
	SigmaCenter = 0.f;
	int cnt = 0;
	for (auto& fp : mesh.allfaces()) {
		auto fhs = mesh.NeighborFh(fp.first);// share common edge
		for (auto fh : fhs) {
			SigmaCenter += (centers[fp.first] - centers[fh]).norm();
			cnt++;
		}
	}
	SigmaCenter /= cnt;
	printf("SigmaCenter %.2f; SigmaNormal: %.2f\n", SigmaCenter, SigmaNormal);
}

void Denoising_BNF::initPositions() {
	int vcnt = mesh.VertexSize();
	positions.resize(vcnt);
	for (auto& vp : mesh.allvertices()) {
		positions[vp.first] = Eigen::Vector3f(vp.second.x(), vp.second.y(), vp.second.z());
	}
}

void Denoising_BNF::initTrianglesData() {
	int fcnt = mesh.FaceSize();
	centers.resize(fcnt, Eigen::Vector3f::Zero());
	normals.resize(fcnt);
	areas.resize(fcnt);
	for (auto& fp : mesh.allfaces()) {
		int f_idx = fp.first;
		auto v_indices = fp.second.getVertexHandle();
		assert(v_indices.size() == 3 && "isn't a triangle");
		Eigen::Vector3f vec23 = positions[v_indices[2]] - positions[v_indices[1]];
		Eigen::Vector3f vec21 = positions[v_indices[0]] - positions[v_indices[1]];

		for (auto& v_idx : v_indices) {
			centers[f_idx] += positions[v_idx];
		}
		centers[f_idx] /= v_indices.size();

		// 计算法向量 注意方向
		normals[f_idx] = (vec23.cross(vec21)).normalized();

		float len1 = vec23.norm();
		float len2 = vec21.norm();
		vec23.normalize();
		vec21.normalize();
		float cos_theta = vec23.dot(vec21);
		float sin_theta = std::sqrt(1 - cos_theta * cos_theta);
		areas[f_idx] = len1 * len2 * sin_theta * 0.5;

	}
	printf("init triangles data success\n");
}