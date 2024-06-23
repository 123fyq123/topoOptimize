#include "Q-MAT.h"

void Q_MAT::simplification(int target_vertices) {

	if (target_vertices >= mesh.vsize()) {
		std::cout << "Q-MAT: It is no need to simplification.\n";
		return;
	}
	int iter_num = mesh.vsize() - target_vertices;

	v_sqe = std::vector<V_SEQ>(mesh.vsize());
	v_collapse = std::vector<std::vector<Skl_CollapseNode*>>(mesh.vsize(), std::vector<Skl_CollapseNode*>{});

	calc_total_sqe();

	for (auto& ep : mesh.alledges()) {
		calc_collapse_cost(ep.first);
	}

	while (iter_num > 0 && !pri_que.empty()) {

		auto* node = pri_que.top();
		pri_que.pop();
		while (node->is_valid == false && !pri_que.empty()) {
			node = pri_que.top();
			pri_que.pop();
		}
		if (node->is_valid == false) break;
		
		auto& edge = mesh.edges(node->eh);
		auto vh_keep = edge.vh1();
		auto vh_remove = edge.vh2();
		assert(mesh.isValid(vh_keep) && mesh.isValid(vh_remove));
		
		//std::cout << "Q-MAT: Cost = " << node->cost << std::endl;
		auto& v_keep = mesh.vertices(vh_keep);
		auto& v_remove = mesh.vertices(vh_remove);
		
		// 为无效的收缩结点做标记
		for (auto* _node : v_collapse[vh_keep]) {
			_node->is_valid = false;
		}
		for (auto* _node : v_collapse[vh_remove]) {
			_node->is_valid = false;
		}
		v_collapse[vh_keep].clear();
		v_collapse[vh_remove].clear();

		// 更新结点几何
		v_keep.setPosition(node->sphere.x(), node->sphere.y(), node->sphere.z());
		mesh.radius[vh_keep] = node->sphere.w();
		v_sqe[vh_keep].A += v_sqe[vh_remove].A;
		v_sqe[vh_keep].b += v_sqe[vh_remove].b;
		v_sqe[vh_keep].c += v_sqe[vh_remove].c;

		std::vector<std::vector<MeshKernel::iGameVertexHandle>> new_triangles;
		for (auto& eh : mesh.NeighborEh(vh_remove)) {
			auto& adje = mesh.edges(eh);
			VH vh_other(adje.vh1() + adje.vh2() - vh_remove);
			if (vh_other == vh_keep || mesh.isConnected(vh_other, vh_keep)) continue;
			mesh.AddEdge(vh_keep, vh_other);
		}
		for (auto& fh : mesh.NeighborFh(vh_remove)) {
			auto vex = (mesh.faces(fh)).getVertexHandle();
			assert(vex.size() > 2);
			if (std::find(vex.begin(), vex.end(), vh_keep) != vex.end()) continue;
			for (int i = 0; i < vex.size(); ++i) {
				if (vex[i] == vh_remove) {
					vex[i] = vh_keep;
					new_triangles.push_back(vex);
					break;
				}
			}
		}
		mesh.DeleteVertex(vh_remove);

		// 添加三角形
		for (auto& tri : new_triangles) {
			mesh.AddFace(tri);
		}
		
		// 重新计算收缩代价
		for (auto& eh : mesh.NeighborEh(vh_keep)) {
			calc_collapse_cost(eh);
		}
		iter_num--;
	}
	mesh.updateAllHandles();
	//mesh.hole_to_triangle();
	std::cout << "Q-MAT: Simplification success. #V = " << mesh.vsize() << ", #E = " << mesh.esize() << ", #F = " << mesh.fsize() << std::endl;

}

void Q_MAT::calc_total_sqe() {
	
	int f_cnt = 0;
	for (auto& fp : mesh.allfaces()) {
		calc_face_sqe(fp.first);
		f_cnt++;
	}

	int e_cnt = 0;
	for (auto& ep : mesh.alledges()) {
		if (mesh.NeighborFh(ep.first).size() < 2) {
			calc_edge_sqe(ep.first);
			e_cnt++;
		}
	}

	int v_cnt = 0;
	for (auto& vp : mesh.allvertices()) {
		if (mesh.NeighborEh(vp.first).size() == 1) {
			calc_vertex_sqe(vp.first);
			v_cnt++;
		}
	}

	std::cout << "Q-MAT: Initial total SQE success. Total calculate #F = " << f_cnt << ", #E = " << e_cnt << ", #V = " << v_cnt << std::endl;

}

void Q_MAT::calc_face_sqe(FH fh) {
	const auto& face = mesh.faces(fh);
	const auto& vhs = face.getVertexHandle();
	auto& v0 = mesh.vertices(vhs[0]);
	auto& v1 = mesh.vertices(vhs[1]);
	auto& v2 = mesh.vertices(vhs[2]);
	//// 方案 1：计算三球的外公切面，可能不正确
	//auto tp = get_tangent_plane(v0, v1, v2);
	//Eigen::Vector4d n0(tp.planes[0].A, tp.planes[0].B, tp.planes[0].C, 1);
	//Eigen::Vector4d n1(tp.planes[1].A, tp.planes[1].B, tp.planes[1].C, 1);
	// 方案 2：拿球心形成的平面通过偏移半径来近似代替
	Vec vec01 = v1 - v0;
	Vec vec02 = v2 - v0;
	Vec normal_triangle = vec01.cross(vec02);
	normal_triangle.normalize();
	Vec p0 = v0 + normal_triangle * mesh.radius[vhs[0]];
	Vec p1 = v1 + normal_triangle * mesh.radius[vhs[1]];
	Vec p2 = v2 + normal_triangle * mesh.radius[vhs[2]];
	Vec _n0 = (p1 - p0).cross(p2 - p0);
	_n0.normalize();
	Eigen::Vector4d n0(_n0.x(), _n0.y(), _n0.z(), 1);
	p0 = v0 - normal_triangle * mesh.radius[vhs[0]];
	p1 = v1 - normal_triangle * mesh.radius[vhs[1]];
	p2 = v2 - normal_triangle * mesh.radius[vhs[2]];
	Vec _n1 = (p1 - p0).cross(p2 - p0);
	_n1.normalize();
	Eigen::Vector4d n1(_n1.x(), _n1.y(), _n1.z(), 1);

	Eigen::Matrix4d A = n0 * n0.transpose() + n1 * n1.transpose();
	for (auto& vh : vhs) {
		auto& v = mesh.vertices(vh);
		Eigen::Vector4d m(v.x(), v.y(), v.z(), 1);
		Eigen::Vector4d b = -2 * A * m;
		double c = m.transpose() * A * m;
		v_sqe[vh].A += A;
		v_sqe[vh].b += b;
		v_sqe[vh].c += c;
	}

}

void Q_MAT::calc_edge_sqe(EH eh) {
	const auto& edge = mesh.edges(eh);
	double factor = soft_k * get_edge_stability(edge);
	const auto& e_fhs = mesh.NeighborFh(eh);
	int degree = e_fhs.size();
	VH vh1 = edge.vh1(), vh2 = edge.vh2();
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);
	if (degree == 0) {

		Vec edge_dir = v2 - v1;// 向外
		edge_dir.normalize();
		// 求法向量为 dir 的平面, 平面方程 ax + by + cz + d = 0，edge_dir = (a, b, c)
		if (std::abs(edge_dir.z()) <= 1E-8F) edge_dir.z() = 1E-8F;
		double from_d = -(edge_dir.dot(v1));
		double point_z = (-from_d) / edge_dir.z();// 平面上的某个点(0, 0, point_z)
		Vec tp_vec0(-v1.x(), -v1.y(), point_z - v1.z());// 切平面上的某个向量
		Vec tp_vec1 = tp_vec0.cross(edge_dir);
		tp_vec0.normalize();
		tp_vec1.normalize();
		std::vector<std::vector<int>> move_dirs = {// 在切平面上移动的方向
			{-1, -1}, {-1, 1}, {1, -1}, {1, 1}
		};
		std::vector<Vec> from_points, to_points;// 在端点切平面上的点，八个顶点形成一个棱柱
		for (auto& move_dir : move_dirs) {
			Vec move_vec0 = tp_vec0 * move_dir[0];
			Vec move_vec1 = tp_vec1 * move_dir[1];
			from_points.push_back(v1 + move_vec0 * mesh.radius[vh1] + move_vec1 * mesh.radius[vh1]);
			to_points.push_back(v2 + move_vec0 * mesh.radius[vh2] + move_vec1 * mesh.radius[vh2]);
		}
		for (int i = 0; i < 4; ++i) {// 外移，形成大棱柱
			from_points[i] -= (edge_dir * mesh.radius[vh1]);
			to_points[i] += (edge_dir * mesh.radius[vh2]);
		}
		for (int i = 0; i < 4; ++i) {// 施加棱柱的四个表面约束
			Vec plane_p0 = from_points[i];
			Vec plane_p1 = to_points[i];
			Vec plane_p2 = to_points[(i + 1) % 4];
			Vec plane_normal = (plane_p1 - plane_p0).cross(plane_p2 - plane_p0);
			plane_normal.normalize();
			Eigen::Vector4d n(plane_normal.x(), plane_normal.y(), plane_normal.z(), 1);
			Eigen::Matrix4d A = n * n.transpose() * factor;
			Eigen::Vector4d m(v1.x(), v1.y(), v1.z(), 1);
			Eigen::Vector4d b = -2 * A * m;
			double c = m.transpose() * A * m;
			v_sqe[vh1].A += A;
			v_sqe[vh1].b += b;
			v_sqe[vh1].c += c;
			m = Eigen::Vector4d(v2.x(), v2.y(), v2.z(), 1);
			b = -2 * A * m;
			c = m.transpose() * A * m;
			v_sqe[vh2].A += A;
			v_sqe[vh2].b += b;
			v_sqe[vh2].c += c;
		}
		
	} else if (degree == 1) {
		// 需增加边界边软约束
		auto fh = *(e_fhs.begin());
		auto& face = mesh.faces(fh);
		const auto& f_vhs = face.getVertexHandle();
		Vec vec_e = v1 - v2;
		vec_e.normalize();
		Vec nor_f = (mesh.vertices(f_vhs[1]) - mesh.vertices(f_vhs[0])).cross(
			mesh.vertices(f_vhs[2]) - mesh.vertices(f_vhs[0]));
		nor_f.normalize();
		Vec move_dir = vec_e.cross(nor_f);
		move_dir.normalize();
		Vec edge_mid = (v1 + v2) * 0.5f;
		Vec centroid = mesh.getFaceCenter(fh);
		Vec toward_c = centroid - edge_mid;
		toward_c.normalize();
		if (move_dir.dot(toward_c) > 0.f) move_dir *= -1;
		Vec p0 = v1 + nor_f * mesh.radius[vh1];
		Vec p1 = v2 + nor_f * mesh.radius[vh2];
		Vec p2 = v2 - nor_f * mesh.radius[vh2];
		p0 += (move_dir * mesh.radius[vh1]);
		p1 += (move_dir * mesh.radius[vh2]);
		p2 += (move_dir * mesh.radius[vh2]);
		Vec plane_normal = (p1 - p0).cross(p2 - p0);
		plane_normal.normalize();
		Eigen::Vector4d n(plane_normal.x(), plane_normal.y(), plane_normal.z(), 1);
		Eigen::Matrix4d A = n * n.transpose() * factor;

		Eigen::Vector4d m(v1.x(), v1.y(), v1.z(), 1);
		Eigen::Vector4d b = -2 * A * m;
		double c = m.transpose() * A * m;
		v_sqe[vh1].A += A;
		v_sqe[vh1].b += b;
		v_sqe[vh1].c += c;
		m = Eigen::Vector4d(v2.x(), v2.y(), v2.z(), 1);
		b = -2 * A * m;
		c = m.transpose() * A * m;
		v_sqe[vh2].A += A;
		v_sqe[vh2].b += b;
		v_sqe[vh2].c += c;
	}
}

void Q_MAT::calc_vertex_sqe(VH vh) {

	EH eh(*(mesh.NeighborEh(vh).begin()));
	auto& edge = mesh.edges(eh);
	VH vh_o(edge.vh1() + edge.vh2() - vh);
	auto& v = mesh.vertices(vh);
	auto& v_o = mesh.vertices(vh_o);
	Vec edge_dir = (v - v_o).normalize();
	double factor = soft_k * get_edge_stability(edge);
	
	Eigen::Vector4d n(edge_dir.x(), edge_dir.y(), edge_dir.z(), 1);
	Eigen::Matrix4d A = n * n.transpose() * factor;
	Eigen::Vector4d m(v.x(), v.y(), v.z(), 1);
	Eigen::Vector4d b = -2 * A * m;
	double c = m.transpose() * A * m;
	v_sqe[vh].A += A;
	v_sqe[vh].b += b;
	v_sqe[vh].c += c;

}


void Q_MAT::calc_collapse_cost(EH eh) {

	auto& edge = mesh.edges(eh);
	auto& sqe0 = v_sqe[edge.vh1()];
	auto& sqe1 = v_sqe[edge.vh2()];
	Eigen::Matrix4d A = sqe0.A + sqe1.A;
	Eigen::Vector4d b = sqe0.b + sqe1.b;
	double c = sqe0.c + sqe1.c;// not be affected by sphere_x
	double cost = 99999999999;
	Eigen::Vector4d X;
	if (A.determinant()) {
		X = (-0.5) * (A.inverse() * b);
		cost = X.transpose() * A * X + b.dot(X);
	} else {
		for (double ratio = 0; ratio < 1.01; ratio += 0.5) {
			Vec _pos = mesh.vertices(edge.vh1()) * ratio + mesh.vertices(edge.vh2()) * (1 - ratio);
			Eigen::Vector3d pos(_pos.x(), _pos.y(), _pos.z());
			double r = mesh.radius[edge.vh1()] * ratio + mesh.radius[edge.vh2()] * (1 - ratio);
			Eigen::Vector4d m = Eigen::Vector4d(pos.x(), pos.y(), pos.z(), r);
			double tmp_cost = m.transpose() * A * m + b.dot(m);
			if (tmp_cost < cost) {
				cost = tmp_cost;
				X = m;
			}
		}
	}
	Skl_CollapseNode* node = new Skl_CollapseNode();
	double stability_weight = get_edge_stability(edge);
	stability_weight *= stability_weight;
	//std::cout << "vh1.r = " << mesh.radius[edge.vh1()] << ", vh2.r = " << mesh.radius[edge.vh2()] << "; Solved r = " << X.w() << ", ";
	X.w() = get_radius_by_interpolation(Vec(X.x(), X.y(), X.z()), edge.vh1(), edge.vh2());
	//std::cout << "Interpolation r = " << X.w() << std::endl;
	node->cost = (cost + c + hard_k) * stability_weight;
	if(std::isnan(node->cost)) return;
	//std::cout << "Geometric weight: " << cost + c << std::endl;
	//std::cout << "Old radius: " << edge->v0->radius << ", " << edge->v1->radius << "; new eadius: " << X.w() << std::endl;
	node->eh = eh;
	node->sphere = X;
	pri_que.push(node);
	v_collapse[edge.vh1()].push_back(node);
	v_collapse[edge.vh2()].push_back(node);

}

double Q_MAT::get_edge_stability(const MeshKernel::iGameEdge& edge) {
	auto& vh1 = edge.vh1();
	auto& vh2 = edge.vh2();
	auto& v1 = mesh.vertices(vh1);
	auto& v2 = mesh.vertices(vh2);
	double diff_radius = std::abs(mesh.radius[vh1] - mesh.radius[vh2]);
	double dis_center = (v1 - v2).norm();
	double dis_sphere = std::max((double)0, dis_center - diff_radius);
	double stability_ratio = dis_sphere / dis_center;
	//std::cout << "Stability: " << stability_ratio << std::endl;
	return stability_ratio;
}
//
////TwoPlane Q_MAT::get_tangent_plane(Skl_Vertex* v1, Skl_Vertex* v2, Skl_Vertex* v3) {
////	// 参考文献：三球的外公切面，文章编号：1673- 260X( 2007) 04- 0001- 02
////	TwoPlane res;
////	double a1 = v1->x(), b1 = v1->y(), c1 = v1->z(), r1 = v1->radius;
////	double a2 = v2->x(), b2 = v2->y(), c2 = v2->z(), r2 = v2->radius;
////	double a3 = v3->x(), b3 = v3->y(), c3 = v3->z(), r3 = v3->radius;
////	Eigen::Matrix3d matA, matB, matC, matA_, matB_, matC_;
////	matA << 1, 1, 1, a1, a2, a3, r1, r2, r3;
////	matB << 1, 1, 1, b1, b2, b3, r1, r2, r3;
////	matC << 1, 1, 1, c1, c2, c3, r1, r2, r3;
////	matA_ << 1, 1, 1, b1, b2, b3, c1, c2, c3;
////	matB_ << 1, 1, 1, a1, a2, a3, c1, c2, c3;
////	matC_ << 1, 1, 1, c1, c2, c3, b1, b2, b3;
////	double detA = matA.determinant(), detB = matB.determinant(), detC = matC.determinant();
////	double detA_ = matA_.determinant(), detB_ = matB_.determinant(), detC_ = matC_.determinant();
////	/*std::cout << "Determinant: \n";
////	std::cout << detA << std::endl;
////	std::cout << detB << std::endl;
////	std::cout << detC << std::endl;
////	std::cout << detA_ << std::endl;
////	std::cout << detB_ << std::endl;
////	std::cout << detC_ << std::endl;*/
////	// calculate lambdas
////	double num = detA_ * detA_ + detB_ * detB_ + detC_ * detC_ + detA * detA - detB * detB - detC * detC;
////	double sqrt_num = std::sqrt(num);
////	double B_C_ = detB_ * detC_;
////	double BC = detB * detC;
////	double denum = detA * detA + detC * detC - detB_ * detB_;
////	std::vector<double> lambdas(2);
////	lambdas[0] = (B_C_ - BC + detA * sqrt_num) / denum;
////	lambdas[1] = (B_C_ - BC - detA * sqrt_num) / denum;
////	// calculate planes
////	for (auto& lambda : lambdas) {
////		double A = detB + lambda * detC;
////		double B = -detA;
////		double C = -detA * lambda;
////		double D = A * (-a1) + B * (-b1) + C * (-c1) - r1 * detC_ - lambda * r1 * detB_;
////		double len = std::sqrt(A * A + B * B + C * C);
////		A /= len;
////		B /= len;
////		C /= len;
////		D /= len;
////		Plane p(A, B, C, D);
////		res.planes.push_back(p);
////	}
////	return res;
////}
//
double Q_MAT::get_radius_by_interpolation(const MeshKernel::iGameVertex& v_n, VH vh0, VH vh1) {
	auto& v0 = mesh.vertices(vh0);
	auto& v1 = mesh.vertices(vh1);
	double radius = 0.f;
	double dis_0 = (v_n - v0).norm();
	double dis_1 = (v_n - v1).norm();
	double weight_1 = dis_0 / (dis_0 + dis_1);
	radius = weight_1 * mesh.radius[vh1] + (1 - weight_1) * mesh.radius[vh0];

	return radius;

}