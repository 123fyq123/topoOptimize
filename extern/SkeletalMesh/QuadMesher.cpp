#include "QuadMesher.h"

bool QuadMesher::meshing(MeshKernel::SurfaceMesh& quadmesh) {

	if (!trimesh.isTriangleMesh() || !trimesh.hasLoopBoundary()) {
		std::cerr << "QuadMesher: Input mesh isn't triangle mesh or no loop boundary." << std::endl;
		return false;
	}

	for (auto& vp : trimesh.allvertices()) {
		if (std::abs(vp.second.z()) > DBL_EPSILON) {// 非二维
			std::cerr << "QuadMesher: Input mesh isn't a 2D plane." << std::endl;
			return false;
		}
	}

	perimeter = 0;
	for (auto& ep : trimesh.alledges()) {
		if (trimesh.isOnBoundary(ep.first)) {
			boundary_ehs.push_back(ep.first);
			auto& edge = trimesh.edges(ep.first);
			auto& v1 = trimesh.vertices(edge.vh1());
			auto& v2 = trimesh.vertices(edge.vh2());
			perimeter += (v1 - v2).norm();
		}
	}

	//smoothing_interior(5);
	calculate_cross_field();
	generate_initial_streamlines();
	simplify_streamlines();


}

void QuadMesher::quad_meshing() {



}

void QuadMesher::simplify_streamlines() {

	std::unordered_map<int, std::vector<int>> out_streamline;
	std::unordered_map<int, std::vector<int>> in_streamline;
	// 1. 参数化每条流线
	for (int i = 0; i < streamlines.size(); ++i) {
		auto& sl = streamlines[i];
		sl.init_length_ratio();
		out_streamline[sl.beg_sid].push_back(i);
		if (!sl.end_on_boundary) in_streamline[sl.end_sid].push_back(i);
	}

	// 2. 通过拓扑和几何来简化流线
	std::vector<bool> need_to_erase(streamlines.size(), false);// 记录旧流线是否可直接作为新流线使用
	std::vector<Streamline> streamlines_new;

	// 2.1 简化终点在边界且意义相同的流线
	for (int sid = 0; sid < singular_points.size(); ++sid) {
		auto& point = singular_points[sid];
		if (point.mode != 0) continue;
		std::vector<std::pair<double, int>> bd_sl;
		for (auto& slid : out_streamline[sid]) {
			if (streamlines[slid].end_on_boundary) {
				bd_sl.push_back(std::pair<double, int> { streamlines[slid].length.back(), slid });
			}
		}
		std::sort(bd_sl.begin(), bd_sl.end());
		for (int i = 1; i < bd_sl.size(); ++i) {
			if (bd_sl[i].first < bd_sl[i - 1].first * 2) {// 选择两条长度基本相等的边界流线简化
				int bid0 = bd_sl[i - 1].second, bid1 = bd_sl[i].second;
				need_to_erase[bid0] = true;
				need_to_erase[bid1] = true;
				// 两条旧流线 线性插值得到新的流线
				Streamline sl_new;
				get_streamline_itrp(streamlines[bd_sl[0].second], streamlines[bd_sl[1].second], sl_new);
				streamlines_new.push_back(sl_new);
				i++;
			}
		}
		
	}

	// 2.2 简化起点终点均为奇异点的流线
	for (int i = 0; i < streamlines.size(); ++i) {
		auto& sl0 = streamlines[i];
		if (need_to_erase[i] || sl0.end_on_boundary) continue;
		int beg_sid = sl0.beg_sid, end_sid = sl0.end_sid;
		
		std::vector<std::pair<double, int>> candiate_sid;

		for (int j = i + 1; j < streamlines.size(); ++j) {
			auto& sl1 = streamlines[j];
			if (need_to_erase[j] || sl1.end_on_boundary) continue;
			if (sl1.beg_sid == end_sid && sl1.end_sid == beg_sid) {// 找到两条意义相同, 方向相反的流线
				candiate_sid.push_back( {sl1.length.back(), j} );
			} else if (sl1.beg_sid == beg_sid && sl1.end_sid == end_sid) {// 找到两条意义相同, 方向相同的流线
				candiate_sid.push_back({ sl1.length.back(), j });
			}
			
		}

		if (!candiate_sid.empty()) {
			int sid1 = -1;
			double lendiff = 0x3f3f3f3f;
			need_to_erase[i] = true;
			for (auto& p : candiate_sid) {
				need_to_erase[p.second] = true;
				double tmp = std::abs(sl0.length.back() - p.first);
				if (tmp < lendiff) {
					lendiff = tmp;
					sid1 = p.second;
				}
			}
			// 两条旧流线 线性插值得到新的流线
			auto& sl1 = streamlines[sid1];

			if (sl1.beg_sid == sl0.end_sid) {// 反向
				auto& points = sl1.points;
				std::reverse(points.begin(), points.end());
				sl1.beg_sid = sl0.beg_sid, sl1.end_sid = sl0.end_sid;
				sl1.init_length_ratio();
			}

			Streamline sl_new;
			get_streamline_itrp(sl0, sl1, sl_new);
			streamlines_new.push_back(sl_new);
			//std::cout << "Itrp success: from " << beg_sid << " to " << end_sid << std::endl;
		}

	}
	
	for (int i = 0; i < streamlines.size(); ++i) {// 检查有无路过过其他奇异点的冗余边界线
		if (streamlines[i].end_on_boundary && !need_to_erase[i]) {
			double len_sl = streamlines[i].length.back();
			for (int sid = 0; sid < singular_points.size(); ++sid) {
				if (sid == streamlines[i].beg_sid) continue;
				double len_ = (singular_points[sid].pos - streamlines[i].points.back()).norm();
				if (len_ < len_sl) {
					need_to_erase[i] = true;
					break;
				}
			}
		}
	}

	for (int i = 0; i < streamlines.size(); ++i) {
		if (!need_to_erase[i]) {
			streamlines_new.push_back(streamlines[i]);
		}
	}

	streamlines = streamlines_new;

	for (int i = 0; i < streamlines.size(); ++i) {
		auto& sl = streamlines[i];
		std::cout << "Streamline #" << i << ": beg_vh=" << sl.beg_sid;
		if (sl.end_on_boundary) std::cout << " to boundary. ";
		else std::cout << " to end_vh=" << sl.end_sid << ". ";
		std::cout << "\tbeg_pos=(" << sl.points[0].x() << ", " << sl.points[0].y() << ", " << sl.points[0].z() << "), ";
		std::cout << "end_pos=(" << sl.points.back().x() << ", " << sl.points.back().y() << ", " << sl.points.back().z() << ").\n";
		
	}

}

void QuadMesher::generate_initial_streamlines() {

	for (auto& singular_point : singular_points) {

		if (singular_point.mode == 0) generate_initial_streamlines_v(singular_point);
		else if (singular_point.mode == 1) generate_initial_streamlines_e(singular_point);
		else if (singular_point.mode == 2) generate_initial_streamlines_f(singular_point);

	}

}

void QuadMesher::morch_on_streamline(const MorchNode& node, Streamline& streamline, std::set<EH>& visited) {

	EH eh = node.eh;
	if (trimesh.isOnBoundary(eh)) {// 终止条件1: 到达边界边
		streamline.points.push_back(node.pos);
		streamline.end_on_boundary = true;
		streamline.end_edgehandle = eh;
		return;
	}
	auto& edge = trimesh.edges(eh);
	if (vh_to_spid.count(edge.vh1()) && streamline.beg_sid != edge.vh1()) {// 终止条件2: 到达非起始奇异点的邻域边
		streamline.points.push_back(trimesh.vertices(edge.vh1()));
		streamline.end_sid = vh_to_spid[edge.vh1()];
		return;
	} else if (vh_to_spid.count(edge.vh2()) && streamline.beg_sid != edge.vh2()) {
		streamline.points.push_back(trimesh.vertices(edge.vh2()));
		streamline.end_sid = vh_to_spid[edge.vh2()];
		return;
	} else if (eh_to_spid.count(eh)) {// 到达奇异边
		int sid = eh_to_spid[eh];
		streamline.points.push_back(singular_points[sid].pos);
		streamline.end_sid = sid;
		return;
	} else if (fh_to_spid.count(node.pre_fh)) {// 经过奇异面
		int sid = fh_to_spid[node.pre_fh];
		streamline.points.push_back(singular_points[sid].pos);
		streamline.end_sid = sid;
		return;
	}

	if (visited.count(eh)) {// 终止条件3: 重复遍历一条边
		streamline.end_on_fault = true;
		return;
	}

	// 将当前交点放入流线中
	streamline.points.push_back(node.pos);
	visited.insert(eh);

	std::vector<Vec> crossdirs_p = get_cross_dirs(node.pos, node.eh, node.pre_fh);// 在前一个三角形面片上插值得到标架场
	Vec dir_fp = get_most_similar_dir(node.dir, crossdirs_p);// f(p_i)
	EH eh_next(-1);
	FH fh_next(-1);
	// 1. 先求得临时交点 p'_i+1
	Vex intersection;// p'_i+1
	if (!get_instersection_polygon(node.pos, dir_fp, node.eh, node.pre_fh, eh_next, fh_next, intersection)) {
		std::cerr << "QuadMesher: Wrong intersection count!" << std::endl;
		streamline.end_on_fault = true;
		return;
	}
	
	// 2. 使用 临时交点 与 源点的标架场分量算出 本次求交的最终方向向量 d_i+1
	Vec tmp_dirsp = (intersection - node.pos).normalize();
	std::vector<Vec> crossdirs_tp = get_cross_dirs(intersection, eh_next, fh_next);
	Vec dir_ftp = get_most_similar_dir(tmp_dirsp, crossdirs_tp);// f(p^'_i+1)
	Vec dir_next = (dir_ftp + dir_fp) * 0.5;
	dir_next.normalize();

	// 3. 算出本次迭代的最终交点 p_i+1
	if (!get_instersection_polygon(node.pos, dir_next, node.eh, node.pre_fh, eh_next, fh_next, intersection)) {
		std::cerr << "QuadMesher: Wrong intersection count!" << std::endl;
		streamline.end_on_fault = true;
		return;
	}


	MorchNode node_next;
	node_next.pos = intersection;
	node_next.dir = dir_next;// 到达交点时的方向向量
	node_next.eh = eh_next;// 交点所在的边
	node_next.pre_fh = fh_next;// 交线经过的面

	morch_on_streamline(node_next, streamline, visited);

}

bool QuadMesher::is_on_streamline(Vec source, EH eh, FH fh, Vec& march_dir, double& weight_vh1) {// 源点, 目标边, 所在面
	auto& edge = trimesh.edges(eh);
	const auto& fvhs = trimesh.faces(fh).getVertexHandle();
	VH vh1 = edge.vh1();
	VH vh2 = edge.vh2();
	auto& v1 = trimesh.vertices(vh1);
	auto& v2 = trimesh.vertices(vh2);
	
	int fvid1 = -1, fvid2 = -1;
	for (int i = 0; i < 3; ++i) {
		if (fvhs[i] == vh1) fvid1 = i;
		else if (fvhs[i] == vh2) fvid2 = i;
	}
	if (fvid1 == -1 || fvid2 == -1) return false;
	std::vector<Vec> cross_edge_dirs(2);
	std::vector<std::vector<Vec>> cross_dirs = {
		get_cross_dirs(crossfield[fh][fvid1]) , get_cross_dirs(crossfield[fh][fvid2])
	};
	for (int i = 0; i < 2; ++i) {
		Vec dire = (v2 - v1).normalize();
		double cos_max = -2;
		for (auto& dirv : cross_dirs[i]) {
			double _cos = dire.dot(dirv);
			if (_cos > cos_max) {
				cos_max = _cos;
				cross_edge_dirs[i] = dirv;
			}
		}
	}

	double cosine_max = -2;
	for (double w = 0.f; w <= 1.f; w += 0.001) {
		Vec point = v1 * w + v2 * (1 - w);
		Vec sp = (point - source).normalize();
		Vec dirp = cross_edge_dirs[0] * w + cross_edge_dirs[1] * (1 - w);
		dirp.normalize();
		std::vector<Vec> cross_p = get_cross_dirs(dirp);
		for (auto& dir : cross_p) {
			double cosine = dir.dot(sp);
			if (cosine > cosine_max) {
				cosine_max = cosine;
				weight_vh1 = w;
				march_dir = dir;
			}
		}
		if (cosine_max > 0.99) break;
	}
	//std::cout << "Cosine_max = " << cosine_max << std::endl;
	return cosine_max > 0.99;

}

bool QuadMesher::is_on_streamline(Vec source, VH vh, EH eh, Vec& march_dir) {
	const auto& adjfhs = trimesh.NeighborFh(eh);
	Vec cross_dirv(0, 0, 0);
	for (auto& fh : adjfhs) {
		const auto& fvhs = trimesh.faces(fh).getVertexHandle();
		int local_vid = std::find(fvhs.begin(), fvhs.end(), vh) - fvhs.begin();
		if (local_vid < 0 || local_vid > 2) return false;
		cross_dirv += crossfield[fh][local_vid];
	}
	cross_dirv /= adjfhs.size();
	cross_dirv.normalize();
	auto cross_dirvs = get_cross_dirs(cross_dirv);
	Vec sp = trimesh.vertices(vh) - source;
	sp.normalize();
	double cosine_max = -2;
	for (auto& dir : cross_dirvs) {
		double tmp_cos = sp.dot(dir);
		if (tmp_cos > cosine_max) {
			cosine_max = tmp_cos;
			march_dir = dir;
		} 
	}
	//std::cout << "Cosine_max = " << cosine_max << std::endl;
	return cosine_max > 0.99;
}


std::vector<Vec> QuadMesher::get_cross_dirs(const Vec& point, EH eh, FH fh) {// 返回 point 在边上插值得到的标架场 
	auto& edge = trimesh.edges(eh);
	auto& face = trimesh.faces(fh);
	auto vh1 = edge.vh1(), vh2 = edge.vh2();
	auto& v1 = trimesh.vertices(vh1);
	auto& v2 = trimesh.vertices(vh2);
	const auto& fvhs = face.getVertexHandle();
	int fvid1 = -1, fvid2 = -1;
	for (int i = 0; i < 3; ++i) {
		if (fvhs[i] == vh1) fvid1 = i;
		else if (fvhs[i] == vh2) fvid2 = i;
	}
	if (fvid1 == -1 || fvid2 == -1) {
		std::cerr << "QuadMesher: EH is not belong to FH.\n";
		return {};
	}
	std::vector<Vec> cross_edge_dirs(2);
	std::vector<std::vector<Vec>> cross_dirs = {
		get_cross_dirs(crossfield[fh][fvid1]) , get_cross_dirs(crossfield[fh][fvid2])
	};
	for (int i = 0; i < 2; ++i) {
		Vec dire = (v2 - v1).normalize();
		double cos_max = -2;
		for (auto& dirv : cross_dirs[i]) {
			double _cos = dire.dot(dirv);
			if (_cos > cos_max) {
				cos_max = _cos;
				cross_edge_dirs[i] = dirv;
			}
		}
	}
	// 理论上使用横坐标之差即可
	double dis1 = (point - v1).norm();
	double dis2 = (point - v2).norm();
	Vec cross_dirp = cross_edge_dirs[0] * dis2 + cross_edge_dirs[1] * dis1;
	cross_dirp.normalize();
	return get_cross_dirs(cross_dirp);
}

std::vector<Vec> QuadMesher::get_cross_dirs(const Vec& vec) {
	std::vector<Vec> res(4);
	Vec dirz(0, 0, 1);
	res[0] = vec;
	res[1] = dirz.cross(vec);
	res[1].normalize();
	res[2] = res[0] * -1;
	res[3] = res[1] * -1;
	return res;
}

void QuadMesher::calculate_cross_field() {

	/* VH is equal to ID */

	std::vector<Vec3> points;
	std::vector<ID3> triangles;
	std::vector<ID2> feature_edges;
	std::vector<ID3> singularities;// 奇异点的 id
	std::vector<std::array<double, 9>> global_triangle_dir;

	points.reserve(trimesh.vsize());
	triangles.reserve(trimesh.fsize());

	for (int i = 0; i < trimesh.vsize(); ++i) {
		VH vh(i);
		auto& v = trimesh.vertices(vh);
		points.push_back(Vec3{ v.x(), v.y(), v.z() });
	}

	for (int i = 0; i < trimesh.fsize(); ++i) {
		FH fh(i);
		auto& f = trimesh.faces(fh);
		const auto& vhs = f.getVertexHandle();
		triangles.push_back(ID3{ (ID)vhs[0], (ID)vhs[1], (ID)vhs[2] });
	}

	for (auto& eh : feature_ehs) {
		auto& e = trimesh.edges(eh);
		feature_edges.push_back(ID2{ (ID)e.vh1(), (ID)e.vh2() });
	}

	/*std::cout << "Initialize mesh success. #V = " << points.size() <<
		", #E = " << feature_edges.size() << ", #F = " << triangles.size() << std::endl;*/

	MeshMath::calc_cross_field(points, triangles, feature_edges, singularities, global_triangle_dir);

	int sv_cnt = 0, se_cnt = 0, sf_cnt = 0;
	SingularPoint point;
	for (auto& singular_id : singularities) {// 所有的奇异点
		if (singular_id[0] != NO_ID) {
			auto& v = trimesh.vertices(VH(singular_id[0]));
			point.pos = v;
			point.mode = 0;
			point.handle = singular_id[0];
			vh_to_spid[VH(singular_id[0])] = singular_points.size();
			sv_cnt++;
		} else if (singular_id[1] != NO_ID) {
			int fidx = singular_id[1] / 3;
			int eidx = singular_id[1] % 3;
			FH fh(fidx);
			EH eh = trimesh.faces(fh).eh(eidx);
			auto& edge = trimesh.edges(eh);
			auto& v1 = trimesh.vertices(edge.vh1());
			auto& v2 = trimesh.vertices(edge.vh2());
			point.pos = (v1 + v2) * 0.5;
			point.mode = 1;
			point.handle = eh;
			eh_to_spid[eh] = singular_points.size();
			se_cnt++;
		} else if (singular_id[2] != NO_ID) {
			FH fh(singular_id[2]);
			const auto& vhs = trimesh.faces(fh).getVertexHandle();
			auto& v0 = trimesh.vertices(vhs[0]);
			auto& v1 = trimesh.vertices(vhs[1]);
			auto& v2 = trimesh.vertices(vhs[2]);
			point.pos = (v0 + v1 + v2) / 3;
			point.mode = 2;
			point.handle = fh;
			fh_to_spid[fh] = singular_points.size();
			sf_cnt++;
		}
		singular_points.push_back(point);
	}

	for (int i = 0; i < trimesh.fsize(); ++i) {
		crossfield.push_back({
			Vec(global_triangle_dir[i][0], global_triangle_dir[i][1], global_triangle_dir[i][2]),
			Vec(global_triangle_dir[i][3], global_triangle_dir[i][4], global_triangle_dir[i][5]),
			Vec(global_triangle_dir[i][6], global_triangle_dir[i][7], global_triangle_dir[i][8])
			});
	}

	std::cout << "QuadMesher: Calculate cross field success. Total singular points = " << (sv_cnt + se_cnt + sf_cnt) << "."<< std::endl;
	std::cout << "\tSingular vertices = " << sv_cnt << ".\n";
	std::cout << "\tSingular edges = " << se_cnt << ".\n";
	std::cout << "\tSingular faces = " << sf_cnt << ".\n";

}

bool QuadMesher::get_intersection_line2(const Vex& src, const Vec& dir, EH eh, Vex& intersection) {

	double k0, b0, k1, b1;// 直线方程 y = kx + b
	bool line0_is_vertical = false, line1_is_vertical = false;
	double line0_x, line1_x;

	auto& edge = trimesh.edges(eh);
	auto& v1 = trimesh.vertices(edge.vh1());
	auto& v2 = trimesh.vertices(edge.vh2());
	if (v2.x() != v1.x()) {
		k0 = (v2.y() - v1.y()) / (v2.x() - v1.x());
		b0 = v2.y() - k0 * v2.x();
	} else {
		line0_is_vertical = true;
		line0_x = v2.x();
	}
	if (dir.x() != 0.f) {
		k1 = dir.y() / dir.x();
		b1 = src.y() - k1 * src.x();
	} else {
		line1_is_vertical = true;
		line1_x = src.x();
	}

	if (!line0_is_vertical && !line1_is_vertical) {

		double diffk = k0 - k1;
		if (diffk == 0) return false;// 平行
		double x = (b1 - b0) / diffk;
		double y = k0 * x + b0;
		if (x < std::min(v1.x(), v2.x()) || x > std::max(v1.x(), v2.x()) || 
			y < std::min(v1.y(), v2.y()) || y > std::max(v1.y(), v2.y())) {
			return false;
		}
		intersection = Vec(x, y, 0);

	} else if (line0_is_vertical && line1_is_vertical) {

		if (line0_x == line1_x) {
			double dis1 = (v1 - src).norm();
			double dis2 = (v2 - src).norm();
			if (dis1 < dis2) intersection = v2;
			else intersection = v1;
		} else return false;

	} else if (line0_is_vertical) {

		double y = line0_x * k1 + b1;
		if (y < std::min(v1.y(), v2.y()) || y > std::max(v1.y(), v2.y())) return false;
		intersection = Vec(line0_x, y, 0);

	} else if (line1_is_vertical) {

		double y = line1_x * k0 + b0;
		if (y < std::min(v1.y(), v2.y()) || y > std::max(v1.y(), v2.y())) return false;
		intersection = Vec(line1_x, y, 0);

	}

	/*std::cout << "Intersection! source = (" << src.x() << ", " << src.y() << "), dir = (" << dir.x() << ", " << dir.y() << "), v1 = ("
		<< v1.x() << ", " << v1.y() << "), v2 = (" << v2.x() << ", " << v2.y() << ")" << std::endl;
	std::cout << "\tIntersection = (" << intersection.x() << ", " << intersection.y() << ")\n";
	std::cout << "\tLine0: " << intersection.y() - (intersection.x() * k0 + b0) << ", Line1: " << intersection.y() - (intersection.x() * k1 + b1) << "\n";*/
	return true;

}

Vec QuadMesher::get_most_similar_dir(const Vec& ref_dir, const std::vector<Vec>& dirs) {
	double cosine_max = -2;
	Vec res;
	for (auto& dir : dirs) {
		double cosine = ref_dir.dot(dir);
		if (cosine > cosine_max) {
			cosine_max = cosine;
			res = dir;
		}
	}
	return res;
}

std::unordered_set<EH> QuadMesher::get_smallest_bounding_ehs(const EH& eh) {
	std::unordered_set<EH> res;
	auto& edge = trimesh.edges(eh);
	VH vh1 = edge.vh1(), vh2 = edge.vh2();
	std::unordered_set<FH> adjfhs = trimesh.NeighborFh(vh1);
	for (auto& adjfh : trimesh.NeighborFh(vh2)) adjfhs.insert(adjfh);
	for (auto& fh : adjfhs) {
		const auto& ehs = trimesh.faces(fh).getEdgeHandle();
		for (auto& eh : ehs) {
			auto& e = trimesh.edges(eh);
			if (e.vh1() == vh1 || e.vh2() == vh2) continue;
			res.insert(eh);
		}
	}
	return res;
}

void QuadMesher::get_streamline_itrp(const Streamline& sl0, const Streamline& sl1, Streamline& sl_new) {
	sl_new.beg_sid = sl0.beg_sid;
	sl_new.end_on_boundary = sl0.end_on_boundary && sl1.end_on_boundary;
	sl_new.end_sid = sl0.end_sid;
	int num = sl0.points.size() + sl1.points.size();
	if (num == 0) return;
	auto& points = sl_new.points;
	points.push_back((sl0.points[0] + sl1.points[0]) * 0.5f);
	double step = 1.f / num;
	int pid0 = 1, pid1 = 1;
	for (double r = step; r <= 1.f; r += step) {
		while (pid0 < sl0.ratio.size() && sl0.ratio[pid0] < r) pid0++;// 找到大于等于当前 r 的点
		while (pid1 < sl1.ratio.size() && sl1.ratio[pid1] < r) pid1++;
		if (pid0 == sl0.ratio.size() || pid1 == sl1.ratio.size()) break;
		double rdiff_00 = r - sl0.ratio[pid0 - 1];
		double rdiff_01 = sl0.ratio[pid0] - r;
		if (rdiff_00 < 0 || rdiff_01 < 0) {
			std::cerr << "QuadMesher: Wrong ratio!" << std::endl;
			break;
		}
		if (rdiff_00 + rdiff_01 == 0) {
			std::cerr << "QuadMesher: Devide zero!" << std::endl;
			break;
		} 
		Vex pos0 = (sl0.points[pid0 - 1] * rdiff_01 + sl0.points[pid0] * rdiff_00) / (rdiff_00 + rdiff_01);
		double rdiff_10 = r - sl1.ratio[pid1 - 1];
		double rdiff_11 = sl1.ratio[pid1] - r;
		if (rdiff_10 < 0 || rdiff_11 < 0) {
			std::cerr << "QuadMesher: Wrong num!" << std::endl;
			break;
		}
		if (rdiff_10 + rdiff_11 == 0) {
			std::cerr << "QuadMesher: Devide zero!" << std::endl;
			break;
		}
		Vex pos1 = (sl1.points[pid1 - 1] * rdiff_11 + sl1.points[pid1] * rdiff_10) / (rdiff_10 + rdiff_11);
		points.push_back((pos0 + pos1) * 0.5f);
	}
	if (sl_new.end_on_boundary) {// 将端点延伸至边界
		int n = sl_new.points.size();
		Vec dir_end = sl_new.points[n-1] - sl_new.points[n-2];
		dir_end.normalize();
		Vex intersection;
		EH eh0(sl0.end_edgehandle), eh1(sl1.end_edgehandle);
		std::set<EH> bdehs;
		get_among_ehs(eh0, eh1, bdehs);
		bool intersection_flag = false;
		for (auto& bdeh : bdehs) {
			if (get_intersection_line2(sl_new.points.back(), dir_end, bdeh, intersection)) {
				intersection_flag = true;
				break;
			}
		}
		if (!intersection_flag) {
			std::cerr << "QuadMesher: Boundary ehs not intersection with the ray.\n";
			return;
		}
		sl_new.points.push_back(intersection);
	} else {// 确保端点在奇异点上
		Vex endpoint = (sl0.points.back() + sl1.points.back()) * 0.5;
		if ((sl_new.points.back() - endpoint).norm() > 1E-6F) {
			sl_new.points.push_back(endpoint);
		}
	}
	sl_new.init_length_ratio();

}

void QuadMesher::get_among_ehs(const EH& beg, const EH& end, std::set<EH>& res) {
	if (beg == end) {
		res = { beg };
		return;
	}
	if (!trimesh.isOnBoundary(beg) || !trimesh.isOnBoundary(end)) {
		std::cerr << "QuadMesher: There exist eh not on the boundary.\n";
		return;
	}
	std::set<EH> ehs0, ehs1;
	ehs0.insert(beg), ehs1.insert(beg);
	EH eh0(-1), eh1(-1);
	for (auto& eh : trimesh.NeighborEh(beg)) {
		if (trimesh.isOnBoundary(eh)) {
			if (eh0 == -1) eh0 = eh;
			else if (eh1 == -1) eh1 = eh;
			else {
				std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
				return;
			}
		}
	}
	while (1) {
		ehs0.insert(eh0);
		if (eh0 == end) {
			res = ehs0;
			break;
		}

		ehs1.insert(eh1);
		if (eh1 == end) {
			res = ehs1;
			break;
		}

		int degree = 0;
		EH eh_tmp(-1);

		for (auto& eh : trimesh.NeighborEh(eh0)) {
			if (trimesh.isOnBoundary(eh)) {
				if (!ehs0.count(eh)) {
					eh_tmp = eh;
					degree++;
				}
			}
		}
		if (degree != 1 || eh_tmp == -1) {
			std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
			return;
		}
		eh0 = eh_tmp;

		degree = 0;
		eh_tmp = EH(-1);
		for (auto& eh : trimesh.NeighborEh(eh1)) {
			if (trimesh.isOnBoundary(eh)) {
				if (!ehs1.count(eh)) {
					eh_tmp = eh;
					degree++;
				}
			}
		}
		if (degree != 1 || eh_tmp == -1) {
			std::cerr << "QuadMesher: It is not a 2D-manifold mesh.\n";
			return;
		}
		eh1 = eh_tmp;

	}
}


void QuadMesher::generate_initial_streamlines_v(const SingularPoint& point) {

	VH vh(point.handle);

	std::vector<FH> fhs;
	std::vector<EH> ehs;// 最小多边形的边
	for (auto& adjfh : trimesh.NeighborFh(vh)) {
		const auto& f_ehs = (trimesh.faces(adjfh)).getEdgeHandle();
		for (auto& eh : f_ehs) {
			auto& edge = trimesh.edges(eh);
			if (edge.vh1() != vh && edge.vh2() != vh) {
				ehs.emplace_back(eh);
				break;
			}
		}
		fhs.emplace_back(adjfh);
	}
	assert(fhs.size() == ehs.size());
	int degree = 0;
	auto& v = trimesh.vertices(vh);
	std::vector<MorchNode> nodes;
	for (int i = 0; i < fhs.size(); ++i) {
		double weight_vh1 = 0;
		Vec march_dir;
		if (is_on_streamline(v, ehs[i], fhs[i], march_dir, weight_vh1)) {
			degree++;
			auto& edge = trimesh.edges(ehs[i]);
			MorchNode tmp_node;
			tmp_node.dir = march_dir;
			tmp_node.pre_fh = fhs[i];
			tmp_node.eh = ehs[i];
			tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1 + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1);
			nodes.push_back(tmp_node);
		}
	}

	for (auto& node : nodes) {
		Streamline sline;
		sline.beg_sid = vh_to_spid[vh];
		sline.points.push_back(v);
		std::set<EH> visited;
		morch_on_streamline(node, sline, visited);
		std::cout << "Streamline_V #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		if (!sline.end_on_fault) streamlines.push_back(sline);
	}

}

void QuadMesher::generate_initial_streamlines_e(const SingularPoint& point) {

	EH eh(point.handle);
	auto& edge = trimesh.edges(eh);
	std::unordered_set<FH> fhs = trimesh.NeighborFh(eh);
	
	Vex pos = point.pos;

	int degree = 0;
	std::vector<MorchNode> nodes;
	for (auto& fh : fhs) {
		for (auto& adjeh : trimesh.faces(fh).getEdgeHandle()) {
			if (adjeh == eh) continue;
			double weight_vh1 = 0;
			Vec march_dir;
			if (is_on_streamline(pos, adjeh, fh, march_dir, weight_vh1)) {
				degree++;
				auto& edge = trimesh.edges(adjeh);
				MorchNode tmp_node;
				tmp_node.dir = march_dir;
				tmp_node.pre_fh = fh;
				tmp_node.eh = adjeh;
				tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1 + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1);
				nodes.push_back(tmp_node);
			}

		}
	}

	for (auto& node : nodes) {
		Streamline sline;
		sline.beg_sid = eh_to_spid[eh];
		sline.points.push_back(point.pos);
		std::set<EH> visited;
		morch_on_streamline(node, sline, visited);
		std::cout << "Streamline_E #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		if (!sline.end_on_fault) streamlines.push_back(sline);
	}
	
}

void QuadMesher::generate_initial_streamlines_f(const SingularPoint& point) {

	FH fh(point.handle);
	Vex pos = point.pos;
	int degree = 0;
	std::vector<MorchNode> nodes;

	for (auto& eh : trimesh.faces(fh).getEdgeHandle()) {

		double weight_vh1 = 0;
		Vec march_dir;
		if (is_on_streamline(pos, eh, fh, march_dir, weight_vh1)) {
			degree++;
			auto& edge = trimesh.edges(eh);
			MorchNode tmp_node;
			tmp_node.dir = march_dir;
			tmp_node.pre_fh = fh;
			tmp_node.eh = eh;
			tmp_node.pos = (trimesh.vertices(edge.vh1())) * weight_vh1 + (trimesh.vertices(edge.vh2())) * (1 - weight_vh1);
			nodes.push_back(tmp_node);
		}

	}

	for (auto& node : nodes) {
		Streamline sline;
		sline.beg_sid = fh_to_spid[fh];
		sline.points.push_back(point.pos);
		std::set<EH> visited;
		morch_on_streamline(node, sline, visited);
		std::cout << "Streamline_F #" << streamlines.size() << ": #P = " << sline.points.size() << std::endl;
		if (!sline.end_on_fault) streamlines.push_back(sline);
	}

}

bool QuadMesher::get_instersection_polygon(const Vex& src, const Vec& dir, const EH& src_eh, const FH& src_fh, 
	EH& intersect_eh, FH& intersect_fh, Vex& intersection) {

	intersect_eh = EH(-1);
	intersect_fh = FH(-1);

	FH adjfh(-1);
	for (auto& fh : trimesh.NeighborFh(src_eh)) {
		if (fh != adjfh) {
			adjfh = fh;
			break;
		}
	}

	int intersection_cnt = 0;

	for (auto& feh : trimesh.faces(adjfh).getEdgeHandle()) {
		if (feh == src_eh) continue;
		Vex _intersection;
		if (get_intersection_line2(src, dir, feh, _intersection)) {
			// 判断是否是反向交点
			Vec dir_si = (_intersection - src).normalize();
			if (dir_si.dot(dir) <= 0.f) continue;

			intersection_cnt++;
			intersect_eh = feh;
			intersect_fh = adjfh;
			intersection = _intersection;
			break;// 只取一个顶点
		}
	}
	if (intersection_cnt == 0) {// 需要遍历这条边的最小多边形的边来求交点
		std::unordered_set<EH> bounding_ehs = get_smallest_bounding_ehs(src_eh);
		auto& edge = trimesh.edges(src_eh);
		for (auto& feh : bounding_ehs) {
			Vex _intersection;
			if (get_intersection_line2(src, dir, feh, _intersection)) {
				// 判断是否是反向交点
				Vec dir_si = (_intersection - src).normalize();
				if (dir_si.dot(dir) <= 0.f) continue;

				intersection_cnt++;
				intersect_eh = feh;
				intersection = _intersection;
				for (auto& fh : trimesh.NeighborFh(feh)) {
					const auto& vhs = trimesh.faces(fh).getVertexHandle();
					if (std::find(vhs.begin(), vhs.end(), edge.vh1()) != vhs.end()
						|| std::find(vhs.begin(), vhs.end(), edge.vh2()) != vhs.end()) {
						intersect_fh = fh;
						break;
					}
				}
				break;// 只取一个顶点
			}
		}
	}

	return intersect_eh != -1 && intersect_fh != -1;

}