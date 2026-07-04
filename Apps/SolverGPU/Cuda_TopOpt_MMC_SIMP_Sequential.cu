#include "Cuda_TopOpt_MMC_SIMP_Sequential.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <vector>

#include "Cuda_TopOpt_MultiMaterial_SIMP.h"
#include "Cuda_TopOpt_Pipeline.h"
#include "Poster/source/Timer.h"
#include "Poster/source/nanoflann.h"
#include "Util/GP_Log.h"

namespace {
using namespace nanoflann;

struct Point {
	double x, y, z;
};

template <typename T>
struct PointCloud {
	using coord_t = T;
	std::vector<Point> pts;

	inline size_t kdtree_get_point_count() const { return pts.size(); }
	inline T kdtree_get_pt(const size_t idx, const size_t dim) const {
		if (dim == 0) return pts[idx].x;
		if (dim == 1) return pts[idx].y;
		return pts[idx].z;
	}
	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const { return false; }
};

using KDTree = KDTreeSingleIndexAdaptor<
	L2_Simple_Adaptor<double, PointCloud<double>>,
	PointCloud<double>,
	3>;

void GetNearestIds(std::vector<Point>& p, std::vector<Point>& q, std::vector<int>& ans) {
	ans.resize(p.size());
	PointCloud<double> A, B;
	A.pts = p;
	B.pts = q;
	KDTree index(3, B, KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();
	for (size_t a_index = 0; a_index < A.pts.size(); ++a_index) {
		size_t ret_index = 0;
		double out_dist_sqr = 0.0;
		nanoflann::KNNResultSet<double> resultSet(1);
		resultSet.init(&ret_index, &out_dist_sqr);
		index.findNeighbors(resultSet, &A.pts[a_index].x);
		ans[a_index] = static_cast<int>(ret_index);
	}
}
}

void Cuda_TopOpt_MMC_SIMP_Sequential(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double mmc_target_rate,
	double mmc_max_iter,
	double material_ratio_f,
	double simp_max_iter,
	double E1,
	double E2,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model,
	int mmcType) {
	Cuda_TopOpt_MMC_Node(allVertex, cell2V, input_model_file, bdy_file, protect_cell_ids, neighbor_info,
		mmc_target_rate, mmc_max_iter, E1, Poissons, density, is_need_update_model, mmcType);
	std::vector<double> mmc_density = density;

	// MATLAB 参考代码中为 Ei = find(denk > alpha^2)。
	// 当前 3D MMC 输出的 mmc_density 已经是单元节点 Heaviside 值的平均值，
	// 没有再平方，因此这里用 alpha 作为 SIMP 第二阶段有效设计域阈值。
	const double alpha = 1e-3;
	std::vector<int> active_indices;
	for (int i = 0; i < static_cast<int>(mmc_density.size()); ++i) {
		if (mmc_density[i] > alpha) {
			active_indices.push_back(i);
		}
	}
	if (active_indices.empty()) {
		is_need_update_model = true;
		return;
	}

	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;
	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml_seq");

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;
	{
		int control_size = 0;
		std::ifstream in(bdy_file);
		in >> control_size;
		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			bdy_points.emplace_back(x, y, z);
			bdy_value.emplace_back(vx, vy, vz);
			bdy_type.push_back(cur_type);
		}
	}

	gpmesh::GPLog::init(spdlog::level::trace);
	CUDA_ERROR(cudaDeviceReset());

	gpmesh::GPSplineSolidV0 spline;
	spline.init_patch(patch);

	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;

	const SparseSystem& temp_sys = spline.m_sparse_system;
	{
		std::vector<Point> A, B;
		std::vector<int> global_ids;
		for (auto& p : bdy_points) {
			A.push_back(Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < static_cast<int>(patch.size()); ++p_id) {
			for (int local_id = 0; local_id < static_cast<int>(patch[p_id].coefs().size()); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.push_back(Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}
		GetNearestIds(A, B, global_ids);
		auto& global_local_map = temp_sys.get_l2g_map();
		for (int i = 0; i < static_cast<int>(bdy_points.size()); ++i) {
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id];
			if (global_id == -1) {
				continue;
			}
			if (bdy_type[i] == 2) {
				neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
			}
			else if (bdy_type[i] == 1) {
				dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
			}
		}
	}

	spline.init_dirichlet_points(dirichlet_points);

	std::vector<Eigen::MatrixXd> unit_local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
	spline.Get_Local_Matrix(unit_local_matrix, 1.0, Poissons);
	spline.m_sparse_system.rhs().setZero();
	spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
	spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);

	MultiMaterialSIMPParams params;
	params.material_ratio_f = material_ratio_f;
	params.E1 = E1;
	params.E2 = E2;
	// 对应论文 / MATLAB 初始化：Svariable(1:nn, 1) = volfrac。
	// 第二阶段材料选择变量 rho_e 从 MMC 体积分数开始。
	params.initial_rho = mmc_target_rate;
	params.max_iter = static_cast<int>(simp_max_iter);

	RunMultiMaterialSIMPStage(spline, neumann_points, dirichlet_points, neighbor_info, allVertex, cell2V,
		mmc_density, active_indices, unit_local_matrix, params, "mmc_simp_result", density, is_need_update_model);
}
