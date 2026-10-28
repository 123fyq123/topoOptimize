#include "Cuda_TopOpt_MMC_SIMP_Sequential.h"

#include <algorithm>
#include <chrono>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <vector>

#include "Cuda_TopOpt_MultiMaterial_SIMP.h"
#include "Cuda_TopOpt_Pipeline.h"
#include "OptimizationPaperRecorder.h"
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
	const std::string output_dir = "mmc_simp_result";
	const optimization_paper_record::RunMetadata run_metadata =
		optimization_paper_record::CreateRunMetadata(
			output_dir,
			"sequential",
			input_model_file,
			bdy_file);
	const std::string record_dir = run_metadata.record_directory;
	const std::string paper_iteration_path = record_dir + "/paper_iterations.tsv";
	const auto sequential_run_start = std::chrono::steady_clock::now();
	std::filesystem::create_directories(output_dir);
	// 顺序版 MMC 阶段复用纯 MMC，Heaviside 平滑宽度全程固定为 0.2。
	auto config_entries = optimization_paper_record::MetadataEntries(run_metadata);
	config_entries.insert(
		config_entries.end(),
		{
			{ "input_model_file", input_model_file },
			{ "mmc_type", optimization_paper_record::ToString(mmcType) },
			{ "mmc_target_rate", optimization_paper_record::ToString(mmc_target_rate) },
			{ "mmc_max_iterations", optimization_paper_record::ToString(mmc_max_iter) },
			{ "material_max_iterations", optimization_paper_record::ToString(simp_max_iter) },
			{ "mmc_youngs_modulus", optimization_paper_record::ToString(E1) },
			{ "material_ratio_f", optimization_paper_record::ToString(material_ratio_f) },
			{ "material_E1", optimization_paper_record::ToString(E1) },
			{ "material_E2", optimization_paper_record::ToString(E2) },
			{ "poissons_ratio", optimization_paper_record::ToString(Poissons) },
			{ "heaviside_epsilon_mode", "fixed" },
			{ "heaviside_epsilon_initial", "0.2" },
			{ "heaviside_epsilon_middle", "0.2" },
			{ "heaviside_epsilon_final", "0.2" },
			{ "material_penalty", "3" },
			{ "material_active_threshold", "0.01" },
			{ "material_rho_min", "0.001" },
			{ "material_move_limit", "0.2" },
			{ "mmc_fe_solves_per_iteration", "1" },
			{ "material_fe_solves_per_iteration", "1" }
		});
	// MMC 阶段沿用 E1 作为基准模量；后续材料阶段直接使用界面传入的 E1/E2。
	if (!optimization_paper_record::WriteConfig(
		record_dir + "/paper_run_config.tsv",
		config_entries)) {
		throw std::runtime_error("Cannot create the sequential paper run configuration.");
	}
	optimization_paper_record::WriteConfig(
		output_dir + "/latest_run.tsv",
		{
			{ "run_id", run_metadata.run_id },
			{ "record_directory", record_dir }
		});
	int mmc_completed_iterations = 0;
	Cuda_TopOpt_MMC_Node(allVertex, cell2V, input_model_file, bdy_file, protect_cell_ids, neighbor_info,
		mmc_target_rate, mmc_max_iter, E1, Poissons, density, is_need_update_model, mmcType,
		paper_iteration_path, &mmc_completed_iterations);
	const double mmc_elapsed_seconds =
		optimization_paper_record::ElapsedSeconds(sequential_run_start);
	std::vector<double> mmc_density = density;

	// MMC 阶段结束后，只截取实体占有率足够高的单元作为 SIMP 第二阶段设计域。
	// 当前 3D MMC 输出的 mmc_density 是八节点 Heaviside 平均值，因此这里直接按密度阈值筛选。
	// 与同步版保持一致：MMC 占有率超过 0.01 的单元均进入材料优化域和最终 VTK。
	const double mmc_active_density_threshold = 0.01;
	std::vector<int> active_indices;
	for (int i = 0; i < static_cast<int>(mmc_density.size()); ++i) {
		if (mmc_density[i] > mmc_active_density_threshold) {
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
	// rho=1 对应 E1，rho 接近 0 对应 E2；两者均来自顺序优化界面。
	params.E1 = E1;
	params.E2 = E2;
	// MATLAB convention: rho represents E1, so its feasible fraction is V1/(V1+V2).
	params.initial_rho = 1.0 / (1.0 + material_ratio_f);
	params.max_iter = static_cast<int>(simp_max_iter);

	RunMultiMaterialSIMPStage(spline, neumann_points, dirichlet_points, neighbor_info, allVertex, cell2V,
		mmc_density, active_indices, unit_local_matrix, params, output_dir, record_dir, density, is_need_update_model,
		mmc_completed_iterations, mmc_completed_iterations, mmc_elapsed_seconds);
	// 固定目录保留最新结果供界面读取，独立运行目录保存与本次参数严格对应的结果。
	const std::filesystem::path latest_vtk =
		std::filesystem::path(output_dir) / "mmc_simp_result.vtk";
	if (std::filesystem::exists(latest_vtk)) {
		std::filesystem::copy_file(
			latest_vtk,
			std::filesystem::path(record_dir) / latest_vtk.filename(),
			std::filesystem::copy_options::overwrite_existing);
	}
}
