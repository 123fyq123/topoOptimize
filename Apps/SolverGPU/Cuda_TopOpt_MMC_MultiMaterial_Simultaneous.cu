#include "Cuda_TopOpt_MMC_MultiMaterial_Simultaneous.h"

#include <algorithm>
#include <array>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <numeric>
#include <set>
#include <stdexcept>
#include <utility>

#include <Eigen/Dense>

#include "GP_Spline_v0.h"
#include "MMASolver.h"
#include "MMC_Trace_Diagnostics.h"
#include "Poster/source/Timer.h"
#include "Poster/source/nanoflann.h"
#include "Util/GP_Log.h"
#include "Util/GP_Macros.h"

namespace mmc_multimaterial_simultaneous {

// ============================== Review map ==============================
// 这份文件只实现新增的同步流程，不进入现有 MMC/SIMP/顺序优化主循环。
// 建议按以下顺序 review：
//   1. calc_Phi_Node()：直接沿用独立 MMC 的单组件 TDF 与 9 个几何导数；
//   2. ComputeGlobalTDF()/EvaluateGeometry()：KS -> Heaviside -> m_e；
//   3. ComputeMmcDerivativeFromElementWeights()：复用的 MMC 链式法则；
//   4. EvaluateSimultaneousState()：a_e=m_e E_e，一次共享结构分析后同时给出两块梯度；
//   5. Cuda_TopOpt_MMC_MultiMaterial_Simultaneous()：MMC-MMA + 材料 OC；
// 全量有限差分帮助函数暂时保留，但已从正常同步优化入口断开。
//
// 核心符号：
//   d_ir : 第 i 个 MMC 组件的第 r 个几何变量
//   phi_i: 第 i 个组件的局部拓扑描述函数
//   Phi  : 所有组件经 KS 聚合后的全局 TDF
//   m_e  : 单元 8 个顶点 Heaviside 值的平均，即 MMC 占有率
//   eta_e: 单元双材料变量
//   E_e  : eta_e 对应的插值材料模量
//   a_e  : m_e E_e，有限元装配实际使用的单元缩放系数
//   q_e  : U_e^T K_e^0 U_e，单位模量单元能量
// =======================================================================

constexpr int kVariablesPerComponent = 9;
constexpr int kTdfPower = 6;
constexpr double kTiny = 1e-12;
constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kFiniteDifferencePassRelativeTolerance = 1e-3;
constexpr double kFiniteDifferenceWarnRelativeTolerance = 1e-2;
constexpr double kFiniteDifferenceNearZeroTolerance = 1e-10;
constexpr double kFiniteDifferencePassAbsoluteTolerance = 1e-8;

struct MMC_Component {
	// 变量顺序必须与独立 3D MMC 保持一致；MMA 展平和灵敏度索引都依赖该顺序。
	double center_x = 0.0;
	double center_y = 0.0;
	double center_z = 0.0;
	double length_x = 0.0;
	double length_y = 0.0;
	double length_z = 0.0;
	double rotation_alpha = 0.0;
	double rotation_beta = 0.0;
	double rotation_gamma = 0.0;

	MMC_Component() = default;

	MMC_Component(
		double x,
		double y,
		double z,
		double length1,
		double length2,
		double length3,
		double rotation_alpha,
		double rotation_beta,
		double rotation_gamma)
		: center_x(x)
		, center_y(y)
		, center_z(z)
		, length_x(length1)
		, length_y(length2)
		, length_z(length3)
		, rotation_alpha(rotation_alpha)
		, rotation_beta(rotation_beta)
		, rotation_gamma(rotation_gamma) {
	}
};

// 新增的同步辅助函数原先使用 MMCComponent；保留别名避免无意义地扩大改动，
// MMC 公共计算段则直接沿用独立流程的 MMC_Component 名称。
using MMCComponent = MMC_Component;

struct Point {
	double x = 0.0;
	double y = 0.0;
	double z = 0.0;
};

template <typename T>
struct PointCloud {
	using coord_t = T;
	std::vector<Point> pts;

	inline size_t kdtree_get_point_count() const { return pts.size(); }

	inline T kdtree_get_pt(size_t index, size_t dimension) const {
		if (dimension == 0) return static_cast<T>(pts[index].x);
		if (dimension == 1) return static_cast<T>(pts[index].y);
		return static_cast<T>(pts[index].z);
	}

	template <class BBOX>
	bool kdtree_get_bbox(BBOX&) const {
		return false;
	}
};

using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<double, PointCloud<double>>,
	PointCloud<double>,
	3>;

struct GeometryEvaluation {
	// 以下名称与 Cuda_TopOpt_MMC_Node 主循环保持一致，便于逐行对照。
	std::vector<std::vector<double>> allPhis;
	std::vector<std::vector<std::vector<double>>> allPhiDrv;
	std::vector<double> globalPhi;
	std::vector<std::vector<double>> phimaxDphi;
	std::vector<double> vertexH;
	std::vector<double> delta_H;
	// 这里只表示 MMC 单元占有率 m_e；材料场和等效刚度另存，不混入 density。
	std::vector<double> density;
};

struct SimultaneousEvaluation {
	// 一次已接受联合状态的全部有限元量、目标/约束和解析梯度。
	// 两组目标梯度必须共享这里同一个 displacement 和 element_energy。
	GeometryEvaluation geometry;
	Eigen::MatrixXd displacement;
	std::vector<double> material_modulus;
	std::vector<double> material_modulus_derivative;
	std::vector<double> projected_occupancy; // 兼容旧输出字段；当前数值路径中 chi_e == m_e
	std::vector<double> projected_occupancy_derivative;
	double projection_beta = 0.0;
	std::vector<double> effective_stiffness; // a_e = m_e E_e
	std::vector<double> element_energy;      // q_e = U_e^T K_e^0 U_e
	std::vector<double> mmc_objective_gradient;
	std::vector<double> material_objective_gradient;
	std::vector<double> geometry_constraint_gradient;
	std::vector<double> ratio_constraint_mmc_gradient;
	std::vector<double> ratio_constraint_material_gradient;
	std::vector<double> mmc_contribution_max_abs;
	std::vector<double> mmc_contribution_sum_abs;
	double compliance = 0.0;
	double compliance_from_load = 0.0;
	double compliance_relative_mismatch = 0.0;
	double mmc_volume = 0.0;
	double V1 = 0.0;
	double V2 = 0.0;
	double ratio_residual = 0.0;
	double geometry_constraint = 0.0;
	double positive_ratio_constraint = 0.0;
	double negative_ratio_constraint = 0.0;
	double uniform_material_gradient_relative_error = 0.0;
};

bool IsFiniteVector(const std::vector<double>& values) {
	return std::all_of(values.begin(), values.end(), [](double value) {
		return std::isfinite(value);
	});
}

double L2Norm(const std::vector<double>& values) {
	double squared_norm = 0.0;
	for (double value : values) squared_norm += value * value;
	return std::sqrt(squared_norm);
}

void SetFailure(MMCMultiMaterialSimultaneousResult& result, const std::string& message) {
	result.success = false;
	result.message = message;
	std::cerr << "[MMC-MultiMaterial-Simultaneous] " << message << std::endl;
}

bool ValidateInputs(
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::string& input_model_file,
	const std::string& boundary_file,
	const std::vector<int>& protected_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	const MMCMultiMaterialSimultaneousParams& params,
	const std::vector<double>& density,
	std::string& error) {
	const int vertex_count = static_cast<int>(all_vertex.size());
	const int cell_count = static_cast<int>(cell_to_vertex.size());
	if (vertex_count == 0 || cell_count == 0) {
		error = "The design mesh is empty.";
		return false;
	}
	if (!std::filesystem::exists(input_model_file) || !std::filesystem::exists(boundary_file)) {
		error = "The Bezier model file or boundary-condition file does not exist.";
		return false;
	}
	if (static_cast<int>(density.size()) != cell_count) {
		error = "The density count does not match the fixed-grid cell count.";
		return false;
	}
	if (!neighbor_info.empty() && static_cast<int>(neighbor_info.size()) != cell_count) {
		error = "The neighbor_info count does not match the fixed-grid cell count.";
		return false;
	}
	for (int vertex_id = 0; vertex_id < vertex_count; ++vertex_id) {
		if (all_vertex.find(vertex_id) == all_vertex.end()) {
			error = "MMC vertex ids must be contiguous and start at zero.";
			return false;
		}
	}
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		auto cell_it = cell_to_vertex.find(cell_id);
		if (cell_it == cell_to_vertex.end() || cell_it->second.size() != 8) {
			error = "Fixed-grid cell ids must be contiguous and each cell must have eight vertices.";
			return false;
		}
		for (const auto& vertex : cell_it->second) {
			if (vertex.first < 0 || vertex.first >= vertex_count) {
				error = "A cell contains an invalid vertex id.";
				return false;
			}
		}
	}
	for (int cell_id : protected_cell_ids) {
		if (cell_id < 0 || cell_id >= cell_count) {
			error = "A protected cell id is out of range.";
			return false;
		}
	}
	if (!(params.mmc_target_rate > 0.0 && params.mmc_target_rate <= 1.0)) {
		error = "The MMC target volume fraction must be in (0, 1].";
		return false;
	}
	if (params.max_iterations <= 0) {
		error = "The maximum iteration count must be positive.";
		return false;
	}
	if (!(params.E1 > 0.0) || !(params.E2 > 0.0) ||
		!std::isfinite(params.E1) || !std::isfinite(params.E2)) {
		error = "Both material moduli must be finite and positive.";
		return false;
	}
	if (!(params.E1 < params.E2)) {
		error =
			"The soft-material convention requires E1 < E2 so that larger eta means softer material.";
		return false;
	}
	if (!(params.material_penalty >= 1.0) || !std::isfinite(params.material_penalty)) {
		error = "The material penalty must be finite and at least one.";
		return false;
	}
	if (!(params.material_ratio_f >= 0.0) || !std::isfinite(params.material_ratio_f)) {
		error = "The material ratio f must be finite and nonnegative.";
		return false;
	}
	if (1.0 / (1.0 + params.material_ratio_f) < params.eta_min) {
		error =
			"The requested material ratio is infeasible for the configured eta_min.";
		return false;
	}
	if (!(params.eta_min > 0.0 && params.eta_min < 1.0)) {
		error = "eta_min must be in (0, 1).";
		return false;
	}
	if (!(params.projection_threshold > 0.0 &&
		params.projection_threshold < 1.0) ||
		!std::isfinite(params.projection_threshold)) {
		error = "The projection threshold must be finite and in (0, 1).";
		return false;
	}
	if (!(params.material_binary_active_threshold >= 0.0 &&
		params.material_binary_active_threshold < 1.0) ||
		!std::isfinite(params.material_binary_active_threshold)) {
		error = "The material binary threshold must be finite and in [0, 1).";
		return false;
	}
	if (!(params.material_active_threshold >= 0.0 &&
		params.material_active_threshold < 1.0) ||
		!std::isfinite(params.material_active_threshold)) {
		error = "The material OC active threshold must be finite and in [0, 1).";
		return false;
	}
	if (!(params.projection_beta_initial > 0.0) ||
		!(params.projection_beta_max >= params.projection_beta_initial) ||
		!std::isfinite(params.projection_beta_initial) ||
		!std::isfinite(params.projection_beta_max) ||
		params.projection_beta_update_interval <= 0) {
		error = "The projection beta continuation parameters are invalid.";
		return false;
	}
	if (!(params.Poissons > -1.0 && params.Poissons < 0.5)) {
		error = "Poisson's ratio must be in (-1, 0.5).";
		return false;
	}
	if (!(params.heaviside_alpha >= 0.0 && params.heaviside_alpha < 1.0) ||
		!std::isfinite(params.heaviside_alpha)) {
		error = "Heaviside alpha must be finite and in [0, 1).";
		return false;
	}
	if (!(params.heaviside_epsilon > 0.0) ||
		!(params.ks_lambda > 0.0) ||
		!std::isfinite(params.heaviside_epsilon) ||
		!std::isfinite(params.ks_lambda)) {
		error = "Heaviside epsilon and KS lambda must be positive.";
		return false;
	}
	if (params.finite_difference_mmc_sample_count < 0 ||
		params.finite_difference_material_sample_count < 0) {
		error = "Finite-difference sample counts must be nonnegative.";
		return false;
	}
	if (!(params.finite_difference_step > 0.0) ||
		!std::isfinite(params.finite_difference_step)) {
		error = "The finite-difference step must be finite and positive.";
		return false;
	}
	return true;
}

void GetNearestIds(
	const std::vector<Point>& source_points,
	const std::vector<Point>& target_points,
	std::vector<int>& nearest_ids) {
	if (source_points.empty() || target_points.empty()) {
		throw std::runtime_error("The boundary-point or control-point set is empty.");
	}
	PointCloud<double> cloud;
	cloud.pts = target_points;
	KDTree index(3, cloud, nanoflann::KDTreeSingleIndexAdaptorParams(10));
	index.buildIndex();

	nearest_ids.resize(source_points.size());
	for (size_t source_id = 0; source_id < source_points.size(); ++source_id) {
		size_t nearest_id = 0;
		double squared_distance = 0.0;
		nanoflann::KNNResultSet<double> result_set(1);
		result_set.init(&nearest_id, &squared_distance);
		index.findNeighbors(result_set, &source_points[source_id].x);
		nearest_ids[source_id] = static_cast<int>(nearest_id);
	}
}

double Heaviside(double phi, double alpha, double epsilon) {
	// 与独立 MMC 一致的正则化 Heaviside：过渡带外分别取 1 和 alpha，
	// 过渡带 [-epsilon, epsilon] 内使用三次多项式保持可导。
	if (phi > epsilon) return 1.0;
	if (phi < -epsilon) return alpha;
	return 3.0 * (1.0 - alpha) / 4.0
		* (phi / epsilon - std::pow(phi, 3.0) / (3.0 * std::pow(epsilon, 3.0)))
		+ (1.0 + alpha) / 2.0;
}

double HeavisideDerivative(double phi, double alpha, double epsilon) {
	if (std::abs(phi) > epsilon) return 0.0;
	return 3.0 * (1.0 - alpha) / (4.0 * epsilon)
		* (1.0 - phi * phi / (epsilon * epsilon));
}

double InterpolateYoungsModulus(
	double rho,
	const MMCMultiMaterialSimultaneousParams& params) {
	// 与顺序版同名函数保持 rho/rho_power_p 命名；同步流程中的 rho 即 eta。
	const double rho_power_p = std::pow(rho, params.material_penalty);
	return rho_power_p * params.E1 + (1.0 - rho_power_p) * params.E2;
}

double InterpolateYoungsModulusDerivative(
	double rho,
	const MMCMultiMaterialSimultaneousParams& params) {
	// dE_e/deta_e = p eta_e^(p-1) (E1-E2)。
	return params.material_penalty
		* std::pow(rho, params.material_penalty - 1.0)
		* (params.E1 - params.E2);
}

// 与顺序 MMC-SIMP 保持一致：材料 OC 更新前先做 Sigmund 灵敏度过滤。
std::vector<double> FilterSensitivities(
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	const std::vector<double>& rho_full,
	const std::vector<double>& dc_full,
	double rho_min) {
	std::vector<double> filtered(dc_full.size(), 0.0);
	for (int i = 0; i < static_cast<int>(dc_full.size()); ++i) {
		if (rho_full[i] <= 0.0) {
			continue;
		}
		double sum_weight = 0.0;
		double accum = 0.0;
		for (const auto& [j, weight] : neighbor_info[i]) {
			if (j < 0 || j >= static_cast<int>(dc_full.size())) {
				continue;
			}
			if (rho_full[j] <= 0.0) {
				continue;
			}
			const double clamped_weight = std::max(0.0, weight);
			sum_weight += clamped_weight;
			accum += clamped_weight * rho_full[j] * dc_full[j];
		}
		// Sigmund 灵敏度过滤：
		//   dc_hat_i = sum_j(w_ij * rho_j * dc_j) / (rho_i * sum_j(w_ij))
		const double denom = std::max(rho_full[i], rho_min) * std::max(sum_weight, 1e-12);
		filtered[i] = accum / denom;
	}
	return filtered;
}

double ProjectionBetaForIteration(
	const MMCMultiMaterialSimultaneousParams& params,
	int evaluation_iteration) {
	const int continuation_stage = std::max(0, evaluation_iteration - 1)
		/ params.projection_beta_update_interval;
	return std::min(
		params.projection_beta_max,
		params.projection_beta_initial * std::pow(2.0, continuation_stage));
}

double ProjectOccupancy(
	double occupancy,
	double beta,
	double threshold) {
	const double denominator =
		std::tanh(beta * threshold) + std::tanh(beta * (1.0 - threshold));
	return (std::tanh(beta * threshold)
		+ std::tanh(beta * (std::clamp(occupancy, 0.0, 1.0) - threshold)))
		/ denominator;
}

double ProjectOccupancyDerivative(
	double occupancy,
	double beta,
	double threshold) {
	const double denominator =
		std::tanh(beta * threshold) + std::tanh(beta * (1.0 - threshold));
	const double cosh_term = std::cosh(beta * (occupancy - threshold));
	return beta / (denominator * cosh_term * cosh_term);
}

// 材料块沿用 OC 二分更新。MMC 只作为当前轮的常量权重进入目标灵敏度
// 和材料资源统计，不在材料子更新中改变组件几何。
void OCUpdate(
	std::vector<double>& x,
	double material_ratio_f,
	const std::vector<double>& dc,
	const std::vector<double>& raw_occupancy,
	const std::vector<double>& material_weight,
	double active_threshold,
	double rho_min,
	double move_limit) {
	const int nn = static_cast<int>(x.size());
	if (static_cast<int>(dc.size()) != nn ||
		static_cast<int>(raw_occupancy.size()) != nn ||
		static_cast<int>(material_weight.size()) != nn) {
		throw std::runtime_error(
			"The material OC input sizes do not match the design vector.");
	}
	std::vector<double> xnew(x);
	std::vector<unsigned char> active(nn, 0);
	double inactive_soft_resource = 0.0;
	double total_material_weight = 0.0;
	for (int i = 0; i < nn; ++i) {
		total_material_weight += material_weight[i];
		if (raw_occupancy[i] > active_threshold) {
			active[i] = 1;
		}
		else {
			inactive_soft_resource += material_weight[i] * x[i];
		}
	}
	const double target_soft_resource = total_material_weight
		/ (1.0 + material_ratio_f);
	const double target_active_soft_resource =
		std::max(0.0, target_soft_resource - inactive_soft_resource);
	double lambda_low = 0.0;
	double lambda_high = 100000.0;
	while ((lambda_high - lambda_low) / (lambda_low + lambda_high) > 1e-4) {
		const double lmid = 0.5 * (lambda_low + lambda_high);
		double active_soft_resource = 0.0;
		for (int i = 0; i < nn; ++i) {
			if (!active[i]) {
				continue;
			}
			const double positive_dc = std::max(
				dc[i],
				std::numeric_limits<double>::min());
			const double candidate = x[i] * std::sqrt(positive_dc / lmid);
			xnew[i] = std::max(rho_min, std::max(x[i] - move_limit, std::min(1.0, std::min(x[i] + move_limit, candidate))));
			active_soft_resource += material_weight[i] * xnew[i];
		}
		if (active_soft_resource > target_active_soft_resource) {
			lambda_low = lmid;
		}
		else {
			lambda_high = lmid;
		}
	}
	x = xnew;
}

std::vector<MMCComponent> MakeInitialComponents(
	const std::vector<double>& xyz_range,
	double x_interval,
	double y_interval,
	double z_interval,
	const std::vector<double>& base_component,
	int type) {
	std::vector<MMCComponent> components;
	if (xyz_range.size() < 6 || base_component.size() < 6 ||
		x_interval <= 0.0 || y_interval <= 0.0 || z_interval <= 0.0) {
		return components;
	}

	std::vector<double> x_samples;
	std::vector<double> y_samples;
	std::vector<double> z_samples;
	const int first_x = static_cast<int>(xyz_range[0] + x_interval);
	const int first_y = static_cast<int>(xyz_range[2] + y_interval);
	const int first_z = static_cast<int>(xyz_range[4] + z_interval);
	for (double x = first_x; x <= xyz_range[1]; x += 2.0 * x_interval) x_samples.push_back(x);
	for (double y = first_y; y <= xyz_range[3]; y += 2.0 * y_interval) y_samples.push_back(y);
	for (double z = first_z; z <= xyz_range[5]; z += 2.0 * z_interval) z_samples.push_back(z);

	constexpr int pattern_count = 8;
	const double alpha_signs[pattern_count] = {
		1.0, 1.0, 1.0, 1.0, -1.0, -1.0, -1.0, -1.0
	};
	const double beta_signs[pattern_count] = {
		1.0, -1.0, 1.0, -1.0, 1.0, -1.0, 1.0, -1.0
	};
	const double gamma_signs[pattern_count] = {
		1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0, -1.0
	};

	const double l1 = base_component[0];
	const double l2 = base_component[1];
	const double l3 = base_component[2];
	const double alpha = base_component[3];
	const double beta = base_component[4];
	const double gamma = base_component[5];

	for (double z0 : z_samples) {
		for (double y0 : y_samples) {
			for (double x0 : x_samples) {
				if (type >= 5 && type <= 9) {
					components.emplace_back(x0, y0, z0, l1, l2, l3, alpha, beta, gamma);
					components.emplace_back(x0, y0, z0, l1, l2, l3, alpha, -beta, gamma);
					components.emplace_back(x0, y0, z0, l1, l2, l3, alpha, beta, gamma);
					components.emplace_back(x0, y0, z0, l1, l2, l3, alpha, -beta, gamma);
				}
				else if (type == 10 || type == 11) {
					for (int pattern = 0; pattern < pattern_count; ++pattern) {
						for (int duplicate = 0; duplicate < 2; ++duplicate) {
							(void)duplicate;
							components.emplace_back(
								x0,
								y0,
								z0,
								l1,
								l2,
								l3,
								alpha * alpha_signs[pattern],
								beta * beta_signs[pattern],
								gamma * gamma_signs[pattern]);
						}
					}
				}
				else if (type != -1) {
					components.emplace_back(x0, y0, z0, l1, l2, l3, alpha, beta, gamma);
					for (int duplicate = 1; duplicate < 8; ++duplicate) {
						(void)duplicate;
						components.emplace_back(x0, y0, z0, l1, l2, l3, -alpha, beta, gamma);
					}
				}
				else {
					for (int pattern = 0; pattern < pattern_count; ++pattern) {
						for (int duplicate = 0; duplicate < 2; ++duplicate) {
							(void)duplicate;
							components.emplace_back(
								x0,
								y0,
								z0,
								l1,
								l2,
								l3,
								alpha * alpha_signs[pattern],
								beta * beta_signs[pattern],
								gamma * gamma_signs[pattern]);
						}
					}
				}
			}
		}
	}
	return components;
}

std::vector<MMCComponent> InitializeComponentsForType(
	const std::vector<double>& xyz_range,
	int mmc_type) {
	std::vector<double> base_component;
	double x_interval = 0.0;
	double y_interval = 0.0;
	double z_interval = 0.0;
	int initialization_type = mmc_type;

	switch (mmc_type) {
	case 1:
		base_component = { 0.6, 7.0, 0.8, std::atan(1.0), 0.0, 0.0 };
		x_interval = 0.5;
		y_interval = 5.0;
		z_interval = 5.0;
		break;
	case 2:
		base_component = { 0.6, 15.0, 1.0, std::atan(0.5), 0.0, 0.0 };
		x_interval = 0.5;
		y_interval = 10.0;
		z_interval = 5.0;
		break;
	case 3:
		base_component = { 0.6, 15.0, 1.0, std::atan(1.0), 0.0, 0.0 };
		x_interval = 0.5;
		y_interval = 10.0;
		z_interval = 10.0;
		break;
	case 4:
		base_component = { 24.0, 5.0, 0.6, 0.0, std::atan(1.0), 0.0 };
		x_interval = 16.0;
		y_interval = 0.5;
		z_interval = 16.0;
		break;
	case 5:
		base_component = { 6.0, 1.25, 1.0, 0.0, std::atan(1.0), 0.0 };
		x_interval = 4.0;
		y_interval = 2.0;
		z_interval = 4.0;
		break;
	case 6:
		base_component = { 12.0, 2.5, 2.0, 0.0, std::atan(1.0), 0.0 };
		x_interval = 8.0;
		y_interval = 2.0;
		z_interval = 8.0;
		break;
	case 7:
		base_component = { 12.0, 2.5, 2.0, 0.0, std::atan(1.0), 0.0 };
		x_interval = 8.0;
		y_interval = 0.5;
		z_interval = 8.0;
		break;
	case 8:
		base_component = { 16.0, 3.5, 3.0, 0.0, std::atan(1.0), 0.0 };
		x_interval = 10.0;
		y_interval = 0.5;
		z_interval = 10.0;
		break;
	case 9:
		base_component = { 24.0, 5.0, 4.0, 0.0, std::atan(1.0), 0.0 };
		x_interval = 16.0;
		y_interval = 0.5;
		z_interval = 16.0;
		break;
	case 10:
		base_component = { 12.0, 2.0, 2.0, std::atan(1.0), std::atan(1.0), std::atan(1.0) };
		x_interval = 8.0;
		y_interval = 8.0;
		z_interval = 8.0;
		break;
	case 11:
		base_component = { 9.0, 2.0, 2.0, std::atan(1.0), std::atan(1.0), std::atan(1.0) };
		x_interval = 6.0;
		y_interval = 6.0;
		z_interval = 6.0;
		break;
	default:
		base_component = { 0.6, 8.0, 1.0, std::atan(1.0), 0.0, 0.0 };
		x_interval = 0.5;
		y_interval = 5.0;
		z_interval = 5.0;
		initialization_type = -1;
		break;
	}

	return MakeInitialComponents(
		xyz_range,
		x_interval,
		y_interval,
		z_interval,
		base_component,
		initialization_type);
}

void calc_Phi_Node(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const double minL,
	std::vector<MMC_Component> components,
	std::vector<double>& componentPhi,
	std::vector<std::vector<double>>& componentPhiDrv,
	int i,
	int p,
	double epsilon) {
	// 复用独立 3D MMC 的组件 TDF 逻辑。对每个顶点同时返回 phi_i
	// 以及对 9 个组件变量的解析导数，后续不再对这一段做数值差分。
	(void)minL;
	(void)epsilon;
	const int nVertex = static_cast<int>(allVertex.size());
	componentPhi.assign(nVertex, 0.0);
	componentPhiDrv.assign(
		nVertex,
		std::vector<double>(kVariablesPerComponent, 0.0));

	const MMC_Component& component = components[i];
	const double x0 = component.center_x;
	const double y0 = component.center_y;
	const double z0 = component.center_z;
	const double l1 = component.length_x + std::numeric_limits<double>::epsilon();
	const double l2 = component.length_y + std::numeric_limits<double>::epsilon();
	const double l3 = component.length_z + std::numeric_limits<double>::epsilon();
	const double alp = component.rotation_alpha;
	const double bet = component.rotation_beta;
	const double gam = component.rotation_gamma;

	const double sa = std::sin(alp), ca = std::cos(alp);
	const double sb = std::sin(bet), cb = std::cos(bet);
	const double sg = std::sin(gam), cg = std::cos(gam);

	Eigen::Matrix3d R;
	R <<
		cb * cg, cb * sg, -sb,
		sa * sb * cg - ca * sg, sa * sb * sg + ca * cg, sa * cb,
		ca * sb * cg + sa * sg, ca * sb * sg - sa * cg, ca * cb;

	std::vector<Eigen::Vector3d> lc_vec(nVertex);
	std::vector<double> x1_vec(nVertex), y1_vec(nVertex), z1_vec(nVertex);
	std::vector<double> temp_vec(nVertex);

	for (const auto& vertexPair : allVertex) {
		const int vertexTag = vertexPair.first;
		const Eigen::Vector3d lc =
			vertexPair.second - Eigen::Vector3d(x0, y0, z0);
		const Eigen::Vector3d xyz = R * lc;
		const double x1 = xyz.x() + kTiny;
		const double y1 = xyz.y() + kTiny;
		const double z1 = xyz.z() + kTiny;
		const double temp =
			std::pow(x1 / l1, p)
			+ std::pow(y1 / l2, p)
			+ std::pow(z1 / l3, p);

		lc_vec[vertexTag] = lc;
		x1_vec[vertexTag] = x1;
		y1_vec[vertexTag] = y1;
		z1_vec[vertexTag] = z1;
		temp_vec[vertexTag] = temp;
		componentPhi[vertexTag] =
			1.0 - std::pow(temp, 1.0 / static_cast<double>(p));
	}

	Eigen::Matrix3d Ra, Rb, Rg;
	Ra <<
		0.0, 0.0, 0.0,
		R(2, 0), R(2, 1), R(2, 2),
		-R(1, 0), -R(1, 1), -R(1, 2);
	Rb <<
		-sb * cg, -sb * sg, -cb,
		sa * cb * cg, sa * cb * sg, -sa * sb,
		ca * cb * cg, ca * cb * sg, -ca * sb;
	Rg <<
		-cb * sg, cb * cg, 0.0,
		-sa * sb * sg - ca * cg, sa * sb * cg - ca * sg,
		0.0,
		-ca * sb * sg + sa * cg, ca * sb * cg + sa * sg,
		0.0;

	for (const auto& vertexPair : allVertex) {
		const int vertexTag = vertexPair.first;
		const double x1 = x1_vec[vertexTag];
		const double y1 = y1_vec[vertexTag];
		const double z1 = z1_vec[vertexTag];
		const double temp = temp_vec[vertexTag];
		const Eigen::Vector3d& lc = lc_vec[vertexTag];

		const double common_power =
			std::pow(temp, 1.0 / static_cast<double>(p) - 1.0);
		const double temp1 =
			common_power * std::pow(x1 / l1, p - 1);
		const double temp2 =
			common_power * std::pow(y1 / l2, p - 1);
		const double temp3 =
			common_power * std::pow(z1 / l3, p - 1);

		const double dpdx1 = temp1 / l1;
		const double dpdy1 = temp2 / l2;
		const double dpdz1 = temp3 / l3;
		const double dpdl1 = temp1 * x1 / (l1 * l1);
		const double dpdl2 = temp2 * y1 / (l2 * l2);
		const double dpdl3 = temp3 * z1 / (l3 * l3);

		const double dx1_dalpha = (Ra.row(0) * lc)(0);
		const double dy1_dalpha = (Ra.row(1) * lc)(0);
		const double dz1_dalpha = (Ra.row(2) * lc)(0);
		const double dx1_dbeta = (Rb.row(0) * lc)(0);
		const double dy1_dbeta = (Rb.row(1) * lc)(0);
		const double dz1_dbeta = (Rb.row(2) * lc)(0);
		const double dx1_dgamma = (Rg.row(0) * lc)(0);
		const double dy1_dgamma = (Rg.row(1) * lc)(0);
		const double dz1_dgamma = (Rg.row(2) * lc)(0);

		std::vector<double>& drv = componentPhiDrv[vertexTag];
		drv[0] = dpdx1 * R(0, 0) + dpdy1 * R(1, 0) + dpdz1 * R(2, 0);
		drv[1] = dpdx1 * R(0, 1) + dpdy1 * R(1, 1) + dpdz1 * R(2, 1);
		drv[2] = dpdx1 * R(0, 2) + dpdy1 * R(1, 2) + dpdz1 * R(2, 2);
		drv[3] = dpdl1;
		drv[4] = dpdl2;
		drv[5] = dpdl3;
		drv[6] = -(dpdx1 * dx1_dalpha + dpdy1 * dy1_dalpha + dpdz1 * dz1_dalpha);
		drv[7] = -(dpdx1 * dx1_dbeta + dpdy1 * dy1_dbeta + dpdz1 * dz1_dbeta);
		drv[8] = -(dpdx1 * dx1_dgamma + dpdy1 * dy1_dgamma + dpdz1 * dz1_dgamma);
	}
}

void ComputeGlobalTDF(
	const std::vector<std::vector<double>>& allPhis,
	double lmd,
	std::vector<double>& globalPhi,
	std::vector<std::vector<double>>& phimaxDphi) {
	// 用 KS 函数平滑近似 max_i(phi_i)，从而表达多个 MMC 组件的并集：
	//   Phi = log(sum_i exp(lambda*phi_i)) / lambda
	//   dPhi/dphi_i = exp(lambda*phi_i) / sum_j exp(lambda*phi_j)
	//
	// 注意：分母的 +1e-9 是从独立 MMC 实际主循环复用的稳定项。
	// 它会使解析导数成为严格 KS 导数的轻微近似；全有限差分日志用于量化该偏差。
	const int nComp = static_cast<int>(allPhis.size());
	if (nComp == 0) {
		throw std::runtime_error("The MMC component set is empty.");
	}
	const int nGrid = static_cast<int>(allPhis.front().size());
	globalPhi.assign(nGrid, 0.0);
	phimaxDphi.assign(
		nGrid,
		std::vector<double>(nComp, 0.0));
	const double tiny = 1e-9;

	for (int g = 0; g < nGrid; ++g) {
		std::vector<double> expValues(nComp);
		double sumExp = 0.0;
		for (int c = 0; c < nComp; ++c) {
			expValues[c] = std::exp(lmd * allPhis[c][g]);
			sumExp += expValues[c];
		}
		globalPhi[g] = std::max(-1e3, std::log(sumExp) / lmd);
		for (int c = 0; c < nComp; ++c) {
			phimaxDphi[g][c] = expValues[c] / (sumExp + tiny);
		}
	}
}

GeometryEvaluation EvaluateGeometry(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::vector<MMCComponent>& mmcComponents,
	const MMCMultiMaterialSimultaneousParams& params) {
	// 完整几何映射：
	//   d_ir -> phi_i -> KS全局Phi -> H(Phi) -> 单元占有率m_e。
	// 本函数只计算几何场和链式因子，不进行有限元求解。
	GeometryEvaluation evaluation;
	const int nComp = static_cast<int>(mmcComponents.size());
	const int nVertex = static_cast<int>(allVertex.size());
	const int nGrid = static_cast<int>(cell2V.size());

	evaluation.allPhis.resize(nComp);
	evaluation.allPhiDrv.resize(nComp);
	for (int i = 0; i < nComp; ++i) {
		std::vector<double> componentPhi;
		std::vector<std::vector<double>> componentPhiDrv;
		calc_Phi_Node(
			allVertex,
			0.0,
			mmcComponents,
			componentPhi,
			componentPhiDrv,
			i,
			kTdfPower,
			params.heaviside_epsilon);
		evaluation.allPhiDrv[i] = componentPhiDrv;
		evaluation.allPhis[i] = componentPhi;
	}

	ComputeGlobalTDF(
		evaluation.allPhis,
		params.ks_lambda,
		evaluation.globalPhi,
		evaluation.phimaxDphi);

	evaluation.vertexH.resize(nVertex, 0.0);
	evaluation.delta_H.resize(nVertex, 0.0);
	for (int i = 0; i < nVertex; ++i) {
		evaluation.vertexH[i] = Heaviside(
			evaluation.globalPhi[i],
			params.heaviside_alpha,
			params.heaviside_epsilon);
		evaluation.delta_H[i] = HeavisideDerivative(
			evaluation.globalPhi[i],
			params.heaviside_alpha,
			params.heaviside_epsilon);
	}

	// 六面体单元占有率取 8 个顶点 Heaviside 值的算术平均：
	//   m_e = (1/8) sum_{v in e} H(Phi_v)。
	evaluation.density.assign(nGrid, 0.0);
	for (int i = 0; i < nGrid; ++i) {
		for (const auto& vertex : cell2V.at(i)) {
			evaluation.density[i] += evaluation.vertexH[vertex.first];
		}
		evaluation.density[i] /= 8.0;
	}
	if (!IsFiniteVector(evaluation.globalPhi) ||
		!IsFiniteVector(evaluation.vertexH) ||
		!IsFiniteVector(evaluation.density)) {
		throw std::runtime_error("MMC TDF, Heaviside, or cell occupancy contains a non-finite value.");
	}
	return evaluation;
}

std::vector<double> ComputeDomainRange(
	const std::map<int, Eigen::Vector3d>& all_vertex) {
	double xmin = std::numeric_limits<double>::max();
	double xmax = std::numeric_limits<double>::lowest();
	double ymin = std::numeric_limits<double>::max();
	double ymax = std::numeric_limits<double>::lowest();
	double zmin = std::numeric_limits<double>::max();
	double zmax = std::numeric_limits<double>::lowest();
	for (const auto& vertex : all_vertex) {
		xmin = std::min(xmin, vertex.second.x());
		xmax = std::max(xmax, vertex.second.x());
		ymin = std::min(ymin, vertex.second.y());
		ymax = std::max(ymax, vertex.second.y());
		zmin = std::min(zmin, vertex.second.z());
		zmax = std::max(zmax, vertex.second.z());
	}
	if (!(xmax > xmin) || !(ymax > ymin) || !(zmax > zmin)) {
		throw std::runtime_error("All three MMC design-domain dimensions must be positive.");
	}
	return { xmin, xmax, ymin, ymax, zmin, zmax };
}

void BuildPhysicalComponentBounds(
	const std::vector<double>& xyz_range,
	int component_count,
	std::vector<double>& lower_bounds,
	std::vector<double>& upper_bounds) {
	const double dx = xyz_range[1] - xyz_range[0];
	const double dy = xyz_range[3] - xyz_range[2];
	const double dz = xyz_range[5] - xyz_range[4];
	const double maximum_length =
		0.5 * std::sqrt(dx * dx + dy * dy + dz * dz);
	const std::vector<double> single_lower = {
		xyz_range[0], xyz_range[2], xyz_range[4],
		0.0, 0.0, 0.0,
		-kPi, -kPi, -kPi
	};
	const std::vector<double> single_upper = {
		xyz_range[1], xyz_range[3], xyz_range[5],
		maximum_length, maximum_length, maximum_length,
		kPi, kPi, kPi
	};
	lower_bounds.resize(component_count * kVariablesPerComponent);
	upper_bounds.resize(component_count * kVariablesPerComponent);
	for (int component_id = 0; component_id < component_count; ++component_id) {
		std::copy(
			single_lower.begin(),
			single_lower.end(),
			lower_bounds.begin() + component_id * kVariablesPerComponent);
		std::copy(
			single_upper.begin(),
			single_upper.end(),
			upper_bounds.begin() + component_id * kVariablesPerComponent);
	}
}

void FlattenPhysicalComponents(
	const std::vector<MMCComponent>& components,
	std::vector<double>& physical_variables) {
	physical_variables.clear();
	physical_variables.reserve(components.size() * kVariablesPerComponent);
	for (const MMCComponent& component : components) {
		physical_variables.push_back(component.center_x);
		physical_variables.push_back(component.center_y);
		physical_variables.push_back(component.center_z);
		physical_variables.push_back(component.length_x);
		physical_variables.push_back(component.length_y);
		physical_variables.push_back(component.length_z);
		physical_variables.push_back(component.rotation_alpha);
		physical_variables.push_back(component.rotation_beta);
		physical_variables.push_back(component.rotation_gamma);
	}
}

// 与独立 MMC 一样直接保存和更新九个物理设计变量。
void UpdateComponentsFromPhysical(
	std::vector<MMCComponent>& components,
	const std::vector<double>& physical_variables) {
	if (physical_variables.size()
		!= components.size() * kVariablesPerComponent) {
		throw std::runtime_error(
			"The physical MMC design-variable count is invalid.");
	}
	int variable = 0;
	for (MMCComponent& component : components) {
		component.center_x = physical_variables[variable++];
		component.center_y = physical_variables[variable++];
		component.center_z = physical_variables[variable++];
		component.length_x = physical_variables[variable++];
		component.length_y = physical_variables[variable++];
		component.length_z = physical_variables[variable++];
		component.rotation_alpha = physical_variables[variable++];
		component.rotation_beta = physical_variables[variable++];
		component.rotation_gamma = physical_variables[variable++];
	}
}

void ComputeMmcDerivativeFromElementWeights(
	const GeometryEvaluation& geometry,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<double>& element_weights,
	std::vector<double>& gradient,
	std::vector<double>* maximum_absolute_contribution = nullptr,
	std::vector<double>* sum_absolute_contribution = nullptr) {
	// 将任意“对单元占有率 m_e 的导数权重”继续沿 MMC 几何链反传。
	// 若目标量写为 J(m)，调用方传入 element_weights[e] = dJ/dm_e，
	// 本函数计算：
	//   dJ/dd_ir = sum_e (dJ/dm_e) (dm_e/dd_ir)
	//             = sum_v w_v H'(Phi_v) KS'_i(v) phi'_{ir}(v)。
	//
	// 三种主要调用的 element_weights 含义：
	//   柔顺度：-E_e q_e
	//   MMC体积：1 / N
	//   比例残差：((1-eta_e)-f eta_e) / N
	const int nComp =
		static_cast<int>(geometry.allPhiDrv.size());
	const int nVertex =
		static_cast<int>(geometry.vertexH.size());
	const int nGrid = static_cast<int>(cell_to_vertex.size());
	if (static_cast<int>(element_weights.size()) != nGrid) {
		throw std::runtime_error("The element-weight count does not match the MMC grid.");
	}

	// m_e 是 8 个顶点值的平均，因此先把单元权重按 1/8 累加回顶点。
	// 这与原 MMC 的 CeNode / volNod 累加方式一致；这里只把两者统一成
	// elementWeightNode，便于同一段代码服务目标、体积和材料比例三个量。
	std::vector<double> elementWeightNode(nVertex, 0.0);
	for (int cell_id = 0; cell_id < nGrid; ++cell_id) {
		for (const auto& vertex : cell_to_vertex.at(cell_id)) {
			elementWeightNode[vertex.first] += element_weights[cell_id] / 8.0;
		}
	}

	// 与独立 MMC 主循环相同：先计算
	// PhimaxDrvAct[c][v][d] = phimaxDphi[v][c] * allPhiDrv[c][v][d]。
	std::vector<std::vector<std::vector<double>>> PhimaxDrvAct(
		nComp,
		std::vector<std::vector<double>>(
			nVertex,
			std::vector<double>(kVariablesPerComponent, 0.0)));
	for (int c = 0; c < nComp; ++c) {
		for (int v = 0; v < nVertex; ++v) {
			for (int d = 0; d < kVariablesPerComponent; ++d) {
				PhimaxDrvAct[c][v][d] =
					geometry.phimaxDphi[v][c] * geometry.allPhiDrv[c][v][d];
			}
		}
	}

	const int variable_count = nComp * kVariablesPerComponent;
	gradient.assign(variable_count, 0.0);
	if (maximum_absolute_contribution != nullptr) {
		maximum_absolute_contribution->assign(variable_count, 0.0);
	}
	if (sum_absolute_contribution != nullptr) {
		sum_absolute_contribution->assign(variable_count, 0.0);
	}

	for (int component_id = 0; component_id < nComp; ++component_id) {
		for (int local_variable = 0;
			local_variable < kVariablesPerComponent;
			++local_variable) {
			const int variable =
				component_id * kVariablesPerComponent + local_variable;
			double derivative = 0.0;
			double max_abs = 0.0;
			double sum_abs = 0.0;
			for (int vertex_id = 0; vertex_id < nVertex; ++vertex_id) {
				// 显式写出链式乘积：dJ/dm -> dm/dH -> dH/dPhi
				// -> dPhi/dphi_i -> dphi_i/dd_ir。
				const double contribution =
					elementWeightNode[vertex_id]
					* geometry.delta_H[vertex_id]
					* PhimaxDrvAct[component_id][vertex_id][local_variable];
				derivative += contribution;
				max_abs = std::max(max_abs, std::abs(contribution));
				sum_abs += std::abs(contribution);
			}
			gradient[variable] = derivative;
			if (maximum_absolute_contribution != nullptr) {
				(*maximum_absolute_contribution)[variable] = max_abs;
			}
			if (sum_absolute_contribution != nullptr) {
				(*sum_absolute_contribution)[variable] = sum_abs;
			}
		}
	}
	if (!IsFiniteVector(gradient)) {
		throw std::runtime_error("An analytical MMC sensitivity is non-finite.");
	}
}

void ReadBoundaryFile(
	const std::string& boundary_file,
	std::vector<CPoint>& boundary_points,
	std::vector<int>& boundary_types,
	std::vector<CVector>& boundary_values) {
	std::ifstream input(boundary_file);
	if (!input.good()) {
		throw std::runtime_error("Cannot open the boundary-condition file.");
	}
	int boundary_count = 0;
	input >> boundary_count;
	if (!input.good() || boundary_count <= 0) {
		throw std::runtime_error("The boundary-condition file is empty or invalid.");
	}
	boundary_points.reserve(boundary_count);
	boundary_types.reserve(boundary_count);
	boundary_values.reserve(boundary_count);
	for (int boundary_id = 0; boundary_id < boundary_count; ++boundary_id) {
		double x = 0.0;
		double y = 0.0;
		double z = 0.0;
		int type = 0;
		double value_x = 0.0;
		double value_y = 0.0;
		double value_z = 0.0;
		input >> x >> y >> z >> type >> value_x >> value_y >> value_z;
		if (!input.good()) {
			throw std::runtime_error("A boundary-condition record is incomplete.");
		}
		boundary_points.emplace_back(x, y, z);
		boundary_types.push_back(type);
		boundary_values.emplace_back(value_x, value_y, value_z);
	}
}

void BuildAnalysisBoundaryConditions(
	const std::vector<CBSplineSolid>& patch,
	gpmesh::GPSplineSolidV0& spline,
	const std::vector<CPoint>& boundary_points,
	const std::vector<int>& boundary_types,
	const std::vector<CVector>& boundary_values,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points) {
	std::vector<Point> source_points;
	std::vector<Point> target_points;
	source_points.reserve(boundary_points.size());
	for (const CPoint& point : boundary_points) {
		source_points.push_back(Point{ point.getX(), point.getY(), point.getZ() });
	}
	for (const CBSplineSolid& solid : patch) {
		for (const CPoint& point : solid.coefs()) {
			target_points.push_back(Point{ point.getX(), point.getY(), point.getZ() });
		}
	}

	std::vector<int> nearest_ids;
	GetNearestIds(source_points, target_points, nearest_ids);
	const auto& local_to_global = spline.m_sparse_system.get_l2g_map();
	for (int boundary_id = 0;
		boundary_id < static_cast<int>(boundary_points.size());
		++boundary_id) {
		const int patch_id = nearest_ids[boundary_id] / 64;
		const int local_id = nearest_ids[boundary_id] % 64;
		if (patch_id < 0 || patch_id >= static_cast<int>(local_to_global.size()) ||
			local_id < 0 || local_id >= static_cast<int>(local_to_global[patch_id].size())) {
			throw std::runtime_error("A boundary point cannot be mapped to a valid control point.");
		}
		const int global_id = local_to_global[patch_id][local_id];
		if (global_id < 0) {
			std::cerr
				<< "[MMC-MultiMaterial-Simultaneous] Skipping an unmapped boundary control point."
				<< std::endl;
			continue;
		}
		if (boundary_types[boundary_id] == 2) {
			neumann_points.emplace_back(
				condition_type::neumann,
				boundary_values[boundary_id],
				global_id);
		}
		else if (boundary_types[boundary_id] == 1) {
			dirichlet_points.emplace_back(
				condition_type::dirichlet,
				CVector(0.0, 0.0, 0.0),
				global_id);
		}
	}
	if (neumann_points.empty() || dirichlet_points.empty()) {
		throw std::runtime_error("Simultaneous optimization requires at least one load point and one constrained point.");
	}
}

void BuildDisplacementByPatch(
	const gpmesh::GPSplineSolidV0& spline,
	const Eigen::MatrixXd& displacement,
	Eigen::MatrixXd& displacement_by_patch) {
	const int patch_count = spline.m_num_geometry;
	if (displacement.cols() != 1 ||
		displacement.rows() != 3 * spline.m_sparse_system.getControlNums()) {
		throw std::runtime_error("The finite-element displacement vector has an invalid size.");
	}
	displacement_by_patch.resize(patch_count, 64 * 3);
	for (int patch_id = 0; patch_id < patch_count; ++patch_id) {
		for (int local_id = 0; local_id < 64; ++local_id) {
			const int global_id =
				spline.m_sparse_system.get_l2g_map()[patch_id][local_id];
			if (global_id < 0) {
				throw std::runtime_error("A local control point has no valid global id.");
			}
			displacement_by_patch(patch_id, local_id) =
				displacement(global_id, 0);
			displacement_by_patch(patch_id, local_id + 64) =
				displacement(
					global_id + spline.m_sparse_system.getControlNums(),
					0);
			displacement_by_patch(patch_id, local_id + 128) =
				displacement(
					global_id + 2 * spline.m_sparse_system.getControlNums(),
					0);
		}
	}
}

double RelativeDifference(double first, double second) {
	return std::abs(first - second)
		/ std::max({ std::abs(first), std::abs(second), 1e-12 });
}

void RunClosedFormSelfChecks(
	const MMCMultiMaterialSimultaneousParams& params) {
	const double modulus_at_zero = InterpolateYoungsModulus(0.0, params);
	const double modulus_at_one = InterpolateYoungsModulus(1.0, params);
	const double middle_eta = 0.5;
	const double expected_middle =
		std::pow(middle_eta, params.material_penalty) * params.E1
		+ (1.0 - std::pow(middle_eta, params.material_penalty)) * params.E2;
	if (RelativeDifference(modulus_at_zero, params.E2) > 1e-13 ||
		RelativeDifference(modulus_at_one, params.E1) > 1e-13 ||
		RelativeDifference(
			InterpolateYoungsModulus(middle_eta, params),
			expected_middle) > 1e-13) {
		throw std::runtime_error("The closed-form endpoint/interior interpolation check failed.");
	}

	for (double sample_eta : { 0.0, 0.37, 1.0 }) {
		const double expected_derivative =
			params.material_penalty
			* std::pow(sample_eta, params.material_penalty - 1.0)
			* (params.E1 - params.E2);
		if (RelativeDifference(
			InterpolateYoungsModulusDerivative(sample_eta, params),
			expected_derivative) > 1e-13) {
			throw std::runtime_error(
				"The closed-form material interpolation derivative check failed.");
		}
	}

	// ProjectOccupancy() 已不在主模型中使用；下面用固定 m_e 样本核对
	// 材料比例约束的原始 MMC 占有率链。
	const double sample_m = 0.37;
	const double normalized_cell_volume = 0.125;
	const double ratio_material_derivative =
		-(1.0 + params.material_ratio_f)
		* sample_m
		* normalized_cell_volume;
	const double independently_expanded_derivative =
		-sample_m
		* normalized_cell_volume
		- params.material_ratio_f
		* sample_m
		* normalized_cell_volume;
	if (RelativeDifference(
		ratio_material_derivative,
		independently_expanded_derivative) > 1e-13) {
		throw std::runtime_error("The closed-form material-ratio derivative check failed.");
	}

	// 用局部中心差分同时核对 R 对 m 与 eta 的两条链，防止遗漏
	// m_e 权重或材料 eta 权重。
	const double sample_eta = 0.41;
	const double ratio_coefficient =
		(1.0 - sample_eta) - params.material_ratio_f * sample_eta;
	const auto local_ratio = [&](double occupancy, double eta) {
		return occupancy * ((1.0 - eta) - params.material_ratio_f * eta)
			* normalized_cell_volume;
	};
	const double finite_difference_step = 1e-6;
	const double ratio_mmc_derivative = ratio_coefficient
		* normalized_cell_volume;
	const double ratio_mmc_finite_difference =
		(local_ratio(sample_m + finite_difference_step, sample_eta)
			- local_ratio(sample_m - finite_difference_step, sample_eta))
		/ (2.0 * finite_difference_step);
	const double ratio_material_finite_difference =
		(local_ratio(sample_m, sample_eta + finite_difference_step)
			- local_ratio(sample_m, sample_eta - finite_difference_step))
		/ (2.0 * finite_difference_step);
	if (RelativeDifference(
		ratio_mmc_derivative,
		ratio_mmc_finite_difference) > 1e-7 ||
		RelativeDifference(
			ratio_material_derivative,
			ratio_material_finite_difference) > 1e-7) {
		throw std::runtime_error(
			"The projected material-ratio finite-difference self-check failed.");
	}

}

// 同步优化核心逻辑
SimultaneousEvaluation EvaluateSimultaneousState(
	gpmesh::GPSplineSolidV0& spline,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points,
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<MMCComponent>& components,
	const std::vector<double>& material_eta,
	const std::vector<Eigen::MatrixXd>& unit_local_matrix,
	const MMCMultiMaterialSimultaneousParams& params,
	int projection_iteration = 1) {
	// 单轮联合状态计算器。给定当前 MMC 组件和完整网格 eta 向量，按以下顺序：
	//   1. 计算 MMC 单元占有率 m_e；
	//   2. 计算 E_e(eta_e) 和 a_e=m_e E_e；
	//   3. 只求解一次 K(a)U=F；
	//   4. 由同一个 U/q_e 计算 MMC 与材料两组解析梯度；
	//   5. 计算几何体积和材料比例约束及其两块梯度。
	// 该函数不更新 MMA，也不修改设计变量，因此可安全用于有限差分重求解。
	const int cell_count = static_cast<int>(cell_to_vertex.size());
	if (static_cast<int>(material_eta.size()) != cell_count ||
		static_cast<int>(unit_local_matrix.size()) != cell_count ||
		spline.m_num_geometry != cell_count) {
		throw std::runtime_error("The simultaneous state, local matrix, and IGA cell counts differ.");
	}

	SimultaneousEvaluation evaluation;
	evaluation.geometry = EvaluateGeometry(
		all_vertex,
		cell_to_vertex,
			components,
	params);
	// 为了与纯 MMC 模板保持一致，有限元装配直接使用原始占有率 m_e。
	// projected_occupancy 仅作为兼容旧导出字段的拷贝，不再参与投影模型。
	evaluation.material_modulus.resize(cell_count, 0.0);
	evaluation.material_modulus_derivative.resize(cell_count, 0.0);
	evaluation.projected_occupancy.resize(cell_count, 0.0);
	evaluation.projected_occupancy_derivative.resize(cell_count, 0.0);
	evaluation.effective_stiffness.resize(cell_count, 0.0);
	evaluation.projection_beta = ProjectionBetaForIteration(
		params, projection_iteration);
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		evaluation.material_modulus[cell_id] = // \(\eta_e\)SIMP 多材料设计变量。决定当前单元偏向材料 E1 还是 E2。
			InterpolateYoungsModulus(material_eta[cell_id], params);
		evaluation.material_modulus_derivative[cell_id] = // 多材料插值导数
			InterpolateYoungsModulusDerivative(material_eta[cell_id], params);
		evaluation.projected_occupancy[cell_id] = // 旧字段，原来是投影
			evaluation.geometry.density[cell_id];
		evaluation.projected_occupancy_derivative[cell_id] = 1.0; // 旧字段，原来是投影
		evaluation.effective_stiffness[cell_id] = // mmc占有率 * simp多材料插值密度
			evaluation.geometry.density[cell_id]
			* evaluation.material_modulus[cell_id];
	}
	if (!IsFiniteVector(evaluation.material_modulus) ||
		!IsFiniteVector(evaluation.projected_occupancy) ||
		!IsFiniteVector(evaluation.projected_occupancy_derivative) ||
		!IsFiniteVector(evaluation.effective_stiffness)) {
		throw std::runtime_error("A material modulus or effective stiffness value is non-finite.");
	}

	// 共享结构分析：当前联合状态只在此处求解一次位移。
	spline.Solver_With_Density_MMC(
		neumann_points,
		dirichlet_points,
		evaluation.effective_stiffness,
		evaluation.displacement);
	if (!evaluation.displacement.allFinite()) {
		throw std::runtime_error("The shared finite-element solve returned a non-finite displacement.");
	}

	// 和纯MMC一样，把全局位移向量 U 重新整理成“每个单元/patch 一行”的局部位移矩阵。
	Eigen::MatrixXd displacement_by_patch;
	BuildDisplacementByPatch(
		spline,
		evaluation.displacement,
		displacement_by_patch);

	// q_e 使用单位模量局部矩阵 K_e^0。目标值累加方式与纯 MMC 保持一致：C = sum_e q_e。
	evaluation.element_energy.resize(cell_count, 0.0);
	evaluation.compliance = 0.0;
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		const Eigen::VectorXd element_displacement =
			displacement_by_patch.row(cell_id);
		const double energy =
			(element_displacement.transpose()
				* unit_local_matrix[cell_id]
				* element_displacement)(0, 0);
		evaluation.element_energy[cell_id] = energy;
		evaluation.compliance +=
			energy;
	}

	// 保留旧诊断字段，但关闭 F^T U 尺度一致性检查；当前目标值按纯 MMC 的 sum_e q_e 记录。
	evaluation.compliance_from_load = evaluation.compliance;
	evaluation.compliance_relative_mismatch = 0.0;
	if (!std::isfinite(evaluation.compliance) ||
		!IsFiniteVector(evaluation.element_energy)) {
		throw std::runtime_error("The compliance or an element energy value is non-finite.");
	}

	// 柔顺度伴随公式：dC/dx = -U^T(dK/dx)U。
	// 对 MMC：dC/dm_e = -E_e q_e；
	// 对材料：dC/deta_e = -m_e q_e dE_e/deta_e。
	std::vector<double> objective_element_weights(cell_count, 0.0);
	std::vector<double> geometry_volume_weights(
		cell_count,
		1.0 / static_cast<double>(cell_count));
	evaluation.material_objective_gradient.resize(cell_count, 0.0);

	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		const double normalized_cell_volume =
			1.0 / static_cast<double>(cell_count);
		objective_element_weights[cell_id] =
			-evaluation.material_modulus[cell_id]
			* evaluation.element_energy[cell_id];
		evaluation.material_objective_gradient[cell_id] =
			-evaluation.geometry.density[cell_id]
			* evaluation.element_energy[cell_id]
			* evaluation.material_modulus_derivative[cell_id];

		// MMC 占有率加权的材料资源：
		//   V1 = sum_e m_e eta_e / N
		//   V2 = sum_e m_e (1-eta_e) / N
		//   R  = V2 - f V1。
		evaluation.V1 +=
			evaluation.geometry.density[cell_id]
			* material_eta[cell_id]
			* normalized_cell_volume;
		evaluation.V2 +=
			evaluation.geometry.density[cell_id]
			* (1.0 - material_eta[cell_id])
			* normalized_cell_volume;
	}

	// MMC 只反传目标梯度和几何体积约束；材料比例由材料 OC 负责，不再反传到 MMC。
	ComputeMmcDerivativeFromElementWeights(
		evaluation.geometry,
		cell_to_vertex,
		objective_element_weights,
		evaluation.mmc_objective_gradient,
		&evaluation.mmc_contribution_max_abs,
		&evaluation.mmc_contribution_sum_abs);
	ComputeMmcDerivativeFromElementWeights(
		evaluation.geometry,
		cell_to_vertex,
		geometry_volume_weights,
		evaluation.geometry_constraint_gradient);
	evaluation.ratio_constraint_mmc_gradient.assign(
		components.size() * kVariablesPerComponent,
		0.0);
	evaluation.ratio_constraint_material_gradient.assign(
		cell_count,
		0.0);

	evaluation.mmc_volume =
		std::accumulate(
			evaluation.geometry.density.begin(),
			evaluation.geometry.density.end(),
			0.0)
		/ static_cast<double>(cell_count);
	// 材料比例由 OC 二分直接维持；ratio 约束不再作为同步优化约束参与流程。
	evaluation.ratio_residual = 0.0;
	evaluation.positive_ratio_constraint = 0.0;
	evaluation.negative_ratio_constraint = 0.0;

	if (!IsFiniteVector(evaluation.material_objective_gradient) ||
		!IsFiniteVector(evaluation.geometry_constraint_gradient)) {
		throw std::runtime_error("An analytical objective or constraint sensitivity is non-finite.");
	}
	return evaluation;
}

std::vector<int> SelectFiniteDifferenceIndices(
	const std::vector<double>& variables,
	const std::vector<double>& derivative_scores,
	double lower_bound,
	double upper_bound,
	double step,
	int requested_count) {
	std::vector<int> candidates;
	if (variables.size() != derivative_scores.size() || requested_count <= 0) {
		return candidates;
	}
	for (int variable = 0;
		variable < static_cast<int>(variables.size());
		++variable) {
		if (variables[variable] - step >= lower_bound &&
			variables[variable] + step <= upper_bound &&
			std::isfinite(derivative_scores[variable])) {
			candidates.push_back(variable);
		}
	}
	std::stable_sort(
		candidates.begin(),
		candidates.end(),
		[&](int lhs, int rhs) {
			return std::abs(derivative_scores[lhs]) >
				std::abs(derivative_scores[rhs]);
		});
	if (static_cast<int>(candidates.size()) > requested_count) {
		candidates.resize(requested_count);
	}
	return candidates;
}

std::vector<int> SelectPhysicalMmcFiniteDifferenceIndices(
	const std::vector<double>& variables,
	const std::vector<double>& derivative_scores,
	const std::vector<double>& lower_bounds,
	const std::vector<double>& upper_bounds,
	double relative_step,
	int requested_count) {
	std::vector<int> candidates;
	if (variables.size() != derivative_scores.size() ||
		variables.size() != lower_bounds.size() ||
		variables.size() != upper_bounds.size() ||
		requested_count <= 0) {
		return candidates;
	}
	for (int variable = 0;
		variable < static_cast<int>(variables.size());
		++variable) {
		const double range = upper_bounds[variable] - lower_bounds[variable];
		const double step = relative_step * range;
		if (range > 0.0 && variables[variable] - step >= lower_bounds[variable] &&
			variables[variable] + step <= upper_bounds[variable] &&
			std::isfinite(derivative_scores[variable])) {
			candidates.push_back(variable);
		}
	}
	std::stable_sort(
		candidates.begin(),
		candidates.end(),
		[&](int lhs, int rhs) {
			return std::abs(derivative_scores[lhs]) >
				std::abs(derivative_scores[rhs]);
		});
	if (static_cast<int>(candidates.size()) > requested_count) {
		candidates.resize(requested_count);
	}
	return candidates;
}

void WriteFiniteDifferenceHeader(std::ofstream& output) {
	output
		<< "accepted_iteration\tvariable_block\tvariable_index"
		<< "\tdifference_scheme\tstep\tquantity"
		<< "\tbase_value\tplus_value\tminus_value"
		<< "\tanalytical_derivative\tfinite_difference_derivative"
		<< "\tabsolute_error\trelative_error\tstatus\tmessage\n";
}

void WriteFiniteDifferenceFailure(
	std::ofstream& output,
	int accepted_iteration,
	double step,
	const std::string& message) {
	std::string sanitized_message = message;
	std::replace(
		sanitized_message.begin(),
		sanitized_message.end(),
		'\t',
		' ');
	std::replace(
		sanitized_message.begin(),
		sanitized_message.end(),
		'\n',
		' ');
	const double not_a_number = std::numeric_limits<double>::quiet_NaN();
	output << accepted_iteration
		<< "\tdiagnostic\t-1\tnone\t" << step
		<< "\terror\t" << not_a_number
		<< '\t' << not_a_number
		<< '\t' << not_a_number
		<< '\t' << not_a_number
		<< '\t' << not_a_number
		<< '\t' << not_a_number
		<< '\t' << not_a_number
		<< "\tfailed\t" << sanitized_message << '\n';
	output.flush();
}

const char* ClassifyFiniteDifferenceError(
	double analytical_derivative,
	double finite_difference_derivative,
	double absolute_error,
	double relative_error) {
	const double derivative_scale = std::max(
		std::abs(analytical_derivative),
		std::abs(finite_difference_derivative));
	if (derivative_scale <= kFiniteDifferenceNearZeroTolerance) {
		return absolute_error <= kFiniteDifferencePassAbsoluteTolerance
			? "PASS"
			: "FAIL";
	}
	if (relative_error <= kFiniteDifferencePassRelativeTolerance) {
		return "PASS";
	}
	if (relative_error <= kFiniteDifferenceWarnRelativeTolerance) {
		return "WARN";
	}
	return "FAIL";
}

void WriteFiniteDifferenceQuantity(
	std::ofstream& output,
	int accepted_iteration,
	const char* variable_block,
	int variable_index,
	double step,
	const char* quantity,
	double base_value,
	double plus_value,
	double minus_value,
	double analytical_derivative) {
	const double finite_difference_derivative =
		(plus_value - minus_value) / (2.0 * step);
	const double absolute_error = std::abs(
		analytical_derivative - finite_difference_derivative);
	const double relative_error = absolute_error / std::max({
		std::abs(analytical_derivative),
		std::abs(finite_difference_derivative),
		1e-12
	});
	const char* status = ClassifyFiniteDifferenceError(
		analytical_derivative,
		finite_difference_derivative,
		absolute_error,
		relative_error);
	output << accepted_iteration << '\t'
		<< variable_block << '\t'
		<< variable_index << "\tcentral\t"
		<< step << '\t'
		<< quantity << '\t'
		<< base_value << '\t'
		<< plus_value << '\t'
		<< minus_value << '\t'
		<< analytical_derivative << '\t'
		<< finite_difference_derivative << '\t'
		<< absolute_error << '\t'
		<< relative_error << '\t'
		<< status << "\t\n";
}

void RunFiniteDifferenceValidation(
	gpmesh::GPSplineSolidV0& spline,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points,
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<MMCComponent>& base_components,
	const std::vector<double>& physical_geometry,
	const std::vector<double>& physical_lower,
	const std::vector<double>& physical_upper,
	const std::vector<double>& base_material_eta,
	const std::vector<Eigen::MatrixXd>& unit_local_matrix,
	const MMCMultiMaterialSimultaneousParams& params,
	int accepted_iteration,
	const SimultaneousEvaluation& base_evaluation,
	std::ofstream& output) {
	const double relative_step = params.finite_difference_step;

	const std::vector<int> mmc_indices = SelectPhysicalMmcFiniteDifferenceIndices(
		physical_geometry,
		base_evaluation.mmc_objective_gradient,
		physical_lower,
		physical_upper,
		relative_step,
		params.finite_difference_mmc_sample_count);
	const std::vector<int> material_indices = SelectFiniteDifferenceIndices(
		base_material_eta,
		base_evaluation.material_objective_gradient,
		params.eta_min,
		1.0,
		relative_step,
		params.finite_difference_material_sample_count);

	std::cout
		<< "[MMC-MultiMaterial-Simultaneous][FiniteDifference] iteration: "
		<< accepted_iteration
		<< " MMC samples: " << mmc_indices.size()
		<< " material samples: " << material_indices.size()
		<< " relative step: " << relative_step << std::endl;

	for (int variable : mmc_indices) {
		const double step = relative_step
			* (physical_upper[variable] - physical_lower[variable]);
		std::vector<double> plus_geometry = physical_geometry;
		std::vector<double> minus_geometry = physical_geometry;
		plus_geometry[variable] += step;
		minus_geometry[variable] -= step;
		std::vector<MMCComponent> plus_components = base_components;
		std::vector<MMCComponent> minus_components = base_components;
		UpdateComponentsFromPhysical(plus_components, plus_geometry);
		UpdateComponentsFromPhysical(minus_components, minus_geometry);
		const SimultaneousEvaluation plus_evaluation = EvaluateSimultaneousState(
			spline,
			neumann_points,
			dirichlet_points,
			all_vertex,
			cell_to_vertex,
			plus_components,
			base_material_eta,
			unit_local_matrix,
			params);
		const SimultaneousEvaluation minus_evaluation = EvaluateSimultaneousState(
			spline,
			neumann_points,
			dirichlet_points,
			all_vertex,
			cell_to_vertex,
			minus_components,
			base_material_eta,
			unit_local_matrix,
			params);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			step,
			"compliance",
			base_evaluation.compliance,
			plus_evaluation.compliance,
			minus_evaluation.compliance,
			base_evaluation.mmc_objective_gradient[variable]);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			step,
			"mmc_volume",
			base_evaluation.mmc_volume,
			plus_evaluation.mmc_volume,
			minus_evaluation.mmc_volume,
			base_evaluation.geometry_constraint_gradient[variable]);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			step,
			"ratio_residual",
			base_evaluation.ratio_residual,
			plus_evaluation.ratio_residual,
			minus_evaluation.ratio_residual,
			base_evaluation.ratio_constraint_mmc_gradient[variable]);
	}

	for (int variable : material_indices) {
		const double step = relative_step;
		std::vector<double> plus_material = base_material_eta;
		std::vector<double> minus_material = base_material_eta;
		plus_material[variable] += step;
		minus_material[variable] -= step;
		const SimultaneousEvaluation plus_evaluation = EvaluateSimultaneousState(
			spline,
			neumann_points,
			dirichlet_points,
			all_vertex,
			cell_to_vertex,
			base_components,
			plus_material,
			unit_local_matrix,
			params);
		const SimultaneousEvaluation minus_evaluation = EvaluateSimultaneousState(
			spline,
			neumann_points,
			dirichlet_points,
			all_vertex,
			cell_to_vertex,
			base_components,
			minus_material,
			unit_local_matrix,
			params);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"material-eta",
			variable,
			step,
			"compliance",
			base_evaluation.compliance,
			plus_evaluation.compliance,
			minus_evaluation.compliance,
			base_evaluation.material_objective_gradient[variable]);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"material-eta",
			variable,
			step,
			"mmc_volume",
			base_evaluation.mmc_volume,
			plus_evaluation.mmc_volume,
			minus_evaluation.mmc_volume,
			0.0);
		WriteFiniteDifferenceQuantity(
			output,
			accepted_iteration,
			"material-eta",
			variable,
			step,
			"ratio_residual",
			base_evaluation.ratio_residual,
			plus_evaluation.ratio_residual,
			minus_evaluation.ratio_residual,
			base_evaluation.ratio_constraint_material_gradient[variable]);
	}
	output.flush();
}

double WeightedOccupancySum(
	const std::vector<double>& weights,
	const std::vector<double>& occupancy) {
	if (weights.size() != occupancy.size()) {
		throw std::runtime_error(
			"The full sensitivity validation received mismatched element arrays.");
	}
	double value = 0.0;
	for (size_t cell_id = 0; cell_id < weights.size(); ++cell_id) {
		value += weights[cell_id] * occupancy[cell_id];
	}
	return value;
}

std::vector<double> EvaluateOccupancyForPerturbedComponent(
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<MMCComponent>& components,
	int perturbed_component,
	const GeometryEvaluation& base_geometry,
	const std::vector<double>& base_exponential_sum,
	const MMCMultiMaterialSimultaneousParams& params) {
	std::vector<double> perturbed_phi;
	std::vector<std::vector<double>> unused_phi_derivative;
	calc_Phi_Node(
		all_vertex,
		0.0,
		components,
		perturbed_phi,
		unused_phi_derivative,
		perturbed_component,
		kTdfPower,
		params.heaviside_epsilon);

	const int vertex_count = static_cast<int>(all_vertex.size());
	std::vector<double> perturbed_heaviside(vertex_count, 0.0);
	for (int vertex_id = 0; vertex_id < vertex_count; ++vertex_id) {
		const double exponential_sum =
			base_exponential_sum[vertex_id]
			- std::exp(
				params.ks_lambda
				* base_geometry.allPhis[perturbed_component][vertex_id])
			+ std::exp(params.ks_lambda * perturbed_phi[vertex_id]);
		const double global_phi = std::max(
			-1e3,
			std::log(std::max(exponential_sum, kTiny)) / params.ks_lambda);
		perturbed_heaviside[vertex_id] = Heaviside(
			global_phi,
			params.heaviside_alpha,
			params.heaviside_epsilon);
	}

	std::vector<double> occupancy(cell_to_vertex.size(), 0.0);
	for (int cell_id = 0;
		cell_id < static_cast<int>(cell_to_vertex.size());
		++cell_id) {
		for (const auto& vertex : cell_to_vertex.at(cell_id)) {
			occupancy[cell_id] += perturbed_heaviside[vertex.first] / 8.0;
		}
	}
	if (!IsFiniteVector(occupancy)) {
		throw std::runtime_error(
			"A perturbed MMC occupancy used by full sensitivity validation is non-finite.");
	}
	return occupancy;
}

void WriteAllFiniteDifferenceValidationHeader(std::ofstream& output) {
	output
		<< "accepted_iteration\tvalidation_kind\tvariable_block"
		<< "\tvariable_index\tcomponent_index\tlocal_variable_index"
		<< "\tvariable_name\tdifference_scheme\tstep\tquantity"
		<< "\tanalytical_derivative\tfinite_difference_derivative"
		<< "\tabsolute_error\trelative_error\tstatus\tmessage\n";
}

void WriteAllFiniteDifferenceValidationFailure(
	std::ofstream& output,
	int accepted_iteration,
	double step,
	const std::string& message) {
	std::string sanitized_message = message;
	std::replace(sanitized_message.begin(), sanitized_message.end(), '\t', ' ');
	std::replace(sanitized_message.begin(), sanitized_message.end(), '\n', ' ');
	const double not_a_number = std::numeric_limits<double>::quiet_NaN();
	output << accepted_iteration
		<< "\tend-to-end-finite-element\tdiagnostic\t-1\t-1\t-1\tnone\tnone\t"
		<< step << "\terror\t"
		<< not_a_number << '\t' << not_a_number << '\t'
		<< not_a_number << '\t' << not_a_number
		<< "\tERROR\t" << sanitized_message << '\n';
	output.flush();
}

void WriteAllFiniteDifferenceValidationRow(
	std::ofstream& output,
	int accepted_iteration,
	const char* variable_block,
	int variable_index,
	int component_index,
	int local_variable_index,
	const char* variable_name,
	const char* difference_scheme,
	double step,
	const char* quantity,
	double analytical_derivative,
	double finite_difference_derivative) {
	const double absolute_error = std::abs(
		analytical_derivative - finite_difference_derivative);
	const double relative_error = absolute_error / std::max({
		std::abs(analytical_derivative),
		std::abs(finite_difference_derivative),
		1e-12
	});
	const char* status = ClassifyFiniteDifferenceError(
		analytical_derivative,
		finite_difference_derivative,
		absolute_error,
		relative_error);
	output << accepted_iteration
		<< "\tend-to-end-finite-element\t" << variable_block
		<< '\t' << variable_index
		<< '\t' << component_index
		<< '\t' << local_variable_index
		<< '\t' << variable_name
		<< '\t' << difference_scheme
		<< '\t' << step
		<< '\t' << quantity
		<< '\t' << analytical_derivative
		<< '\t' << finite_difference_derivative
		<< '\t' << absolute_error
		<< '\t' << relative_error
		<< '\t' << status << "\t\n";
}

void RunAllFiniteDifferenceValidation(
	gpmesh::GPSplineSolidV0& spline,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points,
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<MMCComponent>& base_components,
	const std::vector<double>& physical_geometry,
	const std::vector<double>& physical_lower,
	const std::vector<double>& physical_upper,
	const std::vector<double>& base_material_eta,
	const std::vector<Eigen::MatrixXd>& unit_local_matrix,
	const MMCMultiMaterialSimultaneousParams& params,
	int accepted_iteration,
	const SimultaneousEvaluation& base_evaluation,
	std::ofstream& output) {
	// 严格的逐变量端到端验证，不固定 U 或 q_e：
	//   x_i +/- h -> 重新计算几何/材料 -> 重新装配 K -> 重新求解 U
	//              -> 重新计算 C、MMC体积和比例残差。
	// 内部变量用中心差分；落在上下界的变量用单边差分，保证变量不漏检。
	// 该过程只读取基础设计状态，不推进 MMA 或 OC 历史。
	const int cell_count = static_cast<int>(cell_to_vertex.size());
	const double step = params.finite_difference_step;
	const char* mmc_variable_names[kVariablesPerComponent] = {
		"x0", "y0", "z0", "l1", "l2", "l3",
		"alpha", "beta", "gamma"
	};

	if (static_cast<int>(base_material_eta.size()) != cell_count ||
		static_cast<int>(unit_local_matrix.size()) != cell_count) {
		throw std::runtime_error(
			"The all-variable finite-difference validation received mismatched cell arrays.");
	}
	std::cout
		<< "[MMC-MultiMaterial-Simultaneous][FiniteDifference][All] "
		<< "full finite-element re-solves enabled; MMC variables: "
		<< physical_geometry.size()
		<< " material variables: " << base_material_eta.size()
		<< " expected maximum FE solves: "
		<< 2 * (physical_geometry.size() + base_material_eta.size())
		<< std::endl;

	// 第一块：逐个扰动物理 MMC 变量；步长按对应物理范围缩放。
	for (int variable = 0;
		variable < static_cast<int>(physical_geometry.size());
		++variable) {
		const double physical_range =
			physical_upper[variable] - physical_lower[variable];
		const double mmc_step = step * physical_range;
		const bool can_step_forward =
			physical_geometry[variable] + mmc_step <= physical_upper[variable];
		const bool can_step_backward =
			physical_geometry[variable] - mmc_step >= physical_lower[variable];
		const char* scheme = can_step_forward && can_step_backward
			? "central"
			: (can_step_forward ? "forward" : "backward");
		const int component_index = variable / kVariablesPerComponent;
		const int local_variable_index = variable % kVariablesPerComponent;
		SimultaneousEvaluation plus_evaluation;
		SimultaneousEvaluation minus_evaluation;

		if (can_step_forward) {
			std::vector<double> plus_geometry = physical_geometry;
			plus_geometry[variable] += mmc_step;
			std::vector<MMCComponent> plus_components = base_components;
			UpdateComponentsFromPhysical(plus_components, plus_geometry);
			plus_evaluation = EvaluateSimultaneousState(
				spline,
				neumann_points,
				dirichlet_points,
				all_vertex,
				cell_to_vertex,
				plus_components,
				base_material_eta,
				unit_local_matrix,
				params);
		}
		else {
			plus_evaluation = base_evaluation;
		}
		if (can_step_backward) {
			std::vector<double> minus_geometry = physical_geometry;
			minus_geometry[variable] -= mmc_step;
			std::vector<MMCComponent> minus_components = base_components;
			UpdateComponentsFromPhysical(minus_components, minus_geometry);
			minus_evaluation = EvaluateSimultaneousState(
				spline,
				neumann_points,
				dirichlet_points,
				all_vertex,
				cell_to_vertex,
				minus_components,
				base_material_eta,
				unit_local_matrix,
				params);
		}
		else {
			minus_evaluation = base_evaluation;
		}
		const double denominator =
			can_step_forward && can_step_backward ? 2.0 * mmc_step : mmc_step;
		const double objective_fd =
			(plus_evaluation.compliance - minus_evaluation.compliance)
			/ denominator;
		const double volume_fd =
			(plus_evaluation.mmc_volume - minus_evaluation.mmc_volume)
			/ denominator;
		const double ratio_fd =
			(plus_evaluation.ratio_residual - minus_evaluation.ratio_residual)
			/ denominator;
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			component_index,
			local_variable_index,
			mmc_variable_names[local_variable_index],
			scheme,
			mmc_step,
			"compliance",
			base_evaluation.mmc_objective_gradient[variable],
			objective_fd);
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			component_index,
			local_variable_index,
			mmc_variable_names[local_variable_index],
			scheme,
			mmc_step,
			"mmc_volume",
			base_evaluation.geometry_constraint_gradient[variable],
			volume_fd);
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"MMC-physical",
			variable,
			component_index,
			local_variable_index,
			mmc_variable_names[local_variable_index],
			scheme,
			mmc_step,
			"ratio_residual",
			base_evaluation.ratio_constraint_mmc_gradient[variable],
			ratio_fd);
		std::cout
			<< "[MMC-MultiMaterial-Simultaneous][FiniteDifference][All] MMC "
			<< variable + 1 << '/' << physical_geometry.size()
			<< " variable: " << mmc_variable_names[local_variable_index]
			<< " scheme: " << scheme << std::endl;
	}

	// 第二块：逐个扰动完整网格上的 eta_e；每个 eta_e 都触发独立 FE 重求解。
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		const double eta = base_material_eta[cell_id];
		const bool can_step_forward = eta + step <= 1.0;
		const bool can_step_backward = eta - step >= params.eta_min;
		const char* scheme = can_step_forward && can_step_backward
			? "central"
			: (can_step_forward ? "forward" : "backward");
		const double eta_plus = can_step_forward ? eta + step : eta;
		const double eta_minus = can_step_backward ? eta - step : eta;
		const double denominator =
			can_step_forward && can_step_backward ? 2.0 * step : step;
		SimultaneousEvaluation plus_evaluation;
		SimultaneousEvaluation minus_evaluation;
		if (can_step_forward) {
			std::vector<double> plus_material = base_material_eta;
			plus_material[cell_id] = eta_plus;
			plus_evaluation = EvaluateSimultaneousState(
				spline,
				neumann_points,
				dirichlet_points,
				all_vertex,
				cell_to_vertex,
				base_components,
				plus_material,
				unit_local_matrix,
				params);
		}
		else {
			plus_evaluation = base_evaluation;
		}
		if (can_step_backward) {
			std::vector<double> minus_material = base_material_eta;
			minus_material[cell_id] = eta_minus;
			minus_evaluation = EvaluateSimultaneousState(
				spline,
				neumann_points,
				dirichlet_points,
				all_vertex,
				cell_to_vertex,
				base_components,
				minus_material,
				unit_local_matrix,
				params);
		}
		else {
			minus_evaluation = base_evaluation;
		}
		const double objective_fd =
			(plus_evaluation.compliance - minus_evaluation.compliance)
			/ denominator;
		const double volume_fd =
			(plus_evaluation.mmc_volume - minus_evaluation.mmc_volume)
			/ denominator;
		const double ratio_fd =
			(plus_evaluation.ratio_residual - minus_evaluation.ratio_residual)
			/ denominator;
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"material-eta",
			cell_id,
			-1,
			-1,
			"eta",
			scheme,
			step,
			"compliance",
			base_evaluation.material_objective_gradient[cell_id],
			objective_fd);
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"material-eta",
			cell_id,
			-1,
			-1,
			"eta",
			scheme,
			step,
			"mmc_volume",
			0.0,
			volume_fd);
		WriteAllFiniteDifferenceValidationRow(
			output,
			accepted_iteration,
			"material-eta",
			cell_id,
			-1,
			-1,
			"eta",
			scheme,
			step,
			"ratio_residual",
			base_evaluation.ratio_constraint_material_gradient[cell_id],
			ratio_fd);
		std::cout
			<< "[MMC-MultiMaterial-Simultaneous][FiniteDifference][All] material "
			<< cell_id + 1 << '/' << cell_count
			<< " scheme: " << scheme << std::endl;
	}
	output.flush();
	std::cout
		<< "[MMC-MultiMaterial-Simultaneous][FiniteDifference][All] full finite-element rows: "
		<< physical_geometry.size() * 3 + base_material_eta.size() * 3
		<< std::endl;
}

std::vector<int> BuildSampleIndices(int value_count, int requested_count) {
	std::vector<int> indices;
	if (value_count <= 0 || requested_count <= 0) return indices;
	const int sample_count = std::min(value_count, requested_count);
	if (sample_count == 1) return { value_count / 2 };
	for (int sample = 0; sample < sample_count; ++sample) {
		const int index = static_cast<int>(std::llround(
			static_cast<double>(sample)
			* static_cast<double>(value_count - 1)
			/ static_cast<double>(sample_count - 1)));
		if (indices.empty() || indices.back() != index) indices.push_back(index);
	}
	return indices;
}

bool WriteLegacyHexVtk(
	const std::string& file_path,
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<double>& mmc_occupancy,
	const std::vector<double>& projected_occupancy,
	const std::vector<double>& material_eta,
	const std::vector<double>& material_modulus,
	const std::vector<double>& effective_stiffness,
	double material_binary_active_threshold) {
	const int cell_count = static_cast<int>(cell_to_vertex.size());
	if (all_vertex.empty() || cell_count == 0 ||
		static_cast<int>(mmc_occupancy.size()) != cell_count ||
		static_cast<int>(projected_occupancy.size()) != cell_count ||
		static_cast<int>(material_eta.size()) != cell_count ||
		static_cast<int>(material_modulus.size()) != cell_count ||
		static_cast<int>(effective_stiffness.size()) != cell_count) {
		return false;
	}

	std::map<int, int> vertex_index;
	std::vector<Eigen::Vector3d> ordered_vertices;
	ordered_vertices.reserve(all_vertex.size());
	for (const auto& vertex : all_vertex) {
		vertex_index[vertex.first] = static_cast<int>(ordered_vertices.size());
		ordered_vertices.push_back(vertex.second);
	}

	std::ofstream output(file_path, std::ios::out | std::ios::trunc);
	if (!output.good()) return false;
	output << std::setprecision(17);
	output << "# vtk DataFile Version 2.0\n";
	output << "MMC multi-material simultaneous optimization\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";
	output << "POINTS " << ordered_vertices.size() << " double\n";
	for (const Eigen::Vector3d& vertex : ordered_vertices) {
		output << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
	}

	output << "CELLS " << cell_count << ' ' << cell_count * 9 << "\n";
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		output << 8;
		for (const auto& vertex : cell_to_vertex.at(cell_id)) {
			auto vertex_it = vertex_index.find(vertex.first);
			if (vertex_it == vertex_index.end()) return false;
			output << ' ' << vertex_it->second;
		}
		output << '\n';
	}
	output << "CELL_TYPES " << cell_count << "\n";
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		(void)cell_id;
		output << 12 << '\n';
	}
	output << "CELL_DATA " << cell_count << "\n";

	const auto write_scalar = [&output, cell_count](
		const char* name,
		const std::vector<double>& values) {
		output << "SCALARS " << name << " double 1\n";
		output << "LOOKUP_TABLE default\n";
		for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
			output << values[cell_id] << '\n';
		}
	};
	write_scalar("mmc_occupancy", mmc_occupancy);
	write_scalar("projected_occupancy", projected_occupancy);
	write_scalar("material_eta", material_eta);
	write_scalar("material_modulus", material_modulus);
	write_scalar("effective_stiffness", effective_stiffness);

	// 材料二值只在 m_e 大于传入覆盖域阈值的 MMC 区域内
	// 组件内部有物理含义：-1 表示组件外，0 表示 E2，1 表示 E1。
	// 该字段保留三类数值，避免将组件外潜伏的 eta 错看成真实材料。
	constexpr double material_binary_threshold = 0.5;
	output << "SCALARS material_binary_region int 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		if (mmc_occupancy[cell_id] <= material_binary_active_threshold) {
			output << -1 << '\n';
		}
		else {
			output << (material_eta[cell_id] >= material_binary_threshold ? 1 : 0)
				<< '\n';
		}
	}

	// ParaView 可直接选择该 RGB 字段显示：组件外为浅灰，组件内
	// E2（eta < 0.5）为蓝色，E1（eta >= 0.5）为红色。
	output << "COLOR_SCALARS material_binary_color 3\n";
	for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
		if (mmc_occupancy[cell_id] <= material_binary_active_threshold) {
			output << "0.85 0.85 0.85\n";
		}
		else if (material_eta[cell_id] >= material_binary_threshold) {
			output << "1.0 0.2 0.15\n";
		}
		else {
			output << "0.15 0.35 1.0\n";
		}
	}
	return output.good();
}

bool WriteActiveMaterialBinaryVtkLikeSequential(
	const std::string& file_path,
	const std::map<int, Eigen::Vector3d>& all_vertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell_to_vertex,
	const std::vector<double>& mmc_occupancy,
	const std::vector<double>& material_eta,
	double material_binary_active_threshold) {
	// 最终材料 VTK：先用传入的 MMC 覆盖域阈值构造
	// active_indices，输出网格中只保留这些组件内单元；然后仅在这些
	// 单元上按 eta>=0.5 二值。组件外单元不会混入材料 0 类。
	if (all_vertex.empty() || cell_to_vertex.empty() ||
		mmc_occupancy.size() != cell_to_vertex.size() ||
		material_eta.size() != cell_to_vertex.size()) {
		return false;
	}
	std::vector<int> active_indices;
	active_indices.reserve(cell_to_vertex.size());
	for (int cell_id = 0;
		cell_id < static_cast<int>(cell_to_vertex.size());
		++cell_id) {
		if (mmc_occupancy[cell_id] > material_binary_active_threshold) {
			active_indices.push_back(cell_id);
		}
	}
	if (active_indices.empty()) return false;

	std::map<int, int> vertex_index;
	std::vector<Eigen::Vector3d> ordered_vertices;
	ordered_vertices.reserve(all_vertex.size());
	for (const auto& vertex : all_vertex) {
		vertex_index[vertex.first] = static_cast<int>(ordered_vertices.size());
		ordered_vertices.push_back(vertex.second);
	}

	std::ofstream output(file_path, std::ios::out | std::ios::trunc);
	if (!output.good()) return false;
	output << std::setprecision(17);
	output << "# vtk DataFile Version 2.0\n";
	output << "MMC multi-material simultaneous active material result\n";
	output << "ASCII\n";
	output << "DATASET UNSTRUCTURED_GRID\n";
	output << "POINTS " << ordered_vertices.size() << " double\n";
	for (const Eigen::Vector3d& vertex : ordered_vertices) {
		output << vertex.x() << ' ' << vertex.y() << ' ' << vertex.z() << '\n';
	}

	output << "CELLS " << active_indices.size() << ' '
		<< active_indices.size() * 9 << "\n";
	for (int cell_id : active_indices) {
		auto cell_it = cell_to_vertex.find(cell_id);
		if (cell_it == cell_to_vertex.end() || cell_it->second.size() != 8) {
			return false;
		}
		output << 8;
		for (const auto& vertex : cell_it->second) {
			auto vertex_it = vertex_index.find(vertex.first);
			if (vertex_it == vertex_index.end()) return false;
			output << ' ' << vertex_it->second;
		}
		output << '\n';
	}
	output << "CELL_TYPES " << active_indices.size() << "\n";
	for (size_t index = 0; index < active_indices.size(); ++index) {
		(void)index;
		output << 12 << '\n';
	}

	output << "CELL_DATA " << active_indices.size() << "\n";
	output << "SCALARS material_eta double 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int cell_id : active_indices) {
		output << material_eta[cell_id] << '\n';
	}

	// 顺序版同名二值语义：0 为 E2，1 为 E1。
	output << "SCALARS material_binary int 1\n";
	output << "LOOKUP_TABLE default\n";
	for (int cell_id : active_indices) {
		const double eta = std::clamp(material_eta[cell_id], 0.0, 1.0);
		output << (eta >= 0.5 ? 1 : 0) << '\n';
	}

	output << "COLOR_SCALARS material_color 3\n";
	for (int cell_id : active_indices) {
		const double eta = std::clamp(material_eta[cell_id], 0.0, 1.0);
		output << eta << " 0 " << (1.0 - eta) << '\n';
	}

	output << "COLOR_SCALARS material_binary_color 3\n";
	for (int cell_id : active_indices) {
		const double eta = std::clamp(material_eta[cell_id], 0.0, 1.0);
		const double material_one = eta >= 0.5 ? 1.0 : 0.0;
		output << material_one << " 0 " << (1.0 - material_one) << '\n';
	}
	return output.good();
}

bool WriteFinalFieldTable(
	const std::string& file_path,
	const std::vector<double>& mmc_occupancy,
	const std::vector<double>& projected_occupancy,
	const std::vector<double>& material_eta,
	const std::vector<double>& material_modulus,
	const std::vector<double>& effective_stiffness) {
	const size_t field_count = mmc_occupancy.size();
	if (projected_occupancy.size() != field_count ||
		material_eta.size() != field_count ||
		material_modulus.size() != field_count ||
		effective_stiffness.size() != field_count) {
		return false;
	}
	std::ofstream output(file_path, std::ios::out | std::ios::trunc);
	if (!output.good()) return false;
	output << std::setprecision(17);
	output << "cell_id\tmmc_occupancy\tprojected_occupancy\tmaterial_eta\tmaterial_modulus"
		<< "\teffective_stiffness\n";
	for (size_t cell_id = 0; cell_id < field_count; ++cell_id) {
		output << cell_id << '\t'
			<< mmc_occupancy[cell_id] << '\t'
			<< projected_occupancy[cell_id] << '\t'
			<< material_eta[cell_id] << '\t'
			<< material_modulus[cell_id] << '\t'
			<< effective_stiffness[cell_id] << '\n';
	}
	return output.good();
}

void WriteIterationHeader(std::ofstream& output) {
	output
		<< "iteration\tcompliance\tcompliance_from_load"
		<< "\tcompliance_relative_mismatch\tprojection_beta"
		<< "\tmmc_occupancy_min\tmmc_occupancy_max"
		<< "\tprojected_occupancy_min\tprojected_occupancy_max"
		<< "\tmmc_volume\tV1\tV2"
		<< "\tgeometry_constraint"
		<< "\tgeometry_max_change\tmaterial_max_change"
		<< "\tgeometry_gradient_norm\tmaterial_gradient_norm"
		<< "\tuniform_material_gradient_relative_error\n";
}

void WriteIterationRecord(
	std::ofstream& output,
	const MMCMultiMaterialIterationDiagnostic& diagnostic) {
	output << diagnostic.iteration << '\t'
		<< diagnostic.compliance << '\t'
		<< diagnostic.compliance_from_load << '\t'
		<< diagnostic.compliance_relative_mismatch << '\t'
		<< diagnostic.projection_beta << '\t'
		<< diagnostic.mmc_occupancy_min << '\t'
		<< diagnostic.mmc_occupancy_max << '\t'
		<< diagnostic.projected_occupancy_min << '\t'
		<< diagnostic.projected_occupancy_max << '\t'
		<< diagnostic.mmc_volume << '\t'
		<< diagnostic.V1 << '\t'
		<< diagnostic.V2 << '\t'
		<< diagnostic.geometry_constraint << '\t'
		<< diagnostic.geometry_max_change << '\t'
		<< diagnostic.material_max_change << '\t'
		<< diagnostic.geometry_gradient_norm << '\t'
		<< diagnostic.material_gradient_norm << '\t'
		<< diagnostic.uniform_material_gradient_relative_error << '\n';
	output.flush();
}

void WriteSensitivityHeader(std::ofstream& output) {
	output
		<< "accepted_iteration\tvariable_block\tvariable_index"
		<< "\tanalytical_value\tmax_abs_contribution\tsum_abs_contribution\n";
}

void WriteSensitivityDiagnostics(
	std::ofstream& output,
	int accepted_iteration,
	const SimultaneousEvaluation& evaluation,
	int requested_samples) {
	const std::vector<int> mmc_samples = BuildSampleIndices(
		static_cast<int>(evaluation.mmc_objective_gradient.size()),
		requested_samples);
	for (int variable : mmc_samples) {
		output << accepted_iteration << "\tMMC\t" << variable << '\t'
			<< evaluation.mmc_objective_gradient[variable] << '\t'
			<< evaluation.mmc_contribution_max_abs[variable] << '\t'
			<< evaluation.mmc_contribution_sum_abs[variable] << '\n';
	}

	const std::vector<int> material_samples = BuildSampleIndices(
		static_cast<int>(evaluation.material_objective_gradient.size()),
		requested_samples);
	for (int variable : material_samples) {
		const double absolute_gradient =
			std::abs(evaluation.material_objective_gradient[variable]);
		output << accepted_iteration << "\tmaterial\t" << variable << '\t'
			<< evaluation.material_objective_gradient[variable] << '\t'
			<< absolute_gradient << '\t'
			<< absolute_gradient << '\n';
	}
	output.flush();
}

void CopyEvaluationToResultState(
	const SimultaneousEvaluation& evaluation,
	const std::vector<double>& physical_geometry,
	const std::vector<double>& material_eta,
	int accepted_iteration,
	MMCMultiMaterialSimultaneousState& state) {
	state.accepted_iteration = accepted_iteration;
	state.physical_mmc_variables = physical_geometry;
	state.material_eta = material_eta;
	state.mmc_occupancy = evaluation.geometry.density;
	state.projected_occupancy = evaluation.projected_occupancy;
	state.material_modulus = evaluation.material_modulus;
	state.effective_stiffness = evaluation.effective_stiffness;
	state.element_energy = evaluation.element_energy;
	state.mmc_objective_gradient = evaluation.mmc_objective_gradient;
	state.material_objective_gradient = evaluation.material_objective_gradient;
	state.geometry_constraint_gradient =
		evaluation.geometry_constraint_gradient;
	state.ratio_constraint_mmc_gradient =
		evaluation.ratio_constraint_mmc_gradient;
	state.ratio_constraint_material_gradient =
		evaluation.ratio_constraint_material_gradient;
}

MMCMultiMaterialIterationDiagnostic BuildDiagnostic(
	int iteration,
	const SimultaneousEvaluation& evaluation,
	double geometry_change,
	double material_change) {
	MMCMultiMaterialIterationDiagnostic diagnostic;
	diagnostic.iteration = iteration;
	diagnostic.compliance = evaluation.compliance;
	diagnostic.compliance_from_load = evaluation.compliance_from_load;
	diagnostic.compliance_relative_mismatch =
		evaluation.compliance_relative_mismatch;
	diagnostic.projection_beta = evaluation.projection_beta;
	if (!evaluation.geometry.density.empty()) {
		const auto [minimum, maximum] = std::minmax_element(
			evaluation.geometry.density.begin(),
			evaluation.geometry.density.end());
		diagnostic.mmc_occupancy_min = *minimum;
		diagnostic.mmc_occupancy_max = *maximum;
	}
	if (!evaluation.projected_occupancy.empty()) {
		const auto [minimum, maximum] = std::minmax_element(
			evaluation.projected_occupancy.begin(),
			evaluation.projected_occupancy.end());
		diagnostic.projected_occupancy_min = *minimum;
		diagnostic.projected_occupancy_max = *maximum;
	}
	diagnostic.mmc_volume = evaluation.mmc_volume;
	diagnostic.V1 = evaluation.V1;
	diagnostic.V2 = evaluation.V2;
	diagnostic.geometry_constraint = evaluation.geometry_constraint;
	diagnostic.geometry_max_change = geometry_change;
	diagnostic.material_max_change = material_change;
	diagnostic.geometry_gradient_norm =
		L2Norm(evaluation.mmc_objective_gradient);
	diagnostic.material_gradient_norm =
		L2Norm(evaluation.material_objective_gradient);
	diagnostic.uniform_material_gradient_relative_error =
		evaluation.uniform_material_gradient_relative_error;
	return diagnostic;
}

bool IsFiniteDiagnostic(const MMCMultiMaterialIterationDiagnostic& diagnostic) {
	const std::array<double, 20> fields = {
		diagnostic.compliance,
		diagnostic.compliance_from_load,
		diagnostic.compliance_relative_mismatch,
		diagnostic.projection_beta,
		diagnostic.mmc_occupancy_min,
		diagnostic.mmc_occupancy_max,
		diagnostic.projected_occupancy_min,
		diagnostic.projected_occupancy_max,
		diagnostic.mmc_volume,
		diagnostic.V1,
		diagnostic.V2,
		diagnostic.geometry_constraint,
		diagnostic.geometry_max_change,
		diagnostic.material_max_change,
		diagnostic.geometry_gradient_norm,
		diagnostic.material_gradient_norm,
		diagnostic.uniform_material_gradient_relative_error
	};
	return std::all_of(fields.begin(), fields.end(), [](double value) {
		return std::isfinite(value);
	});
}

} // namespace mmc_multimaterial_simultaneous

bool Cuda_TopOpt_MMC_MultiMaterial_Simultaneous(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	const MMCMultiMaterialSimultaneousParams& params,
	MMCMultiMaterialSimultaneousResult& result,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model) {
	using namespace mmc_multimaterial_simultaneous;

	// -------------------- 0. 输入检查与独立状态初始化 --------------------
	// result 在入口清零；任何异常都只写入本同步流程的 result，不会进入
	// 独立 MMC、SIMP 或顺序 MMC+SIMP 的状态和优化器。
	result = MMCMultiMaterialSimultaneousResult{};
	std::string validation_error;
	if (!ValidateInputs(
		allVertex,
		cell2V,
		input_model_file,
		bdy_file,
		protect_cell_ids,
		neighbor_info,
		params,
		density,
		validation_error)) {
		SetFailure(result, validation_error);
		return false;
	}

	try {
		RunClosedFormSelfChecks(params);
		Timer timer;
		Eigen::setNbThreads(24);
		Eigen::initParallel();

		// -------------------- 1. IGA 模型与边界条件 --------------------
		std::vector<CBSplineSolid> patch;
		ToolFunc::readFileXml(input_model_file, patch);
		if (patch.empty() ||
			static_cast<int>(patch.size()) != static_cast<int>(cell2V.size())) {
			throw std::runtime_error(
				"The Bezier patch count does not match the fixed-grid cell count.");
		}
		timer.out("simultaneous_read_xml");

		std::vector<CPoint> boundary_points;
		std::vector<int> boundary_types;
		std::vector<CVector> boundary_values;
		ReadBoundaryFile(
			bdy_file,
			boundary_points,
			boundary_types,
			boundary_values);

		gpmesh::GPLog::init(spdlog::level::trace);
		CUDA_ERROR(cudaDeviceReset());

		gpmesh::GPSplineSolidV0 spline;
		spline.init_patch(patch);
		std::vector<boundary_point> neumann_points;
		std::vector<boundary_point> dirichlet_points;
		BuildAnalysisBoundaryConditions(
			patch,
			spline,
			boundary_points,
			boundary_types,
			boundary_values,
			neumann_points,
			dirichlet_points);
		spline.init_dirichlet_points(dirichlet_points);

		// K_e^0 只计算一次。每轮通过 a_e=m_e E_e 缩放并装配，避免把材料
		// 插值或 MMC 占有率重复烘焙进基础单元矩阵。
		std::vector<Eigen::MatrixXd> unit_local_matrix(
			patch.size(),
			Eigen::MatrixXd::Zero(192, 192));
		spline.Get_Local_Matrix(
			unit_local_matrix,
			1.0,
			params.Poissons);
		spline.m_sparse_system.rhs().setZero();
		spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
		spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);
		timer.out("simultaneous_analysis_initialization");

		// -------------------- 2. MMC 物理变量与完整网格材料变量 --------------------
		// mmc_type 只影响初始组件数量/姿态；每个组件仍保持固定 9 变量顺序。
		const std::vector<double> xyz_range = ComputeDomainRange(allVertex);
		std::vector<MMCComponent> mmcComponents =
			InitializeComponentsForType(xyz_range, params.mmc_type);
		if (mmcComponents.empty()) {
			throw std::runtime_error(
				"No initial MMC component was generated for this domain and MMC type.");
		}

		std::vector<double> physical_lower;
		std::vector<double> physical_upper;
		BuildPhysicalComponentBounds(
			xyz_range,
			static_cast<int>(mmcComponents.size()),
			physical_lower,
			physical_upper);
		// 与独立 MMC 保持一致：组件的九个变量直接以物理量进入 MMA，
		// 不在外层再做 [0,1] 归一化，也不重复缩放解析梯度。
		std::vector<double> mmcDesignVariables;
		FlattenPhysicalComponents(mmcComponents, mmcDesignVariables);
		const int geometry_variable_count =
			static_cast<int>(mmcDesignVariables.size());
		const int cell_count = static_cast<int>(cell2V.size());
		double initial_eta = params.initial_eta;
		if (!(initial_eta >= params.eta_min && initial_eta <= 1.0)) {
			initial_eta = 1.0 / (1.0 + params.material_ratio_f);
		}
		initial_eta = std::clamp(initial_eta, params.eta_min, 1.0);
		// 始终为全部单元保留一个 eta_e；材料 OC 只更新当前 MMC 活动单元，
		// 但完整网格映射不重建，避免 MMC 移动时丢失单元材料历史。
		std::vector<double> materialEta(cell_count, initial_eta);
		// 与独立 MMC 一致：几何块只用一个体积约束 MMA；材料块沿用
		// 顺序 MMC-SIMP 的过滤 + OC 二分更新。
		MMASolver mmcMmaSolver(geometry_variable_count, 1);
		std::cout
			<< "[MMC-MultiMaterial-Simultaneous] Update scheme: MMC-MMA + material-OC, "
			<< "mmc_variables=" << geometry_variable_count
			<< ", material_variables=" << cell_count << std::endl;

		std::ofstream iteration_output;
		std::ofstream sensitivity_output;
		std::ofstream finite_difference_output;
		std::ofstream mmc_mma_trace_output;
		std::ofstream material_oc_trace_output;
		std::unique_ptr<mmc_trace::Writer> debug_output;
		if (params.export_results) {
			if (params.output_directory.empty()) {
				throw std::runtime_error("The result output directory cannot be empty.");
			}
			std::filesystem::create_directories(params.output_directory);
			iteration_output.open(
				params.output_directory + "/iterations.tsv",
				std::ios::out | std::ios::trunc);
			sensitivity_output.open(
				params.output_directory + "/analytical_sensitivities.tsv",
				std::ios::out | std::ios::trunc);
			if (params.finite_difference_mmc_sample_count > 0 ||
				params.finite_difference_material_sample_count > 0) {
				finite_difference_output.open(
					params.output_directory + "/finite_difference_validation.tsv",
					std::ios::out | std::ios::trunc);
			}
			mmc_mma_trace_output.open(
				params.output_directory + "/mmc_mma_trace.tsv",
				std::ios::out | std::ios::trunc);
			material_oc_trace_output.open(
				params.output_directory + "/material_oc_trace.tsv",
				std::ios::out | std::ios::trunc);
			if (!iteration_output.good() || !sensitivity_output.good()
				|| !mmc_mma_trace_output.good()
				|| !material_oc_trace_output.good()
				|| (finite_difference_output.is_open()
					&& !finite_difference_output.good())) {
				throw std::runtime_error("Cannot create the simultaneous-optimization diagnostics file.");
			}
			iteration_output << std::setprecision(17);
			sensitivity_output << std::setprecision(17);
			if (finite_difference_output.is_open()) {
				finite_difference_output << std::setprecision(17);
				WriteFiniteDifferenceHeader(finite_difference_output);
			}
			mmc_mma_trace_output << std::setprecision(17);
			material_oc_trace_output << std::setprecision(17);
			WriteIterationHeader(iteration_output);
			WriteSensitivityHeader(sensitivity_output);
			mmc_mma_trace_output
				<< "iteration\tvariable_index"
				<< "\tvalue_before\tvalue_after\tactual_step"
				<< "\tobjective_gradient\tvolume_gradient\n";
			material_oc_trace_output
				<< "iteration\tcell_id\tvalue_before\tvalue_after\tactual_step"
				<< "\tobjective_gradient\toc_sensitivity_filtered"
				<< "\tmmc_occupancy\tmaterial_weight\n";
			debug_output = std::make_unique<mmc_trace::Writer>(
				params.output_directory + "/mmc_debug_trace.tsv");
			if (!debug_output->Good()) {
				throw std::runtime_error(
					"Cannot create the simultaneous MMC debug trace file.");
			}
		}

		double previous_geometry_change = 0.0;
		double previous_material_change = 0.0;
		std::vector<double> previous_geometry_step_physical;
		std::vector<double> previous_material_step;
		double previous_objective_linear_geometry = 0.0;
		double previous_objective_linear_material = 0.0;
		double previous_geometry_constraint_linear_geometry = 0.0;
		int accepted_updates = 0;
		int last_recorded_accepted_iteration = -1;
		SimultaneousEvaluation current_evaluation;
		// 只读记录函数：不会更新设计变量或 MMA，仅输出当前已接受状态。
		const auto record_accepted_state =
			[&](const SimultaneousEvaluation& evaluation) {
			const MMCMultiMaterialIterationDiagnostic diagnostic =
				BuildDiagnostic(
					accepted_updates,
					evaluation,
					previous_geometry_change,
					previous_material_change);
			if (!IsFiniteDiagnostic(diagnostic)) {
				throw std::runtime_error(
					"A simultaneous iteration diagnostic value is non-finite.");
			}
			result.diagnostics.push_back(diagnostic);
			if (params.export_results) {
				WriteIterationRecord(iteration_output, diagnostic);
				WriteSensitivityDiagnostics(
					sensitivity_output,
					accepted_updates,
					evaluation,
					params.sensitivity_sample_count);
			}
			density = evaluation.geometry.density;
			is_need_update_model = true;
			std::cout
				<< "[MMC-MultiMaterial-Simultaneous] Accepted: "
				<< accepted_updates
				<< " Obj.: " << evaluation.compliance
				<< " MMC-Vol.: " << evaluation.mmc_volume
				<< " V1: " << evaluation.V1
				<< " V2: " << evaluation.V2
				<< " beta: " << evaluation.projection_beta
				<< " dMMC: " << previous_geometry_change
				<< " dEta: " << previous_material_change
				<< " |gMMC|: " << diagnostic.geometry_gradient_norm
				<< " |gEta|: " << diagnostic.material_gradient_norm
				<< std::endl;
			last_recorded_accepted_iteration = accepted_updates;
		};

		// -------------------- 3. 同步优化主循环 --------------------
		// 一轮只有一个当前状态。MMC-MMA 和材料-OC 都使用该状态的同一次
		// 位移/单元能量；更新两块变量之间不重新求解有限元。
		for (int iteration = 1;
			iteration <= params.max_iterations;
			++iteration) {
			UpdateComponentsFromPhysical(
				mmcComponents,
				mmcDesignVariables);
			// 在更新任意变量前，用当前 m_e 与 eta_e 重新计算一次共享 FE 状态。
			current_evaluation = EvaluateSimultaneousState(
				spline,
				neumann_points,
				dirichlet_points,
				allVertex,
				cell2V,
				mmcComponents,
				materialEta,
				unit_local_matrix,
				params,
				iteration);
			record_accepted_state(current_evaluation);

			// 这些量只写入诊断日志。与独立 MMC 一致，外层不根据变化量或
			// 约束残差提前结束，而是固定完成 params.max_iterations 次更新。
			const double active_constraint_residual = std::max(
				0.0,
				current_evaluation.geometry_constraint);
			if (debug_output) {
				std::vector<double> physical_variables;
				FlattenPhysicalComponents(mmcComponents, physical_variables);
				mmc_trace::Record trace;
				trace.mode = "simultaneous";
				trace.action = "evaluate";
				trace.iteration = accepted_updates;
				trace.component_count = static_cast<int>(mmcComponents.size());
				trace.vertex_count = static_cast<int>(allVertex.size());
				trace.cell_count = cell_count;
				trace.objective = current_evaluation.compliance;
				trace.volume = current_evaluation.mmc_volume;
				trace.volume_constraint = current_evaluation.geometry_constraint;
				trace.material_v1 = current_evaluation.V1;
				trace.material_v2 = current_evaluation.V2;
				trace.active_constraint_residual = active_constraint_residual;
				trace.geometry_change_max = previous_geometry_change;
				trace.material_change_max = previous_material_change;
				trace.previous_objective_linear_geometry =
					previous_objective_linear_geometry;
				trace.previous_objective_linear_material =
					previous_objective_linear_material;
				trace.previous_geometry_constraint_linear_geometry =
					previous_geometry_constraint_linear_geometry;
				trace.global_phi = &current_evaluation.geometry.globalPhi;
				trace.heaviside = &current_evaluation.geometry.vertexH;
				trace.heaviside_derivative =
					&current_evaluation.geometry.delta_H;
				trace.occupancy = &current_evaluation.geometry.density;
				trace.physical_variables = &physical_variables;
				trace.mmc_objective_gradient =
					&current_evaluation.mmc_objective_gradient;
				trace.geometry_constraint_gradient =
					&current_evaluation.geometry_constraint_gradient;
				trace.material_eta = &materialEta;
				trace.material_modulus = &current_evaluation.material_modulus;
				trace.material_objective_gradient =
					&current_evaluation.material_objective_gradient;
				trace.geometry_step = &previous_geometry_step_physical;
				trace.material_step = &previous_material_step;
				debug_output->Write(trace);
			}
			// -------------------- 4. MMC-MMA + 材料-OC 更新 --------------------
			const std::vector<double> oldMmcDesignVariables = mmcDesignVariables;
			const std::vector<double> oldMaterialEta = materialEta;

			std::vector<double> mmc_constraint_values = {
				current_evaluation.geometry_constraint
			};
			if (!IsFiniteVector(current_evaluation.mmc_objective_gradient) ||
				!IsFiniteVector(current_evaluation.geometry_constraint_gradient) ||
				!IsFiniteVector(mmc_constraint_values)) {
				throw std::runtime_error(
					"An MMC-MMA objective or volume-constraint input is non-finite.");
			}

			// 几何块完全按纯 MMC 模板：物理九变量 + 一个几何体积约束。
			mmcMmaSolver.Update(
				mmcDesignVariables.data(),
				current_evaluation.mmc_objective_gradient.data(),
				mmc_constraint_values.data(),
				current_evaluation.geometry_constraint_gradient.data(),
				physical_lower.data(),
				physical_upper.data());
			if (!IsFiniteVector(mmcDesignVariables)) {
				throw std::runtime_error(
					"The MMC MMA update returned a non-finite design variable.");
			}
			for (int variable = 0; variable < geometry_variable_count; ++variable) {
				mmcDesignVariables[variable] = std::clamp(
					mmcDesignVariables[variable],
					physical_lower[variable],
					physical_upper[variable]);
			}

			std::vector<double> oc_sensitivity_full(
				cell_count,
				std::numeric_limits<double>::min());
			std::vector<double> material_filter_eta(cell_count, 0.0);
			for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
				// eta 越大越偏软材料 E1；在当前 m_e 固定时，真实
				// dC/deta 为正。OC 需要正的下降敏感度量，这里直接使用它。
				oc_sensitivity_full[cell_id] =
					std::max(
						current_evaluation.material_objective_gradient[cell_id],
						std::numeric_limits<double>::min());
				if (current_evaluation.geometry.density[cell_id]
					> params.material_active_threshold) {
					material_filter_eta[cell_id] = materialEta[cell_id];
				}
			}
			const std::vector<double> filtered_oc_sensitivity =
				FilterSensitivities(
					neighbor_info,
					material_filter_eta,
					oc_sensitivity_full,
					params.eta_min);
			OCUpdate(
				materialEta,
				params.material_ratio_f,
				filtered_oc_sensitivity,
				current_evaluation.geometry.density,
				current_evaluation.geometry.density,
				params.material_active_threshold,
				params.eta_min,
				params.material_move_limit);
			if (!IsFiniteVector(materialEta)) {
				throw std::runtime_error(
					"The material OC update returned a non-finite eta value.");
			}
			for (double& eta : materialEta) {
				eta = std::clamp(eta, params.eta_min, 1.0);
			}

			if (params.export_results) {
				for (int variable = 0; variable < geometry_variable_count; ++variable) {
					mmc_mma_trace_output
						<< iteration << '\t'
						<< variable << '\t'
						<< oldMmcDesignVariables[variable] << '\t'
						<< mmcDesignVariables[variable] << '\t'
						<< (mmcDesignVariables[variable]
							- oldMmcDesignVariables[variable]) << '\t'
						<< current_evaluation.mmc_objective_gradient[variable] << '\t'
						<< current_evaluation.geometry_constraint_gradient[variable]
						<< '\n';
				}
				mmc_mma_trace_output.flush();
				for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
					if (current_evaluation.geometry.density[cell_id]
						<= params.material_active_threshold) {
						continue;
					}
					material_oc_trace_output
						<< iteration << '\t'
						<< cell_id << '\t'
						<< oldMaterialEta[cell_id] << '\t'
						<< materialEta[cell_id] << '\t'
						<< (materialEta[cell_id] - oldMaterialEta[cell_id])
						<< '\t'
						<< current_evaluation.material_objective_gradient[cell_id]
						<< '\t'
						<< filtered_oc_sensitivity[cell_id]
						<< '\t'
						<< current_evaluation.geometry.density[cell_id]
						<< '\t'
						<< current_evaluation.geometry.density[cell_id]
						<< '\n';
				}
				material_oc_trace_output.flush();
			}

			// 保存两块实际物理步长，用于下一轮只读诊断。
			previous_geometry_change = 0.0;
			previous_geometry_step_physical.assign(
				geometry_variable_count,
				0.0);
			for (int variable = 0;
				variable < geometry_variable_count;
				++variable) {
				previous_geometry_change = std::max(
					previous_geometry_change,
					std::abs(mmcDesignVariables[variable]
						- oldMmcDesignVariables[variable]));
				previous_geometry_step_physical[variable] =
					mmcDesignVariables[variable]
					- oldMmcDesignVariables[variable];
			}
			previous_material_change = 0.0;
			previous_material_step.assign(cell_count, 0.0);
			for (int cell_id = 0; cell_id < cell_count; ++cell_id) {
				previous_material_change = std::max(
					previous_material_change,
					std::abs(materialEta[cell_id] - oldMaterialEta[cell_id]));
				previous_material_step[cell_id] =
					materialEta[cell_id] - oldMaterialEta[cell_id];
			}
			previous_objective_linear_geometry = mmc_trace::Dot(
				current_evaluation.mmc_objective_gradient,
				previous_geometry_step_physical);
			previous_objective_linear_material = mmc_trace::Dot(
				current_evaluation.material_objective_gradient,
				previous_material_step);
			previous_geometry_constraint_linear_geometry = mmc_trace::Dot(
				current_evaluation.geometry_constraint_gradient,
				previous_geometry_step_physical);
			++accepted_updates;

		}

		// -------------------- 5. 最终状态重算与独立字段导出 --------------------
		// 循环退出后必须用最终设计变量再求解一次，不能把更新前状态当最终结果。
		UpdateComponentsFromPhysical(
			mmcComponents,
			mmcDesignVariables);
		const SimultaneousEvaluation final_evaluation = EvaluateSimultaneousState(
			spline,
			neumann_points,
			dirichlet_points,
			allVertex,
			cell2V,
			mmcComponents,
			materialEta,
			unit_local_matrix,
			params,
			params.max_iterations + 1);
		if (last_recorded_accepted_iteration != accepted_updates) {
			record_accepted_state(final_evaluation);
		}
		CopyEvaluationToResultState(
			final_evaluation,
			mmcDesignVariables,
			materialEta,
			accepted_updates,
			result.state);
		density = final_evaluation.geometry.density;
		is_need_update_model = true;

		if (params.export_results) {
			// VTK/TSV 分开导出 m_e、eta_e、E_e 和 a_e，避免 density 语义混淆。
			const bool vtk_written = WriteLegacyHexVtk(
				params.output_directory
					+ "/mmc_multimaterial_simultaneous.vtk",
				allVertex,
				cell2V,
				final_evaluation.geometry.density,
				final_evaluation.projected_occupancy,
				materialEta,
				final_evaluation.material_modulus,
				final_evaluation.effective_stiffness,
				params.material_binary_active_threshold);
			// 另存一个与顺序 MMC-SIMP 最终二值图相同网格语义的文件：
			// 仅包含 m_e>material_binary_active_threshold 的组件内单元，并在其上
			// 按 eta=0.5 二值。
			const bool binary_vtk_written = WriteActiveMaterialBinaryVtkLikeSequential(
				params.output_directory
					+ "/mmc_multimaterial_simultaneous_material_binary.vtk",
				allVertex,
				cell2V,
				final_evaluation.geometry.density,
				materialEta,
				params.material_binary_active_threshold);
			const bool table_written = WriteFinalFieldTable(
				params.output_directory + "/element_fields.tsv",
				final_evaluation.geometry.density,
				final_evaluation.projected_occupancy,
				materialEta,
				final_evaluation.material_modulus,
				final_evaluation.effective_stiffness);
			if (!vtk_written || !binary_vtk_written || !table_written) {
				throw std::runtime_error("Failed to export the final simultaneous-optimization fields.");
			}
		}

		result.success = true;
		result.message =
			"MMC multi-material simultaneous optimization completed the requested iteration count.";
		return true;
	}
	catch (const std::exception& exception) {
		SetFailure(result, exception.what());
		return false;
	}
	catch (...) {
		SetFailure(result, "An unknown error occurred during simultaneous optimization.");
		return false;
	}
}
