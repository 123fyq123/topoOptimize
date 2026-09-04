#pragma once

#include <atomic>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

// MMC-双材料同步优化的独立公开接口。
//
// 设计变量由两块组成：
//   1. 每个 MMC 组件的 9 个几何变量
//      [x0, y0, z0, l1, l2, l3, alpha, beta, gamma]；
//   2. 固定网格上每个单元一个材料变量 eta_e。
//
// 单元耦合关系：
//   E_e = eta_e^p E1 + (1 - eta_e^p) E2
//   a_e = m_e E_e
//   K_e = a_e K_e^0
// 其中 m_e 是原始 MMC 占有率，K_e^0 是单位
// 杨氏模量下的单元刚度矩阵。
struct MMCMultiMaterialSimultaneousParams {
	// MMC 几何体积与组件初始化参数。
	double mmc_target_rate = 0.3;
	int max_iterations = 100;
	int mmc_type = 1;

	// 双材料插值参数。material_ratio_f 表示目标关系 V2 = f * V1。
	double material_ratio_f = 1.5;
	double E1 = 1.0;
	// 与现有顺序 MMC-SIMP 的实际比较工况保持一致。
	double E2 = 2.0;
	double material_penalty = 3.0;
	// initial_eta < 0 时自动使用 1 / (1 + material_ratio_f)。
	double initial_eta = -1.0;
	double eta_min = 1e-3;
	// 材料 OC 每轮允许的最大移动量。
	double material_move_limit = 0.2;
	// 投影参数仅为兼容旧诊断字段保留；当前主模型直接使用 m_e。
	// 若后续重新启用 ProjectOccupancy()，这些参数仍可复用。
	double projection_threshold = 0.3;
	double projection_beta_initial = 4.0;
	double projection_beta_max = 32.0;
	int projection_beta_update_interval = 20;
	// 材料 OC 只在当前 MMC 活动区域内更新，与顺序版一样按 m_e>0.3 筛选。
	double material_active_threshold = 0.3;
	// 只用于结果显示：与顺序版一致，仅 m_e>0.3 的单元进入材料二值图。
	double material_binary_active_threshold = 0.3;
	// 有限元与 MMC 平滑几何参数。
	double Poissons = 0.3;
	double heaviside_alpha = 1e-3;
	double heaviside_epsilon = 0.20;
	double ks_lambda = 100.0;

	// 同步流程固定运行 max_iterations 次，变量变化量与约束残差只作为
	// 诊断，不参与提前终止。MMC 使用持久化 MMA，材料 eta 使用 OC。
	// ratio 约束已从同步流程移除；该字段仅为旧诊断接口兼容保留。
	double compliance_consistency_tolerance = 5e-3;
	double sensitivity_consistency_tolerance = 1e-9;

	// 正常运行只输出解析灵敏度。以下有限差分参数暂时保留给源码中的
	// 诊断帮助函数，但同步优化入口不再调用它们。
	int sensitivity_sample_count = 8;
	int finite_difference_mmc_sample_count = 0;
	int finite_difference_material_sample_count = 0;
	double finite_difference_step = 1e-4;
	bool export_results = true;
	std::string output_directory = "mmc_multimaterial_simultaneous_result";
};

// 每个已接受联合状态的标量诊断。这里的 iteration 是已接受更新次数，
// 不是有限元内部迭代次数。
struct MMCMultiMaterialIterationDiagnostic {
	int iteration = 0;
	double compliance = 0.0;
	double compliance_from_load = 0.0;
	double compliance_relative_mismatch = 0.0;
	double projection_beta = 0.0;
	double mmc_occupancy_min = 0.0;
	double mmc_occupancy_max = 0.0;
	double projected_occupancy_min = 0.0;
	double projected_occupancy_max = 0.0;
	double mmc_volume = 0.0;
	// V1 = sum(m_e eta_e v_e), V2 = sum(m_e (1-eta_e) v_e)。
	double V1 = 0.0;
	double V2 = 0.0;
	// material_ratio_residual = V2 - material_ratio_f * V1。
	double material_ratio_residual = 0.0;
	double geometry_constraint = 0.0;
	double positive_ratio_constraint = 0.0;
	double negative_ratio_constraint = 0.0;
	double geometry_max_change = 0.0;
	double material_max_change = 0.0;
	double geometry_gradient_norm = 0.0;
	double material_gradient_norm = 0.0;
	double uniform_material_gradient_relative_error = 0.0;
};

// 对外保留的最终同步状态。各数组均按固定网格单元或 MMC 设计变量顺序存储。
struct MMCMultiMaterialSimultaneousState {
	int accepted_iteration = 0;
	std::vector<double> physical_mmc_variables;
	std::vector<double> material_eta;
	// m_e、E_e 和 a_e 分开保存，便于检查原始几何场和实际刚度。
	std::vector<double> mmc_occupancy;
	std::vector<double> projected_occupancy;
	std::vector<double> material_modulus;
	std::vector<double> effective_stiffness;
	std::vector<double> element_energy;
	std::vector<double> mmc_objective_gradient;
	std::vector<double> material_objective_gradient;
	std::vector<double> geometry_constraint_gradient;
	std::vector<double> ratio_constraint_mmc_gradient;
	std::vector<double> ratio_constraint_material_gradient;
};

struct MMCMultiMaterialSimultaneousResult {
	bool success = false;
	std::string message;
	MMCMultiMaterialSimultaneousState state;
	std::vector<MMCMultiMaterialIterationDiagnostic> diagnostics;
};

// 运行一条全新的、独立持有状态的 MMC + 双材料同步优化循环。
// 本函数不会调用独立 MMC、独立 SIMP 或顺序 MMC+SIMP 的完整优化入口。
// protect_cell_ids 和 neighbor_info 当前仅用于维持统一入口契约及输入检查；
// 材料块沿用顺序 MMC-SIMP 的过滤 + OC 更新；材料变量始终固定在完整
// 网格上，不按 MMC 阈值增删。
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
	std::atomic_bool& is_need_update_model);
