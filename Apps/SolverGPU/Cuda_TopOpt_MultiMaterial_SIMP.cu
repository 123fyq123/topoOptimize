#include "Cuda_TopOpt_MultiMaterial_SIMP.h"

#include <algorithm>
#include <cassert>
#include <chrono>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <numeric>
#include <stdexcept>

#include "OptimizationPaperRecorder.h"

namespace {

void BuildDisplacementByPatch(
	const gpmesh::GPSplineSolidV0& spline,
	const Eigen::MatrixXd& U,
	Eigen::MatrixXd& UU) {
	const int patch_count = spline.m_num_geometry;
	UU.resize(patch_count, 64 * 3);
	for (int p_id = 0; p_id < patch_count; ++p_id) {
		for (int local_id = 0; local_id < 64; ++local_id) {
			int global_id = spline.m_sparse_system.get_l2g_map()[p_id][local_id];
			UU(p_id, local_id) = U(global_id, 0);
			UU(p_id, local_id + 64) = U(global_id + spline.m_sparse_system.getControlNums(), 0);
			UU(p_id, local_id + 128) = U(global_id + 2 * spline.m_sparse_system.getControlNums(), 0);
		}
	}
}

double InterpolateYoungsModulus(double eta, const MultiMaterialSIMPParams& params) {
	// xi=1-eta 是硬材料 E2 的设计变量。对 xi 使用标准 SIMP，
	// 避免旧公式在高模量反差下奖励 0.4~0.6 的中间材料。
	const double hard_fraction = 1.0 - eta;
	return params.E1
		+ std::pow(hard_fraction, params.penalty) * (params.E2 - params.E1);
}

double EvaluateCompliance(
	gpmesh::GPSplineSolidV0& spline,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points,
	const std::vector<double>& mmc_density,
	const std::vector<int>& active_indices,
	const std::vector<Eigen::MatrixXd>& unit_local_matrix,
	const MultiMaterialSIMPParams& params,
	const std::vector<double>& active_rho) {
	std::vector<double> stiffness_coeff(mmc_density.size(), 0.0);
	for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
		const int cell_id = active_indices[local_id];
		stiffness_coeff[cell_id] = mmc_density[cell_id]
			* InterpolateYoungsModulus(active_rho[local_id], params);
	}

	Eigen::MatrixXd U;
	spline.Solver_With_Density_MMC(neumann_points, dirichlet_points, stiffness_coeff, U);

	Eigen::MatrixXd UU;
	BuildDisplacementByPatch(spline, U, UU);

	double compliance = 0.0;
	for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
		const int cell_id = active_indices[local_id];
		const Eigen::VectorXd Ue = UU.row(cell_id);
		const double energy = (Ue.transpose() * unit_local_matrix[cell_id] * Ue)(0, 0);
		compliance += stiffness_coeff[cell_id] * energy;
	}
	return compliance;
}

std::vector<int> BuildFiniteDifferenceSampleIds(int active_count, int requested_count) {
	std::vector<int> sample_ids;
	if (active_count <= 0 || requested_count <= 0) {
		return sample_ids;
	}

	const int sample_count = std::min(active_count, requested_count);
	sample_ids.reserve(sample_count);
	if (sample_count == 1) {
		sample_ids.push_back(active_count / 2);
		return sample_ids;
	}

	for (int i = 0; i < sample_count; ++i) {
		const int local_id = static_cast<int>(std::llround(
			static_cast<double>(i) * (active_count - 1) / (sample_count - 1)));
		if (sample_ids.empty() || sample_ids.back() != local_id) {
			sample_ids.push_back(local_id);
		}
	}
	return sample_ids;
}

std::vector<double> FilterSensitivities(
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	const std::vector<double>& design_variable_full,
	const std::vector<double>& dc_full,
	double rho_min) {
	std::vector<double> filtered(dc_full.size(), 0.0);
	for (int i = 0; i < static_cast<int>(dc_full.size()); ++i) {
		if (design_variable_full[i] <= 0.0) {
			continue;
		}
		double sum_weight = 0.0;
		double accum = 0.0;
		for (const auto& [j, weight] : neighbor_info[i]) {
			if (j < 0 || j >= static_cast<int>(dc_full.size())) {
				continue;
			}
			if (design_variable_full[j] <= 0.0) {
				continue;
			}
			const double clamped_weight = std::max(0.0, weight);
			sum_weight += clamped_weight;
			accum += clamped_weight * design_variable_full[j] * dc_full[j];
		}
		// Sigmund 灵敏度过滤：
		//   dc_hat_i = sum_j(w_ij * rho_j * dc_j) / (rho_i * sum_j(w_ij))
		const double denom = std::max(design_variable_full[i], rho_min)
			* std::max(sum_weight, 1e-12);
		filtered[i] = accum / denom;
	}
	return filtered;
}

bool WriteLegacyHexVTK(
	const std::string& file_path,
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::vector<int>& active_indices,
	const std::vector<double>& material_rho) {
	if (allVertex.empty() || cell2V.empty() || active_indices.empty()) {
		return false;
	}

	std::map<int, int> vertex_id_map;
	std::vector<Eigen::Vector3d> ordered_vertices;
	ordered_vertices.reserve(allVertex.size());
	int vertex_counter = 0;
	for (const auto& [vertex_id, position] : allVertex) {
		vertex_id_map[vertex_id] = vertex_counter++;
		ordered_vertices.push_back(position);
	}

	std::ofstream out(file_path, std::ios::out);
	if (!out.good()) {
		return false;
	}

	out << "# vtk DataFile Version 2.0\n";
	out << "MMC SIMP sequential result\n";
	out << "ASCII\n";
	out << "DATASET UNSTRUCTURED_GRID\n";
	out << "POINTS " << ordered_vertices.size() << " double\n";
	for (const auto& vertex : ordered_vertices) {
		out << vertex.x() << " " << vertex.y() << " " << vertex.z() << "\n";
	}

	out << "CELLS " << active_indices.size() << " " << active_indices.size() * 9 << "\n";
	for (int cell_id : active_indices) {
		const auto cell_it = cell2V.find(cell_id);
		if (cell_it == cell2V.end() || cell_it->second.size() != 8) {
			return false;
		}
		const auto& vertices = cell_it->second;
		out << 8;
		for (const auto& [vertex_id, position] : vertices) {
			(void)position;
			auto map_it = vertex_id_map.find(vertex_id);
			if (map_it == vertex_id_map.end()) {
				return false;
			}
			out << " " << map_it->second;
		}
		out << "\n";
	}

	out << "CELL_TYPES " << active_indices.size() << "\n";
	for (size_t i = 0; i < active_indices.size(); ++i) {
		out << 12 << "\n";
	}

	out << "CELL_DATA " << active_indices.size() << "\n";
	out << "SCALARS material_rho double 1\n";
	out << "LOOKUP_TABLE default\n";
	for (int cell_id : active_indices) {
		if (cell_id < 0 || cell_id >= static_cast<int>(material_rho.size())) {
			return false;
		}
		out << material_rho[cell_id] << "\n";
	}

	// Keep the binary value aligned with material_rho so one color map has one meaning:
	// 0 selects E2 (hard), while 1 selects E1 (soft).
	out << "SCALARS material_binary int 1\n";
	out << "LOOKUP_TABLE default\n";
	for (int cell_id : active_indices) {
		const double rho = std::clamp(material_rho[cell_id], 0.0, 1.0);
		out << (rho >= 0.5 ? 1 : 0) << "\n";
	}

	// Direct RGB fields use the same convention as ParaView's blue-to-red scalar map:
	// E2/hard is blue and E1/soft is red.
	out << "COLOR_SCALARS material_color 3\n";
	for (int cell_id : active_indices) {
		const double rho = std::clamp(material_rho[cell_id], 0.0, 1.0);
		out << rho << " 0 " << (1.0 - rho) << "\n";
	}

	out << "COLOR_SCALARS material_binary_color 3\n";
	for (int cell_id : active_indices) {
		const double rho = std::clamp(material_rho[cell_id], 0.0, 1.0);
		const double soft_material = rho >= 0.5 ? 1.0 : 0.0;
		out << soft_material << " 0 " << (1.0 - soft_material) << "\n";
	}
	return true;
}

void OCUpdate(
	std::vector<double>& x,
	double material_ratio_f,
	const std::vector<double>& dc,
	const std::vector<double>& material_weight,
	double rho_min,
	double move_limit) {
	double l1 = 0.0;
	double l2 = 100000.0;
	const int nn = static_cast<int>(x.size());
	if (static_cast<int>(dc.size()) != nn ||
		static_cast<int>(material_weight.size()) != nn) {
		throw std::runtime_error(
			"The sequential material OC input sizes do not match.");
	}
	std::vector<double> xnew(x);
	const double total_material_weight = std::accumulate(
		material_weight.begin(),
		material_weight.end(),
		0.0);
	const double target_hard_volume =
		total_material_weight * material_ratio_f
		/ (1.0 + material_ratio_f);
	while ((l2 - l1) / (l1 + l2) > 1e-4) {
		const double lmid = 0.5 * (l1 + l2);
		double hard_volume = 0.0;
		for (int i = 0; i < nn; ++i) {
			// dc=dC/dxi<0。OC 在硬材料变量 xi=1-eta 上更新，
			// 外部仍保存 eta，以保持文件字段和红蓝材料语义不变。
			const double hard_before = 1.0 - x[i];
			const double positive_weight = std::max(
				material_weight[i],
				std::numeric_limits<double>::min());
			// 材料资源对 xi 的导数为 m_e，标准 OC 比值必须除以该权重。
			const double hard_candidate =
				std::max(hard_before, rho_min)
				* std::sqrt(
					std::max(-dc[i], std::numeric_limits<double>::min())
					/ (lmid * positive_weight));
			const double hard_after = std::max(
				0.0,
				std::max(
					hard_before - move_limit,
					std::min(
						1.0 - rho_min,
						std::min(hard_before + move_limit, hard_candidate))));
			xnew[i] = 1.0 - hard_after;
			hard_volume += material_weight[i] * hard_after;
		}
		// f=V2/V1，因此 MMC 加权硬材料目标为 sum(m_e)*f/(1+f)。
		if (hard_volume > target_hard_volume) {
			l1 = lmid;
		}
		else {
			l2 = lmid;
		}
	}
	x = xnew;
}

}

void RunMultiMaterialSIMPStage(
	gpmesh::GPSplineSolidV0& spline,
	std::vector<boundary_point>& neumann_points,
	std::vector<boundary_point>& dirichlet_points,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::vector<double>& mmc_density,
	const std::vector<int>& active_indices,
	const std::vector<Eigen::MatrixXd>& unit_local_matrix,
	const MultiMaterialSIMPParams& params,
	const std::string& output_dir,
	const std::string& record_dir,
	std::vector<double>& material_density_full,
	std::atomic_bool& is_need_update_model,
	int initial_global_iteration,
	int initial_fe_solve_count,
	double elapsed_time_offset_seconds) {
	const int grid_size = static_cast<int>(mmc_density.size());
	material_density_full.assign(grid_size, 0.0);
	if (active_indices.empty()) {
		return;
	}

	std::filesystem::create_directories(output_dir);
	std::filesystem::create_directories(record_dir);
	optimization_paper_record::IterationWriter paper_output(
		record_dir + "/paper_iterations.tsv",
		true);
	if (!paper_output.Good()) {
		throw std::runtime_error("Cannot create the sequential paper iteration record.");
	}
	const auto material_stage_start = std::chrono::steady_clock::now();
	std::ofstream analytical_out(record_dir + "/analytical_sensitivities.txt", std::ios::out | std::ios::trunc);
	std::ofstream finite_difference_out(record_dir + "/finite_difference_sensitivities.txt", std::ios::out | std::ios::trunc);
	analytical_out << std::setprecision(17);
	finite_difference_out << std::setprecision(17);
	analytical_out
		<< "iter\tlocal_id\tcell_id\trho\tmmc_density\telement_energy"
		<< "\tobjective_gradient_raw\toc_sensitivity_raw\toc_sensitivity_filtered\n";
	finite_difference_out
		<< "iter\tlocal_id\tcell_id\trho\tscheme\th\tcompliance_minus\tcompliance_base"
		<< "\tcompliance_plus\tfinite_difference\tobjective_gradient_raw"
		<< "\toc_sensitivity_raw\toc_sensitivity_filtered"
		<< "\tabs_error\trelative_error\n";

	// Match the MATLAB convention: rho is the E1 (basic/soft material) variable.
	// rho close to 1 selects E1; rho close to rho_min selects E2.
	// active_rho:是真正参与 SIMP 优化的材料选择变量，只在 MMC 筛出来的有效单元里存在。
	std::vector<double> active_rho(active_indices.size(), std::clamp(params.initial_rho, params.rho_min, 1.0));
	std::vector<double> active_material_weight(active_indices.size(), 0.0);
	for (int local_id = 0;
		local_id < static_cast<int>(active_indices.size());
		++local_id) {
		active_material_weight[local_id] = mmc_density[active_indices[local_id]];
	}
	std::vector<double> final_rho_full(grid_size, 0.0);

	for (int iter = 0; iter < params.max_iter; ++iter) {
		const auto iteration_start = std::chrono::steady_clock::now();
		std::vector<double> hard_fraction_full(grid_size, 0.0);
		std::vector<double> stiffness_coeff(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			assert(cell_id >= 0 && cell_id < grid_size);
			const double rho = active_rho[local_id];

			// eta=rho 仍表示软材料 E1，xi=1-eta 表示硬材料 E2：
			// E(eta_e) = E1 + (1-eta_e)^p * (E2-E1)。
			const double interpolated_E = InterpolateYoungsModulus(rho, params);
			hard_fraction_full[cell_id] = std::max(1.0 - rho, params.rho_min);

			// 顺序 MMC + SIMP 的单元刚度系数：
			// k_e = rho_e^M * E(rho_e) * k0
			// 之后复用 MMC 求解器的线性装配，把该系数直接乘到单元刚度上。
			stiffness_coeff[cell_id] = mmc_density[cell_id] * interpolated_E;
		}

		Eigen::MatrixXd U;
		spline.Solver_With_Density_MMC(neumann_points, dirichlet_points, stiffness_coeff, U);

		Eigen::MatrixXd UU;
		BuildDisplacementByPatch(spline, U, UU);

		double comp = 0.0;
		// objective_gradient_active=dC/deta>0；由于 xi=1-eta，
		// dc_active=dC/dxi=-dC/deta<0，可直接用于标准 OC。
		std::vector<double> dc_active(active_indices.size(), 0.0);
		std::vector<double> objective_gradient_active(active_indices.size(), 0.0);
		std::vector<double> energy_active(active_indices.size(), 0.0);
		// dc_full 是把 dc_active 映射回完整网格后的数组。
		std::vector<double> dc_full(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			const Eigen::VectorXd Ue = UU.row(cell_id);
			// 这里实际是U^T * K0 * U
			const double energy = (Ue.transpose() * unit_local_matrix[cell_id] * Ue)(0, 0);
			energy_active[local_id] = energy;
			comp += stiffness_coeff[cell_id] * energy;
			// dE/deta=-p(1-eta)^(p-1)(E2-E1)，故柔顺度伴随导数
			// dC/deta=-m*q*dE/deta 为正。
			objective_gradient_active[local_id] = mmc_density[cell_id]
				* params.penalty
				* std::pow(1.0 - active_rho[local_id], params.penalty - 1.0)
				* (params.E2 - params.E1)
				* energy;
			dc_active[local_id] = -objective_gradient_active[local_id];
			dc_full[cell_id] = dc_active[local_id];
		}
		const std::vector<double> raw_dc_active = dc_active;

		const std::vector<double> filtered_full = FilterSensitivities(
			neighbor_info,
			hard_fraction_full,
			dc_full,
			params.rho_min);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			dc_active[local_id] = filtered_full[active_indices[local_id]];
		}

		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			analytical_out << iter + 1 << '\t'
				<< local_id << '\t'
				<< cell_id << '\t'
				<< active_rho[local_id] << '\t'
				<< mmc_density[cell_id] << '\t'
				<< energy_active[local_id] << '\t'
				<< objective_gradient_active[local_id] << '\t'
				<< raw_dc_active[local_id] << '\t'
				<< dc_active[local_id] << '\n';
		}
		analytical_out.flush();

		// 有限差分
		if (0) {
			const std::vector<int> sample_ids = BuildFiniteDifferenceSampleIds(
				static_cast<int>(active_indices.size()), params.finite_difference_sample_count);
			const double requested_h = std::max(params.finite_difference_step, 1e-8);
			const double nan = std::numeric_limits<double>::quiet_NaN();

			for (int local_id : sample_ids) {
				const double rho = active_rho[local_id];
				double compliance_minus = nan;
				double compliance_plus = nan;
				double finite_difference = nan;
				double h = requested_h;
				const char* scheme = "central";

				if (rho - requested_h >= params.rho_min && rho + requested_h <= 1.0) {
					std::vector<double> rho_minus = active_rho;
					std::vector<double> rho_plus = active_rho;
					rho_minus[local_id] -= h;
					rho_plus[local_id] += h;
					compliance_minus = EvaluateCompliance(spline, neumann_points, dirichlet_points,
						mmc_density, active_indices, unit_local_matrix, params, rho_minus);
					compliance_plus = EvaluateCompliance(spline, neumann_points, dirichlet_points,
						mmc_density, active_indices, unit_local_matrix, params, rho_plus);
					finite_difference = (compliance_plus - compliance_minus) / (2.0 * h);
				}
				else if (rho + requested_h <= 1.0) {
					scheme = "forward";
					std::vector<double> rho_plus = active_rho;
					rho_plus[local_id] += h;
					compliance_plus = EvaluateCompliance(spline, neumann_points, dirichlet_points,
						mmc_density, active_indices, unit_local_matrix, params, rho_plus);
					finite_difference = (compliance_plus - comp) / h;
				}
				else {
					scheme = "backward";
					h = std::min(requested_h, rho - params.rho_min);
					std::vector<double> rho_minus = active_rho;
					rho_minus[local_id] -= h;
					compliance_minus = EvaluateCompliance(spline, neumann_points, dirichlet_points,
						mmc_density, active_indices, unit_local_matrix, params, rho_minus);
					finite_difference = (comp - compliance_minus) / h;
				}

				const double abs_error = std::abs(finite_difference - objective_gradient_active[local_id]);
				const double error_scale = std::max({
					std::abs(finite_difference), std::abs(objective_gradient_active[local_id]), 1e-12 });
				const double relative_error = abs_error / error_scale;
				finite_difference_out << iter + 1 << '\t'
					<< local_id << '\t'
					<< active_indices[local_id] << '\t'
					<< rho << '\t'
					<< scheme << '\t'
					<< h << '\t'
					<< compliance_minus << '\t'
					<< comp << '\t'
					<< compliance_plus << '\t'
					<< finite_difference << '\t'
					<< objective_gradient_active[local_id] << '\t'
					<< raw_dc_active[local_id] << '\t'
					<< dc_active[local_id] << '\t'
					<< abs_error << '\t'
					<< relative_error << '\n';
			}
			finite_difference_out.flush();
		}

		const std::vector<double> old_rho = active_rho;
		OCUpdate(
			active_rho,
			params.material_ratio_f,
			dc_active,
			active_material_weight,
			params.rho_min,
			params.move_limit);

		double change = 0.0;
		material_density_full.assign(grid_size, 0.0);
		final_rho_full.assign(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			change = std::max(change, std::abs(active_rho[local_id] - old_rho[local_id]));
			const int cell_id = active_indices[local_id];
			material_density_full[cell_id] = mmc_density[cell_id];
			final_rho_full[cell_id] = active_rho[local_id];
		}

		// 论文记录与 OC 资源约束统一采用 MMC 占有率加权的材料体积，
		// 保证记录的 V2/V1 与优化器实际维持的比例完全同义。
		double material_v1 = 0.0;
		double material_v2 = 0.0;
		double mean_material_variable = 0.0;
		double material_grayness = 0.0;
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			// 柔顺度和梯度均对应 OC 更新前状态，因此材料统计也使用 old_rho。
			const double rho = old_rho[local_id];
			material_v1 += mmc_density[cell_id] * rho / static_cast<double>(grid_size);
			material_v2 += mmc_density[cell_id] * (1.0 - rho) / static_cast<double>(grid_size);
			mean_material_variable += rho;
			material_grayness += 4.0 * rho * (1.0 - rho);
		}
		mean_material_variable /= static_cast<double>(active_indices.size());
		material_grayness /= static_cast<double>(active_indices.size());
		const double mmc_volume = std::accumulate(
			mmc_density.begin(), mmc_density.end(), 0.0)
			/ static_cast<double>(grid_size);

		optimization_paper_record::IterationRecord paper_record;
		paper_record.method = "sequential";
		paper_record.stage = "simp_material";
		paper_record.global_iteration = initial_global_iteration + iter + 1;
		paper_record.stage_iteration = iter + 1;
		paper_record.compliance = comp;
		paper_record.mmc_volume = mmc_volume;
		paper_record.material_v1 = material_v1;
		paper_record.material_v2 = material_v2;
		paper_record.material_ratio = material_v1 > 1e-15
			? material_v2 / material_v1
			: optimization_paper_record::NaN();
		paper_record.ratio_residual =
			material_v2 - params.material_ratio_f * material_v1;
		paper_record.mean_material_variable = mean_material_variable;
		paper_record.material_grayness = material_grayness;
		paper_record.material_change = change;
		paper_record.material_gradient_norm =
			optimization_paper_record::L2Norm(objective_gradient_active);
		paper_record.active_cell_count = static_cast<int>(active_indices.size());
		paper_record.cumulative_fe_solves = initial_fe_solve_count + iter + 1;
		paper_record.iteration_seconds =
			optimization_paper_record::ElapsedSeconds(iteration_start);
		paper_record.elapsed_seconds = elapsed_time_offset_seconds
			+ optimization_paper_record::ElapsedSeconds(material_stage_start);
		paper_output.Write(paper_record);

		std::cout << "[MMC+SIMP] It.: " << iter + 1
			<< " Obj.: " << comp
			<< " Vol.: " << std::accumulate(active_rho.begin(), active_rho.end(), 0.0) / std::max<int>(1, static_cast<int>(active_rho.size()))
			<< " ch.: " << change << std::endl;
		is_need_update_model = true;

		if (change < params.change_tolerance) {
			break;
		}
	}

	const std::string vtk_name = output_dir + "/mmc_simp_result.vtk";
	// Export only the MMC-selected SIMP design domain. material_rho remains continuous:
	// MATLAB convention: rho close to 1 selects E1 (soft), while rho close to 0 selects E2 (hard).
	WriteLegacyHexVTK(vtk_name, allVertex, cell2V, active_indices, final_rho_full);
}
