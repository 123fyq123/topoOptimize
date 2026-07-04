#include "Cuda_TopOpt_MultiMaterial_SIMP.h"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <numeric>
#include <sstream>

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

bool WriteLegacyHexVTK(
	const std::string& file_path,
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::vector<double>& mmc_density,
	const std::vector<double>& material_rho,
	const std::vector<double>& material_rho_full,
	const std::vector<double>& display_density,
	const std::vector<double>& stiffness_coeff,
	const std::vector<double>& effective_coeff,
	const std::vector<double>& material_label) {
	if (allVertex.empty() || cell2V.empty()) {
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

	out << "CELLS " << cell2V.size() << " " << cell2V.size() * 9 << "\n";
	for (const auto& [cell_id, vertices] : cell2V) {
		(void)cell_id;
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

	out << "CELL_TYPES " << cell2V.size() << "\n";
	for (size_t i = 0; i < cell2V.size(); ++i) {
		out << 12 << "\n";
	}

	out << "CELL_DATA " << cell2V.size() << "\n";
	auto write_scalar_field = [&out](const std::string& name, const std::vector<double>& values) {
		out << "SCALARS " << name << " double 1\n";
		out << "LOOKUP_TABLE default\n";
		for (double value : values) {
			out << value << "\n";
		}
	};

	write_scalar_field("mmc_density", mmc_density);
	write_scalar_field("material_rho", material_rho);
	write_scalar_field("material_rho_full", material_rho_full);
	write_scalar_field("display_density", display_density);
	write_scalar_field("stiffness_coeff", stiffness_coeff);
	write_scalar_field("effective_coeff", effective_coeff);
	write_scalar_field("material_label", material_label);
	return true;
}

std::vector<double> BuildBinaryMaterialLabelField(
	const std::vector<double>& mmc_density,
	const std::vector<double>& material_rho_full,
	double mmc_threshold,
	double material_threshold) {
	std::vector<double> material_label(mmc_density.size(), 0.0);
	for (int i = 0; i < static_cast<int>(mmc_density.size()); ++i) {
		if (mmc_density[i] <= mmc_threshold) {
			material_label[i] = 0.0;
			continue;
		}
		material_label[i] = (material_rho_full[i] >= material_threshold) ? 2.0 : 1.0;
	}
	return material_label;
}

void OCUpdate(
	std::vector<double>& x,
	double material_ratio_f,
	const std::vector<double>& dc,
	double rho_min,
	double move_limit) {
	double l1 = 0.0;
	double l2 = 100000.0;
	const int nn = static_cast<int>(x.size());
	std::vector<double> xnew(nn, rho_min);
	while (l2 - l1 > 1e-4) {
		const double lmid = 0.5 * (l1 + l2);
		for (int i = 0; i < nn; ++i) {
			// OC 更新公式：
			//   rho_new = rho * sqrt(-dc / lambda)
			// 随后施加 move limit 以及 [rho_min, 1] 上下界。
			const double candidate = x[i] * std::sqrt(std::max(0.0, -dc[i] / std::max(lmid, 1e-12)));
			xnew[i] = std::max(rho_min, std::max(x[i] - move_limit, std::min(1.0, std::min(x[i] + move_limit, candidate))));
		}
		double sum_x = 0.0;
		for (double value : xnew) {
			sum_x += value;
		}
		// material_ratio_f = V2 / V1。
		// rho 表示两相插值中基体材料 E1 的占比，因此目标总量为 nn / (1 + f)。
		if (sum_x - nn / (1.0 + material_ratio_f) > 0.0) {
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
	std::vector<double>& material_density_full,
	std::atomic_bool& is_need_update_model) {
	const int grid_size = static_cast<int>(mmc_density.size());
	material_density_full.assign(grid_size, 0.0);
	if (active_indices.empty()) {
		return;
	}

	std::filesystem::create_directories(output_dir);

	// rho_e 是固定 MMC 骨架内部的材料选择变量：
	// rho_e 越接近 1 越偏向 E1；rho_e 越接近 rho_min 越偏向 E2。
	// active_rho:是真正参与 SIMP 优化的材料选择变量，只在 MMC 筛出来的有效单元里存在。
	std::vector<double> active_rho(active_indices.size(), std::clamp(params.initial_rho, params.rho_min, 1.0));

	for (int iter = 0; iter < params.max_iter; ++iter) {
		std::vector<double> rho_full(grid_size, 0.0);
		std::vector<double> stiffness_coeff(grid_size, 0.0);
		std::vector<double> effective_coeff(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			assert(cell_id >= 0 && cell_id < grid_size);
			const double rho = active_rho[local_id];

			// 两相材料插值：
			// E(rho_e) = rho_e^p * E1 + (1 - rho_e^p) * E2
			const double rho_power_p = std::pow(rho, params.penalty);
			const double interpolated_E = rho_power_p * params.E1 + (1.0 - rho_power_p) * params.E2;
			rho_full[cell_id] = rho;

			// 顺序 MMC + SIMP 的单元刚度系数：
			// k_e = rho_e^M * E(rho_e) * k0
			// 之后 Solver_With_MultiMaterial_Density(...) 会拿它去装配刚度矩阵。
			stiffness_coeff[cell_id] = mmc_density[cell_id] * interpolated_E;
			// effective_coeff 目前只是 stiffness_coeff 的拷贝，用于导出 VTK 结果查看。
			effective_coeff[cell_id] = stiffness_coeff[cell_id];
		}

		Eigen::MatrixXd U;
		spline.Solver_With_MultiMaterial_Density(neumann_points, dirichlet_points, stiffness_coeff, U);

		Eigen::MatrixXd UU;
		BuildDisplacementByPatch(spline, U, UU);

		double comp = 0.0;
		// dc_active 只对应 MMC 筛出来的有效单元，也就是 active_indices 里的单元。
		std::vector<double> dc_active(active_indices.size(), 0.0);
		// dc_full 是把 dc_active 映射回完整网格后的数组。
		std::vector<double> dc_full(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			const int cell_id = active_indices[local_id];
			const Eigen::VectorXd Ue = UU.row(cell_id);
			const double energy = (Ue.transpose() * unit_local_matrix[cell_id] * Ue)(0, 0);
			comp += energy;
			// MATLAB 参考代码：
			//   dc = penal * rho^(penal-1) * (E1-E2) * Ue' * k0 * Ue
			// 固定的 MMC 因子 rhoM 已在上面的刚度装配中使用；
			// 这里不再乘入 dc，以保持和论文 MATLAB 代码一致。
			dc_active[local_id] = params.penalty * std::pow(active_rho[local_id], params.penalty - 1.0) * (params.E1 - params.E2) * energy;
			dc_full[cell_id] = dc_active[local_id];
		}

		const std::vector<double> filtered_full = FilterSensitivities(neighbor_info, rho_full, dc_full, params.rho_min);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			dc_active[local_id] = filtered_full[active_indices[local_id]];
		}

		const std::vector<double> old_rho = active_rho;
		OCUpdate(active_rho, params.material_ratio_f, dc_active, params.rho_min, params.move_limit);

		double change = 0.0;
		material_density_full.assign(grid_size, 0.0);
		for (int local_id = 0; local_id < static_cast<int>(active_indices.size()); ++local_id) {
			change = std::max(change, std::abs(active_rho[local_id] - old_rho[local_id]));
			const int cell_id = active_indices[local_id];
			material_density_full[cell_id] = mmc_density[cell_id];
		}

		std::ostringstream vtk_name;
		vtk_name << output_dir << "/mmc_simp_iter_" << std::setw(4) << std::setfill('0') << (iter + 1) << ".vtk";
		// VTK 字段与 MATLAB 变量的对应关系：
		// material_rho / material_rho_full 对应 Sxy00 按 Ei 映射回全网格后的 col。
		// effective_coeff 对应用于装配的 Ndenk。
		const std::vector<double> material_label = BuildBinaryMaterialLabelField(mmc_density, rho_full, 1e-3, 0.5);
		WriteLegacyHexVTK(vtk_name.str(), allVertex, cell2V, mmc_density, rho_full, rho_full, material_density_full, stiffness_coeff, effective_coeff, material_label);

		std::cout << "[MMC+SIMP] It.: " << iter + 1
			<< " Obj.: " << comp
			<< " Vol.: " << std::accumulate(active_rho.begin(), active_rho.end(), 0.0) / std::max<int>(1, static_cast<int>(active_rho.size()))
			<< " ch.: " << change << std::endl;
		is_need_update_model = true;

		if (change < params.change_tolerance) {
			break;
		}
	}
}
