#pragma once

#include <atomic>
#include <map>
#include <string>
#include <vector>

#include <Eigen/Core>

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
	int mmcType);
