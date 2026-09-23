#pragma once

#include <atomic>
#include <map>
#include <string>
#include <vector>
#include <Eigen/Core>

#include "GP_Spline_v0.h"

struct MultiMaterialSIMPParams {
	double material_ratio_f = 1.5;
	double E1 = 1.0;
	double E2 = 2.0;
	double initial_rho = 0.4;
	double penalty = 3.0;
	double rho_min = 1e-3;
	double move_limit = 0.2;
	double change_tolerance = 0.01;
	bool enable_finite_difference_check = true;
	int finite_difference_sample_count = 6;
	double finite_difference_step = 1e-3;
	int max_iter = 100;
};

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
	std::atomic_bool& is_need_update_model,
	int initial_global_iteration = 0,
	int initial_fe_solve_count = 0,
	double elapsed_time_offset_seconds = 0.0);
