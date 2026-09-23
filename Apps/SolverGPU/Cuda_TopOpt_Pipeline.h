//
//
// Created by XieJin on 10/27/2023.
//

#pragma once

#include <string>
#include "Cuda_Linear_Elastic_Solver.h"
#include "SolverGPU/Cuda_TopOpt_Pipeline.h"
#include<atomic>
#include <Eigen/Core>
#include <map>
void Cuda_TopOpt_Pipeline(const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model);

void Cuda_TopOpt_MMC_Node(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model,
	int mmcType,
	const std::string& paper_record_path = std::string(),
	int* completed_iterations = nullptr);

//0-1¿Î…¢Õÿ∆À”≈ªØ
void Cuda_TopOpt_Discrete(const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model);
