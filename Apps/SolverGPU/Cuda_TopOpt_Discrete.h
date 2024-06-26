#pragma once

#include <string>
#include "Cuda_Linear_Elastic_Solver.h"
#include<atomic>

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