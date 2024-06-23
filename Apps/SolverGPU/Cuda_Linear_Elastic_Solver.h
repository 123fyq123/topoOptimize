//
// Created by XieJin on 10/23/2023.
//

#pragma once

#include <string>


void Cuda_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
                                const std::string& output_file, double Yangs, double Poissons);