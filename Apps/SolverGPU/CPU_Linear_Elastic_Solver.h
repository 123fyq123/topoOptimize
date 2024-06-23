#pragma once

#include <string>


void CPU_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
                                const std::string& output_file, double Yangs = 1e5, double Poissons = 0.3);

void CPU_Linear_Elastic_Solver_Paridiso(const std::string& input_model_file, const std::string& boundary_file,
                               const std::string& output_file, double Yangs = 1e5, double Poissons = 0.3);

void CPU_Linear_Elastic_Solver_Eigen(const std::string& input_model_file, const std::string& boundary_file,
                                        const std::string& output_file, double Yangs = 1e5, double Poissons = 0.3, bool bUseCG = true);

void CPU_Linear_Elastic_Solver_WN(const std::string& input_model_file, const std::string& boundary_file,
                                  const std::string& output_file, double Yangs = 1e5, double Poissons = 0.3, bool bUseCG = true);