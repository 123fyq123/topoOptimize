#pragma once

#include <string>
#include <vector>
#include <atomic>

void CPU_Hyperelasticity_TimeIntegrator_Solver(const std::string& input_model_file,
                                               const std::string& boundary_file,
                                               const std::string& output_file,
                                               double Yangs, double Poissons,
                                               int frameCount, double timeStep,
                                               HexV3fMesh* mesh,
                                               std::vector<double> analysis_data[3],
                                               std::atomic_bool &is_need_update_model);