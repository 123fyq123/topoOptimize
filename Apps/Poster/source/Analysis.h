#pragma once
#include <cmath>
#include "Subdivision.h"
#include "FunctionBSpline.h"
#include "FileHandler.h"
#include "BezierSolid.h"
#include <iostream>
#include <vector>
#include <string>
#include <set>
#include <omp.h>
#include "ModelDraw.h"
#include"iGamePoster.h"
#include <windows.h>

enum class Solver : int
{
    NO_Solver = 0,
    PCG_GPU = 1,
    PCG_CPU_Parallel = 2,
    PCG_CPU_Serial = 3,
    Pardiso_LDLT = 4,
    Eigen_LDLT = 5,
    Eigen_CG = 6,
    Debug = 7
};

void Analysis(ModelDraw* widget, std::string _filename, std::vector<std::vector<int>>& _pickface, std::vector<vector<int>>& _pickedge,
	std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y, std::vector<std::vector<double>>& z);

void Analysis_V1(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<double>& x,
	std::vector<double>& y,
	std::vector<double>& z);

void Analysis_V2(ModelDraw* widget,
                 std::string _filename,
                 IGA::HexBdy& hexbdy,
                 double Yangs,
                 double Poissions,
                 std::vector<double>& x,
                 std::vector<double>& y,
                 std::vector<double>& z,
                 bool bUsingGPU);

void Analysis_V3(ModelDraw* widget,
                 std::string _filename,
                 IGA::HexBdy& hexbdy,
                 double Yangs,
                 double Poissons,
                 Solver SelectedSolver,
                 std::vector<double>& x,
                 std::vector<double>& y,
                 std::vector<double>& z);

void Hyperelasticity_Analysis_V1(ModelDraw* widget,
                 std::string _filename,
                 IGA::HexBdy& hexbdy,
                 double Yangs,
                 double Poissons,
                 int frameCount,
                 double timestep,
                 Solver SelectedSolver,
                 std::vector<double> analysis_data[3],
                 std::atomic_bool &is_need_update_model);