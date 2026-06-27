#pragma once
#include <cmath>
#include "Subdivision.h"
#include "FunctionBSpline.h"
#include "FileHandler.h"
#include "BezierSolid.h"
#include <iostream>
#include <vector>
#include<string>
#include<set>
#include<omp.h>
#include"ModelDraw.h"
#include"iGamePoster.h"
#include<windows.h>

void TopOpt(ModelDraw* widget, std::string _filename, std::vector<std::vector<int>>& _pickface, std::vector<vector<int>>& _pickedge,
	std::vector<std::vector<int>>& protectcell);

void TopOpt_CC_Bezier_V1(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model);

void TopOpt_CC_Bezier_V2(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model);

void TopOpt_CC_Bezier_GPUSolver(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model);

void TopOpt_GPU_V0(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model);

void TopOpt_GPU_Discrete(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model);

void TopOpt_GPU_MMC(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	int mmcType,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model);