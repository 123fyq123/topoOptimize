//
// Created by XieJin on 10/27/2023.
//

#include <filesystem>
#include"TopOpt.h"
#include "IGA_CC_Helper.h"
#include "CatmullClarkBezierSolver.h"
#include "Timer.h"
#include "SolverGPU/Cuda_TopOpt_Pipeline.h"
#include "SolverGPU/Cuda_TopOpt_MMC_SIMP_Sequential.h"
#include "SolverGPU/Cuda_TopOpt_MMC_MultiMaterial_Simultaneous.h"
#include "SolverGPU/Cuda_TopOpt_MultiRes.h"

void TopOpt_GPU_V0(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poisson,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model) {

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	std::cout << _filename << "\n";

	BSplineSolid m_solids;
	FileHandler file(_filename);

	std::string model_name = IGA_CC_Helper::get_model_name(_filename);
	HexV3fMesh* mesh = new HexV3fMesh();
	file.read2HexMeshObj(mesh, true, true, true);
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	m_solids = AppCCSIMP.bss;

	std::filesystem::path temp_path("./temp/");
	if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

	std::string bezier_file = "./temp/" + model_name + ".xml";
	std::string bdy_file = "./temp/" + model_name + ".bdy";
	std::string neighbor_file = "./temp/" + model_name + ".tinfo";  // <neighbor_id; weight>
	std::string density_file = "./temp/" + model_name + ".density";
	std::string protect_cells_file = "./temp/" + model_name + ".prot";


	IGA_CC_Helper::write2Bezier(bezier_file, m_solids);    // xml

	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert_V2(_filename, hexbdy, ccbdy, AppCCSIMP);
	IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);  //bdy

	// ʩ��filter��Ҫ����Ϣ
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	std::vector<std::vector<std::pair<int, double>>> neighbor_info;

	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = std::numeric_limits<float>::max();
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}

		rmin = rmin * 1.5;
		//std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		std::vector<std::pair<int, double>> temp_neighbor_info;
		for (int j = 0; j < temp_aroundcell.size(); ++j) {
			temp_neighbor_info.emplace_back(std::pair<int, double>{temp_aroundcell[j], temp_aroundcell_weight[j]});
		}
		neighbor_info.emplace_back(temp_neighbor_info);
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		//std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}



	{
		std::cout << "protected_cell size = " << pro.size() << "\n";
		Cuda_TopOpt_Pipeline(bezier_file, bdy_file, pro, neighbor_info,
			target_rate, max_iter, Yangs, Poisson, density, is_need_update_topopt_model);
		/*MultiRes_TopOpt_Pipeline_10(bezier_file, bdy_file, pro, neighbor_info,
			target_rate, max_iter, Yangs, Poisson, 0, density, is_need_update_topopt_model);*/
			//        MultiRes_TopOpt_Pipeline(bezier_file, bdy_file, pro, neighbor_info,
			//                             target_rate, max_iter, Yangs, Poisson, 1, density, is_need_update_topopt_model);
	}

}


void TopOpt_GPU_MMC(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	int mmcType,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poisson,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model) {

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	std::cout << _filename << "\n";

	BSplineSolid m_solids;
	FileHandler file(_filename);

	std::string model_name = IGA_CC_Helper::get_model_name(_filename);
	HexV3fMesh* mesh = new HexV3fMesh();
	file.read2HexMeshObj(mesh, true, true, true);
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	m_solids = AppCCSIMP.bss;

	std::filesystem::path temp_path("./temp/");
	if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

	std::string bezier_file = "./temp/" + model_name + ".xml";
	std::string bdy_file = "./temp/" + model_name + ".bdy";
	std::string neighbor_file = "./temp/" + model_name + ".tinfo";  // <neighbor_id; weight>
	std::string density_file = "./temp/" + model_name + ".density";
	std::string protect_cells_file = "./temp/" + model_name + ".prot";


	IGA_CC_Helper::write2Bezier(bezier_file, m_solids);    // xml

	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert_V2(_filename, hexbdy, ccbdy, AppCCSIMP);
	IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);  //bdy

	// ʩ��filter��Ҫ����Ϣ
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	std::vector<std::vector<std::pair<int, double>>> neighbor_info;
	std::map<int, Eigen::Vector3d> cell2Center;
	std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>> cell2V;
	std::map<int, Eigen::Vector3d> allVertex;
	for (auto i = 0; i < mesh->vertices_.size(); ++i)
	{
		auto v = mesh->vertices_[i];
		allVertex[i] = Eigen::Vector3d{ v.x, v.y, v.z };
	}

	for (auto i = 0; i < mesh->cells_.size(); ++i)
	{
		cell2Center[i] = { mesh->cell_center[i].x ,mesh->cell_center[i].y , mesh->cell_center[i].z };
		cell2V[i].resize(8);
		auto cellV = mesh->cells_[i].vertex_;
		for (int j = 0; j < 8; ++j)
		{
			int idx = cellV[j].idx();
			auto v = mesh->vertices_[idx];
			cell2V[i][j] = std::make_pair(idx, Eigen::Vector3d{ v.x, v.y, v.z });
		}
	}

	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = std::numeric_limits<float>::max();
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}

		rmin = rmin * 1.5;
		//std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		std::vector<std::pair<int, double>> temp_neighbor_info;
		for (int j = 0; j < temp_aroundcell.size(); ++j) {
			temp_neighbor_info.emplace_back(std::pair<int, double>{temp_aroundcell[j], temp_aroundcell_weight[j]});
		}
		neighbor_info.emplace_back(temp_neighbor_info);
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		//std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}



	{
		std::cout << "protected_cell size = " << pro.size() << "\n";
		Cuda_TopOpt_MMC_Node(allVertex, cell2V, bezier_file, bdy_file, pro, neighbor_info, target_rate, max_iter, Yangs, Poisson, density, is_need_update_topopt_model, mmcType);
	}
}

void TopOpt_GPU_Discrete(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poisson,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model) {

	Timer All;

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	std::cout << _filename << "\n";

	BSplineSolid m_solids;
	FileHandler file(_filename);

	std::string model_name = IGA_CC_Helper::get_model_name(_filename);
	HexV3fMesh* mesh = new HexV3fMesh();
	file.read2HexMeshObj(mesh, true, true, true);

	All.out("construct_bss_new start");
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	m_solids = AppCCSIMP.bss;
	All.out("construct_bss_new end");

	std::filesystem::path temp_path("./temp/");
	if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

	std::string bezier_file = "./temp/" + model_name + ".xml";
	std::string bdy_file = "./temp/" + model_name + ".bdy";
	std::string neighbor_file = "./temp/" + model_name + ".tinfo";  // <neighbor_id; weight>
	std::string density_file = "./temp/" + model_name + ".density";
	std::string protect_cells_file = "./temp/" + model_name + ".prot";

	IGA_CC_Helper::write2Bezier(bezier_file, m_solids);    // xml

	IGA::CCBezierBdy ccbdy;
	All.out("BdyConvert start");
	//IGA::BdyConvert(_filename, hexbdy, ccbdy);
	//IGA::BdyConvert_V2(_filename, hexbdy, ccbdy, AppCCSIMP);
	IGA::BdyConvert_V3(mesh, hexbdy, ccbdy, AppCCSIMP);
	All.out("BdyConvert end");

	IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);  //bdy

	// ʩ��filter��Ҫ����Ϣ
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	std::vector<std::vector<std::pair<int, double>>> neighbor_info;

	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = std::numeric_limits<float>::max();
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}

		rmin = rmin * 1.5;
		//std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		std::vector<std::pair<int, double>> temp_neighbor_info;
		for (int j = 0; j < temp_aroundcell.size(); ++j) {
			temp_neighbor_info.emplace_back(std::pair<int, double>{temp_aroundcell[j], temp_aroundcell_weight[j]});
		}
		neighbor_info.emplace_back(temp_neighbor_info);
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		//std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}

	All.out("Cuda_TopOpt_Discrete start");
	{
		std::cout << "protected_cell size = " << pro.size() << "\n";
		Cuda_TopOpt_Discrete(bezier_file, bdy_file, pro, neighbor_info,
			target_rate, max_iter, Yangs, Poisson, density, is_need_update_topopt_model);
	}
	All.out("Cuda_TopOpt_Discrete end");
	std::cout << _filename << "\n";
}

void TopOpt_GPU_MMC_SIMP_Sequential(std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	int mmcType,
	double mmc_target_rate,
	double mmc_max_iter,
	double material_ratio_f,
	double simp_max_iter,
	double E1,
	double E2,
	double Poisson,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model) {
	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	BSplineSolid m_solids;
	FileHandler file(_filename);
	std::string model_name = IGA_CC_Helper::get_model_name(_filename);
	HexV3fMesh* mesh = new HexV3fMesh();
	file.read2HexMeshObj(mesh, true, true, true);
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	m_solids = AppCCSIMP.bss;

	std::filesystem::path temp_path("./temp/");
	if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

	std::string bezier_file = "./temp/" + model_name + ".xml";
	std::string bdy_file = "./temp/" + model_name + ".bdy";
	IGA_CC_Helper::write2Bezier(bezier_file, m_solids);

	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert_V2(_filename, hexbdy, ccbdy, AppCCSIMP);
	IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);

	std::vector<std::vector<std::pair<int, double>>> neighbor_info;
	std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>> cell2V;
	std::map<int, Eigen::Vector3d> allVertex;
	for (int i = 0; i < static_cast<int>(mesh->vertices_.size()); ++i) {
		auto v = mesh->vertices_[i];
		allVertex[i] = Eigen::Vector3d{ v.x, v.y, v.z };
	}
	for (int i = 0; i < static_cast<int>(mesh->cells_.size()); ++i) {
		cell2V[i].resize(8);
		auto cellV = mesh->cells_[i].vertex_;
		for (int j = 0; j < 8; ++j) {
			int idx = cellV[j].idx();
			auto v = mesh->vertices_[idx];
			cell2V[i][j] = std::make_pair(idx, Eigen::Vector3d{ v.x, v.y, v.z });
		}
	}

	for (int i = 0; i < static_cast<int>(m_solids.Local_V3f.size()); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = std::numeric_limits<float>::max();
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (auto c_c_it = c_c.begin(); c_c_it != c_c.end(); ++c_c_it) {
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}
		rmin = rmin * 1.5;
		for (int k = 0; k < static_cast<int>(mesh->cells_.size()); ++k) {
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin) {
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		std::vector<std::pair<int, double>> temp_neighbor_info;
		for (int j = 0; j < static_cast<int>(temp_aroundcell.size()); ++j) {
			temp_neighbor_info.emplace_back(std::pair<int, double>{temp_aroundcell[j], temp_aroundcell_weight[j]});
		}
		neighbor_info.emplace_back(temp_neighbor_info);
	}

	Cuda_TopOpt_MMC_SIMP_Sequential(allVertex, cell2V, bezier_file, bdy_file, pro, neighbor_info,
		mmc_target_rate, mmc_max_iter, material_ratio_f, simp_max_iter, E1, E2, Poisson,
		density, is_need_update_topopt_model, mmcType);

}
// UI/模型层到新增同步求解器的独立适配入口。
// 本函数只负责网格、Bezier、边界条件和参数转换；数值迭代全部位于
// Cuda_TopOpt_MMC_MultiMaterial_Simultaneous.cu，不会转入原 MMC 或顺序流程。
void TopOpt_GPU_MMC_MultiMaterial_Simultaneous(
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::vector<int>& pro,
	int mmcType,
	double mmc_target_rate,
	double max_iter,
	double material_ratio_f,
	double E1,
	double E2,
	double Poisson,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_topopt_model) {
	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	BSplineSolid solids;
	FileHandler file(_filename);
	const std::string model_name = IGA_CC_Helper::get_model_name(_filename);
	HexV3fMesh* mesh = new HexV3fMesh();
	file.read2HexMeshObj(mesh, true, true, true);
	ApproximationCC approximation(mesh, nullptr);
	approximation.construct_bss_new();
	solids = approximation.bss;

	const std::filesystem::path temp_path("./temp/");
	if (!std::filesystem::exists(temp_path)) {
		std::filesystem::create_directory(temp_path);
	}
	const std::string bezier_file = "./temp/" + model_name + ".xml";
	const std::string bdy_file = "./temp/" + model_name + ".bdy";
	IGA_CC_Helper::write2Bezier(bezier_file, solids);

	IGA::CCBezierBdy converted_boundary;
	IGA::BdyConvert_V2(
		_filename,
		hexbdy,
		converted_boundary,
		approximation);
	IGA_CC_Helper::write2bdy(
		bdy_file,
		converted_boundary,
		solids);

	std::map<int, Eigen::Vector3d> all_vertex;
	std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>> cell_to_vertex;
	for (int vertex_id = 0;
		vertex_id < static_cast<int>(mesh->vertices_.size());
		++vertex_id) {
		const auto vertex = mesh->vertices_[vertex_id];
		all_vertex[vertex_id] =
			Eigen::Vector3d{ vertex.x, vertex.y, vertex.z };
	}
	for (int cell_id = 0;
		cell_id < static_cast<int>(mesh->cells_.size());
		++cell_id) {
		cell_to_vertex[cell_id].resize(8);
		const auto cell_vertices = mesh->cells_[cell_id].vertex_;
		for (int local_vertex = 0; local_vertex < 8; ++local_vertex) {
			const int vertex_id = cell_vertices[local_vertex].idx();
			const auto vertex = mesh->vertices_[vertex_id];
			cell_to_vertex[cell_id][local_vertex] = std::make_pair(
				vertex_id,
				Eigen::Vector3d{ vertex.x, vertex.y, vertex.z });
		}
	}

	std::vector<std::vector<std::pair<int, double>>> neighbor_info;
	// 与顺序 MMC+SIMP 入口使用同一套单元邻域和权重构造。同步求解器
	// 的 FilterSensitivities 直接复用顺序材料阶段，因此不能只传空行。
	for (int i = 0;
		i < static_cast<int>(solids.Local_V3f.size());
		++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = std::numeric_limits<float>::max();
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (auto c_c_it = c_c.begin(); c_c_it != c_c.end(); ++c_c_it) {
			rmin = fmin(
				rmin,
				distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}
		rmin = rmin * 1.5;
		for (int k = 0; k < static_cast<int>(mesh->cells_.size()); ++k) {
			const double rdis = distance(
				mesh->cell_center[i],
				mesh->cell_center[k]);
			if (rdis <= rmin) {
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		std::vector<std::pair<int, double>> temp_neighbor_info;
		for (int j = 0;
			j < static_cast<int>(temp_aroundcell.size());
			++j) {
			temp_neighbor_info.emplace_back(std::pair<int, double>{
				temp_aroundcell[j],
				temp_aroundcell_weight[j] });
		}
		neighbor_info.emplace_back(temp_neighbor_info);
	}
	// 将界面参数显式装入同步求解器自己的参数对象，避免复用旧流程中的
	// 隐式默认值或阶段间状态。
	MMCMultiMaterialSimultaneousParams params;
	params.mmc_type = mmcType;
	params.mmc_target_rate = mmc_target_rate;
	params.max_iterations = std::max(1, static_cast<int>(max_iter));
	params.material_ratio_f = material_ratio_f;
	params.E1 = E1;
	params.E2 = E2;
	params.Poissons = Poisson;
	params.initial_eta = 1.0 / (1.0 + material_ratio_f);
	params.output_directory =
		"./mmc_multimaterial_simultaneous_result/" + model_name;

	MMCMultiMaterialSimultaneousResult result;
	const bool success = Cuda_TopOpt_MMC_MultiMaterial_Simultaneous(
		all_vertex,
		cell_to_vertex,
		bezier_file,
		bdy_file,
		pro,
		neighbor_info,
		params,
		result,
		density,
		is_need_update_topopt_model);
	std::cout
		<< "[MMC-MultiMaterial-Simultaneous] solver entry: "
		<< "Cuda_TopOpt_MMC_MultiMaterial_Simultaneous; "
		<< result.message
		<< std::endl;
	if (!success) {
		std::cerr
			<< "[MMC-MultiMaterial-Simultaneous] failed: "
			<< result.message
			<< std::endl;
	}
	delete mesh;
}
