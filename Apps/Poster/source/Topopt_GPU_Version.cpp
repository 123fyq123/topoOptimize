//
// Created by XieJin on 10/27/2023.
//

#include <filesystem>
#include"TopOpt.h"
#include "IGA_CC_Helper.h"
#include "CatmullClarkBezierSolver.h"
#include "Timer.h"
#include "SolverGPU/Cuda_TopOpt_Pipeline.h"

void TopOpt_GPU_V0(      std::string _filename,
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


    IGA_CC_Helper::write2Bezier(bezier_file,m_solids);    // xml

    IGA::CCBezierBdy ccbdy;
    IGA::BdyConvert(_filename, hexbdy, ccbdy);
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
        for(int j = 0; j < temp_aroundcell.size(); ++j){
            temp_neighbor_info.emplace_back(std::pair<int, double>{temp_aroundcell[j], temp_aroundcell_weight[j]});
        }
        neighbor_info.emplace_back(temp_neighbor_info);
        aroundcell.push_back(temp_aroundcell);
        aroundcell_weight.push_back(temp_aroundcell_weight);
        //std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
    }



    {
        std::cout <<"protected_cell size = " << pro.size() <<"\n";
        Cuda_TopOpt_Pipeline(bezier_file, bdy_file, pro, neighbor_info,
                              target_rate, max_iter, Yangs, Poisson, density, is_need_update_topopt_model);
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
    IGA::BdyConvert(_filename, hexbdy, ccbdy);
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
        Cuda_TopOpt_Discrete(bezier_file, bdy_file, pro, neighbor_info,
            target_rate, max_iter, Yangs, Poisson, density, is_need_update_topopt_model);
    }


}