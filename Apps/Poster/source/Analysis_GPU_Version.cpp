//
// Created by XieJin on 10/25/2023.
//
#include "Analysis.h"
#include "Bdy.h"
#include "CatmullClarkBezierSolver.h"
#include "IGA_CC_Helper.h"
#include "Macro_Helper.h"
#include "SolverGPU/Cuda_Linear_Elastic_Solver.h"
#include "SolverGPU/Cuda_Hyperelasticity_Solver.h"
#include "SolverGPU/CPU_Linear_Elastic_Solver.h"
#include "SolverGPU/CPU_Hyperelasticity_Solver.h"
#include "SolverGPU/CPU_TimeIntegrator_Hyperelasticity_Solver.h"

#include <cstdlib>
#include <filesystem>
#include "Timer.h"


//void Analysis_V1(ModelDraw* widget,
//                 std::string _filename,
//                 IGA::HexBdy& hexbdy,
//                 std::vector<double>& x,
//                 std::vector<double>& y,
//                 std::vector<double>& z)
//{
//
//    if (hexbdy.empty()) {
//        std::cerr << "No Boundary Condition!" << std::endl;
//        return;
//    }
//
//
//    HexV3fMesh* mesh = new HexV3fMesh();
//    BSplineSolid m_solids;
//    FileHandler file(_filename);
//    widget->loadHexMesh(_filename);
//    std::string model_name = IGA_CC_Helper::get_model_name(_filename);
//    size_t patch_size = widget->hexMesh->CellSize();
//
//
//    //std::cout << "read2HexMeshObj" << std::endl;
//    file.read2HexMeshObj(mesh, true, true, true);
//    ApproximationCC AppCCSIMP(mesh, nullptr);
//    AppCCSIMP.construct_bss_new();
//    m_solids = AppCCSIMP.bss;
//
//    std::string bezier_file = "./temp/" + model_name + ".xml";
//    std::string bdy_file = "./temp/" + model_name + ".bdy";
//    std::string u_file = "./temp/" + model_name + ".u";
////    std::cout << bezier_file << std::endl;
////    std::cout << bdy_file << std::endl;
////    std::cout << u_file << std::endl;
//
//    if(!std::filesystem::exists(u_file)){
//        // Todo: sharpedge
//
//
//        IGA_CC_Helper::write2Bezier(bezier_file,m_solids);
//
//        IGA::CCBezierBdy ccbdy;
//        IGA::BdyConvert(_filename, hexbdy, ccbdy);
//        IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);
//        //std::cout << "BdyConvert end" << std::endl;
//
//        //std::cout << "cuda solver..." << std::endl;
//        {
//            //Timer timer;
//            Cuda_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
//        }
//
//        //system("C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\core.exe");
//
//    }
//
//    std::vector<std::vector<CPoint>> uu(patch_size, std::vector<CPoint>(64));
//    IGA_CC_Helper::readUU(u_file, uu);
//    /*CatmullClarkBezier model;
//    model.control_vertices.resize(m_solids.Global_size);
//
//    for (int patch_id = 0; patch_id < m_solids.Control_Point.size(); ++patch_id) {
//        std::vector<int> patch_global_ids;
//        for (int w_id = 0; w_id < 4; ++w_id) {
//            for (int v_id = 0; v_id < 4; ++v_id) {
//                for (int u_id = 0; u_id < 4; ++u_id) {
//                    int cur_local_id = m_solids.Control_Point[patch_id][w_id][v_id][u_id];
//                    int cur_global_id = m_solids.Local_to_Global[patch_id][cur_local_id];
//                    patch_global_ids.push_back(cur_global_id);
//                    assert(cur_global_id >= 0);
//                    double cx = m_solids.Local_V3f[patch_id][cur_local_id].x;
//                    double cy = m_solids.Local_V3f[patch_id][cur_local_id].y;
//                    double cz = m_solids.Local_V3f[patch_id][cur_local_id].z;
//                    model.control_vertices[cur_global_id] = Eigen::Vector3d(cx, cy, cz);
//                }
//            }
//        }
//        model.global_id.push_back(patch_global_ids);
//        assert(model.global_id.size() == 64);
//    }
//    std::cout << "model generate" << std::endl;
//    CatmullClarkBezierLinerElasticlSolver solver(model, ccbdy);
//    std::cout << "Solvering" << std::endl;
//    solver.Solving();
//
//    std::vector<std::vector<double>> m_bdycoef;
//    m_bdycoef.resize(m_solids.Local_V3f.size(), std::vector<double>(64*3, 0));
//
//    for (int patch_id = 0; patch_id < m_solids.Local_V3f.size(); patch_id++)
//    for (int cur_local_id = 0; cur_local_id < m_solids.Local_V3f[patch_id].size(); cur_local_id++)
//        for (int k = 0; k < 3; k++)
//            m_bdycoef[patch_id][3 * cur_local_id + k] = fabs(solver.m_u[m_solids.Local_to_Global[patch_id][cur_local_id] + k * m_solids.Global_size]);*/
//
//
//    /*std::vector<std::vector<double>> m_bdycoef;
//    m_bdycoef.resize(m_solids.Local_V3f.size(), std::vector<double>(64 * 3, 0));*/
//
//
//
//    std::vector<double> m_bdycoef_x, m_bdycoef_y, m_bdycoef_z;
//
//    m_bdycoef_x.resize(mesh->vertices_.size());
//    m_bdycoef_y.resize(mesh->vertices_.size());
//    m_bdycoef_z.resize(mesh->vertices_.size());
//
//    std::vector<int> mm = { 51, 3, 15, 63, 48, 0, 12, 60 };
//    for (int i = 0; i < m_solids.Local_V3f.size(); ++i)
//    {
//        for (int j = 0; j < 8; ++j)
//        {
//            auto vh = mesh->cells_[i].vertex_[j];
//
//            /*m_bdycoef_x[vh] = m_bdycoef[i][3 * j];
//            m_bdycoef_y[vh] = m_bdycoef[i][3 * j + 1];
//            m_bdycoef_z[vh] = m_bdycoef[i][3 * j + 2];*/
//            int jj = mm[j];
//            m_bdycoef_x[vh] = fabs(uu[i][jj].getX());
//            m_bdycoef_y[vh] = fabs(uu[i][jj].getY());
//            m_bdycoef_z[vh] = fabs(uu[i][jj].getZ());
////            if (i == 0) {
////                std::cout << "vh = " << vh << std::endl;
////                std::cout << "solution_z[0] = " << m_bdycoef_z[vh] << std::endl;
////                std::cout << "m_solids.Local_V3f[0] = " << m_solids.Local_V3f[0][j].x <<" " << m_solids.Local_V3f[0][j].y << " " << m_solids.Local_V3f[0][j].z << std::endl;
////            }
//
//        }
//
//    }
//    x = m_bdycoef_x;
//    y = m_bdycoef_y;
//    z = m_bdycoef_z;
//
//
//    // ????????????
//    //std::vector<std::vector<double>> local_bigrhs;
//    //std::vector<Matrix> local_matrix;
//    //std::vector<std::vector<double>> m_bdycoef;
//    //std::vector<BdyBase> bdybaseid;
//    //std::vector<std::vector<int>> m_colids;
//    //std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
//    //m_bdycoef.resize(m_solids.Local_V3f.size());
//
//    //std::vector<double> _material_distribution(m_solids.Local_V3f.size(), 1);
//
//
//    //for (int i = 0; i < m_solids.Local_V3f.size(); i++)
//    //{
//    //	m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);    //???????§Ö????????3
//    //}
//    //std::cout << "ccbdy.size() = " <<  ccbdy.size() << std::endl;
//    //for (auto& bdy : ccbdy)
//    //{
//    //	if (bdy.cell_id == 0) {
//    //		std::cout <<"bdy.cell_id = 0:" << bdy.global_base_id <<" " << bdy.local_base_id << " " << bdy.bdy_value << std::endl;
//    //	}
//    //	if (bdy.bdy_type == 2) {
//    //		for (int i = 0; i < 3; ++i)	m_bdycoef[bdy.cell_id][3 * bdy.local_base_id + i] = bdy.bdy_value[i];
//    //		//m_bdycoef[bdy.cell_id][3 * bdy.local_base_id + 1] = -1;
//    //	}
//    //
//    //	bdybaseid.push_back(BdyBase(bdy.local_base_id, bdy.global_base_id, bdy.cell_id, static_cast<int>(bdy.bdy_type)));
//    //}
//
//
//    //FunctionBSpline fb(m_solids, _material_distribution);
//
//    //std::vector<double> m_bigrhs;
//    //std::vector<int> rows, cols;          //???›¥????????????
//    //std::vector<double> elems;             //???›¥????????
//    //std::vector<double> m_coef;
//    //fb.Global_MatrixAndRhs_Sparse_liner_V1(m_bdycoef, bdybaseid, local_matrix,
//    //	m_bigrhs, local_bigrhs, rows, cols, elems, m_colids, m_bdyquick, m_coef);
//
//    /////*-----------------????????????-----------------*/
//    //
//
//    //for (int Local_V3f_id = 0; Local_V3f_id < m_solids.Local_V3f.size(); Local_V3f_id++)
//    //	for (int v_id = 0; v_id < m_solids.Local_V3f[Local_V3f_id].size(); v_id++)
//    //		for (int k = 0; k < 3; k++)
//    //			m_bdycoef[Local_V3f_id][3 * v_id + k] = fabs(m_coef[3 * m_solids.Local_to_Global[Local_V3f_id][v_id] + k]);
//
//    ////std::cout << m_coef << std::endl;
//
//    //std::vector<double> m_bdycoef_x, m_bdycoef_y, m_bdycoef_z;
//
//    //m_bdycoef_x.resize(mesh->vertices_.size());
//    //m_bdycoef_y.resize(mesh->vertices_.size());
//    //m_bdycoef_z.resize(mesh->vertices_.size());
//    //
//
//    //for (int i = 0; i < m_solids.Local_V3f.size(); ++i)
//    //{
//    //	for (int j = 0; j < 8; ++j)
//    //	{
//    //		auto vh = mesh->cells_[i].vertex_[j];
//
//    //		m_bdycoef_x[vh] = m_bdycoef[i][3 * j];
//    //		m_bdycoef_y[vh] = m_bdycoef[i][3 * j + 1];
//    //		m_bdycoef_z[vh] = m_bdycoef[i][3 * j + 2];
//    //		if (i == 0) {
//    //			std::cout << "vh = " << vh << std::endl;
//    //			std::cout << "solution_z[0] = " << m_bdycoef_z[vh] << std::endl;
//    //			std::cout << "m_solids.Local_V3f[0] = " << m_solids.Local_V3f[0][j].x <<" " << m_solids.Local_V3f[0][j].y << " " << m_solids.Local_V3f[0][j].z << std::endl;
//    //		}
//    //
//    //	}
//    //
//    //}
//    //x = m_bdycoef_x;
//    //y = m_bdycoef_y;
//    //z = m_bdycoef_z;
//    //
//    //return;
//}

void Analysis_V2(ModelDraw* widget,
                 std::string _filename,
                 IGA::HexBdy& hexbdy,
                 double Yangs,
                 double Poissons,
                 std::vector<double>& x,
                 std::vector<double>& y,
                 std::vector<double>& z,
                 bool bUsingGPU)
{

    if (hexbdy.empty()) {
        std::cerr << "No Boundary Condition!" << std::endl;
        return;
    }


    HexV3fMesh* mesh = new HexV3fMesh();
    BSplineSolid m_solids;
    FileHandler file(_filename);
    widget->loadHexMesh(_filename);
    std::string model_name = IGA_CC_Helper::get_model_name(_filename);
    size_t patch_size = widget->hexMesh->CellSize();


    std::cout << "read2HexMeshObj" << std::endl;
    file.read2HexMeshObj(mesh, true, true, true);
    ApproximationCC AppCCSIMP(mesh, nullptr);
    AppCCSIMP.construct_bss_new();
    m_solids = AppCCSIMP.bss;

    std::filesystem::path temp_path("./temp/");
    if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

    std::string bezier_file = "./temp/" + model_name + ".xml";
    std::string bdy_file = "./temp/" + model_name + ".bdy";
    std::string u_file = "./temp/" + model_name + ".u";
    std::cout << bezier_file << std::endl;
    std::cout << bdy_file << std::endl;
    std::cout << u_file << std::endl;

    //if(!std::filesystem::exists(u_file)) //for cache calculation results
    {
        // Todo: sharpedge

        IGA_CC_Helper::write2Bezier(bezier_file,m_solids);

        IGA::CCBezierBdy ccbdy;
        IGA::BdyConvert(_filename, hexbdy, ccbdy);
        IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);
        std::cout << "BdyConvert end" << std::endl;

        if(bUsingGPU)
        {
            std::cout << "cuda solver..." << std::endl;
            Timer timer;
            Cuda_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
        }
        else
        {
            std::cout << "cpu solver..." << std::endl;
            Timer timer;
            CPU_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
        }

        //system("C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\core.exe");

    }

    std::vector<std::vector<CPoint>> uu(patch_size, std::vector<CPoint>(64));
    IGA_CC_Helper::readUU(u_file, uu);

    std::vector<double> m_bdycoef_x, m_bdycoef_y, m_bdycoef_z;

    m_bdycoef_x.resize(mesh->vertices_.size());
    m_bdycoef_y.resize(mesh->vertices_.size());
    m_bdycoef_z.resize(mesh->vertices_.size());

    std::vector<int> mm = { 51, 3, 15, 63, 48, 0, 12, 60 };
    for (int i = 0; i < m_solids.Local_V3f.size(); ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            auto vh = mesh->cells_[i].vertex_[j];

            /*m_bdycoef_x[vh] = m_bdycoef[i][3 * j];
            m_bdycoef_y[vh] = m_bdycoef[i][3 * j + 1];
            m_bdycoef_z[vh] = m_bdycoef[i][3 * j + 2];*/
            int jj = mm[j];
            m_bdycoef_x[vh] = fabs(uu[i][jj].getX());
            m_bdycoef_y[vh] = fabs(uu[i][jj].getY());
            m_bdycoef_z[vh] = fabs(uu[i][jj].getZ());
            if (i == 0) {
                std::cout << "vh = " << vh << std::endl;
                std::cout << "solution_z[0] = " << m_bdycoef_z[vh] << std::endl;
                std::cout << "m_solids.Local_V3f[0] = " << m_solids.Local_V3f[0][j].x <<" " << m_solids.Local_V3f[0][j].y << " " << m_solids.Local_V3f[0][j].z << std::endl;
            }

        }

    }
    x = m_bdycoef_x;
    y = m_bdycoef_y;
    z = m_bdycoef_z;
}

void Analysis_V3(ModelDraw *widget, std::string _filename, IGA::HexBdy &hexbdy, double Yangs, double Poissons,
                 Solver SelectedSolver, vector<double> &x, vector<double> &y, vector<double> &z) {


    if (hexbdy.empty()) {
        std::cerr << "No Boundary Condition!" << std::endl;
        return;
    }


    HexV3fMesh* mesh = new HexV3fMesh();
    BSplineSolid m_solids;
    FileHandler file(_filename);
    widget->loadHexMesh(_filename);
    std::string model_name = IGA_CC_Helper::get_model_name(_filename);
    size_t patch_size = widget->hexMesh->CellSize();


    std::cout << "read2HexMeshObj" << std::endl;
    file.read2HexMeshObj(mesh, true, true, true);
    ApproximationCC AppCCSIMP(mesh, nullptr);
    AppCCSIMP.construct_bss_new();
    m_solids = AppCCSIMP.bss;

    std::filesystem::path temp_path("./temp/");
    if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

    std::string bezier_file = "./temp/" + model_name + ".xml";
    std::string bdy_file = "./temp/" + model_name + ".bdy";
    std::string u_file = "./temp/" + model_name + ".u";
    std::cout << bezier_file << std::endl;
    std::cout << bdy_file << std::endl;
    std::cout << u_file << std::endl;

    //if(!std::filesystem::exists(u_file)) //for cache calculation results
    {
        // Todo: sharpedge

        IGA_CC_Helper::write2Bezier(bezier_file,m_solids);

        IGA::CCBezierBdy ccbdy;
        IGA::BdyConvert(_filename, hexbdy, ccbdy);
        IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);
        std::cout << "BdyConvert end" << std::endl;

        switch(SelectedSolver)
        {
            case Solver::PCG_GPU :
            {
                std::cout << "PCG GPU solver..." << std::endl;
                Cuda_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
                break;
            }
            case Solver::PCG_CPU_Parallel :
            {
                std::cout << "PCG Parallel CPU solver..." << std::endl;
                CPU_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
                break;
            }
            case Solver::Pardiso_LDLT :
            {
                std::cout << "Pardiso LDLT solver..." << std::endl;
                CPU_Linear_Elastic_Solver_Paridiso(bezier_file, bdy_file, u_file, Yangs, Poissons);
                break;
            }
            case Solver::Eigen_LDLT :
            {
                CPU_Linear_Elastic_Solver_Eigen(bezier_file, bdy_file, u_file, Yangs, Poissons, false);
                break;
            }
            case Solver::Eigen_CG :
            {
                CPU_Linear_Elastic_Solver_Eigen(bezier_file, bdy_file, u_file, Yangs, Poissons, true);
                break;
            }
            case Solver::NO_Solver:
            {
                throw std::runtime_error("No Solver Selected");
                break;
            }
            default:
            {
                throw std::runtime_error("No Solver Selected");
                break;
            }
        }
        //system("C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\core.exe");

    }

    std::vector<std::vector<CPoint>> uu(patch_size, std::vector<CPoint>(64));
    IGA_CC_Helper::readUU(u_file, uu);

    std::vector<double> m_bdycoef_x, m_bdycoef_y, m_bdycoef_z;

    m_bdycoef_x.resize(mesh->vertices_.size());
    m_bdycoef_y.resize(mesh->vertices_.size());
    m_bdycoef_z.resize(mesh->vertices_.size());

    std::vector<int> mm = { 51, 3, 15, 63, 48, 0, 12, 60 };
    for (int i = 0; i < m_solids.Local_V3f.size(); ++i)
    {
        for (int j = 0; j < 8; ++j)
        {
            auto vh = mesh->cells_[i].vertex_[j];

            /*m_bdycoef_x[vh] = m_bdycoef[i][3 * j];
            m_bdycoef_y[vh] = m_bdycoef[i][3 * j + 1];
            m_bdycoef_z[vh] = m_bdycoef[i][3 * j + 2];*/
            int jj = mm[j];
            m_bdycoef_x[vh] = fabs(uu[i][jj].getX());
            m_bdycoef_y[vh] = fabs(uu[i][jj].getY());
            m_bdycoef_z[vh] = fabs(uu[i][jj].getZ());
            if (i == 0) {
                std::cout << "vh = " << vh << std::endl;
                std::cout << "solution_z[0] = " << m_bdycoef_z[vh] << std::endl;
                std::cout << "m_solids.Local_V3f[0] = " << m_solids.Local_V3f[0][j].x <<" " << m_solids.Local_V3f[0][j].y << " " << m_solids.Local_V3f[0][j].z << std::endl;
            }

        }

    }
    x = m_bdycoef_x;
    y = m_bdycoef_y;
    z = m_bdycoef_z;
}

void Hyperelasticity_Analysis_V1(ModelDraw* widget,
                                 std::string _filename,
                                 IGA::HexBdy& hexbdy,
                                 double Yangs,
                                 double Poissons,
                                 int frameCount,
                                 double timeStep,
                                 Solver SelectedSolver,
                                 std::vector<double> analysis_data[3],
                                 std::atomic_bool &is_need_update_model)
{
    if (hexbdy.empty()) {
        std::cerr << "No Boundary Condition!" << std::endl;
        return;
    }


    HexV3fMesh* mesh = new HexV3fMesh();
    BSplineSolid m_solids;
    FileHandler file(_filename);
    widget->loadHexMesh(_filename);
    std::string model_name = IGA_CC_Helper::get_model_name(_filename);
    size_t patch_size = widget->hexMesh->CellSize();


    std::cout << "read2HexMeshObj" << std::endl;
    file.read2HexMeshObj(mesh, true, true, true);
    ApproximationCC AppCCSIMP(mesh, nullptr);
    AppCCSIMP.construct_bss_new();
    m_solids = AppCCSIMP.bss;

    std::filesystem::path temp_path("./temp/");
    if (!std::filesystem::exists(temp_path)) std::filesystem::create_directory(temp_path);

    std::string bezier_file = "./temp/" + model_name + ".xml";
    std::string bdy_file = "./temp/" + model_name + ".bdy";
    std::string u_file = "./temp/" + model_name + ".u";
    std::cout << bezier_file << std::endl;
    std::cout << bdy_file << std::endl;
    std::cout << u_file << std::endl;

    //if(!std::filesystem::exists(u_file)) //for cache calculation results
    {
        IGA_CC_Helper::write2Bezier(bezier_file,m_solids);

        IGA::CCBezierBdy ccbdy;
        IGA::BdyConvert(_filename, hexbdy, ccbdy);
        IGA_CC_Helper::write2bdy(bdy_file, ccbdy, m_solids);
        std::cout << "BdyConvert end" << std::endl;

        switch(SelectedSolver)
        {
            case Solver::PCG_GPU :
            {
                std::cout << "PCG GPU solver..." << std::endl;
                Cuda_Hyperelasticity_Solver(bezier_file, bdy_file, u_file,
                                            Yangs, Poissons,
                                            frameCount, timeStep,
                                            mesh,
                                            analysis_data,
                                            is_need_update_model);
                break;
            }
            case Solver::PCG_CPU_Parallel :
            {
                std::cout << "PCG Parallel CPU solver..." << std::endl;
                CPU_Linear_Elastic_Solver(bezier_file, bdy_file, u_file, Yangs, Poissons);
                break;
            }
            case Solver::Pardiso_LDLT :
            {
                std::cout << "Pardiso LDLT solver..." << std::endl;
                CPU_Hyperelasticity_Pardiso_Solver(bezier_file, bdy_file, u_file,
                                                   Yangs, Poissons,
                                                   frameCount, timeStep,
                                                   mesh,
                                                   analysis_data,
                                                   is_need_update_model);
                break;
            }
            case Solver::Eigen_LDLT :
            {
                CPU_Hyperelasticity_Pardiso_Solver(bezier_file, bdy_file, u_file,
                                                   Yangs, Poissons,
                                                   frameCount, timeStep,
                                                   mesh,
                                                   analysis_data,
                                                   is_need_update_model);
                break;
            }
            case Solver::Eigen_CG :
            {
                CPU_Linear_Elastic_Solver_Eigen(bezier_file, bdy_file, u_file, Yangs, Poissons, true);
                break;
            }
            case Solver::Debug :
            {
                CPU_Hyperelasticity_TimeIntegrator_Solver(bezier_file, bdy_file, u_file,
                                                   Yangs, Poissons,
                                                   frameCount, timeStep,
                                                   mesh,
                                                   analysis_data,
                                                   is_need_update_model);
                break;
            }
            case Solver::NO_Solver:
            {
                throw std::runtime_error("No Solver Selected");
                break;
            }
            default:
            {
                throw std::runtime_error("No Solver Selected");
                break;
            }
        }
    }
}