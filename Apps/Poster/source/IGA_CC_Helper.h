//
// Created by XieJin on 10/25/2023.
//

#pragma once

#include <cmath>
#include "Bspline3D/BSplineSolid.h"
#include "FileHandler.h"
#include "BezierSolid.h"
#include <iostream>
#include <vector>
#include <string>
#include "Geom/CPoint.h"
#include <windows.h>
#include "Bdy.h"
#include "CatmullClarkBezierSolver.h"
#include <MeshKernel/Mesh.h>
#include <Geom/CBSplineSurface.h>


class IGA_CC_Helper {
public:

    static vector<CBSplineSurface> muti_bezier_batch_;
    static vector<CBSplineSurface> muti_bezier_batch_dir_;

    static void write2Bezier(const std::string &_filepath, BSplineSolid &m_solids);

    static void write2bdy(const std::string &_filepath, IGA::CCBezierBdy &ccbdy, BSplineSolid &m_solids);

    static void readUU(const std::string &_filename, std::vector<std::vector<CPoint>> &uu);

    static void write2prot(const std::string &_filename, const std::vector<int> &protect_ids);

    static void write2neighbor(const std::string &_filename,
                               const std::vector<std::vector<CellHandle>> &aroundcell,
                               const std::vector<std::vector<double>> &aroundcell_weight);

    static std::string get_model_name(std::string &file_path);

    static bool read2BezierOnly(const std::string &_filename);

    static bool read2Bezier(const std::string &_filename);

    static bool
    cad_restruction(vector<vector<vector<Vec>>> &control_points, vector<vector<vector<Vec>>> &control_points_dir);

    static bool cad_restruction(vector<vector<vector<Vec>>> &control_points);

};



