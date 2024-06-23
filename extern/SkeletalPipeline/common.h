#ifndef COMMON_H__
#define COMMON_H__

#include <iostream>
#include <fstream>
#include <iomanip>
#include <cassert>
#include <string>
#define _USE_MATH_DEFINES
#include <cmath>
#include<numeric>

#include <omp.h>

#include <chrono>

#include <algorithm>
#include <iterator>
#include <vector>
#include <map>
#include <unordered_map>
#include <set>
#include <queue>

#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/Eigenvalues>

#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/PolyMesh_ArrayKernelT.hh>
#include <OpenMesh/Core/Utils/PropertyManager.hh>

typedef OpenMesh::PolyMesh_ArrayKernelT<> polyMesh; 


struct Comp_v2d
{
    Eigen::Vector2d ori;
    Comp_v2d(Eigen::Vector2d x2):ori(x2){}
    bool operator()(const Eigen::Vector2d &a , const Eigen::Vector2d &b)
    {
        return a.dot(ori)<b.dot(ori);   
    }

};

//#define GUI_ON
//#define GUI_OPT_ON
//#define GUI_CSA_ON

#endif // !COMMON_H__
