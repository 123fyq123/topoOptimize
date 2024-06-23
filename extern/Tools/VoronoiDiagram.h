#pragma once
#include <vector>
#include <fstream>
#include <iostream>
#include <set>

#include "TinyVector.h"
#include "MeshKernel/Mesh.h"
#include "SkeletalMesh/SkeletalMesh.h"

#undef  Handle(Class)
#include <CGAL/Handle.h>
#define Handle(Class) opencascade::handle<Class>
//#define Handle(Class) CGAL::handle<Class>
#include <Standard_Handle.hxx>

#include <CGAL/Weighted_point_3.h>
#include <CGAL/bounding_box.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Regular_triangulation_3.h>


// 计算 Voronoi Diagram 的顶点并返回

using namespace std;

namespace VoronoiDiagram {

	void get_vertices_of_voronoi_diagram(const vector<Vector3d>& sample_points, vector<Vector3d>& vertices_of_voronoi_diagram);

	void get_medial_mesh(const vector<Vector3d>& sample_points, SkeletalMesh& medial_mesh);

};
