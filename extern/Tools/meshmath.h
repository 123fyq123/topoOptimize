#pragma once
#include "MeshKernel/Mesh.h"
#include "CrossField/CrossField_Slover.h"

/* 将会在此继续更新一些关于网格的数学函数 */

namespace MeshMath {

	double get_volume_surface_mesh(MeshKernel::SurfaceMesh& mesh);

	double get_volume_tetahedral_mesh(MeshKernel::TetMesh& mesh);

	double get_volume_hexahedral_mesh(MeshKernel::VolumeMesh& mesh);

	void calc_cross_field(std::vector<Vec3>& points, std::vector<ID3>& triangles, std::vector<ID2>& feature_edges,
		std::vector<ID3>& singularities, std::vector<std::array<double, 9>>& global_triangle_dir);

	double get_min_scaled_Jacobian(MeshKernel::VolumeMesh& mesh, CH ch);

	double get_avg_scaled_Jacobian(MeshKernel::VolumeMesh& mesh, CH ch);
}

