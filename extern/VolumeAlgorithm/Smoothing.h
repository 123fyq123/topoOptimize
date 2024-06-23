#pragma once
#include "MeshKernel/Mesh.h"
#include "Tools/meshmath.h"
#include <Eigen/Core>
#include <queue>
#include "omp.h"

using namespace std;
using namespace MeshKernel;

class HexMesh_Smoothing {
public:

	HexMesh_Smoothing(MeshKernel::VolumeMesh& _mesh);
	void Laplacian_Smoothing(int max_iter = 1, std::vector<int> feature_cells = {}, std::vector<int> feature_faces = {}, 
		std::vector<int> feature_edges = {}, std::vector<int> feature_vertices = {});
	void Laplacian_Smoothing_Diffusion(int max_iter = 1, std::vector<int> feature_cells = {}, std::vector<int> feature_faces = {},
		std::vector<int> feature_edges = {}, std::vector<int> feature_vertices = {});
	void Laplacian_Smoothing_HC(std::vector<int> feature_vhs, int iter_max = 5);

private:
	MeshKernel::VolumeMesh& mesh;
	double threshold_J = 0.9;
	vector<iGameVertexHandle> interior_vhs;
	vector<iGameVertexHandle> exterior_vhs;
	double CalculateCellQuality(MeshKernel::iGameCellHandle ch);

};
