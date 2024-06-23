#pragma once
#include "MeshKernel/Mesh.h"
#include <Eigen/Core>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseCore>
#include <Eigen/SparseCholesky>
#include <unordered_map>
#include <omp.h>
//#pragma omp parallel for

using namespace std;
using namespace MeshKernel;

struct transform_node {
	iGameCellHandle ch;
	Eigen::Vector3d pos;
	transform_node(): ch(-1), pos(0, 0, 0){}
	transform_node(iGameCellHandle _ch, Eigen::Vector3d _pos): ch(_ch), pos(_pos){}
};

class HexMesh_Regularization {
public:

	HexMesh_Regularization(VolumeMesh& _mesh, bool _use_exterior_only ): mesh(_mesh), use_exterior_only(_use_exterior_only){

	}

	void regularization(uint32_t iter);

private:

	VolumeMesh& mesh;
	bool use_exterior_only = false;
	bool geometric_constraint_flag = true;
	bool check_weight_flag = false;
	bool check_coordinate_flag = false;
	double scaling_factor = 0.7f;// scaling factor, should > 0
	double weight_e = 5;// Local Regularization use, the weight of edge direction( >= 5 is ok)
	double weight_d = 5;// Global Optimization use, the weight of node displacement, default: 5
	double weight_g = 10;// Global Optimization use, the weight of geometric constraints, it is much larger than weight_d, default: 100
	double factor_e = 1.f / 6;
	double KEEPOPSITIVE = 1E-8F;
	double in_min_scaled_jacobian = 2.f, in_avg_scaled_jacobian = 0.f;
	double out_min_scaled_jacobian = 2.f, out_avg_scaled_jacobian = 0.f;

	unordered_map<int, vector<transform_node>> transforms;// save the temporary transformed element nodes
	unordered_map<int, double> weight_cell;// Global Optimization use, shoule be opsitive
	unordered_map<int, double> weight_edge;// Global Optimization use, shoule be opsitive
	unordered_map<int, Eigen::Vector3d> normal_vertex;// vh and its normal
	unordered_map<int, Eigen::Vector3d> ideal_position;// vh and its ideal position

	//vector<vector<int>> hexFaces = {// 6 rows, 4 cols
	//	{1, 2, 3, 4}, {1, 5, 6, 2}, {2, 6, 7, 3},
	//	{3, 7, 8, 4}, {4, 8, 5, 1}, {5, 8, 7, 6}
	//};

	//vector<vector<int>> dualFaces = {// 8 rows, 3 cols
	//	{1, 2, 5}, {1, 3, 2}, {1, 4, 3}, {1, 5, 4},
	//	{6, 5, 2}, {6, 2, 3}, {6, 3, 4}, {6, 4, 5}
	//};

	//vector<vector<int>> hexEdges = {
	//	{1, 2}, {2, 3}, {3, 4}, {4, 1}, {1, 5}, {2, 6},
	//	{3, 7}, {4, 8}, {5, 6}, {6, 7}, {7, 8}, {8, 5}
	//};

	vector<vector<int>> tetFaces = {
		{1, 0, 2, 5}, {3, 2, 0, 7},
		{4, 0, 5, 7}, {6, 7, 5, 2}
	};

	void CalculateCellQuality(vector<iGameCellHandle>& chs, int mode);// mode 0 estemm initial mesh, mode 1 esteem result mesh
	void Preprocessing(vector<iGameEdgeHandle>& ehs, vector<iGameCellHandle>& chs);
	void LocalRegularization(const iGameCellHandle& ch);
	void GlobalOptimization();
	void ComputeTetrahedronLaplacianWeight(vector<iGameVertexHandle> vhs);//  the weights should be positive

	
	void GlobalOptimization_Exterior();
	void LaplacianSmoothing_Interior();

};
