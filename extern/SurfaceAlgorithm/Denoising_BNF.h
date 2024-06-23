#pragma once
#include "MeshKernel/Mesh.h"
#include <Eigen/Dense>
#include <Eigen/Core>
#include <Eigen/Sparse>

// Author: Zhang Sheng
// Time: 2021.8.26
// 
//MeshKernel::SurfaceMesh mesh = io.ReadFile("bunny_random.obj");
//Denoising_BNF test(mesh, 0.5, 5, 15);
//test.Execute();
//io.WriteFile(mesh, "bunny_denoise.obj");


class Denoising_BNF {
public:
	Denoising_BNF(MeshKernel::SurfaceMesh& _mesh, float _SigmaNormal = 0.5f, int it_normal_filtering = 5, int it_vertex_updating = 15) :
		mesh(_mesh), SigmaNormal(_SigmaNormal), it_normal(it_normal_filtering), it_vertex(it_vertex_updating) {
	};
	void Execute();

private:

	MeshKernel::SurfaceMesh& mesh;

	// 算法参数
	float SigmaNormal, SigmaCenter = 0.f;
	int it_normal, it_vertex;

	// 三角形数据
	std::vector<Eigen::Vector3f> centers;
	std::vector<Eigen::Vector3f> normals;
	std::vector<float> areas;

	// 顶点数据
	std::vector<Eigen::Vector3f> positions;

	void initTrianglesData();
	void initSigmaCenter();
	void initPositions();
	void BilateralNormalFiltering();
	void VertexUpdate();
};




