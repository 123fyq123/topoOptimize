#pragma once
#include "MeshKernel/Mesh.h"

#ifndef D_PI
#define D_PI 6.2831852f
#endif // !D_PI

#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

class Analysis_HexMesh {
private:
	MeshKernel::VolumeMesh& mesh;

public:
	Analysis_HexMesh(MeshKernel::VolumeMesh& _mesh) : mesh(_mesh) {

	}

	void calculateCellQuality_Jacobian(std::vector<double>& scaled_jacobian);// the value of Scaled Jacobian is the	determinant of the matrix Jacobian

};
