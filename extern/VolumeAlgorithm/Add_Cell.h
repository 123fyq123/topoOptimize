#pragma once
#include "MeshKernel/Mesh.h"

#ifndef IGMAE_KERNEL_SIMPLIFY
#define IGMAE_KERNEL_SIMPLIFY
typedef MeshKernel::iGameVertex Vex;
typedef MeshKernel::iGameVertex Vec;
typedef MeshKernel::iGameVertexHandle VH;
typedef MeshKernel::iGameEdgeHandle EH;
typedef MeshKernel::iGameFaceHandle FH;
typedef MeshKernel::iGameCellHandle CH;
#endif

using namespace std;

void add_cell_2fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2);

void add_cell_3fh(MeshKernel::VolumeMesh& mesh, FH fh1, FH fh2, FH fh3);

void add_cell_eh_fh(MeshKernel::VolumeMesh& mesh, EH eh, FH fh);

void add_cell_fh(MeshKernel::VolumeMesh& mesh, FH fh);
