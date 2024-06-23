#pragma once
#ifndef IGAME_SurfaceAlgorithm
#define IGAME_SurfaceAlgorithm

#include "SurfaceAlgorithm/Parameterization.h"
#include "SurfaceAlgorithm/Remeshing.h"
#include "SurfaceAlgorithm/Subdivision_Loop.h"
#include "SurfaceAlgorithm/Denoising_BNF.h"
#include "SurfaceAlgorithm/Simplification_QEM.h"
#include "SurfaceAlgorithm/Simplification_Powered_QEM.h"
#include "SurfaceAlgorithm/Deformation.h"
#include "SurfaceAlgorithm/Analysis_TriMesh.h"
#include "../Tools/TriangleKdTree.h"

// OpenMesh
#include "SurfaceAlgorithm/OpenMeshAlgorithm/OpenMesh_Smoothing.h"
#include "SurfaceAlgorithm/OpenMeshAlgorithm/OpenMesh_Remeshing.h"
#include "SurfaceAlgorithm/OpenMeshAlgorithm/OpenMesh_Simplification.h"
#include "SurfaceAlgorithm/OpenMeshAlgorithm/OpenMesh_Denoising.h"

#endif // !IGAME_SurfaceAlgorithm
