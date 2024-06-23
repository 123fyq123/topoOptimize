#pragma once

#include "Common/GeometryKernel.h"

class  SurfaceFitter
{
public:
	SurfaceFitter(HexV3fMesh* _SkelMesh, HexV3fMesh* _OrigMesh) :SkelMesh(_SkelMesh), OrigMesh(_OrigMesh) {};
	void CCSubdivisionFitting();
	HexV3fMesh* OrigMesh;
	HexV3fMesh* SkelMesh;
};
