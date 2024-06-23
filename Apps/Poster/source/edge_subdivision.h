#pragma once
/* iGame Kernel */
#include "MeshKernel/Mesh.h"
#include "MeshKernel/IO.h"
#include "MeshKernel/MeshFileIO.h"

class EdgeSubdivision
{
public:
	void Spliting(MeshKernel::VolumeMesh* hexMesh, int eid, int split_num);

};

