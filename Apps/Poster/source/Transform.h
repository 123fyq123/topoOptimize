#pragma once

#include "Common/GeometryKernel.h"
#include<fstream>
class Transform {
public:
	Transform(HexV3fMesh *_HexModel): HexModel(_HexModel){};
	void Quad2Triangle();
	void Hex2IGEs();
public:
	HexV3fMesh* HexModel;
	std::string OutputPath;
};