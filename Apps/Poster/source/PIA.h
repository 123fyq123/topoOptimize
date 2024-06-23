#pragma once
#include "Common/GeometryKernel.h"
#include <iostream>
#include<unordered_map>

class PIA_FOR_CC_SOLID {

public:
	HexV3fMesh *mesh;
	HexV3fMesh *firstmesh;
	HexV3fMesh *oldmesh;
	double err_now;
	PIA_FOR_CC_SOLID(HexV3fMesh *_firstmesh,HexV3fMesh *_mesh) :firstmesh(_firstmesh),mesh(_mesh),err_now(0.0) {};
public:
	void PIA_FOR_CC_SOLID_ONCE();
	void PIA_FOR_HIGH_ORDERED_MESH();

};

class PartPIA { 
public:
	HexV3fMesh* mesh;                                       // 操作的网格；
	std::vector<VertexHandle> MovablePart;                  // 可移动的部分；
	std::vector<V3f> TargetPoint;                           // 需要插值到的点
	std::vector<V3f> MeshPoint;
	std::vector<EdgeHandle> SharpEdge;                      // 尖锐边
	double CurrentError;                                    // 当前误差
	PartPIA(HexV3fMesh* _mesh, std::vector<V3f>& _TargetPoint, std::vector<V3f>& _MeshPoint, std::vector<EdgeHandle>& _SharpEdge) ;
	void GetTargetPoint(std::vector<V3f>& _TargetPoint);
	void PartUpDate();
	void GetMovablePart();
};


class SurfacePIA {
public:
	HexV3fMesh* mesh;                                       // 操作的网格；
	std::unordered_map<int, int>FixedPart;                  // 固定的部分；
	std::vector<VertexHandle> MovablePart;                  // 可移动的点
	std::vector<V3f> TargetPoint;                           // 需要插值到的点
	std::vector<EdgeHandle> SharpEdge;                      // 尖锐边
	std::unordered_map<int, int> SharpPoints;               // 尖锐边上的点，使用哈希映射便于查找
	double CurrentError;                                    // 当前误差
	SurfacePIA(HexV3fMesh* _mesh, std::vector<V3f>& _TargetPoint, std::vector<EdgeHandle>& _FixedEdge);
	void GetFixedPart(std::vector<EdgeHandle>& _FixedEdge);
	void GetSharpEdge(std::vector<EdgeHandle>& _SharpEdge);
	void GetSharpPoints();
	void GetMovablePart();
	void GetTargetPoint(std::vector<V3f>& _TargetPoint);
	void SurfaceUpDate();
};



