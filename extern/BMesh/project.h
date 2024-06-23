/*
 * @Author: lishaohsuai 1049188593@qq.com
 * @Date: 2022-05-10 16:53:02
 * @LastEditors: lishaohsuai 1049188593@qq.com
 * @LastEditTime: 2022-05-10 23:19:59
 * @FilePath: \ToPoster\BMesh\project.h
 * @Description: 这是默认设置,请设置`customMade`, 打开koroFileHeader查看配置 进行设置: https://github.com/OBKoro1/koro1FileHeader/wiki/%E9%85%8D%E7%BD%AE
 */
#pragma once
#include <iostream>
#include <string>
#include "Common/GeometryKernel.h"
#include "Meshkernel/Mesh.h"
#include <vector>
#include "triangle.hpp"
#include "geometry.hpp"
#include <limits>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>
#ifdef _DEBUG
#pragma comment(lib, "OpenMeshCored.lib")//�������ӿ��ķ�ʽ�ƺ�����
#pragma comment(lib, "OpenMeshToolsd.lib")
#else
#pragma comment(lib, "OpenMeshCore.lib")
#pragma comment(lib, "OpenMeshTools.lib")
#endif

#include "Tools/BVH.h"

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
class Project{
public:
	std::string triFile;
	std::vector<BMesh::Triangle> tv;
	MyMesh triMeshtv;
	std::vector<BMesh::Lvector> triangleVector;
	std::vector<BMesh::Point3D> points;
	//std::set<int> fixPoints;
	HexV3fMesh *hMesh;
	MeshKernel::VolumeMesh* hexMesh;
	BMesh::Point3D center;
	double radius;
	BMesh::Point3D boxMinPoint;
	BMesh::Point3D boxMaxPoint;
	std::vector<std::pair<V3f, V3f>> projectLines;
	std::vector<std::pair<V3f, V3f>> projectLines2;
	MyMesh triMesh;
	std::map<int, int> m;// �����ı��α߽��Ķ���idӳ��Ϊ�µ�����������id
	std::map<int, int> mm; // ������������idӳ��Ϊ�ı���������id

	//Project(std::string triFile_, HexV3fMesh* hMesh_) {
	//	triFile = triFile_;
	//	hMesh = hMesh_;
	//	readTriMesh(triFile);
	//}
public:
	void readTriMesh(std::string triFile);
	void projectVolumeToOff(); 
	void projectVolumeToOff_sec();
	void projectVolumeToOff_sec1();
	void projectVolumeToOffFindClosePoint(); // not to do really
	void dense();// ������������
	void movePoint();
	void movePoint1();
	bool surfaceQuad2Tri();
	bool surfaceQuad2Tri1();
	

public:
	// ��ͼ����
	void calCenter();
};