#ifndef __SOLUTIONSKEL_H__
#define __SOLUTIONSKEL_H__
#include <string>
#include "BMesh/bmesh_solution.hpp"
#include <QObject>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include <glm/gtx/euler_angles.hpp>
#include <glm/gtx/norm.hpp>
#include "MeshKernel/Mesh.h"


typedef enum Estate{
	haveNotGenJN = 0,
	haveGenJN,
}STATE;

class SolutionSkel : public QObject {
	Q_OBJECT

public:
	SolutionSkel();
	~SolutionSkel();
public:
	bool genNodeHex();
	void genViewHex(HexV3fMesh &mesh);
	void genViewHex(MeshKernel::VolumeMesh& mesh);
	void setMesh(HexV3fMesh *mesh_);
	BMesh::Node *selectWhichNode(V3f &p);
	BMesh::Node* selectWhichNode(const V3f& p);
	BMesh::Node *selectWhichNode_(BMesh::Point3D &p);
	BMesh::Node *getRoot();
	std::vector<BMesh::Node*> getNodeVec();
	void rotateCube(int direction, int plus);
	void moveCube(int direction, int plus);
	void clearSkel();
	bool readSkel(std::string filename);
	void saveSkel(std::string name);
	void clearJNCube();
	void getALine(std::vector<BMesh::Node*>& ret);

	bool changeNodeType();
	bool changeSnap(const BMesh::Point3D&pOld, const BMesh::Point3D &pNew, BMesh::Point3D &out);
	double bigStringToZero(double p, int count);

public slots:
	void rscaleNode(double deltaR, int type);
	void rdeltaEditNode(double x, double y, double z, int type);
	void reDrawCubeR();
	void changeToPolylines();
	void simpleWithCgal();
	void pipeline2();


	/*
	void rsetStepDouble(double);
	*/
	void pipeline();
	void changeAngleAuto();
	void generateLinkBody();
	void topologicalSplit();
	void refineJN();

public:
	BMesh::Node * selectNode;
	std::vector<BMesh::Node *> selectNodes;
	BMesh::Node* mirrorNode;
	BMesh::Point3D positionForSnif;
	BMesh::Point3D snifDelta;
	HexV3fMesh *mesh;
	STATE state;
	bool bdragMode;
	bool bshowTriMeshChanged;
	bool bRotateByMouse;
	bool bVertexSnap;
	bool bmirrorMode;

public:
	BMesh::Solution * solutionOldBMesh;
	double step;

	/*******************/
public:
	void reditNode(double x, double y, double z, double angleX, double angleY, double angleZ, double length, double width, double height);
	void raddNode(double x, double y, double z, double angleX, double angleY, double angleZ, double length, double width, double height);
	void rcopyNode();
	void rinsertNode(BMesh::Node* pre);
	void rdelNode();
	void rXyzCubeByMouse(double x, double y, double z, double length, double width, double height);
	void rRotateCubeByMouse(V3f axis, double angle);
	void spiltX();
	void spiltY();
	void spiltZ();
	bool smartConnectNodes(BMesh::Node* a, BMesh::Node* b);
	bool connectNodes(BMesh::Node* a, BMesh::Node* b);
	bool disConnectNodes(BMesh::Node* a, BMesh::Node* b);

private:
	glm::quat SolutionSkel::inAngleAxis_(glm::vec3 RotationAxis, double RotationAngle);
	void xyzCubeByMouse(double x, double y, double z, double length, double width, double height);
	void rotateCubeByMouse(V3f, double);

};


#endif