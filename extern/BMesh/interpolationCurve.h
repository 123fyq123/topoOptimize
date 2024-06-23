#ifndef __INTERPOLATIONCURVE_H__
#define __INTERPOLATIONCURVE_H__

#include <QMouseEvent>
#include <QGLWidget>
#include <Eigen/Dense>
#include <cmath>
#include <BMesh/point.hpp>
#include <OpenMesh/Core/IO/MeshIO.hh>
#include <OpenMesh/Core/Mesh/TriMesh_ArrayKernelT.hh>

typedef OpenMesh::TriMesh_ArrayKernelT<> MyMesh;
typedef struct SkelBoll {
	BMesh::Point3D p;
	double r;
	int index; // to curvePointss[i] index;
	int id;
};

//class Node {
//public:
//	Point3d p;
//	double r;
//	std::vector<int> children;
//	
//};
struct MyCmp {
	bool operator()(SkelBoll l, SkelBoll r) const { 
		if (fabs(l.p.x - r.p.x) > 1e-7) {
			return l.p.x < r.p.x;
		}
		else {
			if (fabs(l.p.y - r.p.y) > 1e-7)
				return l.p.y < r.p.y;
			else
			{
				if (fabs(l.p.z - r.p.z) > 1e-7) {
					return l.p.z < r.p.z;
				}
				else {
					return false;
				}
			}
		}
	}
};
class InterpolateCurve {
public:
	InterpolateCurve();
	~InterpolateCurve();
	double chordLength(const BMesh::Point3D& a, const BMesh::Point3D&b);
	void jsonReader(std::string);
	void genPoints(std::string);
	int findSpan(int numOfCtlPoint, int order, double u, const std::vector<double> &U);
	void basisFuns(int i, double u, int order, const std::vector<double> &U, std::vector<double> &base);
	void curvePoint(int n, int p, const std::vector<double>&U, std::vector<BMesh::Point3D>&P, double u, BMesh::Point3D&C);
	void resizeGL(int width, int height);
	void initializeGL();

	void globalCurveInterp(const std::vector<BMesh::Point3D> &Q, int p, int &m, std::vector<double> &U, std::vector<BMesh::Point3D> &ctlPts);

	void readOff();
	void readOff(std::string filename);
	void simpleSkel();
	void exportSkel();
protected:
	//void paintEvent(QPaintEvent *);
	void paintGL();//绘制opengl图形
	void mousePressEvent(QMouseEvent *event);
	void mouseReleaseEvent(QMouseEvent *event);
	void mouseMoveEvent(QMouseEvent *event);
	void wheelEvent(QWheelEvent *event);
public:
	//SweptVolumeConfig SVCfg;// 两条 nurbs 曲线
	std::string name;

	int p;
	std::vector<BMesh::Point3D> Q;//要拟合的点
	std::vector<BMesh::Point3D> curvePoints;//曲线上的点
	std::vector<BMesh::Point3D> ctlPts;//曲线的控制点
	int curvePointNum;
	//多条曲线
	std::vector<int> ps;
	std::vector<std::vector<BMesh::Point3D>> Qs;
	std::vector<std::vector<BMesh::Point3D>> curvePointss;
	std::vector<std::vector<BMesh::Point3D>> ctlPtss;
	std::vector<int> curvePointNums;
	// 对于多条曲线要每段骨架信息结构
	std::vector<std::vector<SkelBoll>> Qos;
	MyMesh mesh;// 三角网格
	int flag;// 判断选中的是哪一个控制点
	QPointF lastPos;// 因为是平面，所以只能获得（x,y)二维坐标
	double rotationX;
	double rotationY;
	double rotationZ;
	double times;
};

#endif