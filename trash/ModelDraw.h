#pragma once
/* 传入模型对象，处理模型 */
#include "ViewOpenGL.h"
/* std */
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <vector>
#include <cstdio>
/* Qt */
#include <QtCore/QString>
#include <QtCore/qstring.h>
#include <QtGUI/QtGui>
#include <QtWidgets/qfiledialog.h>
#include <QtCore/QtCore>
#include <QtWidgets/QLayout>
#include <QtWidgets/QMessageBox>
#include <QtWidgets/qprogressdialog.h>
#include <QtWidgets/QApplication>
#include <QtWidgets/qinputdialog.h>
#include <QtGUI/Qsurfaceformat.h>
#include <QtGUI/QPixmap>
#include <QtGUI/QPainter>
/* Eigen */
#include <Eigen/Dense>
#include <Eigen/Core>

/* iGame Kernel */
#include "Kernel/Mesh.h"
#include "Kernel/IO.h"
#include "Kernel/MeshFileIO.h"
#include "Algorithm/MeshAlgorithm.h"
/*    Tools    */
#include "Tools/BVH.h"
#include "Tools/TinyVector.h"
#include "Tools/ColorBar.h"
/*    SKEL     */
#include "BMesh/simpleWithCgal.h"
#include "BMesh/project.h"
#include "solutionSkel.h"

/* Geom */
#include "Geom/CGeom.h"
#include "Geom/CGeomSurface.h"
#include "Geom/CPoint.h"
#include "Geom/CLine.h"
#include "Geom/CPlane.h"
#include "Geom/CBezierCurve.h"
#include "Geom/CBezierSurface.h"
#include "Geom/CCylindricalSurface.h"
#include "Geom/CCircle.h"
#include "Geom/CSphere.h"
#include "Geom/SurfaceOfLinearExtrusion.h"
#include "Geom/SurfaceOfRevolution.h"
#include "Geom/CAidFunction.h"
/* Model Reconstruction */
#include "ModelRestruct/BezierSolid.h"
#include "ModelRestruct/FileHandler.h"
#include "ModelRestruct/Subdivision.h"
#include "ModelRestruct/Transform.h"
#include "ModelRestruct/GenHexSkel.h"
#include "ModelRestruct/SurfaceFitter.h"

/* MMC 组件转换 */
#include "Pipeline/include/Pipeline.h"

#include "BMesh/mlca_subdivision.h"

/*    otherz m yes    */
#include "NewFilter/ProDraw.h"
#include "HexMeshTool.h"
#include "CoordinateAxis.h"
#include "edge_subdivision.h"
#include "NewFilter/DotLineCalculation.h"

#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Point_set_3.h>
#include <CGAL/poisson_surface_reconstruction.h>
#include <CGAL/mst_orient_normals.h>
#include <CGAL/Surface_mesh/Surface_mesh.h>
#include <CGAL/remove_outliers.h>
#include <CGAL/grid_simplify_point_set.h>
#include <CGAL/jet_smooth_point_set.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/Scale_space_reconstruction_3/Jet_smoother.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Polygon_mesh_processing/border.h>
#include <boost/iterator/function_output_iterator.hpp>
#include <CGAL/Polygon_mesh_processing/detect_features.h>
#include <CGAL/Polygon_mesh_processing/smooth_mesh.h>
#include <CGAL/Polygon_mesh_processing/detect_features.h>



using namespace std;

// types
typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
typedef Kernel::FT FT;
typedef Kernel::Point_3 Point_3;
typedef Kernel::Vector_3 Vector_3;
typedef Kernel::Sphere_3 Sphere_3;
typedef CGAL::Point_set_3<Point_3, Vector_3> Point_set;
typedef CGAL::Surface_mesh<Kernel::Point_3>   MeshCGAL;

typedef boost::graph_traits<MeshCGAL>::halfedge_descriptor        halfedge_descriptor;
typedef boost::graph_traits<MeshCGAL>::edge_descriptor            edge_descriptor;

namespace PMP = CGAL::Polygon_mesh_processing;

class ModelDraw :public ViewOpenGL {

	Q_OBJECT

public:
	ModelDraw(QWidget* parent = 0);
	virtual ~ModelDraw(void);

	/* 输出信息 */
	void setPickInfo(std::function<void(void)> _func) 
	{
		showPickInfo = _func;
	}
	std::string pick_info;

	bool compare_flag = false;
	std::vector<GLfloat> compared_vertices_data;// Just draw points
	std::unordered_map<int, int> compared_to_current;
	void compared_bind_vertices();
	int used_last_vcnt = 0;
	std::vector<GLfloat> translucence_vertices_data;
	bool show_mesh_quality = false;



signals:
	void Send(QString text); // clicked;
public slots:
	void Receive(QString text); 
	void Receive();

protected:

	std::function<void(void)> showPickInfo;

	GLuint TestVAO = 0, TestVBO = 0, TestaEBO = 0;
	GLuint VAO = 0, VBO = 0, EBO = 0;// 模型数据
	vector<GLuint> VAOs;
	vector<GLuint> VBOs;
	vector<GLuint> EBOs;
	vector<GLuint> VAOs_Cylinder;
	vector<GLuint> VBOs_Cylinder;
	GLuint VAO_Possible_Shpere;
	GLuint VBO_Possible_Shpere;
	MeshKernel::CurveSkel* curve_skel_ForDraw;
	std::unordered_map<int, int> IntVertexHandleToSphereMesh;
	std::unordered_map<int, int> SphereMeshToIntVertexHand;
	std::unordered_map<int, int> IntEdgeHandleToCylinderMesh; 
	std::unordered_map<int, int> CylinderMeshToIntEdgeHandle;


	std::vector<float> Data_Possible_Sphere;

	// Test绘制数据
	std::vector<GLfloat> Testvertices_data;

	// 绘制数据
	std::vector<GLfloat> vertices_data;
	
	size_t tri_cnt = 0, quad_cnt = 0, line_cnt = 0, translucence_tricnt = 0;// 记录绘制数据信息
	

	// 面网格
	MeshKernel::SurfaceMesh* triMesh = nullptr;// 三角形网格
	MeshKernel::SurfaceMesh* quadMesh = nullptr;// 四边形网格
	
	// 体网格
	MeshKernel::VolumeMesh* hexMesh = nullptr;// 六面体网格
	MeshKernel::VolumeMesh* tetMesh = nullptr;// 四面体网格

	//曲线骨架
	MeshKernel::CurveSkel* curveSkel = nullptr; // 曲线骨架

	

	// 拓扑优化约束数据
	std::unordered_map<int, bool> force_id;// 约束顶点ID
	std::unordered_map<int, bool> spc_id;
	std::unordered_map<int, Vex> constrained_point3;

	// 稀疏优化数据
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> axis;
	std::vector<std::vector<Eigen::Vector2d>> bdy;
	polyMesh* poly_mesh = nullptr;
	Pipeline* pipeline = nullptr;

	// 点选数据
	BVH_Tree* bvh_tree = new BVH_Tree();
	bool mouse_move_flag = false;
	std::vector<int> selected_v_indices;
	std::vector<int> selected_e_indices;
	std::vector<int> selected_f_indices;
	std::vector<int> selected_c_indices;
	std::vector<int> selected_v_draw_indices;
	std::vector<QVector3D> selected_add_possible_point;

	std::vector<int> SelectedPossibleVhIndices;
	QVector3D add_possible_point; // 传到绘图中画出来
	
	bool is_change_radius = false;

	double selected_v_radius;

	
	MeshKernel::iGameVertexHandle copy_point; // 是用于复制的节点的handle
	//MeshKernel::iGameVertexHandle copy_point_draw; // 是用于绘制的 点
	MeshKernel::iGameVertexHandle MoveSelectedPoint;
	
	MeshKernel::iGameVertexHandle delet_point_vh;


	//bool Can_send_VAO = false;
	bool is_edit_point = false;

	bool is_line_draw = true;
	bool is_possible_add_point_draw = true;
	bool is_copy_point_draw = false;

	// 选择框数据
	vector<QVector3D> selected_rectangle;
	ColorBar colorbar;


public:

	// wu_yuan
	std::vector<MeshKernel::iGameVertexHandle> Move_Point_vh;

	void TestloadShader();
	void loadShader();

	void ModelClear();// 清除所有模型数据
	void updateDrawData(bool need_update_bvhtree = true);// 更新绘制数据
	void updateDrawData_CutPlane();
	void updateMeshModel();// 根据 bbox 设置中心
	void updateBVH();

	int getVerticesSize();

	
	/* 曲面网格函数 */
	bool loadSurfaceMesh(std::string filename);
	bool saveSurfaceMesh(std::string filename);

	void TestDrawSurfaceMesh();
	void drawSurfaceMesh();
	void drawSurfaceMesh_skel();

	void triangle_meshing();

	void updateDrawData_SurfaceMesh();
	void TestupdateDrawData_SurfaceMesh();
	void updateDrawData_SurfaceMesh_CutPlane();

	void pick_point_with_depth(float x, float y, double(&p)[3], double itemPosition[3]);
	
	/* 六面体函数 */
	int loadMesh(std::string filename);
	bool saveMesh(std::string filename);

	void regularization_hexahedron(bool, int);
	void smoothing_hexahedron(int iter);
	void subdivision(uint16_t);

	void drawHexMesh();// 根据 vertices_data 绘制 HexMesh
	void updateDrawData_HexMesh();
	void updateDrawData_HexMesh_CutPlane();
	int getNumOfModules();
	int delete_islands();
	int delete_redundancy_elements();

	inline bool isVaildHexMesh() { return this->hexMesh; }

	void read_constraint_force(string str);
	void read_constraint_spc(string str);
	void write_constraint_force(string str);
	void write_constraint_spc(string str);
	void remove_all_constraints();
	void smoothing(int iter, bool flat_surface_preserve = false);
	void remeshing();// 重新网格化
	/* 四面体生成 */
	std::string convertMeshToOffString();
	bool loadTetMeshFromStr(const std::string& data);

	void save_as_iges(string filename);// 将六面体转化为样条体

	void get_mesh_quality(string& str);

	Vex get_point_via_vh(const VH&);

	void updateDrawData_TetMesh();
	void updateDrawData_TetMesh_CutPlane();
	void updateSurfaceMesh_skel();
	void drawTetMesh();

	void updateDrawData_QuadMesh();
	void drawQuadMesh();

	void drawQuadMeshfromHexmesh();
	//bool QuadMeshToSpline();
	//六面体拓扑操作
	//bool hexMeshConvertToIGES();

public:
	
	/* 绘制点选对象 */
	void drawSelectedObject();
	void drawSelectedObject_triMesh();
	void drawSelectedObject_hexMesh();
	void drawSelectedObject_quadMesh();
	void drawSelectedObject_tetMesh();
	void drawSelectedObject_SKEL();

	// wuyuan
	void draw_skeleton_selected_point();
	void draw_possible_add_point();
	void DrawDiscretizeSphere(double x, double y, double z, double radius);
	/*void updateDrawDate_Skel();*/
	void DrawPointRadius(MeshKernel::iGameVertexHandle, double radius);
	void DrawPointRadius();

	QMatrix4x4 translate_matrix_ = QMatrix4x4(), rotate_matrix_ = QMatrix4x4();	// 想实现固定管线的glTranslate() glRotatef(),
	QMatrix4x4 model_matrix_ = QMatrix4x4();

	void SetModelMatrix(ModelDraw& modelDraw, QMatrix4x4 matrix);
private:
	void Translated(double x, double y, double z);
	double len = 1.0f;

public:

	void viewCenter();//视点中心
	void EnableLighting(bool b);//光照
	void LoadRotation(void);
	void ResetView(void);//重置视图
	void getPatchStatus(bool status);//设置图元改变状态

	//不同的编辑模式，供鼠标事件选择
	enum EditModel
	{
		DEFULTEDIT, PICKITEM, COORDINATE, PICK_CELLS_SKEL, SET_CELL_SKEL, CON_CELL_SKEL, DEFORMATION, TEST,
		POINT_COORDINATE,
		Skel_Move_Point,
		Skel_Point_add,
	};
	void SetEditModel(const EditModel& editmodel);

	//写错了fivewuyuan
	void PoissonReconstructionUsingPointCloud();

	void loadShader_Possible_Sphere();
	void loadShader_Skel_Sphere(MeshKernel::iGameVertexHandle Node_h); //用于改变VAO VBO EBO
	void loadShader_Skel_Cylinder(MeshKernel::iGameEdgeHandle Edge_h);
	void SendToVAO_Skel();
	void SendChangeToVAO_Skel();

	void InitAllPara();
	void InitOriginMesh();
	void CADContrast();

	//显示点云
	void TurnToCloudPoint();
	//骨架生成的三角网格传入CGAL库中做Union操作
	void GeneratetTriangularMeshforUnionCGAL();

	void write_off_mesh(const MeshCGAL& mesh, const std::string& filename);

protected:
	void DrawScene();//绘制模型
	void DrawSceneModel(void);//传入模型
	void DrawStyle(void);// 绘制样式：点、线、面片、平滑面
	void ColorStyle(void);// 着色样式：模拟光照着色、根据点法向着色
	void updateCurrentSkelSphere(MeshKernel::iGameVertexHandle Node_h); // 更新vertices 数据
	void updatePossibleSphere();

	// five_yuan //用于将控制网格转换为四边形网格
	void ControlMeshChangeToQuadMesh();
	void drawMeshCloudPoint();
	bool IsReadyCloudPoint = false;

	void DrawSkelDiscretize(MeshKernel::iGameVertexHandle Node_h); //绑定VAO 画画

	//void DrawSkelDiscretize_Cylinder(MeshKernel::iGameEdgeHandle Edge_h);
	void InitIndices();
	void InitVaoSphereSize();

	void InitVaoCylinderSize();
	void updateCurrentCylinderSkel(MeshKernel::iGameEdgeHandle Edge_h);
	void DrawSkelDiscretize_Cylinder(MeshKernel::iGameEdgeHandle Edge_h);
	void DrawPossibleSphere();


	//std::vector<GLfloat> vertices;

	std::vector<glm::vec3> positions;
	std::vector<glm::vec2> uv;
	std::vector<glm::vec3> normals;
	std::vector<unsigned int> indices; // 这个能不能公用
	std::vector<float> data;
	std::vector<vector<float>> data_Sphere;
	std::vector<vector<float>> data_Cylinder;
	vector<MeshKernel::SurfaceMesh> SphereMeshForSkeletonPoint;
	vector<MeshKernel::SurfaceMesh> CylinderMeshForSkeletonPoint;

	//std::vector<MeshKernel::SurfaceMesh> Sphere_Mesh;
	Point_set CloudPointForPoissonReconstruction;
	int TestNomal = 1;
	int totalSmapleNum = 0;
	//int NumOfSample = 3;

	

	std::vector<iGameVertex> Vertices_Origin; //用于对比拟合误差的三角形网格构造
	std::vector<std::vector<iGameVertexHandle>> Faces_Origin; //用于对比拟合误差的三角形网格构造
	int PreNumHandle = 0;
	MeshKernel::SurfaceMesh*  SphereMesh_Origin = nullptr;
	std::vector<iGameVertex> VerticesMeshCADForContrast; // 用于记录后面CAD模型的四边形网格顶点
	std::vector<std::vector<iGameVertexHandle>> FacesMeshCADForContrast; // 用于记录后面CAD模型的四边形网格面
	MeshKernel::SurfaceMesh* SphereMesh_MeshCADForContrast = nullptr;

	std::vector<glm::vec3> positions_Cylinder_Zero;
	std::vector<glm::vec3> positions_Cylinder_One;
	std::vector<glm::vec3> normals_Cylinder;
	std::vector<unsigned int> indices_Cylinder;
	//std::vector<float> data_Cylinder;

	const unsigned int X_SEGMENTS = 16; //切的刀数
	const unsigned int Y_SEGMENTS = 16;
	const unsigned int RadiusScaling = 5;

	// 平均半径
	double avg_radius = 0.0;
	double max_radius = -1;
	double min_radius = 9999;
	bool IsCalAvgRadius = false;

	unsigned int indexCount; //这个是要去掉的

	const unsigned int Circle_SEGMENTS = 10;

	void pickEvent(double winX, double winY);
	void pickEvent_TriMesh_BVH(double winX, double winY);// 新方式，快
	void pickEvent_HexMesh_BVH(double winX, double winY);
	void pickEvent_QuadMesh_BVH(double winX, double winY);
	void pickEvent_TetMesh_BVH(double winX, double winY);
	void pickEvent_TriMesh(double winX, double winY);// 旧方式，慢
	void pickEvent_HexMesh(double winX, double winY);

	void pickPoint(double winX, double winY);// 获取二维坐标的三维点 (伟伟)
	void pickPoint(int x, int y, vector<double>& p);// 获取二维坐标的三维点
	void pickPoint_getDepth(int x, int y, vector<double>& p, double& d);// 获取点击位置二维坐标的三维点和深度，用于鼠标选点
	void pickPointDepth(int x, int y, vector<double>& p, double& d);// 使用深度获取三维坐标对应的三维点，用于拖动和绘制点

	void mousePressEvent(QMouseEvent*);
	void mouseMoveEvent(QMouseEvent*);
	void mouseReleaseEvent(QMouseEvent*);
	void wheelEvent(QWheelEvent*);
	void keyPressEvent(QKeyEvent*);
	void keyReleaseEvent(QKeyEvent*);

	// 线与点的距离
	double point2LineDis(QVector3D point1, QVector3D point2, QVector3D point0);
	// 判断线是否穿过三角形
	bool isIntersectTriangle(QVector3D orig, QVector3D end, QVector3D v0, QVector3D v1, QVector3D v2, QVector3D& intersection);

	

protected:

	EditModel editmodel;//编辑模式
	int upatchnum = 10;//渲染时的片数
	int vpatchnum = 10;//渲染时的片数
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool PatchChange = false;//判断是否更新图元绘制数量
	bool isControlLine = false;//是否显示控制曲线
	bool isControlMesh = false;//是否显示控制网格

	bool ControlMeshHasFinished = false; //控制网格是否已经被绘制过了

	double depth = 0.996766;//获取点击位置的深度
	int curvepointIndex[2];//保存更改曲线的控制点索引
	int surfacepointIndex[3];//保存更改曲面的控制点的索引

	// 保存鼠标选点信息
	Vector3d clipoint;

	enum PickType 
	{// 选择类型
		PickPoint, PickEdge, PickFace, PickCell, PickCells
	};

	PickType pickType = PickPoint;

public:

	// Bounding Box
	Vector3d bbMax, bbMin;

	bool drawSelectedAllFlag = true;
	bool drawMeshLine = true;
	double meshline_width = 1.25f;

	// 切片
	bool cutFlag = false;
	//double cutX = 0.f, cutY = 0.f, cutZ = 0.f;// 0 全部显示，1 全部不显示

	bool resetViewFlag = true;

	void setPickType(int mode) 
	{
		switch (mode) {
		case 1:
			pickType = PickEdge;
			break;
		case 2:
			pickType = PickFace;
			break;
		case 3:
			pickType = PickCell;
			break;
		case 4:
			pickType = PickCells;
			break;
		default:// case 0
			pickType = PickPoint;
			break;
		}
	}

	void removeSelectedOne();
	void removeSelectedAll();
	void deleteSelectedAll();
	void add_new_cell();

	bool updateAllHandles();
	
	/*Temp-lee*/
public:
	// 用以区分不同的连通分量
	std::vector<int> module_colors;
	// 随机生成的颜色，用以区分不同面片
	std::vector<Eigen::Vector3d> colors;

	//曲线曲面绘制集合
	vector<CBSplineCurve> curvesInfo;
	vector<CBSplineSurface> surfacesInfo;
	
	//所有曲线曲面转成B样条进行绘制
	vector<CBSplineCurve> curvesDraw;//样条曲线集合
	vector<CBSplineSurface> surfacesDraw;//样条曲面集合
	
	//曲线曲面集合
	vector<CBSplineCurve> bsplineCurves;//样条曲线集合,保存当前操作的曲线集
	vector<CBSplineSurface> bsplineSurfaces;//样条曲面集合，保存当前操作的曲面集

	/* 保存曲线曲面的、控制网格、边界的离散采样数据 */
	// 所有类型的曲线曲面数据
	// 曲线曲面控制点
	vector<vector<CPoint>> controlLine;// 曲线控制点线
	vector<vector<vector<CPoint>>> controlMesh;// 控制点网格

	vector<std::vector<CPoint>> data_curves;// 曲线离散点集合，data_curves[i][j]表示第i条曲线的第j个采样点
	vector<vector<std::vector<CPoint>>> data_surfaces;// 曲面离散点集合,data_surface[i][j][k]表示第i个曲面的第j个面片的第k个点（0<=k<=3）

	vector<vector<vector<CPoint>>> data_bound_surface;// 曲面边界曲线离散点集合,datda_bound[i][j][k]表示第i个曲面的第j个边界（0<=j<=3)的第k个点
	vector<vector<CPoint>> data_controlLine;// 控制曲线离散点
	vector<vector<vector<CPoint>>> data_controlMesh; // 控制点网格离散点 是来自于controlMesh

	vector<vector<vector<CVector>>> surface_normals;

	CPoint curvecenter, surfacecenter;//模型中心点

	


	bool bshowAxis;
	bool bshowGround;
	bool bshowSKEL;
	bool bshowTriMesh;
	bool bdragMode;
	bool bVertexSnap;
	SolutionSkel* solutionSkel;
	Project* projectSkel;
	HexV3fMesh* oldHexMesh;
	BMesh::Node* root;

	std::vector<BMesh::Node*> nodeVec;
	std::unordered_map<int, bool> visited;
	std::queue<BMesh::Node*> Q;
	double px; double py; double pz;	// 点击ray交点
	double oldDepth = 0;
	void setNodeVec(std::vector<BMesh::Node*>& v);
	void setNodeVec(const std::vector<BMesh::Node*>& v);
	void setRoot(BMesh::Node* _root);
	void BFSTraverse(BMesh::Node* root);
	void BFS(BMesh::Node* root);
	
	void drawAxis();
	void drawHorizon();
	void drawSkel();
	void drawTriMesh_skel();
	void drawProjectLine();
	void translationCell(QPoint p);
	bool pickEvent_SKEL(double winX, double winY);
	void Rotation(const QPoint& p);
	int loadSkeletonLee(std::string filename);
	// wr.z
	bool copyTime;

	bool copyTime_point = false;
	MeshKernel::iGameVertexHandle vh_copy_add; // 这个是骨架新生成的点
	MeshKernel::iGameEdgeHandle e_add_h1;

	void Sleep(int msec);
	void initializePick();
	int findNearestCell(double aim[3]);


	
	void rupdateGL();
	void drawNormalSkel();
	void drawNormalSkel_1();
	void rselectALine();
	void rconnectNodes();
	void rdisconnectNodes();

	void project();
	void PIA();
	void subdivisionNormal();
	void edgeSubdivision(int divnum);

	bool hexMeshConvertToIGES();
	bool loadIGESSurface(std::string fileName);
	bool saveIGESFile(std::string fileName);

	//five_yuan
	bool saveBSplineToStep(std::string fileName);
	bool QuadMeshToSpline();
		

	bool loadIGESSurface(std::vector<CBSplineSurface>& geomsurface);
	bool loadIGESSurface(std::vector<CBSplineSurface>& geomsurface, SurfaceMesh& trimesh);
	bool loadSurfaceErrorIGES(std::vector<CBSplineSurface>& geomsurface, SurfaceMesh& trimesh);
	void bsplineSurfaceToDiscrete();
	void controlpointsToDiscrete();
	void updateSurfaceModel();
	void drawControlMesh();
	void drawSurface();
	void DrawPointForOrder();
	void drawControlPoint(double xx, double yy, double zz, double radius, GLfloat M, GLfloat N);
	void updateDrawData_Surface();
	void updateDrawData_Surface_IGES_Error(SurfaceMesh& trimesh);
	void updateDrawData_SurfaceMesh_Error_IGES();
	std::vector<vector<std::vector<CPoint>>> getData_boundcurve(vector<CBSplineSurface> boundsurface);
	HexV3fMesh* convent();
	void genColors(size_t sz = 10007);

	//five_yuan
	void updateControlPointAndLine();//用于传入控制点的坐标到opengl
	void drawControlPointAndLine(); //画beizer的控制点和线

/* 扩展类 CoordinateAxis */
public:
	CoordinateAxis::ModeType select_coordinater_mode;
	CoordinateAxis* local_coordinater_axis;

	friend void CoordinateAxis::Draw3DCoord(ModelDraw&);
	friend void CoordinateAxis::SetModelMatrix(ModelDraw& modelDraw, QMatrix4x4 matrix);
	friend bool CoordinateAxis::PickEvent(ModelDraw& modelDraw, double winX, double winY);
	friend QVector3D CoordinateAxis::MoveEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);
	friend QVector4D CoordinateAxis::RotateEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);
	friend QVector3D CoordinateAxis::FlexibleEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);

	/* 扩展类 ProDraw, HexDraw */
	friend class ProDraw;
	ProDraw* proDraw;


	//class DotLineCalculation;
	friend class DotLineCalculation;
	DotLineCalculation* dotLineCalculation;


	bool use_pro_draw = false;

	/* MMC 组件转换(LING Ran) */
	bool BuildSkeletonToFile(std::string filename, uint32_t width, uint32_t height);
	void SparseOptimization(double max_error = 1.f);
	void GenerateSymAxis();
	void ShapeApproximation(uint32_t max_iter);

	void drawPolyMesh_OPT();
	void updateDrawData_PolyMesh_OPT_1();
	void updateDrawData_PolyMesh_OPT_2();
	void updateDrawData_PolyMesh_OPT_V();


/* HexBuildDraw.cpp */
	HexMeshTool* hex_mesh_tool_;
	void selectedUpToDown();

	void refinement_test();


};