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
#include <QtGUI/Qsurfaceformat.h>
#include <QtGUI/QPixmap>
#include <QPainter>
/* Eigen */
#include <Eigen/Dense>
#include <Eigen/Core>
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
/* Model Reconstruction */
#include "BezierSolid.h"
#include "FileHandler.h"
#include "Subdivision.h"
#include "Transform.h"

#include "GenHexSkel.h"
#include "SurfaceFitter.h"
#include "PIA.h"
#include <Geom/CAidFunction.h>

/* iGame Kernel */
#include "MeshKernel/Mesh.h"
#include "MeshKernel/IO.h"
#include "MeshKernel/MeshFileIO.h"
#include "SurfaceAlgorithm/OpenMeshAlgorithm/OpenMesh_Typedef.h"
#include "VolumeAlgorithm/Volume_Algorithm.h"
#include "SurfaceAlgorithm/Surface_Algorithm.h"
/* Sparse Optimization */
#include "SkeletalPipeline/Pipeline.h"
#include "solutionSkel.h"
#include "BMesh/project.h"
#include "BMesh/mlca_subdivision.h"
/*    Tools    */
#include "Tools/BVH.h"
#include "CoordinateAxis.h"
#include "edge_subdivision.h"

/*	Skeletal Mesh	*/
#include "SkeletalMesh/SkeletalMesh.h"
#include "SkeletalMesh/Q-MAT.h"
#include "SkeletalMesh/Subdivision_Skl.h"

#include "VolumeAlgorithm/Constraint.h"
#include "Tools/meshmath.h"
#include "SkeletalMesh/QuadMesher.h"

// 样条体渲染
#include "CBSplineSolid.h"
#include "Bdy.h"

#include <atomic>

using namespace std;

class ModelDraw :public ViewOpenGL {
	Q_OBJECT
public:
	ModelDraw(QWidget* parent = 0);
	virtual ~ModelDraw(void);

	/* 输出信息 */
	void setPickInfo(std::function<void(void)> _func) {
		showPickInfo = _func;
	}
	std::string pick_info;

	// 拓扑优化
	void updatePickPara();
	bool analysis_flag = false;
	std::atomic_bool is_need_update_topopt_model = false;
    std::atomic_bool is_need_update_spline_simulation_model = false;
    std::mutex update_spline_simulation_model_lock;

	int analysis_mode = 1;
	std::vector<double> analysis_data[3];
	std::vector<std::vector<int>> pickface;
	std::vector<std::vector<int>> pickedge;
	std::vector<std::vector<int>> pickcell;
	IGA::HexBdy hexbdy;
	// IGA::CCBezierBdy ccbdy;
    std::vector<int> topopt_pro;
	void Read_Bdy(std::string& filename);
	/// 
	void savePickItem(string name);
	bool loadPickFile(const std::string& filename);
	void addPickData(vector<double>, int);
	void TobeProtected();

	
	bool show_tangent_plane = false;
	bool compare_flag = false;
	std::vector<GLfloat> compared_vertices_data;// Just draw points
	std::unordered_map<int, int> compared_to_current;
	void compared_bind_vertices();
	int used_last_vcnt = 0;
	std::vector<GLfloat> translucence_vertices_data;

public:

	std::function<void(void)> showPickInfo;

	GLuint VAO = 0, VBO = 0, EBO = 0;// 模型数据

	// 如果拓扑被改变，将会有无效的handle，需要调用 updateAllHandles() 更新所有handle
	// surface_mesh->updateAllHandles()
	bool topology_change_flag = false;


	// 绘制数据
	std::vector<GLfloat> vertices_data;
	std::vector<GLfloat> fixedVertices;
	std::vector<GLfloat> forceVertices;
	std::vector<GLfloat> forceEdges;

	size_t tri_cnt = 0, quad_cnt = 0, line_cnt = 0, translucence_tricnt = 0;// 记录绘制数据信息
	size_t v_cnt = 0, e_cnt = 0;
	size_t pick_cnt = 0, force_cnt = 0;

	std::vector<Streamline> streamlines;
	std::vector<SingularPoint> singularity_vertices;

	/* 骨架网格函数 */
	SkeletalMesh* skeletal_mesh = nullptr;
	bool scale_to_uint_flag = true;
	

	// 稀疏优化
	std::vector<std::pair<Eigen::Vector2d, Eigen::Vector2d>> axis;
	std::vector<std::vector<Eigen::Vector2d>> bdy;
	polyMesh* poly_mesh = nullptr;
	Pipeline* pipeline = nullptr;

	// OpenMesh 三角形网格
	iGame_OpenMesh_TriMesh* tri_mesh = nullptr;

	// 面网格
	MeshKernel::SurfaceMesh* surface_mesh = nullptr;
	std::vector<int> fixed_indices;
	std::vector<pair<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex>> moved_vertices;
	Deformation* tri_deformation = nullptr;

	std::map<int, MeshKernel::SurfaceMesh> remeshing_result;
	std::map<int, MeshKernel::SurfaceMesh> remeshing_hole_result;
	std::map<int, MeshKernel::SurfaceMesh> remeshing_holefilled_result;
	std::map<int, MeshKernel::SurfaceMesh> deletehole;
	std::map<int, MeshKernel::SurfaceMesh> fillpatch;

	int remeshing_show_model = 0;

	// 体网格
	MeshKernel::VolumeMesh* hexMesh = nullptr;// 六面体网格
	MeshKernel::TetMesh* tetMesh = nullptr;// 四面体网格

	TetDeformation* tet_deformation = nullptr;

	// 拓扑优化约束数据
	std::unordered_map<int, bool> force_id;// 约束顶点ID
	std::unordered_map<int, bool> spc_id;

    // material
    double YoungsModulus = 1e5;
    double PoissonsRatio = 0.3;

    //hyperelasticity time-dependent simulation
    double Duration_he = 2;
    int FrameCount_he = 1;
    double TimeStep_he = 0.1;

	// 点选数据
	BVH_Tree* bvh_tree = new BVH_Tree();
	bool mouse_move_flag = false;
	std::vector<int> selected_v_indices;
	std::vector<int> selected_e_indices;
	std::vector<int> selected_f_indices;
	std::vector<int> selected_c_indices;

    std::vector<int> force_v_idx;
    std::vector<int> force_e_idx;
    std::vector<int> force_f_idx;
    std::vector<int> force_c_idx;
    std::vector<int> fix_v_idx;
    std::vector<int> fix_e_idx;
    std::vector<int> fix_f_idx;
    std::vector<int> fix_c_idx;



	// 选择框数据
	vector<QVector3D> selected_rectangle;

public:

	void test_action(double);

	void loadShader();
	void genColors(size_t sz = 10007);
	void initColorbar();
	Vector3d getColor(double c);// c should in [0, 1]
	std::vector<pair<double, Vector3d>> colorbar;

	void ModelClear();// 清除所有模型数据
	void updateDrawData(bool need_update_bvhtree = true);// 更新绘制数据
	void updateDrawData_CutPlane();
    void updateDrawData_Surface_IGES_Color(int dirID = 0);
	void updateMeshModel();// 根据 bbox 设置中心
	void updateBVH();

	double getVolume();// 返回模型的体积
	int getVerticesSize();

	/* 骨架网格函数 */
	void loadSkeletalMesh(std::string);
	void saveSkeletalMesh(std::string);
	void updateDrawData_SkelMesh();
	void drawSkelMesh();

	void simplification_qmat(int target_vsize);
	void subdivision_skl(int iter);
	void calc_cross_field();

	/* 稀疏优化函数 */
	void drawPolyMesh_OPT();
	void updateDrawData_PolyMesh_OPT_1();
	void updateDrawData_PolyMesh_OPT_2();
	void updateDrawData_PolyMesh_OPT_V();

	void updateDrawData_DensityChange();
    void updateDrawData_Topopt();
	/* OpenMesh 三角形网格函数 */
	bool loadTriMesh(std::string filename);
	bool saveTriMesh(std::string filename);
	void drawTriMesh();
	void updateDrawData_TriMesh();

	void smoothing_openmesh(size_t);
	void remeshing_uniform_openmesh(double);
	void simplification_qem_openmesh(size_t);
	void denoising_bnf_openmesh(size_t);

	/* 曲面网格函数 */
	bool loadSurfaceMesh(std::string filename);
	bool saveSurfaceMesh(std::string filename);

	void drawSurfaceMesh();

	void updateDrawData_SurfaceMesh();
	void updateDrawData_SurfaceMesh_Analysis();
	void updateDrawData_SurfaceMesh_CutPlane();

	void remeshing_adaptive();
	void remeshing_uniform();
	void denoising_bnf(int iter = 1);
	void subdivision_loop(int iter = 1);
	void simplification_qem(int target_vertices_size, bool powered_flag = false);
	void parameterization(int mode);
	void deformation_arap();

	void update_fixed_indices();// 将选中的所有点固定

	std::string add_moved_vertex(double x, double y, double z);// 将用户指定点移动信息保存

	void pick_point_with_depth(float x, float y, double(&p)[3], double itemPosition[3]);
	/* 四面体函数 */
	void loadTetMesh(std::string filename);
	void saveTetMesh(std::string filename);

	void drawTetMesh();
	void updateDrawData_TetMesh();
	void updateDrawData_TetMesh_CutPlane();

	/* 四面体生成函数 */
	std::string convertMeshToString();
	bool loadTetMeshFromStr(const std::string&);

	/* 六面体函数 */
	int loadHexMesh(std::string filename, bool reset_view = true);
	bool saveHexMesh(std::string filename);

	void regularization_hexahedron(bool, int);
	void smoothing_hexahedron(int iter);
	void subdivision(uint16_t);
	void add_noise(double ratio);

	void drawHexMesh();// 根据 vertices_data 绘制 HexMesh
	void updateDrawData_HexMesh();
	void updateDrawData_HexMesh_Analysis();
	void updateDrawData_HexMesh_CutPlane();
	int iGameDFS();
	int delete_island_HexMesh();
	int delete_redundancy_cell_HexMesh();
    void generate_coverage_axis();
    void hex_meshing_test();
    void test_det_faces();

	void read_constraint_force(string str);
	void read_constraint_spc(string str);
	void write_constraint_force(string str);
	void write_constraint_spc(string str);
	void remove_all_constraints();

	inline bool isVaildHexMesh() { return this->hexMesh; }


	bool hexMeshConvertToIGES();
    bool hexMeshConvertToIGES_2(std::string model_name, int Smooth_Times);
	bool hexMeshConvertToIGES_skeletonExtraction();
	bool mesh_to_iges(std::string meshFileName, std::string igesFileName);

public:
	/* 样条曲面函数 */
	bool loadIGESSurface(std::string filename);
	bool loadIGESSurface(std::vector<CBSplineSurface>& geomsurface);
    bool loadIGESSurfaceColor(std::vector<CBSplineSurface>& geomsurface, std::vector<CBSplineSurface>& geomsurface_color);

	bool saveIGESFile(std::string filename);
	bool saveIGESFromMesh(std::string filename);

	void bsplineSurfaceToDiscrete();
    void bsplineSurfaceWireLineToDiscrete();
	void updateSurfaceModel();
	void drawSurface();
	void updateDrawData_Surface();
	void updateDrawData_Surface_CutPlane();



	/* 绘制点选对象 */
	void drawSelectedObject();
	void drawSelectedObject_triMesh();
	void drawSelectedObject_tetMesh();
	void drawSelectedObject_hexMesh();
	void drawSelectedObject_sklMesh();
	void drawSelectedObject_SKEL();// Lee

	vector<vector<std::vector<CPoint>>> getData_boundcurve(vector<CBSplineSurface>);//获取样条曲面边界处采样点

	// 样条曲线（未启用）
	bool loadIGESCurve(std::string filename);
	void bsplineCurveToDiscrete();
	void updateCurveModel();
	void drawCurve();
	// 未启用
	void controlpointsToDiscrete();
	void drawControlLine();
	void drawControlMesh();
	void drawControlPoint(double xx, double yy, double zz, double radius, GLfloat M, GLfloat N);

public:

	void viewCenter();//视点中心
	void EnableLighting(bool b);//光照
	void LoadRotation(void);
	void ResetView(void);//重置视图
	void getPatchStatus(bool status);//设置图元改变状态
	//不同的编辑模式，供鼠标事件选择
	enum EditModel {
		DEFULTEDIT, PICKITEM, COORDINATE, PICK_CELLS_SKEL, SET_CELL_SKEL, CON_CELL_SKEL, DEFORMATION, BSPLINESOLIDEDIT
	};
	void SetEditModel(const EditModel& editmodel);

protected:
	void DrawScene();//绘制模型
	void DrawSceneModel(void);//传入模型
	void DrawStyle(void);// 绘制样式：点、线、面片、平滑面
	void ColorStyle(void);// 着色样式：模拟光照着色、根据点法向着色


	void pickEvent(double winX, double winY);
	void pickEvent_TriMesh(double winX, double winY);
	void pickEvent_TriMesh_BVH(double winX, double winY);
	void pickEvent_TetMesh(double winX, double winY);
	void pickEvent_TetMesh_BVH(double winX, double winY);
	void pickEvent_HexMesh(double winX, double winY);
	void pickEvent_HexMesh_BVH(double winX, double winY);
	void pickEvent_Skeleton(double winX, double winY);
    void pickEvent_Bsplinesurface(double winX, double winY);
	

	void pickPoint_getDepth(int x, int y, vector<double>& p, double& d);// 获取点击位置二维坐标的三维点和深度，用于鼠标选点
	void pickPointDepth(int x, int y, vector<double>& p, double& d);// 使用深度获取三维坐标对应的三维点，用于拖动和绘制点

	//ToDo:重写鼠标事件，增加点击操作
	void clickCurveControlPoint(QPoint point);
	void clickSurfaceControlPoint(QPoint point);
	void moveCurveControlPoint(QPoint point);
	void moveSurfaceControlPoint(QPoint point);
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
	int wpatchnum = 10;
	bool isEnableLighting;
	bool isTwoSideLighting;
	bool PatchChange = false;//判断是否更新图元绘制数量
	bool isControlLine = false;//是否显示控制曲线
	bool isControlMesh = false;//是否显示控制网格
	bool isHexMesh = false;//是否显示B样条体
	bool isRealtime = false;//是否实时渲染
	bool isShowinside = false;//是否渲染内部
	bool ifchooseface = false;//是否选面
	double depth = 0.996766;//获取点击位置的深度
	int curvepointIndex[2];//保存更改曲线的控制点索引
	int surfacepointIndex[3];//保存更改曲面的控制点的索引

    bool bshowWireLine = false;
    bool have_light = false;

    vector<vector<vector<Vector3d>>> data_surfaces_color[3];

	// 保存鼠标选点信息
	Vector3d clipoint;

	enum PickType {// 选择类型
		PickPoint, PickEdge, PickFace, PickCell, MoveSKEL, PickCells, SetCell, RotateCell
	};

	PickType pickType = PickPoint;

public:

	// Bounding Box
	Vector3d bbMax, bbMin;

	bool drawSelectedAllFlag = true;
	bool drawMeshLine = true;
	double meshline_width = 1.5f;

	// 切片
	bool cutFlag = false;
	double cutX = 0.f, cutY = 0.f, cutZ = 0.f;// 0 全部显示，1 全部不显示
	double topopt_density = 0.5;
    double angle_range_value = 0.;
    std::vector<iGameFaceHandle> RangeTodoFaces;
	std::vector<double> hex_density;
	bool resetViewFlag = true;

	void enable_realtime_deformation();

	inline void disable_realtime_deformation() {
		SetEditModel(EditModel::DEFULTEDIT);
	}



	void setPickType(int mode) {
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
			pickType = MoveSKEL;
			break;
		case 5:
			pickType = PickCells;
			break;
		case 6:
			pickType = SetCell;
			break;
		case 7:
			pickType = RotateCell;
			break;
		default:// case 0
			pickType = PickPoint;
			break;
		}
	}

	/* 稀疏优化(LING Ran) */
	bool BuildSkeletonToFile(std::string filename, uint32_t width, uint32_t height);
	void SparseOptimization(double max_error = 1.f);
	void GenerateSymAxis();
	void ShapeApproximation(uint32_t max_iter);

	void removeSelectedOne();
	void removeSelectedAll();
	void deleteSelectedAll();
	void add_new_cell();

	bool updateAllHandles();

	/* 测试函数 */
	void test_bspline();



	//保存要绘制的样条体

    Poster_Needed::CBSplineSolid bsplinesolid;
	double ControlPointSize=0.01;
	double ControlLineWidth = 50* ControlPointSize;
	std::vector<Poster_Needed::CBSplineSolid> bsplinesolidmesh;
	int getBsplineSolidFromParamXJ(int& udegree, int& vdegree, int& wdegree, vector<vector<vector<CPoint>>>&controlPoints, vector<double>&uknots, vector<double>&vknots, vector<double>&wknots, vector<int>&umultis, vector<int>&vmultis, vector<int>&wmultis);//获取构造体
	void getBsplineSolidFromParam(int& udegree, int& vdegree, int& wdegree, vector<vector<vector<CPoint>>>& controlPoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis);//获取构造体
	void ModelDraw::trans11(std::vector<double>&modelMatrix_, QMatrix4x4 & modelmatrix);
	//void bsplinesolid_controlpointsToDiscrete();
	void bsplinesolid_drawControlMesh(void);
	void ModelDraw::moveSolidControlPoint(double winX, double winY);
	//void ModelDraw::pickPointDepth(int x, int y, double p[3], double& d);// 使用深度获取三维坐标对应的三维点，用于拖动和绘制点
	//void getPointDepth(int x, int y, double point[3], double& dz);//获取点击位置二维坐标的三维点和深度，用于鼠标选点
	////void bsplineSolidToDiscrete();
	void bsplineSolidSimulation();
	void bsplineSolidToDiscreteXJ();
	void getControlMeshStatus(bool status);
	void draw_BSplineSolid_HexMesh();// 根据 vertices_data 绘制 HexMesh
	void showHexMesh(bool status);
	void set_isRealtime(bool status);//是否实时渲染
	void set_isShowinside(bool status);//是否渲染内部
	void set_chooseface(bool status);//是否进入选面模式

	void aotu_createBSS(int ucnum, int vcnum, int wcnum, int udegree, int vdegree, int wdegree, vector<vector<vector<CPoint>>>& controlpoints, vector<vector<CPoint>>& singlepatch, vector<CPoint>& singlepoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis);
	//void manually_createBSS(int ucnum, int vcnum, int wcnum, int udegree, int vdegree, int wdegree, vector<vector<vector<CPoint>>>& controlpoints, vector<vector<CPoint>>& singlepatch, vector<CPoint>& singlepoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis);
	void pickPoint(double x, double y);
	double pickedItemPosition[3];
	int controlLineOrMesh = 0;// 1代表controlLine，2代表controlMesh
	struct CPFunc //可参考源码模仿
	{
		//重载函数调用操作符，实现CPoint类的排序功能
		constexpr bool operator()(const CPoint& p1, const CPoint& p2) const
		{
			if (p1.m_dX == p2.m_dX && p1.m_dY == p2.m_dY && p1.m_dZ == p2.m_dZ)return 0;
			return p1.m_dX < p2.m_dX || (p1.m_dX == p2.m_dX && p1.m_dY < p2.m_dY) || (p1.m_dX == p2.m_dX && p1.m_dY == p2.m_dY && p1.m_dZ < p2.m_dZ);
		}
	};
	void SaveIndex(vector<vector<vector<CPoint>>> controlmesh,int handle);
	map<CPoint, int, CPFunc> Is_Highlight;//存放高亮的控制顶点
	set<CPoint,CPFunc> AllControlPoint;//存放所有控制顶点
	map<CPoint, vector<int>,CPFunc> PointBelongToSolid;//存放每个控制点对应的样条体
	void ClearAllVector();//清除所有容器和bsplinesolidmesh的数据
	vector<int> singlepatch_vertics_number;

	void pickPoint_chooseface(double winX, double winY);
	int chooseface_index = -1;
	int bp_index=-1;
	int face_handle=-1;
	typedef struct {
		int block_id;
		int face_id;
		Eigen::Vector3f force;
	} Force;
	Force My_Force;
	typedef struct {
		int block_id;
		int face_id;
		int dir;
		Eigen::Vector3f force;
	} Offset;
	Offset My_Offset;
	void set_F(double x,double y,double z);
	void set_V(double x, double y, double z,int dir);
public:
	/*Temp-lee*/
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
	vector<vector<vector<CPoint>>> data_controlMesh;// 控制点网格离散点

	vector<vector<vector<CVector>>> surface_normals;

	CPoint curvecenter, surfacecenter;//模型中心点
    vector<vector<CPoint>> data_surface_line;

	bool ctlPressed;// scale
	bool bshowAxis;
	bool bshowGround;
	bool bshowSKEL;
	bool bdragMode;
	bool bshowTriMesh;
	bool bRotateByMouse;
	bool bVertexSnap;
	SolutionSkel* solutionSkel;
	Project* projectSkel;
	HexV3fMesh* oldHexMesh;
	BMesh::Node* root;

	std::vector<BMesh::Node*> nodeVec;
	std::unordered_map<int, bool> visited;
	std::queue<BMesh::Node*> Q;
	double px; double py; double pz;
	double oldDepth = 0;
	void setNodeVec(std::vector<BMesh::Node*>& v);
	void setNodeVec(const std::vector<BMesh::Node*>& v);
	void setRoot(BMesh::Node* _root);
	void BFSTraverse(BMesh::Node* root);
	void BFS(BMesh::Node* root);
	void drawAxis();
	void drawHorizon();
	void drawSkel();
	void drawSurfaceMesh_skel();
	void translationCell(QPoint p);
	bool pickEvent_SKEL(double winX, double winY);
	void Rotation(const QPoint& p);
	int loadSkeletonLee(std::string filename);
	int loadMesh(std::string filename);
	void subdivisionNormal();
	void edgeSubdivision(int divnum);
	HexV3fMesh* convent();
	// wr.z
	bool copyTime;
	void Sleep(int msec);
	void initializePick();
	int findNearestCell(double aim[3]);

public slots:
	void rupdateGL();
	void drawNormalSkel();
	void drawNormalSkel_1();
	void rselectALine();
	void rconnectNodes();
	void rdisconnectNodes();
	void project();

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
};