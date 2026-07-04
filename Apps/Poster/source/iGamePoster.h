#pragma once
/*
主窗口类，实现有关界面UI和信号&槽的连接
*/
#include <time.h>
#include <QtWidgets/Qtabwidget.h>
#include <QtWidgets/QMainWindow>
#include <QtGUI/QtGui>
#include <QtWidgets/QtWidgets>
#include <QtConcurrent/QtConcurrent>
#include "ui_iGame-Poster.h"
#include "aboutDialog.h"
#include "helpDoc.h"
#include "ModelDraw.h"
#include "singalFather.h"

/* Top Optimization */
#include "TopOpt.h"

#include "Bdy.h"

#include "IsoThermalSolver.h"

// 建模部分
#include "BuildUI.h"
#include "SingalSlot.h"

class ModelView;
class iGamePoster : public QMainWindow
{
	Q_OBJECT
	friend class BuildUI;
	friend class SingalSlot;
private:
	Ui::MainWindow *ui;
	BuildUI* build_ui_;
	SingalSlot* singal_slot_;
	QTimer* topopt_refresh_timer_ = nullptr;

    QButtonGroup buttonGroup_forPick;

	std::string poster_filename;
	
    std::string model_name;
public:
	iGamePoster( QWidget *parent = Q_NULLPTR );
	~iGamePoster(void);
public slots:

	bool savePickItem();//导出点选的信息
	bool openPickItem();//导出点选的信息
	bool addPickForce();//施加力的约束
	bool addPickFixed();//施加固定约束
    bool addHyperelasticityTimeSettings();
	bool TobeProtected();//增加保护区域
    void onCheckBoxClicked(QCheckBox* checkBox1, QCheckBox* checkBox2);

	bool LoadSkeletalMesh(void);
	bool SaveSkeletalMesh(void);

	/**************************** 三角网格函数 ******************************************/
	bool LoadSurfaceMesh();
	bool SaveSurfaceMesh();

	/**************************** 四面体函数函数 ******************************************/
	bool LoadTetrahedronMesh(void);
	bool SaveTetrahedronMesh(void);
	void read_remeshing();
	void set_hole();
	void set_remeshing();
	void set_hole_filled();
	void set_holepatch();
	void set_deletehole();
	/**************************** 六面网格函数 ******************************************/

	bool LoadHexahedralMesh(void);
	bool SaveHexahedralMesh(void);
	
	bool ConvertToIGES(void);
	bool ConvertToIGES_Skel(void);

    /**************************** 后处理 ******************************************/

    bool Line_surface_skel(void);
    bool testDetFaces(void);
    bool testResample(void);
    bool hexmeshingBaseOnMedialMesh(void);
    bool ConvertToIGES_2(void);

	/**************************** 几何函数 ******************************************/

	
	bool OpenIGESSurface(void);
	bool OpenIGESCurve(void);
	bool SaveIGESFile(void);

	/**************************** 信息看板 ******************************************/

	bool InfoClear(void);
	bool InfoCopy(void);

	/**************************** 渲染函数 ******************************************/

	bool ViewStylePoints ();//渲染样式点
	bool ViewStyleLines ();//渲染样式线
	bool ViewStyleFill ();//渲染样式片
	bool ViewNormal();// 默认颜色模式
	bool ViewLight();// 光照模式
	bool ViewClear();// 清除网格

	bool ViewMeshLine();
	bool UpdateMeshLineWidth();

	bool UpdateResetViewFlag();

	/**************************** 截面函数 ******************************************/
	bool ChangeCutFlag();
	bool UpdateCutXYZ();


	/**************************** 点选函数 ******************************************/
	bool ChangePickFlag();
	bool setPickPoint();
	bool setPickEdge();
	bool setPickFace();
	bool setPickCell();

    bool ChangeRangePickFlag();
    bool AddRangePoints();
    bool AddRangeLines();
    bool AddRangeFaces();
    bool AddRangeCells();
    bool GetRangeTempFaces(std::vector<iGameFaceHandle>& faces, double tol);
    bool UpdateAngleRangeValue();
    bool DisablePick();

	bool drawSelectAll();
	bool removeSelectedOne();
	bool removeSelectedAll();
	bool deleteSelectedAll();
	bool addNewCell();

	/**************************** 形变函数 ******************************************/
	bool updateFixedHandle();
	bool changeRealtimeDeformation();

	/**************************** MMC稀疏优化 ******************************************/

	void MMC_GenSkel();
	void MMC_GenAxis();
	void MMC_Opt();
	void MMC_appr();
	void TopOpt_Run();
    void TopOpt_Run_2();
	void TopOpt_Run_4();
	void TopOpt_Run_Seq();
	void Protect_topopt();

	bool UpdateTopoptDensity();
    bool SetFinalTopOptMesh();
	bool isoValueSlider();
	/**************************** IGA Linear ******************************************/
	void Conver_Mesh_To_XML();
	void IGA_Vis_X();
	void IGA_Vis_Y();
	void IGA_Vis_Z();
	void IGA_Linear_Pro();
	void IGA_Thermal();
    void IGA_hyperelasticity_simulation();
    void set_material();
	void Read_Bdy();

	void Poster_Smooth();
	void Poster_DeleteElement();
	void Poster_GenSpline();

	/**************************** 测试函数 ******************************************/
	bool testAction();

	bool simulation();
	bool EditBsplineSolid();//编辑b样条体
	bool auto_ShowBsplineSolid();//绘制B样条体
	bool manually_ShowBsplineSolid();//绘制B样条体
	bool IF_ShowBsplineSolid();//显示b样条体
	bool IF_Real_time_rander();//实时渲染样条体
	bool IF_ShowBsplineSolid_inside();//渲染样条体内部
	bool Choose_Face();//选面模式，修改mousepress鼠标事件
	void generate_F();
	void generate_V();
signals:
	bool editOpen ();//开启编辑模式

};
