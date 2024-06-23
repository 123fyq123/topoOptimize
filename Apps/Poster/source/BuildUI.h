#pragma once
/*
窗口类，实现建模子界面的槽函数
*/

#include <QtWidgets/QtWidgets>
#include "ui_iGame-Poster.h"

//#include "SkeletalModelDraw.h"
class BuildUI : public QWidget
{
	Q_OBJECT
private:
	//SkeletalModelDraw func_;
	Ui::MainWindow* ui;
public:
	BuildUI(Ui::MainWindow* ui, QWidget* parent = Q_NULLPTR);
	~BuildUI(void);
public slots:

	void transformCell();
	void editCell();
	void addCell();
	void connnectCell();

	void editCubeSignal();
	void addCubeSignal();
	
	void selectMutiNodesSignal();
	void selectALineSignal();

	// 六面体模型构建
	void genMesh();
	void changeAngleAutoSignal();
	void generateLinkBodySignal();
	void topoSplitSignal();
	void refineSignal();

	// 渲染
	void groundShowChanged(int arg);
	void rskelShowChanged(int args);
	void axisShowChanged(int args);
	void showTriMeshChanged(int arg);

	//
	void viewFront();
	void eSubdivision();
	bool actSubdivisionNormal();
};
