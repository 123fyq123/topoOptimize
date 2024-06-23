#pragma once
/*
* 此类对应解决singalFather的信号处理
* singalFather与此类的目的：1. 处理涉及UI的功能为其他类提供接口 2. 供solutionSkel去调用ModelDraw对象功能（涉及到opengl的功能）
* 实现信号接收的槽函数与UI更变
*/

#include <QtWidgets/QtWidgets>
#include "ui_iGame-Poster.h"
#include "SingalFather.h"
//#include "SkeletalModelDraw.h"
class SingalSlot : public QWidget
{
	Q_OBJECT
private:
	Ui::MainWindow* ui;
	int process_status_;
public:
	SingalSlot(Ui::MainWindow* ui, QWidget* parent = Q_NULLPTR);
	~SingalSlot(void);
public slots:

	void reNewDataSignal(double x, double y, double z, double angleX,
		double angleY, double angleZ, double length, double width, double height);

	void genCoarseMeshFinishMesh();

	void checkProcess(int pre_status,int now_status, bool& checked);
};

