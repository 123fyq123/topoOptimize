#include "SingalSlot.h"

SingalSlot::SingalSlot(Ui::MainWindow* ui, QWidget* parent) : QWidget(parent) {
	this->ui = ui;
	process_status_ = 0;
	QObject::connect(SingalFather::getInstance(), SIGNAL(supdateGL()), ui->openGLWidget, SLOT(rupdateGL()));
	// 根据solutionSkil网格重建modelModel网格
	QObject::connect(SingalFather::getInstance(), SIGNAL(sdrawNormalSkel()), ui->openGLWidget, SLOT(drawNormalSkel()));
	
	// 
	QObject::connect(SingalFather::getInstance(), SIGNAL(genCoarseMeshFinish()), this, SLOT(genCoarseMeshFinishMesh()));
	// cube info
	QObject::connect(SingalFather::getInstance(), SIGNAL(getData(double, double, double, double, double, double, double, double, double)),
		this, SLOT(reNewDataSignal(double, double, double, double, double, double, double, double, double)));

	QObject::connect(SingalFather::getInstance(), SIGNAL(scheckProcess(int, int, bool&)), this, SLOT(checkProcess(int, int, bool&)));

}

SingalSlot::~SingalSlot(void) {
}


void SingalSlot::reNewDataSignal(double x, double y, double z, double angleX,
	double angleY, double angleZ, double length, double width, double height) {
//	ui->XPosition->setText(QString::number(x, 'f', 4));
//	ui->YPosition->setText(QString::number(y, 'f', 4));
//	ui->ZPosition->setText(QString::number(z, 'f', 4));
//	ui->Length->setText(QString::number(length, 'f', 4));
//	ui->Width->setText(QString::number(width, 'f', 4));
//	ui->Height->setText(QString::number(height, 'f', 4));
//	ui->XRotate->setText(QString::number(angleX, 'f', 4));
//	ui->YRotate->setText(QString::number(angleY, 'f', 4));
//	ui->ZRotate->setText(QString::number(angleZ, 'f', 4));
//	ui->btn_edit_cube->setEnabled(true);
}

void SingalSlot::genCoarseMeshFinishMesh() {
	std::string filePath = "coarse.mesh";
	int module_cnt = ui->openGLWidget->loadMesh(filePath);
	assert(ui->openGLWidget->isVaildHexMesh());
}

void SingalSlot::checkProcess(int pre_status, int now_status, bool& checked) {
	if (pre_status == 0)
		checked = true;
	else {
		if (process_status_ == pre_status) {
			checked = true;		// pass check
		}
		else {
			checked = false;	// dispass check
		}
	}
	if (checked) process_status_ = now_status;
}