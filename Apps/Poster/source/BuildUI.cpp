#include "BuildUI.h"
#include "SingalFather.h"

BuildUI::BuildUI(Ui::MainWindow* ui, QWidget* parent) : QWidget(parent) {
	this->ui = ui;
	//QButtonGroup* group = new QButtonGroup(this);
	//group->addButton(ui->transformMode);
	//group->addButton(ui->editMode);
	//group->addButton(ui->addMode);
	//group->addButton(ui->connectMode);

//	QObject::connect(ui->btn_edit_cube, SIGNAL(clicked()), this, SLOT(editCubeSignal()));
//	QObject::connect(ui->btn_add_cube, SIGNAL(clicked()), this, SLOT(addCubeSignal()));
	QObject::connect(ui->transformMode, SIGNAL(clicked()), this, SLOT(transformCell()));
	QObject::connect(ui->editMode, SIGNAL(clicked()), this, SLOT(editCell()));
	QObject::connect(ui->addMode, SIGNAL(clicked()), this, SLOT(addCell()));
	QObject::connect(ui->connectMode, SIGNAL(clicked()), this, SLOT(connnectCell()));
	
//	// 六面体模型构建
//	QObject::connect(ui->boom, SIGNAL(clicked()), this, SLOT(genMesh()));
//	QObject::connect(ui->boom1, SIGNAL(clicked()), this, SLOT(changeAngleAutoSignal()));
//	QObject::connect(ui->boom2, SIGNAL(clicked()), this, SLOT(generateLinkBodySignal()));
//	QObject::connect(ui->boom3, SIGNAL(clicked()), this, SLOT(topoSplitSignal()));
//	QObject::connect(ui->refine, SIGNAL(clicked()), this, SLOT(refineSignal()));
	
	// 渲染部分
	QObject::connect(ui->groundShow, SIGNAL(stateChanged(int)), this, SLOT(groundShowChanged(int)));
	QObject::connect(ui->axisShow, SIGNAL(stateChanged(int)), this, SLOT(axisShowChanged(int)));
	QObject::connect(ui->skelShow, SIGNAL(stateChanged(int)), this, SLOT(rskelShowChanged(int)));
	QObject::connect(ui->showTriMesh, SIGNAL(stateChanged(int)), this, SLOT(showTriMeshChanged(int)));

	// other
	QObject::connect(ui->viewInit, SIGNAL(clicked()), this, SLOT(viewFront()));
	QObject::connect(ui->subdivedge, SIGNAL(clicked()), this, SLOT(eSubdivision()));
	QObject::connect(ui->actSubdivisionNormal, SIGNAL(clicked()), this, SLOT(actSubdivisionNormal()));
}

BuildUI::~BuildUI(void) {
}

void BuildUI::editCubeSignal() {
//	double x = ui->XPosition->text().toDouble();
//	double y = ui->YPosition->text().toDouble();
//	double z = ui->ZPosition->text().toDouble();
//	double angleX = ui->XRotate->text().toDouble();
//	double angleY = ui->YRotate->text().toDouble();
//	double angleZ = ui->ZRotate->text().toDouble();
//	double length = ui->Length->text().toDouble();
//	double width = ui->Width->text().toDouble();
//	double height = ui->Height->text().toDouble();
//
//	ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
//	ui->openGLWidget->solutionSkel->reditNode(x, y, z, angleX, angleY, angleZ, length, width, height);
}
//
void BuildUI::addCubeSignal() {
//	double x = ui->XPosition->text().toDouble();
//	double y = ui->YPosition->text().toDouble();
//	double z = ui->ZPosition->text().toDouble();
//	double angleX = ui->XRotate->text().toDouble();
//	double angleY = ui->YRotate->text().toDouble();
//	double angleZ = ui->ZRotate->text().toDouble();
//	double length = ui->Length->text().toDouble();
//	double width = ui->Width->text().toDouble();
//	double height = ui->Height->text().toDouble();
//
//	ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
//	ui->openGLWidget->solutionSkel->raddNode(x, y, z, angleX, angleY, angleZ, length, width, height);
}

void BuildUI::transformCell() {
	//ui->openGLWidget->solutionSkel->selectNodes.clear();
	if (ui->transformMode->isChecked()) {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::TRANSFROM;
		if (ui->openGLWidget->solutionSkel->selectNode != nullptr)
			ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::TRANSFROM);
		ui->openGLWidget->SetEditModel(ModelDraw::COORDINATE);
		ui->openGLWidget->setPickType(3);						// 单选cell
	}
	else {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ModelDraw::PICKITEM);
		ui->openGLWidget->makeCurrent();
		update();
	}
	
}

void BuildUI::editCell() {
	//ui->openGLWidget->solutionSkel->selectNodes.clear();
	if (ui->editMode->isChecked()) {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::EDIT;
		if(ui->openGLWidget->solutionSkel->selectNode != nullptr)
			ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::EDIT);
		ui->openGLWidget->SetEditModel(ModelDraw::EditModel::COORDINATE);
		ui->openGLWidget->setPickType(3);						// 单选cell

	}
	else {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ModelDraw::PICKITEM);
		ui->openGLWidget->makeCurrent();
		update();
	}
	
}

void BuildUI::addCell() {
	if (ui->addMode->isChecked()) {
		ui->openGLWidget->initializePick();
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ui->openGLWidget->EditModel::SET_CELL_SKEL);
		ui->openGLWidget->setPickType(3);
	}
	else {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ModelDraw::PICKITEM);
		ui->openGLWidget->makeCurrent();
		update();
	}
}

void BuildUI::connnectCell() {
	if (ui->connectMode->isChecked()) {
		ui->openGLWidget->initializePick();
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ui->openGLWidget->EditModel::CON_CELL_SKEL);
		ui->openGLWidget->setPickType(3);
	}
	else {
		ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
		ui->openGLWidget->SetEditModel(ModelDraw::PICKITEM);
		ui->openGLWidget->makeCurrent();
		update();
	}
}

void BuildUI::selectMutiNodesSignal() {
	//SingalFather::getInstance()->sselectMutiNodes();
	ui->openGLWidget->SetEditModel(ModelDraw::PICK_CELLS_SKEL);
	ui->openGLWidget->setPickType(5);
}

void BuildUI::selectALineSignal() {
	ui->openGLWidget->SetEditModel(ModelDraw::PICK_CELLS_SKEL);
	ui->openGLWidget->setPickType(5);
	//SingalFather::getInstance()->sselectALine();
	ui->openGLWidget->rselectALine();
}



/***************************************************/
void BuildUI::genMesh() {
	/*
	bool checked;
	SingalFather::getInstance()->scheckProcess(1, 5, checked);
	if (!checked)  return;
	*/

	ui->openGLWidget->solutionSkel->pipeline();
}

void BuildUI::changeAngleAutoSignal() {
	/*
	bool checked;
	SingalFather::getInstance()->scheckProcess(1, 2, checked);
	if (checked) {
		return;
	}
	*/
	ui->openGLWidget->solutionSkel->changeAngleAuto();
	/*
	SingalFather::getInstance()->scheckProcess(2, 2, checked);
	if (checked) {
		ui->openGLWidget->solutionSkel->changeAngleAuto();
		return;
	}
	SingalFather::getInstance()->scheckProcess(3, 3, checked);
	if (checked) {
		ui->openGLWidget->solutionSkel->changeAngleAuto();
		return;
	}
	*/
} 

void BuildUI::generateLinkBodySignal() {
	/*
	bool checked;
	SingalFather::getInstance()->scheckProcess(2, 3, checked);
	if (!checked)  return;
	*/
	ui->openGLWidget->solutionSkel->generateLinkBody();
}

void BuildUI::topoSplitSignal() {
	/*
	bool checked;
	SingalFather::getInstance()->scheckProcess(3, 4, checked);
	if (!checked)  return;
	*/

	ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
	if (ui->openGLWidget->solutionSkel->selectNode != nullptr)
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
	ui->openGLWidget->SetEditModel(ModelDraw::EditModel::DEFULTEDIT);
	ui->openGLWidget->solutionSkel->topologicalSplit();
	ui->openGLWidget->SetDrawModule(ModelDraw::DrawModule::HEXMESH);
}

void BuildUI::refineSignal() {
	/*
	bool checked;
	SingalFather::getInstance()->scheckProcess(4, 5, checked);
	if (!checked)  return;
	*/

	ui->openGLWidget->select_coordinater_mode = CoordinateAxis::ModeType::HIDE;
	if (ui->openGLWidget->solutionSkel->selectNode != nullptr)
		ui->openGLWidget->local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
	ui->openGLWidget->SetEditModel(ModelDraw::EditModel::DEFULTEDIT);
	ui->openGLWidget->solutionSkel->refineJN();
}


/********************************************/
void BuildUI::groundShowChanged(int arg) {
	//SingalFather::getInstance()->sgroundShowChanged(arg);
	if (arg) {
		ui->openGLWidget->bshowGround = true;
	}
	else {
		ui->openGLWidget->bshowGround = false;
	}
	ui->openGLWidget->update();
}

void BuildUI::rskelShowChanged(int args) {
	if (args) {
		ui->openGLWidget->bshowSKEL = true;
	}
	else {
		ui->openGLWidget->bshowSKEL = false;
	}
	ui->openGLWidget->update();
}

void BuildUI::axisShowChanged(int args) {
	//SingalFather::getInstance()->saxisShowChanged(args);
	if (args) {
		ui->openGLWidget->bshowAxis = true;
	}
	else {
		ui->openGLWidget->bshowAxis = false;
	}
	ui->openGLWidget->update();
}

void BuildUI::showTriMeshChanged(int arg) {
	if (arg) {
		ui->openGLWidget->bshowTriMesh = true;
	}
	else {
		ui->openGLWidget->bshowTriMesh = false;
	}
	ui->openGLWidget->update();
}

void BuildUI::viewFront() {
	ui->openGLWidget->updateMeshModel();
}

void BuildUI::eSubdivision() {
	ui->openGLWidget->edgeSubdivision(ui->subdivnum->text().toDouble());
}

bool BuildUI::actSubdivisionNormal() {
	ui->openGLWidget->subdivisionNormal();
	//ui->openGLWidget->SetDrawModule(ModelDraw::DrawModule::SKEL);
	return true;
}