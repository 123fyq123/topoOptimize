#include "iGamePoster.h"
#include "Analysis.h"
#include<fstream>
#include<sstream>
#include<iostream>
#include <atomic>
#include <magic_enum.hpp>
#include <QMessageBox>
#pragma execution_character_set("utf-8")
#include <filesystem>
#include "IGA_CC_Helper.h"

namespace {

#define IGAME_STRINGIFY_IMPL(x) #x
#define IGAME_STRINGIFY(x) IGAME_STRINGIFY_IMPL(x)

bool WriteRegularHexMeshFile(const std::filesystem::path& output_path, int nx, int ny, int nz)
{
	if (nx <= 0 || ny <= 0 || nz <= 0) {
		return false;
	}

	std::ofstream out(output_path, std::ios::out);
	if (!out.good()) {
		return false;
	}

	const int vertex_count = (nx + 1) * (ny + 1) * (nz + 1);
	const int cell_count = nx * ny * nz;
	const double half_x = static_cast<double>(nx) / 2.0;
	const double half_y = static_cast<double>(ny) / 2.0;
	const double half_z = static_cast<double>(nz) / 2.0;
	out << "MeshVersionFormatted 1\n";
	out << "Dimension 3\n";
	out << "Vertices\n";
	out << vertex_count << "\n";
	for (int k = 0; k <= nz; ++k) {
		const double z = half_z - static_cast<double>(k);
		for (int j = 0; j <= ny; ++j) {
			const double y = -half_y + static_cast<double>(j);
			for (int i = 0; i <= nx; ++i) {
				const double x = -half_x + static_cast<double>(i);
				out << x << " " << y << " " << z << " -1\n";
			}
		}
	}

	auto vertex_id = [nx, ny](int i, int j, int k) {
		return k * (ny + 1) * (nx + 1) + j * (nx + 1) + i + 1;
	};

	out << "Hexahedra\n";
	out << cell_count << "\n";
	for (int k = 0; k < nz; ++k) {
		for (int j = 0; j < ny; ++j) {
			for (int i = 0; i < nx; ++i) {
				const int v1 = vertex_id(i, j, k);
				const int v2 = vertex_id(i + 1, j, k);
				const int v3 = vertex_id(i + 1, j + 1, k);
				const int v4 = vertex_id(i, j + 1, k);
				const int v5 = vertex_id(i, j, k + 1);
				const int v6 = vertex_id(i + 1, j, k + 1);
				const int v7 = vertex_id(i + 1, j + 1, k + 1);
				const int v8 = vertex_id(i, j + 1, k + 1);
				out << v1 << " " << v2 << " " << v3 << " " << v4 << " "
					<< v5 << " " << v6 << " " << v7 << " " << v8 << " 1\n";
			}
		}
	}
	out << "End\n";

	return out.good();
}

}

iGamePoster::iGamePoster(QWidget* parent) :QMainWindow(parent), ui(new Ui::MainWindow)
{
	Sleep(1000);
	ui->setupUi(this);
	build_ui_ = new BuildUI(ui);
	singal_slot_ = new SingalSlot(ui);
	topopt_refresh_timer_ = new QTimer(this);
	topopt_refresh_timer_->setInterval(80);
	connect(topopt_refresh_timer_, &QTimer::timeout, this, [this]() {
		if (ui == nullptr || ui->openGLWidget == nullptr) return;
		if (ui->openGLWidget->is_need_update_topopt_model) {
			ui->openGLWidget->update();
		}
		});
	topopt_refresh_timer_->start();
	ui->tabWidget->removeTab(5);
	buttonGroup_forPick.setExclusive(true);
	buttonGroup_forPick.addButton(ui->ckbEnablePick);
	//    buttonGroup_forPick.addButton(ui->isRangeSelectRadioButton);
	buttonGroup_forPick.addButton(ui->EndSelectRadioButton);
	// menu bar
	QMenuBar* menu_bar = menuBar();

	ui->toolBox_4->setVisible(true);
	ui->btnTestDetFaces->setVisible(false);
	ui->btnTestResample->setVisible(false);

	{
		QWidget* topopt_page = ui->TopOpt_MMC_type->parentWidget();
		btn_generate_hex_mesh_ = new QPushButton("Generate Mesh", topopt_page);
		btn_generate_hex_mesh_->setObjectName("btnGenerateHexMesh");
		btn_generate_hex_mesh_->setGeometry(20, 305, 231, 27);
	}

	//	ui->tabWidget->removeTab(0);  // 删除Tab，实际效果是隐藏Tab
	//	ui->tabWidget->removeTab(1);  // 删除Tab，实际效果是隐藏Tab
	//    ui->tabWidget->removeTab(ui->tabWidget->indexOf(ui->tab_3));  // 删除Tab，实际效果是隐藏Tab

	{
		QDoubleValidator* validator_2 = new QDoubleValidator();
		validator_2->setRange(0, 10000000); // 用实际的最小值和最大值替换 minValue 和 maxValue
		validator_2->setDecimals(6); // 设置允许的小数位数，例如 2 位
		validator_2->setNotation(QDoubleValidator::StandardNotation);


		QDoubleValidator* validator_0 = new QDoubleValidator();
		validator_0->setRange(0, 1); // 用实际的最小值和最大值替换 minValue 和 maxValue
		validator_0->setDecimals(6); // 设置允许的小数位数，例如 2 位
		validator_0->setNotation(QDoubleValidator::StandardNotation);

		ui->TopOpt_Rate_2->setText("0.3");
		ui->TopOpt_Rate_2->setValidator(validator_0);


		QDoubleValidator* validator_1 = new QDoubleValidator();
		validator_1->setDecimals(6); // 设置允许的小数位数，例如 2 位
		validator_1->setNotation(QDoubleValidator::StandardNotation);
		validator_1->setRange(1, 10000); // 用实际的最小值和最大值替换 minValue 和 maxValue
		ui->TopOpt_Iters_2->setText("20");
		ui->TopOpt_Iters_2->setValidator(validator_2);


		ui->XPosition_fixed->setText("0.0");
		ui->XPosition_fixed->setValidator(validator_2);
		ui->YPosition_fixed->setText("0.0");
		ui->YPosition_fixed->setValidator(validator_2);
		ui->ZPosition_fixed->setText("0.0");
		ui->ZPosition_fixed->setValidator(validator_2);

		ui->XPosition_force->setText("0.0");
		ui->XPosition_force->setValidator(validator_2);
		ui->YPosition_force->setText("0.0");
		ui->YPosition_force->setValidator(validator_2);
		ui->ZPosition_force->setText("0.0");
		ui->ZPosition_force->setValidator(validator_2);

		ui->YangShiMoLiang->setText("100000");
		ui->YangShiMoLiang->setValidator(validator_2);

		ui->PoSongBi->setText("0.3");
		ui->PoSongBi->setValidator(validator_2);

		ui->Thermal_alpha->setText("0.2");
		ui->Thermal_alpha->setValidator(validator_2);

		ui->Thermal_stepTime->setText("0.001");
		ui->Thermal_stepTime->setValidator(validator_2);

		ui->Thermal_SumTime->setText("1.0");
		ui->Thermal_SumTime->setValidator(validator_2);

		QIntValidator* validator_3 = new QIntValidator();
		validator_3->setRange(0, 20);

		ui->Smooth_Times->setText("2");
		ui->Smooth_Times->setValidator(validator_3);
	}
	///************************ skeleton menu ************************/
	/*{
		QMenu* menu_skeleton = menu_bar->addMenu("Skeleton");
		QAction* load_skeleton = menu_skeleton->addAction("Load Skeleton");
		connect(load_skeleton, &QAction::triggered, [=]() {
			LoadSkeleton();
			});
		QAction* load_skeletonLee = menu_skeleton->addAction("Load SkeltonLee");
		connect(load_skeletonLee, &QAction::triggered, [=]() {
			LoadSkeletonLee();
			});

	}*/

	/************************ file ************************/


	/************************ Temp-lee ************************/
	{
		QObject::connect(ui->editbsplinesolid, SIGNAL(clicked(bool)), this, SLOT(EditBsplineSolid()));
		QObject::connect(ui->bssCreator, SIGNAL(clicked(bool)), this, SLOT(auto_ShowBsplineSolid()));
		QObject::connect(ui->bssCreator1, SIGNAL(clicked(bool)), this, SLOT(manually_ShowBsplineSolid()));
		//QObject::connect(ui->bssCreator_2, SIGNAL(clicked(bool)), this, SLOT(manually_ShowBsplineSolid()));
		//QObject::connect(ui->bssCreator_3, SIGNAL(clicked(bool)), this, SLOT(generate_F()));
		//QObject::connect(ui->bssCreator_4, SIGNAL(clicked(bool)), this, SLOT(generate_V()));

		QObject::connect(ui->editbsplinesolid_2, SIGNAL(clicked(bool)), this, SLOT(IF_ShowBsplineSolid()));
		QObject::connect(ui->editbsplinesolid_3, SIGNAL(clicked(bool)), this, SLOT(IF_Real_time_rander()));
		QObject::connect(ui->editbsplinesolid_4, SIGNAL(clicked(bool)), this, SLOT(IF_ShowBsplineSolid_inside()));
		QObject::connect(ui->editbsplinesolid_5, SIGNAL(clicked(bool)), this, SLOT(Choose_Face()));

		///点选部分内容
		QObject::connect(ui->btnSavePickItem, SIGNAL(clicked(bool)), this, SLOT(savePickItem()));
		QObject::connect(ui->btnLoadPickItem, SIGNAL(clicked(bool)), this, SLOT(openPickItem()));
		QObject::connect(ui->btnAddForce, SIGNAL(clicked(bool)), this, SLOT(addPickForce()));
		QObject::connect(ui->btnAddFixed_2, SIGNAL(clicked(bool)), this, SLOT(addPickFixed()));
		QObject::connect(ui->Protect_topopt, SIGNAL(clicked(bool)), this, SLOT(TobeProtected()));
		//QObject::connect(ui->btnAddHyperelasticityTimeSettings, SIGNAL(clicked(bool)), this, SLOT(addHyperelasticityTimeSettings()));
	}


	/************************ file ************************/
	{
		QMenu* menu_file = menu_bar->addMenu("文件");
		menu_file->setTitle("文件");
		QAction* new_model = menu_file->addAction("打开");
		menu_file->addSeparator();
		QAction* save_model = menu_file->addAction("保存");
		menu_file->addSeparator();
		QAction* open_model = menu_file->addAction("新建");
		menu_file->addSeparator();
		QAction* clear_model = menu_file->addAction("清空");
		menu_file->addSeparator();
		QAction* export_model = menu_file->addAction("导出");
		menu_file->addSeparator();
		connect(new_model, &QAction::triggered, [=]() {
			LoadHexahedralMesh();
			ui->rtnCell->setCheckable(true);
			});
		connect(save_model, &QAction::triggered, [=]() {
			SaveHexahedralMesh();
			});
	}

	/************************ 编辑操作 ************************/
	{
		QMenu* menu_edit = menu_bar->addMenu("编辑");
		menu_edit->setTitle("编辑");
		QAction* choose_obj = menu_edit->addAction("选择对象");
		menu_edit->addSeparator();
		QAction* open_model = menu_edit->addAction("删除对象");
		menu_edit->addSeparator();
		QAction* undo_opt = menu_edit->addAction("撤销");
		menu_edit->addSeparator();
	}


	/************************ 视图 ************************/
	{
		QMenu* menu_view = menu_bar->addMenu("视图");
		menu_view->setTitle("视图");
		QAction* zoom_in = menu_view->addAction("放大");
		menu_view->addSeparator();
		QAction* zoom_out = menu_view->addAction("缩小");
		menu_view->addSeparator();
		QAction* reset = menu_view->addAction("默认视角");
		menu_view->addSeparator();
		QAction* cross_section = menu_view->addAction("显示截面");
		menu_view->addSeparator();
		QAction* axial = menu_view->addAction("轴向对齐");
		menu_view->addSeparator();

	}

	/************************ 渲染 ************************/
	{
		QMenu* menu_render = menu_bar->addMenu("渲染");
		menu_render->setTitle("渲染");
		QAction* point_mode = menu_render->addAction("点云模式");
		menu_render->addSeparator();
		QAction* line_mode = menu_render->addAction("线框模式");
		menu_render->addSeparator();
		QAction* face_mode = menu_render->addAction("面片模式");
		menu_render->addSeparator();
		QAction* volume_mode = menu_render->addAction("Volume模式");
		volume_mode->setText("\xE4\xBD\x93\xE6\xA8\xA1\xE5\xBC\x8F");
		menu_render->addSeparator();
		QAction* spline_mode = menu_render->addAction("Spline模式");
		spline_mode->setText("\xE6\xA0\xB7\xE6\x9D\xA1\xE4\xBD\x93\xE6\xA8\xA1\xE5\xBC\x8F");
		menu_render->addSeparator();
		QAction* cross_section = menu_render->addAction("显示截面");
		menu_render->addSeparator();
		QAction* axial = menu_render->addAction("轴向对齐");
		menu_render->addSeparator();

	}

	/************************ 网格************************/
	{
		QMenu* menu_surface_mesh = menu_bar->addMenu("网格");
		//menu_surface_mesh->setTitle("网格");
		menu_surface_mesh->setTitle(QString("\xE6\x9B\xB2\xE9\x9D\xA2\xE7\xBD\x91\xE6\xA0\xBC"));// 曲面网格
		//menu_surface_mesh->setIcon(QIcon("./Resources/triangle.avg"));

		QAction* load_surface_mesh = menu_surface_mesh->addAction("Load Surface Mesh");
		load_surface_mesh->setText(QString("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B"));
		QAction* save_surface_mesh = menu_surface_mesh->addAction("Save Surface Mesh");
		save_surface_mesh->setText(QString("\xE4\xBF\x9D\xE5\xAD\x98\xE6\xA8\xA1\xE5\x9E\x8B"));
		menu_surface_mesh->addSeparator();
		//QAction* remeshing_adaptive_surface_mesh = menu_surface_mesh->addAction("Remeshing(Adaptive Edge Length)");
		//remeshing_adaptive_surface_mesh->setText(QString("\xE9\x87\x8D\xE6\x96\xB0\xE7\xBD\x91\xE6\xA0\xBC\xE5\x8C\x96\xEF\xBC\x88\xE8\x87\xAA\xE9\x80\x82\xE5\xBA\x94\xEF\xBC\x89"));
		QAction* remeshing_uniform_surface_mesh = menu_surface_mesh->addAction("Remeshing(Uniform Edge Length)");
		remeshing_uniform_surface_mesh->setText(QString("\xE9\x87\x8D\xE6\x96\xB0\xE7\xBD\x91\xE6\xA0\xBC\xE5\x8C\x96\xEF\xBC\x88\xE5\x9D\x87\xE5\x8C\x80\xEF\xBC\x89"));
		menu_surface_mesh->addSeparator();
		QAction* simplification_surface_mesh = menu_surface_mesh->addAction("Simplification(QEM)");
		simplification_surface_mesh->setText(QString("\xE7\xBD\x91\xE6\xA0\xBC\xE7\xAE\x80\xE5\x8C\x96"));
		/*QAction* simplification_powered_surface_mesh = menu_surface_mesh->addAction("Simplification(Powered QEM)");
		simplification_powered_surface_mesh->setText(QString("\xE7\xBD\x91\xE6\xA0\xBC\xE7\xAE\x80\xE5\x8C\x96\xEF\xBC\x88\xE8\xB4\xA8\xE9\x87\x8F\xE6\x95\x8F\xE6\x84\x9F\xEF\xBC\x89"));*/
		menu_surface_mesh->addSeparator();
		QAction* denoising_surface_mesh = menu_surface_mesh->addAction("Denoising");
		denoising_surface_mesh->setText(QString("\xE5\x8E\xBB\xE5\x99\xAA"));
		menu_surface_mesh->addSeparator();
		QAction* subdivision_surface_mesh = menu_surface_mesh->addAction("Subdivision");
		subdivision_surface_mesh->setText(QString("\xE7\xBB\x86\xE5\x88\x86"));
		menu_surface_mesh->addSeparator();
		QAction* parameterization_surface_mesh = menu_surface_mesh->addAction("Parameterization");
		parameterization_surface_mesh->setText(QString("\xE5\x8F\x82\xE6\x95\xB0\xE5\x8C\x96"));
		menu_surface_mesh->addSeparator();
		QAction* subdivision_hexahedral_mesh = menu_surface_mesh->addAction("Parameterization");
		subdivision_hexahedral_mesh->setText(QString("六面细分"));


		connect(subdivision_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input subdivision times",
				1, 1, 10, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->subdivision(iter);
			int spend = clock() - t1;
			QString str = "Subdivision success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});


		connect(load_surface_mesh, &QAction::triggered, [=]() {
			LoadSurfaceMesh();
			ui->rtnCell->setCheckable(false);
			});

		connect(save_surface_mesh, &QAction::triggered, [=]() {
			SaveSurfaceMesh();
			});


		connect(remeshing_uniform_surface_mesh, &QAction::triggered, [=]() {
			time_t t1 = clock();
			ui->openGLWidget->remeshing_uniform();
			int spend = clock() - t1;
			QString str = "Remeshing(uniform) success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(denoising_surface_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Surface Mesh", "Please input denoising times",
				1, 1, 15, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->denoising_bnf(iter);
			int spend = clock() - t1;
			QString str = "Denoising success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(subdivision_surface_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Surface Mesh", "Please input subdivision times",
				1, 0, 5, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->subdivision_loop(iter);
			int spend = clock() - t1;
			QString str = "Subdivision success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(simplification_surface_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int vcnt = QInputDialog::getInt(this, "Surface Mesh", "Please input target vertices size",
				ui->openGLWidget->getVerticesSize(), 3, ui->openGLWidget->getVerticesSize(), 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->simplification_qem(vcnt);
			int spend = clock() - t1;
			QString str = "Simplification success, current vertices size: " + QString::number(ui->openGLWidget->getVerticesSize()) +
				"cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(parameterization_surface_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int mode = QInputDialog::getInt(this, "Surface Mesh", "Please select boundary: 1 Squard, 2 Circle",
				2, 1, 2, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->parameterization(mode);
			int spend = clock() - t1;
			QString str = "Parameterization success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});
	}

	/************************ skeletal mesh ************************/
	{
		QMenu* menu_skeletal_mesh = menu_bar->addMenu("骨架网格");
		menu_skeletal_mesh->setTitle("骨架网格");
		//menu_skeletal_mesh->menuAction()->setVisible(false);
		//menu_skeletal_mesh->setTitle("\xE9\xAA\xA8\xE6\x9E\xB6\xE7\xBD\x91\xE6\xA0\xBC");
		QAction* load_skeletal_mesh = menu_skeletal_mesh->addAction("\xE5\x8A\xA0\xE8\xBD\xBD\xE9\xAA\xA8\xE6\x9E\xB6\xE7\xBD\x91\xE6\xA0\xBC");
		QAction* save_skeletal_mesh = menu_skeletal_mesh->addAction("\xE4\xBF\x9D\xE5\xAD\x98\xE9\xAA\xA8\xE6\x9E\xB6\xE7\xBD\x91\xE6\xA0\xBC");
		menu_skeletal_mesh->addSeparator();
		QAction* simplification = menu_skeletal_mesh->addAction("\xE7\xAE\x80\xE5\x8C\x96");
		menu_skeletal_mesh->addSeparator();
		QAction* subdivision = menu_skeletal_mesh->addAction("细分");
		menu_skeletal_mesh->addSeparator();
		QAction* update_all_handles = menu_skeletal_mesh->addAction("更新节点");
		menu_skeletal_mesh->addSeparator();
		QAction* change_show_tangent_plane = menu_skeletal_mesh->addAction("\xE6\x98\xBE\xE7\xA4\xBA\xE5\xA4\x96\xE5\x88\x87\xE7\xAB\x8B\xE6\x96\xB9\xE4\xBD\x93");
		connect(load_skeletal_mesh, &QAction::triggered, [=]() {
			LoadSkeletalMesh();
			});
		connect(save_skeletal_mesh, &QAction::triggered, [=]() {
			SaveSkeletalMesh();
			});
		connect(update_all_handles, &QAction::triggered, [=]() {
			ui->openGLWidget->updateAllHandles();
			});
		connect(simplification, &QAction::triggered, [=]() {
			bool input_ok = false;
			size_t target_vcnt = QInputDialog::getInt(this, "Skeletal Mesh", "Please input target vertices size",
				ui->openGLWidget->getVerticesSize(), 3, ui->openGLWidget->getVerticesSize(), 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->simplification_qmat(target_vcnt);
			time_t t2 = clock();
			QString str = "Simplification success, cost " + QString::number(int(t2 - t1)) + "ms";
			ui->info->printString(str);

			});
		connect(subdivision, &QAction::triggered, [=]() {
			bool input_ok = false;
			size_t iter = QInputDialog::getInt(this, "Skeletal Mesh", "Please input subdivision times",
				2, 1, 5, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->subdivision_skl(iter);
			time_t t2 = clock();
			QString str = "Subdivision success, cost " + QString::number(int(t2 - t1)) + "ms";
			ui->info->printString(str);

			});
		connect(change_show_tangent_plane, &QAction::triggered, [=]() {
			ui->openGLWidget->show_tangent_plane = !ui->openGLWidget->show_tangent_plane;
			if (ui->openGLWidget->show_tangent_plane) {
				ui->openGLWidget->updateDrawData();
				change_show_tangent_plane->setText("\xE9\x9A\x90\xE8\x97\x8F\xE5\xA4\x96\xE5\x88\x87\xE7\xAB\x8B\xE6\x96\xB9\xE4\xBD\x93");

			}
			else {
				ui->openGLWidget->updateDrawData();
				change_show_tangent_plane->setText("\xE6\x98\xBE\xE7\xA4\xBA\xE5\xA4\x96\xE5\x88\x87\xE7\xAB\x8B\xE6\x96\xB9\xE4\xBD\x93");
			}
			});
	}

	//	/************************ surface mesh menu ************************/
	//	{
	//		QMenu* menu_surface_mesh = menu_bar->addMenu("Surface Mesh");
	//		menu_surface_mesh->menuAction()->setVisible(false);
	//		//QTextCodec* codec = QTextCodec::codecForName("UTF-8");
	//		menu_surface_mesh->setTitle(QString("\xE6\x9B\xB2\xE9\x9D\xA2\xE7\xBD\x91\xE6\xA0\xBC"));// 曲面网格
	//		//menu_surface_mesh->setIcon(QIcon("./Resources/triangle.avg"));
	//
	//		QAction* load_surface_mesh = menu_surface_mesh->addAction("Load Surface Mesh");
	//		load_surface_mesh->setText(QString("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B"));
	//		QAction* save_surface_mesh = menu_surface_mesh->addAction("Save Surface Mesh");
	//		save_surface_mesh->setText(QString("\xE4\xBF\x9D\xE5\xAD\x98\xE6\xA8\xA1\xE5\x9E\x8B"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* remeshing_adaptive_surface_mesh = menu_surface_mesh->addAction("Remeshing(Adaptive Edge Length)");
	//		remeshing_adaptive_surface_mesh->setText(QString("\xE9\x87\x8D\xE6\x96\xB0\xE7\xBD\x91\xE6\xA0\xBC\xE5\x8C\x96\xEF\xBC\x88\xE8\x87\xAA\xE9\x80\x82\xE5\xBA\x94\xEF\xBC\x89"));
	//		QAction* remeshing_uniform_surface_mesh = menu_surface_mesh->addAction("Remeshing(Uniform Edge Length)");
	//		remeshing_uniform_surface_mesh->setText(QString("\xE9\x87\x8D\xE6\x96\xB0\xE7\xBD\x91\xE6\xA0\xBC\xE5\x8C\x96\xEF\xBC\x88\xE5\x9D\x87\xE5\x8C\x80\xEF\xBC\x89"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* simplification_surface_mesh = menu_surface_mesh->addAction("Simplification(QEM)");
	//		simplification_surface_mesh->setText(QString("\xE7\xBD\x91\xE6\xA0\xBC\xE7\xAE\x80\xE5\x8C\x96"));
	//		QAction* simplification_powered_surface_mesh = menu_surface_mesh->addAction("Simplification(Powered QEM)");
	//		simplification_powered_surface_mesh->setText(QString("\xE7\xBD\x91\xE6\xA0\xBC\xE7\xAE\x80\xE5\x8C\x96\xEF\xBC\x88\xE8\xB4\xA8\xE9\x87\x8F\xE6\x95\x8F\xE6\x84\x9F\xEF\xBC\x89"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* denoising_surface_mesh = menu_surface_mesh->addAction("Denoising");
	//		denoising_surface_mesh->setText(QString("\xE5\x8E\xBB\xE5\x99\xAA"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* subdivision_surface_mesh = menu_surface_mesh->addAction("Subdivision");
	//		subdivision_surface_mesh->setText(QString("\xE7\xBB\x86\xE5\x88\x86"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* parameterization_surface_mesh = menu_surface_mesh->addAction("Parameterization");
	//		parameterization_surface_mesh->setText(QString("\xE5\x8F\x82\xE6\x95\xB0\xE5\x8C\x96"));
	//		menu_surface_mesh->addSeparator();
	//		QAction* enter_analysis = menu_surface_mesh->addAction("Exter Analysis Mode");
	//		enter_analysis->setText("\xE8\xBF\x9B\xE5\x85\xA5\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");
	//		QAction* exit_analysis = menu_surface_mesh->addAction("Exter Analysis Mode");
	//		exit_analysis->setText("\xE9\x80\x80\xE5\x87\xBA\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");
	//		menu_surface_mesh->addSeparator();
	//
	//		QAction* get_volume = menu_surface_mesh->addAction("Get Volume");
	//		QAction* test_action = menu_surface_mesh->addAction("Test");
	//		QAction* add_noise = menu_surface_mesh->addAction("Add Gauss Noise");
	//
	//
	//
	//		connect(add_noise, &QAction::triggered, [=]() {
	//			bool input_ok = true;
	//			double max_error = QInputDialog::getDouble(this, "Triangle Mesh", "Please input noise ratio",
	//				0.5f, 0.1f, 1.f, 2, &input_ok);
	//			if (!input_ok) return;
	//			ui->openGLWidget->add_noise(max_error);
	//			});
	//
	//		connect(test_action, &QAction::triggered, [=]() {
	//			bool input_ok = true;
	//			double alpha = QInputDialog::getDouble(this, "Alpha", "",
	//				0.f, 0, 1.f, 2, &input_ok);
	//			if (!input_ok) return;
	//			ui->openGLWidget->test_action(alpha);
	//			});
	//
	//		connect(get_volume, &QAction::triggered, [=]() {
	//			double volume = ui->openGLWidget->getVolume();
	//			QString str = "Current volume: " + QString::number(volume) + ".";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(load_surface_mesh, &QAction::triggered, [=]() {
	//			LoadSurfaceMesh();
	//			ui->rtnCell->setCheckable(false);
	//			});
	//
	//		connect(save_surface_mesh, &QAction::triggered, [=]() {
	//			SaveSurfaceMesh();
	//			});
	//
	//		connect(remeshing_adaptive_surface_mesh, &QAction::triggered, [=]() {
	//			time_t t1 = clock();
	//			ui->openGLWidget->remeshing_adaptive();
	//			int spend = clock() - t1;
	//			QString str = "Remeshing(adaptive) success, cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(remeshing_uniform_surface_mesh, &QAction::triggered, [=]() {
	//			time_t t1 = clock();
	//			ui->openGLWidget->remeshing_uniform();
	//			int spend = clock() - t1;
	//			QString str = "Remeshing(uniform) success, cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(denoising_surface_mesh, &QAction::triggered, [=]() {
	//			bool input_ok = false;
	//			int iter = QInputDialog::getInt(this, "Surface Mesh", "Please input denoising times",
	//				1, 1, 15, 1, &input_ok);
	//			if (!input_ok) return;
	//			time_t t1 = clock();
	//			ui->openGLWidget->denoising_bnf(iter);
	//			int spend = clock() - t1;
	//			QString str = "Denoising success, cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(subdivision_surface_mesh, &QAction::triggered, [=]() {
	//			bool input_ok = false;
	//			int iter = QInputDialog::getInt(this, "Surface Mesh", "Please input subdivision times",
	//				1, 0, 5, 1, &input_ok);
	//			if (!input_ok) return;
	//			time_t t1 = clock();
	//			ui->openGLWidget->subdivision_loop(iter);
	//			int spend = clock() - t1;
	//			QString str = "Subdivision success, cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(simplification_surface_mesh, &QAction::triggered, [=]() {
	//			bool input_ok = false;
	//			int vcnt = QInputDialog::getInt(this, "Surface Mesh", "Please input target vertices size",
	//				ui->openGLWidget->getVerticesSize(), 3, ui->openGLWidget->getVerticesSize(), 1, &input_ok);
	//			if (!input_ok) return;
	//			time_t t1 = clock();
	//			ui->openGLWidget->simplification_qem(vcnt);
	//			int spend = clock() - t1;
	//			QString str = "Simplification success, current vertices size: " + QString::number(ui->openGLWidget->getVerticesSize()) +
	//				"cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(simplification_powered_surface_mesh, &QAction::triggered, [=]() {
	//			bool input_ok = false;
	//			int vcnt = QInputDialog::getInt(this, "Surface Mesh", "Please input target vertices size",
	//				ui->openGLWidget->getVerticesSize(), 3, ui->openGLWidget->getVerticesSize(), 1, &input_ok);
	//			if (!input_ok) return;
	//			time_t t1 = clock();
	//			ui->openGLWidget->simplification_qem(vcnt, true);
	//			int spend = clock() - t1;
	//			QString str = "Simplification(Powered) success, current vertices size: " + QString::number(ui->openGLWidget->getVerticesSize()) +
	//				"cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(parameterization_surface_mesh, &QAction::triggered, [=]() {
	//			bool input_ok = false;
	//			int mode = QInputDialog::getInt(this, "Surface Mesh", "Please select boundary: 1 Squard, 2 Circle",
	//				2, 1, 2, 1, &input_ok);
	//			if (!input_ok) return;
	//			time_t t1 = clock();
	//			ui->openGLWidget->parameterization(mode);
	//			int spend = clock() - t1;
	//			QString str = "Parameterization success, cost: " + QString::number(spend) + "ms";
	//			ui->info->printString(str);
	//			});
	//
	//		connect(enter_analysis, &QAction::triggered, [=]() {
	//				bool input_ok = false;
	//				int mode = QInputDialog::getInt(this, "Surface Mesh", "Please select data: 1 Angle, 2 Curvature",
	//					2, 1, 2, 1, &input_ok);
	//				if (!input_ok) return;
	//				ui->openGLWidget->analysis_flag = true;
	//				ui->openGLWidget->analysis_mode = mode;
	//				ui->openGLWidget->updateDrawData();
	//			});
	//
	//		connect(exit_analysis, &QAction::triggered, [=]() {
	//			ui->openGLWidget->analysis_flag = false;
	//			ui->openGLWidget->updateDrawData();
	//			});
	//
	//	}

		/************************ tetrahedron mesh menu ************************/
	{
		QMenu* menu_tetrahedron_mesh = menu_bar->addMenu("Tet Mesh");
		menu_tetrahedron_mesh->menuAction()->setVisible(false);
		menu_tetrahedron_mesh->setTitle("\xE5\x9B\x9B\xE9\x9D\xA2\xE4\xBD\x93\xE7\xBD\x91\xE6\xA0\xBC");
		QAction* load_tetrahedron_mesh = menu_tetrahedron_mesh->addAction("Load Tetrahedral Mesh");
		load_tetrahedron_mesh->setText("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B");
		QAction* save_tetrahedron_mesh = menu_tetrahedron_mesh->addAction("Save Tetrahedral Mesh");
		save_tetrahedron_mesh->setText("\xE4\xBF\x9D\xE5\xAD\x98\xE6\xA8\xA1\xE5\x9E\x8B");
		connect(load_tetrahedron_mesh, &QAction::triggered, [=]() {
			// 读函数还有问题
			LoadTetrahedronMesh();
			ui->rtnCell->setCheckable(true);
			});

		connect(save_tetrahedron_mesh, &QAction::triggered, [=]() {
			SaveTetrahedronMesh();
			});
	}


	/************************ BSpline surface menu ************************/
	{
		QMenu* menu_bspline_surface = menu_bar->addMenu("BSpline Surface");
		menu_bspline_surface->menuAction()->setVisible(false);
		menu_bspline_surface->setTitle("\xE6\xA0\xB7\xE6\x9D\xA1\xE6\x9B\xB2\xE9\x9D\xA2");
		QAction* load_bspline_surface = menu_bspline_surface->addAction("Load BSpline Surface");
		load_bspline_surface->setText("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B");
		QAction* save_bspline_surface = menu_bspline_surface->addAction("Save BSpline Surface");
		save_bspline_surface->setText("\xE4\xBF\x9D\xE5\xAD\x98\xE6\xA8\xA1\xE5\x9E\x8B");
		QAction* test_bspline_surface = menu_bspline_surface->addAction("Load Surface From Patch Test");
		connect(load_bspline_surface, &QAction::triggered, [=]() {
			OpenIGESSurface();
			});
		connect(save_bspline_surface, &QAction::triggered, [=]() {
			SaveIGESFile();
			});
		connect(test_bspline_surface, &QAction::triggered, [=]() {
			ui->openGLWidget->test_bspline();
			});

	}

	/************************ Sparse Optimization menu ************************/


	{
		QMenu* menu_sparse_optimization = menu_bar->addMenu("Sparse Optimization");
		menu_sparse_optimization->menuAction()->setVisible(false);
		menu_sparse_optimization->setTitle("\x4D\x4D\x43\xE8\xBD\xAC\xE6\x8D\xA2\x0D\x0A");
		/*QAction* load_poly_mesh = menu_sparse_optimization->addAction("Load Poly Model");
		load_poly_mesh->setText("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B");*/
		QAction* build_skl = menu_sparse_optimization->addAction("Build Skeleton");
		build_skl->setText("\xE7\x94\x9F\xE6\x88\x90\xE9\xAA\xA8\xE6\x9E\xB6");
		QAction* sparse_opt = menu_sparse_optimization->addAction("Sparse Optimization");
		sparse_opt->setText("\xE7\xA8\x80\xE7\x96\x8F\xE4\xBC\x98\xE5\x8C\x96");
		QAction* generate_axis = menu_sparse_optimization->addAction("Generate Symmetry Axis");
		generate_axis->setText("\xE7\x94\x9F\xE6\x88\x90\xE5\xAF\xB9\xE7\xA7\xB0\xE8\xBD\xB4");
		QAction* shape_approx = menu_sparse_optimization->addAction("Shape Approximation");
		shape_approx->setText("\xE5\xBD\xA2\xE7\x8A\xB6\xE9\x80\xBC\xE8\xBF\x91");

		/*connect(load_poly_mesh, &QAction::triggered, [=]() {
			QString filename = QFileDialog::getOpenFileName(this, tr("Load Obj file"), tr(""), tr("obj Files(*.obj)"));
			poster_filename = std::string((const char*)filename.toLocal8Bit());
			ui->openGLWidget->loadPolyMesh(poster_filename);
			ui->rtnCell->setCheckable(false);
			});*/

		connect(build_skl, &QAction::triggered, [=]() {
			//QString filename = QFileDialog::getOpenFileName(this, tr("Load Obj file"), tr(""), tr("obj Files(*.obj)"));
			QString filename = QFileDialog::getOpenFileName(this, tr("Load Obj file"), tr(""), tr("obj Files(*.obj)"));
			poster_filename = std::string((const char*)filename.toLocal8Bit());
			ui->rtnCell->setCheckable(false);
			bool input_ok = false;
			uint32_t width = QInputDialog::getInt(this, "Poly Mesh", "Please input width",
				722, 0, 2048, 1, &input_ok);
			if (!input_ok) return;
			uint32_t height = QInputDialog::getInt(this, "Poly Mesh", "Please input height",
				122, 0, 2048, 1, &input_ok);
			if (!input_ok) return;
			ui->openGLWidget->BuildSkeletonToFile(poster_filename, width, height);
			});
		connect(sparse_opt, &QAction::triggered, [=]() {
			bool input_ok = false;
			double max_error = QInputDialog::getDouble(this, "Poly Mesh", "Please input max error",
				5.f, 1.f, 10.f, 2, &input_ok);
			if (input_ok) ui->openGLWidget->SparseOptimization(max_error);
			});
		connect(generate_axis, &QAction::triggered, [=]() {
			ui->openGLWidget->GenerateSymAxis();
			});
		connect(shape_approx, &QAction::triggered, [=]() {
			bool input_ok = false;
			uint32_t iter = QInputDialog::getInt(this, "Poly Mesh", "Please input max_iter",
				4000, 0, 1600000, 500, &input_ok);
			if (!input_ok) return;
			ui->openGLWidget->ShapeApproximation(iter);
			});

	}



	/************************ Top Optimization menu ************************/
	{
		QMenu* menu_top_optimization = menu_bar->addMenu("Top Optimization");
		menu_top_optimization->menuAction()->setVisible(false);
		menu_top_optimization->setTitle("\xE6\x8B\x93\xE6\x89\x91\xE4\xBC\x98\xE5\x8C\x96");
		QAction* load_hex_mesh = menu_top_optimization->addAction("Load Mesh");
		load_hex_mesh->setText("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B");
		menu_top_optimization->addSeparator();
		QMenu* menu_analysis = menu_top_optimization->addMenu("Analysis Menu");
		menu_analysis->setTitle("\xE7\xAD\x89\xE5\x87\xA0\xE4\xBD\x95\xE5\x88\x86\xE6\x9E\x90\x0D\x0A");
		QAction* enter_analysis = menu_analysis->addAction("Enter Analysis Mode");
		enter_analysis->setText("\xE8\xBF\x9B\xE5\x85\xA5\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");
		QAction* exit_analysis = menu_analysis->addAction("Exit Analysis Mode");
		exit_analysis->setText("\xE9\x80\x80\xE5\x87\xBA\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");
		menu_analysis->addSeparator();
		QAction* show_x = menu_analysis->addAction("Show Result X");
		show_x->setText("\xE6\x98\xBE\xE7\xA4\xBA\x58\xE8\xBD\xB4\xE7\xBB\x93\xE6\x9E\x9C");
		QAction* show_y = menu_analysis->addAction("Show Result Y");
		show_y->setText("\xE6\x98\xBE\xE7\xA4\xBA\x59\xE8\xBD\xB4\xE7\xBB\x93\xE6\x9E\x9C");
		QAction* show_z = menu_analysis->addAction("Show Result Z");
		show_z->setText("\xE6\x98\xBE\xE7\xA4\xBA\x5A\xE8\xBD\xB4\xE7\xBB\x93\xE6\x9E\x9C");
		menu_top_optimization->addSeparator();
		QAction* top_opt = menu_top_optimization->addAction("Execute");
		top_opt->setText("\xE6\x89\xA7\xE8\xA1\x8C\xE4\xBC\x98\xE5\x8C\x96");


		connect(load_hex_mesh, &QAction::triggered, [=]() {
			LoadHexahedralMesh();
			ui->rtnCell->setCheckable(true);
			});

		connect(enter_analysis, &QAction::triggered, [=]() {
			QtConcurrent::run([&]() {
				/*ui->openGLWidget->updatePickPara();
				std::cout << "Current: " << poster_filename << std::endl;
				std::cout << ui->openGLWidget->pickedge.size() << " " << ui->openGLWidget->pickface.size() << std::endl;
				vector<double> ax, ay, az;
				Analysis(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, ax, ay, az);
				ui->openGLWidget->analysis_data[0] = ax, ui->openGLWidget->analysis_data[1] = ay, ui->openGLWidget->analysis_data[2] = az;
				std::cout << "Analysis success!" << std::endl;
				ui->openGLWidget->analysis_flag = true;
				ui->openGLWidget->analysis_mode = 1;
				ui->openGLWidget->updateDrawData();*/
				});
			});

		connect(exit_analysis, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_flag = false;
			ui->openGLWidget->updateDrawData();
			});

		connect(top_opt, &QAction::triggered, [=]() {
			QtConcurrent::run([&]() {
				//ui->openGLWidget->updatePickPara();
				//std::cout << "current: " << poster_filename << std::endl;
				//std::cout << ui->openGLWidget->pickedge.size() << " " << ui->openGLWidget->pickface.size() << std::endl;
				//TopOpt(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, vv);
				});
			});

		connect(show_x, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_mode = 0;
			ui->openGLWidget->updateDrawData();
			});

		connect(show_y, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_mode = 1;
			ui->openGLWidget->updateDrawData();
			});

		connect(show_z, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_mode = 2;
			ui->openGLWidget->updateDrawData();
			});

	}

	/************************ hexahedron mesh menu ************************/
	{
		QMenu* menu_hexahedral_mesh = menu_bar->addMenu("Hex Mesh");
		menu_hexahedral_mesh->menuAction()->setVisible(false);
		menu_hexahedral_mesh->setTitle(QString("\xE6\x8B\x93\xE6\x89\x91\xE4\xBC\x98\xE5\x8C\x96\xE5\x90\x8E\xE5\xA4\x84\xE7\x90\x86"));
		QAction* load_hexahedral_mesh = menu_hexahedral_mesh->addAction("Load Hexahedral Mesh");
		load_hexahedral_mesh->setText(QString("\xE5\x8A\xA0\xE8\xBD\xBD\xE6\xA8\xA1\xE5\x9E\x8B"));
		QAction* save_hexahedral_mesh = menu_hexahedral_mesh->addAction("Save Hexahedral Mesh");
		save_hexahedral_mesh->setText(QString("\xE4\xBF\x9D\xE5\xAD\x98\xE6\xA8\xA1\xE5\x9E\x8B"));
		menu_hexahedral_mesh->addSeparator();
		QAction* regularization_hexahedral_mesh = menu_hexahedral_mesh->addAction("Regularization");
		regularization_hexahedral_mesh->setText(QString("\xE8\xA7\x84\xE6\x95\xB4\xE5\x8C\x96"));
		QAction* regularization_faster_hexahedral_mesh = menu_hexahedral_mesh->addAction("Regularization(Faster)");
		regularization_faster_hexahedral_mesh->setText(QString("\xE8\xA7\x84\xE6\x95\xB4\xE5\x8C\x96\xEF\xBC\x88\xE5\xA4\x96\xE9\x83\xA8\xEF\xBC\x89"));
		menu_hexahedral_mesh->addSeparator();
		QAction* smoothing_hexahedral_mesh = menu_hexahedral_mesh->addAction("Smoothing");
		smoothing_hexahedral_mesh->setText("\xE5\x85\x89\xE9\xA1\xBA");
		menu_hexahedral_mesh->addSeparator();
		QAction* subdivision_hexahedral_mesh = menu_hexahedral_mesh->addAction("Subdivision");
		subdivision_hexahedral_mesh->setText("\xE7\xBB\x86\xE5\x88\x86");
		menu_hexahedral_mesh->addSeparator();
		QAction* read_force_hexahedral_mesh = menu_hexahedral_mesh->addAction("Read Force Data");
		QAction* read_spc_hexahedral_mesh = menu_hexahedral_mesh->addAction("Read SPC Data");
		QAction* delete_island_hexahedral_mesh = menu_hexahedral_mesh->addAction("Delete island");
		delete_island_hexahedral_mesh->setText("\xE5\x88\xA0\xE9\x99\xA4\xE5\xAD\xA4\xE7\xAB\x8B\xE5\xB2\x9B");
		QAction* delete_redundancy_mesh = menu_hexahedral_mesh->addAction("Delete redundancy cell");
		QAction* write_force_hexahedral_mesh = menu_hexahedral_mesh->addAction("Write Force Data");
		QAction* write_spc_hexahedral_mesh = menu_hexahedral_mesh->addAction("Write SPC Data");
		QAction* remove_constraint_hexahedral_mesh = menu_hexahedral_mesh->addAction("Remove all constraints");
		menu_hexahedral_mesh->addSeparator();
		QAction* reconstruction_hexahedral_mesh = menu_hexahedral_mesh->addAction("Reconstruction");
		reconstruction_hexahedral_mesh->setText("\xE9\x87\x8D\xE5\xBB\xBA");
		//QAction* reconstruction_ske_hexahedral_mesh = menu_hexahedral_mesh->addAction("Reconstruction(Skeleton)");
		menu_hexahedral_mesh->addSeparator();
		QAction* add_noise_hexahedral_mesh = menu_hexahedral_mesh->addAction("Add Noise");
		add_noise_hexahedral_mesh->setText("\xE5\x8A\xA0\xE5\x99\xAA");
		menu_hexahedral_mesh->addSeparator();
		QAction* enter_analysis_hexahedral_mesh = menu_hexahedral_mesh->addAction("Enter Analysis Mode");
		enter_analysis_hexahedral_mesh->setText("\xE8\xBF\x9B\xE5\x85\xA5\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");
		QAction* exit_analysis_hexahedral_mesh = menu_hexahedral_mesh->addAction("Exit Analysis Mode");
		exit_analysis_hexahedral_mesh->setText("\xE9\x80\x80\xE5\x87\xBA\xE5\x88\x86\xE6\x9E\x90\xE6\xA8\xA1\xE5\xBC\x8F");

		menu_hexahedral_mesh->addSeparator();
		QAction* load_compare_mesh = menu_hexahedral_mesh->addAction("Load Compared Mesh");
		QAction* clear_compare_mesh = menu_hexahedral_mesh->addAction("Clear Compared Mesh");

		menu_hexahedral_mesh->addSeparator();
		QAction* update_handles = menu_hexahedral_mesh->addAction("Update All Handles");
		update_handles->setText("\xE6\x9B\xB4\xE6\x96\xB0\xE5\x8F\xA5\xE6\x9F\x84");

		QAction* get_volume = menu_hexahedral_mesh->addAction("Get Volume");

		connect(get_volume, &QAction::triggered, [=]() {
			double volume = ui->openGLWidget->getVolume();
			QString str = "Current volume: " + QString::number(volume) + ".";
			ui->info->printString(str);
			});

		connect(update_handles, &QAction::triggered, [=]() {
			if (ui->openGLWidget->updateAllHandles()) {
				ui->info->printString("update all handles success");
			}
			else {
				ui->info->printString("update all handles fail");
			}
			});

		connect(enter_analysis_hexahedral_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_flag = true;
			ui->openGLWidget->analysis_mode = 3;
			ui->openGLWidget->updateDrawData();
			});

		connect(exit_analysis_hexahedral_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->analysis_flag = false;
			ui->openGLWidget->updateDrawData();
			});

		connect(load_hexahedral_mesh, &QAction::triggered, [=]() {
			LoadHexahedralMesh();
			ui->rtnCell->setCheckable(true);
			});

		connect(load_compare_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->compare_flag = true;
			LoadHexahedralMesh();
			ui->openGLWidget->compare_flag = false;
			QString str = "Read compared vertices data success.";
			ui->info->printString(str);
			});

		connect(clear_compare_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->compared_vertices_data.clear();
			ui->openGLWidget->updateDrawData();
			QString str = "Clear compared vertices data success.";
			ui->info->printString(str);
			});

		connect(save_hexahedral_mesh, &QAction::triggered, [=]() {
			SaveHexahedralMesh();
			});

		connect(regularization_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input iteration times",
				2, 1, 20, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->regularization_hexahedron(false, iter);
			int spend = clock() - t1;
			QString str = "Regularization success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(regularization_faster_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input iteration times",
				2, 1, 20, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->regularization_hexahedron(true, iter);
			int spend = clock() - t1;
			QString str = "Regularization(faster) success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});


		connect(smoothing_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input smoothing times",
				2, 1, 10, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->smoothing_hexahedron(iter);
			int spend = clock() - t1;
			QString str = "Smoothing success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(subdivision_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input subdivision times",
				1, 1, 10, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->subdivision(iter);
			int spend = clock() - t1;
			QString str = "Subdivision success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(delete_island_hexahedral_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->delete_island_HexMesh();
			QString str = "Delete isolated islands success.";
			ui->info->printString(str);
			});

		connect(delete_redundancy_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->delete_redundancy_cell_HexMesh();
			QString str = "Delete redundancy cells success.";
			ui->info->printString(str);
			});

		connect(reconstruction_hexahedral_mesh, &QAction::triggered, [=]() {
			ConvertToIGES();
			});

		connect(add_noise_hexahedral_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			double noise_rate = QInputDialog::getDouble(this, "Hex Mesh", "Please input noise level",
				0.5f, 0.f, 1.f, 2, &input_ok);
			if (input_ok) ui->openGLWidget->add_noise(noise_rate);
			});

		connect(read_force_hexahedral_mesh, &QAction::triggered, [=]() {

			QString filePath = QFileDialog::getOpenFileName(this, tr("Load force_nid file"), tr(""), tr("Files(*.txt)"));
			if (filePath.isEmpty())
				return false;
			else {
				std::string filename = filePath.toStdString();
				ui->openGLWidget->read_constraint_force(filename);
				QString str = "Read force node id success.";
				ui->info->printString(str);
				return true;
			}

			});

		connect(write_force_hexahedral_mesh, &QAction::triggered, [=]() {

			QString filePath = QFileDialog::getSaveFileName(this, tr("Save txt file"), tr(""), tr("txt Files(*.txt)"));
			if (filePath.isEmpty())
				return false;
			else {
				std::string filename = filePath.toStdString();
				ui->openGLWidget->write_constraint_force(filename);
				QString str = "Write force node id success.";
				ui->info->printString(str);
				return true;
			}

			});

		connect(read_spc_hexahedral_mesh, &QAction::triggered, [=]() {

			QString filePath = QFileDialog::getOpenFileName(this, tr("Load spc_nid file"), tr(""), tr("Files(*.txt)"));
			if (filePath.isEmpty())
				return false;
			else {
				std::string filename = filePath.toStdString();
				ui->openGLWidget->read_constraint_spc(filename);
				QString str = "Read spc node id success.";
				ui->info->printString(str);
				return true;
			}

			});

		connect(write_spc_hexahedral_mesh, &QAction::triggered, [=]() {
			QString filePath = QFileDialog::getSaveFileName(this, tr("Save txt file"), tr(""), tr("txt Files(*.txt)"));
			if (filePath.isEmpty())
				return false;
			else {
				std::string filename = filePath.toStdString();
				ui->openGLWidget->write_constraint_spc(filename);
				QString str = "Write spc node id success.";
				ui->info->printString(str);
				return true;
			}

			});

		connect(remove_constraint_hexahedral_mesh, &QAction::triggered, [=]() {
			ui->openGLWidget->remove_all_constraints();
			});

	}
	///************************ Help menu ************************/
	{
		QMenu* menu_help = menu_bar->addMenu("Help");
		menu_help->setTitle("\xE5\xB8\xAE\xE5\x8A\xA9");
		QAction* about = menu_help->addAction("About");
		about->setText("\xE5\x85\xB3\xE4\xBA\x8E");
		QAction* help_doc = menu_help->addAction("Help Doc");
		help_doc->setText("\xE5\xB8\xAE\xE5\x8A\xA9");
		connect(about, &QAction::triggered, [=]() {
			AboutDialog* aboutDialog = new AboutDialog(this);
			aboutDialog->show();
			});

		connect(help_doc, &QAction::triggered, [=]() {
			HelpDoc* helpDoc = new HelpDoc(this);
			helpDoc->show();
			});

	}

	///************************ OpenMesh: poly mesh menu ************************/
	/*{
		QMenu* menu_tri_mesh = menu_bar->addMenu("TriMesh (Beta)");
		QAction* load_tri_mesh = menu_tri_mesh->addAction("Load TriMesh");
		QAction* save_tri_mesh = menu_tri_mesh->addAction("Save TriMesh");
		menu_tri_mesh->addSeparator();
		QAction* smoothing_tri_mesh = menu_tri_mesh->addAction("Smoothing");
		menu_tri_mesh->addSeparator();
		QAction* remeshing_uniform_tri_mesh = menu_tri_mesh->addAction("Remeshing(Uniform)");
		menu_tri_mesh->addSeparator();
		QAction* simplification_qem_tri_mesh = menu_tri_mesh->addAction("Simplification");
		menu_tri_mesh->addSeparator();
		QAction* denoising_bnf_tri_mesh = menu_tri_mesh->addAction("Denoising");
		menu_tri_mesh->addSeparator();

		connect(denoising_bnf_tri_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Triangle Mesh", "Please input denoising times",
				1, 1, 15, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->denoising_bnf_openmesh(iter);
			int spend = clock() - t1;
			QString str = "Denoising success, cost: " + QString::number(spend) + "ms";
			ui->info->printString(str);
			});

		connect(simplification_qem_tri_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			size_t target_vcnt = QInputDialog::getInt(this, "Triangle Mesh", "Please input target vertices size",
				ui->openGLWidget->getVerticesSize(), 3, ui->openGLWidget->getVerticesSize(), 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->simplification_qem_openmesh(target_vcnt);
			time_t t2 = clock();
			QString str = "Simplification success, cost: " + QString::number(int(t2 - t1)) + "ms";
			ui->info->printString(str);
			});

		connect(load_tri_mesh, &QAction::triggered, [=]() {
			QString filename = QFileDialog::getOpenFileName(this, tr("Load Obj file"), tr(""), tr("obj Files(*.obj)"));
			if (!filename.isEmpty()) {
				ui->openGLWidget->loadTriMesh(std::string((const char*)filename.toLocal8Bit()));
				ui->rtnCell->setCheckable(false);
			}
			});

		connect(save_tri_mesh, &QAction::triggered, [=]() {
			QString filename = QFileDialog::getSaveFileName(this, tr("Save Obj file"), tr(""), tr("obj Files(*.obj)"));
			if (!filename.isEmpty())
				ui->openGLWidget->saveTriMesh(std::string((const char*)filename.toLocal8Bit()));
			});
		connect(smoothing_tri_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			int iter = QInputDialog::getInt(this, "Triangle Mesh", "Please input denoising times",
				1, 1, 15, 1, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->smoothing_openmesh(iter);
			time_t t2 = clock();
			QString str = "Smoothing success, cost: " + QString::number(int(t2 - t1)) + "ms";
			ui->info->printString(str);
			});

		connect(remeshing_uniform_tri_mesh, &QAction::triggered, [=]() {
			bool input_ok = false;
			double edge_length_ratio = QInputDialog::getDouble(this, "Triangle Mesh", "Please input length ratio",
				0.5f, 0.f, 1.f, 2, &input_ok);
			if (!input_ok) return;
			time_t t1 = clock();
			ui->openGLWidget->remeshing_uniform_openmesh(edge_length_ratio);
			time_t t2 = clock();
			QString str = "Remeshing(Uniform) success, cost: " + QString::number(int(t2 - t1)) + "ms";
			ui->info->printString(str);
			});


	}*/

	////选择绘制样式（点、边界线、线、填充）
	connect(ui->btnViewPoints, SIGNAL(clicked(bool)), this, SLOT(ViewStylePoints()));
	connect(ui->btnViewLines, SIGNAL(clicked(bool)), this, SLOT(ViewStyleLines()));
	connect(ui->btnViewFaces, SIGNAL(clicked(bool)), this, SLOT(ViewStyleFill()));
	connect(ui->ckbFixView, SIGNAL(clicked(bool)), this, SLOT(UpdateResetViewFlag()));

	connect(ui->btnViewMeshLine, SIGNAL(clicked(bool)), this, SLOT(ViewMeshLine()));
	connect(ui->meshline_width, SIGNAL(valueChanged(int)), this, SLOT(UpdateMeshLineWidth()));

	connect(ui->ckbCutPlane, SIGNAL(clicked(bool)), this, SLOT(ChangeCutFlag()));

	connect(ui->cutPlaneX, SIGNAL(valueChanged(int)), this, SLOT(UpdateCutXYZ()));
	connect(ui->cutPlaneY, SIGNAL(valueChanged(int)), this, SLOT(UpdateCutXYZ()));
	connect(ui->cutPlaneZ, SIGNAL(valueChanged(int)), this, SLOT(UpdateCutXYZ()));
	connect(ui->DensitySlider, SIGNAL(valueChanged(int)), this, SLOT(UpdateTopoptDensity()));
	//	connect(ui->isoValueSlider, SIGNAL(valueChanged(int)), this, SLOT(isoValueSlider()));

	connect(ui->ckbEnablePick, SIGNAL(clicked(bool)), this, SLOT(ChangePickFlag()));
	connect(ui->rtnPoint, SIGNAL(clicked(bool)), this, SLOT(setPickPoint()));
	connect(ui->rtnEdge, SIGNAL(clicked(bool)), this, SLOT(setPickEdge()));
	connect(ui->rtnFace, SIGNAL(clicked(bool)), this, SLOT(setPickFace()));
	connect(ui->rtnCell, SIGNAL(clicked(bool)), this, SLOT(setPickCell()));

	//connect(ui->isRangeSelectRadioButton, SIGNAL(clicked(bool)), this, SLOT());
//    connect(ui->AngleRangeSlider, SIGNAL(valueChanged(int)), this, SLOT(UpdateAngleRangeValue()));
//    connect(ui->AddRangePointPushButton, SIGNAL(clicked(bool)), this, SLOT(AddRangePoints()));
//    connect(ui->AddRangeLinePushButton, SIGNAL(clicked(bool)), this, SLOT(AddRangeLines()));
//    connect(ui->AddRangeFacePushButton, SIGNAL(clicked(bool)), this, SLOT(AddRangeFaces()));
//    connect(ui->AddRangeCellPushButton, SIGNAL(clicked(bool)), this, SLOT(AddRangeCells()));

	connect(ui->EndSelectRadioButton, SIGNAL(clicked(bool)), this, SLOT(DisablePick()));



	connect(ui->btn_pick_remove_one, SIGNAL(clicked(bool)), this, SLOT(removeSelectedOne()));
	connect(ui->btn_pick_remove_all, SIGNAL(clicked(bool)), this, SLOT(removeSelectedAll()));
	connect(ui->btn_delete_pick, SIGNAL(clicked(bool)), this, SLOT(deleteSelectedAll()));
	connect(ui->ckbDrawSelectedAll, SIGNAL(clicked(bool)), this, SLOT(drawSelectAll()));
	//connect(ui->btn_add_cell, SIGNAL(clicked(bool)), this, SLOT(addNewCell()));//添加体


	/*connect(ui->btnUpdateFixedHandle, SIGNAL(clicked(bool)), this, SLOT(updateFixedHandle()));
	connect(ui->ckbRealtimeDeformation, SIGNAL(clicked(bool)), this, SLOT(changeRealtimeDeformation()));*/


	//	connect(ui->MMC_GenSkel, SIGNAL(clicked(bool)), this, SLOT(MMC_GenSkel()));
	//	connect(ui->MMC_GenAxis, SIGNAL(clicked(bool)), this, SLOT(MMC_GenAxis()));
	//	connect(ui->MMC_Opt, SIGNAL(clicked(bool)), this, SLOT(MMC_Opt()));
	//	connect(ui->MMC_appr, SIGNAL(clicked(bool)), this, SLOT(MMC_appr()));
	connect(ui->TopOpt_Run, SIGNAL(clicked(bool)), this, SLOT(TopOpt_Run()));
	connect(ui->TopOpt_Run_2, SIGNAL(clicked(bool)), this, SLOT(TopOpt_Run_2()));
	connect(ui->TopOpt_Run_4, SIGNAL(clicked(bool)), this, SLOT(TopOpt_Run_4()));
	connect(ui->TopOpt_Run_Seq, SIGNAL(clicked(bool)), this, SLOT(TopOpt_Run_Seq()));
	connect(btn_generate_hex_mesh_, SIGNAL(clicked(bool)), this, SLOT(GenerateRegularHexMesh()));
	connect(ui->fix_topopt_result, SIGNAL(clicked(bool)), this, SLOT(SetFinalTopOptMesh()));
	// connect(ui->Protect_topopt, SIGNAL(clicked(bool)), this, SLOT(Protect_topopt()));

	//IGA
	connect(ui->btnConverMeshToXML, SIGNAL(clicked(bool)), this, SLOT(Conver_Mesh_To_XML()));
	connect(ui->btnConverMeshToIGES, SIGNAL(clicked(bool)), this, SLOT(ConvertToIGES_2()));

	connect(ui->line_surface_skel, SIGNAL(clicked(bool)), this, SLOT(Line_surface_skel()));
	connect(ui->save_topopt, SIGNAL(clicked(bool)), this, SLOT(SaveSurfaceMesh()));
	connect(ui->volumetric_sub, SIGNAL(clicked(bool)), this, SLOT(ConvertToIGES()));
	connect(ui->btnTestDetFaces, SIGNAL(clicked(bool)), this, SLOT(testDetFaces()));
	connect(ui->btnTestResample, SIGNAL(clicked(bool)), this, SLOT(testResample()));
	connect(ui->btnReadQuadInfo, SIGNAL(clicked(bool)), this, SLOT(hexmeshingBaseOnMedialMesh()));
	//Build Comobox for IGA Solver
	{
		constexpr auto entries = magic_enum::enum_entries<Solver>();
		for (const auto& e : entries) {
			ui->ComboBoxSelectSolver->addItem(QString(e.second.data()), QVariant(static_cast<int>(e.first)));
		}

		auto debug_name = magic_enum::enum_name(Solver::Debug);
		QString q_name = QString(debug_name.data());
		for (size_t i = 0; i < ui->ComboBoxSelectSolver->count(); ++i)
		{
			if (ui->ComboBoxSelectSolver->itemText(i) == q_name)
			{
				ui->ComboBoxSelectSolver->setCurrentIndex(i);
				break;
			}
		}
	}

	connect(ui->IGA_Vis_X, SIGNAL(clicked(bool)), this, SLOT(IGA_Vis_X()));
	connect(ui->IGA_Vis_Y, SIGNAL(clicked(bool)), this, SLOT(IGA_Vis_Y()));
	connect(ui->IGA_Vis_Z, SIGNAL(clicked(bool)), this, SLOT(IGA_Vis_Z()));
	//connect(ui->Read_Bdy, SIGNAL(clicked(bool)), this, SLOT(Read_Bdy()));
	connect(ui->IGA_Linear_Pro, SIGNAL(clicked(bool)), this, SLOT(IGA_Linear_Pro()));
	//    connect(ui->IGA_hyperelasticity_simulation, SIGNAL(clicked(bool)), this, SLOT(IGA_hyperelasticity_simulation()));
		//connect(ui->pushButton_11, SIGNAL(clicked(bool)), this, SLOT(set_material()));

	connect(ui->IGA_Thermal, SIGNAL(clicked(bool)), this, SLOT(IGA_Thermal()));

	//	connect(ui->Poster_GenSpline, SIGNAL(clicked(bool)), this, SLOT(Poster_GenSpline()));
	connect(ui->Poster_DeleteElement, SIGNAL(clicked(bool)), this, SLOT(Poster_DeleteElement()));
	connect(ui->Poster_Smooth, SIGNAL(clicked(bool)), this, SLOT(Poster_Smooth()));
	//
	//	connect(ui->read_remeshing, SIGNAL(clicked(bool)), this, SLOT(read_remeshing()));
	//	connect(ui->set_remeshing, SIGNAL(clicked(bool)), this, SLOT(set_remeshing()));
	//	connect(ui->set_hole, SIGNAL(clicked(bool)), this, SLOT(set_hole()));
	//	connect(ui->set_hole_filled, SIGNAL(clicked(bool)), this, SLOT(set_hole_filled()));
	//	connect(ui->set_deletehole, SIGNAL(clicked(bool)), this, SLOT(set_deletehole()));
	//	connect(ui->set_holepatch, SIGNAL(clicked(bool)), this, SLOT(set_holepatch()));


		//connect(ui->read_hole, SIGNAL(clicked(bool)), this, SLOT(read_hole()));
		/*connect(ui->btnTest, SIGNAL(clicked(bool)), this, SLOT(testAction()));*/

	ui->openGLWidget->setPickInfo([&]() {
		QString pick_info;
		pick_info = pick_info.fromLocal8Bit(ui->openGLWidget->pick_info.c_str());
		if (!pick_info.isEmpty()) ui->info->printString(pick_info);
		});

	/************************ 设置 UI ************************/
	//ui->btn_delete_pick->setStyleSheet("color: red");
	//ui->btnTest->setVisible(false);
	/*QWidget* tmp = new QWidget();
	tmp->setVisible(true);
	tmp->setWindowTitle("test");*/

}

iGamePoster::~iGamePoster(void) {
}

void iGamePoster::onCheckBoxClicked(QCheckBox* checkBox1, QCheckBox* checkBox2)
{
	QCheckBox* senderCheckBox = qobject_cast<QCheckBox*>(sender());
	QCheckBox* otherCheckBox = (senderCheckBox == checkBox1) ? checkBox2 : checkBox1;

	if (senderCheckBox->isChecked())
	{
		otherCheckBox->setChecked(false);
	}
}


bool iGamePoster::savePickItem()
{
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save Mesh File"), tr(""), tr("Picked Mesh File(*.bdy)"));
	if (!fileName.isEmpty())
	{
		ui->openGLWidget->savePickItem(fileName.toStdString());
		//名字的具体格式规定
		return true;
	}
	return false;
}

bool iGamePoster::openPickItem()
{
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load file"), tr(" "), tr("Picked Mesh File(*.bdy)"));
	if (filePath.isEmpty())
	{
		return false;
	}
	else
	{
		ui->openGLWidget->loadPickFile(filePath.toStdString());
		ui->info->printString(QString("Read bdy file done!"));
	}
	return false;
}

bool iGamePoster::addPickForce()
{
	double force_x = ui->XPosition_force->text().toDouble();
	double force_y = ui->YPosition_force->text().toDouble();
	double force_z = ui->ZPosition_force->text().toDouble();
	std::vector<double> cc{ force_x,force_y,force_z };
	ui->openGLWidget->addPickData(cc, 1);
	return false;
}

bool iGamePoster::addPickFixed()
{
	double fixed_x = ui->XPosition_fixed->text().toDouble();
	double fixed_y = ui->YPosition_fixed->text().toDouble();
	double fixed_z = ui->ZPosition_fixed->text().toDouble();
	std::vector<double> cc{ fixed_x,fixed_y,fixed_z };
	ui->openGLWidget->addPickData(cc, 0);
	return false;
}

bool iGamePoster::addHyperelasticityTimeSettings()
{
	//    ui->openGLWidget->FrameCount_he = ui->SimulationFrameCount_he->text().toInt();
	//    ui->openGLWidget->TimeStep_he = ui->SimulationTimeStep_he->text().toDouble();
	return true;
}

bool iGamePoster::TobeProtected()
{
	ui->openGLWidget->TobeProtected();
	return true;
}

bool iGamePoster::LoadSkeletalMesh(void) {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load ma file"), tr(""), tr("Skeletal Mesh Files(*.ma)"));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		ui->openGLWidget->loadSkeletalMesh(filename);
		int beg = filename.find_last_of('/');
		filePath = QString::fromStdString(filename.substr(++beg, filename.length()));
		ui->info->printString("read " + filePath + " success");
		return true;
	}
}

bool iGamePoster::SaveSkeletalMesh(void) {
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save ma file:"), tr("ma_untitld"), tr("Skeletal Mesh File(*.ma)"));
	if (!fileName.isEmpty()) {
		std::string filename = std::string((const char*)fileName.toLocal8Bit());
		ui->openGLWidget->saveSkeletalMesh(filename);
		int beg = filename.find_last_of('/');
		fileName = fileName.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		QString str = "write " + fileName + " success";
		ui->info->printString(str);
		return true;
	}
	return false;
}

bool iGamePoster::OpenIGESSurface(void) {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load IGES file"), tr(""), tr("IGES Files(*.iges)"));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		ui->openGLWidget->loadIGESSurface(filename);
		int beg = filename.find_last_of('/');
		filePath = QString::fromStdString(filename.substr(++beg, filename.length()));
		ui->info->printString("read " + filePath + " success");
		return true;
	}
}

bool iGamePoster::OpenIGESCurve(void) {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load IGES file"), tr(""), tr("IGES Files(*.iges)"));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		ui->openGLWidget->loadIGESCurve(filename);
		int beg = filename.find_last_of('/');
		filePath = filePath.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		ui->info->printString("read " + filePath + " success");
		return true;
	}
}

bool iGamePoster::SaveIGESFile(void) {
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save IGES file:"), tr("iges_untitld"), tr("IGES File(*.iges)"));
	if (!fileName.isEmpty()) {
		std::string filename = std::string((const char*)fileName.toLocal8Bit());
		ui->openGLWidget->saveIGESFile(filename);
		int beg = filename.find_last_of('/');
		fileName = fileName.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		ui->info->printString("write " + fileName + " success");
		return true;
	}
	return false;
}

bool iGamePoster::LoadHexahedralMesh(void) {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load Mesh file"), tr(""), tr("Mesh Files(*.mesh)"));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		poster_filename = filename;
		ui->openGLWidget->analysis_flag = false;
		int module_cnt = ui->openGLWidget->loadHexMesh(filename);
		int beg = filename.find_last_of('/');
		filePath = filePath.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		auto file_name_without_suffix = filename.substr(0, filename.rfind("."));
		model_name = file_name_without_suffix;

		assert(ui->openGLWidget->isVaildHexMesh());

		{
			ui->openGLWidget->hex_density.clear();
			ui->openGLWidget->hex_density.resize(ui->openGLWidget->hexMesh->CellSize(), 1.0);

			ui->openGLWidget->hexbdy.clear();
			ui->openGLWidget->topopt_pro.clear();
			ui->openGLWidget->topopt_density = 0.5;
			ui->openGLWidget->angle_range_value = 0.;
			ui->openGLWidget->RangeTodoFaces.clear();

			ui->openGLWidget->YoungsModulus = 1e5;
			ui->openGLWidget->PoissonsRatio = 0.3;

			ui->openGLWidget->selected_v_indices.clear();
			ui->openGLWidget->selected_e_indices.clear();
			ui->openGLWidget->selected_f_indices.clear();
			ui->openGLWidget->selected_c_indices.clear();

			ui->openGLWidget->force_v_idx.clear();
			ui->openGLWidget->force_e_idx.clear();
			ui->openGLWidget->force_f_idx.clear();
			ui->openGLWidget->force_c_idx.clear();

			ui->openGLWidget->fix_v_idx.clear();
			ui->openGLWidget->fix_e_idx.clear();
			ui->openGLWidget->fix_f_idx.clear();
			ui->openGLWidget->fix_c_idx.clear();

			ui->openGLWidget->analysis_flag = false;
			ui->openGLWidget->analysis_mode = 0;
		}



		QString str = "read " + filePath + " done!";
		ui->info->printString(str);
		return true;
	}
}

bool iGamePoster::SaveHexahedralMesh(void) {
	if (!ui->openGLWidget->isVaildHexMesh()) {
		ui->info->printString("hexmesh model is null");
		return false;
	}
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save mesh file:"), tr("mesh_untitld"), tr("Mesh File(*.mesh)"));
	if (!fileName.isEmpty()) {
		std::string filename = std::string((const char*)fileName.toLocal8Bit());
		ui->openGLWidget->saveHexMesh(filename);
		int beg = filename.find_last_of('/');
		fileName = fileName.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		QString str = "write " + fileName + " success";
		ui->info->printString(str);
		return true;
	}
	return false;
}

bool iGamePoster::LoadTetrahedronMesh(void) {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load Mesh file"), tr(""), tr("Mesh Files(*.mesh)"));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		ui->openGLWidget->loadTetMesh(filename);
		int beg = filename.find_last_of('/');
		filePath = filePath.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
		QString str = "read " + filePath + " success";
		ui->info->printString(str);
		return true;
	}
}


bool iGamePoster::SaveTetrahedronMesh(void) {
	ui->info->printString("The function is not implemented yet.");
	return false;
}


bool iGamePoster::ConvertToIGES(void) {
	if (!ui->openGLWidget->isVaildHexMesh()) {
		ui->info->printString("hexmesh model is null");
		return false;
	}
	if (ui->openGLWidget->hexMeshConvertToIGES()) ui->info->printString("convert success");
	else ui->info->printString("convert failure");
	return true;
}

bool iGamePoster::ConvertToIGES_2(void) {
	time_t t1 = clock();
	std::string filename = model_name + ".mesh";
	if (filename.empty())
		return false;
	else {
		ui->openGLWidget->loadHexMesh(filename);
	}
	Sleep(1000);
	int spend = clock() - t1;

	if (!ui->openGLWidget->isVaildHexMesh()) {
		ui->info->printString("hexmesh model is null");
		return false;
	}

	int Smooth_Times = ui->Smooth_Times->text().toInt();

	if (ui->openGLWidget->hexMeshConvertToIGES_2(model_name, Smooth_Times)) ui->info->printString("convert success");
	else ui->info->printString("convert failure");
	return true;
}

//bool iGamePoster::Line_surface_skel(void) {
//    time_t t1 = clock();
//    ui->openGLWidget->generate_coverage_axis();
//    int spend = clock() - t1;
//    QString str = "Coverage_axis success, cost: " + QString::number(spend) + "ms";
//    ui->info->printString(str);
//    return true;
//};

bool iGamePoster::Line_surface_skel(void) {
	time_t t1 = clock();
	std::string filename = model_name + ".ma";
	if (filename.empty())
		return false;
	else {
		ui->openGLWidget->loadSkeletalMesh(filename);
	}

	Sleep(1000);
	//ui->openGLWidget->generate_coverage_axis();
	int spend = clock() - t1;
	QString str = "Coverage_axis success, cost: " + QString::number(spend) + "ms";
	//ui->info->printString(str);
	return true;
};

bool iGamePoster::hexmeshingBaseOnMedialMesh(void) {
	time_t t1 = clock();
	//QString filePath = QFileDialog::getOpenFileName(this, tr("Load ma file"), tr(""), tr("Skeletal Mesh Files(*.ma)"));
	//QString filePath =
	std::string filename = model_name + ".mesh";
	if (filename.empty())
		return false;
	else {
		ui->openGLWidget->loadHexMesh(filename);
	}

	Sleep(1000);
	//ui->openGLWidget->generate_coverage_axis();
	int spend = clock() - t1;
	QString str = "Coverage_axis success, cost: " + QString::number(spend) + "ms";
	//ui->info->printString(str);
	return true;
};


bool iGamePoster::LoadSurfaceMesh() {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load file"), tr(""), tr("Obj or Off or Stl or Step or Igs Files(*.obj *.off *.stl *.step *.igs);;OFF Files(*.off);;STL Files(*.stl);;STEP Files(*.step);;IGS Files(*.igs)"

	));
	if (filePath.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)filePath.toLocal8Bit());
		std::string file_name_without_path;
		std::string file_name_without_suffix;

		ui->openGLWidget->loadSurfaceMesh(filename);

		int beg = filename.find_last_of('/');
		file_name_without_path = filename.substr(++beg, filename.length());
		filePath = filePath.fromLocal8Bit(file_name_without_path.c_str());
		QString str = "read " + filePath + " success";

		//QString fileName;
		file_name_without_suffix = filename.substr(0, filename.rfind("."));
		model_name = file_name_without_suffix;
		//fileName = fileName.fromLocal8Bit(file_name_without_suffix.c_str());
		//QString str2 = "fileName: " + fileName;

		ui->info->printString(str);
		//ui->info->printString(str2);
		return true;
	}
}

//bool iGamePoster::hexmeshingBaseOnMedialMesh(void) {
//    ui->openGLWidget->hex_meshing_test();// 基于骨架网格的六面体网格剖分
//    return true;
//}

bool iGamePoster::testResample(void) {

	//    ui->openGLWidget->test_resample();
	//    ui->info->printString("Execute success.");
	return true;
}

bool iGamePoster::testDetFaces(void) {

	ui->openGLWidget->test_det_faces();

	QString str = "Execute success.";
	ui->info->printString(str);
	return true;
}

bool iGamePoster::ConvertToIGES_Skel(void) {// 测试中……
	/*if (!ui->openGLWidget->isVaildHexMesh()) {
		ui->info->printString("hexmesh model is null");
		return false;
	}*/
	ui->openGLWidget->hexMeshConvertToIGES_skeletonExtraction();
	ui->info->printString("convert success");
	return true;
}

bool iGamePoster::InfoClear(void) {
	ui->info->clear();
	return true;
}

bool iGamePoster::InfoCopy(void) {
	ui->info->printString("The function is not implemented yet.");
	return true;
}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

bool iGamePoster::ViewStylePoints() {
	ui->openGLWidget->SetDrawMode(ModelDraw::POINTS);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewStyleLines() {
	ui->openGLWidget->SetDrawMode(ModelDraw::FLATLINES);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewStyleFill() {
	ui->openGLWidget->SetDrawMode(ModelDraw::FILL);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewNormal() {
	ui->openGLWidget->SetColorMode(ModelDraw::NORMAL);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewLight() {
	ui->openGLWidget->SetColorMode(ModelDraw::LIGHT);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewClear() {
	ui->openGLWidget->SetDrawModule(ModelDraw::EMPTY);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ViewMeshLine() {
	ui->openGLWidget->drawMeshLine = !(ui->openGLWidget->drawMeshLine);
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::ChangeCutFlag() {
	if (ui->ckbCutPlane->isChecked()) {
		ui->info->printString("enable cut plane");
		ui->openGLWidget->cutFlag = true;
		UpdateCutXYZ();// 会自动更新 drawdata
	}
	else {
		ui->info->printString("disable cut plane");
		ui->openGLWidget->cutFlag = false;
		ui->openGLWidget->updateDrawData();
	}
	update();
	return true;
}

bool iGamePoster::changeRealtimeDeformation() {
	return false;
	//	if (ui->ckbRealtimeDeformation->isChecked()) {
	//		ui->openGLWidget->enable_realtime_deformation();
	//		//ui->info->printString("enable realtime deformation");
	//		ui->ckbEnablePick->setChecked(false);
	//	} else {
	//		ui->openGLWidget->disable_realtime_deformation();
	//		ui->openGLWidget->updateDrawData();
	//		//ui->info->printString("disable realtime deformation");
	//	}
	//	update();
	//	return true;
}

bool iGamePoster::ChangePickFlag() {
	if (ui->ckbEnablePick->isChecked()) {
		ui->info->printString("change to pick");
		ui->openGLWidget->SetEditModel(ModelDraw::PICKITEM);
		ui->openGLWidget->makeCurrent();
		update();
		/*ui->ckbRealtimeDeformation->setChecked(false);*/
		//ui->info->printString("enable pick");
	}
	else {
		ui->openGLWidget->SetEditModel(ModelDraw::DEFULTEDIT);
		//ui->info->printString("disable pick");
	}
	update();
	return true;
}


bool iGamePoster::ChangeRangePickFlag() {

	return true;
}

bool iGamePoster::AddRangePoints() {
	return true;
}
bool iGamePoster::AddRangeLines() {
	return true;
}
bool iGamePoster::AddRangeFaces() {
	for (auto& fh : ui->openGLWidget->RangeTodoFaces) {
		ui->openGLWidget->selected_f_indices.push_back(fh.idx());
	}
	std::unordered_set<int> seen;
	// 使用 std::unordered_set 去除重复元素
	ui->openGLWidget->selected_f_indices.erase(std::remove_if(ui->openGLWidget->selected_f_indices.begin(), ui->openGLWidget->selected_f_indices.end(), [&](int i) {
		return !seen.insert(i).second;
		}), ui->openGLWidget->selected_f_indices.end());

	ui->openGLWidget->RangeTodoFaces.clear();
	return true;
}
bool iGamePoster::AddRangeCells() {
	return true;
}
bool iGamePoster::GetRangeTempFaces(std::vector<iGameFaceHandle>& faces, double tol) {
	faces.clear();
	std::queue<iGameFaceHandle> q;
	std::unordered_set<iGameFaceHandle> visited;
	for (auto& fid : ui->openGLWidget->selected_f_indices) {
		if (visited.count(iGameFaceHandle(fid)) == 0) {
			q.emplace(fid);
			visited.insert(iGameFaceHandle(fid));
		}
	}
	auto CalcAngle = [&](iGameFaceHandle fh1, iGameFaceHandle fh2)->double {
		auto n1 = ui->openGLWidget->hexMesh->getQuadNormal(fh1);
		auto n2 = ui->openGLWidget->hexMesh->getQuadNormal(fh2);
		double cosAngle = n1.dot(n2) / n1.norm() / n2.norm();
		std::cout << "acos = " << acos(cosAngle) << "\n";

		return acos(cosAngle);
		};
	while (!q.empty()) {
		auto fh = q.front();
		q.pop();
		faces.push_back(fh);

		for (auto& nfh : ui->openGLWidget->hexMesh->NeighborFh(fh)) {
			if (visited.count(nfh)) continue;
			if (!ui->openGLWidget->hexMesh->isOnBoundary(nfh)) continue;
			if (CalcAngle(nfh, fh) > tol) continue;
			visited.insert(nfh);
			q.push(nfh);
		}
	}
	return true;
}

bool iGamePoster::DisablePick() {
	ui->openGLWidget->SetEditModel(ModelDraw::DEFULTEDIT);
	update();
	return true;
}


bool iGamePoster::setPickPoint() {
	ui->openGLWidget->setPickType(0);
	ui->info->printString("pick obj: point");
	ui->openGLWidget->makeCurrent();
	update();
	update();
	return true;
}

bool iGamePoster::setPickEdge() {
	ui->openGLWidget->setPickType(1);
	ui->info->printString("pick obj: edge");
	ui->openGLWidget->makeCurrent();
	update();
	return true;
}

bool iGamePoster::setPickFace() {
	ui->openGLWidget->setPickType(2);
	ui->info->printString("pick obj: face");
	ui->openGLWidget->makeCurrent();
	update();
	update();
	return true;
}

bool iGamePoster::setPickCell() {
	ui->openGLWidget->setPickType(3);
	ui->info->printString("pick obj: cell");
	ui->openGLWidget->makeCurrent();
	update();
	update();
	return true;
}

bool iGamePoster::isoValueSlider() {
	//	std::cout << "ui->isoValueSlider->value() = " << ui->isoValueSlider->value() << std::endl;
	//	int dd = ui->isoValueSlider->value() * 2;
	//	if (dd % 2 == 1) dd + 1;
	//	if (dd <= 6) dd = 6;
	//	if (dd >= 200) dd = 200;
	//	std::cout << "dd = " << dd << std::endl;
	//
	//	if (ui->openGLWidget->remeshing_result.count(dd) == 0) return true;
	//	ui->openGLWidget->surface_mesh = &ui->openGLWidget->remeshing_result[dd];
	//	if (ui->openGLWidget->remeshing_show_model == 2) {
	//		ui->openGLWidget->surface_mesh = &ui->openGLWidget->remeshing_hole_result[dd];
	//	}
	//	else if (ui->openGLWidget->remeshing_show_model == 3) {
	//		ui->openGLWidget->surface_mesh = &ui->openGLWidget->remeshing_holefilled_result[dd];
	//	}
	//	else if (ui->openGLWidget->remeshing_show_model == 4) {
	//		ui->openGLWidget->surface_mesh = &ui->openGLWidget->deletehole[dd];
	//	}
	//	else if (ui->openGLWidget->remeshing_show_model == 5) {
	//		ui->openGLWidget->surface_mesh = &ui->openGLWidget->fillpatch[dd];
	//	}
	//
	//	std::cout << "SetDrawModule\n";
	//	ui->openGLWidget->SetDrawModule(ui->openGLWidget->SURFACEMESH);
	//	std::cout << "updateDrawData\n";
	//	ui->openGLWidget->updateDrawData();
	//	std::cout << "update\n";
	//	update();
	return true;
}

bool iGamePoster::UpdateAngleRangeValue() {
	//    if(!ui->isRangeSelectRadioButton->isChecked()) return true;
	//
	//     auto x = (ui->AngleRangeSlider->value() / 100.f) * M_PI;
	//
	//    GetRangeTempFaces(ui->openGLWidget->RangeTodoFaces, x);
	//    update();
	return true;
}
bool iGamePoster::UpdateTopoptDensity() {

	ui->openGLWidget->topopt_density = (ui->DensitySlider->value() / 100.f);
	ui->openGLWidget->updateDrawData_DensityChange();

	printf("cut rate: %.2f\n", ui->openGLWidget->topopt_density);

	update();
	return true;
}

bool iGamePoster::SetFinalTopOptMesh() {

	double cut_density = ui->openGLWidget->topopt_density;

	std::unordered_set<int> need_to_delete;

	auto& cells = ui->openGLWidget->hexMesh->allcells();

	for (auto& cp : cells) {
		auto vhs = cp.second.getVertexHandle();
		if (ui->openGLWidget->hex_density[cp.first.idx()] <
			ui->openGLWidget->topopt_density) {
			need_to_delete.insert(cp.first.idx());
		}
	}

	for (auto& e : need_to_delete) {
		ui->openGLWidget->hexMesh->DeleteCell(iGameCellHandle(e));
	}

	ui->openGLWidget->saveHexMesh("./temp/result_hex.mesh");
	ui->openGLWidget->loadHexMesh("./temp/result_hex.mesh");

	{
		ui->openGLWidget->hex_density.clear();
		ui->openGLWidget->hex_density.resize(ui->openGLWidget->hexMesh->CellSize(), 1.0);

		ui->openGLWidget->hexbdy.clear();
		ui->openGLWidget->topopt_pro.clear();
		ui->openGLWidget->topopt_density = 0.5;
		ui->openGLWidget->angle_range_value = 0.;
		ui->openGLWidget->RangeTodoFaces.clear();

		ui->openGLWidget->YoungsModulus = 1e5;
		ui->openGLWidget->PoissonsRatio = 0.3;

		ui->openGLWidget->selected_v_indices.clear();
		ui->openGLWidget->selected_e_indices.clear();
		ui->openGLWidget->selected_f_indices.clear();
		ui->openGLWidget->selected_c_indices.clear();

		ui->openGLWidget->force_v_idx.clear();
		ui->openGLWidget->force_e_idx.clear();
		ui->openGLWidget->force_f_idx.clear();
		ui->openGLWidget->force_c_idx.clear();

		ui->openGLWidget->fix_v_idx.clear();
		ui->openGLWidget->fix_e_idx.clear();
		ui->openGLWidget->fix_f_idx.clear();
		ui->openGLWidget->fix_c_idx.clear();

		ui->openGLWidget->analysis_flag = false;
		ui->openGLWidget->analysis_mode = 0;
	}

	return true;
}

bool iGamePoster::UpdateCutXYZ() {
	ui->openGLWidget->cutX = (ui->cutPlaneX->value() / 100.f);
	ui->openGLWidget->cutY = (ui->cutPlaneY->value() / 100.f);
	ui->openGLWidget->cutZ = (ui->cutPlaneZ->value() / 100.f);
	if (ui->ckbCutPlane->isChecked()) ui->openGLWidget->updateDrawData_CutPlane();
	printf("cut rate: %.2f, %.2f, %.2f\n", ui->openGLWidget->cutX, ui->openGLWidget->cutY, ui->openGLWidget->cutZ);
	//ui->info->printString("update cut plane success");
	update();
	return true;
}

bool iGamePoster::UpdateMeshLineWidth() {
	double value = ui->meshline_width->value();
	if (value <= 15) {
		ui->openGLWidget->meshline_width = value / 15 * 1.5f;
	}
	else {
		ui->openGLWidget->meshline_width = 1.5f + (value - 15) * 0.05f;
	}

	ui->openGLWidget->updateDrawData();
	update();
	return true;
}

bool iGamePoster::UpdateResetViewFlag() {
	ui->openGLWidget->resetViewFlag = !ui->ckbFixView->isChecked();
	update();
	return true;
}

bool iGamePoster::removeSelectedOne() {
	ui->openGLWidget->removeSelectedOne();
	return true;
}

bool iGamePoster::removeSelectedAll() {
	ui->openGLWidget->removeSelectedAll();
	return true;
}

bool iGamePoster::deleteSelectedAll() {
	ui->openGLWidget->deleteSelectedAll();
	return true;
}

bool iGamePoster::drawSelectAll() {
	ui->openGLWidget->drawSelectedAllFlag = ui->ckbDrawSelectedAll->isChecked();
	if (ui->openGLWidget->drawSelectedAllFlag) ui->info->printString("enable draw all selected objects");
	else ui->info->printString("disable draw all selected objects");
	return true;
}

bool iGamePoster::addNewCell() {
	ui->openGLWidget->add_new_cell();
	return true;
}

bool iGamePoster::updateFixedHandle() {
	ui->openGLWidget->update_fixed_indices();
	return true;
}

//bool iGamePoster::LoadSurfaceMesh() {
//	QString filePath = QFileDialog::getOpenFileName(this, tr("Load file"), tr(""), tr("Obj or Off or Stl or Step or Igs Files(*.obj *.off *.stl *.step *.igs);;OFF Files(*.off);;STL Files(*.stl);;STEP Files(*.step);;IGS Files(*.igs)"
//
//
//
//
//
//    ));
//	if (filePath.isEmpty())
//		return false;
//	else {
//		std::string filename = std::string((const char*)filePath.toLocal8Bit());
//		ui->openGLWidget->loadSurfaceMesh(filename);
//		int beg = filename.find_last_of('/');
//		filePath = filePath.fromLocal8Bit((filename.substr(++beg, filename.length())).c_str());
//		QString str = "read " + filePath + " success";
//		ui->info->printString(str);
//		return true;
//	}
//}

bool iGamePoster::SaveSurfaceMesh() {
	QString fileName = QFileDialog::getSaveFileName(this, tr("Save model file:"), tr("untitld"), tr("Obj File(*.obj);;IGS File(*.IGS);;Off File(*.off);;Stl File(*.stl);;Step File(*.step);;Vtk File(*.vtk)"));
	if (fileName.isEmpty())
		return false;
	else {
		std::string filename = std::string((const char*)fileName.toLocal8Bit());
		ui->openGLWidget->saveSurfaceMesh(filename);
		int beg = filename.find_last_of('/');
		fileName = fileName.fromLocal8Bit(filename.substr(++beg, filename.length()).c_str());
		QString str = "save " + fileName + " success";
		ui->info->printString(str);
		return true;
	}
}

bool iGamePoster::testAction() {
	/*bool input_ok = true;
	double alpha = QInputDialog::getDouble(this, "Alpha", "",
		0.f, 0, 1.f, 2, &input_ok);
	if (!input_ok) return true;*/
	double alpha = 1;
	ui->openGLWidget->test_action(alpha);
	QString str = "Execute success.";
	ui->info->printString(str);
	return true;
}

void iGamePoster::MMC_GenSkel() {
	QString filename = QFileDialog::getOpenFileName(this, tr("Load Obj file"), tr(""), tr("obj Files(*.obj)"));
	poster_filename = std::string((const char*)filename.toLocal8Bit());
	ui->rtnCell->setCheckable(false);
	bool input_ok = false;
	uint32_t width = QInputDialog::getInt(this, "Poly Mesh", "Please input width",
		722, 0, 2048, 1, &input_ok);
	if (!input_ok) return;
	uint32_t height = QInputDialog::getInt(this, "Poly Mesh", "Please input height",
		122, 0, 2048, 1, &input_ok);
	if (!input_ok) return;
	ui->openGLWidget->BuildSkeletonToFile(poster_filename, width, height);
}

void iGamePoster::MMC_GenAxis() {
	ui->openGLWidget->GenerateSymAxis();
}

void iGamePoster::MMC_Opt() {
	bool input_ok = false;
	double max_error = QInputDialog::getDouble(this, "Poly Mesh", "Please input max error",
		5.f, 1.f, 10.f, 2, &input_ok);
	if (input_ok) ui->openGLWidget->SparseOptimization(max_error);
}

void iGamePoster::MMC_appr() {
	bool input_ok = false;
	uint32_t iter = QInputDialog::getInt(this, "Poly Mesh", "Please input max_iter",
		4000, 0, 1600000, 500, &input_ok);
	if (!input_ok) return;
	ui->openGLWidget->ShapeApproximation(iter);
}
vector<int> vv;
void iGamePoster::TopOpt_Run() {

	QtConcurrent::run([&]() {
		// ui->openGLWidget->updatePickPara();
		ui->openGLWidget->analysis_flag = false;
		//TopOpt(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, ui->openGLWidget->pickcell);
		//TopOpt_CC_Bezier_V2(ui->openGLWidget, poster_filename, ui->openGLWidget->hexbdy, ui->openGLWidget->is_need_update_topopt_model);
		TopOpt_GPU_V0(poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->topopt_pro,
			ui->TopOpt_Rate_2->text().toDouble(),
			ui->TopOpt_Iters_2->text().toDouble(),
			//ui->openGLWidget->YoungsModulus,
			1.0,
			ui->openGLWidget->PoissonsRatio,

			ui->openGLWidget->hex_density,
			ui->openGLWidget->is_need_update_topopt_model);
		});
}

void iGamePoster::TopOpt_Run_2() {

	QtConcurrent::run([&]() {
		// ui->openGLWidget->updatePickPara();
		ui->openGLWidget->analysis_flag = false;
		//TopOpt(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, ui->openGLWidget->pickcell);
		//TopOpt_CC_Bezier_V2(ui->openGLWidget, poster_filename, ui->openGLWidget->hexbdy, ui->openGLWidget->is_need_update_topopt_model);
		TopOpt_GPU_Discrete(poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->topopt_pro,
			ui->TopOpt_Rate_2->text().toDouble(),
			ui->TopOpt_Iters_2->text().toDouble(),
			ui->openGLWidget->YoungsModulus,
			ui->openGLWidget->PoissonsRatio,

			ui->openGLWidget->hex_density,
			ui->openGLWidget->is_need_update_topopt_model);
		});
}

void iGamePoster::TopOpt_Run_4() {

	QtConcurrent::run([&]() {
		// ui->openGLWidget->updatePickPara();
		ui->openGLWidget->analysis_flag = false;
		//TopOpt(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, ui->openGLWidget->pickcell);
		//TopOpt_CC_Bezier_V2(ui->openGLWidget, poster_filename, ui->openGLWidget->hexbdy, ui->openGLWidget->is_need_update_topopt_model);
		TopOpt_GPU_MMC(poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->topopt_pro,
			ui->TopOpt_MMC_type->text().toInt(),
			ui->TopOpt_Rate_2->text().toDouble(),
			ui->TopOpt_Iters_2->text().toDouble(),
			//ui->openGLWidget->YoungsModulus,
			1,
			ui->openGLWidget->PoissonsRatio,
			ui->openGLWidget->hex_density,
			ui->openGLWidget->is_need_update_topopt_model);
		});
}

void iGamePoster::TopOpt_Run_Seq() {

	QtConcurrent::run([&]() {
		ui->openGLWidget->analysis_flag = false;
		TopOpt_GPU_MMC_SIMP_Sequential(poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->topopt_pro,
			ui->TopOpt_MMC_type->text().toInt(),
			ui->TopOpt_Rate_2->text().toDouble(),
			ui->TopOpt_Iters_2->text().toDouble(),
			1.5,
			ui->TopOpt_Iters_2->text().toDouble(),
			1.0,
			2.0,
			ui->openGLWidget->PoissonsRatio,
			ui->openGLWidget->hex_density,
			ui->openGLWidget->is_need_update_topopt_model);
		});
}

void iGamePoster::GenerateRegularHexMesh() {
	QDialog dialog(this);
	dialog.setWindowTitle("Generate regular hex mesh");

	QFormLayout form(&dialog);
	QSpinBox x_cells(&dialog);
	QSpinBox y_cells(&dialog);
	QSpinBox z_cells(&dialog);
	for (QSpinBox* spin_box : { &x_cells, &y_cells, &z_cells }) {
		spin_box->setRange(1, 10000);
		spin_box->setSingleStep(1);
	}
	x_cells.setValue(64);
	y_cells.setValue(32);
	z_cells.setValue(4);
	form.addRow("X cells:", &x_cells);
	form.addRow("Y cells:", &y_cells);
	form.addRow("Z cells:", &z_cells);

	QDialogButtonBox buttons(QDialogButtonBox::Ok | QDialogButtonBox::Cancel, &dialog);
	form.addRow(&buttons);
	connect(&buttons, &QDialogButtonBox::accepted, &dialog, &QDialog::accept);
	connect(&buttons, &QDialogButtonBox::rejected, &dialog, &QDialog::reject);
	if (dialog.exec() != QDialog::Accepted) {
		return;
	}

	const int nx = x_cells.value();
	const int ny = y_cells.value();
	const int nz = z_cells.value();

	const std::filesystem::path output_dir = std::filesystem::path(IGAME_STRINGIFY(INPUT_DIR)) / "topopt";
	std::error_code ec;
	std::filesystem::create_directories(output_dir, ec);
	if (ec) {
		QMessageBox::warning(this, "Mesh generation failed",
			"Could not create the output directory.");
		return;
	}

	const std::string filename = "hex_mesh_" + std::to_string(nx) + "_" +
		std::to_string(ny) + "_" + std::to_string(nz) + ".mesh";
	const std::filesystem::path output_path = output_dir / filename;
	if (!WriteRegularHexMeshFile(output_path, nx, ny, nz)) {
		QMessageBox::warning(this, "Mesh generation failed",
			"Could not write the mesh file.");
		return;
	}

	QMessageBox::information(this, "Mesh generated",
		QString::fromStdString("Generated: " + output_path.string()));
}
void iGamePoster::Protect_topopt() {
	vv = ui->openGLWidget->selected_c_indices;
	ui->openGLWidget->selected_c_indices.clear();
}

void iGamePoster::IGA_Vis_X() {
	ui->openGLWidget->analysis_mode = 0;
	ui->openGLWidget->analysis_flag = true;
	ui->openGLWidget->updateDrawData();
}

void iGamePoster::IGA_Vis_Y() {
	ui->openGLWidget->analysis_mode = 1;
	ui->openGLWidget->analysis_flag = true;
	ui->openGLWidget->updateDrawData();
}
void iGamePoster::IGA_Vis_Z() {
	ui->openGLWidget->analysis_mode = 2;
	ui->openGLWidget->analysis_flag = true;
	ui->openGLWidget->updateDrawData();
}

void iGamePoster::Read_Bdy() {
	QString filePath = QFileDialog::getOpenFileName(this, tr("Load file"), tr(""), tr("Bdy Files(*.bdy)"));
	std::string filename = std::string((const char*)filePath.toLocal8Bit());
	ui->openGLWidget->Read_Bdy(filename);
}

void iGamePoster::read_remeshing() {

	int k;
	MeshKernel::IO io;
	//ui->openGLWidget->loadSurfaceMesh("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\RemesingResult\\TeapotWithPipe_" + std::to_string(200) + ".obj");
	std::cout << "Start\n";
	for (int i = 6; i <= 200; i += 2) {

		/*ui->openGLWidget->remeshing_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\RemshingResult\\TeapotWithPipe_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->remeshing_hole_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\hole\\TeapotWithPipe_with_hole_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->remeshing_holefilled_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\repairhole\\TeapotWithPipe_filled_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->deletehole[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\deletehole\\TeapotWithPipe_deletehole_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->fillpatch[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\HolePatch\\TeapotWithPipe_holepatch_200_vs_" + std::to_string(i) + ".obj", k);*/

		ui->openGLWidget->remeshing_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\RemeshingResult\\TeapotWithPip_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->remeshing_hole_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\hole\\TeapotWithPipe_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->remeshing_holefilled_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\repairhole\\TeapotWithPipe_filled_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->deletehole[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\deletehole\\TeapotWithPipe_deletehole_200_vs_" + std::to_string(i) + ".obj", k);
		ui->openGLWidget->fillpatch[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\Teapot\\HolePatch\\TeapotWithPipe_holepatch_200_vs_" + std::to_string(i) + ".obj", k);


		//ui->openGLWidget->remeshing_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\buddhas\\RemshingResult\\buddhas_" + std::to_string(i) + ".obj", k);
		//ui->openGLWidget->remeshing_hole_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\buddhas\\hole\\buddhas_with_hole_200_vs_" + std::to_string(i) + ".obj", k);
		//ui->openGLWidget->remeshing_holefilled_result[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\buddhas\\repairhole\\buddhas_filled_200_vs_" + std::to_string(i) + ".obj", k);
		//ui->openGLWidget->deletehole[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\buddhas\\deletehole\\buddhas_deletehole_200_vs_" + std::to_string(i) + ".obj", k);
		//ui->openGLWidget->fillpatch[i] = io.ReadObjFile("C:\\Users\\ASUS\\Desktop\\AllPoly\\HoleFillingModel\\output\\buddhas\\HolePatch\\buddhas_holepatch_200_vs_" + std::to_string(i) + ".obj", k);
	}
}
void iGamePoster::set_hole() {
	ui->openGLWidget->remeshing_show_model = 2;
	isoValueSlider();
}
void iGamePoster::set_hole_filled() {
	ui->openGLWidget->remeshing_show_model = 3;
	isoValueSlider();
}

void iGamePoster::set_remeshing() {
	ui->openGLWidget->remeshing_show_model = 1;
	isoValueSlider();
}
void iGamePoster::set_deletehole() {
	ui->openGLWidget->remeshing_show_model = 4;
	isoValueSlider();
}
void iGamePoster::set_holepatch() {
	ui->openGLWidget->remeshing_show_model = 5;
	isoValueSlider();
}
void iGamePoster::set_material() {
	ui->openGLWidget->YoungsModulus = ui->YangShiMoLiang->text().toDouble();
	ui->openGLWidget->PoissonsRatio = ui->PoSongBi->text().toDouble();
}
void iGamePoster::IGA_Linear_Pro() {
	//QtConcurrent::run([&]() {
	ui->info->printString(QString("Start IGA linear elasticity solving..."));
	QCoreApplication::processEvents();
	//ui->openGLWidget->updatePickPara();
	//std::cout << "current: " << poster_filename << std::endl;
	//std::cout << ui->openGLWidget->pickedge.size() << " " << ui->openGLWidget->pickface.size() << std::endl;
	vector<double> ax, ay, az;
	bool bUsingGPU = ui->CheckBoxEnableGPU->isChecked();

	// old
	{
		//Analysis(ui->openGLWidget, poster_filename, ui->openGLWidget->pickface, ui->openGLWidget->pickedge, ax, ay, az);
	}
	// new
	{

		//Analysis_V1(ui->openGLWidget, poster_filename, ui->openGLWidget->hexbdy, ax, ay, az);
	}

	{


		//            Analysis_V2(ui->openGLWidget,
		//                    poster_filename,
		//                    ui->openGLWidget->hexbdy,
		//                    ui->openGLWidget->YoungsModulus,
		//                    ui->openGLWidget->PoissonsRatio,
		//                    ax,
		//                    ay,
		//                    az,
		//                    bUsingGPU);
	}

	{
		Solver SelectedSolver = Solver{ ui->ComboBoxSelectSolver->currentData().toInt() };
		cout << "Selected Solver is " << magic_enum::enum_name(SelectedSolver) << endl;
		Analysis_V3(ui->openGLWidget,
			poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->YoungsModulus,
			ui->openGLWidget->PoissonsRatio,
			SelectedSolver,
			ax,
			ay,
			az);
	}

	ui->openGLWidget->analysis_data[0] = ax, ui->openGLWidget->analysis_data[1] = ay, ui->openGLWidget->analysis_data[2] = az;
	//std::ofstream out1("x.txt");
	//std::ofstream out2("y.txt");
	//std::ofstream out3("z.txt");

	//for (int i = 0; i < ax.size(); ++i) {
	//	for (int j = 0; j < ax[i].size(); ++j) {
	//		out1 << ax[i][j] << " ";
	//		out2 << ay[i][j] << " ";
	//		out3 << az[i][j] << " ";
	//	}

	//}
	//out1.close();
	//out2.close();
	//out3.close();
	//Sleep(2000);
	/*std::cout << ax.size() << std::endl;

	std::ifstream in1("x.txt");
	std::ifstream in2("y.txt");
	std::ifstream in3("z.txt");
	vector<vector<double>> xx, yy, zz;
	for (int i = 0; i < 324; ++i) {
		vector<double> xxx, yyy, zzz;
		for (int j = 0; j < 8; ++j) {
			double xxxx, yyyy, zzzz;
			in1 >> xxxx;
			in2 >> yyyy;
			in3 >> zzzz;
			xxx.push_back(xxxx);
			yyy.push_back(yyyy);
			zzz.push_back(zzzz);
		}
		xx.push_back(xxx);
		yy.push_back(yyy);
		zz.push_back(zzz);
	}
	in1.close();
	in2.close();
	in3.close();*/

	//ui->openGLWidget->analysis_data[0] = xx, ui->openGLWidget->analysis_data[1] = yy, ui->openGLWidget->analysis_data[2] = zz;
	//std::cout << "analysis success!" << std::endl;
	ui->info->printString(QString("Solve done!"));
	ui->openGLWidget->analysis_flag = true;
	ui->openGLWidget->analysis_mode = 1;
	ui->openGLWidget->SetDrawModule(ViewOpenGL::HEXMESH);
	ui->openGLWidget->updateDrawData();
	//});
}

void iGamePoster::IGA_Thermal() {
	std::string filename = poster_filename;

	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(filename);

	std::cout << "read2HexMeshObj" << std::endl;
	file.read2HexMeshObj(mesh, true, true, true);

	double dx = 0.0;
	const auto& d_cell = mesh->cells_[0];
	const auto& d_v_h = d_cell.vertex_[0];
	const auto& d_v_e = d_cell.vertex_[4];

	dx = std::abs(mesh->vertices_[d_v_h].x - mesh->vertices_[d_v_e].x + mesh->vertices_[d_v_h].y - mesh->vertices_[d_v_e].y + mesh->vertices_[d_v_h].z - mesh->vertices_[d_v_e].z);

	std::cout << "dx:" << dx << std::endl;
	IsoThermalSolver* ts = new IsoThermalSolver(ui->Thermal_alpha->text().toDouble(),
		ui->Thermal_stepTime->text().toDouble(),
		ui->Thermal_SumTime->text().toDouble(),
		dx);

	for (const auto& cell : mesh->cells_) {
		CubicBezierSolid cs;

		int Lp_count = 0;
		//		std::vector<VertexHandle> vhs;
		//		for (const auto& v_h : cell.vertex_) {
		//			vhs.push_back(v_h);
		//		}
		//
		//		for (int k = 0; k <= 3; k++) {
		//			V3f v1 = ((mesh->vertices_[vhs[0]] * k) + (mesh->vertices_[vhs[4]] * (3 - k))) / 3;
		//			V3f v2 = ((mesh->vertices_[vhs[1]] * k) + (mesh->vertices_[vhs[5]] * (3 - k))) / 3;
		//			V3f v3 = ((mesh->vertices_[vhs[2]] * k) + (mesh->vertices_[vhs[6]] * (3 - k))) / 3;
		//			V3f v4 = ((mesh->vertices_[vhs[3]] * k) + (mesh->vertices_[vhs[7]] * (3 - k))) / 3;
		//
		//			for (int i = 0; i <= 3; i++) {
		//				V3f v_s = ((v1 * i) + (v4 * (3 - i))) / 3;
		//				V3f v_t = ((v2 * i) + (v3 * (3 - i))) / 3;
		//
		//				for (int j = 0; j <= 3; j++) {
		//					cs.control_points.push_back(Eigen::Vector4d((v_s.x * j + v_t.x * (3 - j)) / 3, (v_s.y * j + v_t.y * (3 - j)) / 3, (v_s.z * j + v_t.z * (3 - j)) / 3, Lp_count++));
		//				}
		//			}
		//
		//			//std::cout << v1.x << v1.y << v1.z << std::endl;
		//		}
		for (const auto& v_h : cell.vertex_) {
			cs.control_points.push_back(Eigen::Vector4d(mesh->vertices_[v_h].x, mesh->vertices_[v_h].y, mesh->vertices_[v_h].z, Lp_count++));
		}

		ts->addBezierSolid(cs);
	}

	std::vector<std::pair<Eigen::Vector4d, double>> TemperatureField;
	ts->run(TemperatureField);

	if (ui && ui->openGLWidget) {
		std::vector<double> sliced_temp(TemperatureField.size(), 0.0);

		for (size_t i = 0; i < TemperatureField.size(); ++i) {
			const auto& coord = TemperatureField[i].first;
			const double temp = TemperatureField[i].second;

			double x2 = coord.x();
			double y2 = coord.y();
			double z2 = coord.z();
			for (const auto& cell : mesh->cells_) {
				bool fl = false;
				for (const auto& v_h : cell.vertex_) {
					double x1 = mesh->vertices_[v_h].x;
					double y1 = mesh->vertices_[v_h].y;
					double z1 = mesh->vertices_[v_h].z;
					double dist = std::abs(x1 - x2) + std::abs(y1 - y2) + std::abs(z1 - z2);
					if (dist < 1e-4) {
						sliced_temp[v_h] = temp;
						fl = true;
						break;
					}
				}
				if (fl) {
					break;
				}
			}
		}

		ui->openGLWidget->analysis_data[0] = sliced_temp;
		ui->openGLWidget->analysis_data[1] = sliced_temp;
		ui->openGLWidget->analysis_data[2] = sliced_temp;
		ui->openGLWidget->analysis_flag = true;
		ui->openGLWidget->analysis_mode = 0;
		ui->openGLWidget->SetDrawModule(ViewOpenGL::HEXMESH);
		ui->openGLWidget->updateDrawData();
		QApplication::processEvents();
	}
}

void iGamePoster::IGA_hyperelasticity_simulation() {
	ui->info->printString(QString("Start IGA hyperelasticity solving..."));

	ui->info->printString(QString("Now Only Support PCG GPU Solver"));

	QCoreApplication::processEvents();

	Solver SelectedSolver = Solver{ ui->ComboBoxSelectSolver->currentData().toInt() };
	cout << "Selected Solver is " << magic_enum::enum_name(SelectedSolver) << endl;

	//for debug
/*
	 Hyperelasticity_Analysis_V1(
				ui->openGLWidget,
				poster_filename,
				ui->openGLWidget->hexbdy,
				ui->openGLWidget->YoungsModulus,
				ui->openGLWidget->PoissonsRatio,
				ui->openGLWidget->FrameCount_he,
				ui->openGLWidget->TimeStep_he,
				SelectedSolver,
				ui->openGLWidget->analysis_data,
				ui->openGLWidget->is_need_update_spline_simulation_model
				);
*/


	QFuture<void> result = QtConcurrent::run([&]() {
		Hyperelasticity_Analysis_V1(
			ui->openGLWidget,
			poster_filename,
			ui->openGLWidget->hexbdy,
			ui->openGLWidget->YoungsModulus,
			ui->openGLWidget->PoissonsRatio,
			ui->openGLWidget->FrameCount_he,
			ui->openGLWidget->TimeStep_he,
			SelectedSolver,
			ui->openGLWidget->analysis_data,
			ui->openGLWidget->is_need_update_spline_simulation_model
		);
		});

	//result.waitForFinished();
	//ui->info->printString(QString("Solve done!"));
}

void iGamePoster::Conver_Mesh_To_XML()
{
	std::string filename = poster_filename;

	std::filesystem::path p(filename);
	p.replace_extension("xml");
	std::string output_filename = p.string();

	HexV3fMesh* mesh = new HexV3fMesh();
	BSplineSolid m_solids;
	FileHandler file(filename);
	//ui->openGLWidget->loadHexMesh(filename);
	size_t patch_size = ui->openGLWidget->hexMesh->CellSize();

	std::cout << "read2HexMeshObj" << std::endl;
	file.read2HexMeshObj(mesh, true, true, true);
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	m_solids = AppCCSIMP.bss;

	IGA_CC_Helper::write2Bezier(output_filename, m_solids);
	std::cout << "XML Convert end" << std::endl;
}



void iGamePoster::Poster_Smooth() {
	bool input_ok = false;
	int iter = QInputDialog::getInt(this, "Hex Mesh", "Please input smoothing times",
		2, 1, 10, 1, &input_ok);
	if (!input_ok) return;
	time_t t1 = clock();
	ui->openGLWidget->smoothing_hexahedron(iter);
	int spend = clock() - t1;
	QString str = "Smoothing success, cost: " + QString::number(spend) + "ms";
	ui->info->printString(str);
}

void iGamePoster::Poster_DeleteElement() {
	ui->openGLWidget->delete_island_HexMesh();
	QString str = "Delete isolated islands success.";
	ui->info->printString(str);
}

void iGamePoster::Poster_GenSpline() {
	ConvertToIGES();
}

bool iGamePoster::auto_ShowBsplineSolid() {
	ui->openGLWidget->ClearAllVector();
	int udegree, vdegree, wdegree, ucnum, vcnum, wcnum, uknum, vknum, wknum;//degree是次数
	vector<vector<vector<CPoint>>>controlpoints;//保存3个参数方向的控制点
	vector<vector<CPoint>> singlepatch;//保存u/v两个参数方向的片
	vector<CPoint>singlepoints;//临时保存u/v/w单个参数方向的控制点
	vector<double>uknots, vknots, wknots;
	vector<int>umultis, vmultis, wmultis;
	CPoint du(0.0, 0.0, 0.0);//初始控制点坐标
	ucnum = ui->bssUcnum->text().toInt();//前端传入控制点个数
	vcnum = ui->bssVcnum->text().toInt();
	wcnum = ui->bssWcnum->text().toInt();
	udegree = ui->bssUdnum->text().toInt();//前端传入次数
	vdegree = ui->bssVdnum->text().toInt();
	wdegree = ui->bssWdnum->text().toInt();
	if (udegree >= ucnum || vdegree >= vcnum || wdegree >= wcnum)
	{
		QMessageBox::warning(this, QString::fromLocal8Bit("warning"), QString::fromLocal8Bit("The number of control points should be greater than the order"));
		return false;
	}
	else
	{
		ui->openGLWidget->aotu_createBSS(ucnum, vcnum, wcnum, udegree, vdegree, wdegree, controlpoints, singlepatch, singlepoints, uknots, vknots, wknots, umultis, vmultis, wmultis);

		std::vector<double> _uknots, _vknots, _wknots;
		for (int i = 0; i < uknots.size(); ++i)
		{
			for (int j = 0; j < umultis[i]; ++j)
				_uknots.push_back(uknots[i]);
		}
		for (int i = 0; i < vknots.size(); ++i)
		{
			for (int j = 0; j < vmultis[i]; ++j)
				_vknots.push_back(vknots[i]);
		}
		for (int i = 0; i < wknots.size(); ++i)
		{
			for (int j = 0; j < wmultis[i]; ++j)
				_wknots.push_back(wknots[i]);
		}


		ui->openGLWidget->showHexMesh(true);
		ui->openGLWidget->getBsplineSolidFromParam(udegree, vdegree, wdegree, controlpoints, _uknots, _vknots, _wknots, umultis, vmultis, wmultis);
		ui->openGLWidget->bsplineSolidToDiscreteXJ();
		ui->openGLWidget->SetDrawModule(ModelDraw::BSPLINESOLID);
		ui->openGLWidget->makeCurrent();
		update();
		//ui->openGLWidget->setInfoFromSurface(ui->openGLWidget->getBsplineSurfaces());
		return true;

	}

}
bool iGamePoster::manually_ShowBsplineSolid() {
	//ui->openGLWidget->ClearAllVector();
	//ifstream in("./C.txt");
	////ifstream in("E:\\poster\\poster\\test1.txt");
	//string line;
	//int count = 1;
	//while (getline(in, line)) {
	//	int cnt = 0;
	//	if (line[0] == 'p') {
	//		int udegree, vdegree, wdegree, ucnum, vcnum, wcnum, uknum, vknum, wknum;   //degree是次数
	//		vector<vector<double>> knots(3);
	//		vector<double>uknots, vknots, wknots;
	//		vector<int> degrees(3);
	//		stringstream ss;
	//		for (int i = 0; i < 3; ++i) {
	//			line.clear();
	//			getline(in, line);
	//			
	//			ss.str("");
	//			ss.clear();
	//			ss.str(line);
	//			double kk;
	//			while (ss >> kk) {
	//				knots[i].push_back(kk);
	//				//std::cout << kk << std::endl;
	//			}
	//		}
	//		for (int i = 0; i < 3; ++i) {
	//			line.clear();
	//			getline(in, line);
	//			ss.str("");
	//			ss.clear();
	//			ss.str(line);
	//			ss >> degrees[i];
	//			//std::cout << degrees[i]<< std::endl;
	//		}
	//		ucnum = knots[0].size() - degrees[0] - 1;
	//		vcnum = knots[1].size() - degrees[1] - 1;
	//		wcnum = knots[2].size() - degrees[2] - 1;
	//		line.clear();
	//		getline(in, line);
	//		ss.str("");
	//		ss.clear();
	//		ss.str(line);
	//
	//
	//		vector<vector<vector<CPoint>>>controlpoints(ucnum, vector<vector<CPoint>>(vcnum, vector<CPoint>(wcnum)));//保存3个参数方向的控制点
	//		for (int i = 0; i < ucnum; ++i) {
	//			for (int j = 0; j < vcnum; ++j) {
	//				for (int s = 0; s < wcnum; ++s) {
	//					double x, y, z;
	//					ss >> x >> y >> z;
	//					controlpoints[i][j][s] = { x*0.5, y*0.5, z*0.5 };
	//				}
	//			}
	//		}
	//		udegree = degrees[0];  vdegree = degrees[1];  wdegree = degrees[2];
	//		uknum = knots[0].size();   vknum = knots[1].size();   wknum = knots[2].size();
	//		uknots = knots[0];  vknots = knots[1]; wknots = knots[2];

	//		vector<int>umultis, vmultis, wmultis;
	//		
	//		int handle=ui->openGLWidget->getBsplineSolidFromParamXJ(udegree, vdegree, wdegree, controlpoints, uknots, vknots, wknots, umultis, vmultis, wmultis);
	//		ui->openGLWidget->SaveIndex(controlpoints, handle);
	//	}
	//	//break;
	//}		

	QString filePath = QFileDialog::getOpenFileName(this, tr("Load file"), tr(" "), tr("Picked Mesh File(*.txt)"));
	if (filePath.isEmpty())
	{
		return false;
	}
	else
	{
		//ui->openGLWidget->loadPickFile(filePath.toStdString());
		//ui->info->printString(QString("Read bdy file done!"));
		ui->openGLWidget->ClearAllVector();
		std::ifstream file(filePath.toStdString(), std::ios::in);
		std::string line;

		std::getline(file, line);
		std::stringstream ss(line);
		int num;
		ss >> num;
		//num = 3;//test
		std::cout << "num=" << num << std::endl;
		for (int linenum = 0; linenum < num; ++linenum) {
			int cnt = 0;
			int udegree, vdegree, wdegree, ucnum, vcnum, wcnum, uknum, vknum, wknum;   //degree是次数
			vector<vector<double>> knots(3);
			vector<double>uknots, vknots, wknots;
			vector<int> degrees(3);
			stringstream ss;
			for (int i = 0; i < 3; i++)
			{
				for (int j = 0; j < 4; ++j)
				{
					knots[i].push_back(0);
				}
				for (int j = 0; j < 4; ++j)
				{
					knots[i].push_back(1);
				}
			}
			degrees[0] = degrees[1] = degrees[2] = 3;

			ucnum = knots[0].size() - degrees[0] - 1;
			vcnum = knots[1].size() - degrees[1] - 1;
			wcnum = knots[2].size() - degrees[2] - 1;

			//std::cout<<"ucnum= "<<ucnum<<"  , vcnum= "<<vcnum<<" , wcnum= "<<wcnum<<std::endl;

			vector<vector<vector<CPoint>>>controlpoints(ucnum, vector<vector<CPoint>>(vcnum, vector<CPoint>(wcnum)));//保存3个参数方向的控制点
			for (int i = 0; i < ucnum; ++i) {
				for (int j = 0; j < vcnum; ++j) {
					for (int s = 0; s < wcnum; ++s) {
						std::getline(file, line);
						ss.clear();
						ss.str(line);
						double x, y, z;
						ss >> x >> y >> z;
						//std::cout<<"x= "<<x<<" , y = "<<y<<" , z = "<<z<<std::endl;
						line.clear();
						controlpoints[i][j][s] = { x, y, z };
						//std::cout<<	controlpoints[i][j][s].getX()<<",  "<<	controlpoints[i][j][s].getY()<<", "<<	controlpoints[i][j][s].getZ()<<std::endl;
					}
				}
			}
			udegree = degrees[0];  vdegree = degrees[1];  wdegree = degrees[2];
			uknum = knots[0].size();   vknum = knots[1].size();   wknum = knots[2].size();
			uknots = knots[0];  vknots = knots[1]; wknots = knots[2];

			vector<int>umultis, vmultis, wmultis;



			int handle = ui->openGLWidget->getBsplineSolidFromParamXJ(udegree, vdegree, wdegree, controlpoints, uknots, vknots, wknots, umultis, vmultis, wmultis);
			ui->openGLWidget->SaveIndex(controlpoints, handle);
		}

	}



	ui->openGLWidget->showHexMesh(true);
	ui->openGLWidget->bsplineSolidToDiscreteXJ();

	ui->openGLWidget->SetDrawModule(ModelDraw::BSPLINESOLID);
	ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
	ui->openGLWidget->makeCurrent();

	update();
	//ui->openGLWidget->setInfoFromSurface(ui->openGLWidget->getBsplineSurfaces());
	return true;


}

bool iGamePoster::simulation()
{
	QtConcurrent::run([&]() {
		ui->openGLWidget->bsplineSolidSimulation();
		});
	return true;
}

bool iGamePoster::EditBsplineSolid() {
	if (ui->editbsplinesolid->isChecked()) {
		ui->openGLWidget->getControlMeshStatus(true);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
	else {
		ui->openGLWidget->getControlMeshStatus(false);
		ui->openGLWidget->SetEditModel(ModelDraw::DEFULTEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
}

bool iGamePoster::IF_ShowBsplineSolid()//显示b样条体
{
	if (ui->editbsplinesolid_2->isChecked()) {
		ui->openGLWidget->showHexMesh(true);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
	else {
		ui->openGLWidget->showHexMesh(false);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
}
bool iGamePoster::IF_Real_time_rander()//实时渲染样条体
{
	if (ui->editbsplinesolid_3->isChecked()) {
		ui->openGLWidget->set_isRealtime(true);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
	else {
		ui->openGLWidget->set_isRealtime(false);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
}
bool iGamePoster::IF_ShowBsplineSolid_inside()//渲染样条体内部
{
	if (ui->editbsplinesolid_4->isChecked()) {
		ui->openGLWidget->set_isShowinside(true);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
	else {
		ui->openGLWidget->set_isShowinside(false);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
}
bool iGamePoster::Choose_Face()//选面模式，修改mousepress鼠标事件
{
	if (ui->editbsplinesolid_5->isChecked()) {
		ui->openGLWidget->set_chooseface(true);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
	else {
		ui->openGLWidget->set_chooseface(false);
		ui->openGLWidget->SetEditModel(ModelDraw::BSPLINESOLIDEDIT);
		ui->openGLWidget->makeCurrent();
		update();
		return true;
	}
}

void iGamePoster::generate_F()
{
	//	double force_x = ui->Force_X->text().toDouble();
	//	double force_y = ui->Force_Y->text().toDouble();
	//	double force_z = ui->Force_Z->text().toDouble();
	//	ui->openGLWidget->set_F(force_x, force_y, force_z);
}

void iGamePoster::generate_V()
{
	//	double offset_x = ui->Offset_X->text().toDouble();
	//	double offset_y = ui->Offset_Y->text().toDouble();
	//	double offset_z = ui->Offset_Z->text().toDouble();
	//	bool flag1 = ui->x_radioButton->isChecked();
	//	bool flag2 = ui->y_radioButton->isChecked();
	//	bool flag3 = ui->z_radioButton->isChecked();
	//	int dir = flag1 * 4 + flag2 * 2 + flag3;
	//	ui->openGLWidget->set_V(offset_x, offset_y, offset_z,dir);
}

