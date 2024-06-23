#include "ModelDraw.h"
#include "UI/SingalFather.h"
#include "Geom/PIA.h"
#include "Geom/FileHandler.h"
#include "NewFilter/ErrorRender.h"
using namespace std;

ModelDraw::ModelDraw(QWidget* parent) :ViewOpenGL(parent), isEnableLighting(false), isTwoSideLighting(false), editmodel(DEFULTEDIT) {
	genColors();
	/*Temp-lee*/
	oldHexMesh = new HexV3fMesh();
	depth = 0;
	bshowAxis = false;
	bshowGround = false;
	bshowSKEL = false;
	bdragMode = false;
	copyTime = false;
	//bRotateByMouse = false;
	bVertexSnap = false;
	bshowTriMesh = false;
	solutionSkel = new SolutionSkel();
	projectSkel = new Project();
	mouse_move_flag = false;
	local_coordinater_axis = new CoordinateAxis();

	proDraw = new ProDraw(this);
	dotLineCalculation = new DotLineCalculation(this);
}
ModelDraw::~ModelDraw(void) {

}

bool ModelDraw::loadSurfaceMesh(std::string filename) {
	ModelClear();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	//std::cout << "file format: " << fileformat << std::endl;

	MeshKernel::IO io;
	if (fileformat == "obj")
	{
		//this->triMesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename));
		this->quadMesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename));
	}
	else if (fileformat == "off") {
		//this->triMesh = new MeshKernel::SurfaceMesh(io.ReadOffFile(filename));
		this->quadMesh = new MeshKernel::SurfaceMesh(io.ReadOffFile(filename));

	}
	else if (fileformat == "stl") {
		std::cout << "support obj and off files only!\n";
		return false;
	}

	///* iGameKernel 三角网格读入以及 KD 树使用示例 */
	//MeshKernel::SurfaceMesh* mesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename, tmp));
	///*if (fileformat == "obj") {
	//	this->triMesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename, tmp));
	//} 	else if (fileformat == "off") {
	//	this->triMesh = new MeshKernel::SurfaceMesh(io.ReadOffFile(filename, tmp));
	//}*/
	//TriangleKdTree* kd_tree = new TriangleKdTree(*triMesh, 1);
	//for (auto& vp : mesh->allvertices()) {
	//	auto pos = vp.second;
	//	std::cout << "Input pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << "; ";
	//	auto face_nearest = kd_tree->nearest(pos);// 输入一个顶点(位置), 函数返回该曲面模型上离这个点最近的点坐标(可能在面上或边上)
	//	pos = face_nearest.nearest;
	//	std::cout << "Nearest pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
	//}
	
	SetDrawModule(SURFACEMESH);
	SetDrawModule(QUADMESH);
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	//updateDrawData_SurfaceMesh();
	update();
	return true;
}

bool ModelDraw::saveSurfaceMesh(std::string filename) {
	if (!this->triMesh) return false;
	triMesh->updateAllHandles();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	MeshKernel::IO io;
	if (fileformat == "obj") {
		io.WriteObjFile(*(this->triMesh), filename);
	}
	else if (fileformat == "off") {
		io.WriteOffFile(*(this->triMesh), filename);
	}
	else if (fileformat == "stl") {
		//io.WriteStlFile(*(this->Mesh), fileName);
		std::cout << "support obj and off files only!\n";
		return false;
	}
	else if (fileformat == "vtk") {
		//io.WriteVtkFile(*(this->Mesh), fileName);
		std::cout << "support obj and off files only!\n";
		return false;
	}
	return true;
}

int ModelDraw::loadMesh(std::string filename) 
{

	MeshFileIO io;
	int d, c;
	io.read_mesh_file(filename, d, c);
	std::vector<Vex> _vertices;
	std::vector<std::vector<VH>> _vhs;
	io.get_mesh_data(_vertices, _vhs);

	ModelClear();
	if (d == 2) 
	{
		if (c == 4) 
		{// 四边形网格
			quadMesh = new MeshKernel::SurfaceMesh(_vertices, _vhs);
			std::cout << "Initialize success. #V = " << quadMesh->vsize() << ", #F = " << quadMesh->fsize() << std::endl;
			SetDrawModule(QUADMESH);
		} else if (c == 3) 
		{// 三角形网格
			triMesh = new MeshKernel::SurfaceMesh(_vertices, _vhs);
			std::cout << "Initialize success. #V = " << triMesh->vsize() << ", #F = " << triMesh->fsize() << std::endl;
			SetDrawModule(SURFACEMESH);
		}
	} else if (d == 3) 
	{
		if (c == 4) 
		{// 四面体网格
			tetMesh = new MeshKernel::VolumeMesh(_vertices, _vhs);
			std::cout << "Initialize success. #V = " << tetMesh->vsize() << ", #C = " << tetMesh->csize() << std::endl;
			SetDrawModule(TETMESH);
		} else if (c == 6)
		{// 六面体网格
			hexMesh = new MeshKernel::VolumeMesh(_vertices, _vhs);
			//hexMesh = new MeshKernel::VolumeMesh();
			std::cout << "Initialize success. #V = " << hexMesh->vsize() << ", #C = " << hexMesh->csize() << std::endl;
			int facenum = 0;
			int allfacenum = 0;
			for (auto f : hexMesh->allfaces())
			{
				allfacenum++;
				if (hexMesh->isOnBoundary(f.first))
					facenum++;
			}
			std::cout << "All facenum: " << allfacenum << std::endl;
			std::cout << "Boundary facenum: " << facenum << std::endl;
			hex_mesh_tool_ = new HexMeshTool(hexMesh);
			//hex_mesh_tool_->CreateCilcle(Vec(0, 0, 0), 2.0, 1.0, 1.0, 8);
			//hex_mesh_tool_->CreateCilcle(Vec(4, 10, 0), 2.0, 1.0, 1.0, 8);
			SetDrawModule(HEXMESH);
		}
	}
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	update();

	return 0;
}

bool ModelDraw::saveMesh(std::string filename) 
{
	ControlMeshChangeToQuadMesh();
	MeshKernel::IO io1;
	int beg = filename.find_last_of('.');
	string fileformat = filename.substr(++beg, filename.length());
	if (fileformat == "off") {
		if (this->triMesh) io1.WriteOffFile(*triMesh, filename);
		return true;
	}

	MeshFileIO io;
	if (this->hexMesh) {
		return io.write_mesh_file(filename, *hexMesh);
	} else if (this->quadMesh) {
		return io.write_mesh_file(filename, *quadMesh);
	} else if (this->tetMesh) {
		return io.write_mesh_file(filename, *tetMesh);
	} else if (this->triMesh) {
		return io.write_mesh_file(filename, *triMesh);
	}
	return false;
}

void ModelDraw::drawHexMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	glDepthRange(0.01, 1.0);
	shaderProgram.setUniformValue("colormode", -1);
	glPointSize(9);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);


	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		// 绘制线框
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width);
		shaderProgram.setUniformValue("colormode", 530);
		core->glDrawArrays(GL_LINES, tri_cnt * 3, line_cnt * 2);

		{
			// 绘制约束

			glDepthRange(0.01, 0.99999);// 放在前面
			glPointSize(9);

			shaderProgram.setUniformValue("colormode", 231);
			glBegin(GL_POINTS);
			for (auto& fp : spc_id) {
				VH vh(fp.first);
				if (!hexMesh->isValid(vh)) continue;
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
			glEnd();

			shaderProgram.setUniformValue("colormode", 232);
			glBegin(GL_POINTS);
			for (auto& fp : force_id) {
				VH vh(fp.first);
				if (!hexMesh->isValid(vh)) continue;
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
			glEnd();

		}

	}

	if (translucence_tricnt) {// 请放在最后绘制
		glDepthRange(0.01, 2.0);// 半透明：可以被其它物体遮挡
		shaderProgram.setUniformValue("colormode", 518);// 半透明：alpha < 1
		core->glDrawArrays(GL_TRIANGLES, tri_cnt * 3 + line_cnt * 2, translucence_tricnt * 3);
	}

	update();
}

void ModelDraw::drawQuadMesh() 
{
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	glDepthRange(0.01, 1.0);
	shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) 
	{

		// 绘制线框
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		if (editmodel == EditModel::PICKITEM && pickType == PickType::PickEdge) 
		{
			glLineWidth(meshline_width * 2.5f);
		} else 
		{
			glLineWidth(meshline_width);
		}

		shaderProgram.setUniformValue("colormode", 530);
		core->glDrawArrays(GL_LINES, tri_cnt * 3, line_cnt * 2);

	}

	{
		// 绘制约束
		shaderProgram.setUniformValue("colormode", 231);
		glDepthRange(0.01, 0.99999);
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& fp : spc_id) {
			VH vh(fp.first);
			if (!quadMesh->isValid(vh)) continue;
			auto& v = quadMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

		shaderProgram.setUniformValue("colormode", 232);
		glBegin(GL_POINTS);
		for (auto& fp : force_id) {
			VH vh(fp.first);
			if (!quadMesh->isValid(vh)) continue;
			auto& v = quadMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

	}

	update();
}

void ModelDraw::drawQuadMeshfromHexmesh()
{
	CAD_Reconstruction app;
	this->quadMesh = app.TestDraw_generate_surface_mesh(*this->hexMesh);
	SetDrawModule(QUADMESH);
	updateDrawData();
}

void ModelDraw::drawTetMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	glDepthRange(0.01, 1.0);
	shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		// 绘制线框
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		if (editmodel == EditModel::PICKITEM && pickType == PickType::PickEdge) {
			glLineWidth(meshline_width * 2.5f);
		} else {
			glLineWidth(meshline_width);
		}

		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < tri_cnt; ++i) {
			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
		}

	}

	{
		// 绘制约束
		shaderProgram.setUniformValue("colormode", 231);
		glDepthRange(0.01, 0.99999);
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& fp : spc_id) {
			VH vh(fp.first);
			if (!tetMesh->isValid(vh)) continue;
			auto& v = tetMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

		shaderProgram.setUniformValue("colormode", 232);
		glBegin(GL_POINTS);
		for (auto& fp : force_id) {
			VH vh(fp.first);
			if (!tetMesh->isValid(vh)) continue;
			auto& v = tetMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

	}

	if (translucence_tricnt) {// 请放在最后绘制
		glDepthRange(0.01, 2.0);// 半透明：可以被其它物体遮挡
		shaderProgram.setUniformValue("colormode", 518);// 半透明：alpha < 1
		core->glDrawArrays(GL_TRIANGLES, tri_cnt * 3 + line_cnt * 2, translucence_tricnt * 3);
	}


	update();
}

void ModelDraw::TestDrawSurfaceMesh()
{
	core->glBindVertexArray(TestVAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	shaderProgram.setUniformValue("colormode", 1);
	if (show_mesh_quality) shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE)
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width);
		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < tri_cnt; ++i)
		{
			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
		}

	}

	{
		// 绘制约束
		shaderProgram.setUniformValue("colormode", 231);
		glDepthRange(0.01, 0.99999);
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& fp : spc_id)
		{
			VH vh(fp.first);
			if (!triMesh->isValid(vh)) continue;
			auto& v = triMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

		shaderProgram.setUniformValue("colormode", 232);
		glBegin(GL_POINTS);
		for (auto& fp : force_id)
		{
			VH vh(fp.first);
			if (!triMesh->isValid(vh)) continue;
			auto& v = triMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

	}

	update();
}

void ModelDraw::drawSurfaceMesh()
{
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	shaderProgram.setUniformValue("colormode", 1);
	if (show_mesh_quality) shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) 
	{
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width);
		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < tri_cnt; ++i)
		{
			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
		}
		
	}

	{
		// 绘制约束
		shaderProgram.setUniformValue("colormode", 231);
		glDepthRange(0.01, 0.99999);
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& fp : spc_id) 
		{
			VH vh(fp.first);
			if (!triMesh->isValid(vh)) continue;
			auto& v = triMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

		shaderProgram.setUniformValue("colormode", 232);
		glBegin(GL_POINTS);
		for (auto& fp : force_id) 
		{
			VH vh(fp.first);
			if (!triMesh->isValid(vh)) continue;
			auto& v = triMesh->vertices(vh);
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

	}

	update();

}

void ModelDraw::updateDrawData_HexMesh_CutPlane() {
	if (!this->hexMesh) return;

	std::vector<GLfloat> translucence_vertices_data;
	std::vector<GLfloat> line_vertices_data;

	/*double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];*/
	//("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	auto& cells = hexMesh->allcells();
	std::unordered_map<int, bool> current_drawed;
	for (auto& cp : cells) {
		bool drawFlag = true;
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& vex = hexMesh->vertices(vh);
		/*	if (vex.x() < draw_min_x || vex.y() < draw_min_y || vex.z() < draw_min_z) {
				drawFlag = false;
				break;
			}*/
		}
		Vector3d color = (hexMesh->isOnBoundary(cp.first)) ? Vector3d(0.349, 0.6, 0.349) : Vector3d(0.6588, 0.6588, 0.0862);
		if (show_mesh_quality) color = colorbar.getColor(1 - MeshMath::get_min_scaled_Jacobian_HexMesh(*hexMesh, cp.first));
		if (drawFlag) {
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				auto f_vhs = hexMesh->faces(fh).getVertexHandle();
				for (int i = 2; i < f_vhs.size(); ++i) {
					std::vector<VH> trivhs = { f_vhs[0], f_vhs[i - 1], f_vhs[i] };
					for (auto& f_vh : trivhs) {
						auto& f_v = hexMesh->vertices(f_vh);
						current_drawed[f_vh] = true;
						vertices_data.push_back(f_v.x());
						vertices_data.push_back(f_v.y());
						vertices_data.push_back(f_v.z());
						vertices_data.push_back(color[0]);
						vertices_data.push_back(color[1]);
						vertices_data.push_back(color[2]);
					}
				}
				for (int i = 0; i < f_vhs.size(); ++i) {
					auto& v0 = hexMesh->vertices(f_vhs[i]);
					auto& v1 = hexMesh->vertices(f_vhs[(i + 1) % f_vhs.size()]);
					line_vertices_data.push_back(v0.x());
					line_vertices_data.push_back(v0.y());
					line_vertices_data.push_back(v0.z());
					line_vertices_data.push_back(0);
					line_vertices_data.push_back(0);
					line_vertices_data.push_back(0);
					line_vertices_data.push_back(v1.x());
					line_vertices_data.push_back(v1.y());
					line_vertices_data.push_back(v1.z());
					line_vertices_data.push_back(0);
					line_vertices_data.push_back(0);
					line_vertices_data.push_back(0);
				}

			}
		} else {// 被截掉，放入透明绘制数据
			if (!hexMesh->isOnBoundary(cp.first)) continue;
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				if (!hexMesh->isOnBoundary(fh)) continue;
				auto f_vhs = hexMesh->faces(fh).getVertexHandle();
				for (int i = 2; i < f_vhs.size(); ++i) {
					std::vector<VH> trivhs = { f_vhs[0], f_vhs[i - 1], f_vhs[i] };
					for (auto& f_vh : trivhs) {
						auto& f_v = hexMesh->vertices(f_vh);
						translucence_vertices_data.push_back(f_v.x());
						translucence_vertices_data.push_back(f_v.y());
						translucence_vertices_data.push_back(f_v.z());
						translucence_vertices_data.push_back(0.7725);
						translucence_vertices_data.push_back(0.7725);
						translucence_vertices_data.push_back(0.7725);
					}
				}
			}

		}

	}

	tri_cnt = vertices_data.size() / 18;
	line_cnt = line_vertices_data.size() / 12;
	translucence_tricnt = translucence_vertices_data.size() / 18;

	vertices_data.insert(vertices_data.end(), line_vertices_data.begin(), line_vertices_data.end());
	vertices_data.insert(vertices_data.end(), translucence_vertices_data.begin(), translucence_vertices_data.end());

	update();

}

void ModelDraw::updateDrawData_SurfaceMesh_CutPlane() {
	if (!this->triMesh || quad_cnt != 0) return;
	vertices_data.clear();
	triMesh->genAllVerticesNormal();
	/*double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];*/
	//printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	auto& faces = triMesh->allfaces();
	for (auto& fp : faces) {
		bool drawFlag = true;
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = triMesh->vertices(vh);
			/*if (v.x() < draw_min_x || v.y() < draw_min_y || v.z() < draw_min_z) {
				drawFlag = false;
				break;
			}*/
		}
		if (drawFlag) {
			for (auto& f_vh : vhs) {
				auto& f_v = triMesh->vertices(f_vh);
				vertices_data.push_back(f_v.x());
				vertices_data.push_back(f_v.y());
				vertices_data.push_back(f_v.z());
				vertices_data.push_back(f_v.getNormalX());
				vertices_data.push_back(f_v.getNormalY());
				vertices_data.push_back(f_v.getNormalZ());
			}
		}

	}
	update();

}

void ModelDraw::updateMeshModel() {// 一般只在重新加载模型的时候使用
	if (drawmodule == SURFACEMESH) {
		assert(triMesh != nullptr);
		this->triMesh->initBBox();
		this->bbMin = Vector3d(triMesh->BBoxMin.x(), triMesh->BBoxMin.y(), triMesh->BBoxMin.z());
		this->bbMax = Vector3d(triMesh->BBoxMax.x(), triMesh->BBoxMax.y(), triMesh->BBoxMax.z());
	} else if (drawmodule == HEXMESH || drawmodule == SKEL) {
		assert(hexMesh != nullptr);
		this->hexMesh->initBBox();
		this->bbMin = Vector3d(hexMesh->BBoxMin.x(), hexMesh->BBoxMin.y(), hexMesh->BBoxMin.z());
		this->bbMax = Vector3d(hexMesh->BBoxMax.x(), hexMesh->BBoxMax.y(), hexMesh->BBoxMax.z());
	} else if (drawmodule == QUADMESH) {
		assert(quadMesh != nullptr);
		this->quadMesh->initBBox();
		this->bbMin = Vector3d(quadMesh->BBoxMin.x(), quadMesh->BBoxMin.y(), quadMesh->BBoxMin.z());
		this->bbMax = Vector3d(quadMesh->BBoxMax.x(), quadMesh->BBoxMax.y(), quadMesh->BBoxMax.z());
	} else if (drawmodule == TETMESH) {
		assert(tetMesh != nullptr);
		this->tetMesh->initBBox();
		this->bbMin = Vector3d(tetMesh->BBoxMin.x(), tetMesh->BBoxMin.y(), tetMesh->BBoxMin.z());
		this->bbMax = Vector3d(tetMesh->BBoxMax.x(), tetMesh->BBoxMax.y(), tetMesh->BBoxMax.z());
	} else if (drawmodule == MMC_Convert) {
		if (axis.empty()) {
			double min_x = 999999999, min_y = 999999999, min_z = 999999999;
			double max_x = -999999999, max_y = -999999999, max_z = -999999999;
			assert(poly_mesh);
			poly_mesh->request_vertex_status();
			for (auto v_it = poly_mesh->vertices_begin(); v_it != poly_mesh->vertices_end(); ++v_it) {
				if (poly_mesh->status(*v_it).deleted()) continue;
				auto V = poly_mesh->point(*v_it);
				auto P = V.data();
				double x = P[0], y = P[1], z = P[2];
				min_x = std::min(x, min_x), min_y = std::min(y, min_y), min_z = std::min(z, min_z);
				max_x = std::max(x, max_z), max_y = std::max(y, max_y), max_z = std::max(z, max_z);
			}
			this->bbMin = Vector3d(min_x, min_y, min_z);
			this->bbMax = Vector3d(max_x, max_y, max_z);
		}
	}
	
	Vector3d center = (bbMin + bbMax) / 2;
	double radius = (bbMax - bbMin).Length();
	SetScenePosition(center, radius);

	/*MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	double Avg_radius = 0;
	for (auto& vmap : curve_skel_->allvertices())
	{
		auto& vh = vmap.first;
		Avg_radius += curve_skel_->vertex_radius_[vh];
	}
	Avg_radius = Avg_radius / curve_skel_->VertexSize();*/

	//local_coordinater_axis->set_len(Avg_radius * 2);
	local_coordinater_axis->set_len(1);
}

void ModelDraw::smoothing_hexahedron(int iter) {
	if (!this->hexMesh) return;
	
	updateDrawData();
	update();
}

void ModelDraw::subdivision(uint16_t iter) {

	if (!this->hexMesh) return;
	Subdivision_CC app(*this->hexMesh);
	app.subdivision_ordinary(iter);

	updateDrawData();
	update();

}

void ModelDraw::regularization_hexahedron(bool use_exterior_only, int max_iter) {
	if (!this->hexMesh) return;

	HexMesh_Regularization app(*this->hexMesh, use_exterior_only);
	app.regularization(max_iter);

	updateDrawData();
	update();
}

int ModelDraw::getNumOfModules() {

	if (hexMesh) return MeshVisiter::get_num_of_modules(hexMesh);
	else if (tetMesh) return MeshVisiter::get_num_of_modules(tetMesh);
	else if (triMesh) return MeshVisiter::get_num_of_modules(triMesh);
	else if (quadMesh) return MeshVisiter::get_num_of_modules(quadMesh);

	return 0;
}

void ModelDraw::remeshing() {

	if (!triMesh) return;


	bool input_ok = false;
	double edge_length_ratio = QInputDialog::getDouble(this, " Remeshing ", "Please input edge length ratio",
		1.f, 0.1f, 2.f, 2, &input_ok);
	if (!input_ok) return;

	std::vector<VH> feature_vhs_processing;

	for (auto& vid : spc_id) feature_vhs_processing.push_back(VH(vid.first));
	for (auto& vid : force_id) feature_vhs_processing.push_back(VH(vid.first));

	Remeshing app(*triMesh, true, feature_vhs_processing, edge_length_ratio);
	app.Execute();

	updateAllHandles_kepp_constraint(*triMesh, force_id, spc_id, constrained_point3);
	updateDrawData();

}


void ModelDraw::drawSelectedObject() 
{
	
	shaderProgram.setUniformValue("colormode", 510);// Blue
	glDepthRange(0.01, 0.02);
	if (!selected_rectangle.empty() && editmodel == EditModel::PICKITEM) 
	{
		// 矩形
		glLineWidth(meshline_width);
		glBegin(GL_LINES);
		for (auto& vex : selected_rectangle)
		{
			glVertex3d(vex.x(), vex.y(), vex.z());
		}
		glEnd();
	}

	shaderProgram.setUniformValue("colormode", 999);// Red
	glDepthRange(0.01, 0.99999);

	switch (drawmodule) 
	{
	
	case DrawModule::SURFACEMESH:
		drawSelectedObject_triMesh();
		break;
	case DrawModule::HEXMESH:
		drawSelectedObject_hexMesh();
		break;
	case DrawModule::TETMESH:
		drawSelectedObject_tetMesh();
		break;
	case DrawModule::QUADMESH:
	//case DrawModule::BSPLINESURFACE:
		drawSelectedObject_quadMesh();
		break;
	case DrawModule::SKEL:
	//case DrawModule::Skeleton_edit_move:
		drawSelectedObject_SKEL();
		break;
	default:
		break;
	}

	glDepthRange(0.01, 1.0);
	update();
}

void ModelDraw::Translated(double x, double y, double z) {
	QMatrix4x4 T = QMatrix4x4(
		1, 0, 0, x,
		0, 1, 0, y,
		0, 0, 1, z,
		0, 0, 0, 1);
	translate_matrix_ = T * translate_matrix_;													// 左乘后操作
	model_matrix_ = translate_matrix_ * rotate_matrix_;
}

void ModelDraw::SetModelMatrix(ModelDraw& modelDraw, QMatrix4x4 matrix) 
{
	shaderProgram.setUniformValue("model", modelDraw.modelmatrix * matrix);	// 右乘先操作
}

void ModelDraw::draw_skeleton_selected_point()
{
	//QVector3D point_pos;
	if ((drawSelectedAllFlag || pickType == Skeleton_edit_move) && selected_v_draw_indices.size())
	{
		shaderProgram.setUniformValue("colormode", 3);
		shaderProgram.setUniformValue("color", QVector4D(1, 0, 0, 1));
		QVector3D point_pos;
		
		MeshKernel::CurveSkel* curve_skel_ = nullptr;
		if (proDraw->skel_tool_) curve_skel_ = proDraw-> skel_tool_-> GetProcessSkel();
		if (!curve_skel_) return;

		glPointSize(17);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_draw_indices)
		{ //绘制点
			auto& v = curve_skel_->vertices(MeshKernel::iGameVertexHandle(v_idx));
			//shaderProgram.setUniformValue("color", QVector4D(0 / 255, 0, 1.0, 1));
			glVertex3d(v.x(), v.y(), v.z());
			point_pos.setX(v.x()); 
			point_pos.setY(v.y());
			point_pos.setZ(v.z());
			//DrawDiscretizeSphere(point_pos.x(), point_pos.y(), point_pos.z(), 5); // 画离散球
		}
		glEnd();


	}
}

void ModelDraw::DrawDiscretizeSphere(double Origin_x, double Origin_y, double Origin_z, double radius)
{
	shaderProgram.setUniformValue("colormode", 3);
	shaderProgram.setUniformValue("color", QVector4D(0 / 255, 0, 1.0, 1));

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	unsigned int sphereVAO = 0;
	unsigned int indexCount;
	
	if (sphereVAO == 0)
	{
		core->glGenVertexArrays(1, &sphereVAO);

		unsigned int vbo, ebo;
		core->glGenBuffers(1, &vbo);
		core->glGenBuffers(1, &ebo);

		std::vector<glm::vec3> positions;
		std::vector<glm::vec2> uv;
		std::vector<glm::vec3> normals;
		std::vector<unsigned int> indices;

		const unsigned int X_SEGMENTS = 64;
		const unsigned int Y_SEGMENTS = 64;
		//float PI = 3.14159265359;
		for (unsigned int y = 0; y <= Y_SEGMENTS; ++y)
		{
			for (unsigned int x = 0; x <= X_SEGMENTS; ++x)
			{
				float xSegment = (float)x / (float)X_SEGMENTS;
				float ySegment = (float)y / (float)Y_SEGMENTS;
				
				//这里涉及到函数的变量
				float xPos = radius * std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin_x;
				float yPos = radius * std::cos(ySegment * PI) + Origin_y;
				float zPos = radius * std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin_z;

				positions.push_back(glm::vec3(xPos, yPos, zPos));
				uv.push_back(glm::vec2(xSegment, ySegment));
				normals.push_back(glm::vec3(xPos, yPos, zPos));
			}
		}

		bool oddRow = false;
		for (unsigned int y = 0; y < Y_SEGMENTS; ++y)
		{
			if (!oddRow) // even rows: y == 0, y == 2; and so on
			{
				for (unsigned int x = 0; x <= X_SEGMENTS; ++x)
				{
					indices.push_back(y * (X_SEGMENTS + 1) + x);
					indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
				}
			}

			// 这里奇偶分开添加是有道理的，奇偶分开添加，就能首位相连，自己可以拿笔画一画
			else
			{
				for (int x = X_SEGMENTS; x >= 0; --x)
				{
					indices.push_back((y + 1) * (X_SEGMENTS + 1) + x);
					indices.push_back(y * (X_SEGMENTS + 1) + x);
				}
			}
			oddRow = !oddRow;
		}
		indexCount = indices.size();

		std::vector<float> data;
		for (unsigned int i = 0; i < positions.size(); ++i)
		{
			data.push_back(positions[i].x);
			data.push_back(positions[i].y);
			data.push_back(positions[i].z);
			if (uv.size() > 0)
			{
				data.push_back(uv[i].x);
				data.push_back(uv[i].y);
			}
			if (normals.size() > 0)
			{
				data.push_back(normals[i].x);
				data.push_back(normals[i].y);
				data.push_back(normals[i].z);
			}
		}
		core->glBindVertexArray(sphereVAO);
		core->glBindBuffer(GL_ARRAY_BUFFER, vbo);
		core->glBufferData(GL_ARRAY_BUFFER, data.size() * sizeof(float), &data[0], GL_STATIC_DRAW);
		core->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ebo);
		core->glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);
		float stride = (3 + 2 + 3) * sizeof(float);
		core->glEnableVertexAttribArray(0);
		core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
		core->glEnableVertexAttribArray(1);
		core->glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
		core->glEnableVertexAttribArray(2);
		core->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, stride, (void*)(5 * sizeof(float)));
	}

	core->glBindVertexArray(sphereVAO);
	core->glDrawElements(GL_TRIANGLE_STRIP, indexCount, GL_UNSIGNED_INT, 0);

}


void ModelDraw::draw_possible_add_point()
{
	if ((drawSelectedAllFlag || pickType == Skeleton_edit_add) )
	{

		MeshKernel::CurveSkel* curve_skel_ = nullptr;
		if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
		if (!curve_skel_) return;

		if (is_possible_add_point_draw && is_copy_point_draw == false && selected_add_possible_point.size() != 0)
		{
			glPointSize(5);
			shaderProgram.setUniformValue("colormode", 500); // 黑色
			glBegin(GL_POINTS);
			//for (auto& v_idx : selected_v_indices) { //绘制点
			//	auto& v = curve_skel_->vertices(MeshKernel::iGameVertexHandle(v_idx));
			//	glVertex3d(v.x(), v.y(), v.z());
			glVertex3d(add_possible_point.x(), add_possible_point.y(), add_possible_point.z());
			glEnd();
		}

		////bug1 这里是找不到边
		if (is_line_draw && is_copy_point_draw == false)
		{
			glLineWidth(2.5);
			shaderProgram.setUniformValue("colormode", 3); // 
			shaderProgram.setUniformValue("color", QVector4D(0 / 255, 0, 1.0, 0.5));
			glBegin(GL_LINES);
			for (auto& e_idx : selected_e_indices) 
			{ //绘制边
				auto& e = curve_skel_->edges(MeshKernel::iGameEdgeHandle(e_idx));
				auto& v1 = curve_skel_->vertices(e.vh1());
				auto& v2 = curve_skel_->vertices(e.vh2());
				glVertex3d(v1.x(), v1.y(), v1.z());
				glVertex3d(v2.x(), v2.y(), v2.z());
			}
			glEnd();
		}

		//bug2 这里是找不到点
		if (is_copy_point_draw && copy_point != -1)
		{
			glPointSize(12);
			glBegin(GL_POINTS);
			//for (auto& v_idx : selected_v_indices) { //绘制点
			//	auto& v = curve_skel_->vertices(MeshKernel::iGameVertexHandle(v_idx));
			//	glVertex3d(v.x(), v.y(), v.z());
			//printf("in the draw copy_point is %d\n", copy_point);
			glVertex3d(curve_skel_->vertices(copy_point).x(), curve_skel_->vertices(copy_point).y(), curve_skel_->vertices(copy_point).z());
			glEnd();
		}
		

	}
}



void ModelDraw::drawSelectedObject_SKEL() 
{
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) 
	{
		glPointSize(17);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) { //绘制点
			auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();

		//glLineWidth(17);
		//glBegin(GL_LINES);
		//for (auto& e_idx : selected_e_indices) { //绘制边
		//	auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
		//	auto& v1 = hexMesh->vertices(e.vh1());
		//	auto& v2 = hexMesh->vertices(e.vh2());
		//	glVertex3d(v1.x(), v1.y(), v1.z());
		//	glVertex3d(v2.x(), v2.y(), v2.z());
		//}
		//glEnd();
	}

	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) 
	{
		glLineWidth(7);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = hexMesh->vertices(e.vh1());
			auto& v2 = hexMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}

	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) 
	{
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& c_idx : selected_c_indices) {
			auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}

}

void ModelDraw::drawSelectedObject_hexMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) 
	{
		glPointSize(17);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) 
		{
			auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}

	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) 
	{
		glLineWidth(7);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = hexMesh->vertices(e.vh1());
			auto& v2 = hexMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size())
	{
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) 
		{
			auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) 
	{
		glBegin(GL_QUADS);
		for (auto& c_idx : selected_c_indices) 
		{
			auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}
}

void ModelDraw::drawSelectedObject_triMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = triMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(1.5);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = triMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = triMesh->vertices(e.vh1());
			auto& v2 = triMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_TRIANGLES);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = triMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = triMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
}

void ModelDraw::drawSelectedObject_quadMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = quadMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(1.5);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = quadMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = quadMesh->vertices(e.vh1());
			auto& v2 = quadMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = quadMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = quadMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
}

void ModelDraw::subdivisionNormal() {
	HexV3fMesh* temp = convent();
	BMesh::MLCASubdivision *mlcaSubdivision = new BMesh::MLCASubdivision();
	mlcaSubdivision->Spliting(*temp);
	free(mlcaSubdivision);
	FileHandler fileOut("crose_meshs.mesh");
	fileOut.writeFile("crose_meshs.mesh", temp);
	// 读取文件
	loadMesh("crose_meshs.mesh");
}

void ModelDraw::edgeSubdivision(int divnum) {
	if (selected_e_indices.size() != 1) return;
	static EdgeSubdivision eageSubdivision = EdgeSubdivision();
	for (int i = 0; i < selected_e_indices.size(); i++) {
		eageSubdivision.Spliting(hexMesh, selected_e_indices[i], divnum);
	}
	//hexMesh->updateAllHandles();
	selected_e_indices.clear();
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	update();
}


void ModelDraw::PIA() {
	HexV3fMesh* temp_hexMea = convent();
	HexV3fMesh* temp_hexMeb = convent();
	temp_hexMea->examine_bdy();
	temp_hexMea->examine_cs();
	temp_hexMea->cal_cen();
	temp_hexMea->examine_bdfv();
	temp_hexMea->cal_cell_cen();
	temp_hexMea->cal_bdnorm();
	temp_hexMea->cal_mami_ed();
	temp_hexMea->cal_oriv();

	temp_hexMeb->examine_bdy();
	temp_hexMeb->examine_cs();
	temp_hexMeb->cal_cen();
	temp_hexMeb->examine_bdfv();
	temp_hexMeb->cal_cell_cen();
	temp_hexMeb->cal_bdnorm();
	temp_hexMeb->cal_mami_ed();
	temp_hexMeb->cal_oriv();
	for (int i = 0; i < 10; ++i) {
		PIA_FOR_CC_SOLID pia(temp_hexMea, temp_hexMea);
		pia.oldmesh = temp_hexMeb;

		std::cout << "stdd" << std::endl;
		//ui.textBrowser->append(QString::number(pia.err_now));
		pia.PIA_FOR_CC_SOLID_ONCE();


		temp_hexMea->examine_bdy();
		temp_hexMea->cal_cen();
		temp_hexMea->examine_bdfv();
		temp_hexMea->cal_cell_cen();
		temp_hexMea->cal_bdnorm();
		temp_hexMea->cal_mami_ed();

	}
	FileHandler fileOut("crose_meshs.mesh");
	fileOut.writeFile("crose_meshs.mesh", temp_hexMea);
	// 读取文件
	loadMesh("crose_meshs.mesh");
	//smoothing(2, false);
	std::cout << "[DEBUG LEE] " << std::endl;
	updateDrawData();
	update();
}

void ModelDraw::project() {
	static int f = 0;
	if (f == 0) {
		f++;
		projectSkel->hexMesh = hexMesh;
		projectSkel->projectVolume();
		//test
		//SetEditModel(EditModel::TEST);
	}
	else {
		projectSkel->moveP();
		f = 0;
	}
	updateDrawData();
	update();
}

void ModelDraw::rselectALine() {
	if (selected_c_indices.size() == 1) {

		iGameVertex tmp(0.0, 0.0, 0.0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices.back());
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			tmp += hexMesh->vertices(vh);
		}
		tmp = tmp / 8.0;
		// check with all node position
		BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(tmp[0], tmp[1], tmp[2]));
		std::vector<BMesh::Node*> ret;
		solutionSkel->getALine(ret);
		for (int i = 0; i < ret.size(); i++) {
			iGameVertex temp_pick(ret[i]->getPosition().x, ret[i]->getPosition().y, ret[i]->getPosition().z);
			int csize = hexMesh->csize();
			double mindis = 999; int mineid = 0;
			for (int i = 0; i < csize; i++) {
				iGameVertex meshv(0, 0, 0);
				auto& cell = hexMesh->cells((iGameCellHandle)i);
				auto vhs = cell.getVertexHandle();
				for (auto& vh : vhs) {
					meshv += hexMesh->vertices(vh);
				}
				meshv = meshv / 8;
				double dis = meshv.dist(temp_pick);
				if (dis < mindis) {
					mindis = dis;
					mineid = i;
				}
			}
			// find it
			std::vector<iGameCellHandle>::iterator it;
			if ((std::find(selected_c_indices.begin(), selected_c_indices.end(), iGameCellHandle(mineid))) != selected_c_indices.end()) {
				// pickcell.erase(it);
			}
			else {
				selected_c_indices.push_back(iGameCellHandle(mineid));
			}
		}
	}
	updateMeshModel();
	updateDrawData();
	update();
}

void ModelDraw::setRoot(BMesh::Node* _root) {
	root = _root;
}

void ModelDraw::setNodeVec(std::vector<BMesh::Node*>& v) {
	nodeVec = v;
}

void ModelDraw::setNodeVec(const std::vector<BMesh::Node*>& v) {
	nodeVec = v;
}

void ModelDraw::rconnectNodes() {
	std::vector<BMesh::Node*> vv;
	for (int i = 0; i < selected_c_indices.size(); i++) {
		iGameVertex meshv(0, 0, 0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[i]);
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = hexMesh->vertices(vh);
			meshv = meshv + v;
		}
		BMesh::Node* node1 = solutionSkel->selectWhichNode(V3f(meshv[0], meshv[1], meshv[2]));
		if (node1 != nullptr) {
			vv.push_back(node1);
		}
	}
	// connect vv
	for (int i = 0; i < vv.size() - 1 && vv.size(); i++) {
		for (int j = i + 1; j < vv.size(); j++) {
			solutionSkel->smartConnectNodes(vv[i], vv[j]);
		}
	}
	setNodeVec(solutionSkel->getNodeVec());
	update();
}

void ModelDraw::rdisconnectNodes() {
	std::vector<BMesh::Node*> vv;
	for (int i = 0; i < selected_c_indices.size(); i++) {
		iGameVertex meshv(0, 0, 0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[i]);
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			meshv = meshv + hexMesh->vertices(vh);
		}
		meshv = meshv / 8.0;
		BMesh::Node* node1 = solutionSkel->selectWhichNode(V3f(meshv[0], meshv[1], meshv[2]));
		if (node1 != nullptr) {
			vv.push_back(node1);
		}
	}
	// connect vv
	for (int i = 0; i < vv.size() - 1 && vv.size(); i++) {
		for (int j = i + 1; j < vv.size(); j++) {
			solutionSkel->smartConnectNodes(vv[i], vv[j]);
		}
	}
	setNodeVec(solutionSkel->getNodeVec());
	update();
}

void ModelDraw::BFS(BMesh::Node* root) {
	Q.push(root);
	V3f v1;
	V3f v2;
	glLineWidth(3);
	shaderProgram.setUniformValue("colormode", 235);
	while (!Q.empty()) {
		BMesh::Node* a = Q.front();
		visited[a->id] = true;
		v1.x = a->root->position.x;
		v1.y = a->root->position.y;
		v1.z = a->root->position.z;
		Q.pop();
		for (int i = 0; i < a->children.size(); i++) {
			if (visited[a->children[i]->id] == false) {
				Q.push(a->children[i]);
				v2.x = a->children[i]->root->position.x;
				v2.y = a->children[i]->root->position.y;
				v2.z = a->children[i]->root->position.z;
				glBegin(GL_LINES);

#if 1
				GLfloat color[] = { 0,0,1,0.5 }; // GLfloat color[] = { 0.5,0.5,0.5,0.1 }
				GLfloat shine[] = { 120.0f };
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
#endif
				glVertex3d(v1.x, v1.y, v1.z); glVertex3d(v2.x, v2.y, v2.z);
				glEnd();
			}
		}
	}
	glLineWidth(1);
}

void ModelDraw::BFSTraverse(BMesh::Node* root) {
	visited.clear();
	std::queue<BMesh::Node*> empty;
	std::swap(empty, Q);

	visited[root->id] = true;
	BFS(root);
	for (int i = 0; i < nodeVec.size(); i++) {
		if (!visited[nodeVec[i]->id]) {
			BFS(nodeVec[i]);
		}
	}
}

void ModelDraw::drawAxis() {
	double len = radius * 10;
	glLineWidth(2);
	shaderProgram.setUniformValue("colormode", 232);
	//glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, len, 0);
	glEnd();
	shaderProgram.setUniformValue("colormode", 233);
	//glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(len, 0, 0);
	glEnd();
	shaderProgram.setUniformValue("colormode", 234);
	//glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, len);
	glEnd();
}

void ModelDraw::drawTriMesh_skel() {
	glLineWidth(2);
	shaderProgram.setUniformValue("colormode", 236);
	glBegin(GL_LINES);
	for (auto e_it = projectSkel->triMeshtv.edges_begin(); e_it != projectSkel->triMeshtv.edges_end(); ++e_it)
	{
		// 得到边所代表的半边
		OpenMesh::HalfedgeHandle heh1 = projectSkel->triMeshtv.halfedge_handle(*e_it, 0); // 默认一个方向的半边
		OpenMesh::Vec3d edgePoint(0, 0, 0);
		int edgePointsNumber = 0;
		OpenMesh::DefaultTraits::Point pointV = projectSkel->triMeshtv.point(projectSkel->triMeshtv.from_vertex_handle(heh1)); // 这条（半）边的起点
		OpenMesh::DefaultTraits::Point pointW = projectSkel->triMeshtv.point(projectSkel->triMeshtv.to_vertex_handle(heh1));   // 这条（半）边的终点
		glVertex3f(pointV[0], pointV[1], pointV[2]);
		glVertex3f(pointW[0], pointW[1], pointW[2]);
	}
	glEnd();
}

void ModelDraw::drawSkel() {
	BFSTraverse(root);
}

void ModelDraw::drawProjectLine() {
	std::vector<std::pair<V3f, V3f>> as = projectSkel->getProjectLine();
	glLineWidth(7);
	shaderProgram.setUniformValue("colormode", 3);
	shaderProgram.setUniformValue("color", QVector4D(0, 0, 0, 1));
	glBegin(GL_LINES);
	for (auto i = as.begin(); i < as.end(); i++) {
		glVertex3f(i->first.x, i->first.y, i->first.z);
		glVertex3f(i->second.x, i->second.y, i->second.z);
	}
	glEnd();
}

void  ModelDraw::drawHorizon() {
	float stepd = radius * 0.05;
	shaderProgram.setUniformValue("colormode", 231);
	//glColor3f(164 / 256.0, 128 / 256.0, 251 / 256.0);
	const int step = 1;
	const int maxLen = 100;
	glLineWidth(1);
	glBegin(GL_LINES);
	for (int i = -maxLen; i < maxLen - step; i += step) {
		for (int j = -maxLen; j < maxLen - step; j += step) {
			glVertex3f(i * stepd, j * stepd, 0);
			glVertex3f(i * stepd, (j + step) * stepd, 0);
			glVertex3f(i * stepd, j * stepd, 0);
			glVertex3f((i + step) * stepd, j * stepd, 0);
		}
	}
	glEnd();
}

int ModelDraw::loadSkeletonLee(std::string filename) {
	ModelClear();
	if (!solutionSkel->readSkel(filename)) {
		std::cout << "[DEBUG] read skel error" << std::endl;
		return -1;
	}
	std::cout << "[DEBUG] read skel success" << std::endl;
	hexMesh = new MeshKernel::VolumeMesh();
	// first gen node Hex
	solutionSkel->genNodeHex();
	// let Bmesh hex to mesh
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	updateMeshModel();
	//updateDrawData_HexMesh();
	updateDrawData();
	update();
	return 0;
}

void ModelDraw::drawSelectedObject_tetMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = tetMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(1.5);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = tetMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = tetMesh->vertices(e.vh1());
			auto& v2 = tetMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_TRIANGLES);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = tetMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = tetMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
		glBegin(GL_TRIANGLES);
		for (auto& c_idx : selected_c_indices) {
			auto fhs = tetMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = tetMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle(); 
				for (auto& vh : vhs) {
					auto& v = tetMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}
}


void ModelDraw::DrawPointRadius()
{
	if (Move_Point_vh.size() <= 0)
	{
		return;
	}
	int n = 1000;
	vector<QVector3D> Points;
	double  Copy = PI / 2 / n;
	//radius = selected_v_radius * 2;

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	//std::cout << "!!!!!" << std::endl;

	auto origin_v = curve_skel_->vertices(Move_Point_vh[0]);
	//asdkjhgaskdjhask;
	//以坐标轴的长度来画 半径
	radius =  curve_skel_->vertex_radius_[Move_Point_vh[0]] ;
	//radius = local_coordinater_axis->ReturnAixsLength() * 0.4 * RadiusScaling * 3 * curve_skel_->vertex_radius_[Move_Point_vh[0]];

	for (int i = 0; i <= n; i++)
	{
		QVector3D temp(origin_v.x() + radius * cos(Copy * i), origin_v.y() + radius * sin(Copy * i), origin_v.z() + 0);
		Points.push_back(temp);
	}

	/*glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(len * k, 0, 0);
	glEnd();*/

	glPointSize(5);
	shaderProgram.setUniformValue("colormode", 3);
	

	if (is_change_radius)
		shaderProgram.setUniformValue("color", QVector4D(0 / 255, 255 / 255, 255 / 255, 1));
	else
		shaderProgram.setUniformValue("color", QVector4D(0 / 255, 255 / 255, 0 / 255, 1));

	glDepthRange(0.0, 0.01);
	//glBegin(GL_POINTS);
	glBegin(GL_LINES);
	for (int i = 1; i <= n; i++)
	{
		//glVertex3f(Points[i].x(), Points[i].y(), Points[i].z());
		//glBegin(GL_LINES);
		glVertex3f(Points[i - 1].x(), Points[i - 1].y(), Points[i - 1].z());
		glVertex3f(Points[i].x(), Points[i].y(), Points[i].z());
		//glEnd();

	}
	glEnd();

}

/**************************** 场景绘制 ******************************************/

void ModelDraw::DrawScene(void) // update() 或者 一帧就会刷新 
{
	DrawStyle();
	ColorStyle();
	DrawSceneModel();
}

void ModelDraw::DrawSceneModel(void) 
{

	//glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

	// 绘制点, 用于剪断网格线框以保证框选准确性
	if (editmodel == EditModel::PICKITEM && mouse_move_flag) 
	{
		shaderProgram.setUniformValue("colormode", 999); //231
		glDepthRange(0.01, 0.9999);
		glPointSize(meshline_width * 5.f);
		core->glDrawArrays(GL_POINTS, 0, vertices_data.size() / 6); // 这个地方不是很懂？
	}

	// 坐标轴
	if (editmodel == EditModel::POINT_COORDINATE) // ?好像是移动按钮
	{
		if (select_coordinater_mode == CoordinateAxis::ModeType::TRANSFROM && !selected_v_indices.empty())
		{
			local_coordinater_axis->set_mode_type(select_coordinater_mode); // 设置模式为 TRANSFROM
			Vec& v = hexMesh->vertices((VH)selected_v_indices.back());
			local_coordinater_axis->set_orign(v.x(), v.y(), v.z()); // 设置坐标轴的原点
		}

		else if (select_coordinater_mode == CoordinateAxis::ModeType::CREATE && !selected_f_indices.empty()) 
		{
			local_coordinater_axis->set_mode_type(select_coordinater_mode);
			iGameFace& f = hexMesh->faces((FH)selected_f_indices.back());
			Vec center = Vec(0, 0, 0);
			for (auto vh : f.getVertexHandle())
				center = (center + hexMesh->vertices(vh)/4.0);// ?
			local_coordinater_axis->set_orign(center.x(), center.y(), center.z());
		}
		else
			local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
	}

	

	DrawPointRadius();
	local_coordinater_axis->Draw3DCoord(*this);

	
	proDraw->DrawSceneModel(); // 用于画骨架的  先画不透明的再画透明的
	if (proDraw->skel_tool_)
	{
		curve_skel_ForDraw = proDraw->skel_tool_->GetProcessSkel();
		//if (!curve_skel_) return;
	}
	//这里导致了渲染的问题先屏蔽掉了
	/*MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;*/

	//for (auto vmap : curve_skel_->allvertices())
	//{
	//	//updateCurrentSkelSphere(vmap.first); // 传入 data
	//	loadShader_Skel_Sphere(vmap.first); // 传入 VAO
	//	//DrawSkelDiscretize(vmap.first);
	//}

	//for (auto emap : curve_skel_->alledges())
	//{
	//	//updateCurrentSkel(vmap.first); // 传入 data
	//	loadShader_Skel_Cylinder(emap.first); // 传入 VAO
	//	DrawSkelDiscretize_Cylinder(emap.first);
	//}
	

	loadShader();
	switch (drawmodule)
	{
	case EMPTY:
		break;

	case SURFACEMESH:// 曲面网格绘制
		drawSurfaceMesh();
		if (IsReadyCloudPoint)
		{
			drawMeshCloudPoint();
		}
		//drawMeshCloudPoint();
		//drawMeshCloudPoint();
		break;

	case HEXMESH:// 六面体网格绘制
		//loadShader();
		TestloadShader();
		TestDrawSurfaceMesh();
		loadShader();
		drawHexMesh();
		//drawSurfaceMesh();
		break;

	case TETMESH:
		drawTetMesh();
		break;

	case QUADMESH:
		drawQuadMesh();
		break;

	case SKEL:
		drawHexMesh();
		drawProjectLine();
		break;

	case BSPLINESURFACE:// B样条曲面绘制
		if (PatchChange)
		{
			bsplineSurfaceToDiscrete();
			PatchChange = false;
		}
		//if (isControlMesh) drawControlMesh(); //  立即渲染模式的控制网格
		/*if (!ControlMeshHasFinished)
		{
			drawControlMesh();
			ControlMeshHasFinished = true;
		}*/
		
		//if (bshowTriMesh) drawSurfaceMesh_skel();
		//drawControlPointAndLine(); //
		
		drawSurface(); //要保留的 控制点和贝塞尔都可以画，去按钮那里调节画哪个
		//DrawPointForOrder();

		break;

	case MMC_Convert:
		drawPolyMesh_OPT();
		break;

	case Skeleton_edit_move:
		//draw_skeleton_selected_point();
		for (auto vmap : curve_skel_ForDraw->allvertices())
		{
			//updateCurrentSkelSphere(vmap.first); // 传入 data
			loadShader_Skel_Sphere(vmap.first); // 传入 VAO
			DrawSkelDiscretize(vmap.first);
		}

		for (auto emap : curve_skel_ForDraw->alledges())
		{
			//updateCurrentSkel(vmap.first); // 传入 data
			loadShader_Skel_Cylinder(emap.first); // 传入 VAO
			DrawSkelDiscretize_Cylinder(emap.first);
		}
		break;


	case Skeleton_edit_add:
		for (auto vmap : curve_skel_ForDraw->allvertices())
		{
			//updateCurrentSkelSphere(vmap.first); // 传入 data
			loadShader_Skel_Sphere(vmap.first); // 传入 VAO
			DrawSkelDiscretize(vmap.first);
		}
		for (auto emap : curve_skel_ForDraw->alledges())
		{
			//updateCurrentSkel(vmap.first); // 传入 data
			loadShader_Skel_Cylinder(emap.first); // 传入 VAO
			DrawSkelDiscretize_Cylinder(emap.first);
		}
		draw_possible_add_point();
		if (selected_add_possible_point.size() > 0)
		{
			loadShader_Possible_Sphere();
			DrawPossibleSphere();
		}
		break;

	default:
		break;

	}

	if (bshowSKEL) drawSkel();
	if (bshowTriMesh) drawSurfaceMesh_skel();

	drawSelectedObject(); //这个为什么绘制不出来了很奇怪

	
}


void ModelDraw::InitVaoSphereSize()
{
	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;
	
	VAOs.resize(curve_skel_->VertexSize());
	VBOs.resize(curve_skel_->VertexSize());
	data_Sphere.resize(curve_skel_->VertexSize());
}

void ModelDraw::InitVaoCylinderSize()
{
	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;
	
	VAOs_Cylinder.resize(curve_skel_->EdgeSize());
	VBOs_Cylinder.resize(curve_skel_->EdgeSize());
	data_Cylinder.resize(curve_skel_->EdgeSize());
}

void ModelDraw::InitIndices()
{
	//初始化最上面的三角形
	for (int i = 0; i < X_SEGMENTS; i++)
	{
		if (i == X_SEGMENTS - 1)
		{
			indices.push_back(0);
			indices.push_back(i + 1);
			indices.push_back(1);
		}
		else
		{
			indices.push_back(0);
			indices.push_back(i + 1);
			indices.push_back(i + 2);
		}
		
	}
	
	for (int i = 0; i < Y_SEGMENTS - 2; i++)
	{
		for (int j = 0; j < X_SEGMENTS; j++)
		{
			if (j == X_SEGMENTS - 1)
			{
				indices.push_back(i * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + 0 + 1);
				
				indices.push_back(i * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + 0 + 1);
				indices.push_back(i * X_SEGMENTS + 0 + 1);
			}
			else
			{
				indices.push_back(i * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + j + 2);

				indices.push_back(i * X_SEGMENTS + j + 1);
				indices.push_back((i + 1) * X_SEGMENTS + j + 2);
				indices.push_back(i * X_SEGMENTS + j + 2);
			}
			

		}
	}
	
	for (int i = 0; i < X_SEGMENTS; i++)
	{
		if (i == X_SEGMENTS - 1)
		{
			indices.push_back((Y_SEGMENTS - 2) * X_SEGMENTS + 1 + i);
			indices.push_back((Y_SEGMENTS - 1) * X_SEGMENTS + 1);
			indices.push_back((Y_SEGMENTS - 2) * X_SEGMENTS + 1);
		}
		else
		{
			indices.push_back((Y_SEGMENTS - 2) * X_SEGMENTS + 1 + i);
			indices.push_back((Y_SEGMENTS - 1) * X_SEGMENTS + 1);
			indices.push_back((Y_SEGMENTS - 2) * X_SEGMENTS + 2 + i);
		}

		
	}

	
}

void ModelDraw::updateCurrentSkelSphere(MeshKernel::iGameVertexHandle Node_h) //主要把 节点 Node_h数据传入到vertices数组中
{
	// 这里有一个小bug 要根据骨架节点的半径去显示球的大小
	// 
	//vertices.clear();
	positions.clear();
	//uv.clear();
	normals.clear();
	data_Sphere[Node_h].clear();
	//std::vector<float> data_pre;
	//indices.clear(); //这个可以公用吗

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	auto& Node_v = curve_skel_->vertices(Node_h);
	QVector3D Origin(Node_v.x(), Node_v.y(), Node_v.z());

	double current_radius = curve_skel_->vertex_radius_[Node_h] ; //每个点的实际半径进行一下缩放
	//double current_radius = curve_skel_->vertex_radius_[Node_h] ; //每个点的实际半径进行一下缩放

	//第一个点
	positions.push_back(glm::vec3(Origin.x(), Origin.y() + current_radius, Origin.z()));
	for (unsigned int y = 1; y < Y_SEGMENTS; ++y)
	{
		for (unsigned int x = 0; x < X_SEGMENTS; ++x)
		{
			float xSegment = (float)x / (float)X_SEGMENTS;
			float ySegment = (float)y / (float)Y_SEGMENTS;


			//这里涉及到函数的变量
			float xPos = current_radius * std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin.x();
			float yPos = current_radius * std::cos(ySegment * PI) + Origin.y();
			float zPos = current_radius * std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin.z();

			positions.push_back(glm::vec3(xPos, yPos, zPos));
		}
	}
	positions.push_back(glm::vec3(Origin.x(), Origin.y() - current_radius, Origin.z()));

	std::vector<std::vector<iGameVertexHandle>>  _faces;//用于构建面的顶点坐标

	for (int i = 0; i < indices.size(); i = i + 3)
	{
		//一个临时的面
		std::vector<iGameVertexHandle> tempface;
		for (int j = 0; j < 3; j++)
		{
			tempface.push_back(VH(indices[i + j]));
		}
		_faces.push_back(tempface);
	}

	std::vector<iGameVertex> _vertices;
	//转换成 vector<iGameVertex>
	for (int i = 0; i < positions.size(); i++)
	{
		iGameVertex tempVertex(positions[i].x, positions[i].y, positions[i].z);
		_vertices.push_back(tempVertex);
	}

	SphereMeshForSkeletonPoint.push_back(MeshKernel::SurfaceMesh(_vertices, _faces));

	
	//生成所有点的法向
	
	//for (int i = 1; i < Y_SEGMENTS; i++)
	//{
	//	//把重复的顶点删除
	//	SphereMesh.DeleteVertex(VH(i));
	//}
	//for (int i = X_SEGMENTS * Y_SEGMENTS - 1; i > (X_SEGMENTS - 1) * Y_SEGMENTS; i--)
	//{
	//	//把重复的顶点删除
	//	SphereMesh.DeleteVertex(VH(i));
	//	
	//}
	//SphereMesh.updateAllHandles();
	// 
	// 

	if (!SphereMesh_Origin)
	{
		//用于构造原始三角形网格用于比较拟合误差
		for (int i = 0; i < _vertices.size(); i++)
		{
			Vertices_Origin.push_back(_vertices[i]);
		}

		for (int i = 0; i < _faces.size(); i++)
		{
			std::vector<iGameVertexHandle> TempFaceAddToOrigin;
			//每个面
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][0]) + PreNumHandle));
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][1]) + PreNumHandle));
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][2]) + PreNumHandle));
			//PreNumHandle
			Faces_Origin.push_back(TempFaceAddToOrigin);
		}

		PreNumHandle = PreNumHandle + _vertices.size();
	}

	SphereMeshForSkeletonPoint.back().genAllVerticesNormal();

	for (auto& fmap : SphereMeshForSkeletonPoint.back().allfaces())
	{
		for (auto& vh : fmap.second.getVertexHandle())
		{
			//把每个面都压进去
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).x());
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).y());
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).z());

			//法向信息 顺时针传进去的
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).getNormalX());
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).getNormalY());
			data_Sphere[Node_h].push_back(SphereMeshForSkeletonPoint.back().vertices(vh).getNormalZ());
		}
	}

	//传入点云数据用于重构出网格
	//修复:可能要进行初始化
	//********什么情况下使用这个函数
	//if (TestNomal)
	//{
	//	CloudPointForPoissonReconstruction.add_normal_map(); // 启用法向量映射
	//	TestNomal--;
	//}
	//CloudPointForPoissonReconstruction.add_normal_map(); // 启用法向量映射
	//for (auto vmap : SphereMesh.allvertices())
	//{
	//	auto vh = vmap.first;
	//	Point_3 first_point(SphereMesh.vertices(vh).x(), SphereMesh.vertices(vh).y(), SphereMesh.vertices(vh).z());
	//	////Vector_3 first_normal(SphereMesh.vertices(vh).getNormalX(), SphereMesh.vertices(vh).getNormalY(), SphereMesh.vertices(vh).getNormalZ()); //法相的坐标
	//	////CloudPointForPoissonReconstruction.insert(first_point, first_normal);      
	//	CloudPointForPoissonReconstruction.insert(first_point);
	//	std::cout << "";
	//}
	//for (auto vmap : SphereMesh.allvertices())
	//{
	//	//思路是这样的：得到顶点的临边以及他的圆柱的数据10个点的，抽取4个点作为距离比较的点
	//	auto vh = vmap.first;
	//	Point_3 first_point(SphereMesh.vertices(vh).x(), SphereMesh.vertices(vh).y(), SphereMesh.vertices(vh).z());
	//	////Vector_3 first_normal(SphereMesh.vertices(vh).getNormalX(), SphereMesh.vertices(vh).getNormalY(), SphereMesh.vertices(vh).getNormalZ()); //法相的坐标
	//	////CloudPointForPoissonReconstruction.insert(first_point, first_normal);      
	//	CloudPointForPoissonReconstruction.insert(first_point);
	//	std::cout << "";
	//}

}

void ModelDraw::updatePossibleSphere() //主要把 节点 Node_h数据传入到vertices数组中
{
	// 这里有一个小bug 要根据骨架节点的半径去显示球的大小
	// 
	//vertices.clear();
	if (selected_add_possible_point.size() == 0)
	{
		return;
	}
	Data_Possible_Sphere.clear();
	positions.clear();
	//uv.clear();
	normals.clear();
	//std::vector<float> Data_Possible_Cylinder;
	//data_Sphere[Node_h].clear();
	//std::vector<float> data_pre;
	//indices.clear(); //这个可以公用吗

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	auto& Node_v = selected_add_possible_point[0];
	QVector3D Origin(Node_v.x(), Node_v.y(), Node_v.z());

	auto& eh_selected = selected_e_indices[0];
	auto& e_selected = curve_skel_->edges((EH)eh_selected);
	
	auto& v1 = curve_skel_->vertices(e_selected.vh1());
	auto& v2 = curve_skel_->vertices(e_selected.vh2());
	QVector3D v1_QVector(v1.x(), v1.y(), v1.z());
	QVector3D v2_QVector(v2.x(), v2.y(), v2.z());

	//double LengthZero = (Origin - v1_QVector).length() 
	double radius_Possible = ((Origin - v1_QVector).length() / (v2_QVector - v1_QVector).length() * curve_skel_->vertex_radius_[e_selected.vh2()] + (Origin - v2_QVector).length() / (v2_QVector - v1_QVector).length() * curve_skel_->vertex_radius_[e_selected.vh1()] ); //这里后续可能需要修改
	//printf("%lf\n", radius_Possible);
	

	//第一个点
	positions.push_back(glm::vec3(Origin.x(), Origin.y() + radius_Possible, Origin.z()));
	for (unsigned int y = 1; y < Y_SEGMENTS; ++y)
	{
		for (unsigned int x = 0; x < X_SEGMENTS; ++x)
		{
			float xSegment = (float)x / (float)X_SEGMENTS;
			float ySegment = (float)y / (float)Y_SEGMENTS;


			//这里涉及到函数的变量
			float xPos = radius_Possible * std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin.x();
			float yPos = radius_Possible * std::cos(ySegment * PI) + Origin.y();
			float zPos = radius_Possible * std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin.z();

			positions.push_back(glm::vec3(xPos, yPos, zPos));
		}
	}
	positions.push_back(glm::vec3(Origin.x(), Origin.y() - radius_Possible, Origin.z()));


	std::vector<std::vector<iGameVertexHandle>>  _faces;//用于构建面的顶点坐标

	for (int i = 0; i < indices.size(); i = i + 3)
	{
		//一个临时的面
		std::vector<iGameVertexHandle> tempface;
		for (int j = 0; j < 3; j++)
		{
			tempface.push_back(VH(indices[i + j]));
		}
		_faces.push_back(tempface);
	}

	std::vector<iGameVertex> _vertices;
	//转换成 vector<iGameVertex>
	for (int i = 0; i < positions.size(); i++)
	{
		iGameVertex tempVertex(positions[i].x, positions[i].y, positions[i].z);
		_vertices.push_back(tempVertex);
	}

	MeshKernel::SurfaceMesh SphereMesh(_vertices, _faces);
	//生成所有点的法向

	SphereMesh.genAllVerticesNormal();

	for (auto& fmap : SphereMesh.allfaces())
	{
		for (auto& vh : fmap.second.getVertexHandle())
		{
			//把每个面都压进去
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).x());
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).y());
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).z());

			//法向信息 顺时针传进去的
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).getNormalX());
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).getNormalY());
			Data_Possible_Sphere.push_back(SphereMesh.vertices(vh).getNormalZ());
		}


	}

}

void ModelDraw::ControlMeshChangeToQuadMesh()
{
	std::vector<std::vector<iGameVertexHandle>>  _faces;//用于构建面的顶点坐标

	int index = 0; //用于记录面的id 
	for (auto EachControlFace : data_controlMesh)
	{
		//每个面
		for (auto EachQuadOFEachControlFace : EachControlFace)
		{
			//每个面的三角形
			//std::vector<iGameVertexHandle> tempface;
			for (int i = 0; i < 2; i++)
			{
				//4边形分成2个三角形
				std::vector<iGameVertexHandle> tempface;
				//index = index + 1;
				if (i == 0)
				{
					tempface.push_back(VH(index));
					tempface.push_back(VH(index+1));
					tempface.push_back(VH(index + 2));
					_faces.push_back(tempface);
				}
				else 
				{
					tempface.push_back(VH(index));
					tempface.push_back(VH(index + 2));
					tempface.push_back(VH(index + 3));
					_faces.push_back(tempface);
				}

				
				//tempface.push_back(VH(index));
				//index = index + 1;
			}
			index = index + 4;
		}
	}

	////生成面的handle ，就是每个面是哪几个点组成的，但是如何判断它按某一个顺序排列的呢？是一个问题
	//for (int i = 0; i < data_controlMesh.size(); i = i + 3)
	//{
	//	//一个临时的面
	//	std::vector<iGameVertexHandle> tempface;
	//	for (int j = 0; j < 3; j++)
	//	{
	//		tempface.push_back(VH(indices[i + j]));
	//	}
	//	_faces.push_back(tempface);
	//}

	std::vector<iGameVertex> _vertices;
	//转换成 vector<iGameVertex>
	/*for (int i = 0; i < positions.size(); i++)
	{
		iGameVertex tempVertex(positions[i].x, positions[i].y, positions[i].z);
		_vertices.push_back(tempVertex);
	}*/
	
	for (auto EachControlFace : data_controlMesh)
	{
		//每个面
		for (auto EachQuadOFEachControlFace : EachControlFace)
		{
			//每个面的四边形
			//std::vector<iGameVertexHandle> tempface;
			for (int i = 0; i < 4; i++)
			{
				iGameVertex tempVertex(EachQuadOFEachControlFace[i].getX(), EachQuadOFEachControlFace[i].getY(), EachQuadOFEachControlFace[i].getZ());
				_vertices.push_back(tempVertex);
			}
		}
	}


	this->triMesh = new MeshKernel::SurfaceMesh(_vertices, _faces);
	//this->quadMesh = &SphereMesh;
	//cout << "";
}

void ModelDraw::drawMeshCloudPoint()
{
	glPointSize(2.5);
	shaderProgram.setUniformValue("colormode", 232);
	glBegin(GL_POINTS);
	for (auto& point : CloudPointForPoissonReconstruction.points())
	{
		// draw_->DrawDiscretizeSphere(vp.second.x(), vp.second.y(), vp.second.z(), 5);
		glVertex3f(point.x(), point.y(), point.z());

	}
	glEnd();
}

void ModelDraw::TurnToCloudPoint()
{
	SetDrawModule(SURFACEMESH);
	IsReadyCloudPoint = true;
}

void ModelDraw::GeneratetTriangularMeshforUnionCGAL()
{
	//用于remeshing的数据结构
	struct halfedge2edge
	{
		halfedge2edge(const MeshCGAL& m, std::vector<edge_descriptor>& edges)
			: m_mesh(m), m_edges(edges)
		{}
		void operator()(const halfedge_descriptor& h) const
		{
			m_edges.push_back(edge(h, m_mesh));
		}
		const MeshCGAL& m_mesh;
		std::vector<edge_descriptor>& m_edges;
	};

	vector<MeshCGAL> Mesh_CGAL;
	std::unordered_map<int, int> CGALIntVertexHandleToSphereMesh;
	std::unordered_map<int, int> CGALIntEdgeHandleToCylinderMesh;
	//这里*****写的不好
	Mesh_CGAL.resize(SphereMeshForSkeletonPoint.size() + CylinderMeshForSkeletonPoint.size());
	std::vector<MeshCGAL::Vertex_index> mesh_CGAL_vertices;
	/*for (const auto& vertex : vertices) 
	{
		Mesh::Vertex_handle v = mesh.add_vertex(vertex);
		mesh_vertices.push_back(v);
	}*/
	int indexForOrderMesh = 0;
	for (int i = 0; i < SphereMeshForSkeletonPoint.size(); i++)
	{
		mesh_CGAL_vertices.clear();
		for (auto vmap : SphereMeshForSkeletonPoint[i].allvertices())
		{
			MeshCGAL::Vertex_index vh_CGAL = Mesh_CGAL[indexForOrderMesh].add_vertex(Kernel::Point_3(SphereMeshForSkeletonPoint[i].vertices(vmap.first).x(),
				SphereMeshForSkeletonPoint[i].vertices(vmap.first).y(), SphereMeshForSkeletonPoint[i].vertices(vmap.first).z()));
			mesh_CGAL_vertices.push_back(vh_CGAL); //int = iGameHandle 和 CGALHandle映射
		}
		for (auto fmap : SphereMeshForSkeletonPoint[i].allfaces())
		{
			std::vector<MeshCGAL::Vertex_index> face_vertices_index;
			auto vectorHandle = SphereMeshForSkeletonPoint[i].faces(fmap.first).getVertexHandle();
			//GetSortedVertexHandle 这个函数有问题
			Mesh_CGAL[indexForOrderMesh].add_face(mesh_CGAL_vertices[vectorHandle[2].idx()], mesh_CGAL_vertices[vectorHandle[1].idx()], mesh_CGAL_vertices[vectorHandle[0].idx()]);
		}
		CGALIntVertexHandleToSphereMesh.emplace(SphereMeshToIntVertexHand[i], indexForOrderMesh);//建立 VextexHandle --->> CGALMeshSphere Array Index
		indexForOrderMesh++;
	}

	for (int i = 0; i < CylinderMeshForSkeletonPoint.size(); i++)
	{
		mesh_CGAL_vertices.clear();
		for (auto vmap : CylinderMeshForSkeletonPoint[i].allvertices())
		{
			MeshCGAL::Vertex_index vh_CGAL = Mesh_CGAL[indexForOrderMesh].add_vertex(Kernel::Point_3(CylinderMeshForSkeletonPoint[i].vertices(vmap.first).x(),
				CylinderMeshForSkeletonPoint[i].vertices(vmap.first).y(), CylinderMeshForSkeletonPoint[i].vertices(vmap.first).z()));
			mesh_CGAL_vertices.push_back(vh_CGAL); //int = iGameHandle 和 CGALHandle映射
		}
		for (auto fmap : CylinderMeshForSkeletonPoint[i].allfaces())
		{
			std::vector<MeshCGAL::Vertex_index> face_vertices_index;
			auto vectorHandle = CylinderMeshForSkeletonPoint[i].faces(fmap.first).getVertexHandle();
			//GetSortedVertexHandle 这个函数有问题
			Mesh_CGAL[indexForOrderMesh].add_face(mesh_CGAL_vertices[vectorHandle[0].idx()], mesh_CGAL_vertices[vectorHandle[1].idx()], mesh_CGAL_vertices[vectorHandle[2].idx()]); 
		}
		CGALIntEdgeHandleToCylinderMesh.emplace(CylinderMeshToIntEdgeHandle[i], indexForOrderMesh);////build EdgeHandle --->> CGALMeshCylinder   Array Index
		indexForOrderMesh++; 

	}

	for (int i = SphereMeshForSkeletonPoint.size(); i < SphereMeshForSkeletonPoint.size() + CylinderMeshForSkeletonPoint.size(); i++)
	{
		double target_edge_length = 0.04;
		unsigned int nb_iter = 3;
		std::cout << "Split border...";
		std::vector<edge_descriptor> border;
		PMP::border_halfedges(faces(Mesh_CGAL[i]),
			Mesh_CGAL[i],
			boost::make_function_output_iterator(halfedge2edge(Mesh_CGAL[i], border)));
		PMP::split_long_edges(border, target_edge_length, Mesh_CGAL[i]);
		std::cout << "done." << std::endl;
		std::cout <<
			" (" << num_faces(Mesh_CGAL[i]) << " faces)..." << std::endl;
		PMP::isotropic_remeshing(
			faces(Mesh_CGAL[i]),
			target_edge_length,
			Mesh_CGAL[i],
			PMP::parameters::number_of_iterations(nb_iter)
			.protect_constraints(true)//i.e. protect border, here
		);
	}
	//double target_edge_length = 0.02;
	//unsigned int nb_iter = 3;
	//std::cout << "Split border...";
	//std::vector<edge_descriptor> border;
	//PMP::border_halfedges(faces(DealMesh),
	//	DealMesh,
	//	boost::make_function_output_iterator(halfedge2edge(DealMesh, border)));
	//PMP::split_long_edges(border, target_edge_length, DealMesh);
	//std::cout << "done." << std::endl;
	//std::cout <<
	//	" (" << num_faces(DealMesh) << " faces)..." << std::endl;
	//PMP::isotropic_remeshing(
	//	faces(DealMesh),
	//	target_edge_length,
	//	DealMesh,
	//	PMP::parameters::number_of_iterations(nb_iter)
	//	.protect_constraints(true)//i.e. protect border, here
	//);


	//现在写一个BFS来进行union操作
	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;
	
	vector<bool> VisitedVertex(curve_skel_->allvertices().size(), false);
	vector<bool> VisitedEdge(curve_skel_->alledges().size(), false);
	
	std::queue<MeshKernel::iGameVertexHandle> q;
	q.push(VH(0));
	MeshCGAL MeshCGALUnion = Mesh_CGAL[CGALIntVertexHandleToSphereMesh[0]];

	while (!q.empty())
	{
		auto CureentVh = q.front();
		q.pop();

		if (visited[CureentVh])
		{
			continue;
		}

		visited[CureentVh] = true;

		//顶点没遍历过需要Union
		//遍历过的不需要
		MeshCGAL output;
		PMP::corefine_and_compute_union(Mesh_CGAL[CGALIntVertexHandleToSphereMesh[CureentVh]], MeshCGALUnion, output);
		MeshCGALUnion = output;

		for (auto EhCureent : curve_skel_->NeighborEh(CureentVh))
		{
			if (VisitedEdge[EhCureent])
			{
				continue;//边遍历过了就不union
			}
			VisitedEdge[EhCureent] = true;
			MeshCGAL output;
			PMP::corefine_and_compute_union(MeshCGALUnion, Mesh_CGAL[CGALIntEdgeHandleToCylinderMesh[EhCureent]], output);
			MeshCGALUnion = output;
		}

		for (auto vh : curve_skel_->NeighborVh(CureentVh))
		{
			if (!visited[vh])
			{
				q.push(vh);
			}
			//把相邻的顶点压入进去
			//之后要去标记哪些已经遍历过了
		} 

		//Mesh_CGAL[CGALIntVertexHandleToSphereMesh[CureentVh]] = output;
	}

	double target_edge_length = 0.035;
		unsigned int nb_iter = 3;
		std::cout << "Split border...";
		std::vector<edge_descriptor> border;
		PMP::border_halfedges(faces(MeshCGALUnion),
			MeshCGALUnion,
			boost::make_function_output_iterator(halfedge2edge(MeshCGALUnion, border)));
		PMP::split_long_edges(border, target_edge_length, MeshCGALUnion);
		std::cout << "done." << std::endl;
		std::cout <<
			" (" << num_faces(MeshCGALUnion) << " faces)..." << std::endl;
		PMP::isotropic_remeshing(
			faces(MeshCGALUnion),
			target_edge_length,
			MeshCGALUnion,
			PMP::parameters::number_of_iterations(nb_iter)
			.protect_constraints(true));//i.e. protect border, here

	// Constrain edges with a dihedral angle over 60°
	//typedef boost::property_map<MeshCGAL, CGAL::edge_is_feature_t>::type EIFMap;
	//EIFMap eif = get(CGAL::edge_is_feature, MeshCGALUnion);
	//PMP::detect_sharp_edges(MeshCGALUnion, 60, eif);

	//int sharp_counter = 0;
	//for (edge_descriptor e : edges(MeshCGALUnion))
	//	if (get(eif, e))
	//		++sharp_counter;
	//std::cout << sharp_counter << " sharp edges" << std::endl;
	//const unsigned int nb_iterations = 10;
	//std::cout << "Smoothing MeshCGALUnion... (" << nb_iterations << " iterations)" << std::endl;
	//// Smooth with both angle and area criteria + Delaunay flips
	//PMP::smooth_mesh(MeshCGALUnion, PMP::parameters::number_of_iterations(nb_iterations)
	//	.use_safety_constraints(false) // authorize all moves
	//	.edge_is_constrained_map(eif));

	const char* out_process_path_ = "D:/GraduateProgram/output.off";

	std::ofstream out(out_process_path_);
	out << MeshCGALUnion;

	
}



void ModelDraw::PoissonReconstructionUsingPointCloud()
{
	
	

	//typename Point_set::iterator rout_it = CGAL::remove_outliers<CGAL::Sequential_tag>
	//	(CloudPointForPoissonReconstruction,24, // Number of neighbors considered for evaluation
	//	CloudPointForPoissonReconstruction.parameters().threshold_percent(5.0)); // Percentage of points to remove
	//CloudPointForPoissonReconstruction.remove(rout_it, CloudPointForPoissonReconstruction.end());
	//std::cout << CloudPointForPoissonReconstruction.number_of_removed_points()<< " CloudPointForPoissonReconstruction are outliers." << std::endl;
	//// Applying point set processing algorithm to a CGAL::Point_set_3
 //// object does not erase the points from memory but place them in
 //// the garbage of the object: memory can be freed by the user.
	//CloudPointForPoissonReconstruction.collect_garbage();

	//IsReadyCloudPoint = true;

	// 计算6个点的平均距离
	double spacing = CGAL::compute_average_spacing<CGAL::Sequential_tag>(CloudPointForPoissonReconstruction, 20);
	
	
	//// Simplify using a grid of size 2 * average spacing
	//typename Point_set::iterator gsim_it = CGAL::grid_simplify_point_set(CloudPointForPoissonReconstruction, 1 * spacing);
	//CloudPointForPoissonReconstruction.remove(gsim_it, CloudPointForPoissonReconstruction.end());
	//std::cout << CloudPointForPoissonReconstruction.number_of_removed_points() << " CloudPointForPoissonReconstruction removed after simplification." << std::endl;
	//CloudPointForPoissonReconstruction.collect_garbage();

	CGAL::jet_smooth_point_set<CGAL::Sequential_tag>(CloudPointForPoissonReconstruction, 24);

	CGAL::jet_estimate_normals<CGAL::Sequential_tag>(CloudPointForPoissonReconstruction, 24); // Use 24 neighbors
	// Orientation of normals, returns iterator to first unoriented point
	//typename Point_set::iterator unoriented_points_begin = CGAL::mst_orient_normals(CloudPointForPoissonReconstruction, 24); // Use 24 neighbors

	//CloudPointForPoissonReconstruction.remove(unoriented_points_begin, CloudPointForPoissonReconstruction.end());
	CGAL::Surface_mesh<Point_3> output_mesh;
	CGAL::poisson_surface_reconstruction_delaunay(CloudPointForPoissonReconstruction.begin(), CloudPointForPoissonReconstruction.end(), CloudPointForPoissonReconstruction.point_map(), CloudPointForPoissonReconstruction.normal_map(),output_mesh, spacing);
	


	std::vector<iGameVertex> _vertices;
	for (auto v : output_mesh.vertices()) 
	{
		Point_3 p = output_mesh.point(v);
		_vertices.push_back(iGameVertex(p.x(), p.y(), p.z()));
	}

	std::cout << "Faces and their vertices: " << std::endl;
	std::vector<std::vector<iGameVertexHandle>>  _faces;//用于构建面的顶点坐标
	for (auto f : output_mesh.faces()) 
	{
		//std::cout << f << " consists of vertices: ";
		CGAL::Vertex_around_face_iterator<CGAL::Surface_mesh<Point_3>> vbegin, vend;
		std::vector<iGameVertexHandle> tempface;
		for (boost::tie(vbegin, vend) = vertices_around_face(output_mesh.halfedge(f), output_mesh); vbegin != vend; ++vbegin)
		{
			int vertex_index_int = get(CGAL::vertex_index, output_mesh, *vbegin);
			tempface.push_back(VH(vertex_index_int));
		}
		_faces.push_back(tempface);
		//std::cout << std::endl;
	}

	triMesh = new MeshKernel::SurfaceMesh(_vertices, _faces);
	std::cout << "Reconstruct success. #V = " << triMesh->vsize() << ", #F = " << triMesh->fsize() << std::endl;
	SetDrawModule(SURFACEMESH);
	updateDrawData();
	update();

}


void ModelDraw::loadShader_Possible_Sphere()
{
	shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	shaderProgram.setUniformValue("light.position", 1.0f, 1.0f, 0.0f);  // 光照方向
	//shaderProgram.setUniformValue("light.ambient", 0.2125f, 0.1275f, 0.054f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.714f, 0.4284f, 0.18144f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.393548f, 0.271906f, 0.166721f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 0.393548f);

	//金色
	shaderProgram.setUniformValue("light.ambient", 1.0f, 0.5f, 0.31f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 1.0f, 0.5f, 0.31f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 64.0f);

	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &VAO_Possible_Shpere);
	core->glDeleteBuffers(1, &VBO_Possible_Shpere);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &VAO_Possible_Shpere);
	core->glGenBuffers(1, &VBO_Possible_Shpere);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(VAO_Possible_Shpere);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, VBO_Possible_Shpere);
	core->glBufferData(GL_ARRAY_BUFFER, Data_Possible_Sphere.size() * sizeof(float), Data_Possible_Sphere.data(), GL_STATIC_DRAW);

	/*core->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	core->glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);*/

	/*float stride = (3 + 2 + 3) * sizeof(float);*/
	float stride = (3 + 3) * sizeof(float);
	core->glEnableVertexAttribArray(0);
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
	core->glEnableVertexAttribArray(1);
	core->glVertexAttribPointer(1, 2, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
	/*core->glEnableVertexAttribArray(1);
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));*/

	//这两行是解绑VAO 和VBO EBO
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ModelDraw::loadShader_Skel_Sphere(MeshKernel::iGameVertexHandle Node_h)
{
	//可能需要初始化一下VAO 和VBO EBO的大小

	 //相关数据传入sahder
	shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	shaderProgram.setUniformValue("light.position", 1.0f, 1.0f, 0.0f);  // 光照方向
	//shaderProgram.setUniformValue("light.ambient", 0.2125f, 0.1275f, 0.054f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.714f, 0.4284f, 0.18144f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.393548f, 0.271906f, 0.166721f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 0.393548f);

	//金色
	shaderProgram.setUniformValue("light.ambient", 1.0f, 0.5f, 0.31f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 1.0f, 0.5f, 0.31f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 64.0f);

	//shaderProgram.setUniformValue("light.ambient", 0.19225f, 0.19225f, 0.19225f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.50754f, 0.50754f, 0.50754f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.508273f, 0.508273f, 0.508273f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 32.0f);

	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &VAOs[Node_h]);
	core->glDeleteBuffers(1, &VBOs[Node_h]);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &VAOs[Node_h]);
	core->glGenBuffers(1, &VBOs[Node_h]);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(VAOs[Node_h]);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, VBOs[Node_h]);
	core->glBufferData(GL_ARRAY_BUFFER, data_Sphere[Node_h].size() * sizeof(float), data_Sphere[Node_h].data(), GL_STATIC_DRAW);

	/*core->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, EBO);
	core->glBufferData(GL_ELEMENT_ARRAY_BUFFER, indices.size() * sizeof(unsigned int), &indices[0], GL_STATIC_DRAW);*/

	/*float stride = (3 + 2 + 3) * sizeof(float);*/
	float stride = (3 + 3) * sizeof(float);
	core->glEnableVertexAttribArray(0);
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
	//传入法向信息
	core->glEnableVertexAttribArray(1);
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));
	/*core->glEnableVertexAttribArray(1);
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));*/
	

	//这两行是解绑VAO 和VBO EBO
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);

	//删除对应的VAO和VBO
	/*core->glDeleteVertexArrays(1, &VAOs[Node_h]);
	core->glDeleteBuffers(1, &VBOs[Node_h]);*/
	/*core->glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, 0);*/
}

void ModelDraw::SendToVAO_Skel()
{
	InitVaoSphereSize();// 和顶点个数有关
	InitIndices(); //绘制球体的下标顶点索引

	InitVaoCylinderSize();// 和边的个数有关
	

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	for (auto vmap : curve_skel_->allvertices())
	{
		//**
		updateCurrentSkelSphere(vmap.first); // 传入data
		//IntVertexHandleToSphereMesh.emplace(vmap.first, SphereMeshForSkeletonPoint.size() - 1); // Vetexhandle -->> IGameMesh Array index
		//SphereMeshToIntVertexHand.emplace(SphereMeshForSkeletonPoint.size() - 1, vmap.first);// IGameMesh Array index -->> Vetexhandle
		//**
		//loadShader_Skel_Sphere(vmap.first); // 把VAO索引好
		//DrawSkelDiscretize(vmap.first); //绘图放到DrawScen
	}


	for (auto emap : curve_skel_->alledges())
	{
		updateCurrentCylinderSkel(emap.first); // 传入data
	//	IntEdgeHandleToCylinderMesh.emplace(emap.first, CylinderMeshForSkeletonPoint.size() - 1); // Edgehandle -->> IGameMesh Array index
		//CylinderMeshToIntEdgeHandle.emplace(CylinderMeshForSkeletonPoint.size() - 1, emap.first);// IGameMesh Array index -->> Edgehandle
		// //**
		//loadShader_Skel_Cylinder(emap.first); // 把VAO索引好
		//DrawSkelDiscretize(emap.first); //绘图放到DrawScen
	}
	//printf("Send to VAO success!\n");
}

void ModelDraw::SendChangeToVAO_Skel()
{   
	//InitVaoSphereSize();// 和顶点个数有关
	//InitIndices(); //绘制球体的下标顶点索引

	//InitVaoCylinderSize();// 和边的个数有关

	// 现在实现点的移动模式下的改变VAO

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	auto& vh = selected_v_indices[0];

	
	updateCurrentSkelSphere((VH)vh);

	//auto eh :curve_skel_->NeighborEh((VH)vh);

	for (auto eh_change : curve_skel_->NeighborEh((VH)vh))
	{
		updateCurrentCylinderSkel(eh_change); // 传入data
		//loadShader_Skel_Cylinder(emap.first); // 把VAO索引好
		//DrawSkelDiscretize(emap.first); //绘图放到DrawScen
	}


	printf("Send Change to VAO success!\n");
}

void ModelDraw::InitAllPara()
{
	avg_radius = 0.0;
	max_radius = -1;
	min_radius = 9999;
	IsCalAvgRadius = false;
	local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::HIDE);
	selected_v_indices.clear();
	Move_Point_vh.clear();
	selected_e_indices.clear();
	selected_add_possible_point.clear();
	is_copy_point_draw == false;
	copy_point.reset();
	selected_add_possible_point.clear();
	selected_v_draw_indices.clear();
	MoveSelectedPoint.reset();
	is_change_radius = false;
	is_edit_point = false;
	//add_possible_point;
	//selected_add_possible_point.clear();
	copyTime_point = false;
	vh_copy_add.reset();
	e_add_h1.reset();

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;
	//计算平均半径 用于选点的距离
	
	avg_radius = 0;
	for (auto& vmap : curve_skel_->allvertices())
	{
		avg_radius += curve_skel_->vertex_radius_[vmap.first];
	}
	avg_radius /= curve_skel_->vertex_radius_.size();
	dotLineCalculation->AveRadiusCalculationTool = avg_radius; //传给计算工具 用于点选
	local_coordinater_axis->set_len(avg_radius * 3);

}

void ModelDraw::InitOriginMesh()
{
	SphereMesh_Origin = new MeshKernel::SurfaceMesh(Vertices_Origin, Faces_Origin);
	//std::cout << "";
}

void ModelDraw::CADContrast()
{
	std::unordered_map<VH, Vector3d> vh2color;
	//ErrorRender::get_Hausdorff_distance(*proDraw->tri_mesh_, *SphereMesh_MeshCADForContrast, vh2color);
	//ErrorRender::get_Hausdorff_distance(*proDraw->tri_mesh_, *SphereMesh_MeshCADForContrast, vh2color);
	//vector<CBSplineSurface> surface = bsplineSurfaces;
	loadIGESSurface(bsplineSurfaces, *proDraw->tri_mesh_);
}

void ModelDraw::DrawSkelDiscretize(MeshKernel::iGameVertexHandle Node_h)
{
	//shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	//shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	//shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	//shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	//shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	//shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
	//shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 32.0f);

	//shaderProgram.setUniformValue("colormode", 222);
	if (Node_h == copy_point || Node_h == MoveSelectedPoint)
	{
		shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
		shaderProgram.setUniformValue("light.ambient", 0.0215f, 0.1745f, 0.0215f);	// 环境光照系数
		shaderProgram.setUniformValue("light.diffuse", 0.07568f, 0.61424f, 0.07568f);	// 漫反射系数
		shaderProgram.setUniformValue("light.specular", 0.633f, 0.727811f, 0.633f);	// 镜面反射系数
		shaderProgram.setUniformValue("shininess", 32.0f);
	}
	/*else
	{
		shaderProgram.setUniformValue("colormode", 222);
	}*/
	shaderProgram.setUniformValue("colormode", 222);
	//shaderProgram.setUniformValue("color", QVector4D(0 / 255, 0, 1.0, 0.5));

	core->glBindVertexArray(VAOs[Node_h]);
	/*core->glDrawElements(GL_TRIANGLE_STRIP, indexCount, GL_UNSIGNED_INT, 0);*/
	//data_Sphere 有法向和顶点的信息
	//glDrawArrays(GL_TRIANGLES 是按顺序绘制三角形的，也就是说
	core->glDrawArrays(GL_TRIANGLES, 0, data_Sphere[Node_h].size() / 3 / 2);

	/*core->glDeleteVertexArrays(1, &VAOs[Node_h]);
	core->glDeleteBuffers(1, &VBOs[Node_h]);*/
}

void ModelDraw::DrawPossibleSphere()
{
	//shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	//shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	//shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	//shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	//shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	//shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
	//shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 32.0f);

	shaderProgram.setUniformValue("colormode", 222);
	glDepthRange(0.0, 0.1);
	//shaderProgram.setUniformValue("color", QVector4D(0, 0 , 1, 0.5));

	core->glBindVertexArray(VAO_Possible_Shpere);
	/*core->glDrawElements(GL_TRIANGLE_STRIP, indexCount, GL_UNSIGNED_INT, 0);*/
	//core->glDrawArrays(GL_TRIANGLES, 0, Circle_SEGMENTS * 2 * 3);
	core->glDrawArrays(GL_TRIANGLES, 0, Data_Possible_Sphere.size() / 6);
}

void ModelDraw::DrawSkelDiscretize_Cylinder(MeshKernel::iGameEdgeHandle Edge_h)
{
	//shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	//shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	//shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	//shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	//shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	//shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
	//shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 32.0f);

	shaderProgram.setUniformValue("colormode", 222);
	//shaderProgram.setUniformValue("color", QVector4D(0 / 255, 1, 0.0, 0.25));

	core->glBindVertexArray(VAOs_Cylinder[Edge_h]);
	/*core->glDrawElements(GL_TRIANGLE_STRIP, indexCount, GL_UNSIGNED_INT, 0);*/
	//core->glDrawArrays(GL_TRIANGLES, 0, Circle_SEGMENTS * 2 * 3);
	core->glDrawArrays(GL_TRIANGLES, 0, data_Cylinder[Edge_h].size() / 6);


	//shaderProgram.setUniformValue("colormode", 222);
	////shaderProgram.setUniformValue("color", QVector4D(0 / 255, 0, 1.0, 0.5));

	//core->glBindVertexArray(VAOs[Node_h]);
	///*core->glDrawElements(GL_TRIANGLE_STRIP, indexCount, GL_UNSIGNED_INT, 0);*/
	////data_Sphere 有法向和顶点的信息
	//core->glDrawArrays(GL_TRIANGLES, 0, data_Sphere[Node_h].size() / 3 / 2);

	/*core->glDeleteVertexArrays(1, &VAOs_Cylinder[Edge_h]);
	core->glDeleteBuffers(1, &VBOs_Cylinder[Edge_h]);*/
}

void ModelDraw::updateCurrentCylinderSkel(MeshKernel::iGameEdgeHandle Edge_h)
{

	positions_Cylinder_Zero.clear(); // 用于记录圆柱底面的分割点
	positions_Cylinder_One.clear();
	data_Cylinder[Edge_h].clear();

	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	auto& Edge = curve_skel_->edges(Edge_h);

	auto& Edge_vh1 = Edge.vh1();
	auto& Edge_vh2 = Edge.vh2();

	auto& Edge_v1 = curve_skel_->vertices(Edge_vh1);
	auto& Edge_v2 = curve_skel_->vertices(Edge_vh2);

	MeshKernel::iGameVertex dir = Edge_v2 - Edge_v1;

	QVector3D Edge_Vector(dir.x(), dir.y(), dir.z());
	QVector3D x_Vector(1, 0, 0);
	QVector3D y_Vector(0, 1, 0);
	QVector3D z_Vector(0, 0, 1);
	QVector3D zero(0, 0, 0);
	QVector3D dir_Circular = QVector3D::crossProduct(Edge_Vector, x_Vector);

	if (dir_Circular == zero)
	{
		dir_Circular = QVector3D::crossProduct(Edge_Vector, y_Vector);
		if (dir_Circular == zero)
		{
			dir_Circular = QVector3D::crossProduct(Edge_Vector, z_Vector);
		}
	}

	dir_Circular.normalize();
	//dir_Circular = dir_Circular * 20;

	double v1_radius = curve_skel_->vertex_radius_[Edge_vh1]; //用于移动到切平面的位置
	double v2_radius = curve_skel_->vertex_radius_[Edge_vh2];

	Edge_Vector.normalize();
	QVector3D Bias_Cylinder_One = Edge_Vector.normalized() * v1_radius; // v1->v2
	QVector3D Bias_Cylinder_Two = -(Edge_Vector.normalized()) * v2_radius; // v2->v1

	//QVector3D point_zero = dir_Circular * v1_radius + QVector3D(Edge_v1.x(), Edge_v1.y(), Edge_v1.z()) + Bias_Cylinder_One * 0.;
	//*******
	QVector3D point_zero = dir_Circular * v1_radius + QVector3D(Edge_v1.x(), Edge_v1.y(), Edge_v1.z()) + Bias_Cylinder_One;

	//QVector3D point_one = dir_Circular * v2_radius + QVector3D(Edge_v2.x(), Edge_v2.y(), Edge_v2.z()) + Bias_Cylinder_Two * 0.;
	//******
	QVector3D point_one = dir_Circular * v2_radius + QVector3D(Edge_v2.x(), Edge_v2.y(), Edge_v2.z()) + Bias_Cylinder_Two;


	for (int i = 0; i < Circle_SEGMENTS; i++)
	{
		positions_Cylinder_Zero.push_back(glm::vec3(point_zero.x(), point_zero.y(), point_zero.z()));
		positions_Cylinder_One.push_back(glm::vec3(point_one.x(), point_one.y(), point_one.z()));
		QMatrix4x4 rotatematrix;
		rotatematrix.setToIdentity();
		rotatematrix.translate(Edge_v1.x(), Edge_v1.y(), Edge_v1.z());
		rotatematrix.rotate(360 / Circle_SEGMENTS, Edge_Vector.normalized());
		rotatematrix.translate(-Edge_v1.x(), -Edge_v1.y(), -Edge_v1.z());
		point_zero = rotatematrix * point_zero; //旋转10pai 后再压入
		point_one = rotatematrix * point_one;
	}

	std::vector<iGameVertex> _vertices;
	//vertices : 圆心1 圆周长上的点1 圆心2 圆周长上的点2   1 + Circle_SEGMENTS + 1 + Circle_SEGMENTS
	//圆心1 ****为了重叠
	//_vertices.push_back(iGameVertex(Edge_v1.x() + Bias_Cylinder_One.x() * 0., Edge_v1.y() + Bias_Cylinder_One.y() * 0., Edge_v1.z() + Bias_Cylinder_One.z() * 0.));
	//*******
	_vertices.push_back(iGameVertex(Edge_v1.x() + Bias_Cylinder_One.x() , Edge_v1.y() + Bias_Cylinder_One.y(), Edge_v1.z() + Bias_Cylinder_One.z()));
	//圆周长上的点1
	for (int i = 0; i < positions_Cylinder_Zero.size(); i++)
	{
		iGameVertex tempVertex(positions_Cylinder_Zero[i].x, positions_Cylinder_Zero[i].y, positions_Cylinder_Zero[i].z);
		_vertices.push_back(tempVertex);
	}
	//圆心2
	//****
	_vertices.push_back(iGameVertex(Edge_v2.x() + Bias_Cylinder_Two.x() , Edge_v2.y() + Bias_Cylinder_Two.y() , Edge_v2.z() + Bias_Cylinder_Two.z()));
	//_vertices.push_back(iGameVertex(Edge_v2.x() + Bias_Cylinder_Two.x() * 0., Edge_v2.y() + Bias_Cylinder_Two.y() * 0., Edge_v2.z() + Bias_Cylinder_Two.z() * 0.));

	//圆周长上的点2
	for (int i = 0; i < positions_Cylinder_One.size(); i++)
	{
		iGameVertex tempVertex(positions_Cylinder_One[i].x, positions_Cylinder_One[i].y, positions_Cylinder_One[i].z);
		_vertices.push_back(tempVertex);
	}

	//PreNumHandle = _vertices.size() + PreNumHandle;

	std::vector<std::vector<iGameVertexHandle>> _faces;
	std::vector<iGameVertexHandle> tempface;
	// 用于生成_faces的索引

	for (int i = 0; i < Circle_SEGMENTS; i++)
		//压入圆的数据 朝外
	{
		if (i != Circle_SEGMENTS - 1)
		{
			tempface.clear();
			tempface.push_back(VH(i + 1 + 1));
			tempface.push_back(VH(i + 1 + 0)); // 下标为0存的是圆心
			//tempface.push_back(VH(i + 1 + 1));
			tempface.push_back(VH(0));
			_faces.push_back(tempface);

			tempface.clear();
			//tempface.push_back(VH(i + 1 + Circle_SEGMENTS + 1 + 1));
			tempface.push_back(VH(i + 1 + Circle_SEGMENTS + 1 + 0));
			tempface.push_back(VH(i + 1 + Circle_SEGMENTS + 1 + 1));
			tempface.push_back(VH(1 + Circle_SEGMENTS));
			_faces.push_back(tempface);

		}

		else
		{
			tempface.clear();
			tempface.push_back(VH(1));
			tempface.push_back(VH(i + 1 + 0)); // 下标为0存的是圆心
			//tempface.push_back(VH(1));
			tempface.push_back(VH(0));
			_faces.push_back(tempface);


			tempface.clear();
			//tempface.push_back(VH(1 + Circle_SEGMENTS + 1));
			tempface.push_back(VH(i + 1 + Circle_SEGMENTS + 1 + 0));
			tempface.push_back(VH(1 + Circle_SEGMENTS + 1));
			tempface.push_back(VH(1 + Circle_SEGMENTS));
			_faces.push_back(tempface);

		}
	}

	for (int i = 0; i < Circle_SEGMENTS; i++)
	{
		// 压入圆柱侧面的数据
		if (i != Circle_SEGMENTS - 1)
		{
			//朝外
			tempface.clear();
			tempface.push_back(VH(1 + i));
			tempface.push_back(VH(1 + i + 1));
			tempface.push_back(VH(1 + i + Circle_SEGMENTS + 2));
			_faces.push_back(tempface);

			//朝外
			tempface.clear();
			tempface.push_back(VH(1 + i));
			tempface.push_back(VH(1 + i + Circle_SEGMENTS + 2));
			tempface.push_back(VH(1 + i + Circle_SEGMENTS + 1));
			_faces.push_back(tempface);


		}
		else
		{
			//朝外
			tempface.clear();
			tempface.push_back(VH(1 + i));
			tempface.push_back(VH(1));
			tempface.push_back(VH(Circle_SEGMENTS + 2));
			_faces.push_back(tempface);

			tempface.clear();
			tempface.push_back(VH(1 + i));
			tempface.push_back(VH(Circle_SEGMENTS + 2));
			tempface.push_back(VH(1 + i + Circle_SEGMENTS + 1));
			_faces.push_back(tempface);


		}



	}

	CylinderMeshForSkeletonPoint.push_back(MeshKernel::SurfaceMesh(_vertices, _faces)); //构造的时候顶点也是按之前输入的顺序构造出来的

	if (!SphereMesh_Origin)
	{
		//用于构造原始三角形网格用于比较拟合误差
		for (int i = 0; i < _vertices.size(); i++)
		{
			Vertices_Origin.push_back(_vertices[i]);
		}

		for (int i = 0; i < _faces.size(); i++)
		{
			std::vector<iGameVertexHandle> TempFaceAddToOrigin;
			//每个面
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][0]) + PreNumHandle));
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][1]) + PreNumHandle));
			TempFaceAddToOrigin.push_back(VH(int(_faces[i][2]) + PreNumHandle));
			//PreNumHandle
			Faces_Origin.push_back(TempFaceAddToOrigin);
		}

		PreNumHandle = _vertices.size() + PreNumHandle;
	}


	//生成所有点的法向
	CylinderMeshForSkeletonPoint.back().genAllVerticesNormal();

	for (auto& fmap : CylinderMeshForSkeletonPoint.back().allfaces())
	{
		for (auto& vh : fmap.second.getVertexHandle())
		{
			//把每个面都压进去
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).x());
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).y());
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).z());

			//法向信息 顺时针传进去的
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).getNormalX());
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).getNormalY());
			data_Cylinder[Edge_h].push_back(CylinderMeshForSkeletonPoint.back().vertices(vh).getNormalZ());
		}


	}
	std::cout << "Clingder render is ok!" << std::endl;



	//传入点云数据用于重构出网格
	//修复:可能要进行初始化
	//********什么情况下使用这个函数
	//int IndexForSample = 0;
	//vector<Vector_3> Normal_Cylinder; //法向信息，用于点云的构建
	//for (auto vmap : SphereMesh.allvertices())
	//{

	//	if (IndexForSample == 0 || IndexForSample == (Circle_SEGMENTS + 1))
	//	{
	//		//不要两个圆的中心点
	//		IndexForSample++;
	//		continue;
	//	}

	//	auto vh = vmap.first;
	//	Point_3 first_point(SphereMesh.vertices(vh).x(), SphereMesh.vertices(vh).y(), SphereMesh.vertices(vh).z());
	//	Vector_3 first_normal(SphereMesh.vertices(vh).getNormalX(), SphereMesh.vertices(vh).getNormalY(), SphereMesh.vertices(vh).getNormalZ()); //法相的坐标
	//	CloudPointForPoissonReconstruction.insert(first_point, first_normal);

	//	//记录一下法向，采样点也可以用
	//	if (Normal_Cylinder.size() < Circle_SEGMENTS)
	//	{
	//		Normal_Cylinder.push_back(Vector_3(SphereMesh.vertices(vh).getNormalX(), SphereMesh.vertices(vh).getNormalY(), SphereMesh.vertices(vh).getNormalZ()));

	//	}
	//	//Normal_Cylinder.push_back(Vector_3(SphereMesh.vertices(vh).getNormalX(), SphereMesh.vertices(vh).getNormalY(), SphereMesh.vertices(vh).getNormalZ()));
	//	IndexForSample++;
	//}

	//int NumOfSample = 4;
	//double TotalNumOfSampleBreak = NumOfSample + 1; //int 会错 在计算分数的时候
	//for (int index = 0; index < NumOfSample; index++)
	//{
	//	for (int k = 0; k < Circle_SEGMENTS; k++)
	//	{
	//		Point_3 first_point(_vertices[1 + k].x() * ((index + 1) / TotalNumOfSampleBreak) + _vertices[1 + Circle_SEGMENTS + 1 + k].x() * ((TotalNumOfSampleBreak - index - 1) / TotalNumOfSampleBreak), _vertices[1 + k].y() * ((index + 1) / TotalNumOfSampleBreak) + _vertices[1 + Circle_SEGMENTS + 1 + k].y() * ((TotalNumOfSampleBreak - index - 1) / TotalNumOfSampleBreak), _vertices[1 + k].z() * ((index + 1) / TotalNumOfSampleBreak) + _vertices[1 + Circle_SEGMENTS + 1 + k].z() * ((TotalNumOfSampleBreak - index - 1) / TotalNumOfSampleBreak)); //点的坐标
	//		Vector_3 first_normal(Normal_Cylinder[k]); //法相的坐标
	//		CloudPointForPoissonReconstruction.insert(first_point, first_normal);
	//	
	//	}
	//}
	//std::cout << "Clingder Cloud Point render is ok!" << std::endl;

	
	  
}

void ModelDraw::loadShader_Skel_Cylinder(MeshKernel::iGameEdgeHandle Edge_h) 
{
	//可能需要初始化一下VAO 和VBO EBO的大小

	 //相关数据传入sahder

	shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向

	//翠绿
	//shaderProgram.setUniformValue("light.ambient", 0.0215f, 0.1745f, 0.0215f);	// 环境光照系数
	//shaderProgram.setUniformValue("light.diffuse", 0.07568f, 0.61424f, 0.07568f);	// 漫反射系数
	//shaderProgram.setUniformValue("light.specular", 0.633f, 0.727811f, 0.633f);	// 镜面反射系数
	//shaderProgram.setUniformValue("shininess", 64.0f);


	shaderProgram.setUniformValue("light.ambient", 0.19225f, 0.19225f, 0.19225f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 0.50754f, 0.50754f, 0.50754f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.508273f, 0.508273f, 0.508273f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 32.0f);


	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &VAOs_Cylinder[Edge_h]);
	core->glDeleteBuffers(1, &VBOs_Cylinder[Edge_h]);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &VAOs_Cylinder[Edge_h]);
	core->glGenBuffers(1, &VBOs_Cylinder[Edge_h]);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(VAOs_Cylinder[Edge_h]);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, VBOs_Cylinder[Edge_h]);
	core->glBufferData(GL_ARRAY_BUFFER, data_Cylinder[Edge_h].size() * sizeof(float), data_Cylinder[Edge_h].data(), GL_STATIC_DRAW);


	float stride = (3 + 3) * sizeof(float);
	//传入顶点信息
	core->glEnableVertexAttribArray(0);
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, stride, (void*)0);
	//传入法向信息
	core->glEnableVertexAttribArray(1);
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, stride, (void*)(3 * sizeof(float)));

	//这两行是解绑VAO 和VBO EBO
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);


}

void ModelDraw::drawSurfaceMesh_skel() {
	glLineWidth(2);
	shaderProgram.setUniformValue("colormode", 236);
	glBegin(GL_LINES);
	for (auto e_it = projectSkel->triMeshtv.edges_begin(); e_it != projectSkel->triMeshtv.edges_end(); ++e_it)
	{
		// 得到边所代表的半边
		OpenMesh::HalfedgeHandle heh1 = projectSkel->triMeshtv.halfedge_handle(*e_it, 0); // 默认一个方向的半边
		OpenMesh::Vec3d edgePoint(0, 0, 0);
		int edgePointsNumber = 0;
		OpenMesh::DefaultTraits::Point pointV = projectSkel->triMeshtv.point(projectSkel->triMeshtv.from_vertex_handle(heh1)); // 这条（半）边的起点
		OpenMesh::DefaultTraits::Point pointW = projectSkel->triMeshtv.point(projectSkel->triMeshtv.to_vertex_handle(heh1));   // 这条（半）边的终点
		glVertex3f(pointV[0], pointV[1], pointV[2]);
		glVertex3f(pointW[0], pointW[1], pointW[2]);
	}
	glEnd();
}

void ModelDraw::updateSurfaceMesh_skel() {
	//wait
	for (auto e_it = projectSkel->triMeshtv.vertices_begin(); e_it != projectSkel->triMeshtv.vertices_end(); ++e_it)
	{
		
	}
}

void ModelDraw::DrawStyle(void) {

	if (this->isEnableLighting) glEnable(GL_LIGHTING);
	else glDisable(GL_LIGHTING);

	switch (drawmode) {
	case POINTS:
		core->glPointSize(4);
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		shaderProgram.setUniformValue("rendermode", 1);
		break;
	case FLATLINES:
		core->glLineWidth(1);// 线宽 
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		shaderProgram.setUniformValue("rendermode", 2);
		break;
	case FILL:
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		shaderProgram.setUniformValue("rendermode", 3);
		break;
	default:
		break;
	}
}

void ModelDraw::ColorStyle(void) {
	switch (colormode) {
	case LIGHT:
		shaderProgram.setUniformValue("colormode", 1);
		break;
	case NORMAL:
		shaderProgram.setUniformValue("colormode", 2);
		break;
	default:
		break;
	}
}

void ModelDraw::pickPoint(double x, double y) {

}

void ModelDraw::pickPoint(int x, int y, vector<double>& p) {

}

void ModelDraw::pickPoint_getDepth(int x, int y, vector<double>& p, double& z) {

}

void ModelDraw::pickPointDepth(int x, int y, vector<double>& p, double& d) {// 用深度信息来 pick 点

}

double ModelDraw::point2LineDis(QVector3D point1, QVector3D point2, QVector3D point0) 
{

	//直线的法向量(p,q,r)
	QVector3D normal = (point1 - point2).normalized();

	//两个向量
	QVector3D normal01 = (point1 - point0).normalized();
	QVector3D normal02 = (point2 - point0).normalized();

	//进行叉乘
	double cross_x = normal01.y() * normal02.z() - normal02.y() * normal01.z();
	double cross_y = normal02.x() * normal01.z() - normal02.z() * normal01.x();
	double cross_z = normal01.x() * normal02.y() - normal01.y() * normal02.x();

	double cross_length = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
	double dx = (point2 - point1).x();
	double dy = (point2 - point1).y();
	double dz = (point2 - point1).z();

	double qb_length = sqrt(dx * dx + dy * dy + dz * dz);

	double ds = cross_length / qb_length;
	return ds;

}

bool ModelDraw::isIntersectTriangle(QVector3D orig, QVector3D end, QVector3D v0, QVector3D v1, QVector3D v2, QVector3D& intersection) {
	// Refer to https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html
	float t, u, v;
	QVector3D dir = end - orig;
	QVector3D E1 = v1 - v0;
	QVector3D E2 = v2 - v0;
	// 叉乘
	QVector3D P(dir.y() * E2.z() - E2.y() * dir.z(),
		E2.x() * dir.z() - E2.z() * dir.x(),
		dir.x() * E2.y() - dir.y() * E2.x());
	// 点乘
	float det = E1.x() * P.x() + E1.y() * P.y() + E1.z() * P.z();
	QVector3D T;
	if (det > 0) {
		T = orig - v0;
	}
	else {
		T = v0 - orig;
		det = -det;
	}
	if (det < 0)
		return false;

	// Calculate u and make sure u <= 1
	u = T.x() * P.x() + T.y() * P.y() + T.z() * P.z();
	if (u < 0.0f || u > det)
		return false;

	// Q
	QVector3D Q(T.y() * E1.z() - E1.y() * T.z(),
		E1.x() * T.z() - E1.z() * T.x(),
		T.x() * E1.y() - T.y() * E1.x());

	// Calculate v and make sure u + v <= 1
	v = dir.x() * Q.x() + dir.y() * Q.y() + dir.z() * Q.z();
	if (v < 0.0f || u + v > det)
		return false;

	//// Calculate t, scale parameters, ray intersects triangle
	//t = E2.x() * Q.x() + E2.y() * Q.y() + E2.z() * Q.z();
	float fInvDet = 1.0f / det;
	//t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;
	intersection = u * E1 + v * E2 + v0;
	return true;
}

void ModelDraw::removeSelectedOne() {
	if (!selected_v_indices.empty() && pickType == PickType::PickPoint) selected_v_indices.resize(selected_v_indices.size() - 1);
	if (!selected_e_indices.empty() && pickType == PickType::PickEdge) selected_e_indices.resize(selected_e_indices.size() - 1);
	if (!selected_f_indices.empty() && pickType == PickType::PickFace) selected_f_indices.resize(selected_f_indices.size() - 1);
	if (!selected_c_indices.empty() && pickType == PickType::PickCell) selected_c_indices.resize(selected_c_indices.size() - 1);
}

void ModelDraw::removeSelectedAll() {
	this->selected_v_indices.clear();
	this->selected_e_indices.clear();
	this->selected_f_indices.clear();
	this->selected_c_indices.clear();
	vector<int>().swap(this->selected_v_indices);
	vector<int>().swap(this->selected_e_indices);
	vector<int>().swap(this->selected_f_indices);
	vector<int>().swap(this->selected_c_indices);
}

void ModelDraw::deleteSelectedAll() {
	switch (pickType) {
	case PickPoint: {
		for (auto& i : selected_v_indices) {
			MeshKernel::iGameVertexHandle vh(i);
			if (hexMesh) hexMesh->DeleteVertex(vh);
			else if (triMesh) triMesh->DeleteVertex(vh);
		}
		selected_v_indices.clear();
		break;
	}
	case PickEdge: {
		for (auto& i : selected_e_indices) {
			MeshKernel::iGameEdgeHandle eh(i);
			if (hexMesh) hexMesh->DeleteEdge(eh);
			else if (triMesh) triMesh->DeleteEdge(eh);
		}
		selected_e_indices.clear();
		break;
	}
	case PickFace: {
		for (auto& i : selected_f_indices) {
			MeshKernel::iGameFaceHandle fh(i);
			if (hexMesh) hexMesh->DeleteFace(fh);
			else if (triMesh) triMesh->DeleteFace(fh);
		}
		selected_f_indices.clear();
		break;
	}
	case PickCell: {
		for (auto& i : selected_c_indices) {
			MeshKernel::iGameCellHandle ch(i);
			hexMesh->DeleteCell(ch);
		}
		selected_c_indices.clear();
		break;
	}
	}
	printf("delete selected all success\n");
	updateAllHandles();
	updateDrawData();
	update();
}

void ModelDraw::add_new_cell() {
	if (!hexMesh) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Support Hex Mesh Only! "));
		return;
	}
	if (selected_f_indices.size() == 2) {
		FH fh1(selected_f_indices[0]), fh2(selected_f_indices[1]);
		add_cell_2fh(*hexMesh, fh1, fh2);
	} else if (selected_e_indices.size() == 1 && selected_f_indices.size() == 1) {
		EH eh(selected_e_indices[0]);
		FH fh(selected_f_indices[0]);
		add_cell_eh_fh(*hexMesh, eh, fh);
	} else if (selected_f_indices.size() == 3) {
		FH fh1(selected_f_indices[0]), fh2(selected_f_indices[1]), fh3(selected_f_indices[2]);
		add_cell_3fh(*hexMesh, fh1, fh2, fh3);
	} else if (selected_f_indices.size() == 1) {
		FH fh(selected_f_indices[0]);
		add_cell_fh(*hexMesh, fh);
	}
	else if (selected_c_indices.size() == 2)
	{
		CH ch1(selected_c_indices[0]);
		CH ch2(selected_c_indices[1]);
		Add_Cell_Neighbour_For_MMC_Bifurcate(*hexMesh, ch1, ch2);
		hexMesh->DeleteCell(ch1);
		hexMesh->DeleteCell(ch2);
		selected_c_indices.clear();
	}
	updateAllHandles();
	updateDrawData();
}

void ModelDraw::initializePick()
{
	selected_c_indices.clear();
	solutionSkel->selectNodes.clear();
	solutionSkel->selectNode = nullptr;
}

void ModelDraw::pickEvent(double winX, double winY) 
{
	switch (drawmodule)
	{
	case ViewOpenGL::EMPTY:
		break;
	case ViewOpenGL::SURFACEMESH: //这个会在读入网格的时候 设置类型
		pickEvent_TriMesh_BVH(winX, winY);
		showPickInfo();
		break;
	case ViewOpenGL::QUADMESH:
		pickEvent_QuadMesh_BVH(winX, winY); // 仿照这个写
		break;
	case ViewOpenGL::TETMESH:
		pickEvent_TetMesh_BVH(winX, winY);
		break;
	case ViewOpenGL::HEXMESH:
		pickEvent_HexMesh_BVH(winX, winY); // 压入到selected_v_indices
		showPickInfo();
		break;
	case ViewOpenGL::SKEL:
		pickEvent_SKEL(winX, winY);
		showPickInfo();
		break;
	default:
		break;
	}
}







void ModelDraw::pickEvent_TriMesh(double winX, double winY)
{
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	double zMax = -numeric_limits<double>::max();

	std::vector<QVector3D> tempPoints;
	MeshKernel::iGameFaceHandle tempFaceHandle(-1);
	MeshKernel::iGameFace tempFace;
	QVector3D tempIntersection;

	for (auto& fp : triMesh->allfaces()) 
	{
		MeshKernel::iGameFaceHandle fh = fp.first;
		auto& f = fp.second;
		auto vhs = f.getVertexHandle();
		vector<QVector3D> face_vertices;
		QVector4D face_depth(0, 0, 0, 0);

		for (auto& vh : vhs) 
		{
			auto& v = triMesh->vertices(vh);
			auto qv = QVector3D(v.x(), v.y(), v.z());
			face_vertices.push_back(qv);
			face_depth += QVector4D(qv, 1);
		}

		face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
		face_depth = modelmatrix * face_depth;
		bool intersection_flag = false;
		QVector3D intersection_point;

		if (face_vertices.size() == 3) 
		{
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
		}
		else if (face_vertices.size() == 4) 
		{
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
				| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
		}
		if (intersection_flag && zMax < face_depth[2]) 
		{
			// 记录了相交的面 和 六面体相交的点
			zMax = face_depth[2];
			tempFaceHandle = fh;
			tempFace = f;
			tempIntersection = intersection_point;
			tempPoints = face_vertices;
		}
	}



	if (tempFaceHandle != -1) 
	{
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		MeshKernel::iGameVertex pos(-1, -1, -1);

		switch (pickType)
		{

		case PickType::PickPoint:
			for (auto& vh : tempFace.getVertexHandle())
			{
				auto& v = triMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				// point1 和 point2 是鼠标的线段点
				double dist = point2LineDis(point1, point2, point0); 
				if (dist < min_distance)
				{
					// 选取面上的点
					min_distance = dist;
					nearest_handle = vh;
				}
			}

			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) 
			{
				selected_v_indices.push_back(nearest_handle);
			}

			break;

		case PickType::PickEdge:

			for (auto& eh : tempFace.getEdgeHandle()) 
			{
				auto& e = triMesh->edges(eh);
				auto& v1 = triMesh->vertices(e.vh1());
				auto& v2 = triMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection); // 面上边 和 交点 的距离
				if (dist < min_distance) 
				{
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end())
			{
				selected_e_indices.push_back(nearest_handle);
			}

			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) 
			{
				selected_f_indices.push_back(nearest_handle);
			}

			
			break;
		default:
			break;
		}
	}
	else
	{

	}


}

bool ModelDraw::pickEvent_SKEL(double winX, double winY) 
{
	// winX winY 是鼠标 2D位置
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);		// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	double zMax = -numeric_limits<double>::max();
	double item[3]{ 0, 0, 0 };

	std::vector<QVector3D> selectPoints;
	MeshKernel::iGameCell selectCell;
	MeshKernel::iGameCellHandle selectCellHandle(-1);
	MeshKernel::iGameFace selectFace;
	MeshKernel::iGameFaceHandle selectFaceHandle(-1);

	QVector3D selectIntersection;

	for (auto& cp : hexMesh->allcells()) 
	{		// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		MeshKernel::iGameCell c = cp.second;
		for (auto& fh : c.getFaceHandle())
		{
			auto& f = hexMesh->faces(fh);
			auto vhs = f.getVertexHandle();
			vector<QVector3D> face_vertices;
			QVector4D face_depth(0, 0, 0, 0);
			for (auto& vh : vhs) 
			{
				auto& v = hexMesh->vertices(vh);
				auto qv = QVector3D(v.x(), v.y(), v.z());
				face_vertices.push_back(qv);
				face_depth += QVector4D(qv, 1);
			}
			face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
			face_depth = modelmatrix * face_depth;
			bool intersection_flag = false;
			QVector3D intersection_point;
			if (face_vertices.size() == 3) 
			{
				intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
			}
			else if (face_vertices.size() == 4) 
			{
				intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
					| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
			}
			if (intersection_flag && zMax < face_depth[2]) 
			{
				zMax = face_depth[2];
				selectCellHandle = cp.first;
				selectCell = cp.second;
				selectFaceHandle = fh;
				selectFace = f;
				selectPoints = face_vertices;
				selectIntersection = intersection_point;
			}
		}
	}
	
	// 如果有选中的体 就返回true
	if (selectCellHandle != -1) 
	{
		MeshKernel::iGameVertex pos;
		switch (pickType) { // 某个按键会触发类型
		case PickType::PickPoint:
		//	
		//	for (auto& vh : tempFace.getVertexHandle()) {
		//		auto& v = hexMesh->vertices(vh);
		//		QVector3D point0(v.x(), v.y(), v.z());
		//		double dist = point2LineDis(point1, point2, point0);
		//		if (dist < min_distance) {
		//			min_distance = dist;
		//			nearest_handle = vh;
		//		}
		//	}
		//	if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
		//		selected_v_indices.push_back(nearest_handle);
		//	}
		//	//printf("the selected vertex handle is %d\n", nearest_handle);
		//	pos = hexMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
		//	this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
		//		std::to_string(pos.z()) + ")";
		//	std::cout << pick_info << std::endl;
		//	break;
			
		case PickType::PickEdge:
			/*
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = hexMesh->edges(eh);
				auto& v1 = hexMesh->vertices(e.vh1());
				auto& v2 = hexMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
			printf("the selected edge handle is %d\n", nearest_handle);
			break;
			*/
		case PickType::PickFace:
			/*
			nearest_handle = (int)tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
			printf("the selected face handle is %d\n", nearest_handle);
			break;
			*/

		case PickType::PickCell:
			px = selectIntersection[0];
			py = selectIntersection[1];
			pz = selectIntersection[2];
			if ((selected_c_indices.size() == 1 && selected_c_indices[0] != selectCellHandle)
				|| selected_c_indices.size() == 0) {
				selected_c_indices.clear();
				selected_c_indices.push_back(selectCellHandle);
				pos = hexMesh->getCellCenter(selectCellHandle);
				// check with all node position
				solutionSkel->selectNodes.clear();
				BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(pos[0], pos[1], pos[2]));
				SingalFather::getInstance()->getData(rlt->getPosition().x, rlt->getPosition().y, rlt->getPosition().z,
					rlt->root->getRotateX(), rlt->root->getRotateY(), rlt->root->getRotateZ(),
					rlt->root->getLength(), rlt->root->getWidth(), rlt->root->getHeight());
				
				// 坐标系设置
				local_coordinater_axis->set_orign(pos.x(), pos.y(), pos.z());
				local_coordinater_axis->set_mode_type(select_coordinater_mode);
				printf("the selected cell mid is %lf %lf %lf\n", pos.x(), pos.y(), pos.z());
			}
			printf("the selected cell handle is %d\n", selectCellHandle);

			break;

		/*
		case PickType::PickCells:
			for (auto& cp : hexMesh->allcells()) {
				if (hexMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			px = tempIntersection[0];
			py = tempIntersection[1];
			pz = tempIntersection[2];

			if ((it = std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle)) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
			else {
				selected_c_indices.erase(it);
			}
			// update node
			solutionSkel->selectNodes.clear();
			for (int index = 0; index < selected_c_indices.size(); index++) {
				iGameVertex tmp(0.0, 0.0, 0.0);
				auto& cell = hexMesh->cells((iGameCellHandle)index);
				auto vhs = cell.getVertexHandle();
				for (auto& vh : vhs) {
					tmp += hexMesh->vertices(vh);
				}
				tmp = tmp / 8.0;
				// check with all node position
				BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(tmp[0], tmp[1], tmp[2]));
				if (rlt) {
					solutionSkel->selectNodes.push_back(rlt);
				}
			}

			printf("the selected cell handle is %d\n", nearest_handle);
			break;
		
		*/
		default:
			break;
		}
		return true;
	}

	else 
{
		printf("no intersection!!!\n");
		return false;
	}
}

void ModelDraw::pickEvent_TriMesh_BVH(double winX, double winY) 
{  // 鼠标点 的射线 和三角形网格的交点
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	auto isec = bvh_tree->getIntersection(ray);
	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index); //得到鼠标点击的那个面
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);


	if (tempFaceHandle != -1)
	{
		double min_distance = numeric_limits<double>::max();
		auto& tempFace = triMesh->faces(tempFaceHandle);
		int nearest_handle(-1);
		MeshKernel::iGameVertex pos(-1, -1, -1);
		switch (pickType) 
		{

			case PickType::PickPoint:
				for (auto& vh : tempFace.getVertexHandle()) 
				{
					auto& v = triMesh->vertices(vh);
					QVector3D point0(v.x(), v.y(), v.z());
					double dist = point2LineDis(point1, point2, point0);
					if (dist < min_distance) {
						min_distance = dist;
						nearest_handle = vh;
					}
				}
				if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) 
				{
					selected_v_indices.push_back(nearest_handle);
				}
#ifdef DEBUG
				pos = triMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
				//printf("the selected vertex handle is %d: (%.2f, %.2f, %.2f)\n", nearest_handle, pos.x(), pos.y(), pos.z());
				this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
					std::to_string(pos.z()) + ")";
				std::cout << pick_info << std::endl;
#endif // DEBUG
				break;

			case PickType::PickEdge:
				for (auto& eh : tempFace.getEdgeHandle())
				{
					auto& e = triMesh->edges(eh);
					auto& v1 = triMesh->vertices(e.vh1());
					auto& v2 = triMesh->vertices(e.vh2());
					QVector3D ev1(v1.x(), v1.y(), v1.z());
					QVector3D ev2(v2.x(), v2.y(), v2.z());
					double dist = point2LineDis(ev1, ev2, tempIntersection);
					if (dist < min_distance) 
					{
						min_distance = dist;
						nearest_handle = eh;
					}
				}
				if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) 
				{
					selected_e_indices.push_back(nearest_handle);
				}
#ifdef DEBUG
				printf("the selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
				break;
			case PickType::PickFace:
				nearest_handle = tempFaceHandle;
				if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
					selected_f_indices.push_back(nearest_handle);
				}
#ifdef DEBUG
				printf("the selected face handle is %d\n", nearest_handle);
#endif // DEBUG


				break;
			default:
				break;
		}
	} 	else {
#ifdef DEBUG
		printf("no intersection!!!\n");
#endif // DEBUG
	}


}

void ModelDraw::pick_point_with_depth(float x, float y, double(&p)[3], double itemPosition[3]) 
{
	//if (!this->triMesh && !this->tetMesh) return;
	QVector4D point(x, y, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(x, y, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointView(0, 0, 1, 1.0);
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector4D tpointView = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointView;
	// point1为射线起点，point2为射线终点
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	QVector3D faceNormal(tpointView.x() / tpointView.w(), tpointView.y() / tpointView.w(), tpointView.z() / tpointView.w());

	// 射线向量
	double e[3]{ point2.x() - point1.x(),point2.y() - point1.y(), point2.z() - point1.z() };
	// 射线起点
	double PL[3]{ point1.x(),point1.y(), point1.z() };
	// 平面向量
	double q[3]{ itemPosition[0] - faceNormal.x(),itemPosition[1] - faceNormal.y(), itemPosition[2] - faceNormal.z() };

	double t = (q[0] * (itemPosition[0] - PL[0]) + q[1] * (itemPosition[1] - PL[1]) + q[2] * (itemPosition[2] - PL[2])) / (q[0] * e[0] + q[1] * e[1] + q[2] * e[2]);
	p[0] = PL[0] + e[0] * t;
	p[1] = PL[1] + e[1] * t;
	p[2] = PL[2] + e[2] * t;
}

void ModelDraw::pickEvent_HexMesh(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	double zMax = -numeric_limits<double>::max();
	std::vector<QVector3D> tempPoints;
	MeshKernel::iGameFaceHandle tempFaceHandle(-1);
	MeshKernel::iGameFace tempFace;
	QVector3D tempIntersection;

	for (auto& fp : hexMesh->allfaces()) {// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(fp.first)) continue;
		MeshKernel::iGameFaceHandle fh = fp.first;
		auto& f = fp.second;
		auto vhs = f.getVertexHandle();
		vector<QVector3D> face_vertices;
		QVector4D face_depth(0, 0, 0, 0);
		for (auto& vh : vhs) {
			auto& v = hexMesh->vertices(vh);
			auto qv = QVector3D(v.x(), v.y(), v.z());
			face_vertices.push_back(qv);
			face_depth += QVector4D(qv, 1);
		}
		face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
		face_depth = modelmatrix * face_depth;
		bool intersection_flag = false;
		QVector3D intersection_point;
		if (face_vertices.size() == 3) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
		}
		else if (face_vertices.size() == 4) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
				| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
		}
		if (intersection_flag && zMax < face_depth[2]) {
			zMax = face_depth[2];
			tempFaceHandle = fh;
			tempFace = f;
			tempIntersection = intersection_point;
			tempPoints = face_vertices;
		}
	}
	
	/*Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	int time_bvh = clock() - time2;

	std::cout << "\nOld-Method fh = " << tempFaceHandle << ", cost: " << time_old << "ms. " << 
		" BVH fh = " << isec.index << ", cost: " << time_bvh << "ms.\n";
	printf("Old-IsecPos = (%.4f, %.4f, %.4f), BVH IsecPos = (%.4f, %.4f, %.4f)\n", tempIntersection.x(), tempIntersection.y(), tempIntersection.z(), 
		isec.pos[0], isec.pos[1], isec.pos[2]);*/

	if (tempFaceHandle != -1) {

		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = hexMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) 
		{
		case PickType::PickPoint:

			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = hexMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end())
			{
				selected_v_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected vertex handle is %d\n", nearest_handle);
			pos = hexMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
#endif // DEBUG

			
			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = hexMesh->edges(eh);
				auto& v1 = hexMesh->vertices(e.vh1());
				auto& v2 = hexMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
				// added for Top Opt
				/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
				pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
			}
#ifdef DEBUG
			printf("The selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				std::cout << "select: face handle is " << nearest_handle << std::endl;
				selected_f_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected face handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickCell:
			for (auto& cp : hexMesh->allcells()) {
				if (hexMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected cell handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		default:
			break;
		}
	}
	else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}

void ModelDraw::rupdateGL() {
	if (selected_c_indices.size() == 1) {
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[0]);
		auto vhs = cell.getVertexHandle();
		for (int i = 0; i < 8; i++) {
			iGameVertex t(solutionSkel->selectNode->root->p[i].x, solutionSkel->selectNode->root->p[i].y,
				solutionSkel->selectNode->root->p[i].z);
			auto& v = hexMesh->vertices(vhs[i]);
			v = t;
		}
	}
	SetDrawModule(SKEL);
	//updateMeshModel();
	updateDrawData();
	return;
}

void ModelDraw::drawNormalSkel() {
	ModelClear();
	hexMesh = new MeshKernel::VolumeMesh();
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) 
	{
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	//updateMeshModel();
	updateDrawData();
	update();
	return;
}

void ModelDraw::drawNormalSkel_1() {
	ModelClear();
	hexMesh = new MeshKernel::VolumeMesh();
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	updateMeshModel();
	updateDrawData();
	update();
	return;
}

void ModelDraw::pickEvent_HexMesh_BVH(double winX, double winY) 
{
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	int time_bvh = clock() - time2;

	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);

	if (tempFaceHandle != -1) 
	{
		auto& tempFace = hexMesh->faces(tempFaceHandle);
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = hexMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) 
		{
			case PickType::PickPoint:

				for (auto& vh : tempFace.getVertexHandle()) 
				{
					auto& v = hexMesh->vertices(vh);
					QVector3D point0(v.x(), v.y(), v.z());
					double dist = point2LineDis(point1, point2, point0);
					if (dist < min_distance) {
						min_distance = dist;
						nearest_handle = vh;
					}
				}
				if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) 
				{
					std::cout << "select: point handle is " << nearest_handle << std::endl;
					selected_v_indices.push_back(nearest_handle);
				}
				break;

			case PickType::PickEdge:
				for (auto& eh : tempFace.getEdgeHandle()) {
					auto& e = hexMesh->edges(eh);
					auto& v1 = hexMesh->vertices(e.vh1());
					auto& v2 = hexMesh->vertices(e.vh2());
					QVector3D ev1(v1.x(), v1.y(), v1.z());
					QVector3D ev2(v2.x(), v2.y(), v2.z());
					double dist = point2LineDis(ev1, ev2, tempIntersection);
					if (dist < min_distance) {
						min_distance = dist;
						nearest_handle = eh;
					}
				}
				if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
					std::cout << "select: edge handle is " << nearest_handle << std::endl;
					selected_e_indices.push_back(nearest_handle);
					// added for Top Opt
					/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
					pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
				}
				break;

			case PickType::PickFace:
				nearest_handle = tempFaceHandle;
				if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
					std::cout << "select: face handle is " << nearest_handle << std::endl;
					selected_f_indices.push_back(nearest_handle);
				}
				break;

			case PickType::PickCell:
				for (auto& cp : hexMesh->allcells()) {
					if (hexMesh->isOnBoundary(cp.first)) {
						auto cell_fhs = cp.second.getFaceHandle();
						if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
							nearest_handle = cp.first;
							break;
						}
					}
				}
				if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
					std::cout << "select: cell handle is " << nearest_handle << std::endl;
					selected_c_indices.push_back(nearest_handle);
				}
				break;

			default:
				break;
		}
	}

	else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}

void ModelDraw::pickEvent_TetMesh_BVH(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	int time_bvh = clock() - time2;

	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);

	if (tempFaceHandle != -1) {
		auto& tempFace = tetMesh->faces(tempFaceHandle);
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = tetMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) {
			case PickType::PickPoint:

				for (auto& vh : tempFace.getVertexHandle()) {
					auto& v = tetMesh->vertices(vh);
					QVector3D point0(v.x(), v.y(), v.z());
					double dist = point2LineDis(point1, point2, point0);
					if (dist < min_distance) {
						min_distance = dist;
						nearest_handle = vh;
					}
				}
				if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
					selected_v_indices.push_back(nearest_handle);
				}
				break;

			case PickType::PickEdge:
				for (auto& eh : tempFace.getEdgeHandle()) {
					auto& e = tetMesh->edges(eh);
					auto& v1 = tetMesh->vertices(e.vh1());
					auto& v2 = tetMesh->vertices(e.vh2());
					QVector3D ev1(v1.x(), v1.y(), v1.z());
					QVector3D ev2(v2.x(), v2.y(), v2.z());
					double dist = point2LineDis(ev1, ev2, tempIntersection);
					if (dist < min_distance) {
						min_distance = dist;
						nearest_handle = eh;
					}
				}
				if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
					selected_e_indices.push_back(nearest_handle);
					// added for Top Opt
					/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
					pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
				}
				break;

			case PickType::PickFace:
				nearest_handle = tempFaceHandle;
				if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
					selected_f_indices.push_back(nearest_handle);
				}
				break;

			case PickType::PickCell:
				for (auto& cp : tetMesh->allcells()) {
					if (tetMesh->isOnBoundary(cp.first)) {
						auto cell_fhs = cp.second.getFaceHandle();
						if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
							nearest_handle = cp.first;
							break;
						}
					}
				}
				if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
					selected_c_indices.push_back(nearest_handle);
				}
				break;

			default:
				break;
		}
	} else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}

void ModelDraw::pickEvent_QuadMesh_BVH(double winX, double winY) 
{
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	int time_bvh = clock() - time2;

	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);

	if (tempFaceHandle != -1) 
	{
		auto& tempFace = quadMesh->faces(tempFaceHandle);
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		MeshKernel::iGameVertex pos;

		switch (pickType) 
		{
			case PickType::PickPoint:

				for (auto& vh : tempFace.getVertexHandle())
				{
					auto& v = quadMesh->vertices(vh);
					QVector3D point0(v.x(), v.y(), v.z());
					double dist = point2LineDis(point1, point2, point0);
					if (dist < min_distance) 
					{
						min_distance = dist;
						nearest_handle = vh;
					}
				}
				if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) 
				{
					selected_v_indices.push_back(nearest_handle);
				}
				break;

			case PickType::PickEdge:
				for (auto& eh : tempFace.getEdgeHandle()) 
				{
					auto& e = quadMesh->edges(eh);
					auto& v1 = quadMesh->vertices(e.vh1());
					auto& v2 = quadMesh->vertices(e.vh2());
					QVector3D ev1(v1.x(), v1.y(), v1.z());
					QVector3D ev2(v2.x(), v2.y(), v2.z());
					double dist = point2LineDis(ev1, ev2, tempIntersection);
					if (dist < min_distance) 
					{
						min_distance = dist;
						nearest_handle = eh;
					}
				}
				if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) 
				{
					selected_e_indices.push_back(nearest_handle);
					// added for Top Opt
					/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
					pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
				}
				break;

			case PickType::PickFace:
				nearest_handle = tempFaceHandle;
				if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) 
				{
					std::cout << "select: face handle is " << nearest_handle << std::endl;
					selected_f_indices.push_back(nearest_handle);
				}
				break;

			default:
				break;
		}
	} else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}

// 槽函数只要写好就行了 按线程来呼应
void ModelDraw::Receive(QString text)
{
	//获取骨架
	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	//根据控件来修改半径
	curve_skel_->vertex_radius_[VH(selected_v_indices[0])] = text.toDouble();
	//update();
	updateCurrentSkelSphere(VH(selected_v_indices[0]));
	
	for (auto eh : curve_skel_->NeighborEh(VH(selected_v_indices[0])))
	{
		updateCurrentCylinderSkel(eh); // 传入data
		//loadShader_Skel_Cylinder(emap.first); // 把VAO索引好
		//DrawSkelDiscretize(emap.first); //绘图放到DrawScen
	}
}

void ModelDraw::Receive()
{
	MeshKernel::CurveSkel* curve_skel_ = nullptr;
	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
	if (!curve_skel_) return;

	// curve_skel_->vertex_radius_[VH(selected_v_indices[0])] = ui->RadiusChange->value();
	
	updateCurrentSkelSphere(VH(selected_v_indices[0]));

	for (auto eh : curve_skel_->NeighborEh(VH(selected_v_indices[0])))
	{
		updateCurrentCylinderSkel(eh); // 传入data
		//loadShader_Skel_Cylinder(emap.first); // 把VAO索引好
		//DrawSkelDiscretize(emap.first); //绘图放到DrawScen
	}
}


/****************************鼠标绘制事件******************************************/
void ModelDraw::mousePressEvent(QMouseEvent* _event) {
	
	/* 【通用】
	* 左键： 设定初始移动点供移动事件使用
	* 右键： 设定初始移动点供移动事件使用
	*/

	Qt::MouseButtons btns = _event->buttons();
	lastpoint2 = _event->pos();
	lastpointok = MapToSphere(lastpoint2, lastpoint3);
	mousemode = _event->button();						// 可以不用这个实现(可改)
	// 我要做的是 Skel_Move_Point
	switch (editmodel)  
	{

		case PICKITEM:  // 启动点选那个按钮 启动
		{// 为框选记录鼠标点击的点
			if (_event->button() == Qt::LeftButton)
			{
				lastpoint2 = _event->pos();
				double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
				pickEvent(winX, winY);
			}
			break;
		}

		case POINT_COORDINATE: 
		{
			/*
			* 左键： 1.pickEvent  2.坐标轴功能 XX 生成
			* 右键：
			*/
			if (btns == Qt::LeftButton) 
			{
				// 没点击中CUBE与AXIS
				if (local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK)  //初始化的时候就是 NO_CLICK
				{
					double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
					double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
					if (select_coordinater_mode == CoordinateAxis::ModeType::CREATE) 
					{
						removeSelectedAll();
						pickEvent(winX, winY); // 用于加速 识别选中吧的 点 边 面
					}
					else 
					{
						pickEvent(winX, winY);
						selectedUpToDown();
					}
					break;
				}
				// 点击中CUBE
				/*
				if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
					switch (local_coordinater_axis->get_click_type()) {
						// XX
					case CoordinateAxis::ClickType::X_CUBE_CLICK:
						
						break;
					case CoordinateAxis::ClickType::Y_CUBE_CLICK:
						
						break;
					case CoordinateAxis::ClickType::Z_CUBE_CLICK:
						
						break;
					default:
						break;
					}
				}
				else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::CREATE) {
					if (selected_v_indices.size() != 1) {
						cout << "the number of selected_v must be only one!" << endl;
						break;
					}
					switch (local_coordinater_axis->get_click_type()) {
						// 生成
					case CoordinateAxis::ClickType::X_CUBE_CLICK:
						//hex_mesh_tool_->CreateCellByAixs(1, 0, 0, selected_v_indices[0]);
						break;
					case CoordinateAxis::ClickType::Y_CUBE_CLICK:
						//hex_mesh_tool_->CreateCellByAixs(0, 1, 0, selected_v_indices[0]);
						break;
					case CoordinateAxis::ClickType::Z_CUBE_CLICK:
						//hex_mesh_tool_->CreateCellByAixs(0, 0, 1, selected_v_indices[0]);
						break;
					default:
						break;
					}
				}
				*/
			}
			break;
		}

		case Skel_Move_Point:// 骨架模式 点的移动
		{
			if (btns == Qt::LeftButton)
			{
				
				//实时的鼠标位置
				double winX = (double)(_event->x() * 2.f - width()) / (width()); // 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height()); // 窗口坐标，范围在[-1,1]之内
				QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
				QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
				QVector4D pointEnd = point + 2 * toward;
				QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
				QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
				QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
				QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

				//获取骨架
				MeshKernel::CurveSkel* curve_skel_ = nullptr;
				if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
				if (!curve_skel_) return;

				double shortest_point_length = std::numeric_limits<double>::max();
				MeshKernel::iGameVertexHandle shortest_point_vh(-1);

				// 找到离鼠标最近的点，以及其距离
				dotLineCalculation->DisPointRay(point1, point2, shortest_point_length, shortest_point_vh);



				if (shortest_point_vh != -1 && shortest_point_length < avg_radius)
				{
					auto& shortest_point = curve_skel_->vertices(shortest_point_vh);
					local_coordinater_axis->set_orign(shortest_point.x(), shortest_point.y(), shortest_point.z());
					local_coordinater_axis->set_mode_type(CoordinateAxis::ModeType::SHOW);
					selected_v_indices.clear();
					selected_v_indices.push_back(shortest_point_vh);

					// 主动出发信号
					//std::cout << curve_skel_->vertex_radius_[shortest_point_vh] << std::endl;
					emit Send(QString::number(curve_skel_->vertex_radius_[shortest_point_vh]));
				//	DrawPointRadius(shortest_point_vh, curve_skel_->vertex_radius_[shortest_point_vh]);
					selected_v_radius = curve_skel_->vertex_radius_[shortest_point_vh];
					Move_Point_vh.clear();
					Move_Point_vh.push_back(shortest_point_vh);
					local_coordinater_axis->SetLen(avg_radius * 5);
				}

			}


			break;
		}

		 //这里实现了 点击 选取线段 和 绘制点的功能
		case Skel_Point_add: // 骨架 点的增加 
		{
			if (btns == Qt::LeftButton)
			{
				//实时的鼠标位置
				double winX = (double)(_event->x() * 2.f - width()) / (width()); // 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height()); // 窗口坐标，范围在[-1,1]之内
				//QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
				QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
				QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
				QVector4D pointEnd = point + 2 * toward;
				QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
				QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
				QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
				QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

				

				MeshKernel::CurveSkel* curve_skel_ = nullptr;
				if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
				if (!curve_skel_) return;



				double shortest_point_length = std::numeric_limits<double>::max();
				MeshKernel::iGameVertexHandle shortest_point_vh(-1);

				// 离射线最近点
				dotLineCalculation->DisPointRay(point1, point2, shortest_point_length, shortest_point_vh);


				if (shortest_point_vh != -1 && shortest_point_length < curve_skel_->vertex_radius_[shortest_point_vh] / 2)
				{
					copy_point = shortest_point_vh;
					
				}


				// 计算 直线和直线之间的距离，并确定距离最短的那个点
				MeshKernel::iGameEdgeHandle shortest_eh(-1);
				double shortest_lines_len = std::numeric_limits<double>::max();
				QVector3D temp_selected_point;



				dotLineCalculation->DisLineRay(point1, point2, shortest_lines_len, shortest_eh, temp_selected_point);

				
				// 把店添加进骨架中&& shortest_lines_len < avg_radius
				if (shortest_eh != -1  && is_copy_point_draw == false)
				{
					//只有在不是点选动画的时候才能增加点
					//点击选取了 点击事件距离最近的边 以及 这条边上的点
					selected_e_indices.clear();
					selected_e_indices.push_back(shortest_eh);
					add_possible_point = temp_selected_point;
					selected_add_possible_point.clear();
					selected_add_possible_point.push_back(add_possible_point);

					// 增加了点，对点进行了修改
					double length_all = curve_skel_->getLength(shortest_eh);
					auto& shortest_e = curve_skel_->edges(shortest_eh);
					auto& shortest_v1 = curve_skel_->vertices(shortest_e.vh1());
					auto& shortest_v2 = curve_skel_->vertices(shortest_e.vh2());
					MeshKernel::iGameVertex v_shortest(temp_selected_point.x(), temp_selected_point.y(), temp_selected_point.z());
					MeshKernel::iGameVertexHandle v_shortest_h = curve_skel_->AddVertex(v_shortest);
					MeshKernel::iGameEdgeHandle e_add_h1 = curve_skel_->AddEdge(shortest_e.vh1(), v_shortest_h);
					MeshKernel::iGameEdgeHandle e_add_h2 = curve_skel_->AddEdge(v_shortest_h, shortest_e.vh2());
					//curve_skel_->DeleteEdge(shortest_eh);
					double length_e1 = curve_skel_->getLength(e_add_h1);
					double length_e2 = curve_skel_->getLength(e_add_h2);
					double radiu = length_e2 / length_all * curve_skel_->vertex_radius_[shortest_e.vh1()] + length_e1 / length_all * curve_skel_->vertex_radius_[shortest_e.vh2()];
					curve_skel_->DeleteEdge(shortest_eh);
					curve_skel_->vertex_radius_.insert(std::make_pair(v_shortest_h, radiu));
					curve_skel_->updateAllHandles();// 清除section 和node
					//curve_skel_->createNode(); 
					proDraw->skel_tool_->createNode();//创造node和section
					SendToVAO_Skel();

					//update()后 选择增加点的邻边
					selected_e_indices.clear();
					selected_e_indices.push_back(e_add_h1);
					
				}


				

				
			}

			if (btns == Qt::RightButton)
			{

				//实时的鼠标位置
				double winX = (double)(_event->x() * 2.f - width()) / (width()); // 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height()); // 窗口坐标，范围在[-1,1]之内
				//QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
				QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
				QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
				QVector4D pointEnd = point + 2 * toward;
				QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
				QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
				QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
				QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());



				MeshKernel::CurveSkel* curve_skel_ = nullptr;
				if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
				if (!curve_skel_) return;


				double shortest_point_length = std::numeric_limits<double>::max();
				MeshKernel::iGameVertexHandle shortest_point_vh(-1);

				dotLineCalculation->DisPointRay(point1, point2, shortest_point_length, shortest_point_vh);

				if (shortest_point_vh != -1 && shortest_point_length < avg_radius / 2)
				{
					//printf("before the pick %d\n", copy_point);
					delet_point_vh = shortest_point_vh;
					//printf("after the pick %d\n", copy_point);
					//is_copy_point_draw = true;

					//删除对应的点
					if (curve_skel_->NeighborEh(delet_point_vh).size() == 1)
					{ //  end node
						//printf("copy_point is %d\n", copy_point);
						curve_skel_->vertex_radius_.erase(delet_point_vh); //删除了对应的半径信息
						curve_skel_->DeleteVertex(delet_point_vh);
						curve_skel_->updateAllHandles();
						selected_e_indices.clear();
						selected_v_indices.clear();
						copy_point.reset();
						selected_add_possible_point.clear();
						//copy_point.reset();
						//printf("copy_point is %d\n", copy_point);

						SendToVAO_Skel();
					}

					else if (curve_skel_->NeighborEh(delet_point_vh).size() == 2)
					{
						//printf("copy_point is %d\n", copy_point);
						auto set_vhs = curve_skel_->NeighborVh(delet_point_vh);
						auto it = set_vhs.begin();
						MeshKernel::iGameVertexHandle vh0 = *it;
						MeshKernel::iGameVertexHandle vh1 = *(++it);
						curve_skel_->AddEdge(vh0, vh1);
						curve_skel_->vertex_radius_.erase(delet_point_vh); //删除了对应的半径信息
						curve_skel_->DeleteVertex(delet_point_vh);
						curve_skel_->updateAllHandles();
						selected_e_indices.clear();
						selected_v_indices.clear();
						copy_point.reset();
						selected_add_possible_point.clear();
						//printf("copy_point is %d\n", copy_point);

						SendToVAO_Skel();
						//copy_point.reset();
					}
					else 
					{// >=3
						printf("refuse delete\n", copy_point);
					}
				}


			}

			break;
		}

		case COORDINATE: { // 变换模式下
			/*
			* 左键： 1.cell pick  2.坐标轴功能 分裂
			* 右键： 
			*/
			if (btns == Qt::LeftButton) 
			{
				if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) 
				{
					switch (local_coordinater_axis->get_click_type()) 
					{
					// 分裂
					case CoordinateAxis::ClickType::X_CUBE_CLICK:
						solutionSkel->spiltX();
						drawNormalSkel(); //   生成骨架 (打开骨架模式可以看到)
						break;

					case CoordinateAxis::ClickType::Y_CUBE_CLICK:
						solutionSkel->spiltY();
						drawNormalSkel();
						break;

					case CoordinateAxis::ClickType::Z_CUBE_CLICK:
						solutionSkel->spiltZ();
						drawNormalSkel();
						break;

					default:
						break;
					}
				}
				if (local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK)
				{
					double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
					double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
					pickEvent(winX, winY);
					update();
				}
			}
			break;
		}

		case PICK_CELLS_SKEL:   //废弃了
		{
			if (_event->button() == Qt::LeftButton)
			{
				double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
				pickEvent(winX, winY);
			}
			else 
			{
				SetEditModel(ModelDraw::DEFULTEDIT);
				setPickType(1);
			}
			break;
		}

		case SET_CELL_SKEL:  // 交互模式中 生成模式
		{
			/*
			* 左键： 1.cell pick(之前没选中，左点选中：选点) 2.cell add(之前选中过，左点到空白：增加点) 3.cell insert(之前选中过，左点选中另一个：增加点)
			* 右键： 1.cell not pick(之前选中过，右键点到空白处：取消选点)  2.cell del(之前点选中过，且再次右键选中：删除它)
			*/
			// 左键
			if (btns == Qt::LeftButton) 
			{
				// 点击了左键 就把左键位置记录下来
				double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内

				BMesh::Node* pre = solutionSkel->selectNode;

				initializePick();

				if (pre) // 之前选中过
				{	
					// 鼠标(winX winY) 与cell相交返回true 否则 false
					if (pickEvent_SKEL(winX, winY)) //3  cell insert 
					{
						solutionSkel->rinsertNode(pre);
						selected_c_indices.clear();
					}
					else
					{ // // 2   cell add，以之前点选的深度作为深度
						double p[3];
						double temp[3] = { pre->position.x, pre->position.y, pre->position.z };
						pick_point_with_depth(winX, winY, p, temp);
						solutionSkel->raddNode(p[0], p[1], p[2], 0, 0, 0, radius*0.05, radius*0.05, radius*0.05);
					}
				}
				else  // 之前没选中
				{
					pickEvent_SKEL(winX, winY);
				}
			}
			// 右键
			else { 
				double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
				BMesh::Node* pre = solutionSkel->selectNode; // 这个是哪里判定的没懂？？
				initializePick();
				if (pre) 
				{	
					if (pickEvent_SKEL(winX, winY)) 
					{	// 之前点选中过，且再次右键选中，删除它
						if(pre->getPosition() == solutionSkel->selectNode->getPosition())
							solutionSkel->rdelNode();
					}
					initializePick();
				}
			}
			break;
		}

		case CON_CELL_SKEL:  // 连接模式
		{ 
			/*
			* 左键： 1.cell pick(点到空白处要取消选点)  2.在已有选中后仍选中，连接两点或取消连接两点
			* 右键：
			*/
			if (btns == Qt::LeftButton) 
			{
				double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内

				BMesh::Node* pre = solutionSkel->selectNode;

				initializePick();

				if (pickEvent_SKEL(winX, winY) && pre)  // pickEvent_SKEL(winX, winY) 会更新 solutionSkel->selectNode
				{
					solutionSkel->smartConnectNodes(pre, solutionSkel->selectNode);
				}
				update();
			}
			break;
		}

		case TEST: { // 没懂这个
			if (_event->button() == Qt::LeftButton) {

				double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
				if (pickEvent_SKEL(winX, winY)) {
					std::cout << "MMMM" << std::endl;
					projectSkel->movePoint(solutionSkel->selectNode->root);
				}
				update();
			}
			break;
		}
		default:
			break;
	}

}



void ModelDraw::mouseMoveEvent(QMouseEvent* _event)
{

	switch (editmodel)
	{
	case DEFULTEDIT:  // 拓扑分裂 不需要知道
	{
		QPoint newPoint2D = _event->pos();
		//makeCurrent();
		if (lastpointok)
		{
			switch (mousemode)
			{
			case Qt::LeftButton:
				if (drawmodule == MMC_Convert) return;
				ModelDraw::Rotation(newPoint2D);
				break;

			case Qt::RightButton:
				Translation(newPoint2D);
				break;

			default:
				break;
			}
			lastpoint2 = newPoint2D;
			lastpointok = MapToSphere(lastpoint2, lastpoint3);
			update();
		}
		break;
	}

	case Skel_Move_Point: // 骨架模式 拖动事件
	{

		//实时的鼠标位置
		double winX = (double)(_event->x() * 2.f - width()) / (width()); // 窗口坐标，范围在[-1,1]之内
		double winY = -(double)(_event->y() * 2.f - height()) / (height()); // 窗口坐标，范围在[-1,1]之内
		//QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
		QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
		QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
		QVector4D pointEnd = point + 2 * toward;
		QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
		QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
		QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
		QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

		MeshKernel::CurveSkel* curve_skel_ = nullptr;
		if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
		if (!curve_skel_) return;


		double shortest_point_length = std::numeric_limits<double>::max();
		MeshKernel::iGameVertexHandle shortest_point_vh(-1);



		//计算平均半径 用于选点的距离
		if (IsCalAvgRadius == false)
		{
			avg_radius = 0;
			for (auto& vmap : curve_skel_->allvertices())
			{
				avg_radius += curve_skel_->vertex_radius_[vmap.first];
				/*if (max_radius < curve_skel_->vertex_radius_[vmap.first])
				{
					max_radius = curve_skel_->vertex_radius_[vmap.first];
				}
				if (min_radius > curve_skel_->vertex_radius_[vmap.first])
				{
					min_radius = curve_skel_->vertex_radius_[vmap.first];
				}*/
			}
			avg_radius /= curve_skel_->vertex_radius_.size();
			dotLineCalculation->AveRadiusCalculationTool = avg_radius; //传给计算工具 用于点选
			IsCalAvgRadius = true;
		}



		//拖动旋转
		QPoint newPoint2D = _event->pos();

		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));

		Qt::MouseButtons btns = _event->buttons();

		switch (btns)
		{
		case Qt::NoButton:
		{
			// coordinater PICK
			local_coordinater_axis->PickEvent(*this, newWin.x(), newWin.y()); //这个函数的作用是什么  1.返回是否和坐标轴相交  2.把与什么相交记录在click_type_中
			

			dotLineCalculation->DisPointRay(point1, point2, shortest_point_length, shortest_point_vh);



			if (shortest_point_vh != -1 && shortest_point_length < avg_radius)
			{
				auto& shortest_point = curve_skel_->vertices(shortest_point_vh);
				selected_v_draw_indices.clear();
				selected_v_draw_indices.push_back(shortest_point_vh); //压入要放大的点 
				MoveSelectedPoint = shortest_point_vh;
			}
			else
			{
				selected_v_draw_indices.clear();
				MoveSelectedPoint.reset();
			}

			double Origin_Cirle_raidus_Dis = -1;
			//实现是否选中了当前点半径，进行修改
			QVector3D PointRadiusXoy = local_coordinater_axis->Point_Xoy; //与XOY面的交点
			if (PointRadiusXoy.x() > 0 && PointRadiusXoy.y() > 0 && Move_Point_vh.size() > 0)
				//半径修改只在正半轴有意义
			{
				// /3 半径变化不然太灵敏
				// 减数是画的半径
				//printf("%lf %lf \n",PointRadiusXoy.x(), PointRadiusXoy.y());
				double dis_mouse_radius = abs(PointRadiusXoy.length()  - curve_skel_->vertex_radius_[Move_Point_vh[0]]);
				if (dis_mouse_radius < avg_radius / 5)
				{
					is_change_radius = true;
				}
				else
				{
					is_change_radius = false;
				}
			}
			else
			{
				is_change_radius = false;
			}

			break;



		}


		case Qt::LeftButton: //拖动 点击左键（不要忘记拖动）
		{
			if (is_change_radius == false)
			{
				switch (local_coordinater_axis->get_click_type())
				{
					// 拖动旋转
				case CoordinateAxis::ClickType::NO_CLICK:
				{
					// rotation
					if (lastpointok)
					{
						Rotation(newPoint2D);
						lastpoint2 = newPoint2D;
						lastpointok = MapToSphere(lastpoint2, lastpoint3);
					}
					break;
				}


				case CoordinateAxis::ClickType::X_AXIS_CLICK:
				case CoordinateAxis::ClickType::Y_AXIS_CLICK:
				case CoordinateAxis::ClickType::Z_AXIS_CLICK:
				{
					is_edit_point = true;
					if (lastpointok)
					{
						double x, y, z;
						local_coordinater_axis->get_orign(x, y, z);
						// 移动
						if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::SHOW) // 变换模式下  TRANSFROM
						{
							//上锁 主要是拖动的时候不要显示别的点
							is_edit_point = true;
							// 返回相对于坐标轴原点的方向
							QVector3D dir = local_coordinater_axis->MoveEvent(*this, oldWin, newWin);
							//printf("the direction is: %lf, %lf, %lf\n", dir.x(), dir.y(), dir.z());
							// 生成cell 以及调整对应的角度
							// 而且是一个不间断的过程
							auto& vh = selected_v_indices[0];
							proDraw->skel_tool_->movePoint((VH)vh, dir.x(), dir.y(), dir.z());
							updateCurrentSkelSphere((VH)vh);

							MeshKernel::CurveSkel* curve_skel_ = nullptr;
							if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
							if (!curve_skel_) return;

							for (auto& eh : curve_skel_->NeighborEh((VH)vh))
							{
								updateCurrentCylinderSkel(eh);
							}

							is_edit_point = false;


						}


						lastpoint2 = newPoint2D; // 是一个不间断的过程
					}
					break;
				}



				default:
					break;
				}
			}

			if (is_change_radius)
			{
				if (lastpointok)
				{
					
					{
						double change_raidus = local_coordinater_axis->ChangeRadius(*this, oldWin, newWin);
						
						curve_skel_->vertex_radius_[Move_Point_vh[0]] += change_raidus ; //这样改变半径比较流畅
						updateCurrentSkelSphere(Move_Point_vh[0]);

						for (auto& eh : curve_skel_->NeighborEh(Move_Point_vh[0]))
						{
							updateCurrentCylinderSkel(eh);
						}
					}

					emit Send(QString::number(curve_skel_->vertex_radius_[Move_Point_vh[0]]));

					lastpoint2 = newPoint2D; // 是一个不间断的过程
				}
			}

			break;
		}

		case Qt::RightButton:
		{
			if (lastpointok && local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				// translation
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}

		break;
	}


	case Skel_Point_add: // 骨架模式 点增加模式
	{
		//实时的鼠标位置
		double winX = (double)(_event->x() * 2.f - width()) / (width()); // 窗口坐标，范围在[-1,1]之内
		double winY = -(double)(_event->y() * 2.f - height()) / (height()); // 窗口坐标，范围在[-1,1]之内
		//QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
		QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
		QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
		QVector4D pointEnd = point + 2 * toward;
		QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
		QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
		QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
		QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

		MeshKernel::CurveSkel* curve_skel_ = nullptr;
		if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
		if (!curve_skel_) return;
		


		//拖动旋转
		QPoint newPoint2D = _event->pos();

		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));

		Qt::MouseButtons btns = _event->buttons();

		switch (btns)
		{
		case Qt::NoButton:
		{
			// coordinater PICK
			local_coordinater_axis->PickEvent(*this, newWin.x(), newWin.y()); //这个函数的作用是什么  1.返回是否和坐标轴相交  2.把与什么相交记录在click_type_中

			if (IsCalAvgRadius == false)
			{
				avg_radius = 0;
				for (auto& vmap : curve_skel_->allvertices())
				{
					avg_radius += curve_skel_->vertex_radius_[vmap.first];
					/*if (max_radius < curve_skel_->vertex_radius_[vmap.first])
					{
						max_radius = curve_skel_->vertex_radius_[vmap.first];
					}
					if (min_radius > curve_skel_->vertex_radius_[vmap.first])
					{
						min_radius = curve_skel_->vertex_radius_[vmap.first];
					}*/
				}
				avg_radius /= curve_skel_->vertex_radius_.size();
				dotLineCalculation->AveRadiusCalculationTool = avg_radius; //传给计算工具 用于点选
				IsCalAvgRadius = true;
			}

			MeshKernel::iGameEdgeHandle shortest_eh(-1);
			double shortest_lines_len = std::numeric_limits<double>::max();
			QVector3D temp_selected_point;


			dotLineCalculation->DisLineRay(point1, point2, shortest_lines_len, shortest_eh, temp_selected_point);

			//&& shortest_lines_len < avg_radius
			if (shortest_eh != -1 )
			{
				//拖动 实时渲染出选取的边和点
				selected_e_indices.clear();
				selected_e_indices.push_back(shortest_eh); 
				add_possible_point = temp_selected_point;
				selected_add_possible_point.clear();
				selected_add_possible_point.push_back(add_possible_point);
				//把要渲染的可能的球体数据传进去
				updatePossibleSphere();
			}
			else
			{
				selected_e_indices.clear();
				selected_add_possible_point.clear();

			}
			


			double shortest_point_length = std::numeric_limits<double>::max();
			MeshKernel::iGameVertexHandle shortest_point_vh(-1);

			
			dotLineCalculation->DisPointRay(point1, point2, shortest_point_length, shortest_point_vh);

			
			if (shortest_point_vh != -1 && shortest_point_length < curve_skel_->vertex_radius_[shortest_point_vh] / 2)
			{

				copy_point = shortest_point_vh;
				is_copy_point_draw = true;
				shortest_eh.reset();
				selected_e_indices.clear();
				selected_add_possible_point.clear();
			}
			else 
			{
				is_copy_point_draw = false;
				copy_point.reset();
			}

			break;
		}

		case Qt::LeftButton: //拖动 点击左键（不要忘记拖动）
		{
			switch (local_coordinater_axis->get_click_type())
			{
				// 拖动旋转
			case CoordinateAxis::ClickType::NO_CLICK:
			{
				// rotation
				if (lastpointok && copy_point == -1) 
				{
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
				break;
			}
		
			


			default:
				break;
			}

			// 1外层还要加上选中点的判断
			if (copy_point != -1)
			{
				auto& copy_v = curve_skel_->vertices(copy_point);
				double temp[3] = { copy_v.x() ,copy_v.y(),  copy_v.z() };

				//MeshKernel::iGameVertexHandle vh(-1); //记录一下新增的 vh
				if (copyTime_point)
				{
					double p[3];//移动的点
					pick_point_with_depth(newWin.x(), newWin.y(), p, temp);
					//这里写上拖动点的操作
					curve_skel_->vertices(vh_copy_add).setX(p[0]);
					curve_skel_->vertices(vh_copy_add).setY(p[1]);
					curve_skel_->vertices(vh_copy_add).setZ(p[2]);
					updateCurrentSkelSphere(vh_copy_add);
					updateCurrentCylinderSkel(e_add_h1);
					//SendToVAO_Skel();

				}

				else
				{
					//这个地方需要1.复制节点 2.对骨架进行修改
					MeshKernel::iGameVertex v_iGame(copy_v.x(), copy_v.y(), copy_v.z() + 0.01);
					vh_copy_add = curve_skel_->AddVertex(v_iGame);
					e_add_h1 = curve_skel_->AddEdge(copy_point, vh_copy_add);
					curve_skel_->vertex_radius_.insert(std::make_pair(vh_copy_add, curve_skel_->vertex_radius_[copy_point]));
					// 这里加了一个点 不需要curve_skel_->updateAllHandles();
					//SendToVAO_Skel();
					InitVaoSphereSize();// 和顶点个数有关
					InitVaoCylinderSize();// 和边的个数有关
					copyTime_point = true;


				}
			}

			break;
		}

		case Qt::RightButton:
		{
			if (lastpointok && local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				// translation
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}




			break;
		} 
		default:
			break;
		}

		break;
	}
	


	case PICKITEM:  // 启动点选
	 {
		if (_event->buttons() == Qt::LeftButton)
		{
			mouse_move_flag = true;

			// TODO: 绘制框选区域
			QVector2D point0 = QVector2D((double)(lastpoint2.x() * 2 - width()) / (width()), -(double)(lastpoint2.y() * 2 - height()) / (height()));
			QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
			double x_min = std::min(point0.x(), point1.x()), x_max = std::max(point0.x(), point1.x());
			double y_min = std::min(point0.y(), point1.y()), y_max = std::max(point0.y(), point1.y());

			// std::cout << "ret: (" << x_min << ", " << y_min << ") --> (" << x_max << ", " << y_max << ")\n";
			// 初始化为设备坐标

			std::vector<QVector4D> ret =
			{
				QVector4D(x_min, y_min, 0.9, 1.0), QVector4D(x_min, y_max, 0.9, 1.0),
				QVector4D(x_max, y_min, 0.9, 1.0), QVector4D(x_max, y_max, 0.9, 1.0)
			};

			// 转化为世界坐标
			for (auto& pos : ret)
			{
				pos = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pos;
				pos /= pos.w();
			}

			std::vector<int> draw_indices =
			{
				0, 1, 0, 2, 3, 1, 3, 2
			};

			selected_rectangle.clear();

			for (auto idx : draw_indices)
			{
				auto& pos4d = ret[idx];
				QVector3D pos3d(pos4d.x(), pos4d.y(), pos4d.z());
				selected_rectangle.push_back(pos3d);
			}

			update();
		}
		break;
	}


	case POINT_COORDINATE: // 启动移动按钮
	{

		/*
		* 移动：坐标轴PICK
		* 移动左键： 整体旋转   坐标轴触发事件 : 移动, 旋转，XX拉伸
		* 移动右键： 整体平移
		*/

		QPoint newPoint2D = _event->pos();

		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));

		Qt::MouseButtons btns = _event->buttons();
		switch (btns) {
		case Qt::NoButton:
		{
			// coordinater PICK
			local_coordinater_axis->PickEvent(*this, newWin.x(), newWin.y());
			break;
		}
		case Qt::LeftButton: {
			switch (local_coordinater_axis->get_click_type()) {
			case CoordinateAxis::ClickType::NO_CLICK: {
				// rotation
				if (lastpointok) {
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
				break;
			}
			case CoordinateAxis::ClickType::X_AXIS_CLICK:
			case CoordinateAxis::ClickType::Y_AXIS_CLICK:
			case CoordinateAxis::ClickType::Z_AXIS_CLICK: {
				if (lastpointok) {
					// 移动
					QVector3D dir = local_coordinater_axis->MoveEvent(*this, oldWin, newWin);
					if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::TRANSFROM) {
						hex_mesh_tool_->MovePoints(dir.x(), dir.y(), dir.z(), selected_v_indices);
						updateDrawData();
					}
					// 拉伸
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
						//QVector3D dir = local_coordinater_axis->FlexibleEvent(*this, oldWin, newWin);
						//solutionSkel->rXyzCubeByMouse(x, y, z, rlt->root->getLength() + dir.y(), rlt->root->getWidth() + dir.x(), rlt->root->getHeight() + dir.z());
					}
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::CREATE) {
						if (copyTime) {
							// 移动复制品
							hex_mesh_tool_->MoveFaces(dir.x(), dir.y(), dir.z(), selected_f_indices);
						}
						else {
							// 生成一个复制品 改变选中面
							hex_mesh_tool_->CreateCellByAixs(dir.x(), dir.y(), dir.z(), selected_f_indices);
							copyTime = true;
						}
						updateDrawData();
					}
					lastpoint2 = newPoint2D;
				}
				break;
			}
			case CoordinateAxis::ClickType::X_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Y_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Z_CIRCLE_CLICK: {
				// 旋转
				if (lastpointok) {
					QVector4D quat = local_coordinater_axis->RotateEvent(*this, oldWin, newWin);
					double center_x, center_y, center_z;
					local_coordinater_axis->get_orign(center_x, center_y, center_z);
					if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::TRANSFROM) {
						hex_mesh_tool_->RotatePoints(Vec(center_x, center_y, center_z), quat[0], quat[1], quat[2], quat[3], selected_v_indices);
					}
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
					}
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::CREATE) {
					}
					updateDrawData();
					lastpoint2 = newPoint2D;
				}
				break;
			}
			default:
				break;
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok && local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				// translation
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	case COORDINATE:
	{ //变换模式
		/*
		* 移动：坐标轴PICK
		* 移动左键： 整体旋转   坐标轴触发事件 : 移动, 旋转，拉伸
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();

		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));

		Qt::MouseButtons btns = _event->buttons();

		switch (btns)
		{
		case Qt::NoButton:
		{
			// coordinater PICK
			local_coordinater_axis->PickEvent(*this, newWin.x(), newWin.y()); //这个函数的作用是什么  1.返回是否和坐标轴相交  2.把与什么相交记录在click_type_中
			break;
		}

		case Qt::LeftButton: //拖动 点击左键（不要忘记拖动）
		{
			switch (local_coordinater_axis->get_click_type())
			{
			case CoordinateAxis::ClickType::NO_CLICK:
			{
				// rotation
				if (lastpointok)
				{
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
				break;
			}
			case CoordinateAxis::ClickType::X_AXIS_CLICK:
			case CoordinateAxis::ClickType::Y_AXIS_CLICK:
			case CoordinateAxis::ClickType::Z_AXIS_CLICK:
			{
				if (lastpointok)
				{
					double x, y, z;
					local_coordinater_axis->get_orign(x, y, z);
					BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(x, y, z)); // solutionSkel是干嘛的
					// 移动
					if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::TRANSFROM) // 变换模式下  TRANSFROM
					{
						// 返回相对于坐标轴原点的方向
						QVector3D dir = local_coordinater_axis->MoveEvent(*this, oldWin, newWin);
						// 生成cell 以及调整对应的角度
						// 而且是一个不间断的过程
						// 这里是生成一个六面体，长度什么的都可以调整 是后面的3个参数
						solutionSkel->rXyzCubeByMouse(x + dir.x(), y + dir.y(), z + dir.z(), rlt->root->getLength(), rlt->root->getWidth(), rlt->root->getHeight());
					}

					// 拉伸
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT)
					{
						QVector3D dir = local_coordinater_axis->FlexibleEvent(*this, oldWin, newWin); // 和 MoveEvent 的区别，这里是不改变原点的，上面是要改变原点的
						solutionSkel->rXyzCubeByMouse(x, y, z, rlt->root->getLength() + dir.y(), rlt->root->getWidth() + dir.x(), rlt->root->getHeight() + dir.z());
					}

					lastpoint2 = newPoint2D; // 是一个不间断的过程
				}
				break;
			}

			case CoordinateAxis::ClickType::X_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Y_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Z_CIRCLE_CLICK:
			{
				// 旋转
				if (lastpointok)
				{
					QVector4D quat = local_coordinater_axis->RotateEvent(*this, oldWin, newWin);
					solutionSkel->rRotateCubeByMouse(V3f(quat[0], quat[1], quat[2]), quat[3]);
					lastpoint2 = newPoint2D;
				}
				break;
			}
			default:
				break;
			}
			break;
		}

		case Qt::RightButton:
		{
			if (lastpointok && local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				// translation
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	case PICK_CELLS_SKEL: { // 废弃了
		if (_event->button() == Qt::LeftButton) {
			double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickEvent(winX, winY);
		}
		else {
			SetEditModel(ModelDraw::DEFULTEDIT);
			setPickType(1);
		}
		break;
	}

	case SET_CELL_SKEL: //生成模式
	{
		/*
		* 移动：
		* 移动左键： （未选中）整体旋转   （选中）复制生成
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();
		switch (btns)
		{
		case Qt::LeftButton:
		{
			if (solutionSkel->selectNode)
			{// 选中复制生成
				BMesh::Node* sNode = solutionSkel->selectNode;
				double temp[3] = { sNode->position.x, sNode->position.y, sNode->position.z };
				if (copyTime) //初始为false，
				{
					// 移动复制品
					double p[3];//p存储了移动的点
					pick_point_with_depth(newWin.x(), newWin.y(), p, temp);
					//
					solutionSkel->rXyzCubeByMouse(p[0], p[1], p[2], sNode->length, sNode->width, sNode->height);
				}
				else
				{
					// 生成一个复制品
					solutionSkel->rcopyNode();
					selected_c_indices.clear();
					selected_c_indices.push_back(findNearestCell(temp));
					copyTime = true;
				}
			}
			else
			{
				if (lastpointok)
				{
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok) {
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	case CON_CELL_SKEL:
	{ // 链接模式
		/* 【基础版】
		* 移动：
		* 移动左键： 整体旋转
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();

		switch (btns)
		{
		case Qt::LeftButton:
		{
			if (lastpointok)
			{
				Rotation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		case Qt::RightButton:
		{
			if (lastpointok)
			{
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	case TEST: {
		/* 【基础版】
		* 移动：
		* 移动左键： 整体旋转
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();
		switch (btns) {
		case Qt::LeftButton: {
			if (lastpointok) {
				Rotation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok) {
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	default:
		break;
	}

	update();
	
}


void ModelDraw::mouseReleaseEvent(QMouseEvent * _event) 
{

	//assert(mousemode < N_MOUSE_MODES);
	switch (editmodel) 
	{
	case PICKITEM:  // 启动点选功能
	{// 框选功能 Beta

		selected_rectangle.clear();
		QPoint currentpoint2 = _event->pos();
		int x_min = std::min(lastpoint2.x(), currentpoint2.x()), x_max = std::max(lastpoint2.x(), currentpoint2.x());
		int y_min = std::min(lastpoint2.y(), currentpoint2.y()), y_max = std::max(lastpoint2.y(), currentpoint2.y());
		int red, green, blue, sum_color_threshold = 75, signal_color_threshold = 30;// threshold 越大, 采样点越少, 速度越快, 但可能出现精度问题
		QScreen* screen = QApplication::primaryScreen();
		QPixmap pixmap = screen->grabWindow(this->winId(), x_min, y_min, x_max - x_min + 1, y_max - y_min + 1);
		QImage image = pixmap.toImage();
		image.save("Image/captureImg.jpg", "jpg"); // 保存截图

		int w = image.width(), h = image.height();

		vector<vector<bool>> visited(w, vector<bool>(h, false));
		vector<QPoint> dirs = 
		{
			QPoint(-1, 0), QPoint(1, 0), QPoint(0, -1), QPoint(0, 1),
			QPoint(-1, -1), QPoint(-1, 1), QPoint(1, -1), QPoint(1, 1)
		};
		vector<QPoint> centroids;
		for (int i = 0; i < w; ++i) 
		{
			for (int j = 0; j < h; ++j) 
			{
				if (visited[i][j]) continue;
				queue<pair<int, int>> que;
				que.emplace(i, j);
				visited[i][j] = true;
				QPoint sum(i, j);
				int cnt = 1;
				while (!que.empty()) {
					auto pospair = que.front();
					que.pop();
					int x = pospair.first, y = pospair.second;
					QPoint pos(x, y);
					const QColor& color = image.pixelColor(pos);
					red = color.red(), green = color.green(), blue = color.blue();
					for (auto& move : dirs) {
						QPoint adj = pos + move;
						int _x = adj.x(), _y = adj.y();
						if (_x < 0 || _x >= w || _y < 0 || _y >= h || visited[_x][_y]) continue;
						const QColor& _c = image.pixelColor(adj);
						int diff_r = abs(_c.red() - red);
						int diff_g = abs(_c.green() - green);
						int diff_b = abs(_c.blue() - blue);
						int diff_signal = std::max(diff_r, std::max(diff_g, diff_b));
						int diff_sum = diff_r + diff_g + diff_b;
						if (diff_signal < signal_color_threshold && diff_sum < sum_color_threshold) {
							que.emplace(_x, _y);
							visited[_x][_y] = true;
							sum += adj;
							cnt++;
						}
					}
				}
				sum /= cnt;
				sum.setX(sum.x() + x_min);
				sum.setY(sum.y() + y_min);
				centroids.push_back(sum);
			}
		}

		time_t time1 = clock();
		for (auto& pos : centroids) 
		{
			double winX = (double)(pos.x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(pos.y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickEvent(winX, winY);
		}
		int pick_cost = clock() - time1;
		std::cout << "Sample points size = " << centroids.size() << ", cost " << pick_cost << "ms.\n";
		lastpoint2 = currentpoint2;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
		break;
	}

	case POINT_COORDINATE: 
	{
		if (copyTime)
		{ // 释放鼠标时的操作
			copyTime = false;
		}
		break;
	}

	case SET_CELL_SKEL: //生成模式
	{
		/*
		* 释放左键： （选中）将  复制时间与选中 重置
		* 释放右键：
		*/
		if (copyTime) 
		{
			copyTime = false;
			initializePick(); // 取消选中的节点
		}

		
		break;

	}

	case Skel_Move_Point:
	{
		/*Qt::MouseButtons btns = _event->buttons();*/
		//switch (btns)
		//{
		//	case Qt::LeftButton: // 拖动结束后改为true;
		//	{
		//		if (Can_send_VAO)
		//		{
		//			SendToVAO_Skel();
		//			Can_send_VAO = false;
		//		}


		//		break;
		//	}

		//	default:
		//		break;

		//}
		//if (Can_send_VAO)
		//{
		//	//SendToVAO_Skel();
		//	SendChangeToVAO_Skel();
		//	Can_send_VAO = false;
		//}


		//SendToVAO_Skel();
		
		break;
	}

	case Skel_Point_add:
	{
		if (copyTime_point)
		{
			copyTime_point = false;
			copy_point.reset();
			selected_add_possible_point.clear();
			selected_e_indices.clear();
			is_copy_point_draw = false;
		}
		/*else
		{
			SendToVAO_Skel();
		}*/

		break;
	}

	default:
		break;
	}
	mouse_move_flag = false;

	lastpointok = false;
	mousemode = Qt::NoButton;
}


void ModelDraw::wheelEvent(QWheelEvent * _event) 
{


	// 滚轮拉伸
	if (editmodel == COORDINATE && local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
		//变换模式下 
		double x, y, z;
		local_coordinater_axis->get_orign(x, y, z);
		BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(x, y, z));
		double d = -_event->delta() / 120.0 * 0.01 * radius; // delta()返回鼠标滚轮的情况
		solutionSkel->rXyzCubeByMouse(x, y, z, rlt->root->getLength() + d, rlt->root->getWidth() + d, rlt->root->getHeight() + d);
	}

	else 
	{
		//整体缩放
		// 没懂
		// 0.1代表缩放速度
		double d = -_event->delta() / 120.0 * 0.1 * radius;
		Translate(Vector3d(0.0, 0.0, d));
		//local_coordinater_axis->set_len(local_coordinater_axis->get_len() + d * 0.05 );
		local_coordinater_axis->set_len(local_coordinater_axis->get_len() );
	}
	update();
}

void ModelDraw::keyPressEvent(QKeyEvent * _event) {
	// 没懂 会出现异常
	switch (_event->key()) {
	case Qt::Key_Escape:
		qApp->quit();
		break;
	case Qt::Key_Delete:
	case Qt::Key_D:
		deleteSelectedAll();
		break;
	case Qt::Key_A:
		add_new_cell();
		break;
	}
	_event->ignore();
}

void ModelDraw::keyReleaseEvent(QKeyEvent * _event) {
	_event->ignore();
}

void ModelDraw::Rotation(const QPoint& p) 
{
	Vector3d  newPoint3D;
	bool newPoint_hitSphere = MapToSphere(p, newPoint3D);
	if (newPoint_hitSphere) {
		Vector3d axis = lastpoint3.Cross(newPoint3D);// corss product
		if (axis.Length() < 1e-7) {
			axis = Vector3d(1, 0, 0);
		}
		else {
			axis.Normalize();
		}
		// find the amount of rotation
		Vector3d d = lastpoint3 - newPoint3D;
		double t = 0.5 * d.Length() / trackballradius;
		if (t < -1.0) t = -1.0;
		else if (t > 1.0) t = 1.0;
		double phi = 2.0 * asin(t);
		double angle = phi * 180.0 / M_PI;

		Rotate(axis, angle);
	}
}

void ModelDraw::translationCell(QPoint p) {
	double winX = (double)(p.x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
	double winY = -(double)(p.y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
	double pLocal[3] = { 0 };
	// gluUnProject(winX, winY, depth, modelview_matrix_, projection_matrix_, viewport, &pLocal[0], &pLocal[1], &pLocal[2]);
	double item[3]{ px, py, pz };
	pick_point_with_depth(winX, winY, pLocal, item);
	double relativePosition[3] = { 0 };
	relativePosition[0] = pLocal[0] - px;
	relativePosition[1] = pLocal[1] - py;
	relativePosition[2] = pLocal[2] - pz;
	px = pLocal[0]; py = pLocal[1]; pz = pLocal[2];
	//SingalFather::getInstance()->sdeltaEditCell(relativePosition[0], relativePosition[1], relativePosition[2], -1);
}

void ModelDraw::SetEditModel(const EditModel & editmodel) 
{
	this->editmodel = editmodel;
}


/**************************** 渲染函数 ******************************************/
void ModelDraw::EnableLighting(bool b) {
	isEnableLighting = b;
	update();
}

void ModelDraw::LoadRotation(void) {
	LoadCopyModelViewMatrix();
	update();
}

void ModelDraw::ResetView(void) {
	ResetModelviewMatrix();
	viewCenter();
	update();
}

void ModelDraw::getPatchStatus(bool status) {
	this->PatchChange = status;
}

void ModelDraw::viewCenter() {
	update();
}

void ModelDraw::ModelClear() {

	// 清除绘制数据
	this->vertices_data.clear();
	vector<GLfloat>().swap(this->vertices_data);
	tri_cnt = 0;
	quad_cnt = 0;

	// 清除模型数据
	delete this->hexMesh;
	this->hexMesh = nullptr;
	delete this->triMesh;
	this->triMesh = nullptr;
	delete this->tetMesh;
	this->tetMesh = nullptr;
	delete this->quadMesh;
	this->quadMesh = nullptr;

	// 清除点选数据
	this->selected_v_indices.clear();
	this->selected_e_indices.clear();
	this->selected_f_indices.clear();
	this->selected_c_indices.clear();
	vector<int>().swap(this->selected_v_indices);
	vector<int>().swap(this->selected_e_indices);
	vector<int>().swap(this->selected_f_indices);
	vector<int>().swap(this->selected_c_indices);
	mouse_move_flag = false;

	force_id.clear();
	spc_id.clear();
}


void ModelDraw::TestloadShader()
{
	// 相关数据传入sahder
	shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
	shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 32.0f);

	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &TestVAO); // 这里的1都是数量
	core->glDeleteBuffers(1, &TestVBO);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &TestVAO);
	core->glGenBuffers(1, &TestVBO);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(TestVAO);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, TestVBO);
	core->glBufferData(GL_ARRAY_BUFFER, Testvertices_data.size() * sizeof(GLfloat), Testvertices_data.data(), GL_STATIC_DRAW); // data()返回指针的地址
	// 4. 设定顶点属性指针
	// 顶点坐标传入VBO，location = 0
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
	core->glEnableVertexAttribArray(0);
	// 顶点 法向/颜色 传入VBO，location=1
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
	core->glEnableVertexAttribArray(1); //用于激活我们设置的数据的步长
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);
}

void ModelDraw::loadShader() {

	// 相关数据传入sahder
	shaderProgram.setUniformValue("model", modelmatrix);	          // 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix);               // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix);   // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix);           // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);      // 视角方向
	shaderProgram.setUniformValue("light.position", 0.0f, 1.0f, 1.0f);  // 光照方向
	shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 32.0f);

	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &VAO); // 这里的1都是数量
	core->glDeleteBuffers(1, &VBO);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &VAO);
	core->glGenBuffers(1, &VBO);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(VAO);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, VBO);
	core->glBufferData(GL_ARRAY_BUFFER, vertices_data.size() * sizeof(GLfloat), vertices_data.data(), GL_STATIC_DRAW); // data()返回指针的地址
	// 4. 设定顶点属性指针
	// 顶点坐标传入VBO，location = 0
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
	core->glEnableVertexAttribArray(0);
	// 顶点 法向/颜色 传入VBO，location=1
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
	core->glEnableVertexAttribArray(1); //用于激活我们设置的数据的步长
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);

}

void ModelDraw::updateDrawData_HexMesh() {
	if (!hexMesh) return;

	auto& cells = hexMesh->allcells();
	auto& edges = hexMesh->alledges();

	
	// 更新绘制数据
	unordered_map<int, bool> isDrawed;
	for (auto& cp : cells) {// 更新绘制面
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		Vector3d color(0.349, 0.6, 0.349);
		//Vector3d color(0.725, 0.725, 0.725);
		if (show_mesh_quality) color = colorbar.getColor(1 - MeshMath::get_min_scaled_Jacobian_HexMesh(*hexMesh, cp.first));
		auto fhs = cp.second.getFaceHandle();
		
		for (auto& fh : fhs) {
			if (isDrawed.count(fh) || !hexMesh->isOnBoundary(fh)) continue;

			//if (hexMesh->NeighborCh(cp.first).size() < 3) color = Vector3d(1, 0, 0);// 标记表面冗余体
			isDrawed[fh] = true;
			auto vhs = hexMesh->faces(fh).getVertexHandle();
			for (int i = 2; i < vhs.size(); ++i) {
				std::vector<VH> trivhs = { vhs[0], vhs[i - 1], vhs[i] };
				for (auto& vh : trivhs) {
					//std::cout << vh << std::endl;
					auto& v = hexMesh->vertices(vh);
					vertices_data.push_back(v.x());
					vertices_data.push_back(v.y());
					vertices_data.push_back(v.z());
					vertices_data.push_back(color[0]);
					vertices_data.push_back(color[1]);
					vertices_data.push_back(color[2]);
				}

			}
		}
	}
	tri_cnt = vertices_data.size() / 18;

	for (auto& ep : edges) {// 更新绘制线
		if (!hexMesh->isOnBoundary(ep.first)) continue;
		auto& v0 = hexMesh->vertices(ep.second.vh(0));
		vertices_data.push_back(v0.x());
		vertices_data.push_back(v0.y());
		vertices_data.push_back(v0.z());
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		auto& v1 = hexMesh->vertices(ep.second.vh(1));
		vertices_data.push_back(v1.x());
		vertices_data.push_back(v1.y());
		vertices_data.push_back(v1.z());
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		vertices_data.push_back(0);
	}
	line_cnt = (vertices_data.size() - tri_cnt * 18) / 12;

	update();

}

void ModelDraw::updateDrawData_QuadMesh() {
	if (!quadMesh) return;

	for (auto& fp : quadMesh->allfaces()) {// 更新绘制面
		const auto& vhs = fp.second.getVertexHandle();
		Vector3d color(0.349, 0.6, 0.349);
		if (show_mesh_quality) color = colorbar.getColor(1 - MeshMath::get_min_angle_quality_SurfaceMesh(*quadMesh, fp.first));
		for (int i = 2; i < vhs.size(); ++i) {
			std::vector<VH> trivhs = { vhs[0], vhs[i - 1], vhs[i] };
			for (auto& vh : trivhs) {
				auto& v = quadMesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}

		}
	}
	tri_cnt = vertices_data.size() / 18;

	for (auto& ep : quadMesh->alledges()) {// 更新绘制线
		auto& v0 = quadMesh->vertices(ep.second.vh(0));
		vertices_data.push_back(v0.x());
		vertices_data.push_back(v0.y());
		vertices_data.push_back(v0.z());
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		auto& v1 = quadMesh->vertices(ep.second.vh(1));
		vertices_data.push_back(v1.x());
		vertices_data.push_back(v1.y());
		vertices_data.push_back(v1.z());
		vertices_data.push_back(0);
		vertices_data.push_back(0);
		vertices_data.push_back(0);
	}
	line_cnt = (vertices_data.size() - tri_cnt * 18) / 12;

	update();

}

void ModelDraw::updateDrawData_TetMesh() {
	if (!tetMesh) return;

	for (auto& cp : tetMesh->allcells()) {
		if (!tetMesh->isOnBoundary(cp.first)) continue;
		Vector3d color(0.349, 0.6, 0.349);
		if (show_mesh_quality) {
			double ar = MeshMath::get_aspect_ratio_TetMesh(*tetMesh, cp.first);
			if (ar < 1.0) color = colorbar.getColor(0);// 质量很好
			else if (ar > 256) color = colorbar.getColor(0);// 质量很差
			else color = colorbar.getColor((std::log(ar) / std::log(2)) / 8);// 质量一般
		}
		for (auto& fh : cp.second.getFaceHandle()) {
			if (!tetMesh->isOnBoundary(fh)) continue;
			auto& face = tetMesh->faces(fh);
			const auto& vhs = face.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = tetMesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}
		}
	}

	tri_cnt = vertices_data.size() / 18;

	update();

}

void ModelDraw::updateDrawData_TetMesh_CutPlane() {
	if (!tetMesh) return;

	std::vector<GLfloat> translucence_vertices_data;

	/*double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];*/
	//printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);
	double ar_max = 0;
	auto& cells = tetMesh->allcells();
	std::unordered_map<int, bool> current_drawed;
	for (auto& cp : cells) {
		bool drawFlag = true;
		Vector3d color = (tetMesh->isOnBoundary(cp.first)) ? Vector3d(0.349, 0.6, 0.349) : Vector3d(0.6588, 0.6588, 0.0862);
		if (show_mesh_quality) {
			double ar = MeshMath::get_aspect_ratio_TetMesh(*tetMesh, cp.first);
			ar_max = std::max(ar_max, ar);
			if (ar < 1.0) color = colorbar.getColor(0);// 质量很好
			else if (ar > 256) color = colorbar.getColor(0);// 质量很差
			else color = colorbar.getColor((std::log(ar) / std::log(2)) / 8);// 质量一般
		}
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& vex = tetMesh->vertices(vh);
			/*if (vex.x() < draw_min_x || vex.y() < draw_min_y || vex.z() < draw_min_z) {
				drawFlag = false;
				break;
			}*/
		}
		if (drawFlag) {
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				auto f_vhs = tetMesh->faces(fh).getVertexHandle();
				for (int i = 2; i < f_vhs.size(); ++i) {
					std::vector<VH> trivhs = { f_vhs[0], f_vhs[i - 1], f_vhs[i] };
					for (auto& f_vh : trivhs) {
						auto& f_v = tetMesh->vertices(f_vh);
						current_drawed[f_vh] = true;
						vertices_data.push_back(f_v.x());
						vertices_data.push_back(f_v.y());
						vertices_data.push_back(f_v.z());
						vertices_data.push_back(color[0]);
						vertices_data.push_back(color[1]);
						vertices_data.push_back(color[2]);
					}
				}

			}
		} else {// 被截掉，放入透明绘制数据
			if (!tetMesh->isOnBoundary(cp.first)) continue;
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				if (!tetMesh->isOnBoundary(fh)) continue;
				auto f_vhs = tetMesh->faces(fh).getVertexHandle();
				for (int i = 2; i < f_vhs.size(); ++i) {
					std::vector<VH> trivhs = { f_vhs[0], f_vhs[i - 1], f_vhs[i] };
					for (auto& f_vh : trivhs) {
						auto& f_v = tetMesh->vertices(f_vh);
						translucence_vertices_data.push_back(f_v.x());
						translucence_vertices_data.push_back(f_v.y());
						translucence_vertices_data.push_back(f_v.z());
						translucence_vertices_data.push_back(0.7725);
						translucence_vertices_data.push_back(0.7725);
						translucence_vertices_data.push_back(0.7725);
					}
				}
			}

		}

	}

	if (ar_max > 1000) std::cout << "\n\nWorning: There exist too bad tetrahedron!!!\n\n";

	tri_cnt = vertices_data.size() / 18;
	translucence_tricnt = translucence_vertices_data.size() / 18;

	vertices_data.insert(vertices_data.end(), translucence_vertices_data.begin(), translucence_vertices_data.end());

	update();

}


void ModelDraw::updateDrawData_SurfaceMesh() 
{
	if (!triMesh) return;
	vertices_data.clear();
	vertices_data.reserve(triMesh->FaceSize() * 18); // 可用于改变vector的大小
	triMesh->genAllVerticesNormal(); // 将顶点的法向计算出来

	for (auto& fp : triMesh->allfaces()) 
	{
		auto vhs = fp.second.getVertexHandle();
		auto face_normal = fp.second.getNormal();
		Vector3d color;
		if (show_mesh_quality) color = colorbar.getColor(1 - MeshMath::get_min_angle_quality_SurfaceMesh(*triMesh, fp.first));
		if (vhs.size() == 3) 
		{
			for (auto& vh : vhs) 
			{
				auto& v = triMesh->vertices(vh);
				if (!show_mesh_quality) color = Vector3d(v.getNormalX(), v.getNormalY(), v.getNormalZ());
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}
		}

	}

	tri_cnt = vertices_data.size() / 18;

	update();

}

void ModelDraw::TestupdateDrawData_SurfaceMesh()
{
	if (!triMesh) return;
	Testvertices_data.clear();
	Testvertices_data.reserve(triMesh->FaceSize() * 18); // 可用于改变vector的大小
	triMesh->genAllVerticesNormal(); // 将顶点的法向计算出来

	for (auto& fp : triMesh->allfaces())
	{
		auto vhs = fp.second.getVertexHandle();
		auto face_normal = fp.second.getNormal();
		Vector3d color;
		if (show_mesh_quality) color = colorbar.getColor(1 - MeshMath::get_min_angle_quality_SurfaceMesh(*triMesh, fp.first));
		if (vhs.size() == 3)
		{
			for (auto& vh : vhs)
			{
				auto& v = triMesh->vertices(vh);
				if (!show_mesh_quality) color = Vector3d(v.getNormalX(), v.getNormalY(), v.getNormalZ());
				Testvertices_data.push_back(v.x());
				Testvertices_data.push_back(v.y());
				Testvertices_data.push_back(v.z());
				Testvertices_data.push_back(color[0]);
				Testvertices_data.push_back(color[1]);
				Testvertices_data.push_back(color[2]);
			}
		}

	}

	tri_cnt = Testvertices_data.size() / 18;

	update();
}

//void ModelDraw::updateDrawDate_Skel()
//{
//	MeshKernel::CurveSkel* curve_skel_ = nullptr;
//	if (proDraw->skel_tool_) curve_skel_ = proDraw->skel_tool_->GetProcessSkel();
//	if (!curve_skel_) return;
//	
//	vertices_data.clear();
//
//	VAO = 0;
//	VBO = 0;
//	EBO = 0;
//
//	if (VAO == 0)
//	{
//		core->glGenVertexArrays(1, &VAO);
//
//		//unsigned int vbo, ebo;
//		core->glGenBuffers(1, &VBO);
//		core->glGenBuffers(1, &EBO);
//
//		std::vector<glm::vec3> positions;
//		std::vector<glm::vec2> uv;
//		std::vector<glm::vec3> normals;
//		std::vector<unsigned int> indices;
//
//		const unsigned int X_SEGMENTS = 64;
//		const unsigned int Y_SEGMENTS = 64;
//		//float PI = 3.14159265359;
//		for (unsigned int y = 0; y <= Y_SEGMENTS; ++y)
//		{
//			for (unsigned int x = 0; x <= X_SEGMENTS; ++x)
//			{
//				float xSegment = (float)x / (float)X_SEGMENTS;
//				float ySegment = (float)y / (float)Y_SEGMENTS;
//
//				//这里涉及到函数的变量
//				float xPos = radius * std::cos(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin_x;
//				float yPos = radius * std::cos(ySegment * PI) + Origin_y;
//				float zPos = radius * std::sin(xSegment * 2.0f * PI) * std::sin(ySegment * PI) + Origin_z;
//
//				positions.push_back(glm::vec3(xPos, yPos, zPos));
//				uv.push_back(glm::vec2(xSegment, ySegment));
//				normals.push_back(glm::vec3(xPos, yPos, zPos));
//			}
//		}
//	}
//}

void ModelDraw::updateDrawData(bool need_update_bvhtree) {

	vertices_data.clear();
	quad_cnt = tri_cnt = line_cnt = translucence_tricnt = 0;

	if (cutFlag)
	{
		updateDrawData_CutPlane();
		return;
	}
	switch (drawmodule) {

		case DrawModule::EMPTY:
			return;

		case DrawModule::SURFACEMESH:
			updateDrawData_SurfaceMesh();
			break;
		case DrawModule::HEXMESH:
			updateDrawData_HexMesh();
			break;
		case DrawModule::TETMESH:
			updateDrawData_TetMesh();
			break;
		case DrawModule::QUADMESH:
			updateDrawData_QuadMesh();
			break;
		case DrawModule::SKEL:
			updateDrawData_HexMesh();
			break;
		case DrawModule::BSPLINESURFACE:
			updateDrawData_Surface();
			break;
			
		case Skeleton_edit_move:
			/*updateDrawDate_Skel();*/
			break;

		default:
			break;
	}
	if (need_update_bvhtree) updateBVH();

}

void ModelDraw::updateDrawData_CutPlane() {
	switch (drawmodule) {

	case DrawModule::EMPTY:
		return;
	case DrawModule::SURFACEMESH:
		//updateDrawData_SurfaceMesh_CutPlane();
		break;

	case DrawModule::TETMESH:
		updateDrawData_TetMesh_CutPlane();
		break;

	case DrawModule::HEXMESH:
		updateDrawData_HexMesh_CutPlane();
		break;

	
	default:
		break;
	}
}

bool ModelDraw::updateAllHandles() {
	switch (drawmodule)
	{
		case ViewOpenGL::EMPTY:
			break;
		case ViewOpenGL::SURFACEMESH:
			if (this->triMesh) {
				updateAllHandles_kepp_constraint(*triMesh, force_id, spc_id, constrained_point3);
				removeSelectedAll();
				updateDrawData();
				return true;
			}
			break;
		case ViewOpenGL::HEXMESH:
			if (this->hexMesh) {
				updateAllHandles_kepp_constraint(*hexMesh, force_id, spc_id, constrained_point3);
				removeSelectedAll();
				updateDrawData();
				return true;
			}
			break;
		case ViewOpenGL::QUADMESH:
			if (this->quadMesh) {
				updateAllHandles_kepp_constraint(*quadMesh, force_id, spc_id, constrained_point3);
				removeSelectedAll();
				updateDrawData();
				return true;
			}
		case ViewOpenGL::TETMESH:
			if (this->tetMesh) {
				updateAllHandles_kepp_constraint(*tetMesh, force_id, spc_id, constrained_point3);
				removeSelectedAll();
				updateDrawData();
				return true;
			}
			break;
		case ViewOpenGL::HYBRIDMESH:
			break;
		case ViewOpenGL::BSPLINECURVE:
			break;
		case ViewOpenGL::BSPLINESURFACE:
			break;
		case ViewOpenGL::MMC_Convert:
			break;
		case ViewOpenGL::TRIMESH:
			break;
		default:
			break;
	}
	QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Support iGame-Model Only! "));
	return false;
}

void ModelDraw::updateBVH() {

	//bvh_tree->destoryBVH_Tree();
	
	switch (drawmodule)
	{
		case ViewOpenGL::EMPTY:
			break;
		case ViewOpenGL::SURFACEMESH:
			bvh_tree->buildBVH_Tree(*triMesh);
			break;
		case ViewOpenGL::HEXMESH:
			bvh_tree->buildBVH_Tree(*hexMesh);
			break;
		case ViewOpenGL::TETMESH:
			bvh_tree->buildBVH_Tree(*tetMesh);
			break;
		case ViewOpenGL::QUADMESH:
			bvh_tree->buildBVH_Tree(*quadMesh);
			break;
		case ViewOpenGL::HYBRIDMESH:
			break;
		case ViewOpenGL::BSPLINECURVE:
			break;
		case ViewOpenGL::BSPLINESURFACE:
			break;
		case ViewOpenGL::SKELMESH:
			break;
		case ViewOpenGL::SKEL:
			bvh_tree->buildBVH_Tree(*hexMesh);
			break;
		case ViewOpenGL::MMC_Convert:
			break;
		case ViewOpenGL::TRIMESH:
			break;
		default:
			break;
	}
	std::cout << "Build BVH success.\n";
}

int ModelDraw::delete_islands() {
	int res = -1;
	std::vector<int> feature_vhs;
	for (auto& vp : force_id) feature_vhs.push_back(vp.first);
	for (auto& vp : spc_id) feature_vhs.push_back(vp.first);

	if (hexMesh) {
		res = MeshDeleter::delete_island_VolumeMesh(hexMesh);
	} else if (quadMesh) {
		res = MeshDeleter::delete_island_SurfaceMesh(quadMesh);
	} else if (tetMesh) {
		res = MeshDeleter::delete_island_VolumeMesh(tetMesh);
	} else if (triMesh) {
		res = MeshDeleter::delete_island_SurfaceMesh(triMesh);
	}

	bool feature_deleted = false;
	for (auto& vid : feature_vhs) {
		if (feature_deleted) break;
		VH vh(vid);
		if (hexMesh) {
			feature_deleted = hexMesh->NeighborCh(vh).empty();
		} else if (triMesh) {
			feature_deleted = triMesh->NeighborFh(vh).empty();
		} else if (tetMesh) {
			feature_deleted = tetMesh->NeighborCh(vh).empty();
		} else if (quadMesh) {
			feature_deleted = quadMesh->NeighborFh(vh).empty();
		}
	}
	if (feature_deleted) {
		QMessageBox::warning(this, QString::fromLocal8Bit(" Error!!! "), QString::fromLocal8Bit("Force/spc id were deleted!"));
	}
	updateAllHandles();
	updateDrawData();

	return res;
}

int ModelDraw::delete_redundancy_elements() {

	std::vector<int> feature_vhs;
	for (auto& vp : force_id) feature_vhs.push_back(vp.first);
	for (auto& vp : spc_id) feature_vhs.push_back(vp.first);

	int res = -1;
	if (hexMesh) {
		res = MeshDeleter::delete_redundancy_VolumeMesh(hexMesh);
	} else if (quadMesh) {
		res = MeshDeleter::delete_redundancy_SurfaceMesh(quadMesh);
	} else if (tetMesh) {
		res = MeshDeleter::delete_redundancy_VolumeMesh(tetMesh);
	} else if (triMesh) {
		res = MeshDeleter::delete_redundancy_SurfaceMesh(triMesh);
	}
	bool feature_deleted = false;
	for (auto& vid : feature_vhs) {
		if (feature_deleted) break;
		VH vh(vid);
		if (hexMesh) {
			feature_deleted = hexMesh->NeighborCh(vh).empty();
		} else if (triMesh) {
			feature_deleted = triMesh->NeighborFh(vh).empty();
		} else if (tetMesh) {
			feature_deleted = tetMesh->NeighborCh(vh).empty();
		} else if (quadMesh) {
			feature_deleted = quadMesh->NeighborFh(vh).empty();
		}
	}
	if (feature_deleted) {
		QMessageBox::warning(this, QString::fromLocal8Bit(" Error!!! "), QString::fromLocal8Bit("Force/spc id were deleted!"));
	}
	updateAllHandles();
	updateDrawData();

	return res;
}


int ModelDraw::getVerticesSize() {
	switch (drawmodule)
	{
		case ViewOpenGL::EMPTY:
			return 0;
			break;
		case ViewOpenGL::SURFACEMESH:
			return triMesh->VertexSize();
			break;
		case ViewOpenGL::HEXMESH:
			return hexMesh->vsize();
			break;
		case ViewOpenGL::TETMESH:
			break;
		case ViewOpenGL::HYBRIDMESH:
			break;
		case ViewOpenGL::BSPLINECURVE:
			break;
		case ViewOpenGL::BSPLINESURFACE:
			break;
		case ViewOpenGL::SKELMESH:
			break;
		case ViewOpenGL::SKEL:
			break;
		case ViewOpenGL::MMC_Convert:
			break;
		case ViewOpenGL::TRIMESH:
			break;
		default:
			break;
	}
	return -1;
}


void ModelDraw::compared_bind_vertices() {
	if (hexMesh == nullptr) return;

	compared_to_current.clear();
	int cmp_vcnt = compared_vertices_data.size() / 6;
	for (int i = 0; i < cmp_vcnt; ++i) {
		int min_vh = -1, idx = i * 6;
		double min_dis = 999999999;
		MeshKernel::iGameVertex cmp_vex(compared_vertices_data[idx], compared_vertices_data[idx+1], compared_vertices_data[idx+2]);
		for (auto& vp : hexMesh->allvertices()) {
			double dis = (cmp_vex - vp.second).norm();
			if (dis < min_dis) {
				min_dis = dis;
				min_vh = vp.first;
				if (min_dis < 1E-6F) break;
			}
		}
		compared_to_current[i] = min_vh;
	}

}

void ModelDraw::read_constraint_force(string str) {
	force_id.clear();
	read_constraint_data(str, force_id);
	for (auto& vid : force_id) {
		VH vh(vid.first);
		Vex origin_pos = get_point_via_vh(vh);
		constrained_point3[vh] = origin_pos;
	}
}

void ModelDraw::read_constraint_spc(string str) {
	spc_id.clear();
	read_constraint_data(str, spc_id);
	for (auto& vid : spc_id) {
		VH vh(vid.first);
		Vex origin_pos = get_point_via_vh(vh);
		constrained_point3[vh] = origin_pos;
	}
}

void ModelDraw::write_constraint_force(string str) {
	write_constraint_data(str, force_id);
}

void ModelDraw::write_constraint_spc(string str) {
	write_constraint_data(str, spc_id);
}

void ModelDraw::remove_all_constraints() {
	force_id.clear();
	spc_id.clear();
}

void ModelDraw::smoothing(int iter, bool flat_surface_preserve) {

	time_t beg = clock();

	if (hexMesh || tetMesh) {
		vector<int> feature_chs = selected_c_indices;
		vector<int> feature_fhs = selected_f_indices;
		vector<int> feature_ehs = selected_e_indices;
		vector<int> feature_vhs = selected_v_indices;
		for (auto& vp : force_id) feature_vhs.push_back(vp.first);
		for (auto& vp : spc_id) feature_vhs.push_back(vp.first);

		if (hexMesh) 
		{
			VolumeMesh_Smoothing app(*hexMesh);
			app.Laplacian_Smoothing(iter, false, feature_chs, feature_fhs, feature_ehs, feature_vhs);
		} 
		else if (tetMesh) 
		{
			/*VolumeMesh_Smoothing app(*tetMesh);
			app.Laplacian_Smoothing(iter, feature_chs, feature_fhs, feature_ehs, feature_vhs);*/
			triMesh = new MeshKernel::SurfaceMesh();

			std::unordered_map<VH, VH> tetvh_to_trivh;
			auto backup_force_id = force_id;
			auto backup_spc_id = spc_id;
			auto backup_points = constrained_point3;

			MeshConvertor::get_trimesh_from_tetmesh(tetMesh, triMesh, force_id, spc_id, constrained_point3, tetvh_to_trivh);

			assert(triMesh);
			SetDrawModule(SURFACEMESH);


			std::vector<VH> feature_vhs_processing;// 特征顶点
			for (auto& vid : spc_id) feature_vhs_processing.push_back(VH(vid.first));
			for (auto& vid : force_id) feature_vhs_processing.push_back(VH(vid.first));

			for (int it = 0; it < iter; ++it) {
				Denoising_BNF app(*triMesh, feature_vhs_processing, flat_surface_preserve);
				app.Execute();
			}

			QMessageBox messageBox1(QMessageBox::NoIcon,
				" Poster ", "Denoising success. Remeshing continue?\nIf not, I recommend you save the file.",
				QMessageBox::Yes | QMessageBox::No, NULL); ;
			int result1 = messageBox1.exec();
			if (result1 == QMessageBox::No) {
				updateDrawData();
				return;
			}
			bool input_ok = false;
			double edge_length_ratio = QInputDialog::getDouble(this, " Remeshing ", "Please input edge length ratio",
				1.f, 0.1f, 2.f, 2, &input_ok);
			if (!input_ok) return;
			Remeshing app(*triMesh, true, feature_vhs_processing, edge_length_ratio);
			app.Execute();
			updateAllHandles_kepp_constraint(*triMesh, force_id, spc_id, constrained_point3);

			delete tetMesh;
			tetMesh = nullptr;

		}

	} else if (triMesh || quadMesh) {
		vector<int> feature_fhs = selected_f_indices;
		vector<int> feature_ehs = selected_e_indices;
		vector<int> feature_vhs = selected_v_indices;
		for (auto& vp : force_id) feature_vhs.push_back(vp.first);
		for (auto& vp : spc_id) feature_vhs.push_back(vp.first);

		if (triMesh) {// BNF 去噪
			std::vector<VH> feature_vhs_processing;
			for (auto& vid : spc_id) feature_vhs_processing.push_back(VH(vid.first));
			for (auto& vid : force_id) feature_vhs_processing.push_back(VH(vid.first));
			for (int it = 0; it < iter; ++it) {
				Denoising_BNF app(*triMesh, feature_vhs_processing, flat_surface_preserve);
				app.Execute();
			}
		} else if (quadMesh) {
			SurfaceMesh_Smoothing app(*quadMesh);
			app.Laplacian_Smoothing(iter, feature_fhs, feature_ehs, feature_vhs);
		}

	}
	int cost = clock() - beg;
	updateDrawData();
	update();

	QString str = "Smoothing success. cost " + QString::number(cost) + "ms.";
	QMessageBox::information(this, QString::fromLocal8Bit(" Message "), str);
	

}

// 建议修改（以后）：将HexMesh存在solutionSkel里比较好
int ModelDraw::findNearestCell(double aim[3]) 
{
	iGameVertex aimV = iGameVertex(aim[0], aim[1], aim[2]);
	double minDistance = 999;
	iGameCellHandle record;
	for (auto& cp : hexMesh->allcells()) 
	{		// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		iGameVertex ccenter = hexMesh->getCellCenter(cp.first);
		double dist = ccenter.dist(aimV);
		if (dist < minDistance) 
		{
			dist = minDistance;
			record = cp.first;
		}
	}
	return record;
}

// 延时函数
void ModelDraw::Sleep(int msec) {
	QTime dieTime = QTime::currentTime().addMSecs(msec);
	while (QTime::currentTime() < dieTime)
		QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}

std::string ModelDraw::convertMeshToOffString() {
	if (!triMesh)	return "";
	MeshKernel::IO io;
	return io.WriteOffString(*triMesh);
}

bool ModelDraw::loadTetMeshFromStr(const std::string& data) {
	clock_t time1 = clock();
	MeshKernel::IO io;

	std::stringstream ss(data);
	std::vector<iGameVertex> vertices;
	std::vector<std::vector<iGameVertexHandle>> surface_faces;
	std::vector<std::vector<iGameVertexHandle>> cells;
	std::vector<iGameVertexHandle> cell;
	std::string str;
	enum State
	{
		USELESS = 1, VERTEX, FACE, TET, EDGE
	}state = USELESS;
	while (std::getline(ss, str)) {
		int len = str.length();
		std::vector<std::string>info;
		std::string s;
		for (int i = 0; i < len; i++) {
			if (str[i] == '#')break;
			if (str[i] == ' ') {
				if (s.length() > 0)
					info.push_back(s);
				s = "";
			} else {
				s.push_back(str[i]);
			}
		}
		if (s.length() > 0)
			info.push_back(s);
		if (info.size() == 0)continue;
		if (info[0] == "MeshVersionFormatted") {
			state = USELESS;
		} else if (info[0] == "Dimension") {
			std::getline(ss, str);
			state = USELESS;
		} else if (info[0] == "Vertices") {
			std::getline(ss, str);
			state = VERTEX;
		} else if (info[0] == "Tetrahedra") {
			std::getline(ss, str);
			state = TET;
		} else if (info[0] == "Triangles") {
			std::getline(ss, str);
			state = FACE;
		} else if (info[0] == "Edges") {
			std::getline(ss, str);
			state = EDGE;
		} else if (info[0] == "End") {
			state = USELESS;
		} else {
			if (state == USELESS) {
				continue;
			} else if (state == VERTEX) {
				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
					, std::stod(info[1]), std::stod(info[2])));
			} else if (state == TET) {
				cell = std::vector<MeshKernel::iGameVertexHandle>{
						(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
				};
				cells.push_back(cell);
				cell.clear();
			} else if (state == FACE) {
				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
					(iGameVertexHandle)(std::stoi(info[0]) - 1)
						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
				});
			}
		}
	}

	tetMesh = new MeshKernel::VolumeMesh(vertices, cells);
	SetDrawModule(TETMESH);

	clock_t time2 = clock();
	std::cout << "TetGen success. cost " << time2 - time1 << " ms" << std::endl;

	update_constrain_nid_after_tetgen(*tetMesh, force_id, spc_id, constrained_point3);
	updateDrawData();
	return true;
}

void ModelDraw::get_mesh_quality(string& str) {

	double det_min = 2, det_max = -1, det_avg = 0;
	double ar_min = 999999, ar_max = 0, ar_avg = 0;
	std::string str0, str1, str2;

	switch (drawmodule)
	{
		case ViewOpenGL::EMPTY:
			break;
		case ViewOpenGL::SURFACEMESH:
			break;
		case ViewOpenGL::HEXMESH:

			for (auto& cp : hexMesh->allcells()) {
				double det = MeshMath::get_min_scaled_Jacobian_HexMesh(*hexMesh, cp.first);
				det_min = std::min(det_min, det);
				det_max = std::max(det_max, det);
				det_avg += det;
			}
			det_avg /= hexMesh->csize();
			str0 = std::to_string(det_min);
			str0 = str0.substr(0, std::find(str0.begin(), str0.end(), '.') - str0.begin() + 3);
			if (det_min < 0.2) str0 += " ( Bad! )";
			str1 = std::to_string(det_max);
			str1 = str1.substr(0, std::find(str1.begin(), str1.end(), '.') - str1.begin() + 3);
			str2 = std::to_string(det_avg);
			str2 = str2.substr(0, std::find(str2.begin(), str2.end(), '.') - str2.begin() + 3);
			str = "HexMesh-Scaled Jacobian:\n\tmin\t=\t" + str0 +
				"\t,\n\tmax\t=\t" + str1 + "\t,\n\tavg\t=\t" + str2 + "\t.\n(A high-scaled Jacobian implies a better shape.)";
			break;
		case ViewOpenGL::TETMESH:

			for (auto& cp : tetMesh->allcells()) {
				double ar = MeshMath::get_aspect_ratio_TetMesh(*tetMesh, cp.first);
				ar_min = std::min(ar_min, ar);
				ar_max = std::max(ar_max, ar);
				ar_avg += ar;
			}
			ar_avg /= tetMesh->csize();
			str0 = std::to_string(ar_min);
			str0 = str0.substr(0, std::find(str0.begin(), str0.end(), '.') - str0.begin() + 3);
			str1 = std::to_string(ar_max);
			str1 = str1.substr(0, std::find(str1.begin(), str1.end(), '.') - str1.begin() + 3);
			if (ar_max > 128) str1 += " ( Bad! )";
			str2 = std::to_string(ar_avg);
			str2 = str2.substr(0, std::find(str2.begin(), str2.end(), '.') - str2.begin() + 3);
			str = "TetMesh-Aspect Ratio:\n\tmin\t=\t" + str0 +
				"\t,\n\tmax\t=\t" + str1 + "\t,\n\tavg\t=\t" + str2 + "\t.\n(A low-aspect ratio implies a better shape.)";
			break;
		case ViewOpenGL::QUADMESH:
			break;
		case ViewOpenGL::HYBRIDMESH:
			break;
		case ViewOpenGL::BSPLINECURVE:
			break;
		case ViewOpenGL::BSPLINESURFACE:
			break;
		case ViewOpenGL::SKELMESH:
			break;
		case ViewOpenGL::SKEL:
			break;
		case ViewOpenGL::MMC_Convert:
			break;
		case ViewOpenGL::TRIMESH:
			break;
		default:
			break;
	}

}

void ModelDraw::save_as_iges(string filename) {

	if (!this->hexMesh) {
		QMessageBox::warning(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Support Hex Mesh Only! "));
		return;
	}

	//QProgressDialog* dialog = new QProgressDialog(QString::fromLocal8Bit("Converting..."), QString::fromLocal8Bit("Cancel"), 0, 12, this);
	//dialog->setMinimumSize(QSize(200, 120));
	//dialog->setWindowTitle(QString::fromLocal8Bit("Save as IGES"));
	//dialog->setWindowModality(Qt::WindowModal);
	//QCoreApplication::processEvents();
	//dialog->setMinimumDuration(1);
	//dialog->show();
	//dialog->setValue(0);
	//dialog->update();

	//HexV3fMesh* _hexMesh = new HexV3fMesh();
	//dialog->setValue(1);

	///* iGame_hexMesh to XIE_hexMesh */
	//{
	//	uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();
	//	for (uint64_t i = 0; i < vcnt; ++i) {
	//		MeshKernel::iGameVertexHandle vh(i);
	//		auto& v = hexMesh->vertices(vh);
	//		V3f _v(v.x(), v.y(), v.z());
	//		_hexMesh->add_vertex(_v);
	//	}
	//	dialog->setValue(2);
	//	for (uint64_t i = 0; i < ccnt; ++i) {
	//		MeshKernel::iGameCellHandle ch(i);
	//		if (hexMesh->isValid(ch)) {
	//			auto& c = hexMesh->cells(ch);
	//			auto vhs = c.getVertexHandle();
	//			vector<VertexHandle> _vhs;
	//			for (auto& vh : vhs) {
	//				_vhs.push_back(VertexHandle(vh));
	//			}
	//			_hexMesh->add_cell(_vhs);
	//		}
	//	}
	//	dialog->setValue(3);
	//}


	//if (_hexMesh->vertices_.size() == 0) {
	//	dialog->close();
	//	QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
	//	return;
	//}

	//_hexMesh->examine_bdy();
	//dialog->setValue(4);
	//_hexMesh->cal_cell_cen();
	//_hexMesh->cal_cen();
	//dialog->setValue(5);
	//_hexMesh->examine_cs();
	//dialog->setValue(6);
	////_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	////_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	//std::cout << "Converting... hexMesh: vertices = " << _hexMesh->vertices_.size() << ", cells = " << _hexMesh->cells_.size() << std::endl;
	//int udegree = 3;
	//int vdegree = 3;
	//vector<double> uknots{ 0,1 };
	//vector<double> vknots{ 0,1 };
	//vector<int> umultis{ 4,4 };
	//vector<int> vmultis{ 4,4 };

	//HexV3fMesh* SkelMesh = _hexMesh;
	//std::vector<EdgeHandle> TempEdge;
	//ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

	//int idx = 0;
	//AppCC.approximationing();
	//dialog->setValue(7);
	//dialog->setValue(8);
	//dialog->setValue(9);
	//dialog->setValue(10);

	//vector<CBSplineSurface> MutiBezierPatch;
	//for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
	//	if (SkelMesh->is_bdy((FaceHandle)i)) {
	//		vector<vector<CPoint>> controlPointsSurface;
	//		for (int j = 0; j < 4; ++j) {
	//			vector<CPoint> controlPointLine;
	//			for (int k = 0; k < 4; ++k) {
	//				controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
	//			}
	//			controlPointsSurface.push_back(controlPointLine);
	//		}
	//		CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
	//		MutiBezierPatch.push_back(splineSurface);
	//		idx++;
	//	}
	//}

	//CAidFunction::getIGESSurface(filename, MutiBezierPatch);
	//dialog->setValue(12);
	//dialog->close();

}

Vex ModelDraw::get_point_via_vh(const VH& vh) {
	Vex res(0, 0, 0);
	switch (drawmodule) {
		case ViewOpenGL::EMPTY:
			break;
		case ViewOpenGL::SURFACEMESH:
			if (triMesh && triMesh->isValid(vh)) {
				res = triMesh->vertices(vh);
			}
			break;
		case ViewOpenGL::HEXMESH:
			if (hexMesh && hexMesh->isValid(vh)) {
				res = hexMesh->vertices(vh);
			}
			break;
		case ViewOpenGL::TETMESH:
			if (tetMesh && tetMesh->isValid(vh)) {
				res = tetMesh->vertices(vh);
			}
			break;
		case ViewOpenGL::QUADMESH:
			if (quadMesh && quadMesh->isValid(vh)) {
				res = quadMesh->vertices(vh);
			}
			break;
		case ViewOpenGL::HYBRIDMESH:
			break;
		case ViewOpenGL::BSPLINECURVE:
			break;
		case ViewOpenGL::BSPLINESURFACE:
			break;
		case ViewOpenGL::SKELMESH:
			break;
		case ViewOpenGL::SKEL:
			break;
		case ViewOpenGL::MMC_Convert:
			break;
		case ViewOpenGL::TRIMESH:
			break;
		default:
			break;
	}

	return res;

}


//Xie CAD 重建  
bool ModelDraw::hexMeshConvertToIGES() 
{
	if (!this->hexMesh) return false;

	//std::string filename = "tmp_convert.mesh";

	//MeshKernel::IO io;
	//io.iGameWriteVolumeFile(*hexMesh, filename);
	HexV3fMesh* _hexMesh = new HexV3fMesh();

	/* iGame_hexMesh to XIE_hexMesh */
	uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();

	for (uint64_t i = 0; i < vcnt; ++i)
	{
		MeshKernel::iGameVertexHandle vh(i);
		auto& v = hexMesh->vertices(vh);
		V3f _v(v.x(), v.y(), v.z());
		_hexMesh->add_vertex(_v);
	}

	for (uint64_t i = 0; i < ccnt; ++i) 
	{
		MeshKernel::iGameCellHandle ch(i);
		if (hexMesh->isValid(ch))
		{
			auto& c = hexMesh->cells(ch);
			auto vhs = c.getVertexHandle();
			vector<VertexHandle> _vhs;
			for (auto& vh : vhs) 
			{
				_vhs.push_back(VertexHandle(vh));
			}
			_hexMesh->add_cell(_vhs);
		}
	}

	/*FileHandler file(filename);
	assert(file.read2HexMeshObj(_hexMesh, true, true, true));*/

	if (_hexMesh->vertices_.size() == 0) 
	{
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
		return false;
	}

	_hexMesh->examine_bdy();  //把边界点记录下来，记录在v_bdy
	_hexMesh->cal_cell_cen(); //计算六面体点的体中心，记录在cell_r
	_hexMesh->cal_cen();      //计算一些bbox 中心点 半径什么的
	_hexMesh->examine_cs();   //将折痕点全部初始化为false 没有计算
	//_hexMesh->smooth_mesh();
	//_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	//_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	std::cout << "converting... hexMesh: vertices = " << _hexMesh->vertices_.size() << ", cells = " << _hexMesh->cells_.size() << std::endl;
	
	int udegree = 3;
	int vdegree = 3;
	vector<double> uknots{ 0,1 };
	vector<double> vknots{ 0,1 };
	vector<int> umultis{ 4,4 };
	vector<int> vmultis{ 4,4 };

	HexV3fMesh* SkelMesh = _hexMesh;
	std::vector<EdgeHandle> TempEdge;
	ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

	int idx = 0;
	AppCC.approximationing();

	vector<CBSplineSurface> MutiBezierPatch;
	for (int i = 0; i < SkelMesh->faces_.size(); ++i) 
	{
		if (SkelMesh->is_bdy((FaceHandle)i)) 
		{
			vector<vector<CPoint>> controlPointsSurface;
			for (int j = 0; j < 4; ++j) 
			{
				vector<CPoint> controlPointLine;
				for (int k = 0; k < 4; ++k) 
				{
					controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
				}
				controlPointsSurface.push_back(controlPointLine);
			}
			CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
			MutiBezierPatch.push_back(splineSurface);
			idx++; // 用于记录当前处理的面
		}
	}
	bsplineSurfaces = MutiBezierPatch;
	saveIGESFile("tmp_convert.iges");
	//测试
	loadIGESSurface(MutiBezierPatch);
	//loadIGESSurface(MutiBezierPatch, *SphereMesh_Origin); //用于比较输出误差
	cout << "Convert success.\n";
	return true;
}

bool ModelDraw::loadIGESSurface(std::string fileName) {
	vector<CBSplineSurface>geomsurface;
	if (CAidFunction::readIGESGeomSurface(fileName, geomsurface)) 
	{
		if (geomsurface.empty()) {
			return false;
		}
		bsplineSurfaces.clear();
		loadIGESSurface(geomsurface);
		return true;
	}
	return false;
}

bool ModelDraw::saveIGESFile(std::string fileName) 
{
	if (!bsplineCurves.empty()) 
	{
		CAidFunction::getIGESCurve(fileName, bsplineCurves);
	}
	else if (!bsplineSurfaces.empty()) 
	{
		CAidFunction::getIGESSurface(fileName, bsplineSurfaces);
	}
	else {
		QMessageBox::information(this, QString::fromLocal8Bit(" 警告 "), QString::fromLocal8Bit(" 没有保存任何的模型 "));
		return false;
	}
	return true;
}

bool ModelDraw::saveBSplineToStep(std::string fileName)
{
	if (!bsplineSurfaces.empty())
	{
		CAidFunction::SaveStepSurface(fileName, bsplineSurfaces);
	}
	else {
		QMessageBox::information(this, QString::fromLocal8Bit(" 警告 "), QString::fromLocal8Bit(" 没有保存任何的模型 "));
		return false;
	}
	return true;
}

bool ModelDraw::QuadMeshToSpline() {
	if (drawmodule != QUADMESH) return false;
	// 检测是否为纯四边形
	bool notQuad = false;
	for (auto fmap : quadMesh->allfaces()) {
		if (fmap.second.size() != 4) {
			notQuad = true;
			selected_f_indices.push_back(fmap.first);

			glBegin(GL_QUADS);
			for (auto& f_idx : selected_f_indices) 
			{
				auto vhs = quadMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
				for (auto& vh : vhs) 
				{
					auto& v = quadMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
			glEnd();

			updateDrawData();
			//break;
		}
	}
	if (notQuad) {
		std::cout << "this mesh is not quad! break" << std::endl;
		return false;
	}

	vector<vector<vector<Vec>>> bezierPatchCtrlPoints;

	CAD_Reconstruction alg;
	alg.get_bezier_control_mesh(*quadMesh, bezierPatchCtrlPoints);

	int udegree = 3;
	int vdegree = 3;
	vector<double> uknots{ 0,1 };
	vector<double> vknots{ 0,1 };
	vector<int> umultis{ 4,4 };
	vector<int> vmultis{ 4,4 };
	vector<CBSplineSurface> MutiBezierPatch;
	for (int i = 0; i < bezierPatchCtrlPoints.size(); ++i) {
		assert(bezierPatchCtrlPoints[i].size() == 4);
		vector<vector<CPoint>> controlPointsSurface;
		for (int j = 0; j < 4; ++j) {
			assert(bezierPatchCtrlPoints[i][j].size() == 4);
			vector<CPoint> controlPointLine;
			for (int k = 0; k < 4; ++k) {
				controlPointLine.push_back(CPoint(bezierPatchCtrlPoints[i][j][k][0],
					bezierPatchCtrlPoints[i][j][k][1],
					bezierPatchCtrlPoints[i][j][k][2]));
			}
			controlPointsSurface.push_back(controlPointLine);
		}
		CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
		MutiBezierPatch.push_back(splineSurface);
	}

	// 建立新模型
	SetDrawModule(BSPLINESURFACE);

	/*createModel();
	currentModel->structure = BSPLINE;*/

	if (MutiBezierPatch.empty()) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}
	//bSplineSurfaces = geomSurface;
	//bsplineSurfaceToDiscrete();
	////controlpointsToDiscrete();
	//updateCurrentMesh();

	//vector<CBSplineSurface> MutiBezierPatch;
	/*for (int i = 0; i < SkelMesh->faces_.size(); ++i)
	{
		if (SkelMesh->is_bdy((FaceHandle)i))
		{
			vector<vector<CPoint>> controlPointsSurface;
			for (int j = 0; j < 4; ++j)
			{
				vector<CPoint> controlPointLine;
				for (int k = 0; k < 4; ++k)
				{
					controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
				}
				controlPointsSurface.push_back(controlPointLine);
			}
			CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
			MutiBezierPatch.push_back(splineSurface);
			idx++;
		}
	}*/
	bsplineSurfaces = MutiBezierPatch;
	saveIGESFile("tmp_convert.iges");
	//测试
	loadIGESSurface(MutiBezierPatch);
	//loadIGESSurface(MutiBezierPatch, *SphereMesh_Origin); //用于比较输出误差
	cout << "Convert success.\n";
	return true;

	update();

	cout << "QuadMesh convert to BSplineSurface success.\n";

	return true;
}

//void ModelDraw::bsplineSurfaceToDiscrete() {
//	int upatchnum = 10;//渲染时的片数
//	int vpatchnum = 10;//渲染时的片数
//
//	time_t beg = clock();
//	data_surfaces.clear();
//	controlMesh.clear();
//	double umin, umax, vmin, vmax;
//	vector<CPoint> singlePatch;// 每个面片的离散点
//	vector<vector<CPoint>> singlelSurface;// 每个曲面的离散点
//	vector<vector<vector<CPoint>>> allSurface;// 所有曲面集合的离散点
//
//	vector<CVector> singlePatch_normal;
//	vector<vector<CVector>> singlelSurface_normal;
//	vector<vector<vector<CVector>>> allSurface_normal;
//	CVector normal;//保存法向量
//
//	int totalVerticesSize = 0;
//	for (int i = 0; i < bSplineSurfaces.size(); i++) {
//		bSplineSurfaces[i].getUVMinMax(umin, umax, vmin, vmax);
//		int unumber, vnumber;
//		unumber = upatchnum * (bSplineSurfaces[i].getSizeUControlPoints());
//		vnumber = vpatchnum * (bSplineSurfaces[i].getSizeVControlPoints());
//		unumber = unumber < 100 ? unumber : 100;
//		vnumber = vnumber < 100 ? vnumber : 100;
//		//std::cout << "u_num: " << unumber << "; v_num: " << vnumber << std::endl;
//		//unumber = 5, vnumber = 5;// for faster rendering
//		unumber *= 0.2, vnumber *= 0.2f;
//		double sizeu = (umax - umin) / unumber - 1e-9;// 步距
//		double sizev = (vmax - vmin) / vnumber - 1e-9;
//		for (double u = umin + 1e-9; u < umax - sizeu - 1e-9; u = u + sizeu) {
//			for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v = v + sizev) {
//				singlePatch.push_back(bSplineSurfaces[i].getPointAtParam(u, v));
//				CAidFunction::geomNormal(bSplineSurfaces[i], u, v, normal);
//				singlePatch_normal.push_back(normal);
//
//				singlePatch.push_back(bSplineSurfaces[i].getPointAtParam(u, v + sizev));
//				CAidFunction::geomNormal(bSplineSurfaces[i], u, v + sizev, normal);
//				singlePatch_normal.push_back(normal);
//
//				singlePatch.push_back(bSplineSurfaces[i].getPointAtParam(u + sizeu, v + sizev));
//				CAidFunction::geomNormal(bSplineSurfaces[i], u + sizeu, v + sizev, normal);
//				singlePatch_normal.push_back(normal);
//
//				singlePatch.push_back(bSplineSurfaces[i].getPointAtParam(u + sizeu, v));
//				CAidFunction::geomNormal(bSplineSurfaces[i], u + sizeu, v, normal);
//				singlePatch_normal.push_back(normal);
//
//				singlelSurface.push_back(singlePatch);
//				singlelSurface_normal.push_back(singlePatch_normal);
//				singlePatch.clear();
//				singlePatch_normal.clear();
//			}
//
//		}
//
//		if (!singlelSurface.empty()) totalVerticesSize += singlelSurface.size() * singlelSurface[0].size();
//
//		allSurface.push_back(singlelSurface);
//		allSurface_normal.push_back(singlelSurface_normal);
//		singlelSurface.clear();
//		singlelSurface_normal.clear();
//	}
//	this->data_surfaces = allSurface;// 曲面离散点集合,data_surface[i][j][k]表示第i个曲面的第j个面片的第k个点（0<=k<=3）
//	this->surface_normals = allSurface_normal;
//	// 边界点先不用
//	//this->data_bound_surface = getData_boundcurve(bSplineSurfaces);// 曲面边界曲线离散点集合,datda_bound[i][j][k]表示第i个曲面的第j个边界（0<=j<=3)的第k个点
//	//获取控制点
//	vector<vector<CPoint>> singleSurfaceControlPoints;
//	vector<vector<vector<CPoint>>> allSurfaceControlPoints;
//	for (int i = 0; i < bSplineSurfaces.size(); i++) {
//		singleSurfaceControlPoints = bSplineSurfaces[i].getControlPoints();
//		allSurfaceControlPoints.push_back(singleSurfaceControlPoints);
//	}
//	this->controlMesh = allSurfaceControlPoints;
//	time_t end = clock();
//	printf("bsplineSurfaceToDiscrete cost: %zdms, vertices size: %zd, surface size: %zd\n", int(end - beg), totalVerticesSize, allSurface.size());
//}


bool ModelDraw::loadIGESSurface(std::vector<CBSplineSurface>& geomsurface)
{
	ModelClear();
	if (geomsurface.empty()) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}
	bsplineSurfaces = geomsurface;
	surfacesDraw = bsplineSurfaces;
	SetDrawModule(BSPLINESURFACE);
	bsplineSurfaceToDiscrete();
	controlpointsToDiscrete();
	if (resetViewFlag) updateSurfaceModel();
	updateDrawData();
	update();
	return true;
}

bool ModelDraw::loadIGESSurface(std::vector<CBSplineSurface>& geomsurface, SurfaceMesh& trimesh) {
	ModelClear();
	if (geomsurface.empty())
 {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}
	bsplineSurfaces = geomsurface;
	surfacesDraw = bsplineSurfaces;
	SetDrawModule(BSPLINESURFACE);
	bsplineSurfaceToDiscrete(); //这个地方把 贝塞尔面的绘制数据穿进去了 所以updateDrawData_Surface_IGES_Error(trimesh);我在里面加了一下clear()-
	controlpointsToDiscrete();
	if (resetViewFlag) updateSurfaceModel();
	updateDrawData_Surface_IGES_Error(trimesh);
	update();
	return true;
}

bool ModelDraw::loadSurfaceErrorIGES(std::vector<CBSplineSurface>& geomsurface, SurfaceMesh& trimesh) {
	ModelClear();
	if (geomsurface.empty()) 
	{
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}
	bsplineSurfaces = geomsurface;
	surfacesDraw = bsplineSurfaces;
	bsplineSurfaceToDiscrete();
	controlpointsToDiscrete();
	SetDrawModule(SURFACEMESH);
	triMesh = new SurfaceMesh(trimesh);
	if (resetViewFlag) updateMeshModel();
	updateDrawData_SurfaceMesh_Error_IGES();
	update();
	return true;
}

void ModelDraw::bsplineSurfaceToDiscrete() 
{
	time_t beg = clock();
	data_surfaces.clear();
	controlMesh.clear();
	double umin, umax, vmin, vmax;
	vector<CPoint> singlePatch;// 每个面片的离散点
	vector<vector<CPoint>> singlelSurface;// 每个曲面的离散点
	vector<vector<vector<CPoint>>> allSurface;// 所有曲面集合的离散点

	vector<CVector>singlePatch_normal;
	vector<vector<CVector>>singlelSurface_normal;
	vector<vector<vector<CVector>>>allSurface_normal;
	CVector normal;//保存法向量

	int totalVerticesSize = 0;
	
	//surfacesDraw 应该记录着多个贝塞尔面的函数表达形式，数据格式是来自OCC
	for (int i = 0; i < surfacesDraw.size(); i++) 
	{
		surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
		int unumber, vnumber;
		unumber = upatchnum * (surfacesDraw[i].getSizeUControlPoints()); // 渲染时的片数 * 后者是m_ControlPoints的size
		vnumber = vpatchnum * (surfacesDraw[i].getSizeVControlPoints()); // 渲染时的片数 * 后者是m_ControlPoints[0]的size
		unumber = unumber < 100 ? unumber : 100;
		vnumber = vnumber < 100 ? vnumber : 100;
		//std::cout << "u_num: " << unumber << "; v_num: " << vnumber << std::endl;
		//unumber = 5, vnumber = 5;// for faster rendering
		unumber *= 0.2, vnumber *= 0.2f;
		double sizeu = (umax - umin) / unumber - 1e-9;// 步距
		double sizev = (vmax - vmin) / vnumber - 1e-9;
		for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u = u + sizeu)  // u 和 v 应该是生成面时候的参数 生成的事贝塞尔曲面
		{
			for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v = v + sizev) 
			{
				singlePatch.push_back(surfacesDraw[i].getPointAtParam(u, v)); // 可以从贝塞尔曲线推广，对于每一个参数t，可以对应n次贝塞尔的1个点
				CAidFunction::geomNormal(surfacesDraw[i], u, v, normal); //传入B样条曲面的handle 参数 u v 得到法向normal 是OCC库中的函数
				singlePatch_normal.push_back(normal);

				singlePatch.push_back(surfacesDraw[i].getPointAtParam(u, v + sizev));
				CAidFunction::geomNormal(surfacesDraw[i], u, v + sizev, normal);
				singlePatch_normal.push_back(normal);

				singlePatch.push_back(surfacesDraw[i].getPointAtParam(u + sizeu, v + sizev));
				CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v + sizev, normal);
				singlePatch_normal.push_back(normal);

				singlePatch.push_back(surfacesDraw[i].getPointAtParam(u + sizeu, v));
				CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v, normal);
				singlePatch_normal.push_back(normal);

				singlelSurface.push_back(singlePatch); //整 个曲面包含了众多的patch
				singlelSurface_normal.push_back(singlePatch_normal); 
				singlePatch.clear(); 
				singlePatch_normal.clear();
			}

		}

		allSurface.push_back(singlelSurface); //所有曲面 和 单个曲面有什么区别呢？
		if (!singlelSurface.empty()) totalVerticesSize += singlelSurface.size() * singlelSurface[0].size(); //
		singlelSurface.clear();
		allSurface_normal.push_back(singlelSurface_normal);
		singlelSurface_normal.clear();
	}

	// data_surfaces是渲染相关的
	this->data_surfaces = allSurface;// 曲面离散点集合,data_surface[i][j][k]表示第i个曲面的第j个面片的第k个点（0<=k<=3） 第j个面片是那个四边形
	this->surface_normals = allSurface_normal;
	this->data_bound_surface = getData_boundcurve(surfacesDraw);// 曲面边界曲线离散点集合,datda_bound[i][j][k]表示第i个曲面的第j个边界（0<=j<=3)的第k个点？？//没用下面这个函数会很奇怪，就是会到原点去，感觉是改变相对坐标的

	////用于对比误差 hzk 记录生成的四边形网格
	////std::vector<iGameVertex> VerticesMeshCADForContrast; // 用于记录后面CAD模型的四边形网格顶点
	////std::vector<std::vector<iGameVertexHandle>> FacesMeshCADForContrast; // 用于记录后面CAD模型的四边形网格面
	////MeshKernel::SurfaceMesh* SphereMesh_MeshCADForContrast = nullptr;
	////PreNumHandle = 0;
	//int indexHandleForMeshCADContrast = 0;
	//int indexFaceHandleForMeshCADContrast = 0;
	//for (int i = 0; i < allSurface.size(); i++)
	//{
	//	auto& singleSurface = data_surfaces[i]; //每个贝塞尔面
	//	for (int j = 0; j < singleSurface.size(); j++) 
	//	{
	//		//每个贝塞尔面中64个四边形
	//		std::vector<iGameVertexHandle> tempface;
	//		for (int k = 0; k < singleSurface[j].size(); k++) 
	//		{//每个四边形的4个点
	//			
	//			VerticesMeshCADForContrast.push_back(iGameVertex(singleSurface[j][k].getX(), singleSurface[j][k].getY(), singleSurface[j][k].getZ()));
	//			tempface.push_back(VH(indexHandleForMeshCADContrast++));
	//			/*std::vector<iGameVertexHandle> tempface;
	//			for (int j = 0; j < 3; j++)
	//			{
	//				tempface.push_back(VH(indices[i + j]));
	//			}
	//			_faces.push_back(tempface);
	//			FacesMeshCADForContrast[indexFaceHandleForMeshCADContrast].push_back(VH(indexFaceHandleForMeshCADContrast++));*/
	//			//vertices_data.push_back(normal[i][j].getX());
	//			//vertices_data.push_back(normal[i][j].getY());
	//			//vertices_data.push_back(normal[i][j].getZ());
	//		}
	//		FacesMeshCADForContrast.push_back(tempface);
	//		//indexFaceHandleForMeshCADContrast++;
	//	}
	//}
	//SphereMesh_MeshCADForContrast = new MeshKernel::SurfaceMesh(VerticesMeshCADForContrast, FacesMeshCADForContrast);

	//获取控制点
	vector<vector<CPoint>> singleSurfaceControlPoints;
	vector<vector<vector<CPoint>>> allSurfaceControlPoints;
	for (int i = 0; i < surfacesDraw.size(); i++) 
	{
		singleSurfaceControlPoints = surfacesDraw[i].getControlPoints();
		allSurfaceControlPoints.push_back(singleSurfaceControlPoints);
	}
	this->controlMesh = allSurfaceControlPoints;
	time_t end = clock();
	printf("bsplineSurfaceToDiscrete cost: %zdms, vertices size: %zd, surface size: %zd\n", int(end - beg), totalVerticesSize, allSurface.size());
	updateDrawData();

}

void ModelDraw::controlpointsToDiscrete()
{

	data_controlMesh.clear();
	vector<CPoint>singlePatch;//每个面片的离散点
	vector<vector<CPoint>>singlelSurface;//每个曲面的离散点
	vector<vector<vector<CPoint>>>allSurface;//所有曲面集合的离散点
	// 这里可以看章的毕业论文，一个控制网格面， 也需要多个面片生成
	for (int i = 0; i <  controlMesh.size(); i++) 
	{//controlMesh[i]是第i个曲面
		for (int j = 0; j < controlMesh[i].size() - 1; j++) 
		{//controlMesh[i][j]是第i个曲面控制点集合的第j行
			for (int k = 0; k < controlMesh[i][j].size() - 1; k++)
			{//controlMesh[i][j][k]是第i个曲面的第j行第k列的控制点
				singlePatch.push_back(controlMesh[i][j][k]);
				singlePatch.push_back(controlMesh[i][j][k + 1]);
				singlePatch.push_back(controlMesh[i][j + 1][k + 1]);
				//singlePatch.push_back(controlMesh[i][j + 1][k + 1]);
				singlePatch.push_back(controlMesh[i][j + 1][k]);
				singlelSurface.push_back(singlePatch);
				singlePatch.clear();
			}
		}
		allSurface.push_back(singlelSurface);
		singlelSurface.clear();
	}
	this->data_controlMesh = allSurface;
}

void ModelDraw::updateSurfaceModel() 
{
	if (surfacesDraw.empty()) 
	{
		return;
	}
	bbMin = Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	bbMax = Vector3d(-99999999, -99999999, -99999999);
	for (auto& vv_p : data_bound_surface) {
		for (auto& v_p : vv_p) {
			for (auto& p : v_p) {
				double x = p.getX(), y = p.getY(), z = p.getZ();
				bbMin[0] = std::min(bbMin[0], x), bbMin[1] = std::min(bbMin[1], y), bbMin[2] = std::min(bbMin[2], z);
				bbMax[0] = std::max(bbMax[0], x), bbMax[1] = std::max(bbMax[1], y), bbMax[2] = std::max(bbMax[2], z);
			}
		}
	}
	double u1, u2, v1, v2;
	surfacesDraw[0].getUVMinMax(u1, u2, v1, v2);
	double r = 1.5 * sqrt(
		pow((surfacecenter.getX() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getX()), 2) +
		pow((surfacecenter.getY() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getY()), 2) +
		pow((surfacecenter.getZ() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getZ()), 2)
	);
	SetScenePosition(Vector3d(surfacecenter.getX(), surfacecenter.getY(), surfacecenter.getZ()), r);
	update();
}

std::vector<vector<std::vector<CPoint>>> ModelDraw::getData_boundcurve(vector<CBSplineSurface> boundsurface) 
{
	vector<CPoint>singleBoundPoints;//一条边界上的离散点；
	vector<vector<CPoint>>singleSurfacePoints;//一张曲面上的离散点；
	vector<vector<vector<CPoint>>>allSurfacePoints;//所有边界离散点集合；
	CPoint p(0, 0, 0);
	int number = 0;
	for (int i = 0; i < boundsurface.size(); i++) 
	{
		double umin, umax, vmin, vmax;
		boundsurface[i].getUVMinMax(umin, umax, vmin, vmax); //?
		int unumber, vnumber;
		unumber = upatchnum * (boundsurface[i].getSizeUControlPoints());
		vnumber = vpatchnum * (boundsurface[i].getSizeVControlPoints());
		number += 2 * (unumber + vnumber);
		double sizeu = (umax - umin) / unumber - 1.e-9;
		double sizev = (vmax - vmin) / vnumber - 1.e-9;
		for (double u = umin + 1.e-9; u <= umax - 1.e-9; u = u + sizeu)
		{//第一条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(u, vmin + 1.e-9));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
			
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		for (double u = umin + 1.e-9; u <= umax - 1.e-9; u = u + sizeu) 
		{//第二条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(u, vmax - 1.e-9));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		for (double v = vmin + 1.e-9; v <= vmax - 1.e-9; v = v + sizev) 
		{//第三条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(umin + 1.e-9, v));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		//for (double v = vmin + 1.e-9; v <= vmax - 1.e-9; v = v + sizev)
		//{//第四条边界
		//	singleBoundPoints.push_back(boundsurface[i].getPointAtParam(umax - 1.e-9, v));
		//	p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		//}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		allSurfacePoints.push_back(singleSurfacePoints);
		singleSurfacePoints.clear();
		
		//break; // hzk测试
	}
	surfacecenter.setXYZ(p.getX() / number, p.getY() / number, p.getZ() / number);
	return allSurfacePoints;
}

void ModelDraw::drawControlMesh(void)
{
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	for (int i = 0; i < data_controlMesh.size(); i++)
	{
		for (int j = 0; j < data_controlMesh[i].size(); j++)
		{
			for (int k = 0; k < data_controlMesh[i][j].size(); k++) 
			{
				drawControlPoint(data_controlMesh[i][j][k].getX(), data_controlMesh[i][j][k].getY(), data_controlMesh[i][j][k].getZ(), 0.01, 16, 16);
			}
		}
	}
	glLineWidth(1.0);
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int i = 0; i < data_controlMesh.size(); i++) 
	{
		glBegin(GL_QUADS);
		for (int j = 0; j < data_controlMesh[i].size(); j++) 
		{
			for (int k = 0; k < data_controlMesh[i][j].size(); k++) 
			{
				glVertex3d(data_controlMesh[i][j][k].getX(), data_controlMesh[i][j][k].getY(), data_controlMesh[i][j][k].getZ()); //这个我要做的
			}
		}
		glEnd();
	}
	//编辑模式下，使用填充样式
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void ModelDraw::drawSurface() 
{
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
	shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_QUADS, 0, vertices_data.size() / 6);
}

void ModelDraw::DrawPointForOrder()
{
	////模仿updateDrawData_Surface写的
	//vertices_data.clear();
	////colors[0] = { 227.0 /255, 121.0/255, 77.0/255 };
	//int index = 0;
	//for (int i = 0; i < data_controlMesh.size(); ++i) //data_controlMesh已经把压入的顺序整理好了，直接压入就可以了
	//{
	//	//对于每一个面
	//	auto& singleSurface = data_controlMesh[i];
	//	auto& color = colors[i % colors.size()];
	//	auto& normal = surface_normals[i];
	//	for (int j = 0; j < singleSurface.size(); j++) {
	//		//每一个面的9个面片

	//		//index++;
	//		//auto& color = colors[index % colors.size()];
	//		for (int k = 0; k < singleSurface[j].size(); k++) {
	//			//每个面片有4个点
	//			vertices_data.push_back(singleSurface[j][k].getX());
	//			vertices_data.push_back(singleSurface[j][k].getY());
	//			vertices_data.push_back(singleSurface[j][k].getZ());
	//			vertices_data.push_back(color.x());
	//			vertices_data.push_back(color.y());
	//			vertices_data.push_back(color.z());
	//			//vertices_data.push_back(normal[i][j].getX());
	//			//vertices_data.push_back(normal[i][j].getY());
	//			//vertices_data.push_back(normal[i][j].getZ());
	//		}

	//	}
	//}
	//update();

	//double size = 30.0;
	//glPointSize(2.5);
	shaderProgram.setUniformValue("colormode", 232);
	//glBegin(GL_POINTS);
	for (int i = 0; i < data_controlMesh.size(); ++i) //data_controlMesh已经把压入的顺序整理好了，直接压入就可以了
	{
		//对于每一个面
		double size = 161.0;
		auto& singleSurface = data_controlMesh[i];
		for (int j = 0; j < singleSurface.size(); j++) 
		{
			for (int k = 0; k < singleSurface[j].size(); k++) {
				//每个面片有4个点
				glPointSize(size);
				glBegin(GL_POINTS);
				glVertex3f(singleSurface[j][k].getX(), singleSurface[j][k].getY(), singleSurface[j][k].getZ());
				//vertices_data.push_back(normal[i][j].getX());
				//vertices_data.push_back(normal[i][j].getY());
				//vertices_data.push_back(normal[i][j].getZ());
				glEnd();
				size = size - 10;
			}

		}
	}
	//for (const auto& vp : curve_skel_->allvertices())
	//{
	//	// draw_->DrawDiscretizeSphere(vp.second.x(), vp.second.y(), vp.second.z(), 5);
	//	glVertex3f(vp.second.x(), vp.second.y(), vp.second.z());

	//}
	glEnd();
}

void ModelDraw::updateControlPointAndLine()
{
	//模仿updateDrawData_Surface写的
	vertices_data.clear();
	//colors[0] = { 227.0 /255, 121.0/255, 77.0/255 };
	int index = 0;
	for (int i = 0; i < data_controlMesh.size(); ++i) //data_controlMesh已经把压入的顺序整理好了，直接压入就可以了
	{
		//对于每一个面
		auto& singleSurface = data_controlMesh[i];
		auto& color = colors[i % colors.size()];
		auto& normal = surface_normals[i];
		for (int j = 0; j < singleSurface.size(); j++) {
			//每一个面的9个面片
		
			//index++;
			//auto& color = colors[index % colors.size()];
			
			for (int k = 0; k < singleSurface[j].size(); k++) 
			{
					//每个面片有4个点
					vertices_data.push_back(singleSurface[j][k].getX());
					vertices_data.push_back(singleSurface[j][k].getY());
					vertices_data.push_back(singleSurface[j][k].getZ());
					vertices_data.push_back(color.x());
					vertices_data.push_back(color.y());
					vertices_data.push_back(color.z());
					//vertices_data.push_back(normal[i][j].getX());
					//vertices_data.push_back(normal[i][j].getY());
					//vertices_data.push_back(normal[i][j].getZ());
			}
			
			//if (j == 0)
			//{
			//	for (int k = 0; k < singleSurface[j].size(); k++) {
			//		//每个面片有4个点
			//		vertices_data.push_back(singleSurface[j][k].getX());
			//		vertices_data.push_back(singleSurface[j][k].getY());
			//		vertices_data.push_back(singleSurface[j][k].getZ());
			//		vertices_data.push_back(color.x());
			//		vertices_data.push_back(color.y());
			//		vertices_data.push_back(color.z());
			//		//vertices_data.push_back(normal[i][j].getX());
			//		//vertices_data.push_back(normal[i][j].getY());
			//		//vertices_data.push_back(normal[i][j].getZ());
			//	}
			//}
			//else if ( j == 8)
			//{
			//	for (int k = 0; k < singleSurface[j].size(); k++) {
			//		//每个面片有4个点
			//		vertices_data.push_back(singleSurface[j][k].getX());
			//		vertices_data.push_back(singleSurface[j][k].getY());
			//		vertices_data.push_back(singleSurface[j][k].getZ());
			//		vertices_data.push_back(color.x());
			//		vertices_data.push_back(color.y());
			//		vertices_data.push_back(color.z());
			//		//vertices_data.push_back(normal[i][j].getX());
			//		//vertices_data.push_back(normal[i][j].getY());
			//		//vertices_data.push_back(normal[i][j].getZ());
			//	}
			//}
			//else
			//{
			//	for (int k = 0; k < singleSurface[j].size(); k++) {
			//		//每个面片有4个点
			//		vertices_data.push_back(singleSurface[0][k].getX());
			//		vertices_data.push_back(singleSurface[0][k].getY());
			//		vertices_data.push_back(singleSurface[0][k].getZ());
			//		vertices_data.push_back(color.x());
			//		vertices_data.push_back(color.y());
			//		vertices_data.push_back(color.z());
			//		//vertices_data.push_back(normal[i][j].getX());
			//		//vertices_data.push_back(normal[i][j].getY());
			//		//vertices_data.push_back(normal[i][j].getZ());
			//	}
			//}
			
		}
	}
	update();
}

void ModelDraw::drawControlPointAndLine()
{
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
	shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_QUADS, 0, vertices_data.size() / 6);

	/*core->glDrawArrays(GL_POINTS, 0, vertices_data.size() / 6);
	core->glDrawArrays(GL_LINES, 0, vertices_data.size() / 6);*/
}

void ModelDraw::drawControlPoint(double xx, double yy, double zz, double radius, GLfloat M, GLfloat N) {
	float step_z = PI / M;
	float step_xy = 2 * PI / N;
	double x[4], y[4], z[4];
	float angle_z = 0.0;
	float angle_xy = 0.0;
	int i = 0, j = 0;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);
	for (i = 0; i < M; i++) {
		angle_z = i * step_z;

		for (j = 0; j < N; j++) {
			angle_xy = j * step_xy;

			x[0] = radius * sin(angle_z) * cos(angle_xy);
			y[0] = radius * sin(angle_z) * sin(angle_xy);
			z[0] = radius * cos(angle_z);

			x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);
			y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);
			z[1] = radius * cos(angle_z + step_z);

			x[2] = radius * sin(angle_z + step_z) * cos(angle_xy + step_xy);
			y[2] = radius * sin(angle_z + step_z) * sin(angle_xy + step_xy);
			z[2] = radius * cos(angle_z + step_z);

			x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);
			y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);
			z[3] = radius * cos(angle_z);
			for (int k = 0; k < 4; k++) {
				double distance = sqrt(pow(x[k], 2) + pow(y[k], 2) + pow(z[k], 2));
				glNormal3d(x[k] / distance, y[k] / distance, z[k] / distance);
				glVertex3d(xx + x[k], yy + y[k], zz + z[k]);
			}
		}
	}
	glEnd();

}

void ModelDraw::updateDrawData_Surface() {
	vertices_data.clear();
	//colors[0] = { 227.0 /255, 121.0/255, 77.0/255 };
	for (int i = 0; i < data_surfaces.size(); ++i) {//102个贝塞尔体
		auto& singleSurface = data_surfaces[i]; // 每个贝塞尔面 每个贝塞尔面里有64个四边形，每个四边形有4个点
		auto& color = colors[i % colors.size()];
		auto& normal = surface_normals[i];
		//int index = 0;
		for (int i = 0; i < singleSurface.size(); i++) {
			//auto& color = colors[index % colors.size()];
			//index++;
			for (int j = 0; j < singleSurface[i].size(); j++) {
				//auto& color = colors[j % colors.size()];
				vertices_data.push_back(singleSurface[i][j].getX());
				vertices_data.push_back(singleSurface[i][j].getY());
				vertices_data.push_back(singleSurface[i][j].getZ());
				vertices_data.push_back(color.x());
				vertices_data.push_back(color.y());
				vertices_data.push_back(color.z());
				//vertices_data.push_back(normal[i][j].getX());
				//vertices_data.push_back(normal[i][j].getY());
				//vertices_data.push_back(normal[i][j].getZ());
			}
		}
	}
	update();
}


void ModelDraw::updateDrawData_Surface_IGES_Error(SurfaceMesh& trimesh) 
{
	//hzk测试 可以删掉
	vertices_data.clear();

	SurfaceMesh scrmesh;
	vector<vector<vector<VH>>> to_scrmesh_vh;
	for (int i = 0; i < data_surfaces.size(); ++i) 
	{
		auto& singleSurface = data_surfaces[i];
		vector<vector<VH>> one_surfaces;
		for (int j = 0; j < singleSurface.size(); j++)
		{
			vector<VH> one_patch;
			for (int k = 0; k < singleSurface[j].size(); k++) 
			{
				CPoint cvec = singleSurface[j][k];
				one_patch.push_back(scrmesh.AddVertex(Vec(cvec.getX(),cvec.getY(),cvec.getZ())));
			}
			one_surfaces.push_back(one_patch);
			scrmesh.AddFace(one_patch);
		}
		to_scrmesh_vh.push_back(one_surfaces);
	}
	MeshKernel::IO file;
	file.WriteOffFile(scrmesh, "C:\\Users\\swords\\Desktop\\k.off");

	std::unordered_map<VH, Vector3d> vh2color;
	ErrorRender::get_Hausdorff_distance(scrmesh, trimesh, vh2color);

	for (int i = 0; i < data_surfaces.size(); ++i) 
	{
		auto& singleSurface = data_surfaces[i];
		//auto& color = colors[i % colors.size()];
		auto& normal = surface_normals[i];
		for (int j = 0; j < singleSurface.size(); j++) 
		{
			for (int k = 0; k < singleSurface[j].size(); k++) 
			{
				auto color = vh2color[to_scrmesh_vh[i][j][k]];
				//Vector3d color(1, 1, 1);
				vertices_data.push_back(singleSurface[j][k].getX());
				vertices_data.push_back(singleSurface[j][k].getY());
				vertices_data.push_back(singleSurface[j][k].getZ());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
				//vertices_data.push_back(normal[i][j].getX());
				//vertices_data.push_back(normal[i][j].getY());
				//vertices_data.push_back(normal[i][j].getZ());
			}
		}
	}
	update();
}

void ModelDraw::updateDrawData_SurfaceMesh_Error_IGES() {
	if (!triMesh) return;
	vertices_data.clear();
	vertices_data.reserve(triMesh->FaceSize() * 18);

	std::vector<Vector3f> discrete_vs;
	for (int i = 0; i < data_surfaces.size(); ++i) {
		auto& singleSurface = data_surfaces[i];
		for (int j = 0; j < singleSurface.size(); j++) {
			vector<Vector3f> one_patch;
			assert(singleSurface[j].size() == 4);
			for (int k = 0; k < singleSurface[j].size(); k++) {
				CPoint cvec = singleSurface[j][k];
				one_patch.push_back(Vector3f(cvec.getX(), cvec.getY(), cvec.getZ()));
			}
			discrete_vs.emplace_back(one_patch[0]);
			discrete_vs.emplace_back(one_patch[1]);
			discrete_vs.emplace_back(one_patch[2]);
			discrete_vs.emplace_back(one_patch[0]);
			discrete_vs.emplace_back(one_patch[2]);
			discrete_vs.emplace_back(one_patch[3]);
		}
	}

	AABB_Tree aabb_tree(discrete_vs);
	std::unordered_map<VH, Vector3d> vh2color;
	ErrorRender::get_Hausdorff_distance(*triMesh, aabb_tree, vh2color);

	for (auto& fp : triMesh->allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		if (vhs.size() == 3) {
			for (auto& vh : vhs) {
				Vector3d color = vh2color[vh];
				auto& v = triMesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}
		}

	}
	tri_cnt = vertices_data.size() / 18;

	update();
}

void ModelDraw::genColors(size_t sz) {
	colors.resize(sz);
	srand((unsigned)time(NULL));
	for (int i = 0; i < sz; ++i) {
		double R = ((std::rand() % 9999) + 1) / 10000.f;
		double G = ((std::rand() % 9999) + 1) / 10000.f;
		double B = ((std::rand() % 9999) + 1) / 10000.f;
		colors[i] = Eigen::Vector3d(R, G, B);
	}
}

void ModelDraw::triangle_meshing() {

	if (!hexMesh) return;

	triMesh = new MeshKernel::SurfaceMesh();

	TriMesh_Generator app;
	app.get_triangle_mesh(*hexMesh, *triMesh);

	delete hexMesh;
	hexMesh = nullptr;

	SetDrawModule(SURFACEMESH);
	updateDrawData();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_1() {
	vertices_data.clear();

	auto propsSke = OpenMesh::FProp<int>(*poly_mesh, "ske_flag");
	auto propsBdy = OpenMesh::FProp<int>(*poly_mesh, "bdy_flag");

	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		int s_flag = propsSke[f_it], b_flag = propsBdy[f_it];
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			if (s_flag == 1) {
				vertices_data.push_back(0.4745);
				vertices_data.push_back(0.6902);
				vertices_data.push_back(0.3608);
			} else if (b_flag == -1) {
				vertices_data.push_back(0.2667);
				vertices_data.push_back(0.001f);
				vertices_data.push_back(0.3294);
			} else {
				vertices_data.push_back(0.1882);
				vertices_data.push_back(0.4039);
				vertices_data.push_back(0.4353);
			}

		}

	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_2() {
	vertices_data.clear();

	auto propsBdy = OpenMesh::FProp<int>(*poly_mesh, "ske_sp_flag");
	int max_value = -9999999, min_value = 9999999;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		int value = propsBdy[f_it];
		min_value = std::min(min_value, value);
		max_value = std::max(max_value, value);
	}
	double range_value = max_value - min_value;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		int value = propsBdy[f_it];
		Vector3d color = colorbar.getColor((value - min_value) * 1.f / range_value);
		//OpenMesh::Vec3f N = poly_mesh->normal(*f_it).normalized();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			vertices_data.push_back(color[0]);
			vertices_data.push_back(color[1]);
			vertices_data.push_back(color[2]);
		}

	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_V() {
	vertices_data.clear();

	auto propsBdy = OpenMesh::VProp<double>(*poly_mesh, "mmc_p_val");
	int cnt_val = 0;
	double max_value = -9999999, min_value = 9999999;
	for (auto v_it = poly_mesh->vertices_begin(); v_it != poly_mesh->vertices_end(); ++v_it) {
		double value = propsBdy[v_it];
		if (value < 0.f) value = 0.f;
		min_value = std::min(min_value, value);
		max_value = std::max(max_value, value);
		//if (cnt_val++ % 50 == 0) std::cout << "vertex value: " << value << std::endl;
	}
	//std::cout << "\n\n--------------------\n\n" << std::endl;
	double range_value = max_value - min_value;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		//OpenMesh::Vec3f N = poly_mesh->normal(*f_it).normalized();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			double value = propsBdy[fv_it];
			//if (cnt_val++ % 50 == 0) std::cout << "vertex value: " << value << std::endl;
			double k = (value - min_value) * 1.f / range_value;
			// 将 0 ~ 0.1 拉长到 0 ~ 0.5, 0.1 ~ 1 压缩到 0.5 ~ 1
			if (k <= 0.1) k *= 5;
			else {
				k = 0.5 + 0.5 * (k - 0.1) / 0.9;
			}
			Vector3d color = colorbar.getColor(k);
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			vertices_data.push_back(color[0]);
			vertices_data.push_back(color[1]);
			vertices_data.push_back(color[2]);
		}

	}

	for (int i = 0; i < bdy.size(); ++i) {
		auto& line = bdy[i];
		auto& color = colors[i % colors.size()];
		for (int j = 1; j < line.size(); ++j) {
			auto& p1 = line[j - 1], & p2 = line[j];
			vertices_data.push_back(p1.x()), vertices_data.push_back(p1.y()), vertices_data.push_back(0.f);
			vertices_data.push_back(color.x()), vertices_data.push_back(color.y()), vertices_data.push_back(color.z());
			vertices_data.push_back(p2.x()), vertices_data.push_back(p2.y()), vertices_data.push_back(0.f);
			vertices_data.push_back(color.x()), vertices_data.push_back(color.y()), vertices_data.push_back(color.z());
		}
	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::drawPolyMesh_OPT() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	if (axis.empty()) {
		shaderProgram.setUniformValue("colormode", -1);
		auto total_vcnt = vertices_data.size() / 6;
		auto f_vcnt = poly_mesh->n_faces();
		core->glDrawArrays(GL_QUADS, 0, f_vcnt * 4);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width * 2.f);
		GLsizei l_vcnt = total_vcnt - f_vcnt;
		core->glDrawArrays(GL_LINES, f_vcnt * 4, l_vcnt);// 绘制折线

		/*if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0, 1.0);
			glLineWidth(meshline_width);
			shaderProgram.setUniformValue("colormode", 530);
			for (int i = 0; i < f_vcnt; ++i) {
				core->glDrawArrays(GL_LINE_LOOP, i * 4, 4);
			}
		}*/

	} else {
		// draw 
		shaderProgram.setUniformValue("colormode", 530);
		glLineWidth(meshline_width);
		glBegin(GL_LINES);
		for (auto& pos_pair : axis) {
			auto& p1 = pos_pair.first;
			auto& p2 = pos_pair.second;
			glVertex3f(p1.x(), p1.y(), 0);
			glVertex3f(p2.x(), p2.y(), 0);
		}
		glEnd();

	}

	update();
}


//bool ModelDraw::BuildSkeletonToFile(std::string filename, uint32_t width, uint32_t height) {
//	ModelClear();
//	int dot_idx = filename.find_last_of('.');
//	string fileformat = filename.substr(++dot_idx, filename.length());
//	//std::cout << "file format: " << fileformat << std::endl;
//	std::cout << "filename: " << filename << std::endl;
//	std::cout << "width = " << width << ", height = " << height << std::endl;
//
//	MeshKernel::IO io;
//	assert(fileformat == "obj");
//	time_t t1 = clock();
//	//poly_mesh = pipeline->BuildSkeletonToFile(filename, 722, 242);
//	pipeline = new Pipeline();
//	poly_mesh = new polyMesh();
//	*poly_mesh = pipeline->BuildSkeletonToFile(filename, width, height);
//	std::cout << "Build Skeleton Cost: " << int(clock() - t1) << " ms" << std::endl;
//
//	SetDrawModule(MMC_Convert);
//	updateMeshModel();
//	updateDrawData_PolyMesh_OPT_1();
//	update();
//	return true;
//
//}

//void ModelDraw::SparseOptimization(double max_error) {
//
//	time_t t1 = clock();
//	*poly_mesh = pipeline->SparseOptimization(1, max_error);
//	std::cout << "Sparse Optimization Cost: " << int(clock() - t1) << " ms" << std::endl;
//	updateDrawData_PolyMesh_OPT_2();
//	update();
//
//}

//void ModelDraw::GenerateSymAxis() {
//	time_t t1 = clock();
//	pipeline->GenerateSymAxis(axis);
//	std::cout << "Generate SymAxis Cost: " << int(clock() - t1) << " ms" << std::endl;
//	update();
//}

//void ModelDraw::ShapeApproximation(uint32_t max_iter) {
//	time_t t1 = clock();
//	*poly_mesh = pipeline->ShapeApproximation(max_iter);
//
//	pipeline->ExtractMMCBoundary(bdy);
//	axis.clear();
//	std::cout << "Shpae Approximation Cost: " << int(clock() - t1) << " ms" << std::endl;
//	updateDrawData_PolyMesh_OPT_V();
//	update();
//}



HexV3fMesh* ModelDraw::convent() {
	HexV3fMesh* _hexMesh = new HexV3fMesh();

	/* iGame_hexMesh to XIE_hexMesh */
	{
		uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();
		for (uint64_t i = 0; i < vcnt; ++i) {
			MeshKernel::iGameVertexHandle vh(i);
			auto& v = hexMesh->vertices(vh);
			V3f _v(v.x(), v.y(), v.z());
			_hexMesh->add_vertex(_v);
		}

		for (uint64_t i = 0; i < ccnt; ++i) {
			MeshKernel::iGameCellHandle ch(i);
			if (hexMesh->isValid(ch)) {
				auto& c = hexMesh->cells(ch);
				auto vhs = c.getVertexHandle();
				vector<VertexHandle> _vhs;
				for (auto& vh : vhs) {
					_vhs.push_back(VertexHandle(vh));
				}
				_hexMesh->add_cell(_vhs);
			}
		}

		/*std::cout << "iGameHexMesh convert to OldHexMesh success. OldHexMesh: #V = " << _hexMesh->vertices_.size() <<
		 ", #C = " << _hexMesh->cells_.size() << std::endl;*/
	}


	/*FileHandler file(filename);
	assert(file.read2HexMeshObj(_hexMesh, true, true, true));*/

	if (_hexMesh->vertices_.size() == 0) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
		return nullptr;
	}
	return _hexMesh;
}

