#include "coordinateAxis.h"
#include "ModelDraw.h"
#include "GL/glut.h"
#include <gl/GLU.h>
#include <gl/GL.h>

void CoordinateAxis::SetModelMatrix(ModelDraw& modelDraw, QMatrix4x4 matrix) {
	modelDraw.shaderProgram.setUniformValue("model", modelDraw.modelmatrix * matrix);	// 右乘先操作
}

void CoordinateAxis::Translated(double x, double y, double z) {
	QMatrix4x4 T = QMatrix4x4(
		1, 0, 0, x,
		0, 1, 0, y,
		0, 0, 1, z,
		0, 0, 0, 1);
	translate_matrix_ = T * translate_matrix_;													// 左乘后操作
	model_matrix_ = translate_matrix_ * rotate_matrix_;
}

void CoordinateAxis::Rotatef(double angle, double x, double y, double z) {
	double theta = angle / 180 * PI;
	double s = sin(theta), c = cos(theta), t = 1 - cos(theta);
	QMatrix4x4 T = QMatrix4x4(
		t * x * x + c, t * x * y + s * z, t * x * z - s * y, 0,
		t * x * y - s * z, t * y * y + c, t * y * z + s * x, 0,
		t * x * z + s * y, t * y * z - s * x, t * z * z + x, 0,
		0, 0, 0, 1);
	rotate_matrix_ = T * rotate_matrix_;														// 左乘后操作
	model_matrix_ = translate_matrix_ * rotate_matrix_;
}

void CoordinateAxis::PushMatrix() {
	temp_matrixs_.push(translate_matrix_);
	temp_matrixs_.push(rotate_matrix_);
}

void CoordinateAxis::PopMatrix() {
	if (!temp_matrixs_.empty()) {
		rotate_matrix_ = temp_matrixs_.top();
		temp_matrixs_.pop();
		translate_matrix_ = temp_matrixs_.top();
		temp_matrixs_.pop();
		model_matrix_ = translate_matrix_ * rotate_matrix_;
	}
}

void CoordinateAxis::Draw3DCoord(ModelDraw& modelDraw)
{
	if (mode_type_ == ModeType::HIDE) return;

	modelDraw.shaderProgram.setUniformValue("colormode", 3);
	modelDraw.shaderProgram.setUniformValue("color", QVector4D(1, 1, 0, 1));
	glDepthRange(0.0, 0.01);
	GLUquadric* obj = gluNewQuadric();
	model_matrix_.setToIdentity(), rotate_matrix_.setToIdentity(), translate_matrix_.setToIdentity();
	Translated(x_, y_, z_);
	SetModelMatrix(modelDraw, model_matrix_);
	/*
	glPointSize(7);
	glBegin(GL_POINTS);
	glVertex3f(t_point1.x(), t_point1.y(), t_point1.z());
	glEnd();
	modelDraw.shaderProgram.setUniformValue("color", QVector4D(0, 1, 0, 1));
	glBegin(GL_POINTS);
	glVertex3f(t_point2.x(), t_point2.y(), t_point2.z());
	glEnd();
	*/
	gluSphere(obj, len * 0.05, 20, 20);
	glLineWidth(3.0);
	{
		// x axis
		if (click_type_ == ClickType::X_AXIS_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 48.0/255, 48.0/255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 105.0/255, 180.0/255, 1));
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(len*k, 0, 0);
		glEnd();
		//QVector3D point = QVector3D(len - vlen, 0.0, 0.0);
		//QVector3D dir = QVector3D(1.0, 0.0, 0.0);
		//DrawCone(point, dir, 0.01, vlen);
		PushMatrix();
		Rotatef(-90, 0, 1, 0);
		Translated(len - vlen, 0.0, 0.0);
		SetModelMatrix(modelDraw, model_matrix_);
		// 坐标为obj中心 底座半径、顶座半径、长度
		gluCylinder(obj, pick_threshold, 0.00, vlen, 30, 30);
		PopMatrix();

		// y axis
		if (click_type_ == ClickType::Y_AXIS_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 48.0 / 255, 48.0 / 255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 105.0 / 255, 180.0 / 255, 1));
		SetModelMatrix(modelDraw, model_matrix_);
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, len * k, 0);
		glEnd();
		PushMatrix();
		//model_matrix_.setToIdentity(), rotate_matrix_.setToIdentity(), translate_matrix_.setToIdentity();
		Rotatef(-90, 0, 1, 0);
		Rotatef(-90, 0, 0, 1);
		Translated(0.0, len - vlen, 0.0);
		SetModelMatrix(modelDraw, model_matrix_);
		// 坐标为obj中心 底座半径、顶座半径、长度
		gluCylinder(obj, pick_threshold, 0.00, vlen, 30, 30);
		PopMatrix();

		// z axis
		if (click_type_ == ClickType::Z_AXIS_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 48.0 / 255, 48.0 / 255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 105.0 / 255, 180.0 / 255, 1));
		SetModelMatrix(modelDraw, model_matrix_);
		glBegin(GL_LINES);
		glVertex3f(0, 0, 0);
		glVertex3f(0, 0, len * k);
		glEnd();
		PushMatrix();
		Translated(0.0, 0.0, len - vlen);
		SetModelMatrix(modelDraw, model_matrix_);
		// 坐标为obj中心 底座半径、顶座半径、长度
		gluCylinder(obj, pick_threshold, 0.00, vlen, 30, 30);
		PopMatrix();
	}
	
	if (mode_type_ == ModeType::TRANSFROM)
	{
		// normal z circle
		if (click_type_ == ClickType::Z_CIRCLE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(1.0, 1.0, 0.0, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(173.0/255, 1.0, 47.0/255, 1));
		SetModelMatrix(modelDraw, model_matrix_);
		glBegin(GL_LINE_LOOP);
		for (double angle = 0.0; angle <= 2 * PI; angle += 0.05) {
			glVertex3d(cos(angle) * radios, sin(angle) * radios, 0.0);
		}
		glEnd();

		// normal x circle
		if (click_type_ == ClickType::X_CIRCLE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(255, 255, 0, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(173.0 / 255, 1.0, 47.0 / 255, 1));
		glBegin(GL_LINE_LOOP);
		for (double angle = 0.0; angle <= 2 * PI; angle += 0.05) {
			glVertex3d(0.0, sin(angle) * radios, cos(angle) * radios);
		}
		glEnd();

		// normal y circle
		if (click_type_ == ClickType::Y_CIRCLE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(255, 255, 0, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(173.0 / 255, 1.0, 47.0 / 255, 1));
		glBegin(GL_LINE_LOOP);
		for (double angle = 0.0; angle <= 2 * PI; angle += 0.05) {
			glVertex3d(cos(angle) * radios, 0.0, sin(angle) * radios);
		}
		glEnd();
	}
	else if (mode_type_ == ModeType::EDIT) {
		// x
		if (click_type_ == ClickType::X_CUBE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(46.0 / 255, 139.0 / 255, 87.0 / 255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(84.0 / 255, 1.0, 159.0 / 255, 1));
		PushMatrix();
		Translated(len + clen, 0.0, 0.0);
		SetModelMatrix(modelDraw, model_matrix_);
		DrawCube();
		PopMatrix();
		// y
		if (click_type_ == ClickType::Y_CUBE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(46.0 / 255, 139.0 / 255, 87.0 / 255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(84.0 / 255, 1.0, 159.0 / 255, 1));
		PushMatrix();
		Translated(0.0, len + clen, 0.0);
		SetModelMatrix(modelDraw, model_matrix_);
		DrawCube();
		PopMatrix();
		// z
		if (click_type_ == ClickType::Z_CUBE_CLICK)
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(46.0 / 255, 139.0 / 255, 87.0 / 255, 1));
		else
			modelDraw.shaderProgram.setUniformValue("color", QVector4D(84.0 / 255, 1.0, 159.0 / 255, 1));
		PushMatrix();
		Translated(0.0, 0.0, len + clen);
		SetModelMatrix(modelDraw, model_matrix_);
		DrawCube();
		PopMatrix();
	}
	glDepthRange(0.01, 0.99999);
	modelDraw.shaderProgram.setUniformValue("model", modelDraw.modelmatrix);
}

void CoordinateAxis::DrawCone(QVector3D point, QVector3D dir, double radios, double len) {
	auto GetReal = [=](double angle)-> QVector3D {
		QVector3D local_circle_point;
		if (dir.z() == 1.0) //z
			local_circle_point = QVector3D(cos(angle) * radios, sin(angle) * radios, 0.0);
		else if (dir.x() == 1.0) //x
			local_circle_point = QVector3D(0.0, sin(angle) * radios, cos(angle) * radios);
		else //y
			local_circle_point = QVector3D(cos(angle) * radios, 0.0, sin(angle) * radios);
		return point + local_circle_point;
	};
	//double PI = acos(-1.0);
	glBegin(GL_LINE_LOOP);
	for (double angle = 0.0; angle <= 2 * PI; angle += 0.1) {
		QVector3D real_circle_point = GetReal(angle);
		glVertex3d(real_circle_point.x(), real_circle_point.y(), real_circle_point.z());
	}
	glEnd();
	QVector3D point_to = point + dir * len;
	glBegin(GL_TRIANGLES);
	QVector3D pre_circle_point = GetReal(2 * PI);
	for (double angle = 0.0; angle < 2 * PI; angle += 0.1) {
		QVector3D real_circle_point = GetReal(angle);
		glVertex3d(point_to.x(), point_to.y(), point_to.z());
		glVertex3d(real_circle_point.x(), real_circle_point.y(), real_circle_point.z());
		glVertex3d(pre_circle_point.x(), pre_circle_point.y(), pre_circle_point.z());
	}
	glEnd();
}

void CoordinateAxis::DrawCube() {
	QVector3D tmp[8] = { {1,1,1},{1,-1,1},{-1,1,1},{-1,-1,1},{1,1,-1},{1,-1,-1},{-1,1,-1},{-1,-1,-1} };
	for (int i = 0; i < 8; i++)
		tmp[i] = tmp[i] * clen;
	auto DrawQuat = [&](int* order) {
		glBegin(GL_TRIANGLES);
		for (int i = 0; i < 2; i++) {
			glVertex3d(tmp[order[i]].x(), tmp[order[i]].y(), tmp[order[i]].z());
			glVertex3d(tmp[order[i + 1]].x(), tmp[order[i + 1]].y(), tmp[order[i + 1]].z());
			glVertex3d(tmp[order[i + 2]].x(), tmp[order[i + 2]].y(), tmp[order[i + 2]].z());
		}
		glEnd();
	};
	int order[6][4] = { {0, 1, 2, 3}, {0, 1, 4, 5}, {1, 2, 5, 6}, {2, 3, 6, 7}, {4, 5, 6, 7} };
	for(int i = 0; i < 6; i++)
		DrawQuat(order[i]);
}

/*
* 返回坐标轴的3D位置变化
*/
QVector3D CoordinateAxis::MoveEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP) {
	std::cout << "move event" << std::endl;
	QMatrix4x4 inv_matrix = modelDraw.modelmatrix.inverted() * modelDraw.viewmatrix.inverted() * modelDraw.projectionmatrix.inverted();
	vector<QVector3D> intersect_pt_old = IntersectPlanePoint(inv_matrix, oldP.x(), oldP.y());
	vector<QVector3D> intersect_pt_new = IntersectPlanePoint(inv_matrix, newP.x(), newP.y());
	QVector3D ans = QVector3D(0, 0, 0);
	switch (click_type_) {
	case ClickType::X_AXIS_CLICK:
		ans.setX(intersect_pt_new[1][0] - intersect_pt_old[1][0]);
		break;
	case ClickType::Y_AXIS_CLICK:
		ans.setY(intersect_pt_new[2][1] - intersect_pt_old[2][1]);
		break;
	case ClickType::Z_AXIS_CLICK:
		ans.setZ(intersect_pt_new[0][2] - intersect_pt_old[0][2]);
		break;
	default:
		break;
	}
	x_ += ans.x();
	y_ += ans.y();
	z_ += ans.z();
	std::cout << "dir: " << ans.x() << " " << ans.y() << " " << ans.z() << std::endl;
	return ans;
}

/*
* 返回对应物体的3D旋转变化 前三维为轴，第四维为角度
*/
QVector4D CoordinateAxis::RotateEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP)
{
	std::cout << "rotate event" << std::endl;
	QMatrix4x4 inv_matrix = modelDraw.modelmatrix.inverted() * modelDraw.viewmatrix.inverted() * modelDraw.projectionmatrix.inverted();
	vector<QVector3D> intersect_pt_old = IntersectPlanePoint(inv_matrix, oldP.x(), oldP.y());
	vector<QVector3D> intersect_pt_new = IntersectPlanePoint(inv_matrix, newP.x(), newP.y());
	QVector3D old_vector, new_vector, judge_vector;
	switch (click_type_) {
	case ClickType::X_CIRCLE_CLICK:
		old_vector = intersect_pt_old[0];
		new_vector = intersect_pt_new[0];
		judge_vector = QVector3D(1, 0, 0);
		break;
	case ClickType::Y_CIRCLE_CLICK:
		old_vector = intersect_pt_old[1];
		new_vector = intersect_pt_new[1];
		judge_vector = QVector3D(0, 1, 0);
		break;
	case ClickType::Z_CIRCLE_CLICK:
		old_vector = intersect_pt_old[2];
		new_vector = intersect_pt_new[2];
		judge_vector = QVector3D(0, 0, 1);
		break;
	default:
		break;
	}
	old_vector.normalize();
	new_vector.normalize();
	std::cout << "old: " << old_vector.x() << " " << old_vector.y() << " " << old_vector.z() << std::endl;
	std::cout << "new: " << new_vector.x() << " " << new_vector.y() << " " << new_vector.z() << std::endl;
	double angle = acos(QVector3D::dotProduct(old_vector, new_vector));
	angle = angle / PI * 180.0;
	if (QVector3D::dotProduct(QVector3D::crossProduct(old_vector, new_vector), judge_vector) < 0)
		angle = -angle;
	QVector4D ans = QVector4D(judge_vector, angle);
	std::cout << "quat: " << ans.x() << " " << ans.y() << " " << ans.z() << " " << ans.w() << std::endl;
	return ans;
}

/*
* 返回坐标轴的3D位置变化的两倍，伸缩变化为两倍坐标轴位置变化
*/
QVector3D CoordinateAxis::FlexibleEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP) {
	std::cout << "flexible event" << std::endl;
	QMatrix4x4 inv_matrix = modelDraw.modelmatrix.inverted() * modelDraw.viewmatrix.inverted() * modelDraw.projectionmatrix.inverted();
	vector<QVector3D> intersect_pt_old = IntersectPlanePoint(inv_matrix, oldP.x(), oldP.y());
	vector<QVector3D> intersect_pt_new = IntersectPlanePoint(inv_matrix, newP.x(), newP.y());
	QVector3D ans = QVector3D(0, 0, 0);
	switch (click_type_) {
	case ClickType::X_AXIS_CLICK:
		ans.setX(intersect_pt_new[1][0] - intersect_pt_old[1][0]);
		break;
	case ClickType::Y_AXIS_CLICK:
		ans.setY(intersect_pt_new[2][1] - intersect_pt_old[2][1]);
		break;
	case ClickType::Z_AXIS_CLICK:
		ans.setZ(intersect_pt_new[0][2] - intersect_pt_old[0][2]);
		break;
	default:
		break;
	}
	//x_ += ans.x()/2.0;
	//y_ += ans.y()/2.0;
	//z_ += ans.z()/2.0;
	std::cout << "dir: " << ans.x() << " " << ans.y() << " " << ans.z() << std::endl;
	return ans;
}

bool CoordinateAxis::PickEvent(ModelDraw& modelDraw, double winX, double winY) {
	if (mode_type_ == ModeType::HIDE) return false;
	QMatrix4x4 inv_matrix = modelDraw.modelmatrix.inverted() * modelDraw.viewmatrix.inverted() * modelDraw.projectionmatrix.inverted();
	vector<QVector3D> intersect_pt = IntersectPlanePoint(inv_matrix, winX, winY);

	double min_aix[3] = { pick_threshold, pick_threshold, pick_threshold };		// 与轴的最小长度
	double min_circle[3] = { pick_threshold, pick_threshold, pick_threshold };  // 与圈的最小长度
	double min_cube[3] = { pick_threshold, pick_threshold, pick_threshold };	// 与体的最小长度
	// whether intersect_pt on any plane (没搞懂games104的想法，自己写了个，有问题再改)
	for (int i = 0; i < 3; i++) {
		int j = (i + 1) % 3;
		int k = (i + 2) % 3;
		double dx = fabs(intersect_pt[i][j]);
		double dy = fabs(intersect_pt[i][k]);
		if (dx < pick_threshold && dy < len) 
			min_aix[k] = min_aix[k] < dx ? min_aix[k] : dx;
		if (dy < pick_threshold && dx < len)
			min_aix[j] = min_aix[j] < dy ? min_aix[j] : dy;
		if (mode_type_ == ModeType::TRANSFROM) {
			double dxy = sqrt(dx * dx + dy * dy);
			if (fabs(dxy-len) < pick_threshold)
				min_circle[i] = min_circle[j] < fabs(dxy - len) ? min_circle[j] : fabs(dxy - len);
		}
		else if (mode_type_ == ModeType::EDIT) {
			if (dx < pick_threshold + clen && dy >= len && dy < len + 2 * clen)
				min_cube[k] = min_cube[k] < dx-clen ? min_aix[k] : dx-clen;
			if (dy < pick_threshold + clen && dx >= len && dx < len + 2 * clen)
				min_cube[j] = min_cube[j] < dy-clen ? min_aix[j] : dy-clen;
		}
	}
	if (mode_type_ == ModeType::TRANSFROM) {
		double min_dist = pick_threshold;
		int m_selected = (int)ClickType::NO_CLICK;
		for (int i = 0; i < 3; i++) {
			if (min_dist > min_aix[i]) {
				min_dist = min_aix[i];
				m_selected = i;
			}
			if (min_dist > min_circle[i]) {
				min_dist = min_circle[i];
				m_selected = i+3;
			}
		}
		//std::cout << m_selected << std::endl;
		click_type_ = (ClickType)m_selected;
		if (click_type_ == ClickType::NO_CLICK)
			return false;
		else
			return true;
	}
	else if (mode_type_ == ModeType::EDIT) {
		double min_dist = pick_threshold;
		int m_selected = (int)ClickType::NO_CLICK;
		for (int i = 0; i < 3; i++) {
			if (min_dist > min_aix[i]) {
				min_dist = min_aix[i];
				m_selected = i;
			}
			if (min_dist > min_cube[i]) {
				min_dist = min_cube[i];
				m_selected = i + 6;
			}
		}
		//std::cout << m_selected << std::endl;
		click_type_ = (ClickType)m_selected;
		if (click_type_ == ClickType::NO_CLICK)
			return false;
		else
			return true;
	}
}

double CoordinateAxis::IntersectPlaneRay(QVector3D normal, QVector3D origin, QVector3D dir) {
	double deno = QVector3D::dotProduct(normal, dir);
	if (fabs(deno) < 0.0001)
		deno = 0.0001;
	return -QVector3D::dotProduct(normal, origin) / deno;
}

vector<QVector3D> CoordinateAxis::IntersectPlanePoint(QMatrix4x4 inv_matrix, double winX, double winY) {
	QVector4D point(winX, winY, 1.0, 1.0);			// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1.0, 1.0);		// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = inv_matrix * point;
	QVector4D tpointEnd = inv_matrix * pointEnd;
	QVector3D world_ray_origin(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D world_ray_to(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	QVector3D world_coordiante_orgin(x_, y_, z_);

	QVector3D local_ray_origin = world_ray_origin - world_coordiante_orgin;
	QVector3D local_ray_dir = world_ray_origin - world_ray_to;
	local_ray_dir.normalize();
	QVector3D plane_normals[3] = {
		QVector3D(1, 0, 0), QVector3D(0, 1, 0), QVector3D(0, 0, 1)
	};
	double plane_view_depth[3] = {
		IntersectPlaneRay(plane_normals[0], local_ray_origin, local_ray_dir),
		IntersectPlaneRay(plane_normals[1], local_ray_origin, local_ray_dir),
		IntersectPlaneRay(plane_normals[2], local_ray_origin, local_ray_dir)
	};
	vector<QVector3D> intersect_pt;
	intersect_pt.push_back(local_ray_origin + plane_view_depth[0] * local_ray_dir); // yoz
	intersect_pt.push_back(local_ray_origin + plane_view_depth[1] * local_ray_dir); // xoz
	intersect_pt.push_back(local_ray_origin + plane_view_depth[2] * local_ray_dir);	// xoy
	
	t_point1 = intersect_pt[2];
	//std::cout << "xoy " << t_point1.x() << " " << t_point1.y() << " " << t_point1.z() << std::endl;
	t_point2 = intersect_pt[1];
	//std::cout << "xoz " << t_point2.x() << " " << t_point2.y() << " " << t_point1.z() << std::endl;
	return intersect_pt;
}
