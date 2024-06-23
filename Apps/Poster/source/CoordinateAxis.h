/*
* class: 坐标轴
* 实现坐标轴的绘制与点选交互
* 绘制与点选 需要写成友元函数（因为需要用到原来世界的矩阵，后面有想法给去耦合）
* 坐标轴实现两种模式：TRANSFORM（移动，旋转） EDIT（分裂，扩展，收缩）；以及 HIDE模式 使绘制与点选失效
* 只是实现坐标轴在触发功能下的变化，对应的其他变化不写在这里
* 
*/
#pragma once
#include <cmath>
#include <stack>
#include <QPoint>
#include <QVector>
#include <QMatrix4x4>

class ModelDraw;
class CoordinateAxis {
public:
	enum class ModeType {
		HIDE,
		TRANSFROM,
		EDIT,
	};
	enum class ClickType {
		X_AXIS_CLICK,
		Y_AXIS_CLICK,
		Z_AXIS_CLICK,
		X_CIRCLE_CLICK,
		Y_CIRCLE_CLICK,
		Z_CIRCLE_CLICK,
		X_CUBE_CLICK,
		Y_CUBE_CLICK,
		Z_CUBE_CLICK,
		NO_CLICK
	};
	// basis idea: get screen point, plot the coordinate axis in right position;
	CoordinateAxis() {
		mode_type_ = ModeType::HIDE;
		click_type_ = ClickType::NO_CLICK;
		x_ = 0, y_ = 0, z_ = 0;
		object_id_ = 0;
		model_matrix_ = QMatrix4x4();
		translate_matrix_ = QMatrix4x4();
		rotate_matrix_ = QMatrix4x4();
		t_point1 = QVector3D();
		t_point2 = QVector3D();
	};
	~CoordinateAxis() {};

	double get_len() { return len; };
	void set_len(double d) {
		len = d;  
		vlen = len * (1 - k);
		radios = len;
		step = radios / 25;
		clen = vlen / 2.0;
		pick_threshold = len * 0.05;
	}
	void get_orign(double& x, double& y, double& z) { x = x_; y = y_; z = z_; }
	void set_orign(double x, double y, double z) { x_ = x; y_ = y; z_ = z; }
	unsigned int get_object_id() { return object_id_; }
	void set_object_id(unsigned int id) { object_id_ = id; }
	void set_click_type(ClickType _clicktype) { click_type_ = _clicktype; }
	ClickType get_click_type() { return click_type_; }
	void set_mode_type(ModeType _modetype) { mode_type_ = _modetype; }
	ModeType get_mode_type() { return mode_type_; }
	
	// 绘制
	void Draw3DCoord(ModelDraw& modelDraw);										// 绘制坐标的功能
	void SetModelMatrix(ModelDraw& modelDraw, QMatrix4x4 matrix);				// 在调用了仿固定管线函数后，建议绘制前调用一次
	// 交互判断
	bool PickEvent(ModelDraw& modelDraw, double winX, double winY);
	QVector3D MoveEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);
	QVector4D RotateEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);
	QVector3D FlexibleEvent(ModelDraw& modelDraw, QVector2D oldP, QVector2D newP);
	
private:
	//bool MoveCoordianteAxis(double x, double y, double z);						// 移动轴 移动中心坐标
	//bool RotateCoordianteAxis(double angle, double x, double y, double z);		// 旋转轴 旋转后不变

	void DrawCone(QVector3D point, QVector3D dir ,double radios, double len);	// 绘制圆锥 dir: 目前只允许x,y,z轴
	void DrawCube();															// 绘制交互结点分裂的正方体
	double IntersectPlaneRay(QVector3D normal, QVector3D origin, QVector3D dir);
	std::vector<QVector3D> IntersectPlanePoint(QMatrix4x4 inv_matrix, double winX, double winY);
	// 仿固定管线
	void Translated(double x, double y, double z);
	void Rotatef(double angle, double x, double y, double z);
	void PushMatrix();
	void PopMatrix();

private:
	//const double PI = acos(-1.0);
	const double k = 0.85;				// 坐标轴轴身
	double len = 1.0f;					// 坐标轴的长度
	double vlen = len * (1 - k);		// 坐标轴轴帽
	double radios = len;				// 旋转弧线半径
	double step = radios / 25;			// 画弧半径
	double clen = vlen / 2.0;			// clip_cube 边长的一半 
	double pick_threshold = len * 0.05;	// 选取距离门槛

	ModeType mode_type_;
	ClickType click_type_;
	double x_, y_, z_;					// local coordinate orgin point
	unsigned int object_id_;			// face to object id
	
	std::stack <QMatrix4x4> temp_matrixs_;
	QMatrix4x4 translate_matrix_, rotate_matrix_;	// 想实现固定管线的glTranslate() glRotatef(),
	QMatrix4x4 model_matrix_;

	//test
	QVector3D t_point1, t_point2;
};