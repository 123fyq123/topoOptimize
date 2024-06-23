#pragma once
/*
OpenGL绘图窗口类，包括图形的渲染和视图显示
*/
#include <QtWidgets/qopenglwidget.h>
#include "Tools/TinyVector.h"
#include <QtGUI/qopenglfunctions.h>
#include <QOpenGLFunctions_3_3_Core>
#include <QOpenGLShaderProgram>
#include <GL/freeglut.h>
#include <Eigen/Dense>
#include <Eigen/Core>
//class QOpenGLTexture;


class ViewOpenGL :public QOpenGLWidget {
	Q_OBJECT
public:
	ViewOpenGL(QWidget* _parent = 0);
	virtual ~ViewOpenGL(void);
private:
	void Init(void);// OpenGL初始化
public:
	QSize minimumSizeHint(void) const override;
	QSize sizeHint(void) const override;
	/* 视角半径 */
	const double& Radius(void) const;
	/* 视角中心点 */
	const Vector3d& Center(void) const;
	//const double *GetModelviewMatrix (void) const;
	const QMatrix4x4* GetModelviewMatrix1(void) const;
	void ResetModelviewMatrix(void);
	void CopyModelViewMatrix(void);
	void LoadCopyModelViewMatrix(void);
	//const double *GetProjectionMatrix (void) const;
	const QMatrix4x4* GetProjectionMatrix1(void) const;

	enum ProjectionMode {
		PERSPECTIVE, ORTHOGRAPHIC
	};
	void SetProjectionMode(const ProjectionMode& pm);
	const ProjectionMode& GetProjectionMode(void) const;
	/* 绘制样式---点、线框、线、填充 */
	enum DrawMode {// BOUNDLINE 暂未实现 
		POINTS, BOUNDLINE, FLATLINES, FILL, FLATFILL
	};
	enum ColorMode {
		LIGHT, NORMAL, DIY
	};
	enum DrawModule {
		EMPTY, SURFACEMESH, HEXMESH, TETMESH, HYBRIDMESH, BSPLINECURVE, BSPLINESURFACE, SKELMESH, SKEL, SPARSEOPTIMIZATION, TRIMESH, BSPLINESOLID, BSPLINESOLIDS
	};
	// TEMP-LEE
	/* 鼠标按压模式 */
	Qt::MouseButton mousemode;

	void SetDrawMode(const DrawMode& dm);
	void SetColorMode(const ColorMode& dm);
	void SetDrawModule(const DrawModule& module);

	const DrawMode& GetDrawMode(void) const;
	const ColorMode& GetColorMode(void) const;
	const DrawModule& GetDrawModule(void)const;

protected:

	/* OpenGL绘制重写函数 */
	void initializeGL(void) override;
	void resizeGL(int w, int h) override;
	void paintGL(void) override;
	/* 设置绘制样式 */
	virtual void DrawScene();
	/* 鼠标键盘事件 */
	virtual void mousePressEvent(QMouseEvent*) override;
	virtual void mouseMoveEvent(QMouseEvent*) override;
	virtual void mouseReleaseEvent(QMouseEvent*) override;
	virtual void wheelEvent(QWheelEvent*) override;
	virtual void keyPressEvent(QKeyEvent*) override;
	virtual void keyReleaseEvent(QKeyEvent*) override;
protected:
	/* 旋转、缩放、平移等变换 */
	void Translation(const QPoint& p);
	void Translate(const Vector3d& trans);
	void Rotation(const QPoint& p);
	void Rotate(const Vector3d& axis, const double& angle);
	bool MapToSphere(const QPoint& point, Vector3d& result);
	void UpdateProjectionMatrix(bool is_load = true);
	void UpdatePickMoveMatrix(const QPoint& p);

public:
	/* 设置OpenGL绘制窗口的显示区域---显示中心和显示半径 */
	void SetScenePosition(const Vector3d& c, const double& r);
	void ViewAll(void);//更新模型矩阵
protected:
	DrawMode drawmode;// 绘制样式枚举值 
	ColorMode colormode;// 渲染样式枚举值 
	DrawModule drawmodule;// 绘制模块 
	ProjectionMode projectionmode;// 变换矩阵 
	double windowleft;
	double windowright;
	double windowtop;
	double windowbottom;

	/* 视角中心和半径值，可以通过鼠标操作改变 */
	Vector3d center;
	double radius;
	// 声明时已为单位矩阵
	QMatrix4x4 projectionmatrix;
	QMatrix4x4 modelmatrix;
	QMatrix4x4 viewmatrix;

	QMatrix4x4 tempmatrix;
	QMatrix4x4 tempmatrix2;
	QMatrix4x4 rotatematrix;
	QMatrix4x4 pickmovematrix;
	/* 鼠标移动时记录坐标点的位置 */
	QPoint lastpoint2;
	Vector3d lastpoint3;
	bool lastpointok;

	QOpenGLFunctions_3_3_Core* core;	// OpenGL核心态
	QOpenGLShaderProgram shaderProgram;

protected:
	static const double trackballradius;
};