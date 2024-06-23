#include "interpolationCurve.h"
#include <QFileDialog>
#include <fstream>
#include <QGLWidget>
#include <gl/freeglut.h>

#include <json/config.h>
#include <json/json.h>
/**************************************************
@brief   : 实现曲线的差值
@author  : lee
@step    : 1. json文件读取需要插值的控制点
           2. 调用算法实现样条曲线
		   3. 绘制样条曲线
@output  : none
@time    : none
**************************************************/

/**************************************************
@brief   : InterpolateCurve构造函数
@author  : lee
@input   : none
@output  : none
@time    : none
**************************************************/
InterpolateCurve::InterpolateCurve()  {
	flag = -1;// 没有任何的鼠标事件
	rotationX = 0;
	rotationY = 0;
	rotationZ = 0;
	p = 0;
	curvePointNum = 0;
	times = 1.0;
};



/**************************************************
@brief   : 析构函数释放相关资源
@author  : lee
@input   : none
@output  : none
@time    : none
**************************************************/
InterpolateCurve::~InterpolateCurve() {}

/**************************************************
@brief   : 全局曲线插值	ALGORITHM A9.1
@author  : lee
@input   : Q 要拟合的点
		   n // n+1个要拟合的点
		   r 坐标的维度
		   p 拟合的曲线的阶次
@output  : m // U_0 - U_m 中的m
		   U //
		   P // ctrlPoints
@time    : none
**************************************************/
void InterpolateCurve::globalCurveInterp(const std::vector<BMesh::Point3D> &Q,
	int p, int &m, std::vector<double> &U, std::vector<BMesh::Point3D> &ctlPts) {
	int n = Q.size() - 1;
	m = n + p + 1;
	// 计算u_k 通过 Eq.9.5
	double sumChordLength = 0;
	for (int i = 1; i < n + 1; i++) {
		sumChordLength += chordLength(Q[i], Q[i - 1]);
	}
	std::vector<double> u_k;// Q点对应的u
	u_k.resize(Q.size());
	u_k[0] = 0.0;
	for (int i = 1; i < n + 1; i++) {
		u_k[i] = u_k[i - 1] + chordLength(Q[i], Q[i - 1]) / sumChordLength;
	}//插值点对应的u_k???   然后形成U？？  Eq.(9.5)
	// 计算U 通过 Eq. 9.8
	U.resize(m + 1);
	for (int i = 0; i <= p; i++) {
		U[i] = 0;
	}
	for (int j = 1; j <= n - p; j++) {
		double tmp = 0;
		for (int i = j; i <= j + p - 1; i++) {
			tmp += u_k[i];
		}
		U[j + p] = 1.0 / p * tmp;
	}
	for (int i = m - p; i <= m; i++) {
		U[i] = 1;
	}
	Eigen::MatrixXd N = Eigen::MatrixXd::Zero(n + 1, n + 1);
	
	for (int i = 0; i <= n; i++) {
		// 构造N矩阵
		int span = findSpan(n, p, u_k[i], U);
		std::vector<double> base;
		basisFuns(span, u_k[i], p, U, base);// N[i][span - p]
		for (int j = 0; j < base.size(); j++) {
			N(i, span - p + j) = base[j];
		}
	}
	// test N
	/*for (int i = 0; i < n + 1; i++) {
		for (int j = 0; j < n + 1; j++) {
			std::cout << N(i, j) << " ";
		}
		std::cout << std::endl;
	}*/
	// LUDecomposition(A,n+1,p-1):  LU分解  NP=F  求P
	Eigen::MatrixXd F = Eigen::MatrixXd::Zero(n+1, 3);
	for (int i = 0; i < n+1; ++i)
	{
		F(i, 0) = Q[i].x;
		F(i, 1) = Q[i].y;
		F(i, 2) = Q[i].z;
	}
	Eigen::MatrixXd P(n + 1, 3);
	P = N.fullPivLu().solve(F);
#ifdef __TT__
	std::cout << "N--------------" << std::endl << N << std::endl;
	std::cout << "F--------------" << std::endl << F << std::endl;
	std::cout << "P--------------" << std::endl << P << std::endl;
#endif
	ctlPts.clear();
	for (int i = 0; i < n + 1; i++) {
		BMesh::Point3D point(P(i, 0), P(i, 1), P(i, 2));
		ctlPts.push_back(point);
	}
}


/**************************************************
@brief   : 计算弦长 chord length 单纯两点之间的距离？？
@author  : lee
@input   : vecFitPoints 要拟合的点
@output  : none
@time    : none
**************************************************/
double InterpolateCurve::chordLength(const BMesh::Point3D& a, const BMesh::Point3D&b) {
	double length = sqrt((a.x - b.x)*(a.x - b.x) +
		(a.y - b.y)*(a.y - b.y) + (a.z - b.z)*(a.z - b.z));
	return length;
}


/**************************************************
@brief   : 找到 u 所属于的区间序列
@author  : lee
@input   : n n = m - p -1
		   p 基函数的阶
		   u 递增的一般从0到1上的点
		   U 节点向量 U_0 - U_m
@output  : 节点区间序列
@time    : none
**************************************************/
int InterpolateCurve::findSpan(int n, int p, double u, const std::vector<double> &U) {
	if (u >= U[n + 1]) return n;
	for (int i = p; i < n + 1; i++) {
		if (U[i] <= u && u < U[i + 1]) {
			return i;
		}
	}
}

/**************************************************
@brief   : 计算非空的基函数
		   第一次，p=0时，只有N[0]=1.0,即N0,0=1.0;p=1时，求N[0],N[1],即N0,1和N1,1;p=2时，求N[0],N[1],N[2],即N0,2  N1,2和N2,2
@author  : lee
@input   ：i 节点向量的第i个序列
		   u 生成曲线上的点
		   order 阶次
		   U 节点向量
@output  ：base 从第0阶到后面的阶辐射的三角形所对应的阶函数的值，一整套
@time    : none
**************************************************/
void InterpolateCurve::basisFuns(int i, double u, int order, const std::vector<double> &U, std::vector<double> &base) {
	base.resize(order + 1);
	base[0] = 1.0;
	std::vector<double> left, right;
	left.resize(order + 1);
	right.resize(order + 1);
	for (int j = 1; j <= order; j++) {
		left[j] = u - U[i + 1 - j];
		right[j] = U[i + j] - u;
		double saved = 0.0;
		for (int r = 0; r < j; r++) {
			double temp = base[r] / (right[r + 1] + left[j - r]);
			base[r] = saved + right[r + 1] * temp;
			saved = left[j - r] * temp;
		}
		base[j] = saved;
	}
}


/**************************************************
@brief   : 生成曲线上的点
		   关键公式  m = n + p + 1
@author  : lee
@input   : n+1 控制点的个数
		   p 阶次
		   m+1 节点的个数
		   P  控制点数组
		   u 从小到大，从而绘制出整个曲线
		   C 生成的点的坐标
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::curvePoint(int n, int p, const std::vector<double>&U, std::vector<BMesh::Point3D>&P, double u, BMesh::Point3D&C) {
	int span = findSpan(n, p, u, U);
	std::vector<double> base;
	basisFuns(span, u, p, U, base);
	double x, y, z;
	x = y = z = 0;
	for (int i = 0; i <= p; i++) {
		x = x + P[span - p + i].x * base[i];
		y = y + base[i] * P[span - p + i].y;
		z = z + base[i] * P[span - p + i].z;
	}
	BMesh::Point3D temp;
	temp.x = x;
	temp.y = y;
	temp.z = z;
	C = temp;
}


/**************************************************
@brief   : 简单的讲生成的点加入到曲线上
@author  : lee
@input   ：str  文件路径
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::genPoints(std::string str) {
	jsonReader(str);
	std::cout << "[DEBUG]genPoint" << std::endl;
	if (name == "interpolateCurve") {
		int m;
		std::vector<double> U;
		globalCurveInterp(Q, p, m, U, ctlPts);
		//BSplineCurveAbstract bs(p, ctlPts.size() - 1, U, ctlPts);

		for (int i = 0; i <= curvePointNum; i++) {
			double u = i * 1.0 * (U[U.size() - 1] - U[0]) / curvePointNum;//均匀分成多少个点再所有的节点上
			BMesh::Point3D C;
			curvePoint(ctlPts.size() - 1, p, U, ctlPts, u, C);
			curvePoints.push_back(C);
		}
	}
	else {
		int size = ps.size();
		for (int i = 0; i < size; i++) {
			int m;
			std::vector<double> U;
			std::vector<BMesh::Point3D> ctlPtsTmp;
			globalCurveInterp(Qs[i], ps[i], m, U, ctlPtsTmp);
			//BSplineCurveAbstract bs(ps[i], ctlPtsTmp.size() - 1, U, ctlPtsTmp);
			std::vector<BMesh::Point3D> curvePointsTmp;
			for (int j = 0; j <= curvePointNums[i]; j++) {
				double u = j * 1.0 * (U[U.size() - 1] - U[0]) / curvePointNums[i];//均匀分成多少个点再所有的节点上
				BMesh::Point3D C;
				curvePoint(ctlPtsTmp.size() - 1, ps[i], U, ctlPtsTmp, u, C);
				
				curvePointsTmp.push_back(C);
			}
			curvePointss.push_back(curvePointsTmp);
			ctlPtss.push_back(ctlPtsTmp);
		}
	}
	// check ps endpoint in curvePointss ??
	/*TEST*/
	/*
	int size = Qs.size();
	for (int i = 0; i < size; i++) {
		Point3d a = Qs[i][0];
		Point3d b = Qs[i].back();
		std::cout << "-- TIME -- \n";
		for (int j = 0; j < curvePointss[i].size(); j++) {
			if (a == curvePointss[i][j]) {
				std::cout << "a:\n";
				a.print();
				curvePointss[i][j].print();
				std::cout << "B-:\n";
				b.print();
				curvePointss[i].back().print();
			}
			else if (b == curvePointss[i][j]) {
				std::cout << "b:\n";
				b.print();
				curvePointss[i][j].print();
			}
		}
	}
	*/
}



/**************************************************
@brief   : 从文件中读取需要经过的点
@author  : lee
@input   ：文件名(完整路径)/或空 方便调试
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::jsonReader(std::string fileName) {
	if (fileName == "") {
		QString QfileName = QFileDialog::getOpenFileName(NULL,
			QObject::tr("Open json file"),
			QObject::tr("./config"),
			QObject::tr("mesh (*.json);;"
				"All Files (*)"));
		fileName = QfileName.toStdString();
	}
	Json::Value root;
	Json::Reader reader;
	std::ifstream ifs(fileName.c_str());//open file example.json

	if (!reader.parse(ifs, root)) {
		// fail to parse
		std::cout << "fail to parse json \n";
	}
	else {
		// success
		std::cout << "[DEBUG] config file have found!!\n";
		std::cout << "name: " << root["name"].asString() << std::endl;
		name = root["name"].asString();
		if (name == "interpolateCurve") {
			p = root["p"].asInt();
			curvePointNum = root["curvePointNum"].asInt();

			Json::Value points = root["points"];
			for (int i = 0; i < points.size(); i++) {
				BMesh::Point3D p(0, 0, 0);
				p.x = (points[i][0].asDouble());
				p.y = (points[i][1].asDouble());
				p.z = (points[i][2].asDouble());
				Q.push_back(p);
			}
		}
		else if(name == "interpolateCurves"){
			Json::Value lines = root["lines"];
			for (int i = 0; i < lines.size(); i++) {
				curvePointNums.push_back(lines[i]["curvePointNum"].asInt());
				ps.push_back(lines[i]["p"].asInt());
				std::vector<BMesh::Point3D> Qindex;
				Json::Value points = lines[i]["points"];
				for (int j = 0; j < points.size(); j++) {
					BMesh::Point3D p(0, 0, 0);
					p.x = (points[j][0].asDouble());
					p.y = (points[j][1].asDouble());
					p.z = (points[j][2].asDouble());
					Qindex.push_back(p);
				}
				Qs.push_back(Qindex);
			}
		}
	}
}


/**************************************************
@brief   : 初始化opengl环境
@author  : lee
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::initializeGL() {
	glShadeModel(GL_SMOOTH);
	glClearColor(0.1, 0.1, 0.4, 1.0);
	glClearDepth(1.0);
	glEnable(GL_DEPTH_TEST);
	glDepthFunc(GL_LEQUAL);
	glHint(GL_PERSPECTIVE_CORRECTION_HINT, GL_NICEST);
}


/**************************************************
@brief   : 改变窗口的大小
@author  : lee
@input   ：width 宽度
		   height 高度
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::resizeGL(int width, int height)
{
	// 防止height为0
	if (height == 0) {
		height = 1;
	}
	// 重置当前的视口
	glViewport(0, 0, (GLint)width, (GLint)height);
	// 选择投影矩阵
	glMatrixMode(GL_PROJECTION);
	// 重置观察矩阵/投影矩阵 当调用次函数，实际将当前点移到了屏幕中心
	glLoadIdentity();
	// 建立透视投影矩阵,需要<GL/glu.h>头文件
	gluPerspective(45.0, (GLdouble)width / (GLdouble)height, 0.1, 100.0);//角度，长宽比，远近
	// 选择模型观察矩阵
	glMatrixMode(GL_MODELVIEW);
	// 重置观察矩阵/投影矩阵
	glLoadIdentity();
}


/**************************************************
@brief   : 自动和被动 调用更新整个图形
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::paintGL() {
	//genPoints("normal");
	glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);// 清除屏幕和深度缓存
	glMatrixMode(GL_MODELVIEW);
	glLoadIdentity();
	glTranslatef(1.0, -4, -25.0);//平移

	glRotated(rotationX, 1.0, 0.0, 0.0);
	glRotated(rotationY, 0.0, 1.0, 0.0);
	glRotated(rotationZ, 0.0, 0.0, 1.0);//旋转
	glScalef(times, times, times);//缩放
	// 绘制坐标轴
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 0.0, 1.0); // 蓝色表示x
	glVertex3d(0, 0, 0);
	glVertex3d(1000, 0, 0);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glColor3f(0.0, 1.0, 0.0);// 绿色  表示 y
	glVertex3d(0, 0, 0);
	glVertex3d(0, 1000, 0);
	glEnd();
	glBegin(GL_LINE_STRIP);
	glColor3f(1.0, 0, 0.0);// 红色表示的是 Z
	glVertex3d(0, 0, 0);
	glVertex3d(0, 0, 1000);
	glEnd();
	// 画控制点
	
	glPointSize(8.0);
	glBegin(GL_POINTS);
	glColor3f(1.0, 0.0, 0.0);
	if (name == "interpolateCurve") {
		for (int i = 0; i < ctlPts.size(); i++) {
			glVertex3d(ctlPts[i].x, ctlPts[i].y, ctlPts[i].z);
		}
	}
	else {
		for (int i = 0; i < ctlPtss.size(); i++) {
			for (int j = 0; j < ctlPtss[i].size(); j++) {
				glVertex3d(ctlPtss[i][j].x, ctlPtss[i][j].y, ctlPtss[i][j].z);
			}
		}
	}
	

	glEnd();

	// 画出控制多边形
	
	glColor3f(0.5, 1.0, 0.3);
	if (name == "interpolateCurve") {
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < ctlPts.size(); i++) {
			glVertex3d(ctlPts[i].x, ctlPts[i].y, ctlPts[i].z);
		}
		glEnd();
	}
	else {
		for (int i = 0; i < ctlPtss.size(); i++) {
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < ctlPtss[i].size(); j++) {
				glVertex3d(ctlPtss[i][j].x, ctlPtss[i][j].y, ctlPtss[i][j].z);
			}
			glEnd();
		}
	}

	// 画出要拟合的点
	
	glColor3f(1.0, 1.0, 0.0);
	if (name == "interpolateCurve") {
		glBegin(GL_POINTS);
		for (int i = 0; i < Q.size(); i++) {
			glVertex3d(Q[i].x, Q[i].y, Q[i].z);
		}
		glEnd();
	}
	else {
		for (int i = 0; i < Qs.size(); i++) {
			glBegin(GL_POINTS);
			for (int j = 0; j < Qs[i].size(); j++) {
				glVertex3d(Qs[i][j].x, Qs[i][j].y, Qs[i][j].z);
			}
			glEnd();
		}
	}
	

	// 画出曲线
	
	glColor3f(1.0, 1.0, 1.0);
	if (name == "interpolateCurve") {
		glBegin(GL_LINE_STRIP);
		for (int i = 0; i < curvePoints.size(); i++) {
			glVertex3f(curvePoints[i].x, curvePoints[i].y, curvePoints[i].z);
		}
		glEnd();
	}
	else {
		
		for (int i = 0; i < curvePointss.size(); i++) {
			glBegin(GL_LINE_STRIP);
			for (int j = 0; j < curvePointss[i].size(); j++) {
				glVertex3f(curvePointss[i][j].x, curvePointss[i][j].y, curvePointss[i][j].z);
			}
			glEnd();
		}
		
	}
}



/**************************************************
@brief   : 按键按下事件
@author  : lee
@input   ：e  事件
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::mousePressEvent(QMouseEvent *e) {
	lastPos = e->pos();
	flag = 1;
}


/**************************************************
@brief   : none
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::mouseReleaseEvent(QMouseEvent *e) {
	flag = -1;
}


/**************************************************
@brief   : 鼠标滚轮事件
@author  : lee
@input   ：event
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::wheelEvent(QWheelEvent *event) {
	if (event->delta() > 0) {// 当滚轮远离使用者时
		times += 0.1f;
		//update();
	}
	else {//当滚轮向使用者方向旋转时
		times -= 0.1f;
		//update();
	}
}



/**************************************************
@brief   : none
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::mouseMoveEvent(QMouseEvent *e) {
	//if (flag) {
	//	GLdouble dx = GLdouble(e->x() - lastPos.x()) / width();//QWidght 的 宽度
	//	GLdouble dy = GLdouble(e->y() - lastPos.y()) / height();
	//	if (e->buttons() & Qt::LeftButton) {
	//		rotationX -= 180 * dy;
	//		rotationY -= 180 * dx;
	//		update();
	//	}
	//	else if (e->buttons() & Qt::RightButton) {
	//		rotationX -= 180 * dy;
	//		rotationZ -= 180 * dx;
	//		update();
	//	}
	//	lastPos = e->pos();
	//}
}

/**************************************************
@brief   : read off file then skel to simple
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::readOff() {
	std::string filename;
	QString QfileName = QFileDialog::getOpenFileName(NULL,
		QObject::tr("Open off file"),
		QObject::tr("./config"),
		QObject::tr("OFF (*.off);;"
			"All Files (*)"));
	filename = QfileName.toStdString();

	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::read_mesh(mesh, filename))
		{
			std::cerr << "Cannot write mesh to file " << filename << std::endl;
		}
		else {
			std::cout << "[DEBUG] read mesh success, lucky Guy\n";
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}


void InterpolateCurve::readOff(std::string filename) {
	// write mesh to output.obj
	try
	{
		if (!OpenMesh::IO::read_mesh(mesh, filename))
		{
			std::cerr << "Cannot write mesh to file " << filename << std::endl;
		}
		else {
			std::cout << "[DEBUG] read mesh success, lucky Guy\n";
		}
	}
	catch (std::exception& x)
	{
		std::cerr << x.what() << std::endl;
	}
}


/**************************************************
@brief   : simple SKELD
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/


void InterpolateCurve::simpleSkel() {
	Qos.resize(Qs.size());
	for (int i = 0; i < Qs.size(); i++) {
		for (int j = 0; j < Qs[i].size(); j++) {
			if (j == 0 || j == Qs[i].size() - 1) {
				SkelBoll a;
				a.p = Qs[i][j];
				MyMesh::Point opa(a.p.x, a.p.y, a.p.z);
				double maxR = INT_MAX;
				for (auto it : mesh.all_vertices()) {
					double len = (mesh.point(it) - opa).length();
					if (maxR > len) {
						maxR = len;
					}
				}
				a.r = maxR;
				if (j == 0)
					a.index = 0;
				else
					a.index = curvePointss[i].size() - 1;
				Qos[i].push_back(a);
			}
		}
		bool check = true;
		while (check) {
			std::vector<SkelBoll> pss;
			for (int j = 0; j < Qos[i].size()-1; j++) {
				int indexCurve = (Qos[i][j].index + Qos[i][j + 1].index) / 2;
				SkelBoll a;
				a.p = curvePointss[i][indexCurve];
				a.index = indexCurve;
				double maxR = INT_MAX;
				MyMesh::Point opa(a.p.x, a.p.y, a.p.z);
				for (auto it : mesh.all_vertices()) {
					double len = (mesh.point(it) - opa).length();
					if (maxR > len) {
						maxR = len;
					}
				}
				a.r = maxR;
				MyMesh::Point opb(Qos[i][j].p.x, Qos[i][j].p.y, Qos[i][j].p.z);
				MyMesh::Point opc(Qos[i][j+1].p.x, Qos[i][j+1].p.y, Qos[i][j+1].p.z);
				if (a.r + Qos[i][j].r > (opb - opa).length() || a.r + Qos[i][j+1].r > (opc - opa).length()) {
					check = false;
				}
				pss.push_back(a);
			}
			if (check) {
				int size = Qos[i].size();
				int indexx = pss.size()- 1;
				for (int j = size - 1; j >= 1; j--) {
					Qos[i].insert(Qos[i].begin() + j, pss[indexx]);
					indexx--;
				}
			}
		}
	}
}


/**************************************************
@brief   : export SKEL
@author  : none
@input   ：none
@output  ：none
@time    : none
**************************************************/
void InterpolateCurve::exportSkel() {

	std::vector<SkelBoll> rlt;
	std::set<SkelBoll, MyCmp> s;
	std::map<SkelBoll, int, MyCmp> m;
	std::vector<std::vector<int>> childrens;
	for (int i = 0; i < Qos.size(); i++) {
		for (int j = 0; j < Qos[i].size(); j++) {
			if (m.count(Qos[i][j]) == 0) {
				s.insert(Qos[i][j]);
				rlt.push_back(Qos[i][j]);
				m[Qos[i][j]] = rlt.size() - 1;
			}
		}
	}
	childrens.resize(rlt.size());
	for (int i = 0; i < Qos.size(); i++) {
		for (int j = 0; j < Qos[i].size(); j++) {
			if (j == 0) {
				continue;
			}
			int index = m[Qos[i][j]];
			childrens[index].push_back(m[Qos[i][j - 1]]);
			childrens[m[Qos[i][j - 1]]].push_back(index);
		}
	}

	std::ofstream off = std::ofstream("result.skel", std::ios::out);
	if (!off.good()) {
		std::cout << "[Error] Could not open file " << name << " for writing!" << std::endl;
		off.close();
		return;
	}
	// Write header
	off << "ID Cx Cy Cz Length Width Height #NEIGHBORS NEIGHBORS_LIST" << std::endl;
	off << rlt.size() << std::endl;
	// renew id
	for (int i = 0; i < rlt.size(); i++) {
		rlt[i].id = i;
	}
	for (int i = 0; i < rlt.size(); i++) {
		off << rlt[i].id << " " << rlt[i].p.x << " " << rlt[i].p.y << " " << rlt[i].p.z
			<< " " << rlt[i].r*2/sqrt(3) << " " << rlt[i].r * 2 / sqrt(3) << " " << rlt[i].r * 2 / sqrt(3) << " " << childrens[i].size();
		for (int j = 0; j < childrens[i].size(); j++) {
			off << " " << rlt[childrens[i][j]].id;
		}
		off << std::endl;
	}
	off.close();
}
