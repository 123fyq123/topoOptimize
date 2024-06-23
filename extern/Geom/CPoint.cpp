#include <cmath>
#include "CPoint.h"
#include "CAidFunction.h"
#include "CLine.h"
#include "CPlane.h"


CPoint::CPoint()
{
	m_dX = 0;
	m_dY = 0;
	m_dZ = 0;
	m_dR = 0;
	m_dTheta = 0;
	//m_dFi = 0;
	m_iTypeCoord = 1;
	m_iTypeAngle = 1;
}

CPoint::CPoint(double p1, double p2, double p3, int typeCoord, int typeAngle)
{
	m_iTypeCoord = typeCoord;
	m_iTypeAngle = typeAngle;
	if (typeCoord == 1)
	{
		m_dX = p1;
		m_dY = p2;
		m_dZ = p3;
		//极坐标初始化
		computePolar();
	}
	else if(typeCoord == 0)
	{
		m_dR = p1;
		m_dTheta = p2;
		m_dZ = p3;
		//m_dFi = p3;
		//笛卡尔坐标初始化
		computeCartesian();
	}
	else
	{
		//错误处理，待补充
	}

}


CPoint::CPoint(const CPoint & point)
{
	const_cast<CPoint &>(point).getXYZ(m_dX, m_dY, m_dZ);
	const_cast<CPoint &>(point).getPolar(m_dR, m_dTheta, m_dZ, m_iTypeAngle);
	m_iTypeCoord = point.getTypeCoord();
}

CPoint::CPoint(const gp_Pnt & pnt)
{
	m_dX = pnt.X();
	m_dY = pnt.Y();
	m_dZ = pnt.Z();
	m_iTypeCoord = 1;
	m_iTypeAngle = 1;
	computePolar();
}

CPoint::CPoint(const gp_Vec & vec)
{
	m_dX = vec.X();
	m_dY = vec.Y();
	m_dZ = vec.Z();
	m_iTypeCoord = 1;
	m_iTypeAngle = 1;
}


CPoint::~CPoint()
{
}

bool CPoint::setTypeCoord(int type)
{
	m_iTypeCoord = type;
	return true;
}

int CPoint::getTypeCoord() const
{
	return m_iTypeCoord;
}

bool CPoint::setX(double x)
{
	m_dX = x;
	computePolar();
	return true;
}

double CPoint::getX() const {
	return m_dX;
}

bool CPoint::setY(double y)
{
	m_dY = y;
	computePolar();
	return true;
}

double CPoint::getY() const {
	return m_dY;
}

bool CPoint::setZ(double z)
{
	m_dZ = z;
	//computePolar();
	return true;
}

double CPoint::getZ() const {
	return m_dZ;
}

bool CPoint::setXYZ(double x, double y, double z, int iTypeCoord)
{
	m_dX = x;
	m_dY = y;
	m_dZ = z;
	m_iTypeCoord = iTypeCoord;
	computePolar();
	return true;
}

bool CPoint::getXYZ(double & x, double & y, double & z)
{
	x = m_dX;
	y = m_dY;
	z = m_dZ;
	return true;
}

bool CPoint::setR(double r)
{
	m_dR = r;
	computeCartesian();
	return true;
}

double CPoint::getR(void)
{
	return m_dR;
}

int CPoint::getTypeAngle(void)
{
	return m_iTypeAngle;
}

bool CPoint::setTheta(double theta, int typeAngle)
{
	if (typeAngle != m_iTypeAngle)
	{
		//computeAngle();
		m_iTypeAngle = typeAngle;
	}
	m_dTheta = theta;
	computeCartesian();
	return true;
}

double CPoint::getTheta(void)
{
	return m_dTheta;
}

//bool CPoint::setFi(double fi, int typeAngle)
//{
//	if (typeAngle != m_iTypeAngle)
//	{
//		computeAngle();
//		m_iTypeAngle = typeAngle;
//		m_dFi = fi;
//
//	}
//	else
//	{
//		m_dFi = fi;
//	}
//	computeCartesian();
//	return true;
//}

//double CPoint::getFi(void) const
//{
//	return m_dFi;
//}

bool CPoint::setPolar(double r, double theta, double z, int typeAngle)
{
	if (typeAngle != m_iTypeAngle)
	{
		//computeAngle();
		m_iTypeAngle = typeAngle;
	}

	m_dR = r;
	m_dTheta = theta;
	m_dZ = z;
	computeCartesian();
	return true;
}

bool CPoint::getPolar(double & r, double & theta, double & z, int & typeAngle)
{
	r = m_dR;
	theta = m_dTheta;
	z = m_dZ;
	typeAngle = m_iTypeAngle;
	return true;
}

CPoint CPoint::operator+(const CPoint & point)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX + const_cast<CPoint &>(point).getX();
	double y = m_dY + const_cast<CPoint &>(point).getY();
	double z = m_dZ + const_cast<CPoint &>(point).getZ();
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}

CPoint CPoint::operator-(const CPoint & point)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX - const_cast<CPoint &>(point).getX();
	double y = m_dY - const_cast<CPoint &>(point).getY();
	double z = m_dZ - const_cast<CPoint &>(point).getZ();
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}

CPoint CPoint::operator*(const CPoint & point)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX * const_cast<CPoint &>(point).getX();
	double y = m_dY * const_cast<CPoint &>(point).getY();
	double z = m_dZ * const_cast<CPoint &>(point).getZ();
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}

CPoint CPoint::operator/(const CPoint & point)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX / const_cast<CPoint &>(point).getX();
	double y = m_dY / const_cast<CPoint &>(point).getY();
	double z = m_dZ / const_cast<CPoint &>(point).getZ();
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}

CPoint CPoint::operator+(const CVector & vec)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX + const_cast<CVector &>(vec).getX();
	double y = m_dY + const_cast<CVector &>(vec).getY();
	double z = m_dZ + const_cast<CVector &>(vec).getZ();
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}

CPoint CPoint::operator*(const double d)
{
	int typeAngle = m_iTypeAngle;
	double x = m_dX * d;
	double y = m_dY * d;
	double z = m_dZ * d;
	CPoint finalPoint(x, y, z, 1, typeAngle);
	return finalPoint;
}
bool CPoint::operator ==(const CPoint& point) {
	return getX() == point.m_dX && getY() == point.m_dY && getZ() == point.m_dZ;
}
double CPoint::norm(void) const
{
	double d = sqrt(m_dX * m_dX + m_dY * m_dY + m_dZ * m_dZ);
	return d;
}

bool CPoint::neg(void)
{
	m_dX = -m_dX;
	m_dY = -m_dY;
	m_dZ = -m_dZ;
	computePolar();
	return true;
}

void CPoint::normaliz()
{
	double d = norm();
	if (d < 1e-10)
	{
		//不做处理
	}
	else 
	{
		m_dX = m_dX / d;
		m_dY = m_dY / d;
		m_dZ = m_dZ / d;
	}
	
	
}

void CPoint::rotate(CGeom & line, double angle, int typeAngle, CGeomPoint &point)
{
	//绕直线旋转
	CPoint ps = *this;
	double tmpAngle;
	if (typeAngle == 1)
	{
		tmpAngle = angle * PI / 180.0;
	}
	else
	{
		tmpAngle = angle;
	}
		
	CLine line0 = dynamic_cast<CLine &>(line);
	CVector vec;
	line0.getNormal(vec);
	CPlane plane(ps, vec);
	double t = CAidFunction::intersectLinePlane(line0, plane);
	CPoint pnt = line0.getPointAtParam(t);
	CPoint pv1 = ps - pnt;
	double r = sqrt(CAidFunction::dot3D(pv1, pv1));
	pv1.normaliz();
	CPoint pv2;
	CAidFunction::crossProduct3D(pnt, pv1, pv2);
	pv2.normaliz();
	//CLine line1(&pnt, &pv1);
	CPoint p1 = pnt + pv1 * r * cos(tmpAngle);
	CPoint pp = p1 + pv2 * r * sin(tmpAngle);
	point.setXYZ(pp.getX(), pp.getY(), pp.getZ());
}

void CPoint::mirror(CGeom & plane, CGeomPoint & point)
{
	//绕平面对称
	CPoint ps = *this;
	CPlane plane0 = dynamic_cast<CPlane &>(plane);
	CVector normal;
	plane0.getNormal(normal);
	CLine line(ps, normal);
	double t = CAidFunction::intersectLinePlane(line, plane0);
	CPoint pnt = line.getPointAtParam(t);
	CPoint p = pnt + (pnt - *this);
	point.setXYZ(p.getX(), p.getY(), p.getZ());
}

 void CPoint::getPerpendicular(CGeomPoint &point)
{
	int typeAngle = m_iTypeAngle;
	if (fabs(m_dX) > 1e-8 || fabs(m_dY) > 1e-8)
	{
		//CPoint finalPoint(-m_dY, m_dX, 0, 1, typeAngle);
		point.setXYZ(-m_dY, m_dX, 0);
	}
	else
	{
		CPoint finalPoint(-m_dZ, 0, m_dX, 1, typeAngle);
		point.setXYZ(-m_dZ, 0, m_dX);
	}
		
}

 void CPoint::getPnt(gp_Pnt & pnt) const
 {
	 pnt.SetX(m_dX);
	 pnt.SetY(m_dY);
	 pnt.SetZ(m_dZ);
 }

 gp_Pnt CPoint::getPnt() const
 {
	 gp_Pnt pnt(m_dX, m_dY, m_dZ);
	 return pnt;
 }

void CPoint::computePolar()
{
	//根据笛卡尔坐标求极坐标
	double angle = atan2(m_dY, m_dX); //返回弧度
	if (m_iTypeAngle == 1)
	{
		m_dTheta = (double) angle * 180 / PI;
	}
	else
	{
		m_dTheta = angle;
	}

	m_dR = sqrt(m_dX * m_dX + m_dY * m_dY);
}

void CPoint::computeCartesian()
{
	//根据极坐标求笛卡尔坐标
	double angle = 0;
	if (m_iTypeAngle == 1)
	{
		angle = (double)m_dTheta * PI / 180.0;
	}
	else
	{
		angle = m_dTheta;
	}
	m_dX = m_dR * cos(angle);
	m_dY = m_dR * sin(angle);
	
}

//inline void CPoint::computeAngle()
//{
	//角度弧度转换，待补充
//}
void CPoint::getPointAtParam (CPoint &p1,CPoint &p2,CPoint &point,double &u){
	if(u <= 0 || u >= 1)
		cout << "参数位置不合理" << endl;
	else{
		point.setX (p1.getX () + u * (p2.getX () - p1.getX ()));
		point.setY (p1.getY () + u * (p2.getY () - p1.getY ()));
		point.setZ (p1.getZ () + u * (p2.getZ () - p1.getZ ()));
	}
}

void CPoint::getMidPoint (CPoint &pa,CPoint &pb,CPoint &mpoint){
	mpoint.setX ((pa.getX () + pb.getX ()) / 2.0);
	mpoint.setY ((pa.getY () + pb.getY ()) / 2.0);
	mpoint.setZ ((pa.getZ () + pb.getZ ()) / 2.0);
}

void CPoint::getSymmetryPoint (CPoint &p,CPoint &sp,CPoint &newpoint){
	newpoint.setX (2 * sp.getX () - p.getX ());
	newpoint.setY (2 * sp.getY () - p.getY ());
	newpoint.setZ (2 * sp.getZ () - p.getZ ());
}

double CPoint::getDistance(const CPoint &p) {
    double res = (m_dX - p.m_dX) * (m_dX - p.m_dX) + (m_dY - p.m_dY)* (m_dY - p.m_dY) + (m_dZ - p.m_dZ)* (m_dZ - p.m_dZ);
    return sqrt(res);
}
