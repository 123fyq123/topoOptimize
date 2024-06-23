#include "CLine.h"
#include "CAidFunction.h"
#include "DefineType.h"

CLine::CLine()
{
	m_cVector.setXYZ(0, 0, 0);
	m_cPoint.setXYZ(0, 0, 0);
	m_dStart = -1;
	m_dEnd = -1;
}

CLine::CLine(CGeomPoint & point, CGeomPoint & vec)
{
	m_cPoint = dynamic_cast<CPoint &> (point);
	m_cVector = dynamic_cast<CVector &> (vec);
	m_pCurve = new Geom_Line(m_cPoint.getPnt(), m_cVector.getVec());
	m_dStart = 0;
	m_dEnd = sqrt(m_cVector.getX() * m_cVector.getX() + m_cVector.getY() * m_cVector.getY() + m_cVector.getZ() * m_cVector.getZ());
}

CLine::CLine(CGeomSurface &plane1, CGeomSurface &plane2)
{
	CVector planeNormal1, planeNormal2;
	CPoint p, normal;
	CVector direction, v;
	CVector direc;
	CPoint pnt;

	plane1.getNormal(planeNormal1);
	plane2.getNormal(planeNormal2);
	
	CAidFunction::crossProduct3D(planeNormal1, planeNormal2, direction);
	direction.normalize();

	direc.setXYZ(direction.getX(), direction.getY(), direction.getZ());
	m_cVector = direc;

	plane1.getPoint(p);
	plane1.getNormal(normal);
	normal.getPerpendicular(v);
	CLine line2(p, v);
	double t = CAidFunction::intersectLinePlane(line2, plane2);
	CPoint tmpPnt = line2.getPointAtParam(t);
	pnt.setXYZ(tmpPnt.getX(), tmpPnt.getY(), tmpPnt.getZ());
	m_cPoint = pnt;

	m_pCurve = new Geom_Line(m_cPoint.getPnt(), m_cVector.getVec());
	m_dStart = 0;
	m_dEnd = sqrt(m_cVector.getX() * m_cVector.getX() + m_cVector.getY() * m_cVector.getY() + m_cVector.getZ() * m_cVector.getZ());
}


CLine::~CLine()
{
	if (!m_pCurve.IsNull())
	{
		m_pCurve.Nullify();
	}
	/*if (!m_pTrimmedCurve.IsNull())
	{
		m_pTrimmedCurve.Nullify();
	}*/
}

int CLine::getCurveType() const
{
	return LINE;
}

void CLine::setStartPoint(CPoint & point)
{
	m_cPoint = point;
	m_dStart = 0;

	if (!m_pCurve.IsNull())
	{
		m_pCurve->SetLocation(m_cPoint.getPnt());
	}
	else if (m_dStart != -1 && m_dEnd != -1)
	{
		m_pCurve = new Geom_Line(m_cPoint.getPnt(), m_cVector.getVec());
	}
}

void CLine::setDirection(CVector & vec)
{
	m_cVector = vec;
	m_dEnd = sqrt(m_cVector.getX() * m_cVector.getX() + m_cVector.getY() * m_cVector.getY() + m_cVector.getZ() * m_cVector.getZ());
	if (!m_pCurve.IsNull())
	{
		m_pCurve->SetDirection(m_cVector.getVec());
	}
	else if (m_dStart != -1 && m_dEnd != -1)
	{
		m_pCurve = new Geom_Line(m_cPoint.getPnt(), m_cVector.getVec());
	}
}

CPoint CLine::getPointAtParam(double u)
{
	CVector vec(m_cVector);
	vec.normalize();
	return (m_cPoint) + ((vec) * u);
}

CPoint CLine::getStartPoint(void)
{
	return m_cPoint;
}

CPoint CLine::getEndPoint(void)
{
	return m_cPoint + m_cVector;
}

double CLine::getUmin(void)
{
	return m_dStart;
}

double CLine::getUmax(void)
{
	return m_dEnd;
}

void CLine::getDirection(CGeomPoint &vec)
{
	vec.setXYZ(m_cVector.getX(), m_cVector.getY(), m_cVector.getZ());
}

int CLine::getCurveType()
{
	return LINE;
}

void CLine::rotate(CGeom & line, double angle, int typeAngle, CLine & result)
{
	CPoint point(m_cPoint.getX() + m_cVector.getX(), m_cPoint.getY() + m_cVector.getY(), m_cPoint.getZ() + m_cVector.getZ(), 1, 0);
	CPoint point1;//记录t=0的变换
	CPoint point2;//记录vector的变换

	m_cPoint.rotate(line, angle, typeAngle, point1);
	point.rotate(line, angle, typeAngle, point2);

	CVector vec(point2.getX() - point1.getX(), point2.getY() - point1.getY(), point2.getX() - point1.getZ());
	result.setStartPoint(point1);
	result.setDirection(vec);
}

void CLine::mirror(CGeom & plane, CLine & result)
{
	CPoint point(m_cPoint.getX() + m_cVector.getX(), m_cPoint.getY() + m_cVector.getY(), m_cPoint.getZ() + m_cVector.getZ(), 1, 0);
	CPoint point1;//记录t=0的变换
	CPoint point2;//记录vector的变换

	m_cPoint.mirror(plane, point1);
	point.mirror(plane, point2);

	CVector vec(point2.getX() - point1.getX(), point2.getY() - point1.getY(), point2.getX() - point1.getZ());
	result.setStartPoint(point1);
	result.setDirection(vec);
}

void CLine::translate(CGeomPoint & vec, CLine & result)
{
	CPoint point(m_cPoint.getX() + m_cVector.getX(), m_cPoint.getY() + m_cVector.getY(), m_cPoint.getZ() + m_cVector.getZ(), 1, 0);
	CPoint point1 = m_cPoint + dynamic_cast<CVector &> (vec);//记录t=0的变换
	CPoint point2 = point + dynamic_cast<CVector &> (vec);//记录vector的变换

	CVector vect(point2.getX() - point1.getX(), point2.getY() - point1.getY(), point2.getX() - point1.getZ());
	result.setStartPoint(point1);
	result.setDirection(vect);
}

double CLine::computeLength(double u1, double u2)
{
	CVector vec0 = m_cVector;
	vec0.normalize();
	CVector vec = (vec0)* (u2 - u1);
	double len = 0;
	len = sqrt(vec.getX() * vec.getX() + vec.getY() * vec.getY() + vec.getZ() * vec.getZ());
	return len;
}

void CLine::setTrimmed(double u1, double u2, CLine & line)
{
	if (m_pCurve.IsNull())
		return;
	gp_Pnt pnt1 = m_pCurve->Value(u1);
	gp_Pnt pnt2 = m_pCurve->Value(u2);
	CPoint start(pnt1);
	CVector vec(pnt2.X() - pnt1.X(), pnt2.Y() - pnt1.Y(), pnt2.Z() - pnt1.Z());
	line.setStartPoint(start);
	line.setDirection(vec);
}

Handle_Geom_Curve CLine::getGeomCurve()
{
	if (!m_pCurve.IsNull())
		return m_pCurve;
	else
		return Handle_Geom_Line();
}

//Handle_Geom_TrimmedCurve CLine::getTrimmedCurve()
//{
//	if(m_pCurve.IsNull())
//		return Handle_Geom_TrimmedCurve();
//
//	if (!m_pTrimmedCurve.IsNull())
//		m_pTrimmedCurve.Nullify();
//
//	m_pTrimmedCurve = new Geom_TrimmedCurve(m_pCurve, m_dStart, m_dEnd);
//	return m_pTrimmedCurve;
//}

void CLine::setGeomCurve(Handle_Geom_Line & c)
{
	// 单位线矢量
	m_pCurve = c;
	CPoint p1(m_pCurve->Value(0));
	CPoint p2(m_pCurve->Value(1));
	m_cPoint = p1;

	m_cVector = CVector(p2.getX() - p1.getX(), p2.getY() - p1.getY(), p2.getZ() - p1.getZ());
	m_dStart = 0;
	m_dEnd = sqrt(m_cVector.getX() * m_cVector.getX() + m_cVector.getY() * m_cVector.getY() + m_cVector.getZ() * m_cVector.getZ());
}
