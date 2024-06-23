#include "CCircle.h"
#include "DefineType.h"
#include <Eigen/Dense>
#include <GC_MakeArcOfCircle.hxx>
#include <GeomAdaptor_Curve.hxx>
#include <Extrema_ExtPC.hxx>
CCircle::CCircle()
{
}

CCircle::~CCircle()
{
	if (!m_pCircle.IsNull())
	{
		m_pCircle.Nullify();
	}

	if (!m_pArc.IsNull())
	{
		m_pArc.Nullify();
	}
}

CCircle::CCircle(CPoint & p1, CPoint & p2, CPoint & p3, int type)
{
	//´íÎó¼ì²â ´ý²¹³ä
	if (0 == type)
	{
		bool res = compute(p1, p2, p3);
		if (!res)
			return;
	}
	else if (1 == type)
	{
		double x1 = p1.getX();
		double y1 = p1.getY();
		double z1 = p1.getZ();

		double x2 = p2.getX();
		double y2 = p2.getY();
		double z2 = p2.getZ();

		double x3 = p3.getX();
		double y3 = p3.getY();
		double z3 = p3.getZ();
		if (check(p1, p2, p3))
		{
			m_cCenter = p1;

			double x = (y2 - y1) * (z3 - z1) - (y3 - y1)*(z2 - z1);
			double y = (z2 - z1)*(x3 - x1) - (x2 - x1) * (z3 - z1);
			double z = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
			m_cAxis.setX(x);
			m_cAxis.setY(y);
			m_cAxis.setZ(z);

			m_dRadius = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));
		}
		else
		{
			return ;
		}
	}
	gp_Dir vec(m_cAxis.getX(), m_cAxis.getY(), m_cAxis.getZ());
	gp_Pnt tmpP(m_cCenter.getX(), m_cCenter.getY(), m_cCenter.getZ());
	m_pCircle = new Geom_Circle(gp_Ax2(tmpP, vec), m_dRadius);
}

void CCircle::buildArc(CPoint & p1, CPoint & p2, CPoint & p3, int type)
{
	//´íÎó¼ì²â ´ý²¹³ä
	if (0 == type)
	{
		m_pArc = GC_MakeArcOfCircle(p1.getPnt(), p2.getPnt(), p3.getPnt());
		bool res = compute(p1, p2, p3);
		if (!res)
			return;
	}
	else if (1 == type)
	{
		double x1 = p1.getX();
		double y1 = p1.getY();
		double z1 = p1.getZ();

		double x2 = p2.getX();
		double y2 = p2.getY();
		double z2 = p2.getZ();

		double x3 = p3.getX();
		double y3 = p3.getY();
		double z3 = p3.getZ();
		if (check(p1, p2, p3))
		{
			m_cCenter = p1;

			double x = (y2 - y1) * (z3 - z1) - (y3 - y1)*(z2 - z1);
			double y = (z2 - z1)*(x3 - x1) - (x2 - x1) * (z3 - z1);
			double z = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
			m_cAxis.setX(x);
			m_cAxis.setY(y);
			m_cAxis.setZ(z);

			m_dRadius = sqrt(pow(x2 - x1, 2) + pow(y2 - y1, 2) + pow(z2 - z1, 2));

			gp_Dir vec(m_cAxis.getX(), m_cAxis.getY(), m_cAxis.getZ());
			gp_Pnt tmpP(m_cCenter.getX(), m_cCenter.getY(), m_cCenter.getZ());
			Handle_Geom_Circle tmpCircle = new Geom_Circle(gp_Ax2(tmpP, vec), m_dRadius);
			double u1 = INT_MIN;
			double u2 = INT_MIN;

			gp_Pnt pnt;
			pnt = p2.getPnt();

			GeomAdaptor_Curve apdatorCurve;
			apdatorCurve.Load(tmpCircle);

			Extrema_ExtPC extPC;
			extPC.Initialize(apdatorCurve, apdatorCurve.FirstParameter(), apdatorCurve.LastParameter());
			extPC.Perform(pnt);

			bool isDone = extPC.IsDone() && (extPC.NbExt() > 0);
			int index;
			if (isDone)
			{
				double Dist2, Dist2Min = extPC.SquareDistance(1);
				index = 1;
				for (int i = 2; i <= extPC.NbExt(); i++)
				{
					Dist2 = extPC.SquareDistance(i);
					if (Dist2 < Dist2Min)
					{
						Dist2Min = Dist2;
						index = i;
					}
				}

				u1 = (extPC.Point(index)).Parameter();
			}

			pnt = p3.getPnt();
			extPC.Perform(pnt);
			isDone = extPC.IsDone() && (extPC.NbExt() > 0);
			if (isDone)
			{
				double Dist2, Dist2Min = extPC.SquareDistance(1);
				index = 1;
				for (int i = 2; i <= extPC.NbExt(); i++)
				{
					Dist2 = extPC.SquareDistance(i);
					if (Dist2 < Dist2Min)
					{
						Dist2Min = Dist2;
						index = i;
					}
				}

				u2 = (extPC.Point(index)).Parameter();
			}

			if (u2 != INT_MIN && u1 != INT_MIN)
			{
				pnt = tmpCircle->Value(0.5*(u2 - u1));
				CPoint p4(pnt);
				m_pArc = GC_MakeArcOfCircle(p2.getPnt(), p4.getPnt(), p3.getPnt());
			}
			
			tmpCircle.Nullify();
		}
		else
		{
			return;
		}
		
	}
	
}

int CCircle::getCurveType()
{
	if (!m_pCircle.IsNull())
	{
		return CIRCLE;
	}
	if (!m_pArc.IsNull())
	{
		return ARC;
	}
	return CIRCLE;
}

void CCircle::getCenter(CPoint & center)
{
	if (!m_pCircle.IsNull() || !m_pArc.IsNull())
	{
		center.setXYZ(m_cCenter.getX(), m_cCenter.getY(), m_cCenter.getZ());
	}
	else
	{
		//error
		center.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
	}
}

CPoint CCircle::getCenter()
{
	if (!m_pCircle.IsNull() || !m_pArc.IsNull())
	{
		return m_cCenter;
	}
	else
	{
		//error
		return CPoint(INT_MAX, INT_MAX, INT_MAX, -1);
	}
	
}

void CCircle::getRadius(double & radius)
{
	if (!m_pCircle.IsNull() || !m_pArc.IsNull())
	{
		radius = m_dRadius;
	}
	else
	{
		//error
		radius = 0.0;
	}
}

double CCircle::getRadius()
{
	if (!m_pCircle.IsNull() || !m_pArc.IsNull())
	{
		return m_dRadius;
	}
	else
	{
		//error
		return 0.0;
	}
}

bool CCircle::getDerivation1(double u, CVector & vector)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Vec vec;
	gp_Pnt pnt;
	if (!m_pCircle.IsNull())
	{
		m_pCircle->D1(u, pnt, vec);
	}
	else
	{
		m_pArc->D1(u, pnt, vec);
	}
	
	vector.setXYZ(vec.X(), vec.Y(), vec.Z());
	return true;
}

bool CCircle::getDerivation2(double u, CVector & vector)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}

	int con = m_pCircle.IsNull() ? m_pArc->Continuity() : m_pCircle->Continuity();
	if (con < GeomAbs_C2)
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Pnt p;
	gp_Vec v1, v2;
	m_pCircle.IsNull() ? m_pArc->D2(u, p, v1, v2) : m_pCircle->D2(u, p, v1, v2);
	vector.setXYZ(v2.X(), v2.Y(), v2.Z());
	return true;
}

bool CCircle::getDerivation3(double u, CVector & vector)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	int con = m_pCircle.IsNull() ? m_pArc->Continuity() : m_pCircle->Continuity();
	if (con < GeomAbs_C3)
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Pnt p;
	gp_Vec v1, v2, v3;
	m_pCircle.IsNull() ? m_pArc->D3(u, p, v1, v2, v3) : m_pCircle->D3(u, p, v1, v2, v3);
	vector.setXYZ(v3.X(), v3.Y(), v3.Z());
	return true;
}

CPoint CCircle::getPointAtParam(double u)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}
	gp_Pnt pnt = m_pCircle.IsNull() ? m_pArc->Value(u) : m_pCircle->Value(u);
	CPoint point(pnt);
	return point;
}

bool CCircle::getPointAtParam(double u, CPoint & point)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		point = tmp;
		return false;
	}
	gp_Pnt pnt = m_pCircle.IsNull() ? m_pArc->Value(u) : m_pCircle->Value(u);
	point.setXYZ(pnt.X(), pnt.Y(), pnt.Z());
	return true;
}

bool CCircle::getParamAtPoint(const CPoint & point, double & u)
{
	if (m_pCircle.IsNull() && m_pArc.IsNull())
	{
		return false;
	}
	gp_Pnt pnt;
	pnt = point.getPnt();

	GeomAdaptor_Curve apdatorCurve;
	if (m_pArc.IsNull() && !m_pCircle.IsNull())
	{
		apdatorCurve.Load(m_pCircle);
	}
	else if (!m_pArc.IsNull() && m_pCircle.IsNull())
	{
		apdatorCurve.Load(m_pArc);
	}

	Extrema_ExtPC extPC;
	extPC.Initialize(apdatorCurve, apdatorCurve.FirstParameter(), apdatorCurve.LastParameter());
	extPC.Perform(pnt);

	bool isDone = extPC.IsDone() && (extPC.NbExt() > 0);
	int index;
	if (isDone)
	{
		double Dist2, Dist2Min = extPC.SquareDistance(1);
		index = 1;
		for (int i = 2; i <= extPC.NbExt(); i++)
		{
			Dist2 = extPC.SquareDistance(i);
			if (Dist2 < Dist2Min)
			{
				Dist2Min = Dist2;
				index = i;
			}
		}

		u = (extPC.Point(index)).Parameter();
		return true;
	}

	return false;
}

double CCircle::getUmin()
{
	double umin = -1;
	if (!m_pCircle.IsNull())
	{
		umin = m_pCircle->FirstParameter();
	}
	else if (!m_pArc.IsNull())
	{
		umin = m_pArc->FirstParameter();
	}
	return umin;
}

double CCircle::getUmax()
{
	double umax = -1;
	if (!m_pCircle.IsNull())
	{
		umax = m_pCircle->LastParameter();
	}
	else if (!m_pArc.IsNull())
	{
		umax = m_pArc->LastParameter();
	}
	return umax;
}

double CCircle::getDeltaU()
{
	double umin = 0, umax = 0;
	if (!m_pCircle.IsNull())
	{
		umin = m_pCircle->FirstParameter();
		umax = m_pCircle->LastParameter();
	}
	else if (!m_pArc.IsNull())
	{
		umin = m_pArc->FirstParameter();
		umax = m_pArc->LastParameter();
	}

	return umax - umin;
}

CPoint CCircle::getStartPoint()
{
	if (!m_pArc.IsNull())
	{
		gp_Pnt pnt = m_pArc->StartPoint();
		CPoint point(pnt);
		return point;
	}
	else
	{
		CPoint point(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return point;
	}
}

bool CCircle::getStartPoint(CPoint & point)
{
	if (!m_pArc.IsNull())
	{
		gp_Pnt pnt = m_pArc->StartPoint();
		point.setXYZ(pnt.X(), pnt.Y(), pnt.Z());
		return true;
	}
	else
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		point = tmp;
		return false;
	}
}

CPoint CCircle::getEndPoint()
{
	if (!m_pArc.IsNull())
	{
		gp_Pnt pnt = m_pArc->EndPoint();
		CPoint point(pnt);
		return point;
	}
	else
	{
		CPoint point(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return point;
	}
}

bool CCircle::getEndPoint(CPoint & point)
{
	if (!m_pArc.IsNull())
	{
		gp_Pnt pnt = m_pArc->EndPoint();
		point.setXYZ(pnt.X(), pnt.Y(), pnt.Z());
		return true;
	}
	else
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		point = tmp;
		return false;
	}
}

Handle_Geom_Curve CCircle::getGeomCurve()
{
	if (!m_pCircle.IsNull())
		return m_pCircle;
	if (!m_pArc.IsNull())
		return m_pArc;
	return m_pCircle;
}
//
//void CCircle::setGeomCurve(Handle_Geom_Curve & c)
//{
//	if (!m_pCircle.IsNull())
//	{
//		c = m_pCircle;
//		return;
//	}
//		
//	if (!m_pArc.IsNull())
//	{
//		c = m_pArc;
//		return;
//	}
//		
//}

bool CCircle::compute(CPoint & p1, CPoint & p2, CPoint & p3)
{
	Eigen::MatrixXd matrix(3, 3);
	Eigen::VectorXd vecB(3);
	Eigen::VectorXd vecX(3);
	double x1 = p1.getX();
	double y1 = p1.getY();
	double z1 = p1.getZ();

	double x2 = p2.getX();
	double y2 = p2.getY();
	double z2 = p2.getZ();
	
	double x3 = p3.getX();
	double y3 = p3.getY();
	double z3 = p3.getZ();
	if (check(p1, p2, p3))
	{
		matrix(0, 0) = y1 * z2 - y1 * z3 - z1 * y2 + z1 * y3 + y2 * z3 - y3 * z2;
		matrix(0, 1) = -x1 * z2 + x1 * z3 + z1 * x2 - z1 * x3 - x2 * z3 + x3 * z2;
		matrix(0, 2) = x1 * y2 - x1 * y3 - y1 * x2 + y1 * x3 + x2 * y3 - x3 * y2;
		vecB(0) = -x1 * y2*z3 + x1 * y3*z2 + x2 * y1*z3 - x3 * y1*z2 - x2 * y3*z1 + x3 * y2*z1;

		matrix(1, 0) = 2 * (x2 - x1);
		matrix(1, 1) = 2 * (y2 - y1);
		matrix(1, 2) = 2 * (z2 - z1);
		vecB(1) = pow(x1, 2) + pow(y1, 2) + pow(z1, 2) - pow(x2, 2) - pow(y2, 2) - pow(z2, 2);

		matrix(2, 0) = 2 * (x3 - x1);
		matrix(2, 1) = 2 * (y3 - y1);
		matrix(2, 2) = 2 * (z3 - z1);
		vecB(2) = pow(x1, 2) + pow(y1, 2) + pow(z1, 2) - pow(x3, 2) - pow(y3, 2) - pow(z3, 2);

		vecX = matrix.colPivHouseholderQr().solve(vecB);

		double r1 = sqrt(pow(x1 + vecX(0), 2) + pow(y1 + vecX(1), 2) + pow(z1 + vecX(2), 2));
		double r2 = sqrt(pow(x2 + vecX(0), 2) + pow(y2 + vecX(1), 2) + pow(z2 + vecX(2), 2));
		double r3 = sqrt(pow(x3 + vecX(0), 2) + pow(y3 + vecX(1), 2) + pow(z3 + vecX(2), 2));
		if (abs(r1 - r2) < 1e-15 && abs(r3 - r1) < 1e-15)
		{
			m_cCenter.setX(-vecX(0));
			m_cCenter.setY(-vecX(1));
			m_cCenter.setZ(-vecX(2));
			m_dRadius = r1;
			double x = (y2 - y1) * (z3 - z1) - (y3 - y1)*(z2 - z1);
			double y = (z2 - z1)*(x3 - x1) - (x2 - x1) * (z3 - z1);
			double z = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);
			m_cAxis.setX(x);
			m_cAxis.setY(y);
			m_cAxis.setZ(z);

			return true;
		}
		else
		{
			return false;
		}
	}
	else
		return false;
}

bool CCircle::check(CPoint & p1, CPoint & p2, CPoint & p3)
{
	double x1 = p1.getX();
	double y1 = p1.getY();
	double z1 = p1.getZ();

	double x2 = p2.getX();
	double y2 = p2.getY();
	double z2 = p2.getZ();

	double x3 = p3.getX();
	double y3 = p3.getY();
	double z3 = p3.getZ();

	double x = (y2 - y1) * (z3 - z1) - (y3 - y1)*(z2 - z1);
	double y = (z2 - z1)*(x3 - x1) - (x2 - x1) * (z3 - z1);
	double z = (x2 - x1)*(y3 - y1) - (x3 - x1)*(y2 - y1);

	if ((x < 0 || x > 0) && (y < 0 || y > 0) && (z < 0 || z > 0))
	{
		return true;
	}
	return false;
}
