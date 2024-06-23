#include "CCylindricalSurface.h"
#include "DefineType.h"
#include <GeomAdaptor_Surface.hxx>
#include <Extrema_ExtPS.hxx>
#include<gp_Ax1.hxx>
#include<gp_Ax3.hxx>
#include<gp_Dir.hxx>

CCylindricalSurface::CCylindricalSurface()
{
	m_dRadius = -1;
	m_Direction.setTypeCoord(-1);
	m_Direction.setXYZ(0, 0, 0);
	m_Location.setTypeCoord(-1);
}

CCylindricalSurface::CCylindricalSurface(const CCylindricalSurface & surface)
{
	m_Location = surface.m_Location;
	m_Direction = surface.m_Direction;
	m_dRadius = surface.m_dRadius;
	createCylindricalSurface();
	
}

CCylindricalSurface::CCylindricalSurface(const CPoint & point, const CVector & vec, const double & radius)
{
	m_Location = point;
	m_Direction = vec;
	m_dRadius = radius;
	createCylindricalSurface();
}

CCylindricalSurface::~CCylindricalSurface()
{
	if (!m_pSurface.IsNull())
	{
		m_pSurface.Nullify();
	}
}

int CCylindricalSurface::getSurfaceType()
{
	return CYLINDRICALSURFACE;
}

bool CCylindricalSurface::setRadius(double radius)
{
	if (radius > 0)
	{
		m_dRadius = radius;
		return true;
	}
	return false;
}

bool CCylindricalSurface::setLocation(CPoint & location)
{
	m_Location.setXYZ(location.getX(), location.getY(), location.getZ());
	m_Location.setTypeCoord(location.getTypeCoord());
	return true;
}

bool CCylindricalSurface::setDirection(CVector & direction)
{
	double len = sqrt(pow(direction.getX(), 2) + pow(direction.getY(), 2) + pow(direction.getZ(), 2));
	if (len <= 0)
	{
		return false;
	}
	m_Direction.setXYZ(direction.getX(), direction.getY(), direction.getZ());
	m_Direction.setTypeCoord(direction.getTypeCoord());
	return true;
}

bool CCylindricalSurface::buildSurface()
{
	bool result = createCylindricalSurface();
	return result;
}

bool CCylindricalSurface::getRadius(double & radius)
{
	radius = m_dRadius;
	return true;
}

double CCylindricalSurface::getRadius()
{
	return m_dRadius;
}

bool CCylindricalSurface::getLocation(CPoint & point)
{
	point.setXYZ(m_Location.getX(), m_Location.getY(), m_Location.getZ());
	return true;
}

bool CCylindricalSurface::getDirection(CVector & vec)
{
	vec.setXYZ(m_Direction.getX(), m_Direction.getY(), m_Direction.getZ());
	return true;
}

bool CCylindricalSurface::getDerivation1(double u, double v, vector<CVector>& vectors)
{
	if (m_pSurface.IsNull())
	{
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		return false;
	}

	gp_Vec vecU;
	gp_Vec vecV;
	gp_Pnt pnt;
	m_pSurface->D1(u, v, pnt, vecU, vecV);
	CVector vectorU(vecU);
	CVector vectorV(vecV);
	vectors.push_back(vectorU);
	vectors.push_back(vectorV);
	return true;
}

bool CCylindricalSurface::getDerivation2(double u, double v, vector<CVector>& vectors)
{
	if (m_pSurface.IsNull())
	{
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		return false;
	}
	gp_Vec vecU1, vecU2;
	gp_Vec vecV1, vecV2;
	gp_Vec vecUV;
	gp_Pnt pnt;
	m_pSurface->D2(u, v, pnt, vecU1, vecV1, vecU2, vecV2, vecUV);
	CVector vectorU2(vecU2);
	CVector vectorV2(vecV2);
	CVector vectorUV(vecUV);
	vectors.push_back(vectorU2);
	vectors.push_back(vectorV2);
	vectors.push_back(vectorUV);
	return true;
}

bool CCylindricalSurface::getDerivation3(double u, double v, vector<CVector>& vectors)
{
	if (m_pSurface.IsNull())
	{
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		vectors.push_back(tmpV);
		return false;
	}

	gp_Vec vecU1, vecU2, vecU3;
	gp_Vec vecV1, vecV2, vecV3;
	gp_Vec vecUV, vecUUV, vecUVV;
	gp_Pnt pnt;
	m_pSurface->D3(u, v, pnt, vecU1, vecV1, vecU2, vecV2, vecUV, vecU3, vecV3, vecUUV, vecUVV);
	CVector vectorU3(vecU3);
	CVector vectorV3(vecV3);
	CVector vectorUUV(vecUUV);
	CVector vectorUVV(vecUVV);
	vectors.push_back(vectorU3);
	vectors.push_back(vectorV3);
	vectors.push_back(vectorUUV);
	vectors.push_back(vectorUVV);
	return true;
}

bool CCylindricalSurface::getPointAtParam(double u, double v, CPoint & point)
{
	if (m_pSurface.IsNull())
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		point = tmp;
		return false;
	}
	double umin, umax, vmin, vmax;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	if ((u - umin) >= 0 && (umax - u) >= 0 && (v - vmin) >= 0 && (vmax - v) >= 0)
	{
		gp_Pnt pnt = m_pSurface->Value(u, v);
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

CPoint CCylindricalSurface::getPointAtParam(double u,double v){
	CPoint point;
	getPointAtParam(u,v,point);
	return point;
}

bool CCylindricalSurface::getParamAtPoint(const CPoint & point, double & u, double & v)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	gp_Pnt pnt;
	pnt = point.getPnt();

	double umin, umax, vmin, vmax;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	GeomAdaptor_Surface apdatorSurface;
	apdatorSurface.Load(m_pSurface);

	Extrema_ExtPS extPS;
	extPS.Initialize(apdatorSurface, umin, umax, vmin, vmax, 1.e-10, 1.e-10);
	extPS.Perform(pnt);

	bool isDone = extPS.IsDone() && (extPS.NbExt() > 0);
	int index;
	if (isDone)
	{
		double Dist2, Dist2Min = extPS.SquareDistance(1);
		index = 1;
		for (int i = 2; i <= extPS.NbExt(); i++)
		{
			Dist2 = extPS.SquareDistance(i);
			if (Dist2 < Dist2Min)
			{
				Dist2Min = Dist2;
				index = i;
			}
		}

		(extPS.Point(index)).Parameter(u, v);
		return true;
	}

	return false;
}

bool CCylindricalSurface::getUVMinMax(double & umin, double & umax, double & vmin, double & vmax)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	double u1, u2;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	CPoint pnt1(0, m_dRadius, m_Location.getZ(), 1, 1);
	getParamAtPoint(pnt1, u1, vmin);
	CPoint pnt2(0, m_dRadius, m_Location.getZ() + m_Direction.getZ(), 1, 1);
	getParamAtPoint(pnt2, u2, vmax);
	return true;
}

bool CCylindricalSurface::setTranslation(CGeomPoint & vec, CCylindricalSurface & surface)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	gp_Trsf translation;
	gp_Vec vecDir(vec.getX(), vec.getY(), vec.getZ());
	translation.SetTranslation(vecDir);
	CPoint p;
	CPoint p0;
	p0 = m_Location;
	gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
	translation.Transforms(coord);
	p.setXYZ(coord.X(), coord.Y(), coord.Z());

	surface.setLocation(p);
	surface.setDirection(m_Direction);
	surface.setRadius(m_dRadius);
	bool result = surface.buildSurface();
	return result;
}

bool CCylindricalSurface::setRotation(CGeomPoint & axisOrigin, CGeomPoint & axisDirection, double angle, CCylindricalSurface & surface)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}

	gp_Trsf rotation;
	gp_Ax1 axis;
	axis.SetDirection(gp_Dir(gp_Vec(axisDirection.getX(), axisDirection.getY(), axisDirection.getZ())));
	axis.SetLocation(gp_Pnt(axisOrigin.getX(), axisOrigin.getY(), axisOrigin.getZ()));
	rotation.SetRotation(axis, angle);
	CPoint p;
	CPoint p0;
	p0 = m_Location;
	gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
	rotation.Transforms(coord);
	p.setXYZ(coord.X(), coord.Y(), coord.Z());
	
	surface.setLocation(p);
	surface.setDirection(m_Direction);
	surface.setRadius(m_dRadius);
	bool result = surface.buildSurface();
	return result;
}

Handle_Geom_RectangularTrimmedSurface CCylindricalSurface::getGeomSurface() const
{
	return m_pSurface;

}

bool CCylindricalSurface::createCylindricalSurface()
{
	if (m_dRadius <= 0)
	{
		cerr << "CCylindricalSurface construct failed - radius nonvalid";
		return false;
	}
	double len = pow(m_Direction.getX(), 2) + pow(m_Direction.getY(), 2) + pow(m_Direction.getZ(), 2);
	if (len <= 0)
	{
		cerr << "CCylindricalSurface construct failed - direction lenth nonvalid";
		return false;
	}
	
	if (m_Direction.getTypeCoord() != -1 && m_Location.getTypeCoord() != -1)
	{
		gp_Dir dir(m_Direction.getVec());
		gp_Ax3 axis(m_Location.getPnt(), dir);
		try
		{
			if (!m_pSurface.IsNull())
			{
				m_pSurface.Nullify();
			}
			if (!m_pSurfaceBasis.IsNull())
			{
				m_pSurfaceBasis.Nullify();
			}
			
			m_pSurfaceBasis = new Geom_CylindricalSurface(axis, m_dRadius);
			double umin, umax, vmin, vmax;
			double u1, u2;
			m_pSurfaceBasis->Bounds(umin, umax, vmin, vmax);
			CVector normal;
			m_Direction.getPerpendicular(normal);
			normal.normalize();
			CPoint pnt1(m_Location.getX() + m_dRadius *normal.getX(), m_Location.getY() + m_dRadius * normal.getY(), m_Location.getZ() + m_dRadius *normal.getZ(), 1, 1);
			getParamAtPointBasis(pnt1, u1, vmin);
			CPoint pnt2(m_Location.getX() - m_dRadius * normal.getX() + m_Direction.getX(), m_Location.getY() - m_dRadius * normal.getY() + m_Direction.getY(), m_Location.getZ() - m_dRadius * normal.getZ() + m_Direction.getZ(), 1, 1);
			getParamAtPointBasis(pnt2, u2, vmax);
			m_pSurface = new Geom_RectangularTrimmedSurface(m_pSurfaceBasis, umin, umax, vmin, vmax);
			return true;
		}
		catch (const bad_alloc &)
		{
			cerr << "CylindricalSurface Contruct new fail" << endl;
			return false;
		}
	}
	else
	{
		return false;
	}

}

bool CCylindricalSurface::getParamAtPointBasis(const CPoint & point, double & u, double & v)
{
	if (m_pSurfaceBasis.IsNull())
	{
		return false;
	}
	gp_Pnt pnt;
	pnt = point.getPnt();

	double umin, umax, vmin, vmax;
	m_pSurfaceBasis->Bounds(umin, umax, vmin, vmax);
	GeomAdaptor_Surface apdatorSurface;
	apdatorSurface.Load(m_pSurfaceBasis);

	Extrema_ExtPS extPS;
	extPS.Initialize(apdatorSurface, umin, umax, vmin, vmax, 1.e-10, 1.e-10);
	extPS.Perform(pnt);

	bool isDone = extPS.IsDone() && (extPS.NbExt() > 0);
	int index;
	if (isDone)
	{
		double Dist2, Dist2Min = extPS.SquareDistance(1);
		index = 1;
		for (int i = 2; i <= extPS.NbExt(); i++)
		{
			Dist2 = extPS.SquareDistance(i);
			if (Dist2 < Dist2Min)
			{
				Dist2Min = Dist2;
				index = i;
			}
		}

		(extPS.Point(index)).Parameter(u, v);
		return true;
	}

	return false;
}
