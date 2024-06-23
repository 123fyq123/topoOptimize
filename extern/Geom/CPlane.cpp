#include "CPlane.h"
#include "CAidFunction.h"
#include "DefineType.h"
#include<GeomAdaptor_Surface.hxx>
#include<Extrema_ExtPS.hxx>
CPlane::CPlane()
{
}

CPlane::CPlane( CPoint & point,  CVector & vec)
{
	double x = vec.getX();
	double y = vec.getY();
	double z = vec.getZ();
	if ((x < 0 || x > 0) && (y < 0 || y > 0) && (z < 0 || z > 0))
	{
		m_cPoint = point;
		m_cNormal = vec;
		//计算任意平面上的三个点,待补充
		gp_Dir v(vec.getX(), vec.getY(), vec.getZ());
		/*v.SetX (Standard_Real(vec.getX ()));
		v.SetY (Standard_Real(vec.getY ()));
		v.SetZ (Standard_Real(vec.getZ ()));*/
		this->Plane = new Geom_Plane(point.getPnt(), v);
		double umax, umin, vmax, vmin;
		Plane->Bounds(umin, umax, vmin, vmax);
		CVector oppvec(m_cNormal.getX() - m_cPoint.getX(), m_cNormal.getY() - m_cPoint.getY(), m_cNormal.getZ() - m_cPoint.getZ());
		gp_Pnt pntmin, pntmax;
		Plane->D0(-100, -100, pntmin);
		Plane->D0(100, 100, pntmax);
		CPoint pmin(pntmin);
		getParamAtPointBasis(pmin, umin, vmin);
		CPoint pmax(pntmax);
		getParamAtPointBasis(pmax, umax, vmax);
		mRPlane = new Geom_RectangularTrimmedSurface(Plane, umin, umax, vmin, vmax);
	}
}

CPlane::CPlane( CPoint & point1,  CPoint & point2,  CPoint & point3)
{
	CVector vec;
	CPoint p1(point2.getX() - point1.getX(), point2.getY() - point1.getY(), point2.getY() - point1.getY());
	CPoint p2(point3.getX() - point1.getX(), point3.getY() - point1.getY(), point3.getY() - point1.getY());
	CAidFunction::crossProduct3D(p1, p2, vec);
	double x = vec.getX();
	double y = vec.getY();
	double z = vec.getZ();
	if ((x < 0 || x > 0) && (y < 0 || y > 0) && (z < 0 || z > 0))
	{
		m_cPoint1 = point1;
		m_cPoint2 = point2;
		m_cPoint3 = point3;
		m_cPoint = point1;
		m_cNormal.setXYZ(x, y, z);
		gp_Dir v(x, y, z);
		Plane = new Geom_Plane(point1.getPnt(), v);
		double umax, umin, vmax, vmin;
		Plane->Bounds(umin, umax, vmin, vmax);
		CVector oppvec(m_cNormal.getX() - m_cPoint.getX(), m_cNormal.getY() - m_cPoint.getY(), m_cNormal.getZ() - m_cPoint.getZ());
		gp_Pnt pntmin, pntmax;
		Plane->D0(-100, -100, pntmin);
		Plane->D0(100, 100, pntmax);
		CPoint pmin(pntmin);
		getParamAtPointBasis(pmin, umin, vmin);
		CPoint pmax(pntmax);
		getParamAtPointBasis(pmax, umax, vmax);
		mRPlane = new Geom_RectangularTrimmedSurface(Plane, umin, umax, vmin, vmax);
	}
	
}


CPlane::~CPlane()
{
}

void CPlane::getNormal(CGeomPoint &normal)
{
	//normal = m_cNormal;
	normal.setXYZ(m_cNormal.getX(), m_cNormal.getY(), m_cNormal.getZ());
}

void CPlane::getPoint(CGeomPoint &point)
{
	//point = m_cPoint;
	point.setXYZ(m_cPoint.getX(), m_cPoint.getY(), m_cPoint.getZ());
}

void CPlane::getBounds (double &u1,double &u2,double &v1,double &v2){
	Plane->Bounds (u1,u2,v1,v2);
}

void CPlane::getPointAtParam (CPoint &point,double u,double v){
	gp_Pnt p;
	Plane->D0 (u,v,p);
	point = CPoint (p);
}

void CPlane::getBoundPoint (vector<CPoint> &boundpoint){
	if (Plane.IsNull())
		return;
	double U1,U2,V1,V2;
	gp_Pnt P1,P2,P3,P4;
	Plane->Bounds (U1,U2,V1,V2);
	Plane->D0 (U1,V1,P1);
	Plane->D0 (U1,V2,P2);
	Plane->D0 (U2,V1,P3);
	Plane->D0 (U2,V2,P4);
	boundpoint.push_back (CPoint (P1));
	boundpoint.push_back (CPoint (P2));
	boundpoint.push_back (CPoint (P3));
	boundpoint.push_back (CPoint (P4));

}
bool CPlane::getParamAtPointBasis (const CPoint &point,double &u,double &v){
	if(Plane.IsNull ()){
		return false;
	}
	gp_Pnt pnt;
	pnt = point.getPnt ();

	double umin,umax,vmin,vmax;
	Plane->Bounds (umin,umax,vmin,vmax);
	GeomAdaptor_Surface apdatorSurface;
	apdatorSurface.Load (Plane);
	Extrema_ExtPS extPS;
	extPS.Initialize (apdatorSurface,umin,umax,vmin,vmax,1.e-10,1.e-10);
	extPS.Perform (pnt);

	bool isDone = extPS.IsDone () && (extPS.NbExt () > 0);
	int index;
	if(isDone){
		double Dist2,Dist2Min = extPS.SquareDistance (1);
		index = 1;
		for(int i = 2; i <= extPS.NbExt (); i++){
			Dist2 = extPS.SquareDistance (i);
			if(Dist2 < Dist2Min){
				Dist2Min = Dist2;
				index = i;
			}
		}

		(extPS.Point (index)).Parameter (u,v);
		return true;
	}

	return false;
}
Handle_Geom_Plane CPlane::getGeomPlane (){
	return Plane;
}

Handle_Geom_RectangularTrimmedSurface CPlane::getRGeomPlane (){
	return mRPlane;
}

int CPlane::getSurfaceType()
{
	return PLANE;
}
