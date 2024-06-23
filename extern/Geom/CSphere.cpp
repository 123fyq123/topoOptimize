#include "CSphere.h"
#include<gp_Ax3.hxx>
#include<gp_Dir.hxx>
#include<gp_Vec.hxx>
#include<GeomAdaptor_Surface.hxx>
#include<Extrema_ExtPS.hxx>
CSphere::CSphere (){
}

CSphere::~CSphere (){
	if(!mSurface.IsNull())
	   mSurface.Nullify ();
}

CSphere::CSphere (CPoint point,double radius){
	gp_Pnt pnt;
	point.getPnt (pnt);
	CVector vec (point.getX (),point.getY (),point.getZ ());
	gp_Dir dir(vec.getVec ());
	gp_Ax3 ax3 (pnt,dir);
	if(radius <= 0)
		cerr << "°ë¾¶²ÎÊý´íÎó";
	else{
		mSurface = new Geom_SphericalSurface (ax3,radius);
		upCSpheredata ();
	}

}

void CSphere::setLocation (CPoint point){
	mSurface->SetLocation (point.getPnt ());
	mPoint = CPoint(mSurface->Location ());
}

void CSphere::setRadius (double radius){
	if(radius <= 0)
		cerr << "°ë¾¶²ÎÊý´íÎó";
	else{
		mSurface->SetRadius (radius);
		mRadius = mSurface->Radius ();
	}
}

double CSphere::computeArea (){
	return mSurface->Area();
}

double CSphere::computeVolume (){
	return mSurface->Volume();
}

void CSphere::getUVmin (double &umin,double &vmin){
	double u,v;
	mSurface->Bounds (umin,u,vmin,v);
}

void CSphere::getUVmax (double &umax,double &vmax){
	double u,v;
	mSurface->Bounds (u,umax,v,vmax);
}

void CSphere::computeCoefficients (double &A1,double &A2,double &A3,
	double &B1,double &B2,double &B3,
	double &C1,double &C2,double &C3,
	double &D){
	mSurface->Coefficients (A1,A2,A3,B1,B2,B3,C1,C2,C3,D);
}

bool CSphere::IsUClosed (){
	return mSurface->IsUClosed();
}

bool CSphere::IsVClosed (){
	return mSurface->IsVClosed();
}

bool CSphere::IsUPeriodic (){
	return mSurface->IsUPeriodic();
}

bool CSphere::IsVPeriodic (){
	return mSurface->IsVPeriodic();
}

void CSphere::getPointAtParam (double u,double v, CPoint &point){
	gp_Pnt pnt;
	mSurface->D0 (u,v,pnt);
	point = CPoint (pnt);
}

bool CSphere::getParamAtPoint (const CPoint &point,double &u,double &v){
	if(mSurface.IsNull ()){
		return false;
	}
	gp_Pnt pnt;
	pnt = point.getPnt ();

	double umin,umax,vmin,vmax;
	mSurface->Bounds (umin,umax,vmin,vmax);
	GeomAdaptor_Surface apdatorSurface;
	apdatorSurface.Load (mSurface);

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

void CSphere::getDerivation1 (double &u, double &v, vector<CVector> &vectors){
	gp_Pnt pnt;
	gp_Vec vdu1,vdv1;
	mSurface->D1 (u,v,pnt,vdu1,vdv1);
	CVector u1(vdu1);
	CVector v1(vdv1);
	vectors.push_back(u1);
	vectors.push_back(v1);
}

void CSphere::getDerivation2 (double &u,double &v,vector<CVector> &vectors){
	gp_Pnt pnt;
	gp_Vec vdu1,vdv1,vdu2,vdv2,vduv2;
	mSurface->D2 (u,v,pnt,vdu1,vdv1,vdu2,vdv2,vduv2);
	CVector u2(vdu2);
	CVector v2(vdv2);
	CVector uv2(vduv2);
	vectors.push_back(u2);
	vectors.push_back(v2);
	vectors.push_back(uv2);
}

void CSphere::getDerivation3 (double &u,double &v,vector<CVector> &vectors){
	gp_Pnt pnt;
	gp_Vec vdu1,vdv1,vdu2,vdv2,vduv2,vdu3,vdv3,vduuv3,vduvv3;
	mSurface->D3 (u,v,pnt,vdu1,vdv1,vdu2,vdv2,vduv2,vdu3,vdv3,vduuv3,vduvv3);
	CVector u3(vdu3);
	CVector v3(vdv3);
	CVector uuv3(vduuv3);
	CVector uvv3(vduvv3);
	vectors.push_back(u3);
	vectors.push_back(v3);
	vectors.push_back(uuv3);
	vectors.push_back(uvv3);
}
bool CSphere::setTranslation (CVector &vec,CSphere &sphere){
	if(mSurface.IsNull ()){
		return false;
	}
	gp_Trsf translation;
	gp_Vec vecDir (vec.getX (),vec.getY (),vec.getZ ());
	translation.SetTranslation (vecDir);
	mSurface->Transform (translation);
	mPoint =CPoint(mSurface->Location ());
	return true;
}

int CSphere::getSurfaceType (){
	return SPHERE;
}

Handle_Geom_SphericalSurface CSphere::getGeomSurface (){
	return mSurface;
}

void CSphere::upCSpheredata (){
	mRadius=mSurface->Radius ();
	mPoint = mSurface->Location ();
}
