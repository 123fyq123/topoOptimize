#include "SurfaceOfRevolution.h"
#include<GeomAdaptor_Surface.hxx>
#include<Extrema_ExtPS.hxx>
#include"DefineType.h"
SurfaceOfRevolution::SurfaceOfRevolution (){
}

SurfaceOfRevolution::~SurfaceOfRevolution (){
	if(!mSurface.IsNull ())
		mSurface.Nullify ();
}

SurfaceOfRevolution::SurfaceOfRevolution (CBSplineCurve &curve,CPoint &point,CVector &axl){
	if(curve.getGeomCurve ().IsNull ())
		cerr << "基准曲线参数错误";
	else if(axl.getX () == 0 && axl.getY () == 0 && axl.getZ () == 0){
		cerr << "拉伸方向参数错误";
	}
	else{
		mCurve = curve;
		mPoint = point;
		gp_Dir v (axl.getVec ());
		/*v.SetX (axl.getX ());
		v.SetY (axl.getY ());
		v.SetZ (axl.getZ ());*/
		gp_Ax1 ax1 (point.getPnt (),v);
		mSurface = new Geom_SurfaceOfRevolution (curve.getGeomCurve (),ax1);
	}
}

void SurfaceOfRevolution::setDirection (CPoint &point,CVector &axl){
	 if(axl.getX () == 0 && axl.getY () == 0 && axl.getZ () == 0){
		cerr << "拉伸方向参数错误";
	}
	 else{
		 gp_Dir v;
		 v.SetX (axl.getX ());
		 v.SetY (axl.getY ());
		 v.SetZ (axl.getZ ());
		 mSurface->SetDirection (v);
		 mSurface->SetLocation (point.getPnt ());
		 mPoint = point;
	 }
}

void SurfaceOfRevolution::setBasisBSplineCurve (CBSplineCurve &curve){
	if(curve.getGeomCurve ().IsNull ())
		cerr << "基准曲线参数错误";
	else{
		mSurface->SetBasisCurve (curve.getGeomCurve ());
		mCurve = curve;
	}
}

CPoint SurfaceOfRevolution::getLocation (){
	return mPoint;
}

CVector SurfaceOfRevolution::getAxis (){
	gp_Ax1 ax1;
	ax1=mSurface->Axis ();
	gp_Dir dir;
	dir=ax1.Direction ();
	CVector vec;
	vec.setX(dir.X ());
	vec.setY (dir.Y ());
	vec.setZ(dir.Z());
	return vec;
}

void SurfaceOfRevolution::getBounds (double &U1,double &U2,double &V1,double &V2){
	mSurface->Bounds (U1,U2,V1,V2);
}

bool SurfaceOfRevolution::isClosed (){
	if(mSurface->IsUClosed ())
		return true;
	else
		return false;
}

bool SurfaceOfRevolution::isCNc (int N){
	if(mSurface->IsCNu (N))
		return true;
	else
		return false;
}

void SurfaceOfRevolution::getPointAtParam (double &u,double &v, CPoint &point){
	double umin,umax,vmin,vmax;
	this->getBounds (umin,umax,vmin,vmax);
	if(u<umin || u>umax || v<vmin || v>vmax){
		cerr << "参数范围溢出";
	}
	else{
		gp_Pnt pnt;
		mSurface->D0 (u,v,pnt);
		point = CPoint (pnt);
	}
}

bool SurfaceOfRevolution::getParamAtPoint (const CPoint &point,double &u,double &v){
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

void SurfaceOfRevolution::getDerivation1 (CPoint &point,double &u,double &v,vector<CVector> &vectors){
	double umin,umax,vmin,vmax;
	this->getBounds (umin,umax,vmin,vmax);
	if(u<umin || u>umax || v<vmin || v>vmax){
		cerr << "参数范围溢出";
	}
	else{
		gp_Pnt pnt;
		gp_Vec vdu1,vdv1;
		mSurface->D1 (u,v,pnt,vdu1,vdv1);
		CVector u1(vdu1);
		CVector v1(vdv1);
		vectors.push_back(u1);
		vectors.push_back(v1);
	}
}

void SurfaceOfRevolution::getDerivation2 (CPoint &point,double &u,double &v,vector<CVector> &vectors){
	double umin,umax,vmin,vmax;
	this->getBounds (umin,umax,vmin,vmax);
	if(u<umin || u>umax || v<vmin || v>vmax){
		cerr << "参数范围溢出";
	}
	else{
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
}

void SurfaceOfRevolution::getDerivation3 (CPoint &point,double &u,double &v,vector<CVector> &vectors){
	double umin,umax,vmin,vmax;
	this->getBounds (umin,umax,vmin,vmax);
	if(u<umin || u>umax || v<vmin || v>vmax){
		cerr << "参数范围溢出";
	}
	else{
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
}
void SurfaceOfRevolution::setTranslation (CVector vec){
	gp_Vec v;
	gp_Trsf T;
	v.SetX (vec.getX ());
	v.SetY (vec.getY ());
	v.SetZ (vec.getZ ());
	T.SetTranslation (v);
	mSurface->Transform (T);
	mCurve.setTranslation (vec,mCurve);
	mPoint = CPoint (mSurface->Location ());
}

void SurfaceOfRevolution::setRotation (CPoint p,CVector ax,double ang){
	gp_Pnt pt;
	p.getPnt (pt);
	gp_Dir d (ax.getVec ());
	gp_Ax1 ax1 (pt,d);
	gp_Trsf T;
	T.SetRotation (ax1,ang);
	mSurface->Transform (T);
	mCurve.setRotation (p,ax,ang,mCurve);
	mPoint = CPoint (mSurface->Location ());
}
int SurfaceOfRevolution::getSurfaceType (){
	return REVOSURFACE;
}
Handle_Geom_SurfaceOfRevolution SurfaceOfRevolution::getGeomSurface (){
	return mSurface;
}
