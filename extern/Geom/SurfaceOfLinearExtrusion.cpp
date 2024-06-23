
#include "SurfaceOfLinearExtrusion.h"
#include<gp_Ax1.hxx>
#include<GeomAdaptor_Curve.hxx>
#include<GeomAdaptor_Surface.hxx>
#include<Extrema_ExtPS.hxx>
#include"DefineType.h"
SurfaceOfLinearExtrusion::SurfaceOfLinearExtrusion (){

}
SurfaceOfLinearExtrusion::~SurfaceOfLinearExtrusion (){
	if(!mSurface.IsNull())
		mSurface.Nullify();
}

SurfaceOfLinearExtrusion::SurfaceOfLinearExtrusion (CBSplineCurve &curve,  CVector &axl){
	if(curve.getGeomCurve ().IsNull ())
		cerr << "基准曲线参数错误";
	else if(axl.getX () == 0 && axl.getY () == 0 && axl.getZ () == 0){
		cerr << "拉伸方向参数错误";
	}
	else{
		mcurve = curve;
		maxl = axl;
		gp_Dir v (axl.getVec ());
		/*v.SetX (axl.getX ());
		v.SetY (axl.getY ());
		v.SetZ (axl.getZ ());*/
		mSurface = new Geom_SurfaceOfLinearExtrusion (curve.getGeomCurve (),v);
		double u1,u2,v1,v2;
		double umax, umin;//,vmax,vmin;
		double up,vp;
		mSurface->Bounds (u1,u2,v1,v2);
		umax = mcurve.getUmax ();
		umin = mcurve.getUmin ();
		CPoint pnt (mcurve.getPointAtParam (umax).getX () + maxl.getX (),
			mcurve.getPointAtParam (umax).getY () + maxl.getY (),
			mcurve.getPointAtParam (umax).getZ () + maxl.getZ (),
			1,1);
		getParamAtPointBasis (mcurve.getPointAtParam (umax),up,v1);
		getParamAtPointBasis (pnt,vp,v2);
		mRSurface = new Geom_RectangularTrimmedSurface (mSurface,u1,u2,v1,v2);
	}
}

bool SurfaceOfLinearExtrusion::getParamAtPointBasis (const CPoint &point,double &u,double &v){
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

void SurfaceOfLinearExtrusion::setDirection ( CVector &axl){
	if(axl.getX () == 0 && axl.getY () == 0 && axl.getZ () == 0){
		cerr << "拉伸方向参数错误";
	}
	else{
		maxl = axl;
		gp_Dir v;
		v.SetX (axl.getX ());
		v.SetY (axl.getY ());
		v.SetZ (axl.getZ ());
		mSurface->SetDirection (v);
	}
}

void SurfaceOfLinearExtrusion::setBasisBSplineCurve (CBSplineCurve &curve){
	if(curve.getGeomCurve ().IsNull ()){
		cerr << "曲线参数错误";
	}
	else{
		mcurve = curve;
		mSurface->SetBasisCurve (curve.getGeomCurve ());
	}
}

void SurfaceOfLinearExtrusion::getBounds (double &U1,double &U2,double &V1,double &V2){
	if(mSurface.IsNull ()){
		cerr << "曲面参数错误";
	}
	else{
		mSurface->Bounds (U1,U2,V1,V2);
	}
}

void SurfaceOfLinearExtrusion::uReverse (){
	if(mSurface.IsNull ()){
		cerr << "曲面参数错误";
	}
	else{
		mSurface->UReverse ();
	}
}

void SurfaceOfLinearExtrusion::vReverse (){
	if(mSurface.IsNull ()){
		cerr << "曲面参数错误";
	}
	else{
		mSurface->VReverse ();
	}
}

bool SurfaceOfLinearExtrusion::isClosed (){
	if(mSurface->IsUClosed ())
		return true;
	else
		return false;
}

bool SurfaceOfLinearExtrusion::isCNc (int N){
	if(mSurface->IsCNu (N))
		return true;
	else
		return false;
}

void SurfaceOfLinearExtrusion::getPointAtParam (double &u,double &v, CPoint &point){
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

bool SurfaceOfLinearExtrusion::getParamAtPoint (const CPoint &point,double &u,double &v){
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

void SurfaceOfLinearExtrusion::getDerivation1 (double &u, double &v, vector<CVector> &vectors){
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

void SurfaceOfLinearExtrusion::getDerivation2 (double &u,double &v, vector<CVector> &vectors){
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

void SurfaceOfLinearExtrusion::getDerivation3 (double &u,double &v, vector<CVector> &vectors){
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

void SurfaceOfLinearExtrusion::setTranslation (CVector vec){
	gp_Vec v;
	gp_Trsf T;
	v.SetX (vec.getX ());
	v.SetY (vec.getY ());
	v.SetZ (vec.getZ ());
	T.SetTranslation (v);
	mSurface->Transform (T);
	double u1,u2,v1,v2;
	double umax, umin;// , vmax, vmin;
	double up,vp;
	mSurface->Bounds (u1,u2,v1,v2);
	umax = mcurve.getUmax ();
	umin = mcurve.getUmin ();
	CPoint pnt (mcurve.getPointAtParam (umax).getX () + maxl.getX (),
		mcurve.getPointAtParam (umax).getY () + maxl.getY (),
		mcurve.getPointAtParam (umax).getZ () + maxl.getZ (),
		1,1);
	getParamAtPointBasis (mcurve.getPointAtParam (umax),up,v1);
	getParamAtPointBasis (pnt,vp,v2);
	mRSurface = new Geom_RectangularTrimmedSurface (mSurface,u1,u2,v1,v2);
	mcurve.setTranslation (vec,mcurve);
	gp_Dir dir;
	dir=mSurface->Direction ();
	maxl.setX (dir.X ());
	maxl.setY(dir.Y ());
	maxl.setZ (dir.Z ());
}

void SurfaceOfLinearExtrusion::setRotation (CPoint p,CVector ax,double ang){
	gp_Pnt pt;
	p.getPnt (pt);
	gp_Dir d (ax.getVec ());
	gp_Ax1 ax1 (pt,d);
	gp_Trsf T;
	T.SetRotation (ax1,ang);
	mRSurface->Transform (T);
	double u1,u2,v1,v2;
	double umax, umin;// , vmax, vmin;
	double up,vp;
	mSurface->Bounds (u1,u2,v1,v2);
	umax = mcurve.getUmax ();
	umin = mcurve.getUmin ();
	CPoint pnt (mcurve.getPointAtParam (umax).getX () + maxl.getX (),
		mcurve.getPointAtParam (umax).getY () + maxl.getY (),
		mcurve.getPointAtParam (umax).getZ () + maxl.getZ (),
		1,1);
	getParamAtPointBasis (mcurve.getPointAtParam (umax),up,v1);
	getParamAtPointBasis (pnt,vp,v2);
	mRSurface = new Geom_RectangularTrimmedSurface (mSurface,u1,u2,v1,v2);
	mcurve.setRotation (p,ax,ang,mcurve);
	gp_Dir dir;
	dir = mSurface->Direction ();
	maxl.setX (dir.X ());
	maxl.setY (dir.Y ());
	maxl.setZ (dir.Z ());
}
int SurfaceOfLinearExtrusion::getSurfaceType (){
	return LINESURFACE;
}
Handle_Geom_RectangularTrimmedSurface SurfaceOfLinearExtrusion::getGeomSurface (){
	return mRSurface;
}
