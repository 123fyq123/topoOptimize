#include "CAidFunction.h"
#include <cmath>
#include "CLine.h"
#include "CPlane.h"
#include <TopTools_ListOfShape.hxx>
#include <BRepBuilderAPI_MakeEdge.hxx>
#include <BRepBuilderAPI_MakeFace.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <IGESControl_Controller.hxx>
#include <IGESControl_Writer.hxx>
#include <GeomAPI_ExtremaCurveCurve.hxx>
#include <Geom_Line.hxx>
#include <IntTools_EdgeEdge.hxx>
#include <BRepExtrema_ExtCC.hxx>
#include <GeomTools.hxx>
#include <GeomTools_CurveSet.hxx>
#include <GeomAPI_ExtremaCurveSurface.hxx>
#include <GeomAPI_IntSS.hxx>
#include <GeomAdaptor_Curve.hxx>
#include<GeomAdaptor_Surface.hxx>
#include <GeomAPI_IntCS.hxx>
#include<GeomAPI_ProjectPointOnSurf.hxx>
#include <BRepBuilderAPI_MakeVertex.hxx>
#include <TopoDS_Vertex.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>
#include <GeomConvert.hxx>
#include <IGESControl_Reader.hxx>
#include <TopTools_HSequenceOfShape.hxx>
#include <TopExp_Explorer.hxx>
#include <TopoDS.hxx>
#include <GeomLib.hxx>
#include <BRep_Tool.hxx>
#include<BndLib_Add3dCurve.hxx>
#include<BndLib_AddSurface.hxx>
CAidFunction::CAidFunction (){
}


CAidFunction::~CAidFunction (){
}

double CAidFunction::dot2D (CPoint &point1,CPoint &point2){
	double x = point1.getX ();
	double y = point1.getY ();
	x *= point2.getX ();
	y *= point2.getY ();
	return  x + y;
}

double CAidFunction::angle2D (CPoint &point1,CPoint &point2,int typeAngle){
	double dot = dot2D (point1,point2);
	double len1 = point1.norm ();
	double len2 = point2.norm ();
	double value = dot / len1 / len2;
	if(typeAngle == 1){
		return acos (value) / PI * 180; //degree
	}
	else if(typeAngle == 0){
		return acos (value);
	}
	else{
		return -1; //错误处理
	}
}

double CAidFunction::distance2D (CPoint &point1,CPoint &point2){
	double x = point1.getX ();
	double y = point1.getY ();

	x -= point2.getX ();
	y -= point2.getY ();

	double len = sqrt (x * x + y * y);
	return len;
}

double CAidFunction::crossProduct2D (CPoint &point1,CPoint &point2){
	double x = point1.getX ();
	double y = point1.getY ();

	x *= point2.getY ();
	y *= point1.getX ();

	return x - y;
}

double CAidFunction::dot3D (CGeomPoint &point1,CGeomPoint &point2){
	double x,y,z;
	point1.getXYZ (x,y,z);
	x *= point2.getX ();
	y *= point2.getY ();
	z *= point2.getZ ();

	return x + y + z;
}

double CAidFunction::angle3D (CPoint &point1,CPoint &point2,int typeAngle){
	double dot = dot3D (point1,point2);
	double len1 = point1.norm ();
	double len2 = point2.norm ();
	double value = dot / len1 / len2;
	if(typeAngle == 1){
		return acos (value) / PI * 180.0; //degree
	}
	else if(typeAngle == 0){
		return acos (value);
	}
	else{
		return -1; //错误处理
	}
}

double CAidFunction::distance3D (CPoint &point1,CPoint &point2){
	double x = point1.getX ();
	double y = point1.getY ();
	double z = point1.getZ ();

	x -= point2.getX ();
	y -= point2.getY ();
	z -= point2.getZ ();

	double len = sqrt (x * x + y * y + z * z);
	return len;
}

void CAidFunction::crossProduct3D (CGeomPoint &point1,CGeomPoint &point2,CGeomPoint &vec){
	double x = point1.getY () * point2.getZ () - point1.getZ () * point2.getY ();
	double y = point1.getZ () * point2.getX () - point1.getX () * point2.getZ ();
	double z = point1.getX () * point2.getY () - point1.getY () * point2.getX ();
	vec.setXYZ (x,y,z);
}

double CAidFunction::intersectLinePlane (CGeomCurve &line,CGeomSurface &plane){
	//if (line.getCurveType() != 1 || plane.getSurfaceType() != 1)
	//{
	//	return -1;
	//}
	//CLine line0 = dynamic_cast<CLine &>(line);
	//CPlane plane0 = dynamic_cast<CPlane &>(plane);
	//CGeomPoint planePoint, planeNormal, linePoint, lineVector;
	CPoint planePoint,linePoint;
	CVector planeNormal,lineVector;
	plane.getPoint (planePoint);
	plane.getNormal (planeNormal);
	line.getPoint (linePoint);
	line.getNormal (lineVector);
	//CGeomPoint planeNormal = plane.getNormal();
	//CGeomPoint linePoint = line.getPoint();
	//CGeomPoint lineVector = line.getNormal();
	double d = dot3D (planePoint,planeNormal);
	double div = dot3D (planeNormal,lineVector);

	if(div == 0){
		div = 1e-8;
	}

	double t = (double) (d - (dot3D (planeNormal,linePoint))) / div; //错误处理，待补充
	return t;
}

bool CAidFunction::getCCIntersection (CGeomCurve &curve1,CGeomCurve &curve2,vector<CPoint> &points){
	int type1 = curve1.getCurveType ();
	int type2 = curve2.getCurveType ();
	if(type1 == 1 && type2 == 1){
		//两直线求交
		CLine tmpLine1 = dynamic_cast<CLine &>(curve1);
		Handle_Geom_Curve line1 = tmpLine1.getGeomCurve ();

		CLine tmpLine2 = dynamic_cast<CLine &>(curve2);
		Handle_Geom_Curve line2 = tmpLine2.getGeomCurve ();

		double umin1 = tmpLine1.getUmin ();
		double umax1 = tmpLine1.getUmax ();
		double umin2 = tmpLine2.getUmin ();
		double umax2 = tmpLine2.getUmax ();

		getCCIntersect (line1,line2,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 1 && type2 == 2) || (type1 == 2 && type2 == 1)){
		//直线和b样条求交
		CLine tmpLine = dynamic_cast<CLine &>(type1 == 1 ? curve1 : curve2);
		Handle_Geom_Curve line = tmpLine.getGeomCurve ();

		CBSplineCurve spline = dynamic_cast<CBSplineCurve &>(type1 == 1 ? curve2 : curve1);
		Handle_Geom_Curve  c1 = spline.getGeomCurve ();
		double umin1 = tmpLine.getUmin ();
		double umax1 = tmpLine.getUmax ();
		double umin2 = spline.getUmin ();
		double umax2 = spline.getUmax ();
		getCCIntersect (line,c1,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 1 && type2 == 3) || (type1 == 3 && type2 == 1)){
		//直线和bezier求交
		CLine tmpLine = dynamic_cast<CLine &>(type1 == 1 ? curve1 : curve2);
		Handle_Geom_Curve line = tmpLine.getGeomCurve ();

		CBezierCurve bezier = dynamic_cast<CBezierCurve &>(type1 == 1 ? curve2 : curve1);
		Handle_Geom_Curve  c1 = bezier.getBezierCurve ();
		double umin1 = tmpLine.getUmin ();
		double umax1 = tmpLine.getUmax ();
		double umin2 = 0;
		double umax2 = 1;

		getCCIntersect (line,c1,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 2 && type2 == 2)){
		//两个b样条求交
		CBSplineCurve spline1 = dynamic_cast<CBSplineCurve &>(curve1);
		CBSplineCurve spline2 = dynamic_cast<CBSplineCurve &>(curve2);
		Handle_Geom_Curve c1 = spline1.getGeomCurve ();
		Handle_Geom_Curve c2 = spline2.getGeomCurve ();
		double umin1 = spline1.getUmin ();
		double umax1 = spline1.getUmax ();
		double umin2 = spline2.getUmin ();
		double umax2 = spline2.getUmax ();
		getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 2 && type2 == 3) || (type1 == 3 && type2 == 2)){
		//b样条和bezier求交
		CBSplineCurve spline = dynamic_cast<CBSplineCurve &>(curve1);
		CBezierCurve bezier = dynamic_cast<CBezierCurve &>(curve2);
		Handle_Geom_Curve c1 = spline.getGeomCurve ();
		Handle_Geom_Curve c2 = bezier.getBezierCurve ();
		double umin1 = spline.getUmin ();
		double umax1 = spline.getUmax ();
		double umin2 = 0;
		double umax2 = 1;
		getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
	}
	else if(type1 == 3 && type2 == 3){
		//两条bezier求交
		CBezierCurve bezier1 = dynamic_cast<CBezierCurve &>(curve1);
		CBezierCurve bezier2 = dynamic_cast<CBezierCurve &>(curve2);
		Handle_Geom_Curve c1 = bezier1.getBezierCurve ();
		Handle_Geom_Curve c2 = bezier2.getBezierCurve ();
		int umin1 = 0;
		int umax1 = 1;
		int umin2 = 0;
		int umax2 = 1;
		getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);

	}

	return true;
}

bool CAidFunction::getCCIntersection (CGeomCurve &curve1,CGeomCurve &curve2,vector<CPoint> &points,double umin1,double umax1,double umin2,double umax2){
	int type1 = curve1.getCurveType ();
	int type2 = curve2.getCurveType ();
	bool res = false;
	if(type1 == 1 && type2 == 1){
		//两直线求交
		CLine tmpLine1 = dynamic_cast<CLine &>(curve1);
		Handle_Geom_Curve line1 = tmpLine1.getGeomCurve ();

		CLine tmpLine2 = dynamic_cast<CLine &>(curve2);
		Handle_Geom_Curve line2 = tmpLine2.getGeomCurve ();

		res = getCCIntersect (line1,line2,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 1 && type2 == 2) || (type1 == 2 && type2 == 1)){
		//直线和b样条求交
		CLine tmpLine = dynamic_cast<CLine &>(type1 == 1 ? curve1 : curve2);
		Handle_Geom_Curve line = tmpLine.getGeomCurve ();

		CBSplineCurve spline = dynamic_cast<CBSplineCurve &>(type1 == 1 ? curve2 : curve1);
		Handle_Geom_Curve  c1 = spline.getGeomCurve ();

		res = getCCIntersect (line,c1,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 1 && type2 == 3) || (type1 == 3 && type2 == 1)){
		//直线和bezier求交
		CLine tmpLine = dynamic_cast<CLine &>(type1 == 1 ? curve1 : curve2);
		Handle_Geom_Curve line = tmpLine.getGeomCurve ();

		CBezierCurve bezier = dynamic_cast<CBezierCurve &>(type1 == 1 ? curve2 : curve1);
		Handle_Geom_Curve  c1 = bezier.getBezierCurve ();

		res = getCCIntersect (line,c1,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 2 && type2 == 2)){
		//两个b样条求交
		CBSplineCurve spline1 = dynamic_cast<CBSplineCurve &>(curve1);
		CBSplineCurve spline2 = dynamic_cast<CBSplineCurve &>(curve2);
		Handle_Geom_Curve c1 = spline1.getGeomCurve ();
		Handle_Geom_Curve c2 = spline2.getGeomCurve ();
		res = getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
	}
	else if((type1 == 2 && type2 == 3) || (type1 == 3 && type2 == 2)){
		if(type1 == 2 && type2 == 3){
			//b样条和bezier求交
			CBSplineCurve spline = dynamic_cast<CBSplineCurve &>(curve1);
			CBezierCurve bezier = dynamic_cast<CBezierCurve &>(curve2);
			Handle_Geom_Curve c1 = spline.getGeomCurve ();
			Handle_Geom_Curve c2 = bezier.getBezierCurve ();
			double umin1 = spline.getUmin ();
			double umax1 = spline.getUmax ();
			double umin2 = 0;
			double umax2 = 1;
			res = getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
		}
		else{
			CBSplineCurve spline = dynamic_cast<CBSplineCurve &>(curve2);
			CBezierCurve bezier = dynamic_cast<CBezierCurve &>(curve1);
			Handle_Geom_Curve c2 = spline.getGeomCurve ();
			Handle_Geom_Curve c1 = bezier.getBezierCurve ();
			double umin2 = spline.getUmin ();
			double umax2 = spline.getUmax ();
			double umin1 = 0;
			double umax1 = 1;
			res= getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
		}
	}
	else if(type1 == 3 && type2 == 3){
		//两条bezier求交
		CBezierCurve bezier1 = dynamic_cast<CBezierCurve &>(curve1);
		CBezierCurve bezier2 = dynamic_cast<CBezierCurve &>(curve2);
		Handle_Geom_Curve c1 = bezier1.getBezierCurve ();
		Handle_Geom_Curve c2 = bezier2.getBezierCurve ();
		int umin1 = 0;
		int umax1 = 1;
		int umin2 = 0;
		int umax2 = 1;
		res = getCCIntersect (c1,c2,points,umin1,umax1,umin2,umax2);
	}
	return res;
}

bool CAidFunction::getCSIntersection (CGeomCurve &curve,CGeomSurface &surface,vector<CPoint> &points){
	int type1 = curve.getCurveType ();
	int type2 = surface.getSurfaceType ();
	bool res = false;
	if(type1 == 1 && type2 == 1){
		//直线和平面求交
		CLine tmpLine1 = dynamic_cast<CLine &>(curve);
		Handle_Geom_Curve line = tmpLine1.getGeomCurve ();

		CPlane plane = dynamic_cast<CPlane &>(surface);
		Handle_Geom_Surface planeSurface = plane.getGeomPlane ();

		double u1,u2,umin,umax,vmin,vmax;
		u1 = tmpLine1.getUmin ();
		u2 = tmpLine1.getUmax ();
		plane.getBounds (umin,umax,vmin,vmax);

		res = getCSIntersect (line,planeSurface,points,u1,u2,umin,umax,vmin,vmax);
	}
	else if((type1 == 1 && type2 == 2)){
		//直线和b样条曲面求交
		CLine tmpLine1 = dynamic_cast<CLine &>(curve);
		Handle_Geom_Curve line = tmpLine1.getGeomCurve ();

		CBSplineSurface splineSurface = dynamic_cast<CBSplineSurface &>(surface);
		Handle_Geom_Surface s = splineSurface.getGeomSurface ();

		double umin,umax,vmin,vmax,u1,u2;
		splineSurface.getUVMinMax (umin,umax,vmin,vmax);
		u1 = tmpLine1.getUmin ();
		u2 = tmpLine1.getUmax ();
		res = getCSIntersect (line,s,points,u1,u2,umin,umax,vmin,vmax);

	}
	else if((type1 == 1 && type2 == 3)){
		//直线和bezier曲面求交
		CLine tmpLine1 = dynamic_cast<CLine &>(curve);
		Handle_Geom_Curve line = tmpLine1.getGeomCurve ();

		CBezierSurface bezierSurface = dynamic_cast<CBezierSurface &>(surface);
		Handle_Geom_Surface s = bezierSurface.getGeomSurface ();

		double u1,u2;
		u1 = tmpLine1.getUmin ();
		u2 = tmpLine1.getUmax ();
		res = getCSIntersect (line,s,points,u1,u2,0,1,0,1);

	}
	else if((type1 == 2 && type2 == 2) || (type1 == 2 && type2 == 3)){
		//b样条曲线和b样条/bezier曲面求交
		CBSplineCurve splineCurve = dynamic_cast<CBSplineCurve &>(curve);
		CBSplineSurface splineSurface;
		CBezierSurface bezierSurface;
		if((type1 == 2 && type2 == 2)){
			splineSurface = dynamic_cast<CBSplineSurface &>(surface);
			Handle_Geom_Curve c = splineCurve.getGeomCurve ();
			Handle_Geom_Surface s = splineSurface.getGeomSurface ();
			double umin,umax,vmin,vmax,u1,u2;
			splineSurface.getUVMinMax (umin,umax,vmin,vmax);
			u1 = splineCurve.getUmin ();
			u2 = splineCurve.getUmax ();
			res = getCSIntersect (c,s,points,u1,u2,umin,umax,vmin,vmax);
		}
		if(type1 == 2 && type2 == 3){
			bezierSurface = dynamic_cast<CBezierSurface &>(surface);
			Handle_Geom_Curve c = splineCurve.getGeomCurve ();
			Handle_Geom_Surface s = bezierSurface.getGeomSurface ();
			double u1,u2;
			u1 = splineCurve.getUmin ();
			u2 = splineCurve.getUmax ();
			res = getCSIntersect (c,s,points,u1,u2,0,1,0,1);
		}

	}
	else if((type1 == 2 && type2 == 1)){
		//b样条和平面求交
		CBSplineCurve splineCurve = dynamic_cast<CBSplineCurve &>(curve);
		CPlane plane;
		plane = dynamic_cast<CPlane &>(surface);
		Handle_Geom_Curve c = splineCurve.getGeomCurve ();
		Handle_Geom_Surface s = plane.getGeomPlane ();
		double u1,u2,umin,umax,vmin,vmax;
		u1 = splineCurve.getUmin ();
		u2 = splineCurve.getUmax ();
		plane.getBounds (umin,umax,vmin,vmax);
		res = getCSIntersect (c,s,points,u1,u2,umin,umax,vmin,vmax);
	}
	else if(type1 == 3 && type2 == 3 || (type1 == 3 && type2 == 2)){
		//bezier和b样条/bezier曲面求交  要分开  有问题
		CBezierCurve bezierCurve = dynamic_cast<CBezierCurve &>(curve);
		CBSplineSurface splineSurface;
		CBezierSurface bezierSurface;
		if((type1 == 3 && type2 == 2)){
			splineSurface = dynamic_cast<CBSplineSurface &>(surface);
			Handle_Geom_Curve c = bezierCurve.getBezierCurve ();
			Handle_Geom_Surface s = splineSurface.getGeomSurface ();
			double umin,umax,vmin,vmax,u1,u2;
			splineSurface.getUVMinMax (umin,umax,vmin,vmax);
			u1 = 0;
			u2 = 1;
			res = getCSIntersect (c,s,points,u1,u2,umin,umax,vmin,vmax);
		}
		if(type1 == 3 && type2 == 3){
			bezierSurface = dynamic_cast<CBezierSurface &>(surface);
			Handle_Geom_Curve c = bezierCurve.getBezierCurve ();
			Handle_Geom_Surface s = bezierSurface.getGeomSurface ();
			res = getCSIntersect (c,s,points,0,1,0,1,0,1);
		}
	}
	else if(type1 == 3 && type2 == 1){
		//bezier和平面求交
		CBezierCurve bezierCurve = dynamic_cast<CBezierCurve &>(curve);
		CPlane plane;
		plane = dynamic_cast<CPlane &>(surface);
		Handle_Geom_Curve c = bezierCurve.getBezierCurve ();
		Handle_Geom_Surface s = plane.getGeomPlane ();
		double u1,u2,umin,umax,vmin,vmax;
		u1 = 0;
		u2 = 1;
		plane.getBounds (umin,umax,vmin,vmax);
		res = getCSIntersect (c,s,points,u1,u2,umin,umax,vmin,vmax);
	}

	return res;
}

bool CAidFunction::getSSIntersection (CGeomSurface &surface1,CGeomSurface &surface2,vector<CBSplineCurve> &curves,double tol){
	int type1 = surface1.getSurfaceType ();
	int type2 = surface2.getSurfaceType ();
	Handle_Geom_Surface s1;
	Handle_Geom_Surface s2;
	if(type1 == 1 && type2 == 1){
		//平面和平面求交
		CPlane plane1 = dynamic_cast<CPlane &>(surface1);
		CPlane plane2 = dynamic_cast<CPlane &>(surface2);
		s1 = plane1.getGeomPlane ();
		s2 = plane2.getGeomPlane ();
		//getSSIntersect(s1, s2, curves, tol);
	}
	else if((type1 == 1 && type2 == 2) || (type1 == 2 && type2 == 1)){
		//平面和b样条曲面求交
		CBSplineSurface spline1 = dynamic_cast<CBSplineSurface &>(type1 == 1 ? surface2 : surface1);
		CPlane plane = dynamic_cast<CPlane &>(type1 == 1 ? surface1 : surface2);
		s1 = spline1.getGeomSurface ();
		s2 = plane.getGeomPlane ();
		//getSSIntersect (s1,s2,curves,tol);
	}
	else if((type1 == 1 && type2 == 3) || (type1 == 3 && type2 == 1)){
		//平面和bezier曲面求交
		CBezierSurface bezierface = dynamic_cast<CBezierSurface &>(type1 == 1 ? surface2 : surface1);
		CPlane plane = dynamic_cast<CPlane &>(type1 == 1 ? surface1 : surface2);
		s1 = bezierface.getGeomSurface ();
		s2 = plane.getGeomPlane ();
		//getSSIntersect (s1,s2,curves,tol);
	}
	else if((type1 == 2 && type2 == 2) || (type1 == 2 && type2 == 3)){
		//b样条曲面和b样条/bezier曲面求交
		CBSplineSurface spline1 = dynamic_cast<CBSplineSurface &>(surface1);

		if((type1 == 2 && type2 == 2)){
			CBSplineSurface spline2;
			spline2 = dynamic_cast<CBSplineSurface &>(surface2);
			s1 = spline1.getGeomSurface ();
			s2 = spline2.getGeomSurface ();
			//getSSIntersect (s1,s2,curves,tol);
		}
		else{
			CBezierSurface bezierface;
			bezierface = dynamic_cast<CBezierSurface &>(surface2);
			s1 = spline1.getGeomSurface ();
			s2 = bezierface.getGeomSurface ();
			//getSSIntersect (s1,s2,curves,tol);
		}
	}
	else if(type1 == 3 && type2 == 3 || (type1 == 3 && type2 == 2)){
		CBezierSurface bs1 = dynamic_cast<CBezierSurface &>(surface1);
		if(type1 == 3 && type2 == 3){
			CBezierSurface bezierface;
			bezierface = dynamic_cast<CBezierSurface &>(surface2);
			s1 = bs1.getGeomSurface ();
			s2 = bezierface.getGeomSurface ();
			//getSSIntersect (s1,s2,curves,tol);
		}
		//bezier和b样条/bezier曲面求交
		else{
			CBSplineSurface spline2;
			spline2 = dynamic_cast<CBSplineSurface &>(surface2);
			s1 = bs1.getGeomSurface ();
			s2 = spline2.getGeomSurface ();
			//getSSIntersect (s1,s2,curves,tol);
		}
	}
	else if (type1 == 4 && type2 == 4) {
		//圆柱和圆柱求交
		CCylindricalSurface ccy1 = dynamic_cast<CCylindricalSurface &>(surface1);
		CCylindricalSurface ccy2 = dynamic_cast<CCylindricalSurface &>(surface2);
		s1 = ccy1.getGeomSurface();
		s2 = ccy2.getGeomSurface();
		//getSSIntersect(s1, s2, curves, tol);
	}
	else if ((type1 == 4 && type2 == 2) || (type1 == 2 && type2 == 4)) {
		//圆柱和b样条曲面求交
		CBSplineSurface spline1 = dynamic_cast<CBSplineSurface &>(type1 == 4 ? surface2 : surface1);
		CCylindricalSurface ccy = dynamic_cast<CCylindricalSurface &>(type1 == 4 ? surface1 : surface2);
		s1 = spline1.getGeomSurface();
		s2 = ccy.getGeomSurface();
		//getSSIntersect(s1, s2, curves, tol);

	}
	else if ((type1 == 4 && type2 == 3) || (type1 == 3 && type2 == 4)) {
		//圆柱和bezier曲面求交
		CBezierSurface bezierface = dynamic_cast<CBezierSurface &>(type1 == 4 ? surface2 : surface1);
		CCylindricalSurface ccy = dynamic_cast<CCylindricalSurface &>(type1 == 4 ? surface1 : surface2);
		s1 = bezierface.getGeomSurface();
		s2 = ccy.getGeomSurface();
		//getSSIntersect(s1, s2, curves, tol);
	}

	bool res = getSSIntersect(s1, s2, curves, tol);
	return res;
}
//支持b样条曲面
void CAidFunction::geomNormal (CGeomSurface &surface,double u,double v,CVector &normal){
	CBSplineSurface bsplinesurface = dynamic_cast<CBSplineSurface &>(surface);
	Handle_Geom_BSplineSurface geomsurface = bsplinesurface.getGeomSurface ();
	gp_Pnt2d point2d (u,v);
	gp_Dir dir;
	GeomLib::NormEstim (geomsurface,point2d,1.e-10,dir);
	normal.setX (dir.X ());
	normal.setY (dir.Y ());
	normal.setZ (dir.Z ());

}

void CAidFunction::addCurveToBox (CGeomCurve& curve,double tol,CBndBox& box) {
	int type = curve.getCurveType ();
	Handle_Geom_Curve pCurve;
	if(type == 1){
		pCurve = (dynamic_cast<CLine &>(curve)).getGeomCurve ();
	}
	else if(type == 2){
		pCurve = (dynamic_cast<CBSplineCurve &>(curve)).getGeomCurve ();
	}
	else if(type == 3){
		pCurve = (dynamic_cast<CBezierCurve &>(curve)).getBezierCurve ();
	}
	else if(type == 4||type==5){
		pCurve = (dynamic_cast<CCircle &>(curve)).getGeomCurve ();
	}
	GeomAdaptor_Curve adcurve (pCurve);
	BndLib_Add3dCurve addcurve;
	Bnd_Box mbox;
	addcurve.Add (adcurve,tol,mbox);
	box.setmBndBox (mbox);
}

void CAidFunction::addSurfaceToBox (CGeomSurface& surface,double tol,CBndBox& box) {

	int type = surface.getSurfaceType ();
	Handle_Geom_Surface pSurface;
	if(type == 1){
		pSurface = (dynamic_cast<CPlane &>(surface)).getGeomPlane ();
		
	}
	else if(type == 2){
		pSurface = (dynamic_cast<CBSplineSurface &>(surface)).getGeomSurface ();
		
	}
	else if(type == 3){
		pSurface = (dynamic_cast<CBezierSurface &>(surface)).getGeomSurface ();
		
	}
	else if(type == 4){
		pSurface = (dynamic_cast<CCylindricalSurface &>(surface)).getGeomSurface ();
		
	}
	else if(type == 5){
		pSurface = (dynamic_cast<SurfaceOfLinearExtrusion &>(surface)).getGeomSurface ();
	}
	else if(type == 6){
		pSurface = (dynamic_cast<SurfaceOfRevolution &>(surface)).getGeomSurface ();
		
	}
	else if(type == 7){
		pSurface = (dynamic_cast<CSphere &>(surface)).getGeomSurface ();	
	}

	GeomAdaptor_Surface adsurface (pSurface);
	BndLib_AddSurface addsurface;
	Bnd_Box mbox;
	addsurface.Add (adsurface,tol,mbox);
	box.setmBndBox (mbox);
}

void CAidFunction::getIGESCurve (string fileName,CGeomCurve&curve){
	vector<CBSplineCurve>bscurves;
	CBSplineCurve bscurve;
	adaptorCurve (curve,bscurve);
	
	/*int type = curve.getCurveType ();
	Handle_Geom_Curve pCurve;
	if(type == 1){
		pCurve = (dynamic_cast<CLine &>(curve)).getGeomCurve ();
	}
	else if(type == 2){
		pCurve = (dynamic_cast<CBSplineCurve &>(curve)).getGeomCurve ();
	}
	else if(type == 3){
		pCurve = (dynamic_cast<CBezierCurve &>(curve)).getBezierCurve ();
	}
	else if(type == 4){
		pCurve = (dynamic_cast<CCircle &>(curve)).getGeomCurve ();
	}
	else if(type == 5){
		pCurve = (dynamic_cast<CCircle &>(curve)).getGeomCurve ();
	}*/
	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;
	igesWriter.AddGeom (bscurve.getGeomCurve());
	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

void CAidFunction::getIGESCurve (string fileName,vector<CBSplineCurve> curves){
	
	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;

	size_t size = curves.size ();

	for(int i = 0; i < size; ++i){
		/*int type = curves[i]->getCurveType ();
		if(type == 1){
			Handle_Geom_Curve pCurve = (dynamic_cast<CLine *>(curves[i]))->getGeomCurve ();
			igesWriter.AddGeom (pCurve);
		}
		else if(type == 2){
			Handle_Geom_Curve pCurve = (dynamic_cast<CBSplineCurve *>(curves[i]))->getGeomCurve ();
			igesWriter.AddGeom (pCurve);
		}
		else if(type == 3){
			Handle_Geom_Curve pCurve = (dynamic_cast<CBezierCurve *>(curves[i]))->getBezierCurve ();
			igesWriter.AddGeom (pCurve);
		}
		else if(type == 4 || type == 5){
			Handle_Geom_Curve pCurve = (dynamic_cast<CCircle *>(curves[i]))->getGeomCurve ();
			igesWriter.AddGeom (pCurve);
		}*/
		igesWriter.AddGeom(curves[i].getGeomCurve());
	}
	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

void CAidFunction::getIGESSurface (string fileName,CGeomSurface &surface){
	/*int type = surface.getSurfaceType ();
	Handle_Geom_Surface pSurface;*/
	CBSplineSurface bssurface;
	adaptorSurface (surface,bssurface);
	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;
	/*if(type == 1){
		pSurface = (dynamic_cast<CPlane &>(surface)).getGeomPlane ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 2){
		pSurface = (dynamic_cast<CBSplineSurface &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 3){
		pSurface = (dynamic_cast<CBezierSurface &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 4){
		pSurface = (dynamic_cast<CCylindricalSurface &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 5){
		pSurface = (dynamic_cast<SurfaceOfLinearExtrusion &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 6){
		pSurface = (dynamic_cast<SurfaceOfRevolution &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}
	else if(type == 7){
		pSurface = (dynamic_cast<CSphere &>(surface)).getGeomSurface ();
		igesWriter.AddGeom (pSurface);
	}*/
	igesWriter.AddGeom(bssurface.getGeomSurface());
	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

void CAidFunction::getIGESSurface (string fileName,vector<CBSplineSurface > &surfaces){
	//int type = surfaces[0]->getSurfaceType();

	//int type = surfaces[0]->getSurfaceType();

	IGESControl_Controller::Init();
	IGESControl_Writer igesWriter;

	size_t size = surfaces.size();

	for (int i = 0; i < size; ++i) {
		/*int type = surfaces[i]->getSurfaceType ();
		if(type == 1){
			Handle_Geom_Surface pSurface = (dynamic_cast<CPlane *>(surfaces[i]))->getGeomPlane ();
			igesWriter.AddGeom (pSurface);
		}
		else if(type == 2){
			Handle_Geom_Surface pSurface = (dynamic_cast<CBSplineSurface *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}

		else if(type == 3){
			Handle_Geom_Surface pSurface = (dynamic_cast<CBezierSurface *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}

		else if(type == 4){
			Handle_Geom_Surface pSurface = (dynamic_cast<CCylindricalSurface *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}
		else if(type == 5){
			Handle_Geom_Surface pSurface = (dynamic_cast<SurfaceOfLinearExtrusion *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}
		else if(type == 6){
			Handle_Geom_Surface pSurface = (dynamic_cast<SurfaceOfRevolution *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}
		else if(type == 7){
			Handle_Geom_Surface pSurface = (dynamic_cast<CSphere *>(surfaces[i]))->getGeomSurface ();
			igesWriter.AddGeom (pSurface);
		}*/
		igesWriter.AddGeom(surfaces[i].getGeomSurface());
	}

	igesWriter.ComputeModel();
	igesWriter.Write(fileName.c_str());

}

void CAidFunction::getIGESPoint (string fileName,CPoint &point){
	//gp_Pnt p = dynamic_cast<CPoint &>(point).getPnt ();
	gp_Pnt p = point.getPnt();
	TopoDS_Vertex v = BRepBuilderAPI_MakeVertex (p);
	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;
	igesWriter.AddShape (v);

	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

void CAidFunction::getIGESPoint (string fileName,vector<CPoint > &points){
	size_t size = points.size ();
	vector<TopoDS_Vertex> vs;
	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;

	for(int i = 0; i < size; ++i){
		//gp_Pnt p = (dynamic_cast<CPoint *>(points[i]))->getPnt ();
		gp_Pnt p =points[i].getPnt();
		TopoDS_Vertex v = BRepBuilderAPI_MakeVertex (p);
		igesWriter.AddShape (v);
	}

	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

void CAidFunction::getIGESPntCurSur (string fileName,vector<CPoint > &points,vector<CBSplineCurve > &curves,vector<CBSplineSurface > &surfaces){
	/*对于vector集合，必须使用指针类型才能转换成功，但是这样会导致每次导出文件时需要转换成指针数组*/
	/*vector<CBSplineCurve>bscurves;
	vector<CBSplineSurface>bssurfaces;
	for (int i = 0; i < curves.size (); i++) {
		CBSplineCurve curve;
		adaptorCurve (curves [i],curve);
		bscurves.push_back (curve);
	}
	for (int j = 0; j < surfaces.size (); j++) {
		CBSplineSurface surface;
		adaptorSurface (surfaces [j],surface);
		bssurfaces.push_back (surface);
	}*/

	IGESControl_Controller::Init ();
	IGESControl_Writer igesWriter;
	size_t size = points.size ();
	for(int i = 0; i < size; ++i){
		//gp_Pnt p = (dynamic_cast<CPoint *>(points[i]))->getPnt ();
		gp_Pnt p =points[i].getPnt();
		TopoDS_Vertex v = BRepBuilderAPI_MakeVertex (p);
		igesWriter.AddShape (v);
	}

	size = curves.size ();
	for(int i = 0; i < size; ++i){
		//int type = curves[i]->getCurveType ();
		TopoDS_Edge edge;
		/*if(type == 1){
			Handle_Geom_Curve pCurve = (dynamic_cast<CLine *>(curves[i]))->getGeomCurve ();
			edge = BRepBuilderAPI_MakeEdge (pCurve);
		}
		else if(type == 2){
			Handle_Geom_Curve pCurve = (dynamic_cast<CBSplineCurve *>(curves[i]))->getGeomCurve ();
			edge = BRepBuilderAPI_MakeEdge (pCurve);
		}

		else if(type == 3){
			Handle_Geom_Curve pCurve = (dynamic_cast<CBezierCurve *>(curves[i]))->getBezierCurve ();
			edge = BRepBuilderAPI_MakeEdge (pCurve);
		}*/
		edge = BRepBuilderAPI_MakeEdge(curves[i].getGeomCurve());
		igesWriter.AddShape (edge);
	}

	size = surfaces.size ();
	for(int i = 0; i < size; ++i){
		//int type = surfaces[i]->getSurfaceType ();
		TopoDS_Face face;
		/*if(type == 1){
			Handle_Geom_Surface pSurface = (dynamic_cast<CPlane *>(surfaces[i]))->getGeomPlane ();
			face = BRepBuilderAPI_MakeFace (pSurface,1.e-15);
		}
		else if(type == 2){
			Handle_Geom_Surface pSurface = (dynamic_cast<CBSplineSurface *>(surfaces[i]))->getGeomSurface ();
			face = BRepBuilderAPI_MakeFace (pSurface,1.e-15);
		}

		else if(type == 3){
			Handle_Geom_Surface pSurface = (dynamic_cast<CBezierSurface *>(surfaces[i]))->getGeomSurface ();
			face = BRepBuilderAPI_MakeFace (pSurface,1.e-15);
		}
		else if(type == 4){
			Handle_Geom_Surface pSurface = (dynamic_cast<CCylindricalSurface *>(surfaces[i]))->getGeomSurface ();
			face = BRepBuilderAPI_MakeFace (pSurface,1.e-15);
		}*/
		face = BRepBuilderAPI_MakeFace(surfaces[i].getGeomSurface(),1.e-15);
		igesWriter.AddShape (face);
	}

	igesWriter.ComputeModel ();
	igesWriter.Write (fileName.c_str ());
}

bool CAidFunction::readIGESGeomPoint(string & fileName, vector<CPoint>& points)
{
	vector<CBSplineCurve > curves;
	vector<CBSplineSurface > surfaces;
	bool res = getIGESGeom(0, fileName, points, curves, surfaces);
	return res;
}

bool CAidFunction::readIGESGeomCurve(string & fileName, vector<CBSplineCurve>& curves)
{
	vector<CPoint > pnts;
	vector<CBSplineSurface > surfaces;
	bool res = getIGESGeom(1, fileName, pnts, curves, surfaces);
	return res;
}

bool CAidFunction::readIGESGeomSurface(string & fileName, vector<CBSplineSurface>& surfaces)
{
	vector<CPoint > pnts;
	vector<CBSplineCurve > curves;
	bool res = getIGESGeom(2, fileName, pnts, curves, surfaces);
	return res;
}

void CAidFunction::getDumpCurve (CGeomCurve &curve){
	int type = curve.getCurveType ();
	Handle_Geom_Curve pCurve;
	if(type == 1){
		pCurve = (dynamic_cast<CLine &>(curve)).getGeomCurve ();
	}
	else if(type == 2){
		pCurve = (dynamic_cast<CBSplineCurve &>(curve)).getGeomCurve ();
	}
	else if(type == 3){
		pCurve = (dynamic_cast<CBezierCurve &>(curve)).getBezierCurve ();
	}
	else if(type == 4 || type == 5){
		pCurve = (dynamic_cast<CCircle &>(curve)).getGeomCurve ();
	}


	//std::ofstream dumpFile(fileName);
	//GeomTools::Write(pCurve, dumpFile);
	//GeomTools::Dump(pCurve, dumpFile);
	GeomTools::Dump (pCurve,std::cout);
}

void CAidFunction::getDumpSurface (CGeomSurface &surface){
	int type = surface.getSurfaceType ();
	Handle_Geom_Surface pSurface;
	if(type == 1){
		pSurface = (dynamic_cast<CPlane &>(surface)).getGeomPlane ();
	}
	else if(type == 2){
		pSurface = (dynamic_cast<CBSplineSurface &>(surface)).getGeomSurface ();
	}
	else if(type == 3){
		pSurface = (dynamic_cast<CBezierSurface &>(surface)).getGeomSurface ();
	}

	GeomTools::Dump (pSurface,std::cout);
}

//void CAidFunction::readGeomTurbo (string filePath,string fileName,int typeSurface,vector<CBSplineSurface *> &surfaces){
//	std::ifstream infile (filePath);
//	vector<vector<CPoint>> points1,points2;
//	vector<CPoint> cyPoints1,cyPoints2;
//	int num;
//	string s0,s1,s2,s3;
//	double n0;
//	int n1;
//	double n2,n3;
//	infile >> s0 >> s1;
//	infile >> s0 >> n0;
//	infile >> s0 >> n0;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	//ÇúÏßhubÊý¾Ý´¦Àí
//	infile >> s0 >> s1;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; i++){
//		infile >> n2 >> n3;
//		CPoint pnt (n2,n3,0,1,1);
//		cyPoints1.push_back (pnt);
//	}
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	//ÇúÏßshroudÊý¾Ý´¦Àí
//	infile >> s0 >> s1;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; i++){
//		infile >> n2 >> n3;
//		CPoint pnt (n2,n3,0,1,1);
//		cyPoints2.push_back (pnt);
//	}
//
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	//channel
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//
//	//row begin
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> n1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//
//	//niblade
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//
//	infile >> s0 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//
//	infile >> s0 >> n2;
//	infile >> s0 >> n2;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n2;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; ++i){
//		infile >> s0 >> s1 >> n1;
//		infile >> s0;
//		infile >> n1;
//		double x,y,z;
//		vector<CPoint> ps;
//		for(int j = 1; j <= n1; ++j){
//			infile >> x >> y >> z;
//			CPoint pnt (x,y,z,1,1);
//			ps.push_back (pnt);
//		}
//
//		points1.push_back (ps);
//	}
//
//	infile >> s0;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; ++i){
//		infile >> s0 >> s1 >> n1;
//		infile >> s0;
//		infile >> n1;
//		double x,y,z;
//		vector<CPoint> ps;
//		for(int j = 1; j <= n1; ++j){
//			infile >> x >> y >> z;
//			CPoint pnt (x,y,z,1,1);
//			ps.push_back (pnt);
//		}
//
//		points2.push_back (ps);
//	}
//
//	infile.close ();
//	//¹¹ÔìÔ²ÖùÃæ
//	double zmin1 = 0,zmax1 = 0,zmin2 = 0,zmax2 = 0;
//	for(vector<CPoint>::iterator it = cyPoints1.begin (); it != cyPoints1.end (); ++it){
//		if((zmin1 - (*it).getX ()) > 1.e-10){
//			zmin1 = (*it).getX ();
//		}
//		if(((*it).getX () - zmax1) > 1.e-10){
//			zmax1 = (*it).getX ();
//		}
//	}
//
//	for(vector<CPoint>::iterator it = cyPoints2.begin (); it != cyPoints2.end (); ++it){
//		if((zmin2 - (*it).getX ()) > 1.e-10){
//			zmin2 = (*it).getX ();
//		}
//		if(((*it).getX () - zmax2) > 1.e-10){
//			zmax2 = (*it).getX ();
//		}
//	}
//
//	CCylindricalSurface *pSurface1 = new CCylindricalSurface (CPoint (0,0,zmin1,1,1),CVector (0,0,zmax1 - zmin1),cyPoints1[0].getY ());
//	CCylindricalSurface *pSurface2 = new CCylindricalSurface (CPoint (0,0,zmin2,1,1),CVector (0,0,zmax2 - zmin2),cyPoints2[0].getY ());
//
//	if(!(pSurface1->getGeomSurface ().IsNull ())){
//		surfaces.push_back (pSurface1);
//	}
//	if(!(pSurface2->getGeomSurface ().IsNull ())){
//		surfaces.push_back (pSurface2);
//	}
//
//	if(typeSurface == 2){
//		CBSplineSurface *pSurface3 = new(CBSplineSurface)(points1);
//		CBSplineSurface *pSurface4 = new(CBSplineSurface)(points2);
//
//		if(!(pSurface3->getGeomSurface ().IsNull ())){
//			surfaces.push_back (pSurface3);
//		}
//		if(!(pSurface4->getGeomSurface ().IsNull ())){
//			surfaces.push_back (pSurface4);
//		}
//	}
//
//	else if(typeSurface == 3){
//		CBezierSurface *pSurface3 = new(CBezierSurface)(points1);
//		CBezierSurface *pSurface4 = new(CBezierSurface)(points2);
//
//		if(!(pSurface3->getGeomSurface ().IsNull ())){
//			surfaces.push_back (pSurface3);
//		}
//		if(!(pSurface4->getGeomSurface ().IsNull ())){
//			surfaces.push_back (pSurface4);
//		}
//	}
//
//	getIGESSurface (fileName,surfaces);
//
//}

//void CAidFunction::readGeomTurbo (string filePath,string fileName,CCylindricalSurface &surface1,CCylindricalSurface &surface2,CBSplineSurface &surface3,CBSplineSurface &surface4){
//	std::ifstream infile (filePath);
//	vector<vector<CPoint>> points1,points2;
//	vector<CPoint> cyPoints1,cyPoints2;
//	int num;
//	string s0,s1,s2,s3;
//	double n0;
//	int n1;
//	double n2,n3;
//	infile >> s0 >> s1;
//	infile >> s0 >> n0;
//	infile >> s0 >> n0;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	//ÇúÏßhubÊý¾Ý´¦Àí
//	infile >> s0 >> s1;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; i++){
//		infile >> n2 >> n3;
//		CPoint pnt (n2,n3,0,1,1);
//		cyPoints1.push_back (pnt);
//	}
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	//ÇúÏßshroudÊý¾Ý´¦Àí
//	infile >> s0 >> s1;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; i++){
//		infile >> n2 >> n3;
//		CPoint pnt (n2,n3,0,1,1);
//		cyPoints2.push_back (pnt);
//	}
//
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//	//channel
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//
//	//row begin
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> n1;
//	infile >> s0 >> s1;
//	infile >> s0 >> n1;
//
//	//niblade
//	infile >> s0 >> s1;
//	infile >> s0 >> s1 >> s2;
//	infile >> s0 >> s1;
//	infile >> s0 >> s1;
//
//	infile >> s0 >> n1;
//	infile >> s0 >> s1 >> s2 >> n1;
//
//	infile >> s0 >> n2;
//	infile >> s0 >> n2;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n2;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0 >> n1;
//	infile >> s0;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; ++i){
//		infile >> s0 >> s1 >> n1;
//		infile >> s0;
//		infile >> n1;
//		double x,y,z;
//		vector<CPoint> ps;
//		for(int j = 1; j <= n1; ++j){
//			infile >> x >> y >> z;
//			CPoint pnt (x,y,z,1,1);
//			ps.push_back (pnt);
//		}
//
//		points1.push_back (ps);
//	}
//
//	infile >> s0;
//	infile >> s0;
//	infile >> num;
//
//	for(int i = 1; i <= num; ++i){
//		infile >> s0 >> s1 >> n1;
//		infile >> s0;
//		infile >> n1;
//		double x,y,z;
//		vector<CPoint> ps;
//		for(int j = 1; j <= n1; ++j){
//			infile >> x >> y >> z;
//			CPoint pnt (x,y,z,1,1);
//			ps.push_back (pnt);
//		}
//
//		points2.push_back (ps);
//	}
//
//	infile.close ();
//	//¹¹ÔìÔ²ÖùÃæ
//	double zmin1 = 0,zmax1 = 0,zmin2 = 0,zmax2 = 0;
//	for(vector<CPoint>::iterator it = cyPoints1.begin (); it != cyPoints1.end (); ++it){
//		if((zmin1 - (*it).getX ()) > 1.e-10){
//			zmin1 = (*it).getX ();
//		}
//		if(((*it).getX () - zmax1) > 1.e-10){
//			zmax1 = (*it).getX ();
//		}
//	}
//
//	for(vector<CPoint>::iterator it = cyPoints2.begin (); it != cyPoints2.end (); ++it){
//		if((zmin2 - (*it).getX ()) > 1.e-10){
//			zmin2 = (*it).getX ();
//		}
//		if(((*it).getX () - zmax2) > 1.e-10){
//			zmax2 = (*it).getX ();
//		}
//	}
//
//
//	CPoint p1 (0,0,zmin1,1,1);
//	CPoint p2 (0,0,zmin2,1,1);
//	CVector v1 (0,0,zmax1 - zmin1);
//	CVector v2 (0,0,zmax2 - zmin1);
//	surface1.setLocation (p1);
//	surface1.setDirection (v1);
//	surface1.setRadius (cyPoints1[0].getY ());
//	surface1.buildSurface ();
//
//	surface2.setLocation (p2);
//	surface2.setDirection (v2);
//	surface2.setRadius (cyPoints2[0].getY ());
//	surface2.buildSurface ();
//
//	surface3.setInterpolatePoints (points1);
//	surface4.setInterpolatePoints (points2);
//
//	vector<CGeomSurface *> surfaces;
//	surfaces.push_back (&surface1);
//	surfaces.push_back (&surface2);
//	surfaces.push_back (&surface3);
//	surfaces.push_back (&surface4);
//	getIGESSurface (fileName,surfaces);
//}

void CAidFunction::mergeCurves (CGeomCurve &curve1,CGeomCurve &curve2,CBSplineCurve &finalCurve,double tol){
	int type1 = curve1.getCurveType ();
	int type2 = curve2.getCurveType ();
	CBSplineCurve bscurve1,bscurve2;
	adaptorCurve (curve1,bscurve1);
	adaptorCurve (curve2,bscurve2);
	if(type1 == 2 && type2 == 2){
		CBSplineCurve spline1 = dynamic_cast<CBSplineCurve &>(curve1);
		CBSplineCurve spline2 = dynamic_cast<CBSplineCurve &>(curve2);
		Handle_Geom_BoundedCurve c1 = spline1.getGeomCurve ();
		Handle_Geom_BoundedCurve c2 = spline2.getGeomCurve ();
		GeomConvert_CompCurveToBSplineCurve res;
		res.Add (c1,tol);
		res.Add (c2,tol);
		Handle_Geom_BSplineCurve p = res.BSplineCurve ();
		finalCurve.setGeomCurve (p);
	}

}

bool CAidFunction::extractCurve (CGeomSurface &surface,int ORIENTATIONC,double Param,CBSplineCurve &curve){
	CBSplineSurface bssurface;
	adaptorSurface (surface,bssurface);
	vector<CPoint>parampoints;//抽取的曲线采样点
	double umin,umax,vmin,vmax;
	bssurface.getUVMinMax (umin,umax,vmin,vmax);
	double usize = (umax - umin) / 1000;
	double vsize = (vmax - vmin) / 1000;
	if(ORIENTATIONC == 0){


		if(Param >= umin && Param <= umax){
			for(double v = vmin + 1.e-9;v <= vmax;){//在u方向上进行采样
				parampoints.push_back (bssurface.getPointAtParam (Param,v));
				v = v + vsize;
			}
			CBSplineCurve *sc = new CBSplineCurve (parampoints);//插值构造B样条曲线
			curve = *sc;
			return true;

		}
		else
			return false;
	}
	else if(ORIENTATIONC == 1){
		if(Param >= vmin && Param <= vmax){
			for(double u = umin + 1.e-9;u <= umax;){//在v方向上进行采样
				parampoints.push_back (bssurface.getPointAtParam (u,Param));
				u = u + usize;
			}
			CBSplineCurve *sc = new CBSplineCurve (parampoints);//插值构造B样条曲线
			curve = *sc;
			return true;

		}
		else
			return false;
	}
	else
		return false;
}

void CAidFunction::extractmidCurve (CGeomCurve &curve1,CGeomCurve &curve2,CBSplineCurve& midcurve){
	CBSplineCurve bscurve1,bscurve2;
	adaptorCurve (curve1,bscurve1);
	adaptorCurve (curve2,bscurve2);
	double umin1,umax1;
	umax1 = bscurve1.getUmax ();
	umin1 = bscurve1.getUmin ();
	double size1 = (umax1 - umin1) / 100;
	vector<CPoint>points1;
	for(double i = umin1 + 1.e-9;i < umax1;i = i + size1){
		
		points1.push_back (bscurve1.getPointAtParam (i));
	}
	double umin2,umax2;
	umax2 = bscurve2.getUmax ();
	umin2 = bscurve2.getUmin ();
	double size2 = (umax2 - umin2) / 100;
	vector<CPoint>points2;
	for(double i = umin2 + 1.e-9;i < umax2;i = i + size2){

		points2.push_back (bscurve2.getPointAtParam (i));
	}

	vector<CPoint>midpoints;
	size_t sizeu = points1.size () > points2.size () ? points2.size () : points1.size ();
	for(int i = 0;i < sizeu;i++){
		
		double x = (points1[i].getX () + points2[i].getX ()) / 2;
		double y = (points1[i].getY () + points2[i].getY ()) / 2;
		double z = (points1[i].getZ () + points2[i].getZ ()) / 2;
		midpoints.push_back (CPoint (x,y,z));
	}
	CBSplineCurve *mid=new CBSplineCurve(midpoints);
	midcurve = *mid;
	delete mid;
	mid = nullptr;
}

bool CAidFunction::getCCIntersect (Handle_Geom_Curve &curve1,Handle_Geom_Curve &curve2,vector<CPoint> &points,double umin1,double umax1,double umin2,double umax2){
	if (curve1.IsNull() || curve2.IsNull())
		return false;
	GeomAPI_ExtremaCurveCurve intersector (curve1,curve2,umin1,umax1,umin2,umax2);
	int number = intersector.NbExtrema ();
	for(int i = 1; i <= number; ++i){
		double distance;
		distance = intersector.Distance (i);
		gp_Pnt p1,p2;
		intersector.Points (i,p1,p2);
		if(distance < 1e-10){
			points.push_back (CPoint (p1.X (),p1.Y (),p1.Z (),1,1));
		}
	}
	return true;
}

bool CAidFunction::getCSIntersect (Handle_Geom_Curve &curve,Handle_Geom_Surface &surface,vector<CPoint> &points,double u1,double u2,double umin1,double umax1,double umin2,double umax2){
	if (curve.IsNull() || surface.IsNull())
		return false;
	/*GeomAPI_ExtremaCurveSurface intersector(curve, surface, u1, u2, umin1, umax1, umin2, umax2);
	int number = intersector.NbExtrema();
	for (int i = 1; i <= number; ++i)
	{
		double distance;
		distance = intersector.Distance(i);
		gp_Pnt p1, p2;
		intersector.Points(i, p1, p2);
		if (distance < 1e-10)
		{
			points.push_back(CPoint(p1.X(), p1.Y(), p1.Z(), 1, 1));
		}
	}*/
	GeomAPI_IntCS intersector (curve,surface);
	if(intersector.IsDone ()){
		int num = intersector.NbPoints ();
		for(int i = 1; i <= num; ++i){
			gp_Pnt p = intersector.Point (i);
			points.push_back (CPoint (p.X (),p.Y (),p.Z (),1,1));
		}
	}
	return true;
}
//
//void CAidFunction::getCPIntersect (Handle_Geom_Curve &curve,Handle_Geom_Plane &plane,vector<CPoint> &points,double u1,double u2,double umin1,double umax1,double umin2,double umax2){
//	GeomAPI_IntCS intersector (curve,plane);
//	if(intersector.IsDone ()){
//		int num = intersector.NbPoints ();
//		for(int i = 1; i <= num; ++i){
//			gp_Pnt p = intersector.Point (i);
//			points.push_back (CPoint (p.X (),p.Y (),p.Z (),1,1));
//		}
//	}
//}

bool CAidFunction::getSSIntersect (Handle_Geom_Surface &surface1,Handle_Geom_Surface &surface2,vector<CBSplineCurve> &curves,double tol){
	if (surface1.IsNull() || surface2.IsNull())
		return false;
	GeomAPI_IntSS intersector;
	intersector.Perform (surface1,surface2,tol);
	if(intersector.IsDone ()){
		int n = intersector.NbLines ();

		for(int i = 1; i <= n; ++i){
			Handle_Geom_Curve c = intersector.Line (i);
			GeomAdaptor_Curve c1 (c);
			GeomAbs_CurveType type = c1.GetType();
			Handle_Geom_BSplineCurve spline;
			if (type != GeomAbs_BSplineCurve)
			{
				if (type == GeomAbs_Line)
				{
					Handle_Geom_TrimmedCurve tline = new Geom_TrimmedCurve(c, 0, 100);
					spline = GeomConvert::CurveToBSplineCurve(tline);
				}
				else
				{
					spline = GeomConvert::CurveToBSplineCurve(c);
				}
				
			}			
			else
			{
				spline = c1.BSpline();
			}
			CBSplineCurve curve;
			curve.setGeomCurve (spline);
			curves.push_back (curve);
		}
	}
	return true;
}

bool CAidFunction::getPorjectCurve (CGeomSurface &surface,CGeomCurve &curve,CBSplineCurve &projectCurve){
	int type1 = surface.getSurfaceType ();
	int type2 = curve.getCurveType ();
	Handle_Geom_Surface sur;
	vector<CPoint> points;
	vector<CPoint> proPnts;
	if(type1 == 1){
		CPlane plane = dynamic_cast<CPlane &>(surface);
		sur = plane.getGeomPlane ();
	}
	else if(type1 == 2){
		CBSplineSurface spline;
		spline = dynamic_cast<CBSplineSurface &>(surface);
		sur = spline.getGeomSurface ();
	}
	else if(type1 == 3){
		CBezierSurface bs = dynamic_cast<CBezierSurface &>(surface);
		sur = bs.getGeomSurface ();
	}
	if(type2 == 1){
		CLine line;
		line = dynamic_cast<CLine &>(curve);
		double umin = line.getUmin ();
		double umax = line.getUmax ();
		double delta = (double) (umax - umin) / 100;
		for(int i = 0; i <= 100; ++i){
			CPoint p = line.getPointAtParam (umin + i * delta);
			points.push_back (p);
		}
	}
	else if(type2 == 2){
		CBSplineCurve spline;
		spline = dynamic_cast<CBSplineCurve &>(curve);
		double umin = spline.getUmin ();
		double umax = spline.getUmax ();
		double delta = (double) (umax - umin) / 100;
		for(int i = 0; i <= 100; ++i){
			CPoint p = spline.getPointAtParam (umin + i * delta);
			points.push_back (p);
		}
	}
	else if(type2 == 3){
		CBezierCurve bc;
		bc = dynamic_cast<CBezierCurve &>(curve);

		double delta = (double) (1.0 - 0.0) / 100;
		for(int i = 0; i <= 100; ++i){
			CPoint p;
			bc.getPointAtParam (i * delta,p);
			points.push_back (p);
		}
	}
	else if(type2 == 4 || type2 == 5){
		CCircle circle;
		circle = dynamic_cast<CCircle &>(curve);
		double umin = circle.getUmin ();
		double umax = circle.getUmax ();
		double delta = (double) (umax - umin) / 100;
		for(int i = 0; i <= 100; ++i){
			CPoint p = circle.getPointAtParam (umin + i * delta);
			points.push_back (p);
		}
	}

	getPSProject (sur,points,proPnts);

	if(proPnts.size () == 0){
		return false;
	}

	CBSplineCurve spline (proPnts);
	projectCurve = spline;
	return true;
}


void CAidFunction::getPSProject (Handle_Geom_Surface &surface,vector<CPoint> &points,vector<CPoint> &pnts){
	GeomAPI_ProjectPointOnSurf project;
	double u1,u2,v1,v2;
	surface->Bounds (u1,u2,v1,v2);
	project.Init (surface,u1,u2,v1,v2);
	for(auto it : points){
		project.Perform (it.getPnt ());
		if(project.IsDone ()){
			int num = project.NbPoints ();
			if(num != 0){
				CPoint p (project.Point (1));
				pnts.push_back (p);
			}
		}

	}
}

void CAidFunction::mergeSurfaces (CGeomSurface &surface1,CGeomSurface &surface2,CBSplineSurface& mergesurface){
	CBSplineSurface bssurface1,bssurface2;
	adaptorSurface (surface1,bssurface1);
	adaptorSurface (surface2,bssurface2);
	double umin1,umax1,vmin1,vmax1;
	bssurface1.getUVMinMax (umin1,umax1,vmin1,vmax1);
	double sizeu1 = (umax1 - umin1) / 100;
	double sizev1 = (vmax1 - vmin1) / 100;

	double umin2,umax2,vmin2,vmax2;
	bssurface2.getUVMinMax (umin2,umax2,vmin2,vmax2);
	double sizeu2 = (umax2 - umin2) / 100;
	double sizev2 = (vmax2 - vmin2) / 100;
	vector<vector<CPoint>>pointss;
	vector<CPoint>point;
	for(double i = umin1 +1.e-9;i < umax1;i = i + sizeu1){
		for(double j = vmax1 - 1.e-9;j > vmin1;j = j - sizev1){
			point.push_back (bssurface1.getPointAtParam (i,j));
		}
		for(double j = vmin2 + 1.e-9;j < vmax2;j = j + sizev2){
			point.push_back (bssurface2.getPointAtParam (i,j));
		}
		pointss.push_back (point);
		point.clear ();

	}
	CBSplineSurface *merge=new CBSplineSurface(pointss);
	mergesurface = *merge;
	delete merge;
	merge = nullptr;
}

void CAidFunction::adaptorCurve (CGeomCurve &curve,CBSplineCurve&bsplinecurve){
	int type = curve.getCurveType();
	if(type==1){
		CLine line= dynamic_cast<CLine&>(curve);
		Handle_Geom_TrimmedCurve tline = new Geom_TrimmedCurve(line.getGeomCurve(),line.getUmin(),line.getUmax());
		Handle_Geom_BSplineCurve geombscurve = GeomConvert::CurveToBSplineCurve(tline);
		bsplinecurve.setGeomCurve(geombscurve);
	}
	if (type == 2) {
		CBSplineCurve bcurve= dynamic_cast<CBSplineCurve&>( curve );
		bsplinecurve = bcurve;
	}
	if(type == 3){
		CBezierCurve beziercurve = dynamic_cast<CBezierCurve &>(curve);
		Handle_Geom_BSplineCurve geombscurve=GeomConvert::CurveToBSplineCurve (beziercurve.getBezierCurve ());
		bsplinecurve.setGeomCurve(geombscurve);
		/*CBSplineCurve bscurve;
		bscurve.setGeomCurve (geombscurve);
		return bscurve;*/
	}
	if(type == 4||type==5){
		CCircle circle = dynamic_cast<CCircle &>(curve);
		Handle_Geom_BSplineCurve geombscurve = GeomConvert::CurveToBSplineCurve (circle.getGeomCurve());
		bsplinecurve.setGeomCurve(geombscurve);
		/*CBSplineCurve bscurve;
		bscurve.setGeomCurve (geombscurve);
		return bscurve;*/
	}
	
}

 void CAidFunction::adaptorSurface (CGeomSurface&surface,CBSplineSurface& bsplinesurface){
	int type = surface.getSurfaceType();
	if(type==1){
		CPlane plane = dynamic_cast<CPlane&>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface(plane.getRGeomPlane());
		bsplinesurface.setGeomSurface(geombssurface);
	}
	if (type == 2) {
		CBSplineSurface bsurface= dynamic_cast<CBSplineSurface&>( surface );
		bsplinesurface = bsurface;
	}
	if(type == 3){
		CBezierSurface beziersurface = dynamic_cast<CBezierSurface &>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface (beziersurface.getGeomSurface ());
		bsplinesurface.setGeomSurface(geombssurface);
		/*CBSplineSurface bssurface;
		bssurface.setGeomSurface (geombssurface);
		return bssurface;*/
	}
	if(type == 4){
		CCylindricalSurface cylindricalsurface = dynamic_cast<CCylindricalSurface &>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface (cylindricalsurface.getGeomSurface ());
		bool r=geombssurface->IsVRational();
		bsplinesurface.setGeomSurface(geombssurface);
		/*CBSplineSurface bssurface;
		bssurface.setGeomSurface (geombssurface);
		return bssurface;*/
	}
	if(type == 5){
		SurfaceOfLinearExtrusion linesurface = dynamic_cast<SurfaceOfLinearExtrusion &>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface (linesurface.getGeomSurface ());
		bsplinesurface.setGeomSurface(geombssurface);
		/*CBSplineSurface *bssurface = new CBSplineSurface ();
		bssurface->setGeomSurface (geombssurface);
		return *bssurface;*/
	}

	if(type == 6){
		SurfaceOfRevolution revosurface = dynamic_cast<SurfaceOfRevolution &>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface (revosurface.getGeomSurface ());
		bsplinesurface.setGeomSurface(geombssurface);
		/*CBSplineSurface bssurface;
		bssurface.setGeomSurface (geombssurface);
		return bssurface;*/
	}
	if(type == 7){
		CSphere spheresurface = dynamic_cast<CSphere&>(surface);
		Handle_Geom_BSplineSurface geombssurface = GeomConvert::SurfaceToBSplineSurface (spheresurface.getGeomSurface ());
		bsplinesurface.setGeomSurface(geombssurface);
		/*CBSplineSurface bssurface;
		bssurface.setGeomSurface (geombssurface);
		return bssurface;*/
	}
}

void CAidFunction::extractmidSurface (CGeomSurface &surface1,CGeomSurface &surface2,CBSplineSurface& midsurface){
	CBSplineSurface bssurface1,bssurface2;
	adaptorSurface (surface1,bssurface1);
	adaptorSurface (surface2,bssurface2);
	double umin1,umax1,vmin1,vmax1;
	bssurface1.getUVMinMax (umin1,umax1,vmin1,vmax1);
	double sizeu1 = (umax1 - umin1) / 100;
	double sizev1 = (vmax1 - vmin1) / 100;
	vector<CPoint>point1;
	vector<vector<CPoint>>points1;
	for(double i = umin1 + 0.001;i < umax1;i = i + sizeu1){
		for(double j = vmin1 + 0.001;j < vmax1;j = j + sizev1){
			point1.push_back (bssurface1.getPointAtParam (i,j));
		}
		points1.push_back (point1);
		point1.clear ();
	}
	double umin2,umax2,vmin2,vmax2;
	bssurface2.getUVMinMax (umin2,umax2,vmin2,vmax2);
	double sizeu2 = (umax2 - umin2) / 100;
	double sizev2 = (vmax2 - vmin2) / 100;
	vector<CPoint>point2;
	vector<vector<CPoint>>points2;
	for(double i = umin2 + 1.e-9;i < umax2;i = i + sizeu2){
		for(double j = vmin2 + 1.e-9;j < vmax2;j = j + sizev2){
			point2.push_back (bssurface2.getPointAtParam (i,j));
		}
		points2.push_back (point2);
		point2.clear ();
	}
	vector<CPoint>midpoint;
	vector<vector<CPoint>>midpoints;
	size_t sizeu = points1.size () > points2.size () ? points2.size () : points1.size ();
	size_t sizev = points1[0].size () > points2[0].size () ? points2[0].size () : points1[0].size ();
	for(int i = 0;i < sizeu;i++){
		for(int j = 0;j < sizev;j++){
			double x = (points1[i][j].getX () + points2[i][j].getX ()) / 2;
			double y = (points1[i][j].getY () + points2[i][j].getY ()) / 2;
			double z = (points1[i][j].getZ () + points2[i][j].getZ ()) / 2;
			midpoint.push_back (CPoint (x,y,z,1,1));
		}
		midpoints.push_back (midpoint);
		midpoint.clear ();
	}
	CBSplineSurface *mid=new CBSplineSurface(midpoints);
	midsurface = *mid;
	delete mid;
	mid = nullptr;
}

bool CAidFunction::getIGESGeom(int type, string & fileName, vector<CPoint> points, vector<CBSplineCurve>& curves, vector<CBSplineSurface>& surfaces)
{

	IGESControl_Controller::Init();
	IGESControl_Writer igesWriter;
	string tempFile = fileName;
	IGESControl_Reader reader;
	Standard_Integer status = reader.ReadFile(tempFile.c_str());

	if (status != IFSelect_RetDone)
	{
		return false;
	}

	reader.TransferRoots();

	TopoDS_Shape aShape = reader.OneShape();
	
	if (type == 0)//读取点
	{
		TopExp_Explorer explorer(aShape, TopAbs_VERTEX);

		for (; explorer.More();)
		{
			TopoDS_Shape shape = explorer.Current();
			TopoDS_Vertex pnt = TopoDS::Vertex(shape);
			gp_Pnt point = BRep_Tool::Pnt(pnt);
			CPoint *p = new CPoint(point);
			points.push_back(*p);
			explorer.Next();
		}		
	}

	else if (type == 1)//读取线
	{
		TopExp_Explorer explorer(aShape, TopAbs_FACE);
		for (;explorer.More();)
		{
			return false;
		}
		TopExp_Explorer explorer1(aShape, TopAbs_EDGE);

		for (; explorer1.More();)
		{
			TopoDS_Shape shape = explorer1.Current();
			TopoDS_Edge edge = TopoDS::Edge(shape);
			TopLoc_Location location;
			double pFirst=0, pLast=1;
			Handle_Geom_Curve curve = BRep_Tool::Curve(edge,pFirst, pLast);
			GeomAdaptor_Curve adapter(curve);
			Handle_Geom_BSplineCurve pCurve = adapter.BSpline();
			CBSplineCurve * spline = new CBSplineCurve();
			spline->setGeomCurve(pCurve);
			curves.push_back(*spline);
			explorer1.Next();
		}
	}

	else if (type == 2)//读取面
	{
		TopExp_Explorer explorer(aShape, TopAbs_FACE);

		for (; explorer.More();)
		{
			TopoDS_Shape shape = explorer.Current();
			TopoDS_Face face = TopoDS::Face(shape);
			TopLoc_Location location;
			Handle_Geom_Surface surface = BRep_Tool::Surface(face, location);
			GeomAdaptor_Surface adapter(surface);
			Handle_Geom_BSplineSurface pSurface = adapter.BSpline();
			CBSplineSurface * spline = new CBSplineSurface();
			spline->setGeomSurface(pSurface);
			surfaces.push_back(*spline);
			explorer.Next();
		}
	}
	
	return true;
}
