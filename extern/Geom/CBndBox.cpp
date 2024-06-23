#include "CBndBox.h"

CBndBox::CBndBox ():Xmin (0.),Xmax (0.),Ymin (0.),Ymax (0.),Zmin (0.),Zmax (0.) {
	mBox.SetVoid ();
}

void CBndBox::updateBox (const double aXmin,const double aYmin,const double aZmin,const double aXmax,const double aYmax,const double aZmax) {
	mBox.Update (aXmin,aYmin,aZmin,aXmax,aYmax,aZmax);
	updateBoxParam ();
}

void CBndBox::getBound (double& theXmin,double& theYmin,double& theZmin,double& theXmax,double& theYmax,double& theZmax) const {
	mBox.Get (theXmin,theYmin,theZmin,theXmax,theYmax,theZmax);
}

CPoint CBndBox::getCornerMin () const {
	gp_Pnt pnt;
	pnt=mBox.CornerMin ();
	CPoint point (pnt);
	return point;
}

CPoint CBndBox::getCornerMax () const {
	gp_Pnt pnt;
	pnt = mBox.CornerMax ();
	CPoint point (pnt);
	return point;
}

bool CBndBox::isOut (const CPoint& point) const {
	gp_Pnt pnt;
	point.getPnt (pnt);
	bool out=mBox.IsOut (pnt);
	return out;
}

bool CBndBox::isOut (CLine line) const {
	CPoint point;
	CVector vector;
	line.getPoint (point);
	line.getDirection (vector);
	gp_Pnt pnt = point.getPnt ();
	gp_Dir dir;
	dir.SetX (vector.getX ());
	dir.SetY (vector.getY ());
	dir.SetZ (vector.getZ ());
	gp_Lin lin (pnt,dir);
	bool out=mBox.IsOut (lin);
	return out;
}

bool CBndBox::isOut (CPlane plane) const {
	CPoint point;
	CVector vector;
	plane.getPoint (point);
	plane.getNormal (vector);
	gp_Pnt pnt = point.getPnt ();
	gp_Dir dir;
	dir.SetX (vector.getX ());
	dir.SetY (vector.getY ());
	dir.SetZ (vector.getZ ());
	gp_Pln pln(pnt,dir);
	bool out = mBox.IsOut (pln);
	return out;
}

bool CBndBox::isOut (const CBndBox& other) const {
	bool out=mBox.IsOut (other.mBox);
	return out;
}

double CBndBox::getDistance (const CBndBox& other) const {
	double dis = mBox.Distance (other.mBox);
	return dis;
}

double CBndBox::getDiagonal () const {
	double dia = sqrt (mBox.SquareExtent ());
	return dia;
}

Bnd_Box CBndBox::getmBndBox () {
	return mBox;
}

void CBndBox::setmBndBox (Bnd_Box Box) {
	mBox = Box;
	updateBoxParam ();
}

void CBndBox::updateBoxParam () {
	mBox.Get (Xmin,Xmax,Ymin,Ymax,Zmin,Zmax);
}
