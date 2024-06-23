
#include "CBezierCurve.h"
#include <Geom_BezierCurve.hxx>
#include "DefineType.h"
#include<GeomAdaptor_Curve.hxx>
#include<Extrema_ExtPC.hxx>
void CBezierCurve::TransStruct() {
	m_ControlPoints.clear();
	m_Rational = false;
	m_Closed = true;
	m_Weights.clear();
	for (int i = 1; i <= BezierCurve->Poles().Size(); i++) {
		gp_Pnt pt = BezierCurve->Poles().Value(i);
		m_ControlPoints.push_back(CPoint(pt));
	}
	if (BezierCurve->IsRational()) {
		m_Rational = true;
		for (int j = 1; j <= BezierCurve->Weights()->Size(); j++) {
			m_Weights.push_back(BezierCurve->Weights()->Value(j));
		}
	}

}
Handle_TColgp_HArray1OfPnt CBezierCurve::transStruct(vector<CPoint>P)
{
	gp_Pnt pt;
	Handle_TColgp_HArray1OfPnt points=new TColgp_HArray1OfPnt(1, (int)P.size());
	for (int i = 1; i <= P.size(); i++) {
		P[i - 1].getPnt(pt);
		points->SetValue(i, pt);
	}
	return points;
}
TColStd_HArray1OfReal CBezierCurve::transStruct(vector<double>W)
{
	TColStd_HArray1OfReal Wd(1, (int)W.size());
	for (int i = 1; i <= W.size(); i++) {
		Wd.SetValue(i, W[i - 1]);
	}
	return Wd;
}
void CBezierCurve::setBezierCurve(Handle_Geom_BezierCurve beziercurve)
{
	BezierCurve = beziercurve;
}
CBezierCurve::CBezierCurve()
{
	m_ControlPoints .clear();
	m_Weights.clear();
	m_Closed = true;
	m_Rational = false;
	BezierCurve.Nullify();
}

CBezierCurve::~CBezierCurve()
{

}
CBezierCurve::CBezierCurve(vector<CPoint>& controlPoints)
{
	if (controlPoints.size() < 2 || controlPoints.size() > (BezierCurve->MaxDegree() + 1)) {
		cerr<<"The number of control points is abnormal, and the curve cannot be generated";
		return;
	}
	vector<double>weight;
	weight.clear();
	for (int i = 0; i < controlPoints.size(); i++) {
		weight.push_back(1.0);
	}
	CBezierCurve::createBezierCurve(controlPoints,weight);

}

CBezierCurve::CBezierCurve(vector<CPoint>& controlPoints, vector<double>& weight)
{
	if (controlPoints.size() < 2 || controlPoints.size() > (BezierCurve->MaxDegree() + 1)) {
		cerr<<"The number of control points is abnormal, and the curve cannot be generated";
		return;
	}
	if (weight.size() != controlPoints.size()) {
		cerr << "The number of control points does not match the number of weights";
         return;	
	}
	for (int i = 0; i <controlPoints.size(); i++) {
		if (weight[i]<= gp::Resolution()) {
			cerr<< "The weight is less than the precision";
			return;
		}
	}

	CBezierCurve::createBezierCurve(controlPoints, weight);
}

vector<CBezierCurve> CBezierCurve:: interPolate(vector<CPoint>points, bool periodic, double tolerance)
{
	GeomAPI_Interpolate BS(transStruct(points), periodic, tolerance);
	BS.Perform();
	if(BS.IsDone()){
		GeomConvert_BSplineCurveToBezierCurve BS_BZ(BS.Curve());
		vector<CBezierCurve>bezierarray;
	   for (int i = 1; i <= BS_BZ.NbArcs(); i++) {
		     CBezierCurve bz;
		     bz.setBezierCurve(BS_BZ.Arc(i));
			 bz.TransStruct();
	         bezierarray.push_back(bz);
	   }
	   return bezierarray;
	}
	return vector<CBezierCurve>{};
}

vector<CPoint> CBezierCurve::getControlPoints() const
{
	return m_ControlPoints;
}

vector<double> CBezierCurve::getWeights() const
{
	return m_Weights;
}

int CBezierCurve::getPointNumber()
{
	return (int)m_ControlPoints.size();
}

bool CBezierCurve::getRational() const
{
	return m_Rational;
}

int CBezierCurve::getDegree() const
{
	return (int)((m_ControlPoints).size()-1);
}

bool CBezierCurve::getClosed() const
{
	return m_Closed;
}

Handle_Geom_BezierCurve CBezierCurve::getBezierCurve() const
{
	return BezierCurve;
}

void CBezierCurve::increase(const int Degree)
{
	if (Degree < getDegree()) {
		cerr<<"The order after the promotion cannot be less than the original order";
		return;
	}
	BezierCurve->Increase(Degree);
	this->TransStruct();
}

void CBezierCurve::insertControlPoint(const int Index, CPoint & P, const int type)
{
	if (Index < 0 || Index >= m_ControlPoints.size()) {
		cerr << "The insertion position is out of bounds";
		return;
	}

	gp_Pnt pt;
	P.getPnt(pt);
	if (type == 0) 
		BezierCurve->InsertPoleBefore(Index, pt, 1.0);
	else
	{
		BezierCurve->InsertPoleAfter(Index, pt, 1.0);
	}
	this->TransStruct();
}

void CBezierCurve::insertControlPoint(const int Index,  CPoint & P, const double Weight, const int type)
{
	if (Index < 0 || Index >= m_ControlPoints.size()) {
		cerr << "The insertion position is out of bounds";
		return;
	}
	if (Weight < gp::Resolution()) {
		cerr << "The weight is less than the precision";
		return;
	}
	gp_Pnt pt;
	P.getPnt(pt);
	if (type == 0) {
		BezierCurve->InsertPoleBefore(Index, pt, Weight);
	}
	else
	{
		BezierCurve->InsertPoleAfter(Index, pt, Weight);
	}
	this->TransStruct();
}

void CBezierCurve::removeControlPoint(const int Index)
{
	if (Index < 0 || Index >= m_ControlPoints.size()) {
		cerr << "The remove position is out of bounds";
		return;
	}
	BezierCurve->RemovePole(Index);
	this->TransStruct();
}

void CBezierCurve::setTrimmed(const double U1, const double U2)
{
	BezierCurve->Segment(U1, U2);
	this->TransStruct();
}

void CBezierCurve::setControlPoint(const int Index, CPoint & P)
{
	gp_Pnt pt;
	P.getPnt(pt);
	BezierCurve->SetPole(Index, pt);
	this->TransStruct();
}

void CBezierCurve::setControlPoint(const int Index, CPoint & P, double weight)
{
	gp_Pnt pt;
	P.getPnt(pt);
	BezierCurve->SetPole(Index, pt, weight);
	TransStruct();
}

void CBezierCurve::setWeight(const int Index, const double Weight)
{
	BezierCurve->SetWeight(Index, Weight);
	TransStruct();
}

int CBezierCurve::getCurveType()
{
	return BEZIERCURVE;
}

void CBezierCurve::reverse (){
	BezierCurve->Reverse ();
	this->TransStruct ();
}

bool CBezierCurve::isPeriodic()
{
	return false;
}

void CBezierCurve::getPointAtParam(const double U, CPoint & P) const
{
	gp_Pnt pt;
	if (U < 0.0 || U > 1.0)
	{
		P.setXYZ(INT_MAX, INT_MAX, INT_MAX);
		return;
	}
	BezierCurve->D0(U, pt);
	CPoint Pt=CPoint(pt);
	P = Pt;
}

void CBezierCurve::getParamAtPoint(CPoint P, double & U)
{
	gp_Pnt pnt;
	P.getPnt(pnt);

	GeomAdaptor_Curve apdatorCurve;
	apdatorCurve.Load(BezierCurve);

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

		U= (extPC.Point(index)).Parameter();
		
	}
}

void CBezierCurve::getDerivation1(const double U, CVector & V1) const
{
	gp_Pnt pt;
	gp_Vec vt1;
	if (U < 0.0 || U > 1.0)
	{
		V1.setXYZ(INT_MAX, INT_MAX, INT_MAX);
		return;
	}
	BezierCurve->D1(U, pt, vt1);
	CVector v1(vt1);
	V1.setXYZ(v1.getX(), v1.getY(), v1.getZ());

}

void CBezierCurve::getDerivation2(const double U, CVector & V2) const
{
	gp_Pnt pt;
	gp_Vec vt1, vt2;
	if (U < 0.0 || U > 1.0)
	{
		
		V2.setXYZ(INT_MAX, INT_MAX, INT_MAX);
		return;
	}
	BezierCurve->D2(U, pt, vt1, vt2);
	CVector v2(vt2);
	V2.setXYZ(v2.getX(), v2.getY(), v2.getZ());
}

void CBezierCurve::getDerivation3(const double U, CVector & V3) const
{
	gp_Pnt pt;
	gp_Vec vt1, vt2, vt3;
	if (U < 0.0 || U > 1.0)
	{
		V3.setXYZ(INT_MAX, INT_MAX, INT_MAX);
		return;
	}
	BezierCurve->D3(U, pt, vt1, vt2, vt3);
	CVector v3(vt3);
	V3.setXYZ(v3.getX(), v3.getY(), v3.getZ());
}

void CBezierCurve::setTranslation(CVector vec,CBezierCurve&beziercurve)
{
	gp_Vec v;
	gp_Trsf T;
	v.SetX(vec.getX());
	v.SetY(vec.getY());
	v.SetZ(vec.getZ());
	T.SetTranslation(v);
	BezierCurve->Transform(T);
	/*gp_XYZ coord;
	for (int i = 1; i <= m_ControlPoints.size(); i++) {
		coord.SetX(m_ControlPoints[i-1].getX());
		coord.SetY(m_ControlPoints[i - 1].getY());
		coord.SetY(m_ControlPoints[i - 1].getY());
		T.Transforms(coord);
		m_ControlPoints[i - 1].setX(coord.X());
		m_ControlPoints[i - 1].setY(coord.Y());
		m_ControlPoints[i - 1].setZ(coord.Z());
	}*/
	TransStruct();
	if (!m_Rational)
		beziercurve = CBezierCurve(m_ControlPoints);
	else
		beziercurve = CBezierCurve(m_ControlPoints, m_Weights);
	
}

void CBezierCurve::setRotation(CPoint p, CVector ax, double ang,CBezierCurve&beziercurve)
{
	gp_Pnt pt;
	p.getPnt(pt);
	gp_Dir d(ax.getVec());
	gp_Ax1 ax1(pt,d);
	gp_Trsf T;
	T.SetRotation(ax1, ang);
	BezierCurve->Transform(T);
	/*gp_XYZ coord;
	for (int i = 1; i <= m_ControlPoints.size(); i++) {
		coord.SetX(m_ControlPoints[i - 1].getX());
		coord.SetY(m_ControlPoints[i - 1].getY());
		coord.SetY(m_ControlPoints[i - 1].getY());
		T.Transforms(coord);
		m_ControlPoints[i - 1].setX(coord.X());
		m_ControlPoints[i - 1].setY(coord.Y());
		m_ControlPoints[i - 1].setZ(coord.Z());
	}*/
	TransStruct();
	if (!m_Rational)
		beziercurve = CBezierCurve(m_ControlPoints);
	else
		beziercurve = CBezierCurve(m_ControlPoints, m_Weights);

}

CPoint CBezierCurve::getStartPoint()
{
	gp_Pnt pt;
	pt=BezierCurve->StartPoint();
	return CPoint(pt);
}

CPoint CBezierCurve::getEndPoint()
{
	gp_Pnt pt;
	pt = BezierCurve->EndPoint();
	return CPoint(pt);
}

int CBezierCurve::getNumPoint()
{
	return (int)m_ControlPoints.size();
}

double CBezierCurve::getCurveLength()
{
	double u_first = BezierCurve->FirstParameter();
	double u_last = BezierCurve->LastParameter();
	GeomAdaptor_Curve curve(BezierCurve);
	double length=GCPnts_AbscissaPoint::Length(curve,u_first,u_last);
	return length;
}


void CBezierCurve::createBezierCurve(vector<CPoint>&controlPoints, vector<double>&weight)
{
	
    m_ControlPoints =controlPoints;
	m_Weights = weight;
	size_t size_points = controlPoints.size();
	size_t size_weights = weight.size();
	gp_Pnt pt;
	TColgp_Array1OfPnt points(1, (int)size_points);
	TColStd_Array1OfReal weights(1, (int)size_weights);
	for (int i = 1; i <= size_points; i++) {
		controlPoints[i - 1].getPnt(pt);
		points.SetValue(i, pt);
	}
	if (weight.empty()) {
		m_Rational = false;
		BezierCurve = new Geom_BezierCurve(points);
	}
	else {
		m_Rational = true;
		for (int i = 1; i <= size_weights; i++) {
			weights.SetValue(i, weight[i - 1]);
		}
		BezierCurve = new Geom_BezierCurve(points, weights);
	}
	TransStruct();
}


