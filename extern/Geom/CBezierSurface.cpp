
#include "CBezierSurface.h"


CBezierSurface::CBezierSurface () {
}


CBezierSurface::~CBezierSurface () {
}

CBezierSurface::CBezierSurface (vector<vector<CPoint>>& controlPoints) {
	vector<vector<double>>weights;
	weights.clear ();
	/*for (int i = 0; i < controlPoints.size(); i++) {
		vector<double>w;
		for (int j = 0; j <controlPoints[i].size(); j++) {

			w.push_back(1.0);

		}
		weights.push_back(w);
		w.clear();
	}*/
	CBezierSurface::createBezierSurface (controlPoints,weights);
}

CBezierSurface::CBezierSurface (vector<vector<CPoint>>& controlPoints,vector<vector<double>>& weights) {
	CBezierSurface::createBezierSurface (controlPoints,weights);
}

vector<vector<CBezierSurface>> CBezierSurface::interPolate (vector<vector<CPoint>>& points) {

	GeomAPI_PointsToBSplineSurface pBS;
	pBS.Interpolate (transStruct (points));
	GeomConvert_BSplineSurfaceToBezierSurface  BsToBzier (pBS.Surface ());
	vector<vector<CBezierSurface>> Bziersurface_array;
	for (int i = 1; i <= BsToBzier.NbUPatches (); i++) {
		vector<CBezierSurface>bsf_array;
		for (int j = 1; j <= BsToBzier.NbVPatches (); j++) {
			CBezierSurface bsf;
			bsf.setBezierSurface (BsToBzier.Patch (i,j));
			bsf.TransStruct ();
			bsf_array.push_back (bsf);
		}
		Bziersurface_array.push_back (bsf_array);
	}
	return Bziersurface_array;
}

vector<vector<CPoint>> CBezierSurface::getControlPoints () const {
	return m_ControlPoints;
}

vector<vector<double>> CBezierSurface::getWeights () const {
	return m_Weights;
}

int CBezierSurface::getSizeUControlPoints () const {
	auto it = m_ControlPoints.begin ();
	return (int) m_ControlPoints.size ();
}

int CBezierSurface::getSizeVControlPoints () const {
	auto it = m_ControlPoints.begin ();
	return (int) ( ( *it ).size () );
}

CPoint CBezierSurface::getPoint (const int Uindex,const int Vindex) const {
	return m_ControlPoints [Uindex][Vindex];
}

int CBezierSurface::getUDegree () const {
	return (int) m_ControlPoints.size () - 1;
}

int CBezierSurface::getVDegree () const {
	auto it = m_ControlPoints.begin ();
	return (int) ( ( *it ).size () - 1 );
}

Handle_Geom_BezierSurface CBezierSurface::getGeomSurface () {
	return BezierSurface;
}

//CBezierCurve CBezierSurface::getUIso(const double U)
//{
//	BezierSurface->UIso(U);

//}
void CBezierSurface::TransStruct () {
	m_Urational = BezierSurface->IsURational ();
	m_Vrational = BezierSurface->IsVRational ();
	m_ControlPoints.clear ();
	for (int i = 1; i <= BezierSurface->Poles ().RowLength (); i++) {
		vector<CPoint>pole;
		for (int j = 1; j <= BezierSurface->Poles ().ColLength (); j++) {
			pole.push_back (CPoint (BezierSurface->Poles ().Value (j,i)));
		}
		m_ControlPoints.push_back (pole);
	}
	if (BezierSurface->IsURational () && BezierSurface->IsVRational ()) {
		m_Urational = true;
		m_Vrational = true;
		m_Weights.clear ();
		for (int i = 1; i <= BezierSurface->Weights ()->RowLength (); i++) {
			vector<double>we;
			for (int j = 1; j <= BezierSurface->Weights ()->ColLength (); j++) {
				we.push_back (BezierSurface->Weights ()->Value (j,i));
			}
			m_Weights.push_back (we);
		}
	}
}
TColgp_Array1OfPnt CBezierSurface::transStruct (vector<CPoint>P) {
	gp_Pnt pt;
	TColgp_Array1OfPnt points (1,(int) P.size ());
	for (int i = 1; i <= P.size (); i++) {
		P [i - 1].getPnt (pt);
		points.SetValue (i,pt);
	}
	return points;
}
TColStd_Array1OfReal CBezierSurface::transStruct (vector<double>W) {
	TColStd_Array1OfReal Wd (1,(int) W.size ());
	for (int i = 1; i <= W.size (); i++) {
		Wd.SetValue (i,W [i - 1]);
	}
	return Wd;
}
TColgp_Array2OfPnt CBezierSurface::transStruct (vector<vector<CPoint>>P) {
	vector<vector<CPoint>>::iterator it_v = P.begin ();
	size_t size_upoints = P.size ();
	size_t size_vpoints = ( *it_v ).size ();
	gp_Pnt pt;
	int i = 0;
	TColgp_Array2OfPnt points (1,(int) size_upoints,1,(int) size_vpoints);

	for (it_v; it_v != P.end (); it_v++) {
		i++;
		int j = 0;
		for (auto it_vv = ( *it_v ).begin (); it_vv != ( *it_v ).end (); it_vv++) {
			j++;
			( *it_vv ).getPnt (pt);
			points.SetValue (j,i,pt);

		}

	}
	return  points;
}
void CBezierSurface::setBezierSurface (Handle_Geom_BezierSurface BS) {
	BezierSurface = BS;
}
void CBezierSurface::createBezierSurface (vector<vector<CPoint>>& controlPoints,vector<vector<double>>& weights) {

	m_ControlPoints = controlPoints;
	m_Weights = weights;
	vector<vector<CPoint>>::iterator it_v = controlPoints.begin ();
	vector<vector<double>>::iterator itt_v = weights.begin ();
	size_t size_upoints = controlPoints.size ();
	size_t size_vpoints = ( controlPoints [0] ).size ();
	size_t size_uweights = 0;
	size_t size_vweights = 0;
	if (!weights.empty ()) {
		size_t size_uweights = weights.size ();
		size_t size_vweights = ( weights [0] ).size ();
	}
	gp_Pnt pt;
	int i = 0;
	TColgp_Array2OfPnt points (1,(int) size_upoints,1,(int) size_vpoints);
	//ת
	for (int i = 1; i <= size_upoints; ++i) {
		for (int j = 1; j <= size_vpoints; ++j) {
			gp_Pnt pnt;
			( m_ControlPoints ) [i - 1][j - 1].getPnt (pnt);
			points.SetValue (i,j,pnt);
		}
	}
	/*for (it_v; it_v != controlPoints.end(); it_v++) {
		i++;
		int j = 0;
		for (auto it_vv = (*it_v).begin(); it_vv != (*it_v).end(); it_vv++) {
			j++;
			(*it_vv).getPnt(pt);
			points.SetValue(j, i, pt);

		}

	}*/

	if (weights.size () <= 1) {
		if (weights.empty ())
			m_Urational = false;
		if (weights.size () == 1)
			m_Vrational = false;
		BezierSurface = new Geom_BezierSurface (points);
	}
	else {
		m_Urational = m_Vrational = true;
		TColStd_Array2OfReal weight (1,(int) size_uweights,1,(int) size_vweights);
		//ת
		for (int i = 1; i <= size_uweights; ++i) {
			for (int j = 1; j <= size_vweights; ++j) {
				double w;
				w = ( m_Weights ) [i - 1][j - 1];
				weight.SetValue (i,j,w);
			}
		}

		BezierSurface = new Geom_BezierSurface (points,weight);
	}
	TransStruct ();
}

double CBezierSurface::getSingleWeight (const int Uindex,const int Vindex) const {
	return m_Weights [Uindex][Vindex];
}

bool CBezierSurface::getURational () const {
	return m_Urational;
}

bool CBezierSurface::getVRational () const {
	return m_Vrational;
}

int CBezierSurface::getSurfaceType () {
	return BEZIERSURFACE;
}

void CBezierSurface::exchangeUV () {
	BezierSurface->ExchangeUV ();
	TransStruct ();
}

void CBezierSurface::increaseDgree (const int Udegree,const int Vdegree) {
	BezierSurface->Increase (Udegree,Vdegree);
	TransStruct ();
}

void CBezierSurface::insertVControlPoint (const int Vindex,vector<CPoint>point,int type) {
	gp_Pnt pt;
	TColgp_Array1OfPnt points (1,(int) point.size ());
	for (int i = 1; i <= point.size (); i++) {
		point [i - 1].getPnt (pt);
		points.SetValue (i,pt);

	}
	if (type == 0) {
		BezierSurface->InsertPoleColBefore (Vindex,points);
	}
	else
		BezierSurface->InsertPoleColAfter (Vindex,points);
	CBezierSurface::TransStruct ();
}

void CBezierSurface::insertVControlPoint (const int Vindex,vector<CPoint>& P,const vector<double> Weight,const int type) {
	gp_Pnt pt;
	TColgp_Array1OfPnt points (1,(int) P.size ());
	TColStd_Array1OfReal weights (1,(int) Weight.size ());
	for (int i = 1; i <= P.size (); i++) {
		P [i - 1].getPnt (pt);
		points.SetValue (i,pt);

	}
	for (int j = 1; j <= Weight.size (); j++) {
		weights.SetValue (j,Weight [j - 1]);
	}
	if (type == 0) {
		BezierSurface->InsertPoleColBefore (Vindex,points,weights);
	}
	else
		BezierSurface->InsertPoleColAfter (Vindex,points,weights);
	TransStruct ();
}

void CBezierSurface::insertUControlPoint (const int Uindex,vector<CPoint> P,int type) {
	gp_Pnt pt;
	TColgp_Array1OfPnt points (1,(int) P.size ());
	//TColStd_Array1OfReal weights(1, Weight.size());
	for (int i = 1; i <= P.size (); i++) {
		P [i - 1].getPnt (pt);
		points.SetValue (i,pt);

	}
	//for (int j = 1; j <= Weight.size(); j++) {
	//	weights.SetValue(j, Weight[j - 1]);
	//}
	if (type == 0) {
		BezierSurface->InsertPoleRowBefore (Uindex,points);
	}
	else
		BezierSurface->InsertPoleRowAfter (Uindex,points);
	TransStruct ();
}

void CBezierSurface::insertUControlPoint (const int Uindex,vector<CPoint> P,const vector<double> Weight,const int type) {
	gp_Pnt pt;
	TColgp_Array1OfPnt points (1,(int) P.size ());
	TColStd_Array1OfReal weights (1,(int) Weight.size ());
	for (int i = 1; i <= P.size (); i++) {
		P [i - 1].getPnt (pt);
		points.SetValue (i,pt);

	}
	for (int j = 1; j <= Weight.size (); j++) {
		weights.SetValue (j,Weight [j - 1]);
	}
	if (type == 0) {
		BezierSurface->InsertPoleRowBefore (Uindex,points,weights);
	}
	else
		BezierSurface->InsertPoleRowAfter (Uindex,points,weights);
	TransStruct ();
}

void CBezierSurface::removeVControlPoint (const int Vindex) {
	BezierSurface->RemovePoleCol (Vindex);
	TransStruct ();
}

void CBezierSurface::removeUControlPoint (const int Uindex) {
	BezierSurface->RemovePoleRow (Uindex);
	TransStruct ();
}

void CBezierSurface::setTrimmed (const double U1,const double U2,const double V1,const double V2) {
	BezierSurface->Segment (U1,U2,V1,V2);
	TransStruct ();
}

void CBezierSurface::setSingleControlPoint (const int Uindex,const int Vindex,const CPoint& point) {
	m_ControlPoints [Uindex][Vindex] = point;
	TransStruct ();
}

void CBezierSurface::setSingleControlPoint (const int Uindex,const int Vindex,const CPoint& point,const double weight) {
	m_ControlPoints [Uindex][Vindex] = point;
	m_Weights [Uindex][Vindex] = weight;
	TransStruct ();
}

void CBezierSurface::setVControlPoint (const int Vindex,const vector<CPoint> P) {
	BezierSurface->SetPoleCol (Vindex,transStruct (P));
	TransStruct ();
}

void CBezierSurface::setVControlPoint (const int Vindex,const vector<CPoint> P,const vector<double> weights) {
	BezierSurface->SetPoleCol (Vindex,transStruct (P),transStruct (weights));
	TransStruct ();
}

void CBezierSurface::setUControlPoint (const int Uindex,const vector<CPoint> P) {
	BezierSurface->SetPoleRow (Uindex,transStruct (P));
	TransStruct ();
}

void CBezierSurface::setUControlPoint (const int Uindex,const vector<CPoint> P,const vector<double> weights) {
	BezierSurface->SetPoleRow (Uindex,transStruct (P),transStruct (weights));
	TransStruct ();
}

void CBezierSurface::setSingleWeight (const int Uindex,const int Vindex,const double weight) {
	m_Weights [Uindex][Vindex] = weight;
}

void CBezierSurface::setVWeight (const int Vindex,vector<double> weights) {
	BezierSurface->SetWeightCol (Vindex,transStruct (weights));
	TransStruct ();
}

void CBezierSurface::setUWeight (const int Uindex,vector<double> weights) {
	BezierSurface->SetWeightRow (Uindex,transStruct (weights));
	TransStruct ();
}

void CBezierSurface::getPointAtParam (const double U,const double V,CPoint& P) const {
	gp_Pnt pt;
	if (U < 0 || U > 1 || V < 0 || V > 1) {
		P.setXYZ (INT_MAX,INT_MAX,INT_MAX);
		return;
	}

	BezierSurface->D0 (U,V,pt);
	P.setX (pt.X ());
	P.setX (pt.Y ());
	P.setX (pt.Y ());
}

void CBezierSurface::getDerivation1 (const double U,const double V,vector<CVector>& vectors) const {
	gp_Pnt pt;
	gp_Vec vt1,ut1;
	if (U < 0 || U > 1 || V < 0 || V > 1) {
		CVector tmpV (INT_MAX,INT_MAX,INT_MAX,-1);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		return;
	}
	BezierSurface->D1 (U,V,pt,ut1,vt1);
	CVector v1 (vt1);
	CVector u1 (ut1);
	vectors.push_back (u1);
	vectors.push_back (v1);
}

void CBezierSurface::getDerivation2 (const double U,const double V,vector<CVector>& vectors) const {
	gp_Pnt pt;
	gp_Vec vt1,ut1,ut2,vt2,uvt2;
	if (U < 0 || U > 1 || V < 0 || V > 1) {
		CVector tmpV (INT_MAX,INT_MAX,INT_MAX,-1);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		return;
	}
	BezierSurface->D2 (U,V,pt,ut1,vt1,ut2,vt2,uvt2);
	CVector u2 (ut2);
	CVector v2 (vt2);
	CVector uv2 (uvt2);
	vectors.push_back (u2);
	vectors.push_back (v2);
	vectors.push_back (uv2);
}

void CBezierSurface::getDerivation3 (const double U,const double V,vector<CVector>& vectors) const {
	gp_Pnt pt;
	gp_Vec vt1,ut1,ut2,vt2,uvt2,ut3,vt3,uuvt3,uvvt3;
	if (U < 0 || U > 1 || V < 0 || V > 1) {
		CVector tmpV (INT_MAX,INT_MAX,INT_MAX,-1);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		vectors.push_back (tmpV);
		return;
	}
	BezierSurface->D3 (U,V,pt,ut1,vt1,ut2,vt2,uvt2,ut3,vt3,uuvt3,uvvt3);

	CVector u3 (ut3);
	CVector v3 (vt3);
	CVector uuv3 (uuvt3);
	CVector uvv3 (uvvt3);
	vectors.push_back (u3);
	vectors.push_back (v3);
	vectors.push_back (uuv3);
	vectors.push_back (uvv3);
}

void CBezierSurface::setTranslation (CVector vec) {
	gp_Vec v;
	gp_Trsf T;
	v.SetX (vec.getX ());
	v.SetY (vec.getY ());
	v.SetZ (vec.getZ ());
	T.SetTranslation (v);
	BezierSurface->Transform (T);
	TransStruct ();
	createBezierSurface (m_ControlPoints,m_Weights);
}

void CBezierSurface::setRotation (CPoint p,CVector ax,double ang) {
	gp_Pnt pt;
	p.getPnt (pt);
	gp_Dir d (ax.getVec ());
	gp_Ax1 ax1 (pt,d);
	gp_Trsf T;
	T.SetRotation (ax1,ang);
	BezierSurface->Transform (T);
	TransStruct ();
	createBezierSurface (m_ControlPoints,m_Weights);
}

