#include "CBSplineCurve.h"
#include "CLine.h"
#include "DefineType.h"
#include <GeomAdaptor_Curve.hxx>
#include <GCPnts_AbscissaPoint.hxx>
#include <GeomAPI_ProjectPointOnCurve.hxx>
#include <GeomAPI.hxx>
#include <gp_Pln.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <GeomConvert_CompCurveToBSplineCurve.hxx>

CBSplineCurve::CBSplineCurve()
{
	m_iDegree = 0;
	m_iTypeCoord = 1;
	m_dLength = 0;
	
	m_ControlPoints.clear();
	m_Knots.clear();
	m_Multis.clear();
	m_Weights.clear();
	m_bPeriodic = false;
	m_bRational = false;
	//m_pCurve = nullptr;
}

CBSplineCurve::~CBSplineCurve()
{
	if (!m_pCurve.IsNull())
	{
		//m_pCurve->Delete();
		m_pCurve.Nullify();
	}

	m_ControlPoints.clear();
	m_Knots.clear();
	m_Multis.clear();
	m_Weights.clear();

}

CBSplineCurve::CBSplineCurve(const CBSplineCurve & splineCurve)
{
	m_iDegree = splineCurve.m_iDegree;
	m_iTypeCoord = 1;
	m_ControlPoints = splineCurve.m_ControlPoints;
	m_Knots = splineCurve.m_Knots;
	m_Multis = splineCurve.m_Multis;
	m_bRational = splineCurve.m_bRational;
	m_bPeriodic = splineCurve.m_bPeriodic;
	m_Weights = splineCurve.m_Weights;

	createBSplineCurve();
	if (!m_pCurve.IsNull())
	{
		computeLength();
	}
}

CBSplineCurve::CBSplineCurve(const int degree, const vector<CPoint>& controlPoints, const vector<double>& knots, const vector<int>& multis, bool periodic)
{
	m_iDegree = degree;
	m_iTypeCoord = 1;
	m_ControlPoints = controlPoints;
	m_Knots = knots;
	m_Multis = multis;
	m_bRational = false;
	m_bPeriodic = periodic;

	createBSplineCurve();
	if (!m_pCurve.IsNull())
	{
		computeLength();
	}
}

CBSplineCurve::CBSplineCurve(int degree, const vector<CPoint>& controlPoints, const vector<double>& knots, const vector<int>& multis, const vector<double>& weights/*, bool rational*/, bool periodic/*, int typeCoord*/)
{
	m_iDegree = degree;
	m_iTypeCoord = 1;
	m_ControlPoints = controlPoints;
	m_Knots = knots;
	m_Multis = multis;
	m_Weights = weights;
	m_bRational = true;
	m_bPeriodic = periodic;

	createBSplineCurve();
	if (!m_pCurve.IsNull())
	{
		computeLength();
	}
}

CBSplineCurve::CBSplineCurve(const vector<CPoint>& curvePoints, double tol)
{
	size_t size = curvePoints.size();
	if (size <= 1)
	{
		return;
	}

	Handle_TColgp_HArray1OfPnt interpolatePoints = new TColgp_HArray1OfPnt(1, (int)size);
	for (int i = 1; i <= size; ++i)
	{
		gp_Pnt pnt;
		(curvePoints)[i - 1].getPnt(pnt);
		interpolatePoints->SetValue(i, pnt);
	}

	GeomAPI_Interpolate aInterpolater(interpolatePoints, Standard_False, tol);
	aInterpolater.Perform();
	if (aInterpolater.IsDone())
	{
		m_pCurve = aInterpolater.Curve();
		computeLength();
		updateCurveDate();
	}
	else
	{
		m_pCurve = nullptr;
	}

	interpolatePoints.Nullify();
}

CBSplineCurve::CBSplineCurve(const vector<CPoint>& curvePoints, const vector<double>& weights, double tol)
{
	size_t size = curvePoints.size();
	size_t weightSize = weights.size();
	if (size <= 1)
	{
		return;
	}

	Handle_TColgp_HArray1OfPnt interpolatePoints = new TColgp_HArray1OfPnt(1, (int)size);
	Handle_TColStd_HArray1OfReal parameters = new TColStd_HArray1OfReal(1, (int)weightSize);

	for (int i = 1; i <= size; ++i)
	{
		gp_Pnt pnt;
		(curvePoints)[i - 1].getPnt(pnt);
		interpolatePoints->SetValue(i, pnt);
	}

	for (size_t i = 1; i <= weightSize; ++i)
	{
		parameters->SetValue((int)i, weights[i- 1]);
	}

	GeomAPI_Interpolate aInterpolater(interpolatePoints, parameters, Standard_False, tol);
	aInterpolater.Perform();
	if (aInterpolater.IsDone())
	{
		m_pCurve = aInterpolater.Curve();
		computeLength();
		updateCurveDate();
	}
	else
	{
		m_pCurve = nullptr;
	}

	interpolatePoints.Nullify();
	parameters.Nullify();
}

CBSplineCurve::CBSplineCurve(const vector<CBSplineCurve> &curves, double tol)
{
	GeomConvert_CompCurveToBSplineCurve compCurve;
	size_t size = curves.size();
	for (size_t i = 0; i < size; ++i)
	{
		compCurve.Add(curves[i].getGeomCurve(), tol);
	}
	m_pCurve = compCurve.BSplineCurve();
	computeLength();
	updateCurveDate();

	//释放合并的曲线  假设由底层解决
}

int CBSplineCurve::getCurveType()
{
	return BSPLINECURVE;
}

bool CBSplineCurve::setDegree(int degree) 
{	
	if (m_pCurve.IsNull())
	{
		m_iDegree = degree;
		return true;
	}
	else
	{
		return false;
		//暂时屏蔽升阶
		if (degree > m_iDegree)  //升阶
		{
			m_pCurve->IncreaseDegree(degree);
			//m_iDegree = degree;
			computeLength();
			//获取数据
			updateCurveDate();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineCurve::setControlPoints(const vector<CPoint>& controlPoints)
{
	if (m_pCurve.IsNull())
	{
		m_ControlPoints = controlPoints;
		return true;
	}
	else
	{
		size_t size = m_ControlPoints.size();
		if (size == controlPoints.size())
		{
			m_ControlPoints = controlPoints;
			for (int i = 1; i <= size; ++i)
			{
				gp_Pnt p = controlPoints[i - 1].getPnt();
				m_pCurve->SetPole(i, p);
			}
			computeLength();
			//bool state = createBSplineCurve();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineCurve::setControlPoints(const vector<CPoint>& controlPoints, const vector<double>& weights)
{
	size_t sizew = weights.size();
	for (size_t i = 0; i < sizew; ++i)
	{
		if (weights[i] < 1e-10)
		{
			return false;
		}
	}
	if (m_pCurve.IsNull())
	{
		m_ControlPoints = controlPoints;
		m_Weights = weights;
		m_bRational = true;
		return true;
	}
	else
	{
		size_t size = m_ControlPoints.size();
		if (size == controlPoints.size())
		{
			m_ControlPoints = controlPoints;
			m_Weights = weights;
			m_bRational = true;
			for (int i = 1; i <= size; ++i)
			{
				gp_Pnt p = controlPoints[i - 1].getPnt();
				m_pCurve->SetPole(i, p, weights[i - 1]);
			}
			computeLength();
			//bool state = createBSplineCurve();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineCurve::setSingleControlPoint(int index, const CPoint & controlPoint)
{
	size_t size = m_ControlPoints.size();
	if (m_pCurve.IsNull() && (index < 1 || index > size) )
	{
		return false;
	}
	m_ControlPoints[index - 1] = controlPoint;
	gp_Pnt p;
	controlPoint.getPnt(p);
	m_pCurve->SetPole(index, p);
	//长度更新
	computeLength();
	return true;
}

bool CBSplineCurve::setSingleControlPoint(int index, const CPoint & controlPoint, double weight)
{
	size_t size = m_ControlPoints.size();
	if (m_pCurve.IsNull() && (index < 1 || index > size) && !m_bRational)
	{
		return false;
	}
	if (weight < 1e-10)
	{
		return false;
	}
	m_ControlPoints[index - 1] = controlPoint;
	m_Weights[index - 1] = weight;
	gp_Pnt p;
	controlPoint.getPnt(p);
	m_pCurve->SetPole(index, p, weight);
	computeLength();
	return true;
}

bool CBSplineCurve::setWeights(const vector<double>& weights)
{
	size_t sizew = weights.size();
	if (sizew == 0)
	{
		return false;
	}
	for (size_t i = 0; i < sizew; ++i)
	{
		if (weights[i] < 1e-10)
		{
			return false;
		}
	}
	if (m_pCurve.IsNull())
	{
		m_Weights = weights;
		m_bRational = true;
		return true;
	}
	else
	{
		if (m_bRational == false)
		{
			m_Weights = weights;
			m_bRational = true;
		}
		size_t size = m_Weights.size();
		if (size == sizew)
		{
			m_Weights = weights;
			for (int i = 1; i <= size; ++i)
			{
				m_pCurve->SetWeight(i, weights[i]);
			}
			computeLength();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineCurve::setSingleWeight(int index, double weight)
{
	size_t size = m_Weights.size();
	if (m_pCurve.IsNull() && (index < 1 || index > size) && !m_bRational)
	{
		return false;
	}
	
	if (weight < 1e-10)
	{
		return false;
	}

	m_Weights[index - 1] = weight;
	m_pCurve->SetWeight(index, weight);
	computeLength();
	return true;
}

bool CBSplineCurve::setKnots(const vector<double>& knots)
{
	//递增序列检查
	if (m_pCurve.IsNull())
	{
		m_Knots = knots;
		return true;
	}
	size_t size = m_Knots.size();
	if (size == knots.size())
	{
		m_Knots = knots;
		for (int i = 1; i <= size; ++i)
		{
			m_pCurve->SetKnot(i, knots[i - 1]);
		}
		computeLength();
		return true;
	}
	return false;

}

bool CBSplineCurve::setSingleKnot(int index, double knot)
{
	size_t size = m_Knots.size();
	if (!m_pCurve.IsNull() && (index < 1 || index > size))
	{
		return false;
	}
	//递增序列检查
	m_Knots[index - 1] = knot;
	m_pCurve->SetKnot(index, knot);
	computeLength();
	return true;
}

bool CBSplineCurve::setMultis(const vector<int>& multis)
{
	if (m_pCurve.IsNull())
	{
		m_Multis = multis;
		return true;
	}
	size_t size = m_Multis.size();
	if (size == multis.size())
	{
		bool flag = false;
		for (int i = 1; i <= size; ++i)
		{
			if (m_Multis[i - 1] < multis[i - 1])
			{
				m_pCurve->IncreaseMultiplicity(i, multis[i - 1]);  //类似degree  只能增加至degree  不能减少
				m_Multis[i - 1] = multis[i - 1];
				flag = true;
			}
			
		}
		if (flag)
		{
			//数据更新
			updateCurveDate();
			computeLength();
		}
		
		return true;
	}
	return false;
}

bool CBSplineCurve::setSingleMultis(int index, int multi)
{
	size_t size = m_Knots.size();
	if (m_pCurve.IsNull() && (index < 1 || index > size))
	{
		return false;
	}
	if (m_Multis[index - 1] >= multi)
	{
		return false;
	}
	m_Multis[index - 1] = multi;
	m_pCurve->IncreaseMultiplicity(index, multi);
	//数据更新
	updateCurveDate();
	computeLength();

	return true;
}

bool CBSplineCurve::setPeriodic(bool periodic)
{
	if (m_bPeriodic != periodic)
	{
		if (periodic)
		{
			m_pCurve->SetPeriodic();
		}
		else
		{
			m_pCurve->SetNotPeriodic();
		}
		//数据更新，待补充
		updateCurveDate();
		computeLength();
		m_bPeriodic = periodic;
		return true;
	}
	return false;
}

bool CBSplineCurve::insertKnot(double u, int multi)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	/*vector<double> cd;
	vector<int> ci;
	int size = m_Knots.size();
	cd.resize(size + 1);
	ci.resize(size + 1);
	for (int i = 0; i <size; ++i)
	{
		double d = m_Knots[i];
		if ((d - u) >= 1e-10)
		{
			cd.push_back(u);
			ci.push_back(multi);
		}
		cd.push_back(d);
		ci.push_back(m_Multis[i]);
	}

	m_Multis = ci;
	m_Knots = cd;*/
	double umin = m_pCurve->FirstParameter();
	double umax = m_pCurve->LastParameter();
	if ((u - umin) >= 0.0 && (umax - u) >= 0.0)
	{
		m_pCurve->InsertKnot(u, multi);
		//数据更新
		updateCurveDate();
		computeLength();
		return true;
	}
	
	return false;
}

bool CBSplineCurve::removeKnot(int index, int multi)
{
	size_t size = m_Knots.size();
	if (m_pCurve.IsNull() || index < 1 || index > size)
	{
		return false;
	}
	/*vector<double> cd;
	vector<int> ci;
	double u = m_pCurve->Knot(index);
	cd.resize(size - 1);
	ci.resize(size - 1);
	for (int i = 0; i < size; ++i)
	{
		double d = m_Knots[i];
		if (fabs(d - u) < 1e-10)
		{
			continue;
		}
		cd.push_back(d);
		ci.push_back(m_Multis[i]);
	}

	m_Multis = ci;
	m_Knots = cd;*/
	bool state = m_pCurve->RemoveKnot(index, multi, 0.0);
	//数据更新，待验证  是否需要进行更新
	if (state)
	{
		updateCurveDate();
		computeLength();
	}
	return state;
}

bool CBSplineCurve::buildCurve()
{
	bool state = createBSplineCurve();
	return state;
}

vector<CPoint> CBSplineCurve::getControlPoints() const
{
	return m_ControlPoints;
}

int CBSplineCurve::getSizeControlPoints() const
{
	if (m_pCurve.IsNull())
	{
		return -1;
	}
	return (int)m_ControlPoints.size();
}

vector<double> CBSplineCurve::getKnots() const
{
	return m_Knots;
}

int CBSplineCurve::getSizeKnots() const
{
	if (m_pCurve.IsNull())
	{
		return -1;
	}
	return (int)m_Knots.size();
}

vector<double> CBSplineCurve::getWeights() const
{
	return m_Weights;
}

int CBSplineCurve::getSizeWeight() const
{
	if (m_pCurve.IsNull())
	{
		return -1;
	}
	return (int)m_Weights.size();
}

vector<int> CBSplineCurve::getMultis() const
{
	return m_Multis;
}

int CBSplineCurve::getSizeMultis() const
{
	if (m_pCurve.IsNull())
	{
		return -1;
	}
	return (int)m_Multis.size();
}

CPoint CBSplineCurve::getSingleControlPoint(int index)
{
	size_t size = m_ControlPoints.size();
	if (m_pCurve.IsNull() || (index < 1 || index > size))
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}
	return m_ControlPoints[index - 1];
}

double CBSplineCurve::getSingleWeight(int index)
{
	size_t size = m_ControlPoints.size();
	if (m_pCurve.IsNull() || (index < 1 || index > size))
	{
		return -1;
	}
	return m_Weights[index - 1];
}

double CBSplineCurve::getSingleKnot(int index)
{
	size_t size = m_Knots.size();
	if (m_pCurve.IsNull() || (index < 1 || index > size))
	{
		return -1;
	}
	return m_Knots[index - 1];
}

int CBSplineCurve::getSingleMulti(int index)
{
	size_t size = m_Multis.size();
	if (m_pCurve.IsNull() || (index < 1 || index > size))
	{
		return -1;
	}
	return m_Multis[index - 1];
}

double CBSplineCurve::getLength() const
{
	return m_dLength;
}

bool CBSplineCurve::getRational() const
{
	return m_bRational;
}

int CBSplineCurve::getDegree() const
{
	return m_iDegree;
}

bool CBSplineCurve::getClosed() const
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	bool state = m_pCurve->IsClosed();
	return state;
}

bool CBSplineCurve::getPeriodic() const
{
	return m_bPeriodic;
}


/*CVector CBSplineCurve::getDerivation1(double u)
{
	if (m_pCurve.IsNull())
	{
		CVector tmp(INT_MAX, INT_MAX, INT_MAX, -1);
		return tmp;
	}
	gp_Vec vec;
	gp_Pnt pnt;
	m_pCurve->D1(u, pnt, vec);
	CVector vector(vec);
	return vector;
}*/

bool CBSplineCurve::getDerivation1(double u, CVector & vector)
{
	if (m_pCurve.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Vec vec;
	gp_Pnt pnt;
	m_pCurve->D1(u, pnt, vec);
	vector.setXYZ(vec.X(), vec.Y(), vec.Z());
	return true;
}

/*CVector CBSplineCurve::getDerivation2(double u)
{
	if (m_pCurve.IsNull())
	{
		CVector vec(INT_MAX, INT_MAX, INT_MAX, -1);
		return vec;
	}
	int con = m_pCurve->Continuity();
	if (con < GeomAbs_C2)
	{
		CVector vec(INT_MAX, INT_MAX, INT_MAX, -1); 
		return vec;
	}
	gp_Pnt p;
	gp_Vec v1, v2;
	m_pCurve->D2(u, p, v1, v2);
	CVector vec(v2);
	return vec;
}*/

bool CBSplineCurve::getDerivation2(double u, CVector & vector)
{
	if (m_pCurve.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}

	int con = m_pCurve->Continuity();
	if (con < GeomAbs_C2)
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Pnt p;
	gp_Vec v1, v2;
	m_pCurve->D2(u, p, v1, v2);
	vector.setXYZ(v2.X(), v2.Y(), v2.Z());
	return true;
}

/*CVector CBSplineCurve::getDerivation3(double u)
{
	if (m_pCurve.IsNull())
	{
		CVector vec(INT_MAX, INT_MAX, INT_MAX, -1);
		return vec;
	}
	int con = m_pCurve->Continuity();
	if (con < GeomAbs_C3)
	{
		CVector vec(INT_MAX, INT_MAX, INT_MAX, -1); 
		return vec;
	}
	gp_Pnt p;
	gp_Vec v1, v2, v3;
	m_pCurve->D3(u, p, v1, v2, v3);
	CVector vec(v2);
	return vec;
}*/

bool CBSplineCurve::getDerivation3(double u, CVector & vector)
{
	if (m_pCurve.IsNull())
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	int con = m_pCurve->Continuity();
	if (con < GeomAbs_C3)
	{
		vector.setXYZ(INT_MAX, INT_MAX, INT_MAX, -1);
		return false;
	}
	gp_Pnt p;
	gp_Vec v1, v2, v3;
	m_pCurve->D3(u, p, v1, v2, v3);
	vector.setXYZ(v3.X(), v3.Y(), v3.Z());
	return true;
}

CPoint CBSplineCurve::getPointAtParam(double u)
{
	if (m_pCurve.IsNull())
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}
	gp_Pnt pnt = m_pCurve->Value(u);
	CPoint point(pnt);
	return point;
}

bool CBSplineCurve::getPointAtParam(double u, CPoint & point)
{
	if (m_pCurve.IsNull())
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		point = tmp;
		return false;
	}
	gp_Pnt pnt = m_pCurve->Value(u);
	point.setXYZ(pnt.X(), pnt.Y(), pnt.Z());
	return true;
}

bool CBSplineCurve::getParamAtPoint(const CPoint & point, double & u)
{
	gp_Pnt pnt;
	pnt = point.getPnt();

	GeomAdaptor_Curve apdatorCurve;
	apdatorCurve.Load(m_pCurve);

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

double CBSplineCurve::getUmin()
{
	double umin = -INT_MAX;
	if (!m_pCurve.IsNull())
	{
		umin = m_pCurve->FirstParameter();
	}
	return umin;
}

double CBSplineCurve::getUmax()
{
	double umax = INT_MAX;
	if (!m_pCurve.IsNull())
	{
		umax = m_pCurve->LastParameter();
	}
	return umax;
}

double CBSplineCurve::getDeltaU()
{
	double umin = -INT_MAX, umax = INT_MAX;
	if (!m_pCurve.IsNull())
	{
		umin = m_pCurve->FirstParameter();
		umax = m_pCurve->LastParameter();
	}

	return umax - umin;
}

CPoint CBSplineCurve::getStartPoint()
{
	if (!m_pCurve.IsNull())
	{
		gp_Pnt pnt = m_pCurve->StartPoint();
		CPoint point(pnt);
		return point;
	}
	else
	{
		CPoint point(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return point;
	}	
}

bool CBSplineCurve::getStartPoint(CPoint & point)
{
	if (!m_pCurve.IsNull())
	{
		gp_Pnt pnt = m_pCurve->StartPoint();
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

CPoint CBSplineCurve::getEndPoint()
{
	if (!m_pCurve.IsNull())
	{
		gp_Pnt pnt = m_pCurve->EndPoint();
		CPoint point(pnt);
		return point;
	}
	else
	{
		CPoint point(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return point;
	}
}

bool CBSplineCurve::getEndPoint(CPoint & point)
{
	if (!m_pCurve.IsNull())
	{
		gp_Pnt pnt = m_pCurve->EndPoint();
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

//bool CBSplineCurve::setTranslation(CGeomPoint & vec)
//{
//	if (m_pCurve.IsNull())
//	{
//		return false;
//	}
//	gp_Trsf translation;
//	gp_Vec vecDir(vec.getX(), vec.getY(), vec.getZ());
//	translation.SetTranslation(vecDir);
//	m_pCurve->Transform(translation);
//	//数据更新
//}

bool CBSplineCurve::setTranslation(CGeomPoint & vec, CBSplineCurve & spline)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	gp_Trsf translation;
	gp_Vec vecDir(vec.getX(), vec.getY(), vec.getZ());
	translation.SetTranslation(vecDir);
	size_t sizeControlPoints = (m_ControlPoints).size();
	vector<CPoint> cp;
	for (int i = 1; i <= sizeControlPoints; ++i)
	{
		CPoint p;
		CPoint p0;
		p0 = (m_ControlPoints)[i - 1];
		gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
		translation.Transforms(coord);
		p.setXYZ(coord.X(), coord.Y(), coord.Z());
		//(*m_ControlPoints)[i - 1].rotate(line, angle, 0, p);
		cp.push_back(p);
	}
	spline.setDegree(m_iDegree);
	spline.setControlPoints(cp);
	spline.setKnots(m_Knots);
	spline.setMultis(m_Multis);
	spline.setWeights(m_Weights);
	spline.setPeriodic(m_bPeriodic);
	spline.buildCurve();
	//if (m_bRational)
	//{
	//	CBSplineCurve cc(m_iDegree, cp, m_Knots, m_Multis, m_Weights, m_bPeriodic);
	//	spline = cc;
	//}
	//else
	//{
	//	CBSplineCurve cc(m_iDegree, cp, m_Knots, m_Multis, m_bPeriodic);
	//	spline = cc;
	//}
	
	//还有一种是直接自己实现对点进行平移 待补充
	return true;
}
//
//void CBSplineCurve::setRotation(CGeomPoint & axisOrigin, CGeomPoint & axisDirection, double angle)
//{
//	//CLine line(axisOrigin, axisDirection);
//	//int sizeControlPoints = (*m_ControlPoints).size();
//	//vector<CPoint> cp;
//
//
//	//for (int i = 1; i <= sizeControlPoints; ++i)
//	//{
//	//	CPoint p;
//	//	gp_Pnt pnt;
//	//	CPoint p0;
//	//	p0 = (*m_ControlPoints)[i - 1];
//	//	(*m_ControlPoints)[i - 1].rotate(line, angle, 0, p);
//	//	cp.push_back(p);
//	//}
//	//m_ControlPoints.reset();
//	//m_ControlPoints = make_shared<vector<CPoint>>(cp);
//	//createBSplineCurve();
//	gp_Trsf rotation;
//	gp_Ax1 axis;
//	axis.SetDirection(gp_Dir(gp_Vec(axisDirection.getX(), axisDirection.getY(), axisDirection.getZ())));
//	axis.SetLocation(gp_Pnt(axisOrigin.getX(), axisOrigin.getY(), axisOrigin.getZ()));
//	rotation.SetRotation(axis, angle);
//	m_pCurve->Transform(rotation);
//}

bool CBSplineCurve::setRotation(CGeomPoint & axisOrigin, CGeomPoint & axisDirection, double angle, CBSplineCurve & spline)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	size_t sizeControlPoints = (m_ControlPoints).size();
	vector<CPoint> cp;
	gp_Trsf rotation;
	gp_Ax1 axis;
	axis.SetDirection(gp_Dir(gp_Vec(axisDirection.getX(), axisDirection.getY(), axisDirection.getZ())));
	axis.SetLocation(gp_Pnt(axisOrigin.getX(), axisOrigin.getY(), axisOrigin.getZ()));
	rotation.SetRotation(axis, angle);

	for (int i = 1; i <= sizeControlPoints; ++i)
	{
		CPoint p;
		CPoint p0;
		p0 = (m_ControlPoints)[i - 1];
		gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
		rotation.Transforms(coord);
		p.setXYZ(coord.X(), coord.Y(), coord.Z());
		//(*m_ControlPoints)[i - 1].rotate(line, angle, 0, p);
		cp.push_back(p);
	}

	spline.setDegree(m_iDegree);
	spline.setControlPoints(cp);
	spline.setKnots(m_Knots);
	spline.setMultis(m_Multis);
	spline.setWeights(m_Weights);
	spline.setPeriodic(m_bPeriodic);
	spline.buildCurve();

	/*if (m_bRational)
	{
		CBSplineCurve cc(m_iDegree, cp, m_Knots, m_Multis, m_Weights, m_bPeriodic);
		spline = cc;
	}
	else
	{
		CBSplineCurve cc(m_iDegree, cp, m_Knots, m_Multis, m_bPeriodic);
		spline = cc;
	}
*/
	//还有一种是直接自己实现对点进行旋转
	return true;
}

bool CBSplineCurve::setTrimmed(double u1, double u2)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	double umin = m_pCurve->FirstParameter();
	double umax = m_pCurve->LastParameter();
	if ((u1 - umin) >= 0.0 && (umax - u2) >= 0.0 &&(u2 - u1) >= 0.0)
	{
		m_pCurve->Segment(u1, u2);
		updateCurveDate();
		return true;
	}

	return false;
}

bool CBSplineCurve::setMirror(CGeom & plane, CBSplineCurve & spline)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}

	size_t sizeControlPoints = (m_ControlPoints).size();
	vector<CPoint> cp;
	for (int i = 1; i <= sizeControlPoints; ++i)
	{
		CPoint p;
		
		(m_ControlPoints)[i - 1].mirror(plane, p);
		cp.push_back(p);
	}
	spline.setDegree(m_iDegree);
	spline.setControlPoints(cp);
	spline.setKnots(m_Knots);
	spline.setMultis(m_Multis);
	spline.setWeights(m_Weights);
	spline.setPeriodic(m_bPeriodic);
	bool flag = spline.buildCurve();
	return flag;
}

void CBSplineCurve::reverse (){
	m_pCurve->Reverse ();
	updateCurveDate ();
}

bool CBSplineCurve::setExtendLength(double length, CBSplineCurve & spline)
{
	if (m_pCurve.IsNull())
	{
		return false;
	}
	double umax = m_pCurve->LastParameter();
	vector<CPoint> points;
	gp_Vec vec;
	gp_Pnt pnt;
	gp_Pnt pp = m_pCurve->Value(umax);
	m_pCurve->D1(umax, pnt, vec);
	CVector vector;
	vector.setXYZ(vec.X(), vec.Y(), vec.Z());
	vector.normalize();

	CPoint p1(pp);
	CPoint p2(vector.getX()*length + pnt.X(), vector.getY()*length + pnt.Y(), vector.getZ()*length + pnt.Z());
	
	points.push_back(p1);
	points.push_back(p2);
	CBSplineCurve lineCurve(points);

	Handle_Geom_BSplineCurve c2 = lineCurve.getGeomCurve();
	GeomConvert_CompCurveToBSplineCurve res;

	bool flag1 = res.Add(m_pCurve, 0.1);
	bool flag2 = res.Add(c2, 0.1);
	Handle_Geom_BSplineCurve p = res.BSplineCurve();
	spline.setGeomCurve(p);
	return flag1 & flag2;
}

Handle_Geom_BSplineCurve CBSplineCurve::getGeomCurve() const
{
	return m_pCurve;
}

void CBSplineCurve::setGeomCurve(Handle_Geom_BSplineCurve &c) 
{
	if (!m_pCurve.IsNull())
	{
		m_pCurve.Nullify();
	}
	m_pCurve = c;
	updateCurveDate();
}

bool CBSplineCurve::createBSplineCurve()
{
	if (m_ControlPoints.size() == 0 || m_Knots.size() == 0 || m_Multis.size() == 0)
	{
		cerr << "BSplineCurve Contruct Error!" << endl;
		return false;
	}
	size_t sizeControlPoints = (m_ControlPoints).size();
	size_t sizeKonts = (m_Knots).size();
	size_t sizeMultis = (m_Multis).size();
	size_t sizeWeights = 0;

	TColgp_Array1OfPnt controlPoints(1, (int)sizeControlPoints);
	TColStd_Array1OfReal knots(1, (int)sizeKonts);
	TColStd_HArray1OfInteger mults(1, (int)sizeMultis);

	if (m_Weights.size() != 0)
	{
		m_bRational = true;
		sizeWeights = (m_Weights).size();
		if (sizeWeights != sizeControlPoints)
		{
			cerr << "BSplineCurve Contruct Error!" << endl;
			return false;
		}
	}


	if (m_iDegree < 1 || m_iDegree > 25)
	{
		cerr << "BSplineCurve Contruct Error! Degree Out of Range" << endl;
		return false;
	}

	if( sizeControlPoints < 2 )
	{
		cerr << "BSplineCurve Contruct Error!" << endl;
		return false;
	}

	if (sizeKonts != sizeMultis)
	{
		cerr << "BSplineCurve Contruct Error!" << endl;
		return false;
	}

	int sum = 0;
	for (size_t i = 0; i < sizeMultis; ++i)
	{
		sum += m_Multis[i];
	}

	if (sum != 1 + m_iDegree + sizeControlPoints)
	{
		cerr << "BSplineCurve Contruct Error!" << endl;
		return false;
	}

	//数据转换
	for (int i = 1; i <= sizeControlPoints; ++i)
	{
		gp_Pnt pnt;
		(m_ControlPoints)[i - 1].getPnt(pnt);
		controlPoints.SetValue(i, pnt);
	}

	for (int i = 1; i <= sizeKonts; ++i)
	{
		knots.SetValue(i, (m_Knots)[i - 1]);
	}

	for (int i = 1; i <= sizeMultis; ++i)
	{
		mults.SetValue(i, (m_Multis)[i - 1]);
	}


	if (m_bRational)
	{
		TColStd_Array1OfReal weights(1, (int)sizeWeights);
		for (int i = 1; i <= sizeWeights; ++i)
		{
			weights.SetValue(i, (m_Weights)[i - 1]);
		}

		try
		{
			if (!m_pCurve.IsNull())
			{
				m_pCurve->Delete();
				m_pCurve.Nullify();
			}
			m_pCurve = new Geom_BSplineCurve(controlPoints, weights, knots, mults, m_iDegree, m_bPeriodic);
			computeLength();
			return true;
		}
		catch (const bad_alloc &)
		{
			cerr << "BSplineCurve Contruct new fail" << endl;
			return false;
		}
		
	}
	else
	{
		try
		{
			if (!m_pCurve.IsNull())
			{
				m_pCurve.Nullify();
			}
			m_pCurve = new Geom_BSplineCurve(controlPoints, knots, mults, m_iDegree, m_bPeriodic);
			computeLength();
			return true;
		}
		catch (const bad_alloc &)
		{
			cerr << "BSplineCurve Contruct new fail" << endl;
			return false;
		}
		
	}	
}

void CBSplineCurve::computeLength()
{
	double umin = m_pCurve->FirstParameter();
	double umax = m_pCurve->LastParameter();
	GeomAdaptor_Curve curve(m_pCurve);
	m_dLength = GCPnts_AbscissaPoint::Length(curve, umin, umax);
}

void CBSplineCurve::updateCurveDate()
{
	if (m_pCurve.IsNull())
	{
		return;
	}

	m_iDegree = m_pCurve->Degree();
	m_bPeriodic = m_pCurve->IsPeriodic();
	updateControlPoints();
	updateKnots();
	updateMultis();
	updateWeights();
	computeLength();
}

void CBSplineCurve::updateControlPoints()
{
	if (m_pCurve.IsNull())
	{
		return;
	}

	int size = m_pCurve->Poles().Size();
	m_ControlPoints.clear();
	m_ControlPoints.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		gp_Pnt pnt;
		pnt = m_pCurve->Poles().Value(i);
		CPoint point(pnt);
		m_ControlPoints[i - 1] = point;

	}
}

void CBSplineCurve::updateWeights()
{
	if (m_pCurve.IsNull())
	{
		return;
	}
	bool isRational = m_pCurve->IsRational();
	
	if (!isRational)
	{
		m_Weights.clear();
		m_bRational = false;
		return;
	}

	int size = m_pCurve->Weights()->Size();
	m_Weights.clear();
	m_Weights.resize(size);
	m_bRational = true;
	for (int i = 1; i <= size; ++i)
	{
		double weight;
		weight = m_pCurve->Weights()->Value(i);
		m_Weights[i - 1] = weight;

	}
}

void CBSplineCurve::updateKnots()
{
	if (m_pCurve.IsNull())
	{
		return;
	}
	int size = m_pCurve->Knots().Size();
	m_Knots.clear();
	m_Knots.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		double knot;
		knot = m_pCurve->Knots().Value(i);
		m_Knots[i - 1] = knot;
	}
}

void CBSplineCurve::updateMultis()
{
	if (m_pCurve.IsNull())
	{
		return;
	}

	int size = m_pCurve->Multiplicities().Size();
	m_Multis.clear();
	m_Multis.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		int multi;
		multi = m_pCurve->Multiplicities().Value(i);
		m_Multis[i - 1] = multi;
	}
}
