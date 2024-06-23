#include "CBSplineSurface.h"
#include "DefineType.h"
#include <GeomAPI_PointsToBSplineSurface.hxx>
#include <gp_Ax1.hxx>
#include <GeomAdaptor_Surface.hxx>
#include <Extrema_ExtPS.hxx>
#include <GeomFill_FillingStyle.hxx>
#include <GeomFill_BSplineCurves.hxx>
#include<GeomFill_Pipe.hxx>
#include<GeomPlate_BuildPlateSurface.hxx>
#include<GeomPlate_MakeApprox.hxx>
#include<GeomPlate_Surface.hxx>
#include<Geom_Curve.hxx>
#include<Geom_BSplineSurface.hxx>
#include<Geom_SurfaceOfRevolution.hxx>
#include<Geom_SurfaceOfRevolution.hxx>
CBSplineSurface::CBSplineSurface()
{
	m_iUDegree = 0;
	m_iVDegree = 0;
	m_iTypeCoord = 1;
	m_dArea = 0;

	m_ControlPoints.clear();
	m_UKnots.clear();
	m_VKnots.clear();
	m_UMultis.clear();
	m_VMultis.clear();
	m_Weights.clear();
	m_bUPeriodic = false;
	m_bVPeriodic = false;
	m_bRational = false;
	m_pSurface = nullptr;
}


CBSplineSurface::~CBSplineSurface()
{
	if (!m_pSurface.IsNull())
	{
		m_pSurface.Nullify();
	}

	m_ControlPoints.clear();
	m_UKnots.clear();
	m_VKnots.clear();
	m_UMultis.clear();
	m_VMultis.clear();
	m_Weights.clear();
}

CBSplineSurface::CBSplineSurface(const CBSplineSurface & splineSurface)
{
	m_iUDegree = splineSurface.m_iUDegree;
	m_iVDegree = splineSurface.m_iVDegree;
	m_iTypeCoord = 1;
	m_ControlPoints = splineSurface.m_ControlPoints;
	m_UKnots = splineSurface.m_UKnots;
	m_VKnots = splineSurface.m_VKnots;
	m_UMultis = splineSurface.m_UMultis;
	m_VMultis = splineSurface.m_VMultis;
	m_bRational = splineSurface.m_bRational;
	m_bUPeriodic = splineSurface.m_bUPeriodic;
	m_bVPeriodic = splineSurface.m_bVPeriodic;
	m_Weights = splineSurface.m_Weights;

	createBSplineSurface();
	if (!m_pSurface.IsNull())
	{
		//computeArea();
	}
}

CBSplineSurface::CBSplineSurface(const int & udegree, const int & vdegree, const vector<vector<CPoint>>& controlPoints, const vector<double>& uknots, const vector<double>& vknots, const vector<int>& umultis, const vector<int>& vmultis, bool uperiodic, bool vperiodic)
{
	m_iUDegree = udegree;
	m_iVDegree = vdegree;
	m_iTypeCoord = 1;
	m_ControlPoints = controlPoints;
	m_UKnots = uknots;
	m_VKnots = vknots;
	m_UMultis = umultis;
	m_VMultis = vmultis;
	m_bRational = false;
	m_bUPeriodic = uperiodic;
	m_bVPeriodic = vperiodic;

	createBSplineSurface();
	if (!m_pSurface.IsNull())
	{
		//computeArea();
	}
}

CBSplineSurface::CBSplineSurface(const int & udegree, const int & vdegree, const vector<vector<CPoint>>& controlPoints, const vector<double>& uknots, const vector<double>& vknots, const vector<int>& umultis, const vector<int>& vmultis, const vector<vector<double>>& weights, bool uperiodic, bool vperiodic)
{
	m_iUDegree = udegree;
	m_iVDegree = vdegree;
	m_iTypeCoord = 1;
	m_ControlPoints = controlPoints;
	m_UKnots = uknots;
	m_VKnots = vknots;
	m_UMultis = umultis;
	m_VMultis = vmultis;
	m_Weights = weights;
	m_bRational = true;
	m_bUPeriodic = uperiodic;
	m_bVPeriodic = vperiodic;

	createBSplineSurface();
	if (!m_pSurface.IsNull())
	{
		//computeArea();
	}
}

CBSplineSurface::CBSplineSurface(const vector<vector<CPoint>>& surfacePoints)
{
	size_t sizeu = surfacePoints.size();
	if (sizeu <= 1)
	{
		return;
	}

	size_t sizev = surfacePoints[0].size();
	if (sizev <= 1)
	{
		return;
	}

	TColgp_Array2OfPnt points = TColgp_Array2OfPnt(1, (int)sizeu, 1, (int)sizev);

	for (int i = 1; i <= sizeu; ++i)
	{
		for (int j = 1; j <= sizev; ++j)
		{
			gp_Pnt p;
			surfacePoints[i - 1][j - 1].getPnt(p);
			points.SetValue(i, j, p);
		}
	}

	GeomAPI_PointsToBSplineSurface aSurface;
	aSurface.Interpolate(points);
	m_pSurface = aSurface.Surface();
	updateSurfaceDate();
}

CBSplineSurface::CBSplineSurface(const vector<vector<CPoint>>& surfacePoints, const vector<vector<double>>& weights)
{
	size_t sizeu = surfacePoints.size();
	if (sizeu <= 1)
	{
		return;
	}

	size_t sizev = surfacePoints[0].size();
	if (sizev <= 1)
	{
		return;
	}

	TColgp_Array2OfPnt points = TColgp_Array2OfPnt(1, (int)sizeu, 1, (int)sizev);

	for (int i = 1; i <= sizeu; ++i)
	{
		for (int j = 1; j <= sizev; ++j)
		{
			gp_Pnt p;
			surfacePoints[i - 1][j - 1].getPnt(p);
			points.SetValue(i, j, p);
		}
	}

	GeomAPI_PointsToBSplineSurface aSurface;
	aSurface.Interpolate(points);
	m_pSurface = aSurface.Surface();

	setWeights(weights);
	updateSurfaceDate();
}

CBSplineSurface::CBSplineSurface(const CBSplineCurve & curve1, const CBSplineCurve & curve2)
{
	GeomFill_FillingStyle type = GeomFill_StretchStyle;
	GeomFill_BSplineCurves mergeCurve(curve1.getGeomCurve(), curve2.getGeomCurve(), type);
	m_pSurface = mergeCurve.Surface();
	updateSurfaceDate();
}

CBSplineSurface::CBSplineSurface (const CBSplineCurve &curve1,const CBSplineCurve &curve2,const CBSplineCurve &curve3,const CBSplineCurve &curve4){
	GeomFill_FillingStyle type = GeomFill_StretchStyle;
	try{
		GeomFill_BSplineCurves mergeCurve (curve1.getGeomCurve (),curve2.getGeomCurve (),curve3.getGeomCurve (),curve4.getGeomCurve (),type);
		m_pSurface = mergeCurve.Surface ();
	}
	catch(...){
		cerr << "曲线不封闭" << endl;
	}
	
	updateSurfaceDate ();
}

CBSplineSurface::CBSplineSurface (CBSplineCurve &pathcurve, double radius){
	
	GeomFill_Pipe sweepSurface (pathcurve.getGeomCurve(),radius);
	sweepSurface.Perform ();
	GeomAdaptor_Surface surface(sweepSurface.Surface());
	m_pSurface=surface.BSpline ();
	updateSurfaceDate ();
}

CBSplineSurface::CBSplineSurface ( CBSplineCurve &pathcurve, CBSplineCurve &sectioncurve,bool sweep){
	if(sweep){
		GeomFill_Pipe sweepSurface (pathcurve.getGeomCurve (),sectioncurve.getGeomCurve ());
		sweepSurface.Perform ();
		GeomAdaptor_Surface surface (sweepSurface.Surface ());
		m_pSurface = surface.BSpline ();
		updateSurfaceDate ();
	}
	else
		cout << "sweep参数错误，构造失败";

}

int CBSplineSurface::getSurfaceType()
{
	return BSPLINESURFACE;
}

bool CBSplineSurface::setInterpolatePoints(const vector<vector<CPoint>>& surfacePoints)
{
	size_t sizeu = surfacePoints.size();
	if (sizeu <= 1)
	{
		return false;
	}

	size_t sizev = surfacePoints[0].size();
	if (sizev <= 1)
	{
		return false;
	}

	TColgp_Array2OfPnt points = TColgp_Array2OfPnt(1, (int)sizeu, 1, (int)sizev);

	for (int i = 1; i <= sizeu; ++i)
	{
		for (int j = 1; j <= sizev; ++j)
		{
			gp_Pnt p;
			surfacePoints[i - 1][j - 1].getPnt(p);
			points.SetValue(i, j, p);
		}
	}

	GeomAPI_PointsToBSplineSurface aSurface;
	aSurface.Interpolate(points);
	m_pSurface = aSurface.Surface();
	updateSurfaceDate();
	return true;
}

//bool CBSplineSurface::PlateSurface (CBSplineSurface *bsplinesurface,int degree){
//	GeomPlate_BuildPlateSurface plateSurface (bsplinesurface->getGeomSurface (),degree);
//	Handle(GeomPlate_Surface)surface=plateSurface.Surface ();
//	GeomPlate_MakeApprox appsurface (surface,1.e-15,10,4,5);
//	m_pSurface = appsurface.Surface ();
//	updateSurfaceDate ();
//	return true;
//}

bool CBSplineSurface::setUDegree(int udegree)
{
	if (m_pSurface.IsNull())
	{
		m_iUDegree = udegree;
		return true;
	}
	else
	{
		return false;
		//暂时屏蔽升阶
	}
}

bool CBSplineSurface::setVDegree(int vdegree)
{
	if (m_pSurface.IsNull())
	{
		m_iVDegree = vdegree;
		return false;
	}
	else
	{
		return false;
		//暂时屏蔽升阶
	}
}

bool CBSplineSurface::setControlPoints(const vector<vector<CPoint>>& controlPoints)
{
	if (m_pSurface.IsNull())
	{
		m_ControlPoints = controlPoints;
		return true;
	}
	else
	{
		size_t usize = m_ControlPoints.size();
		size_t vsize = m_ControlPoints[0].size();
		if (usize == controlPoints.size() && vsize == controlPoints[0].size())
		{
			m_ControlPoints = controlPoints;
			for (int i = 1; i <= vsize; ++i)
			{
				for (int j = 1; j <= usize; ++j)
				{
					gp_Pnt p = controlPoints[j - 1][i - 1].getPnt();
					m_pSurface->SetPole(j, i, p);
				}
				
			}
			//computeArea();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineSurface::setControlPoints(const vector<vector<CPoint>>& controlPoints, const vector<vector<double>>& weights)
{
	size_t sizewu = weights.size();
	size_t sizewv = weights[0].size();
	size_t sizepu = controlPoints.size();
	size_t sizepv = controlPoints[0].size();
	if (sizewu != sizepu || sizewv != sizepv)
	{
		return false;
	}
	for (int i = 0; i < sizewu; ++i)
	{
		for (int j = 0; j < sizewv; ++j)
		{
			if (weights[i][j] < 1e-10)
			{
				return false;
			}
		}
	}

	if (m_pSurface.IsNull())
	{
		m_ControlPoints = controlPoints;
		m_Weights = weights;
		m_bRational = true;
		return true;
	}
	else
	{
		size_t sizeu = m_ControlPoints.size();
		size_t sizev = m_ControlPoints[0].size();
		if (sizeu != sizepu || sizev != sizepv)
		{
			return false;
		}
		
		m_ControlPoints = controlPoints;
		m_Weights = weights;
		m_bRational = true;
		for (int i = 1; i <= sizeu; ++i)
		{
			for (int j = 1; j <= sizev; ++j)
			{
				gp_Pnt p = controlPoints[i - 1][j - 1].getPnt();
				m_pSurface->SetPole(i, j, p, weights[i - 1][j - 1]);
			}
			
		}
		//computeArea();
		return true;
	}
}

bool CBSplineSurface::setSingleControlPoint(int & uindex, int & vindex, const CPoint & controlPoint)
{
	size_t sizeu = m_ControlPoints.size();
	size_t sizev = m_ControlPoints[0].size();
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev))
	{
		return false;
	}
	m_ControlPoints[uindex - 1][vindex - 1] = controlPoint;
	gp_Pnt p;
	controlPoint.getPnt(p);
	m_pSurface->SetPole(uindex, vindex, p);
	//面积更新
	//computeArea();
	return true;
}

bool CBSplineSurface::setSingleControlPoint(int & uindex, int & vindex, const CPoint & controlPoint, double & weight)
{
	size_t sizeu = m_ControlPoints.size();
	size_t sizev = m_ControlPoints[0].size();
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev))
	{
		return false;
	}
	if (weight < 1e-10)
	{
		return false;
	}
	m_ControlPoints[uindex - 1][vindex - 1] = controlPoint;
	m_Weights[uindex - 1][vindex - 1] = weight;
	gp_Pnt p;
	controlPoint.getPnt(p);
	m_pSurface->SetPole(uindex, vindex, p, weight);
	//面积更新
	//computeArea();
	return true;
}

bool CBSplineSurface::setWeights(const vector<vector<double>>& weights)
{
	size_t sizeu = weights.size();
	if (sizeu == 0)
	{
		return false;
	}
	size_t sizev = weights[0].size();
	if (sizev == 0)
	{
		return false;
	}
	for (int i = 0; i < sizeu; ++i)
	{
		for (int j = 0; j < sizev; ++j)
		{
			if (weights[i][j] < 1e-10)
			{
				return false;
			}
		}		
	}
	if (m_pSurface.IsNull())
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
		size_t size2 = m_Weights[0].size();
		if (size == sizeu && size2 == sizev)
		{
			m_Weights = weights;
			for (int i = 1; i <= size; ++i)
			{
				for (int j = 1; j <= sizev; ++j)
				{
					m_pSurface->SetWeight(i,j, weights[i-1][j - 1]);
				}
				
			}
			//computeArea();
			return true;
		}
		else
		{
			return false;
		}
	}
}

bool CBSplineSurface::setSingleWeight(int uindex, int vindex, double weight)
{
	size_t sizeu = m_Weights.size();
	size_t sizev = m_Weights[0].size();
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev) && !m_bRational)
	{
		return false;
	}

	if (weight < 1e-10)
	{
		return false;
	}

	m_Weights[uindex - 1][vindex - 1] = weight;
	m_pSurface->SetWeight(uindex, vindex, weight);
	//computeArea();
	return true;
}

bool CBSplineSurface::setUKnots(const vector<double>& uknots)
{
	//递增序列检查
	size_t size = uknots.size();
	for (int i = 0; i < size - 1; ++i)
	{
		if (uknots[i] > uknots[i + 1])
			return false;
	}
	if (m_pSurface.IsNull())
	{
		m_UKnots = uknots;
		return true;
	}
	size = m_UKnots.size();
	if (size == uknots.size())
	{
		m_UKnots = uknots;
		for (int i = 1; i <= size; ++i)
		{
			m_pSurface->SetUKnot(i, uknots[i - 1]);
		}
		//computeArea();
		return true;
	}
	return false;
}

bool CBSplineSurface::setVKnots(const vector<double>& vknots)
{
	//递增序列检查
	size_t size = vknots.size();
	for (int i = 0; i < size - 1; ++i)
	{
		if (vknots[i] > vknots[i + 1])
			return false;
	}
	if (m_pSurface.IsNull())
	{
		m_VKnots = vknots;
		return true;
	}
	size = m_VKnots.size();
	if (size == vknots.size())
	{
		m_VKnots = vknots;
		for (int i = 1; i <= size; ++i)
		{
			m_pSurface->SetVKnot(i, vknots[i - 1]);
		}
		//computeArea();
		return true;
	}
	return false;
}

bool CBSplineSurface::setSingleUKnot(int uindex, double uknot)
{
	size_t size = m_UKnots.size();
	if (!m_pSurface.IsNull() && (uindex < 1 || uindex > size))
	{
		return false;
	}
	//递增序列检查
	if (uindex != size && uknot > m_UKnots[uindex])
		return false;
	if (uindex != 1 && uknot < m_UKnots[uindex - 2])
		return false;

	m_UKnots[uindex - 1] = uknot;
	m_pSurface->SetUKnot(uindex, uknot);
	//computeArea()
	return true;
}

bool CBSplineSurface::setSingleVKnot(int vindex, double vknot)
{
	size_t size = m_VKnots.size();
	if (!m_pSurface.IsNull() && (vindex < 1 || vindex > size))
	{
		return false;
	}
	//递增序列检查
	if (vindex != size && vknot > m_UKnots[vindex])
		return false;
	if (vindex != 1 && vknot < m_UKnots[vindex - 2])
		return false;

	m_VKnots[vindex - 1] = vknot;
	m_pSurface->SetVKnot(vindex, vknot);
	//computeArea()
	return true;;
}

bool CBSplineSurface::setUMultis(const vector<int>& umultis)
{
	if (m_pSurface.IsNull())
	{
		m_UMultis = umultis;
		return true;
	}
	size_t size = m_UMultis.size();
	if (size == umultis.size())
	{
		bool flag = false;
		for (int i = 1; i <= size; ++i)
		{
			if (m_UMultis[i - 1] < umultis[i - 1])
			{
				m_pSurface->IncreaseUMultiplicity(i, umultis[i - 1]);  //类似degree  只能增加至degree  不能减少
				m_UMultis[i - 1] = umultis[i - 1];
				flag = true;
			}

		}
		if (flag)
		{
			//数据更新
			updateSurfaceDate();
			//computeArea()
		}

		return true;
	}
	return false;
}

bool CBSplineSurface::setVMultis(const vector<int>& vmultis)
{
	if (m_pSurface.IsNull())
	{
		m_VMultis = vmultis;
		return true;
	}
	size_t size = m_VMultis.size();
	if (size == vmultis.size())
	{
		bool flag = false;
		for (int i = 1; i <= size; ++i)
		{
			if (m_VMultis[i - 1] < vmultis[i - 1])
			{
				m_pSurface->IncreaseVMultiplicity(i, vmultis[i - 1]);  //类似degree  只能增加至degree  不能减少
				m_VMultis[i - 1] = vmultis[i - 1];
				flag = true;
			}

		}
		if (flag)
		{
			//数据更新
			updateSurfaceDate();
			//computeArea()
		}

		return true;
	}
	return false;
}

bool CBSplineSurface::setSingleUMultis(int uindex, int multi)
{
	size_t size = m_UKnots.size();
	if (m_pSurface.IsNull() && (uindex < 1 || uindex > size))
	{
		return false;
	}
	if (m_UMultis[uindex - 1] >= multi  || multi > m_iUDegree)
	{
		return false;
	}
	m_UMultis[uindex - 1] = multi;
	m_pSurface->IncreaseUMultiplicity(uindex, multi);
	//数据更新
	updateSurfaceDate();
	//computeArea()

	return true;
}

bool CBSplineSurface::setSingleVMultis(int vindex, int multi)
{
	size_t size = m_VKnots.size();
	if (m_pSurface.IsNull() && (vindex < 1 || vindex > size))
	{
		return false;
	}
	if (m_VMultis[vindex - 1] >= multi || multi > m_iVDegree)
	{
		return false;
	}
	m_VMultis[vindex - 1] = multi;
	m_pSurface->IncreaseVMultiplicity(vindex, multi);
	//数据更新
	updateSurfaceDate();
	//computeArea()

	return true;
}

bool CBSplineSurface::setUPeriodic(bool uperiodic)
{
	if (m_bUPeriodic != uperiodic)
	{
		if (uperiodic)
		{
			m_pSurface->SetUPeriodic();
		}
		else
		{
			m_pSurface->SetUNotPeriodic();
		}
		//数据更新，待补充
		updateSurfaceDate();
		m_bUPeriodic = uperiodic;
		//computeArea();
		return true;
	}
	return false;
}

bool CBSplineSurface::setVPeriodic(bool vperiodic)
{
	if (m_bVPeriodic != vperiodic)
	{
		if (vperiodic)
		{
			m_pSurface->SetVPeriodic();
		}
		else
		{
			m_pSurface->SetVNotPeriodic();
		}
		//数据更新，待补充
		updateSurfaceDate();
		m_bVPeriodic = vperiodic;
		//computeArea();
		return true;
	}
	return false;
}

bool CBSplineSurface::insertUKnot(double u, int umulti)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	double umin, umax, vmin, vmax;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	if ((u - umin) >= 0.0 && (umax - u) >= 0.0)
	{
		m_pSurface->InsertUKnot(u, umulti, 1.e-7);
		//数据更新
		updateSurfaceDate();
		//computeArea();
		return true;
	}

	return false;
}

bool CBSplineSurface::removeUKnot(int uindex, int umulti)
{
	size_t size = m_UKnots.size();
	if (m_pSurface.IsNull() || uindex < 1 || uindex > size)
	{
		return false;
	}
	bool state = m_pSurface->RemoveUKnot(uindex, umulti, 0.0);
	//数据更新，待验证  是否需要进行更新
	if (state)
	{
		updateSurfaceDate();
		//computeArea();
	}
	return state;
}

bool CBSplineSurface::insertVKnot(double v, int vmulti)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	double umin, umax, vmin, vmax;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	if ((v - vmin) >= 0.0 && (vmax - v) >= 0.0)
	{
		m_pSurface->InsertVKnot(v, vmulti, 1.e-7);
		//数据更新
		updateSurfaceDate();
		//computeArea();
		return true;
	}

	return false;
}

bool CBSplineSurface::removeVKnot(int vindex, int vmulti)
{
	size_t size = m_VKnots.size();
	if (m_pSurface.IsNull() || vindex < 1 || vindex > size)
	{
		return false;
	}
	bool state = m_pSurface->RemoveVKnot(vindex, vmulti, 0.0);
	//数据更新，待验证  是否需要进行更新
	if (state)
	{
		updateSurfaceDate();
		//computeArea();
	}
	return state;
}

bool CBSplineSurface::buildSurface()
{
	bool state = createBSplineSurface();
	return state;
}

vector<vector<CPoint>> CBSplineSurface::getControlPoints() const
{
	return m_ControlPoints;
}

int CBSplineSurface::getSizeUControlPoints() const
{
	if (m_ControlPoints.empty())
	{
		return -1;
	}
	size_t size =  m_ControlPoints.size();
	return (int)size;
}

int CBSplineSurface::getSizeVControlPoints() const
{
	if (m_ControlPoints.empty())
	{
		return -1;
	}
	size_t size = m_ControlPoints[0].size();
	return (int)size;
}

vector<double> CBSplineSurface::getUKnots() const
{
	return m_UKnots;
}

int CBSplineSurface::getSizeUKnots() const
{
	if (m_UKnots.empty())
	{
		return -1;
	}
	size_t size =  m_UKnots.size();
	return (int)size;
}

int CBSplineSurface::getSizeVKnots() const
{
	if (m_VKnots.empty())
	{
		return -1;
	}
	size_t size = m_VKnots.size();
	return (int)size;
}

vector<double> CBSplineSurface::getVKnots() const
{
	return m_VKnots;
}

vector<vector<double>> CBSplineSurface::getWeights() const
{
	return m_Weights;
}

int CBSplineSurface::getSizeUWeights() const
{
	if (m_Weights.empty())
	{
		return -1;
	}
	size_t size = m_Weights.size();
	return (int)size;
}

int CBSplineSurface::getSizeVWeights() const
{
	if (m_Weights.empty())
	{
		return -1;
	}
	size_t size = m_Weights[0].size();
	return (int)size;
}

vector<int> CBSplineSurface::getUMultis() const
{
	return m_UMultis;
}

int CBSplineSurface::getSizeUMultis() const
{
	if (m_UMultis.empty())
	{
		return -1;
	}
	size_t size = m_UMultis.size();
	return (int)size;
}

vector<int> CBSplineSurface::getVMultis() const
{
	return m_VMultis;
}

int CBSplineSurface::getSizeVMultis() const
{
	if (m_VMultis.empty())
	{
		return -1;
	}
	size_t size = m_VMultis.size();
	return (int)size;
}

CPoint CBSplineSurface::getSingleControlPoint(int uindex, int vindex)
{
	size_t size = m_ControlPoints.size();	
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > size))
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}

	size_t vsize = m_ControlPoints[0].size();
	if ((vindex < 1 || vindex > vsize))
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}
	return m_ControlPoints[uindex - 1][vindex - 1];
}

double CBSplineSurface::getSingleWeight(int uindex, int vindex)
{
	size_t size = m_ControlPoints.size();
	
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > size) )
	{
		return -1;
	}

	size_t vsize = m_ControlPoints[0].size();
	if ((vindex < 1 || vindex > vsize))
	{
		return -1;
	}
	return m_Weights[uindex - 1][vindex - 1];
}

double CBSplineSurface::getSingleUKnot(int uindex)
{
	size_t size = m_UKnots.size();
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > size))
	{
		return -1;
	}
	return m_UKnots[uindex - 1];
}

double CBSplineSurface::getSingleVKnot(int vindex)
{
	size_t size = m_VKnots.size();
	if (m_pSurface.IsNull() || (vindex < 1 || vindex > size))
	{
		return -1;
	}
	return m_VKnots[vindex - 1];
}

int CBSplineSurface::getSingleUMulti(int uindex)
{
	size_t size = m_UMultis.size();
	if (m_pSurface.IsNull() || (uindex < 1 || uindex > size))
	{
		return -1;
	}
	return m_UMultis[uindex - 1];
}

int CBSplineSurface::getSingleVMulti(int vindex)
{
	size_t size = m_VMultis.size();
	if (m_pSurface.IsNull() || (vindex < 1 || vindex > size))
	{
		return -1;
	}
	return m_VMultis[vindex - 1];
}

double CBSplineSurface::getArea() const
{
	return m_dArea;
}

bool CBSplineSurface::getRational() const
{
	return m_bRational;
}

int CBSplineSurface::getUDegree() const
{
	return m_iUDegree;
}

int CBSplineSurface::getVDegree() const
{
	return m_iVDegree;
}

bool CBSplineSurface::getUClosed() const
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	bool state = m_pSurface->IsUClosed();
	return state;
}

bool CBSplineSurface::getVClosed() const
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	bool state = m_pSurface->IsVClosed();
	return state;
}

bool CBSplineSurface::getUPeriodic() const
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	return m_bUPeriodic;
}

bool CBSplineSurface::getVPeriodic() const
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	return m_bVPeriodic;
}

vector<CVector> CBSplineSurface::getDerivation1(double u, double v)
{
	vector<CVector> tmp;
	if (m_pSurface.IsNull())
	{		
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		return tmp;
	}
	gp_Vec vecU;
	gp_Vec vecV;
	gp_Pnt pnt;
	m_pSurface->D1(u, v, pnt, vecU, vecV);
	CVector vectorU(vecU);
	CVector vectorV(vecV);
	tmp.push_back(vectorU);
	tmp.push_back(vectorV);
	return tmp;
}

bool CBSplineSurface::getDerivation1(double u, double v, vector<CVector>& vectors)
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

vector<CVector> CBSplineSurface::getDerivation2(double u, double v)
{
	vector<CVector> tmp;
	if (m_pSurface.IsNull())
	{
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		return tmp;
	}
	gp_Vec vecU1, vecU2;
	gp_Vec vecV1, vecV2;
	gp_Vec vecUV;
	gp_Pnt pnt;
	m_pSurface->D2(u, v, pnt, vecU1, vecV1, vecU2, vecV2, vecUV);
	CVector vectorU2(vecU2);
	CVector vectorV2(vecV2);
	CVector vectorUV(vecUV);
	tmp.push_back(vectorU2);
	tmp.push_back(vectorV2);
	tmp.push_back(vectorUV);
	return tmp;
}

bool CBSplineSurface::getDerivation2(double u, double v, vector<CVector>& vectors)
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

vector<CVector> CBSplineSurface::getDerivation3(double u, double v)
{
	vector<CVector> tmp;
	if (m_pSurface.IsNull())
	{
		CVector tmpV(INT_MAX, INT_MAX, INT_MAX, -1);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		tmp.push_back(tmpV);
		return tmp;
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
	tmp.push_back(vectorU3);
	tmp.push_back(vectorV3);
	tmp.push_back(vectorUUV);
	tmp.push_back(vectorUVV);
	return tmp;
}

bool CBSplineSurface::getDerivation3(double u, double v, vector<CVector>& vectors)
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

CPoint CBSplineSurface::getPointAtParam(double u, double v)
{
	CPoint pnt;
	getPointAtParam(u, v, pnt);
	return pnt;
}

bool CBSplineSurface::getPointAtParam(double u, double v, CPoint & point)
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

bool CBSplineSurface::getParamAtPoint(const CPoint & point, double & u, double & v)
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
	extPS.Initialize(apdatorSurface, umin, umax, vmin, vmax, 1.e-7, 1.e-7);
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

bool CBSplineSurface::getUVMinMax(double & umin, double & umax, double & vmin, double & vmax)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}

	m_pSurface->Bounds(umin,umax,vmin,vmax);
	return true;
}

bool CBSplineSurface::setTranslation(CGeomPoint & vec, CBSplineSurface & surface)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	gp_Trsf translation;
	gp_Vec vecDir(vec.getX(), vec.getY(), vec.getZ());
	translation.SetTranslation(vecDir);
	size_t sizeUControlPoints = (m_ControlPoints).size();
	size_t sizeVControlPoints = (m_ControlPoints[0]).size();
	vector<vector<CPoint>> cp;
	for (int i = 1; i <= sizeUControlPoints; ++i)
	{
		vector<CPoint> tmp;
		for (int j = 1; j <= sizeVControlPoints; ++j)
		{
			CPoint p;
			CPoint p0;
			p0 = (m_ControlPoints)[i - 1][j - 1];
			gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
			translation.Transforms(coord);
			p.setXYZ(coord.X(), coord.Y(), coord.Z());
			//(*m_ControlPoints)[i - 1].rotate(line, angle, 0, p);
			tmp.push_back(p);
		}
		cp.push_back(tmp);
		
	}
	surface.setUDegree(m_iUDegree);
	surface.setVDegree(m_iVDegree);
	surface.setControlPoints(cp);
	surface.setUKnots(m_UKnots);
	surface.setVKnots(m_VKnots);
	surface.setUMultis(m_UMultis);
	surface.setVMultis(m_VMultis);
	surface.setWeights(m_Weights);
	surface.setUPeriodic(m_bUPeriodic);
	surface.setVPeriodic(m_bVPeriodic);
	surface.buildSurface();

	return true;
}

bool CBSplineSurface::setRotation(CGeomPoint & axisOrigin, CGeomPoint & axisDirection, double angle, CBSplineSurface & surface)
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
	size_t sizeUControlPoints = (m_ControlPoints).size();
	size_t sizeVControlPoints = (m_ControlPoints[0]).size();
	vector<vector<CPoint>> cp;
	for (int i = 1; i <= sizeUControlPoints; ++i)
	{
		vector<CPoint> tmp;
		for (int j = 1; j <= sizeVControlPoints; ++j)
		{
			CPoint p;
			CPoint p0;
			p0 = (m_ControlPoints)[i - 1][j - 1];
			gp_XYZ coord(p0.getX(), p0.getY(), p0.getZ());
			rotation.Transforms(coord);
			p.setXYZ(coord.X(), coord.Y(), coord.Z());
			tmp.push_back(p);
		}
		cp.push_back(tmp);

	}
	surface.setUDegree(m_iUDegree);
	surface.setVDegree(m_iVDegree);
	surface.setControlPoints(cp);
	surface.setUKnots(m_UKnots);
	surface.setVKnots(m_VKnots);
	surface.setUMultis(m_UMultis);
	surface.setVMultis(m_VMultis);
	surface.setWeights(m_Weights);
	surface.setUPeriodic(m_bUPeriodic);
	surface.setVPeriodic(m_bVPeriodic);
	surface.buildSurface();

	return true;
}

bool CBSplineSurface::setTrimmed(double u1, double u2, double v1, double v2)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}
	double umin, umax, vmin, vmax;
	m_pSurface->Bounds(umin, umax, vmin, vmax);
	if ((u1 - umin) >= 0.0 && (umax - u2) >= 0.0 && (u2 - u1) >= 0.0)
	{
		if ((v1 - vmin) >= 0.0 && (vmax - v2) >= 0.0 && (v2 - v1) >= 0.0)
		{
			m_pSurface->Segment(u1, u2, v1, v2);
			updateSurfaceDate();
			return true;
		}		
	}
	return false;
}

bool CBSplineSurface::setMirror(CGeom & plane, CBSplineSurface & surface)
{
	if (m_pSurface.IsNull())
	{
		return false;
	}

	size_t sizeUControlPoints = (m_ControlPoints).size();
	size_t sizeVControlPoints = (m_ControlPoints[0]).size();
	vector<vector<CPoint>> cp;
	for (int i = 1; i <= sizeUControlPoints; ++i)
	{
		vector<CPoint> tmp;
		for (int j = 1; j <= sizeVControlPoints; ++j)
		{
			CPoint p;
			(m_ControlPoints)[i - 1][j - 1].mirror(plane, p);
			
			tmp.push_back(p);
		}
		cp.push_back(tmp);

	}
	surface.setUDegree(m_iUDegree);
	surface.setVDegree(m_iVDegree);
	surface.setControlPoints(cp);
	surface.setUKnots(m_UKnots);
	surface.setVKnots(m_VKnots);
	surface.setUMultis(m_UMultis);
	surface.setVMultis(m_VMultis);
	surface.setWeights(m_Weights);
	surface.setUPeriodic(m_bUPeriodic);
	surface.setVPeriodic(m_bVPeriodic);
	bool flag = surface.buildSurface();
	return flag;
}

Handle_Geom_BSplineSurface CBSplineSurface::getGeomSurface() const
{
	return m_pSurface;
}

void CBSplineSurface::setGeomSurface (Handle_Geom_BSplineSurface surface){
	m_pSurface = surface;
	updateSurfaceDate ();
}

bool CBSplineSurface::createBSplineSurface()
{
	if (m_ControlPoints.size() == 0 || m_ControlPoints[0].size() == 0 || m_UKnots.size() == 0 || m_UMultis.size() == 0 || m_VKnots.size() ==0 || m_VMultis.size() == 0)
	{
		cerr << "BSplineSurface Contruct Error!1" << endl;
		return false;
	}
	size_t sizeUControlPoints = (m_ControlPoints).size();
	size_t sizeVControlPoints = (m_ControlPoints[0]).size();
	size_t sizeUKonts = (m_UKnots).size();
	size_t sizeUMultis = (m_UMultis).size();
	size_t sizeVKonts = (m_VKnots).size();
	size_t sizeVMultis = (m_VMultis).size();
	size_t sizeUWeights = 0, sizeVWeights = 0;


	TColgp_Array2OfPnt controlPoints(1, (int)sizeUControlPoints, 1, (int)sizeVControlPoints);
	TColStd_Array1OfReal uknots(1, (int)sizeUKonts);
	TColStd_HArray1OfInteger umults(1, (int)sizeUMultis);
	TColStd_Array1OfReal vknots(1, (int)sizeVKonts);
	TColStd_HArray1OfInteger vmults(1, (int)sizeVMultis);

	if (m_Weights.size() != 0 && m_Weights[0].size() != 0)
	{
		m_bRational = true;
		sizeUWeights = (m_Weights).size();
		sizeVWeights = (m_Weights[0]).size();
		if (sizeUWeights != sizeUControlPoints || sizeVWeights != sizeVControlPoints)
		{
			cerr << "BSplineSurface Contruct Error!2" << endl;
			return false;
		}
	}


	if (m_iUDegree < 1 || m_iUDegree > 25 || m_iVDegree < 1 || m_iVDegree > 25)
	{
		cerr << "BSplineSurface Contruct Error! Degree Out of Range" << endl;
		return false;
	}

	if (sizeUControlPoints < 2 || sizeVControlPoints < 2)
	{
		cerr << "BSplineSurface Contruct Error!3" << endl;
		return false;
	}

	if (sizeUKonts != sizeUMultis || sizeVKonts != sizeVMultis)
	{
		cerr << "BSplineSurface Contruct Error!4" << endl;
		return false;
	}

	//int sumu = 0, sumv = 0;
	//for (int i = 0; i < sizeUMultis; ++i)
	//{
	//	sumu += m_UMultis[i];
	//}

	/*if (sumu != 1 + m_iUDegree + sizeUControlPoints)
	{
		cerr << "BSplineSurface Contruct Error!" << endl;
		return false;
	}

	for (int i = 0; i < sizeVMultis; ++i)
	{
		sumv += m_VMultis[i];
	}

	if (sumv != 1 + m_iVDegree + sizeVControlPoints)
	{
		cerr << "BSplineSurface Contruct Error!" << endl;
		return false;
	}*/

	//数据转换
	for (int i = 1; i <= sizeUControlPoints; ++i)
	{
		for (int j = 1; j <= sizeVControlPoints; ++j)
		{
			gp_Pnt pnt;
			(m_ControlPoints)[i - 1][j - 1].getPnt(pnt);
			controlPoints.SetValue(i, j, pnt);
		}		
	}

	for (int i = 1; i <= sizeUKonts; ++i)
	{
		uknots.SetValue(i, (m_UKnots)[i - 1]);
	}

	for (int i = 1; i <= sizeUMultis; ++i)
	{
		umults.SetValue(i, (m_UMultis)[i - 1]);
	}

	for (int i = 1; i <= sizeVKonts; ++i)
	{
		vknots.SetValue(i, (m_VKnots)[i - 1]);
	}

	for (int i = 1; i <= sizeVMultis; ++i)
	{
		vmults.SetValue(i, (m_VMultis)[i - 1]);
	}


	if (m_bRational)
	{
		TColStd_Array2OfReal  weights(1, (int)sizeUWeights, 1, (int)sizeVWeights);
		for (int i = 1; i <= sizeUWeights; ++i)
		{
			for (int j = 1; j <= sizeVWeights; ++j)
			{
				weights.SetValue(i, j, (m_Weights)[i - 1][j - 1]);
			}			
		}

		try
		{
			if (!m_pSurface.IsNull())
			{
				m_pSurface.Nullify();
			}
			m_pSurface = new Geom_BSplineSurface(controlPoints, weights, uknots, vknots, umults, vmults, m_iUDegree, m_iVDegree, m_bUPeriodic, m_bVPeriodic);
			//computeArea();
			return true;
		}
		catch (const bad_alloc &)
		{
			cerr << "BSplineSurface Contruct new fail" << endl;
			return false;
		}

	}
	else
	{
		try
		{
			if (!m_pSurface.IsNull())
			{
				m_pSurface.Nullify();
			}
			m_pSurface = new Geom_BSplineSurface(controlPoints, uknots, vknots, umults, vmults, m_iUDegree, m_iVDegree, m_bUPeriodic, m_bVPeriodic);
			//computeArea();
			return true;
		}
		catch (const bad_alloc &)
		{
			cerr << "BSplineSurface Contruct new fail" << endl;
			return false;
		}

	}
}

void CBSplineSurface::getBoundCurves (vector<CBSplineCurve> &curves){
	if(m_pSurface.IsNull ()){
		return;
	}

	double umin,umax,vmin,vmax;
	m_pSurface->Bounds (umin,umax,vmin,vmax);
	vector<vector<CPoint>> points;
	points.resize (4);
	for(int i = 0; i < 4; ++i){
		points[i].resize (101);
	}
	double deltaU = (double) (umax - umin) / 100.0;
	double deltaV = (double) (vmax - vmin) / 100.0;
	for(int i = 0; i <= 100; ++i){
		points[0][i] = CPoint (m_pSurface->Value (umin + deltaU * i,vmin));
		points[1][i] = CPoint (m_pSurface->Value (umin + deltaU * i,vmax));

		points[2][i] = CPoint (m_pSurface->Value (umin,vmin + deltaV * i));
		points[3][i] = CPoint (m_pSurface->Value (umax,vmin + deltaV * i));
	}

	for(int i = 0; i < 4; ++i){
		CBSplineCurve tmpCurve (points[i]);
		curves.push_back (tmpCurve);
	}

}

void CBSplineSurface::updateSurfaceDate()
{
	if (m_pSurface.IsNull())
	{
		return;
	}

	m_iUDegree = m_pSurface->UDegree();
	m_iVDegree = m_pSurface->VDegree();
	m_bUPeriodic = m_pSurface->IsUPeriodic();
	m_bVPeriodic = m_pSurface->IsVPeriodic();
	updateControlPoints();
	updateWeights();
	updateUKnots();
	updateVKnots();
	updateUMultis();
	updateVMultis();
}

void CBSplineSurface::updateControlPoints()
{
	if (m_pSurface.IsNull())
	{
		return;
	}

	int usize = m_pSurface->NbUPoles();
	int vsize = m_pSurface->NbVPoles();
	m_ControlPoints.clear();
	m_ControlPoints.resize(usize);

	for (int i = 1; i <= usize; ++i)
	{
		m_ControlPoints[i - 1].resize(vsize);
		for (int j = 1; j <= vsize; ++j)
		{
			gp_Pnt pnt;
			pnt = m_pSurface->Pole(i, j);
			CPoint point(pnt);
			m_ControlPoints[i - 1][j - 1] = point;
		}
	}
}

void CBSplineSurface::updateWeights()
{
	if (m_pSurface.IsNull())
	{
		return;
	}

	bool isURational = m_pSurface->IsURational();
	bool isVRational = m_pSurface->IsVRational();
	if (!isURational && !isVRational)
	{
		m_Weights.clear();
		m_bRational = false;
		return;
	}
	//有问题
	int usize = m_pSurface->Weights()->ColLength();
	int vsize = m_pSurface->Weights()->RowLength();
	m_Weights.clear();
	m_Weights.resize(usize);
	m_bRational = true;
	for (int i = 1; i <= usize; ++i)
	{
		m_Weights[i - 1].resize(vsize);
		for (int j = 1; j <= vsize; ++j)
		{
			double weight;
			weight = m_pSurface->Weight(i, j);
			m_Weights[i - 1][j - 1] = weight;
		}
	
	}
}

void CBSplineSurface::updateUKnots()
{
	if (m_pSurface.IsNull())
	{
		return;
	}
	int size = m_pSurface->NbUKnots();
	m_UKnots.clear();
	m_UKnots.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		double knot;
		knot = m_pSurface->UKnot(i);
		m_UKnots[i - 1] = knot;
	}
}

void CBSplineSurface::updateVKnots()
{
	if (m_pSurface.IsNull())
	{
		return;
	}
	int size = m_pSurface->NbVKnots();
	m_VKnots.clear();
	m_VKnots.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		double knot;
		knot = m_pSurface->VKnot(i);
		m_VKnots[i - 1] = knot;
	}
}

void CBSplineSurface::updateUMultis()
{
	if (m_pSurface.IsNull())
	{
		return;
	}

	int size = m_pSurface->UMultiplicities().Size();
	m_UMultis.clear();
	m_UMultis.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		int multi;
		multi = m_pSurface->UMultiplicities().Value(i);
		m_UMultis[i - 1] = multi;
	}
}

void CBSplineSurface::updateVMultis()
{
	if (m_pSurface.IsNull())
	{
		return;
	}

	int size = m_pSurface->VMultiplicities().Size();
	m_VMultis.clear();
	m_VMultis.resize(size);

	for (int i = 1; i <= size; ++i)
	{
		int multi;
		multi = m_pSurface->VMultiplicities().Value(i);
		m_VMultis[i - 1] = multi;
	}
}
