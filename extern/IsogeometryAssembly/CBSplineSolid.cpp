#include "CBSplineSolid.h"
#include "Geom/DefineType.h"
#include "GeomAPI_PointsToBSplineSurface.hxx"
#include "opencascade/opencascade770-install/inc/gp_Ax1.hxx"
#include "opencascade/opencascade770-install/inc/GeomAdaptor_Surface.hxx"
#include "opencascade/opencascade770-install/inc/Extrema_ExtPS.hxx"
#include "opencascade/opencascade770-install/inc/GeomFill_FillingStyle.hxx"
#include "opencascade/opencascade770-install/inc/GeomFill_BSplineCurves.hxx"
#include"opencascade/opencascade770-install/inc/GeomFill_Pipe.hxx"
#include"opencascade/opencascade770-install/inc/GeomPlate_BuildPlateSurface.hxx"
#include"opencascade/opencascade770-install/inc/GeomPlate_MakeApprox.hxx"
#include"opencascade/opencascade770-install/inc/GeomPlate_Surface.hxx"
#include"opencascade/opencascade770-install/inc/Geom_Curve.hxx"
#include"opencascade/opencascade770-install/inc/Geom_BSplineSurface.hxx"
#include"opencascade/opencascade770-install/inc/Geom_SurfaceOfRevolution.hxx"
#include"opencascade/opencascade770-install/inc/Geom_SurfaceOfRevolution.hxx"

CBSplineSolid::CBSplineSolid()
{ 
	m_iUDegree = 0;
	m_iVDegree = 0;
	m_iWDegree = 0;

	m_ControlPoints.clear();
	m_Coefs.clear();
	m_UKnots.clear();
	m_VKnots.clear();
	m_WKnots.clear();
	m_UMultis.clear();
	m_VMultis.clear();
	m_WMultis.clear();
	m_Weights.clear();
	m_bUPeriodic = false;
	m_bVPeriodic = false;
	m_bWPeriodic = false;
	m_bRational = false;

	u_knots = nullptr;
	v_knots = nullptr;
	w_knots = nullptr;
}


CBSplineSolid::~CBSplineSolid()
{
	m_ControlPoints.clear();
	m_UKnots.clear();
	m_VKnots.clear();
	m_WKnots.clear();
	m_UMultis.clear();
	m_VMultis.clear();
	m_Weights.clear();
}

CBSplineSolid::CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints, const vector<double>& uknots, const vector<double>& vknots,
	const vector<double>& wknots, const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis, bool uperiodic, bool vperiodic , bool wperiodic )
{
	m_iUDegree = udegree;
	m_iVDegree = vdegree;
	m_iWDegree = wdegree;

	m_ControlPoints = controlPoints;

	std::vector<double> _uknots;
	std::vector<double> _vknots;
	std::vector<double> _wknots;
	std::vector<int> _uMknots;
	std::vector<int> _vMknots;
	std::vector<int> _wMknots;

	for (int i = 0; i < uknots.size(); ++i)
	{
		if (i == 0 || fabs(uknots[i] - uknots[i - 1]) > 1e-8) {
			_uknots.emplace_back(uknots[i]);
			_uMknots.emplace_back(1);
		}
		else
			_uMknots.back()++;
	}

	for (int i = 0; i < vknots.size(); ++i)
	{
		if (i == 0 || fabs(vknots[i] - vknots[i - 1]) > 1e-8) {
			_vknots.emplace_back(vknots[i]);
			_vMknots.emplace_back(1);
		}
		else
			_vMknots.back()++;
	}

	for (int i = 0; i < wknots.size(); ++i)
	{
		if (i == 0 || fabs(wknots[i] - wknots[i - 1]) > 1e-8) {
			_wknots.emplace_back(wknots[i]);
			_wMknots.emplace_back(1);
		}
		else
			_wMknots.back()++;
	}

	m_UKnots = _uknots;
	m_VKnots = _vknots;
	m_WKnots = _wknots;
	m_UMultis = _uMknots;
	m_VMultis = _vMknots;
	m_WMultis = _wMknots;
	m_bRational = false;
	m_bUPeriodic = uperiodic;
	m_bVPeriodic = vperiodic;
	m_bWPeriodic = wperiodic;
	createBSplineSolid();
}

CBSplineSolid::CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints, const vector<double>& uknots,
	const vector<double>& vknots, const vector<double>& wknots, const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis,
	const vector<vector<vector<double>>>& weights, bool uperiodic, bool vperiodic , bool wperiodic )
{
	m_iUDegree = udegree;
	m_iVDegree = vdegree;
	m_iWDegree = wdegree;
	m_iTypeCoord = 1;
	m_ControlPoints = controlPoints;
	
	std::vector<double> _uknots;
	std::vector<double> _vknots;
	std::vector<double> _wknots;
	std::vector<double> _uMknots;
	std::vector<double> _vMknots;
	std::vector<double> _wMknots;

	int pre = 0;
	for (int i = 0; i < uknots.size(); ++i)
	{
		if (i == 0 || fabs(uknots[i] - uknots[i - 1]) > 1e-8) {
			_uknots.emplace_back(uknots[i]);
			_uMknots.emplace_back(i - pre + 1);
			pre = i + 1;
		}
	}
	pre = 0;
	for (int i = 0; i < vknots.size(); ++i)
	{
		if (i == 0 || fabs(vknots[i] - vknots[i - 1]) > 1e-8) {
			_vknots.emplace_back(vknots[i]);
			_vMknots.emplace_back(i - pre + 1);
			pre = i + 1;
		}
	}
	pre = 0;
	for (int i = 0; i < wknots.size(); ++i)
	{
		if (i == 0 || fabs(wknots[i] - wknots[i - 1]) > 1e-8) {
			_wknots.emplace_back(wknots[i]);
			_wMknots.emplace_back(i - pre + 1);
			pre = i + 1;
		}
	}

	m_UKnots = uknots;
	m_VKnots = vknots;
	m_WKnots = wknots;



	m_UMultis = umultis;
	m_VMultis = vmultis;
	m_Weights = weights;
	m_bRational = true;
	m_bUPeriodic = uperiodic;
	m_bVPeriodic = vperiodic;
	m_bVPeriodic = wperiodic;

	//createBSplineSolid();
}

int CBSplineSolid::getSizeUControlPoints() const
{
	if (m_ControlPoints.empty())
	{
		return -1;
	}
	size_t size = m_ControlPoints.size();
	return (int)size;
}

int CBSplineSolid::getSizeVControlPoints() const
{
	if (m_ControlPoints.empty())
	{
		return -1;
	}
	size_t size = m_ControlPoints[0].size();
	return (int)size;
}

int CBSplineSolid::getSizeWControlPoints() const
{
	if (m_ControlPoints.empty())
	{
		return -1;
	}
	size_t size = m_ControlPoints[0][0].size();
	return (int)size;
}

int CBSplineSolid::getSizeControlPoints() const
{
	return m_Coefs.size();
}

bool CBSplineSolid::setUDegree(int udegree)
{
		m_iUDegree = udegree;
		return true;

}

bool CBSplineSolid::setVDegree(int vdegree)
{
		m_iVDegree = vdegree;
		return true;

}

bool CBSplineSolid::setWDegree(int wdegree)
{
	    m_iWDegree = wdegree;
		return true;

}
int CBSplineSolid::getWDegree() const
{
	return m_iWDegree;
}
bool CBSplineSolid::getControlPointsByIndex(const std::vector<int>& localIndex, Eigen::MatrixXd& res)
{
	res.resize(localIndex.size(), 3);
	for (size_t i = 0; i < localIndex.size(); ++i) {
		res(i, 0) = m_Coefs[localIndex[i]].getX();
		res(i, 1) = m_Coefs[localIndex[i]].getY();
		res(i, 2) = m_Coefs[localIndex[i]].getZ();
	}
	return true;
}
bool CBSplineSolid::getControlPointsByCoef(const std::vector<int>& localIndex, const Eigen::MatrixXd& u_coef, Eigen::MatrixXd& res)
{
	res.resize(localIndex.size(), 3);
	if (u_coef.rows() != m_Coefs.size()) {
		cout << "error input coef" << endl;
		exit(1);
	}
	for (size_t i = 0; i < localIndex.size(); ++i) {
		res(i, 0) = u_coef(localIndex[i], 0);
		res(i, 1) = u_coef(localIndex[i], 1);
		res(i, 2) = u_coef(localIndex[i], 2);
	}
	return true;
}
int CBSplineSolid::gethandle()const {
	return SolidHandle;
}
void CBSplineSolid::sethandle(int handle) {
	SolidHandle = handle;
}
bool CBSplineSolid::setControlPoints(const vector<vector<vector<CPoint>>>& controlPoints)
{
		m_ControlPoints = controlPoints;
		return true;
}

bool CBSplineSolid::setControlPoints(const vector<vector<vector<CPoint>>>& controlPoints, const vector<vector<vector<double>>>& weights)
{
	m_ControlPoints = controlPoints;
	m_Weights = weights;
	return true;
	
}

bool CBSplineSolid::setSingleControlPoint(int& uindex, int& vindex,int& windex, const CPoint& controlPoint)
{
	size_t sizeu = m_ControlPoints.size();
	size_t sizev = m_ControlPoints[0].size();
	size_t sizew = m_ControlPoints[0][0].size();
	if ( (uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev)|| (windex < 1 || windex > sizew))
	{
		return false;
	}
	m_ControlPoints[uindex - 1][vindex - 1][windex] = controlPoint;
	gp_Pnt p;
	controlPoint.getPnt(p);
	return true;
}

bool CBSplineSolid::setSingleControlPoint(int& uindex, int& vindex,int& windex, const CPoint& controlPoint, double& weight)
{
	size_t sizeu = m_ControlPoints.size();
	size_t sizev = m_ControlPoints[0].size();
	size_t sizew = m_ControlPoints[0][0].size();
	if ((uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev) || (windex < 1 || windex > sizew))
	{
		return false;
	}
	m_ControlPoints[uindex - 1][vindex - 1][windex - 1] = controlPoint;
	m_Weights[uindex - 1][vindex - 1][windex - 1] = weight;
	gp_Pnt p;
	controlPoint.getPnt(p);
	return true;
}

bool CBSplineSolid::setWeights(const vector<vector<vector<double>>>& weights)
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
	size_t sizew = weights[0][0].size();
	if (sizew == 0)
	{
		return false;
	}
	for (int i = 0; i < sizeu; ++i)
	{
		for (int j = 0; j < sizev; ++j)
		{
			for (int k = 0; k < sizew; k++)
			{
				if (weights[i][j][k] < 1e-10)
				{
					return false;
				}
			}
		}
	}
		m_Weights = weights;
		m_bRational = true;
		return true;
}

bool CBSplineSolid::setSingleWeight(int uindex, int vindex,int windex, double weight)
{
	size_t sizeu = m_Weights.size();
	size_t sizev = m_Weights[0].size();
	size_t sizew = m_Weights[0][0].size();
	if ( (uindex < 1 || uindex > sizeu) || (vindex < 1 || vindex > sizev) || (windex < 1 || windex > sizew) && !m_bRational)
	{
		return false;
	}

	if (weight < 1e-10)
	{
		return false;
	}

	m_Weights[uindex - 1][vindex - 1][windex - 1] = weight;
	return true;
}

const vector<double>& CBSplineSolid::getUKnots() const
{
	return m_UKnots;
}

const vector<double>& CBSplineSolid::getVKnots() const
{
	return m_VKnots;
}

const vector<double>& CBSplineSolid::getWKnots() const
{
	return m_WKnots;
}

vector<double> CBSplineSolid::getUMKnots() const
{
	vector<double> res = { m_UKnots[0] };
	for (auto& it : m_UKnots) res.push_back(it);
	res.push_back(m_UKnots.back());
	return res;
}

vector<double> CBSplineSolid::getVMKnots() const
{
	vector<double> res = { m_VKnots[0] };
	for (auto& it : m_VKnots) res.push_back(it);
	res.push_back(m_VKnots.back());
	return res;
}

vector<double> CBSplineSolid::getWMKnots() const
{
	vector<double> res = { m_WKnots[0] };
	for (auto& it : m_WKnots) res.push_back(it);
	res.push_back(m_WKnots.back());
	return res;
}


vector<vector<vector<double>>>CBSplineSolid::getWeights() const {
	return m_Weights;
}


vector<vector<vector<CPoint>>> CBSplineSolid::getControlPoints() const
{
	return m_ControlPoints;
}

const vector<CPoint>& CBSplineSolid::coefs() const
{
	return m_Coefs;
}






CPoint CBSplineSolid::getSingleControlPoint(int uindex, int vindex,int windex)
{
	size_t size = m_ControlPoints.size();
	if ((uindex < 1 || uindex > size))
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

	size_t wsize = m_ControlPoints[0][0].size();
	if ((windex < 1 || windex > wsize))
	{
		CPoint tmp(INT_MAX, INT_MAX, INT_MAX, -1, 1);
		return tmp;
	}
	return m_ControlPoints[uindex - 1][vindex - 1][windex - 1];
}

double CBSplineSolid::getSingleWeight(int uindex, int vindex,int windex)
{
	size_t size = m_ControlPoints.size();

	if ((uindex < 1 || uindex > size))
	{
		return -1;
	}

	size_t vsize = m_ControlPoints[0].size();
	if ((vindex < 1 || vindex > vsize))
	{
		return -1;
	}

	size_t wsize = m_ControlPoints[0][0].size();
	if ((windex < 1 || windex > wsize))
	{
		return -1;
	}
	return m_Weights[uindex - 1][vindex - 1][windex - 1];
}



bool CBSplineSolid::getRational() const
{
	return m_bRational;
}

int CBSplineSolid::getUDegree() const
{
	return m_iUDegree;
}

int CBSplineSolid::getVDegree() const
{
	return m_iVDegree;
}


CPoint CBSplineSolid::getPointAtParam(double u, double v,double w)
{
	CPoint pnt;
	getPointAtParam(u, v,w, pnt);
	//std::cout << "getPointAtParam" << std::endl;
	return pnt;
}

bool CBSplineSolid::getPointAtParam(double u, double v,double w ,CPoint& point)

{
	//std::cout << "eval begin" << std::endl;
	eval(u, v, w, basisValue);
	//std::cout << "eval finish" << std::endl;
	getLocation(u, v, w, P);
	//std::cout << "getLocation" << std::endl;
	Location = basisValue.transpose() * P;
	point = CPoint(Location(0, 0), Location(0, 1), Location(0, 2));
	return true;
}

void CBSplineSolid::eval(double u, double v, double w, Eigen::VectorXd& basisValue)
{
	int tol = 1e-15;
	int u_ctrlPtsNum = u_size - 1 - m_iUDegree - 1;
	int v_ctrlPtsNum = v_size - 1 - m_iVDegree - 1;
	int w_ctrlPtsNum = w_size - 1 - m_iWDegree - 1;

	int uspan = nrbFindSpan(u_ctrlPtsNum, m_iUDegree, u, u_knots);
	int vspan = nrbFindSpan(v_ctrlPtsNum, m_iVDegree, v, v_knots);
	int wspan = nrbFindSpan(w_ctrlPtsNum, m_iWDegree, w, w_knots);

	if (fabs(u - u_knots[u_ctrlPtsNum - 1]) < tol)
		u = u_knots[u_size - 1] - tol;

	if (fabs(v - v_knots[v_ctrlPtsNum - 1]) < tol)
		v = v_knots[v_size - 1] - tol;

	if (fabs(w - w_knots[w_ctrlPtsNum - 1]) < tol)
		w = w_knots[w_size - 1] - tol;

	double* Nu = (double*)malloc(sizeof(double) * (m_iUDegree + 1));
	double* Nv = (double*)malloc(sizeof(double) * (m_iVDegree + 1));
	double* Nw = (double*)malloc(sizeof(double) * (m_iWDegree + 1));

	nrbAllBsplineBasisFuns(uspan, u, m_iUDegree, u_knots, Nu);
	nrbAllBsplineBasisFuns(vspan, v, m_iVDegree, v_knots, Nv);
	nrbAllBsplineBasisFuns(wspan, w, m_iWDegree, w_knots, Nw);

	int vind, uind, wind;
	double wgt = 1, we = 0.0, fac, dwdxi = 0.0, dwdeta = 0.0, dwdgemma = 0.0;
	vector<double> res;
	for (int k = 0; k <= m_iWDegree; k++) {
		wind = wspan - m_iWDegree + k;
		for (int j = 0; j <= m_iVDegree; j++)
		{
			vind = vspan - m_iVDegree + j;
			for (int i = 0; i <= m_iUDegree; i++)
			{
				uind = uspan - m_iUDegree + i;
				we += Nu[i] * Nv[j] * Nw[k] * wgt;
			}
		}
	}

	basisValue.resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1));
	int k = 0;
	for (int t = 0; t <= m_iWDegree; t++) {
		wind = wspan - m_iWDegree + t;
		for (int j = 0; j <= m_iVDegree; j++)
		{
			vind = vspan - m_iVDegree + j;
			for (int i = 0; i <= m_iUDegree; i++)
			{
				uind = uspan - m_iUDegree + i;
				fac = wgt / (we * we);
				basisValue[k] = Nu[i] * Nv[j] * Nw[t] * fac * we;
				++k;
			}
		}
	}
	free(Nu);
	free(Nv);
	free(Nw);
}

void CBSplineSolid::getLocation(double u, double v, double w, Eigen::MatrixXd& P)
{
	int u_interval = 0, v_interval = 0, w_interval = 0;
	for (int i = 1; i < m_UKnots.size(); ++i) {
		if (m_UKnots[i] < u)
			u_interval = i;
	}
	for (int i = 1; i < m_VKnots.size(); ++i) {
		if (m_VKnots[i] < v)
			v_interval = i;
	}
	for (int i = 1; i < m_WKnots.size(); ++i) {
		if (m_WKnots[i] < w)
			w_interval = i;
	}
	P.resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1), 3);
	int count = 0;
	for (int k = w_interval; k < w_interval + m_iWDegree + 1; ++k) {
		for (int j = v_interval; j < v_interval + m_iVDegree + 1; ++j) {
			for (int i = u_interval; i < u_interval + m_iUDegree + 1; ++i) {
				P(count, 0) = m_ControlPoints[i][j][k].getX();
				P(count, 1) = m_ControlPoints[i][j][k].getY();
				P(count, 2) = m_ControlPoints[i][j][k].getZ();
				count++;
			}
		}
	}
}

bool CBSplineSolid::getParamAtPoint(const CPoint& point, double& u,double &v, double& w)
{
	return true;
}

bool CBSplineSolid::getUVMinMax(double& umin, double& umax, double& vmin, double& vmax, double& wmin, double& wmax)
{
	umin = 0;
	umax = 1;
	vmin = 0;
	vmax = 1;
	wmin = 0;
	wmax = 1;
	return true;
}


void CBSplineSolid::updateUKnots()
{
	m_UKnots_nonredundant = m_UKnots;
	vector<double > tmp;
	tmp.push_back(0);
	for (int i = 0; i < m_UKnots.size(); ++i)
	{
		for (int j = 0; j < m_UMultis[i]; j++)
		{
			tmp.push_back(m_UKnots[i]);
		}
	}
	m_UKnots = tmp;
	u_size = m_UKnots.size();
}

void CBSplineSolid::updateVKnots()
{
	m_VKnots_nonredundant = m_VKnots;
	vector<double > tmp;
	tmp.push_back(0);
	for (int i = 0; i < m_VKnots.size(); ++i)
	{
		for (int j = 0; j < m_VMultis[i]; j++)
		{
			tmp.push_back(m_VKnots[i]);
		}
	}
	m_VKnots = tmp;
	v_size = m_VKnots.size();
}


void CBSplineSolid::updateWKnots()
{
	m_WKnots_nonredundant = m_WKnots;
	vector<double > tmp;
	tmp.push_back(0);
	for (int i = 0; i < m_WKnots.size(); ++i)
	{
		for (int j = 0; j < m_WMultis[i]; j++)
		{
			tmp.push_back(m_WKnots[i]);
		}
	}
	m_WKnots = tmp;
	w_size = m_WKnots.size();
}

bool CBSplineSolid::getJacobi(double u, double v, double w, std::vector<int>& localIndex, Eigen::MatrixXd& jac)
{
	vector<Eigen::VectorXd> basisValue;
	eval_2Ders(u, v, w, basisValue);
	int len = localIndex.size();
	if ((int)basisValue[1].rows() != 3 * len) return false;

	const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue[1].data(), 3, len);

	Eigen::MatrixXd X;
	getControlPointsByIndex(localIndex, X);
	jac = Nder * X;
	return true;
}

bool CBSplineSolid::getJacobi(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, Eigen::MatrixXd& jac)
{
	int len = localIndex.size();
	if ((int)basisValue.rows() != 3 * len) return false;

	const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue.data(), 3, len);
	Eigen::MatrixXd X;
	getControlPointsByIndex(localIndex, X);
	jac = Nder * X;
	return false;
}

bool CBSplineSolid::getJacobiU(double u, double v, double w, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement, Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU)
{
	vector<Eigen::VectorXd> basisValue;
	eval_2Ders(u, v, w, basisValue);
	int len = localIndex.size();
	if ((int)basisValue[1].rows() != 3 * len) return false;
	//cout << basisValue[1] << "\n\n\n";
	const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue[1].data(), 3, len);

	Eigen::MatrixXd X, U;
	getControlPointsByIndex(localIndex, X);
	getControlPointsByCoef(localIndex, displacement, U);
	jac = Nder * X;
	jacU = Nder * U;
	return true;
}

bool CBSplineSolid::getJacobiU(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement, Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU)
{
	int len = localIndex.size();
	if ((int)basisValue.rows() != 3 * len) {
		cout << "error" << endl;
		return false;
	}
	//cout << basisValue << "\n\n\n";
	const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue.data(), 3, len);

	Eigen::MatrixXd X, U;
	getControlPointsByIndex(localIndex, X);
	getControlPointsByCoef(localIndex, displacement, U);
	jac = Nder * X;
	jacU = Nder * U;
	return true;
}

bool CBSplineSolid::updateCoefs(const Eigen::MatrixXd& coefs)
{
	if (coefs.rows() != m_Coefs.size()) return false;

	for (size_t i = 0; i < coefs.rows(); ++i) {
		double x = m_Coefs[i].getX() + coefs(i, 0);
		double y = m_Coefs[i].getY() + coefs(i, 1);
		double z = m_Coefs[i].getZ() + coefs(i, 2);
		m_Coefs[i].setXYZ(x, y, z);
	}

	int cnt = 0;
	int u_p = u_size - m_iUDegree - 1;
	int v_p = v_size - m_iVDegree - 1;
	int w_p = w_size - m_iWDegree - 1;
	for (int k = 0; k < w_p; ++k)
	{
		for (int j = 0; j < v_p; ++j)
		{
			for (int i = 0; i < u_p; ++i)
			{
				m_ControlPoints[i][j][k] = m_Coefs[cnt++];
			}
		}
	}


	return true;
}

void CBSplineSolid::eval_AllDers(std::vector<double> u, std::vector<double> v, std::vector<double> w, vector<Eigen::MatrixXd>& basisValue)
{
	int tol = 1e-15;
	int u_ctrlPtsNum = u_size - 1 - m_iUDegree - 1;
	int v_ctrlPtsNum = v_size - 1 - m_iVDegree - 1;
	int w_ctrlPtsNum = w_size - 1 - m_iWDegree - 1;
	basisValue.resize(2);
	basisValue[0].resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1), u.size() * v.size() * w.size());
	basisValue[1].resize(3 * (m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1), u.size() * v.size() * w.size());
	vector<Eigen::VectorXd> singleValue;
	for (int k = 0; k < w.size(); ++k) {
		for (int j = 0; j < v.size(); ++j) {
			for (int i = 0; i < u.size(); ++i) {
				eval_2Ders(u[i], v[j], w[k], singleValue);
				basisValue[0].col(k * u.size() * v.size() + j * u.size() + i) = singleValue[0];
				basisValue[1].col(k * u.size() * v.size() + j * u.size() + i) = singleValue[1];
			}
		}
	}
}
void CBSplineSolid::eval_into(std::vector<double> u, std::vector<double> v, std::vector<double> w, Eigen::MatrixXd& basisValue)
{
	int tol = 1e-15;
	basisValue.resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1), u.size() * v.size() * w.size());
	Eigen::VectorXd singleValue;
	for (int k = 0; k < w.size(); ++k) {
		for (int j = 0; j < v.size(); ++j) {
			for (int i = 0; i < u.size(); ++i) {
				eval(u[i], v[j], w[k], singleValue);
				basisValue.col(k * v.size() * u.size() + j * u.size() + i) = singleValue;
			}
		}
	}

}



void CBSplineSolid::eval_2Ders(double u, double v, double w, vector<Eigen::VectorXd>& basisValue)
{
	int tol = 1e-15;
	int u_ctrlPtsNum = u_size - 1 - m_iUDegree - 1;
	int v_ctrlPtsNum = v_size - 1 - m_iVDegree - 1;
	int w_ctrlPtsNum = w_size - 1 - m_iWDegree - 1;
	basisValue.resize(2);

	int uspan = nrbFindSpan(u_ctrlPtsNum, m_iUDegree, u, u_knots);
	int vspan = nrbFindSpan(v_ctrlPtsNum, m_iVDegree, v, v_knots);
	int wspan = nrbFindSpan(w_ctrlPtsNum, m_iWDegree, w, w_knots);

	if (fabs(u - u_knots[u_ctrlPtsNum - 1]) < tol)
		u = u_knots[u_size - 1] - tol;

	if (fabs(v - v_knots[v_ctrlPtsNum - 1]) < tol)
		v = v_knots[v_size - 1] - tol;

	if (fabs(w - w_knots[w_ctrlPtsNum - 1]) < tol)
		w = w_knots[w_size - 1] - tol;

	double* Nu = (double*)malloc(sizeof(double) * (m_iUDegree + 1));
	double* Nv = (double*)malloc(sizeof(double) * (m_iVDegree + 1));
	double* Nw = (double*)malloc(sizeof(double) * (m_iWDegree + 1));
	double** Nu_ders = init2DArray(u_ctrlPtsNum + 1, m_iUDegree + 1);
	double** Nv_ders = init2DArray(v_ctrlPtsNum + 1, m_iVDegree + 1);
	double** Nw_ders = init2DArray(w_ctrlPtsNum + 1, m_iWDegree + 1);

	nrbAllBsplineBasisFuns(uspan, u, m_iUDegree, u_knots, Nu);
	nrbAllBsplineBasisFuns(vspan, v, m_iVDegree, v_knots, Nv);
	nrbAllBsplineBasisFuns(wspan, w, m_iWDegree, w_knots, Nw);
	nrbAllBsplineBasisFunsDers(uspan, u, m_iUDegree, u_knots, u_ctrlPtsNum, Nu_ders);
	nrbAllBsplineBasisFunsDers(vspan, v, m_iVDegree, v_knots, v_ctrlPtsNum, Nv_ders);
	nrbAllBsplineBasisFunsDers(wspan, w, m_iWDegree, w_knots, w_ctrlPtsNum, Nw_ders);

	int vind, uind, wind;
	double wgt = 1, we = 0.0, fac, dwdxi = 0.0, dwdeta = 0.0, dwdgemma = 0.0;
	vector<double> res;
	for (int k = 0; k <= m_iWDegree; k++) {
		wind = wspan - m_iWDegree + k;
		for (int j = 0; j <= m_iVDegree; j++)
		{
			vind = vspan - m_iVDegree + j;
			for (int i = 0; i <= m_iUDegree; i++)
			{
				uind = uspan - m_iUDegree + i;
				we += Nu[i] * Nv[j] * Nw[k] * wgt;
				dwdxi += Nu_ders[1][i] * Nv[j] * Nw[k] * wgt;
				dwdeta += Nu[i] * Nv_ders[1][j] * Nw[k] * wgt;
				dwdgemma += Nu[i] * Nv[j] * Nw_ders[1][k] * wgt;
			}
		}
	}
	basisValue[0].resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1));
	basisValue[1].resize(3 * (m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree + 1));
	int k = 0;
	for (int t = 0; t <= m_iWDegree; t++) {
		wind = wspan - m_iWDegree + t;
		for (int j = 0; j <= m_iVDegree; j++)
		{
			vind = vspan - m_iVDegree + j;
			for (int i = 0; i <= m_iUDegree; i++)
			{
				uind = uspan - m_iUDegree + i;
				fac = wgt / (we * we);
				//cout << Nu[i] << " " << Nv[j] << " " << Nw[t] << "\n";
				basisValue[0][k] = Nu[i] * Nv[j] * Nw[t] * fac * we;
				basisValue[1][3 * k] = (Nu_ders[1][i] * Nv[j] * Nw[t] * we - Nu[i] * Nv[j] * Nw[t] * dwdxi) * fac;
				basisValue[1][3 * k + 1] = (Nu[i] * Nv_ders[1][j] * Nw[t] * we - Nu[i] * Nv[j] * Nw[t] * dwdeta) * fac;
				basisValue[1][3 * k + 2] = (Nu[i] * Nv[j] * Nw_ders[1][t] * we - Nu[i] * Nv[j] * Nw[t] * dwdgemma) * fac;
				//cout << basisValue[1][3 * k] << " " << basisValue[1][3 * k + 1] << " " << basisValue[1][3 * k + 2] << "\n";
				++k;
			}
		}
	}
	free(Nu);
	free(Nv);
	free(Nw);
	free2Darray(Nu_ders, (u_ctrlPtsNum + 1));
	free2Darray(Nv_ders, (v_ctrlPtsNum + 1));
	free2Darray(Nw_ders, (w_ctrlPtsNum + 1));
}
void CBSplineSolid::active(int u_interval, int v_interval, int w_interval, std::vector<int>& localIndex) const
{
	int count = 0;
	for (int k = 0; k < m_iWDegree + 1; ++k) {
		for (int j = 0; j < m_iVDegree + 1; ++j) {
			for (int i = 0; i < m_iUDegree + 1; ++i)
				localIndex[count++] = (w_interval + k) * (v_size - m_iVDegree - 1) * (u_size - m_iUDegree - 1) +
				(v_interval + j) * (u_size - m_iUDegree - 1) + (u_interval + i);
		}
	}
}


void CBSplineSolid::activeFromParam(double u, double v, double w, std::vector<int>& localIndex) const
{
	int u_interval = 0, v_interval = 0, w_interval = 0;
	for (int i = 1; i < m_UKnots.size(); ++i) {
		if (m_UKnots[i] < u)
			u_interval = i;
	}
	for (int i = 1; i < m_VKnots.size(); ++i) {
		if (m_VKnots[i] < v)
			v_interval = i;
	}
	for (int i = 1; i < m_WKnots.size(); ++i) {
		if (m_WKnots[i] < w)
			w_interval = i;
	}
	active(u_interval, v_interval, w_interval, localIndex);
}

CPoint CBSplineSolid::getPointFromParam(double u, double v, double w)
{
	Eigen::VectorXd basisValue;
	Eigen::MatrixXd P;
	eval(u, v, w, basisValue);
	getLocation(u, v, w, P);
	auto it = basisValue.transpose() * P;
	CPoint res{ it(0,0),it(0,1),it(0,2) };
	return res;
}

bool CBSplineSolid::createBSplineSolid()
{
	u_size = m_ControlPoints.size() + m_iUDegree + 1;
	v_size = m_ControlPoints[0].size() + m_iVDegree + 1;
	w_size = m_ControlPoints[0][0].size() + m_iWDegree + 1;

	u_knots = (double*)malloc(sizeof(double) * u_size);
	v_knots = (double*)malloc(sizeof(double) * v_size);
	w_knots = (double*)malloc(sizeof(double) * w_size);
	int count = 0;
	for (int i = 0; i < m_UKnots.size(); ++i) {
		for (int j = 0; j < m_UMultis[i]; ++j)
			u_knots[count++] = m_UKnots[i];
	}

	count = 0;
	for (int i = 0; i < m_VKnots.size(); ++i) {
		for (int j = 0; j < m_VMultis[i]; ++j)
			v_knots[count++] = m_VKnots[i];
	}
	count = 0;
	for (int i = 0; i < m_WKnots.size(); ++i) {
		for (int j = 0; j < m_WMultis[i]; ++j)
			w_knots[count++] = m_WKnots[i];
	}

	if (m_iUDegree < 1 || m_iUDegree > 25 || m_iVDegree < 1 || m_iVDegree > 25)
	{
		cerr << "BSplineSurface Contruct Error! Degree Out of Range" << endl;
		return false;
	}

	int cnt = 0;
	int u_p = u_size - m_iUDegree - 1;
	int v_p = v_size - m_iVDegree - 1;
	int w_p = w_size - m_iWDegree - 1;
	m_Coefs.resize(u_p * v_p * w_p);
	for (int k = 0; k < w_p; ++k)
	{
		for (int j = 0; j < v_p; ++j)
		{
			for (int i = 0; i < u_p; ++i)
			{
				m_Coefs[cnt++] = m_ControlPoints[i][j][k];
			}
		}
	}


	facePointIndex.resize(6);
	int usize = m_ControlPoints.size();
	int vsize = m_ControlPoints[0].size();
	int wsize = m_ControlPoints[0][0].size();
	auto& coefs = m_ControlPoints;
	cnt = 0;
	for (int w = 0; w < wsize; ++w)
	{
		for (int j = 0; j < vsize; ++j)
		{
			for (int i = 0; i < usize; ++i)
			{
				if (i == 0) facePointIndex[0].insert(cnt);
				if (i == usize - 1) facePointIndex[1].insert(cnt);
				if (j == 0) facePointIndex[2].insert(cnt);
				if (j == vsize - 1) facePointIndex[3].insert(cnt);
				if (w == 0) facePointIndex[4].insert(cnt);
				if (w == wsize - 1) facePointIndex[5].insert(cnt);
				cnt++;
			}
		}
	}
	return true;
}
