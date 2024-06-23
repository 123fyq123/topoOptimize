#include "CBSplineSolid.h"

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

CBSplineSolid::CBSplineSolid(const CBSplineSolid& solid_)
{
    m_iUDegree = solid_.m_iUDegree;
    m_iVDegree = solid_.m_iVDegree;
    m_iWDegree = solid_.m_iWDegree;

    m_ControlPoints = solid_.m_ControlPoints;
    m_Coefs = solid_.m_Coefs;
    m_UKnots = solid_.m_UKnots;
    m_VKnots = solid_.m_VKnots;
    m_WKnots = solid_.m_WKnots;
    m_UMultis = solid_.m_UMultis;
    m_VMultis = solid_.m_VMultis;
    m_WMultis = solid_.m_WMultis;
    m_bRational = solid_.m_bRational;
    m_bUPeriodic = solid_.m_bUPeriodic;
    m_bVPeriodic = solid_.m_bVPeriodic;
    m_bWPeriodic = solid_.m_bWPeriodic;
    m_Weights = solid_.m_Weights;

    createBSplineSolid();
}

CBSplineSolid::CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints, const vector<double>& uknots, const vector<double>& vknots, const vector<double>& wknots,
                             const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis, bool uperiodic, bool vperiodic, bool wperiodic)
{
    m_iUDegree = udegree;
    m_iVDegree = vdegree;
    m_iWDegree = wdegree;

    m_ControlPoints = controlPoints;
    m_UKnots = uknots;
    m_VKnots = vknots;
    m_WKnots = wknots;
    m_UMultis = umultis;
    m_VMultis = vmultis;
    m_WMultis = wmultis;
    m_bRational = false;
    m_bUPeriodic = uperiodic;
    m_bVPeriodic = vperiodic;
    m_bWPeriodic = wperiodic;
    createBSplineSolid();
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

vector<double> CBSplineSolid::getUAllKnots() const
{
    vector<double> res;
    for (int i = 0; i < m_UKnots.size(); ++i) {
        for (int j = 0; j < m_UMultis[i]; ++j)
            res.push_back(m_UKnots[i]);
    }
    return res;
}

vector<double> CBSplineSolid::getVAllKnots() const
{
    vector<double> res;
    for (int i = 0; i < m_VKnots.size(); ++i) {
        for (int j = 0; j < m_VMultis[i]; ++j)
            res.push_back(m_VKnots[i]);
    }
    return res;
}

vector<double> CBSplineSolid::getWAllKnots() const
{
    vector<double> res;
    for (int i = 0; i < m_WKnots.size(); ++i) {
        for (int j = 0; j < m_WMultis[i]; ++j)
            res.push_back(m_WKnots[i]);
    }
    return res;
}

int CBSplineSolid::getUDegree() const
{
    return m_iUDegree;
}

int CBSplineSolid::getVDegree() const
{
    return m_iVDegree;
}

int CBSplineSolid::getWDegree() const
{
    return m_iWDegree;
}

bool CBSplineSolid::getControlPointsByIndex(const std::vector<int>& localIndex, Eigen::MatrixXd& res) const
{
    res.resize(localIndex.size(), 3);
    for (size_t i = 0; i < localIndex.size(); ++i) {
        res(i,0) = m_Coefs[localIndex[i]].getX();
        res(i,1) = m_Coefs[localIndex[i]].getY();
        res(i,2) = m_Coefs[localIndex[i]].getZ();
    }

    //int m = m_ControlPoints.size(), n = m_ControlPoints[0].size();
    //res.resize(3); res[0].resize(localIndex.size());
    //res[1].resize(localIndex.size()); res[2].resize(localIndex.size());
    //int count = 0;
    //for (auto index : localIndex) {
    //	int z = index / (m * n); index %= (m * n);
    //	int y = index / m, x = index % m;
    //	res[0](count) = m_ControlPoints[x][y][z].getX();
    //	res[1](count) = m_ControlPoints[x][y][z].getY();
    //	res[2](count++) = m_ControlPoints[x][y][z].getZ();
    //}
    return true;
}

bool CBSplineSolid::getControlPointsByCoef(const std::vector<int>& localIndex, const Eigen::MatrixXd& u_coef, Eigen::MatrixXd& res) const
{
    res.resize(localIndex.size(), 3);
    if (u_coef.rows() != m_Coefs.size()) {
        cout << "error input coef" << endl;
        exit(1);
    }
    for (size_t i = 0; i < localIndex.size(); ++i) {
        res(i, 0) = u_coef(localIndex[i],0);
        res(i, 1) = u_coef(localIndex[i],1);
        res(i, 2) = u_coef(localIndex[i],2);
    }
    return true;
}

vector<vector<vector<CPoint>>> CBSplineSolid::getControlPoints() const
{
    return m_ControlPoints;
}

const vector<CPoint>& CBSplineSolid::coefs() const
{
    return m_Coefs;
}

bool CBSplineSolid::getJacobi(double u, double v, double w, std::vector<int> &localIndex, Eigen::MatrixXd& jac)
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

bool CBSplineSolid::getJacobi(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, Eigen::MatrixXd& jac) const
{
    int len = localIndex.size();
    if ((int)basisValue.rows() != 3 * len) return false;

    const Eigen::Map<const Eigen::MatrixXd> Nder(basisValue.data(), 3, len);
    Eigen::MatrixXd X;
    getControlPointsByIndex(localIndex, X);
    jac = Nder * X;
    return false;
}

bool CBSplineSolid::getJacobiU(double u, double v, double w, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement,
                               Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU) const
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

    //cout << jac<<endl;
    jacU = Nder * U;
    return true;
}

bool CBSplineSolid::getJacobiU(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement, Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU) const
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

bool CBSplineSolid::setCoefs(const Eigen::MatrixXd& coefs)
{
    if (coefs.rows() != m_Coefs.size()) return false;

    for (size_t i = 0; i < coefs.rows(); ++i) {
        double x = coefs(i, 0);
        double y = coefs(i, 1);
        double z = coefs(i, 2);
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

void CBSplineSolid::eval_AllDers(std::vector<double> u, std::vector<double> v, std::vector<double> w, vector<Eigen::MatrixXd>& basisValue) const
{
    int tol = 1e-15;
    int u_ctrlPtsNum = u_size - 1 - m_iUDegree - 1;
    int v_ctrlPtsNum = v_size - 1 - m_iVDegree - 1;
    int w_ctrlPtsNum = w_size - 1 - m_iWDegree - 1;
    basisValue.resize(2);
    basisValue[0].resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree+1), u.size() * v.size() * w.size());
    basisValue[1].resize(3 * (m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree+1), u.size() * v.size() * w.size());
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
    basisValue.resize((m_iUDegree + 1) * (m_iVDegree + 1) * (m_iWDegree+1), u.size() * v.size()*w.size());
    Eigen::VectorXd singleValue;
    for (int k = 0; k < w.size(); ++k) {
        for (int j = 0; j < v.size(); ++j) {
            for (int i = 0; i < u.size(); ++i) {
                eval(u[i], v[j], w[k], singleValue);
                basisValue.col(k * v.size()*u.size() + j * u.size() + i) = singleValue;
            }
        }
    }

}



void CBSplineSolid::eval_2Ders(double u, double v, double w, vector<Eigen::VectorXd>& basisValue) const
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

    assert(Nu[0] = Nu_ders[0][0]);
    assert(Nu[1] = Nu_ders[0][1]);
    assert(Nu[2] = Nu_ders[0][2]);
    assert(Nu[3] = Nu_ders[0][3]);

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
                we += Nu[i] * Nv[j]* Nw[k] * wgt;
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

void CBSplineSolid::active(int u_interval, int v_interval, int w_interval, std::vector<int>& localIndex) const
{
    int count = 0;
    for (int k = 0; k < m_iWDegree + 1; ++k) {
        for (int j = 0; j < m_iVDegree + 1; ++j) {
            for (int i = 0; i < m_iUDegree + 1; ++i)
                localIndex[count++] = (w_interval + k)*(v_size - m_iVDegree-1)* (u_size - m_iUDegree - 1) +
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

CPoint CBSplineSolid::getPointFromParam(double u, double v, double w)
{
    Eigen::VectorXd basisValue;
    Eigen::MatrixXd P;
    eval(u, v, w, basisValue);
    //cout << basisValue << endl;
    getLocation(u, v, w, P);
    //cout << "\n";
    //cout << P << endl;
    auto it = basisValue.transpose()* P;
    //cout << it << endl;
    CPoint res{ it(0,0),it(0,1),it(0,2) };
    return res;
}

Eigen::MatrixXd CBSplineSolid::eval(const Eigen::MatrixXd& u)
{
    Eigen::MatrixXd points(3, u.cols());
    for (unsigned i = 0; i < u.cols(); ++i)
    {
        auto point = getPointFromParam(u(0, i), u(1, i), u(2, i));
        points(0, i) = point.getX();
        points(1, i) = point.getY();
        points(2, i) = point.getZ();
    }
    return points;
}
int CBSplineSolid::getUSize() const
{
    return u_size;
}

int CBSplineSolid::getVSize() const
{
    return v_size;
}

int CBSplineSolid::getWSize() const
{
    return w_size;
}
const vector<CPoint>& CBSplineSolid::getCoefs() const
{
    return m_Coefs;
}
