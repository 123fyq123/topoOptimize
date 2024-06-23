#pragma once
#include "Geom/CBSplineSurface.h"
#include<unordered_set>
#include <Eigen/Core>

class CBSplineSolid {


    bool m_bRational;
    bool m_bUPeriodic;
    bool m_bVPeriodic;
    bool m_bWPeriodic;
    vector<vector<vector<CPoint>>> m_ControlPoints;
    vector<CPoint> m_Coefs;
    vector<double> m_UKnots;
    vector<double> m_VKnots;
    vector<double> m_WKnots;
    vector<vector<double>> m_Weights;
    vector<int>    m_UMultis;
    vector<int>    m_VMultis;
    vector<int>    m_WMultis;

    double* u_knots;
    double* v_knots;
    double* w_knots;
    int u_size = 0;
    int v_size = 0;
    int w_size = 0;

    int m_iUDegree;
    int m_iVDegree;
    int m_iWDegree;

    std::vector<unordered_set<int>> facePointIndex;

public:
    CBSplineSolid();


    CBSplineSolid(const CBSplineSolid& solid_);


    CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints, const vector<double>& uknots,
                  const vector<double>& vknots, const vector<double>& wknots, const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis,
                  bool uperiodic = false, bool vperiodic = false, bool wperiodic = false);

    bool createBSplineSolid();

    int getSizeUControlPoints() const;

    int getSizeVControlPoints() const;

    int getSizeWControlPoints() const;

    int getSizeControlPoints() const;

    const vector<double>& getUKnots() const;

    const vector<double>& getVKnots() const;

    const vector<double>& getWKnots() const;

    vector<double> getUMKnots() const;

    vector<double> getVMKnots() const;

    vector<double> getWMKnots() const;

    vector<double> getUAllKnots() const;

    vector<double> getVAllKnots() const;

    vector<double> getWAllKnots() const;

    int getUSize() const;
    int getVSize() const;
    int getWSize() const;

    int  getUDegree() const;

    int  getVDegree() const;

    int  getWDegree() const;

    bool getControlPointsByIndex(const std::vector<int>& localIndex, Eigen::MatrixXd& res) const;

    bool getControlPointsByCoef(const std::vector<int>& localIndex, const Eigen::MatrixXd& u_coef, Eigen::MatrixXd& res) const;

    vector<vector<vector<CPoint>>> getControlPoints() const;

    const vector<CPoint>& getCoefs() const;
    const vector<CPoint>& coefs() const;

    bool getJacobi(double u, double v, double w, std::vector<int>& localIndex, Eigen::MatrixXd& jac);

    bool getJacobi(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, Eigen::MatrixXd& jac) const;

    bool getJacobiU(double u, double v, double w, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement,
                    Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU) const;

    bool getJacobiU(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement,
                    Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU) const;



    bool updateCoefs(const Eigen::MatrixXd& coefs);
    bool setCoefs(const Eigen::MatrixXd& coefs);

    void eval_AllDers(std::vector<double> u, std::vector<double> v, std::vector<double> w, vector<Eigen::MatrixXd>& basisValue) const;

    void eval_into(std::vector<double> u, std::vector<double> v, std::vector<double> w, Eigen::MatrixXd& basisValue);

    void eval_2Ders(double u, double v, double w, vector<Eigen::VectorXd>& basisValue) const;

    void eval(double u, double v, double w, Eigen::VectorXd& basisValue);

    void active(int u_interval, int v_interval, int w_interval, std::vector<int>& localIndex) const;
    void activeFromParam(double u, double v, double w, std::vector<int>& localIndex) const;
    void getLocation(double u, double v, double w, Eigen::MatrixXd& P);

    CPoint getPointFromParam(double u, double v, double w);

    bool isBelongToTheFace(int dir, int id) { return facePointIndex[dir - 1].count(id) > 0; }

    inline Eigen::MatrixXd support() const
    {
        Eigen::MatrixXd supp;
        supp.resize(3, 2);
        supp << 0, 1, 0, 1, 0, 1;
        return supp;
    }

    Eigen::MatrixXd eval(const Eigen::MatrixXd& u);


};