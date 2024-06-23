#pragma once

#include <unordered_set>
#include "Geom/CBSplineCurve.h"
#include "nrbNurbs.h"
#include "Eigen/Dense"
#include "Eigen/Core"
/// \brief b样条曲面类.
///
/// 提供一些基本b样条曲面操作.
/// 当前点的处理面向XYZ坐标系.
///
class CBSplineSolid
{
	friend class CAidFunction;
public:
	CBSplineSolid();
	virtual ~CBSplineSolid();

	/// \brief 构造.
	/// 
	/// \param splineSolid 拷贝构造.
	//CBSplineSolid(const CBSplineSolid& splineSolid);

	/// \brief 构造一个无理b样条体.
	/// 
	/// \param udegree u次数.
	/// \param vdegree v次数.
	/// \param wdegree w次数.
	/// \param controlPoints 控制点.
	/// \param uknots 节点序列（不重复的非递减序列）.
	/// \param vknots 节点序列（不重复的非递减序列）.
	/// \param wknots 节点序列（不重复的非递减序列）.
	/// \param umultis 记录每个uknot的重复度.
	/// \param vmultis 记录每个vknot的重复度.
	/// \param wmultis 记录每个wknot的重复度.
	/// \param uperiodic u向周期.
	/// \param vperiodic v向周期.
	/// \param vperiodic w向周期.
	CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints, const vector<double>& uknots, const vector<double>& vknots,
		const vector<double>& wknots,const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis, bool uperiodic = false, bool vperiodic = false,bool wperiodic = false);

	/// \brief 构造一个有理b样条体.
	/// 
	/// \param udegree u次数.
	/// \param vdegree v次数.
	/// \param vdegree w次数.
	/// \param controlPoints 控制点.
	/// \param uknots 节点序列（不重复的非递减序列）.
	/// \param vknots 节点序列（不重复的非递减序列）.
	/// \param wknots 节点序列（不重复的非递减序列）.
	/// \param umultis 记录每个uknot的重复度.
	/// \param vmultis 记录每个vknot的重复度.
	/// \param vmultis 记录每个wknot的重复度.
	/// \param weights 权重. 
	/// \param uperiodic u向周期.
	/// \param vperiodic v向周期.
	/// \param wperiodic w向周期.
	CBSplineSolid(const int& udegree, const int& vdegree, const int& wdegree, const vector<vector<vector<CPoint>>>& controlPoints,const vector<double>& uknots,
		const vector<double>& vknots,const vector<double>& wknots, const vector<int>& umultis, const vector<int>& vmultis, const vector<int>& wmultis,
		const vector<vector<vector<double>>>& weights, bool uperiodic = false, bool vperiodic = false, bool wperiodic = false);

public:
	/// \brief 获取曲面u控制点数目.
	///
	int getSizeUControlPoints() const;

	/// \brief 获取曲面v控制点数目.
	///
	int getSizeVControlPoints() const;

	/// \brief 获取曲面w控制点数目.
	///
	int getSizeWControlPoints() const;
	/// \brief 设置曲面u向次数.
	/// 
	/// \param udegree u向次数.
	/// \return true if the spline surface is not constructed, then the spline surface's udegree will be set.
	/// \note increase degree is not allowed.
	/// 
	
	int getSizeControlPoints() const;
	
	bool setUDegree(int udegree);

	/// \brief 设置曲面v向次数.
	/// 
	/// \param vdegree v向次数.
	/// \return true if the spline surface is not constructed, then the spline surface's vdegree will be set.
	/// \note increase degree is not allowed.
	bool setVDegree(int vdegree);
	/// \brief 设置曲面w向次数.
	/// 
	/// \param wdegree w向次数.
	/// \return true if the spline surface is not constructed, then the spline surface's vdegree will be set.
	/// \note increase degree is not allowed.
	bool setWDegree(int wdegree);

	/// \brief 设置/修改曲面控制点.
	/// 
	/// \param controlPoints 控制点.
	bool setControlPoints(const vector<vector<vector<CPoint>>>& controlPoints);

	/// \brief 设置/修改曲面控制点+权重.
	/// 
	/// \param controlPoints 控制点.
	/// \param weights 权重.
	bool setControlPoints(const vector<vector<vector<CPoint>>>& controlPoints, const vector<vector<vector<double>>>& weights);

	/// \brief 修改曲面单个控制点 controlPoints[u][v][w].
	/// \param uindex controlPoints 中的u索引,从1开始.
	/// \param vindex controlPoints 中的v索引，从1开始.
	/// \param vindex controlPoints 中的w索引，从1开始.
	/// \param controlPoint 新控制点值.
	bool setSingleControlPoint(int& uindex, int& vindex, int& windex, const CPoint& controlPoint);

	/// \brief 修改曲面单个控制点 controlPoints[u][v][w].
	/// \param uindex controlPoints 中的u索引,从1开始.
	/// \param vindex controlPoints 中的v索引，从1开始.
	/// \param vindex controlPoints 中的w索引，从1开始.
	/// \param controlPoint 新控制点值.
	/// \param weight 权重.
	bool setSingleControlPoint(int& uindex, int& vindex, int& windex, const CPoint& controlPoint, double& weight);

	/// \brief 设置/修改曲面权重.
	/// 
	/// \param weights 权重, 可以给无理b样条添加权重.
	/// \return true if the spline surface's weights are set.
	bool setWeights(const vector<vector<vector<double>>>& weights);

	/// \brief 修改曲线权重.
	/// 
	bool setSingleWeight(int uindex, int vindex, int windex, double weight);


	const vector<double>& getUKnots() const;

	const vector<double>& getVKnots() const;

	const vector<double>& getWKnots() const;

	vector<double> getUMKnots() const;

	vector<double> getVMKnots() const;

	vector<double> getWMKnots() const;

public:
	/// \brief 获取曲面控制点.
	///
	vector<vector<vector<CPoint>>> getControlPoints() const;

	const vector<CPoint>& coefs() const;

	/// \brief 获取曲面权重.
	///
	vector<vector<vector<double>>> getWeights() const;

	

	/// \brief 获取单个控制点 索引从1开始.
	///
	CPoint getSingleControlPoint(int uindex, int vindex,int windex);

	/// \brief 获取单个权重 索引从1开始.
	/// 
	double getSingleWeight(int uindex, int vindex,int windex);


	
	bool   getRational() const;

	/// \brief 获取udegree.
	///
	int    getUDegree() const;

	/// \brief 获取vdegree.
	///
	int    getVDegree() const;

	/// \brief 获取wdegree.
	///
	int    getWDegree() const;

	bool getControlPointsByIndex(const std::vector<int>& localIndex, Eigen::MatrixXd& res);

	bool getControlPointsByCoef(const std::vector<int>& localIndex, const Eigen::MatrixXd& u_coef, Eigen::MatrixXd& res);


	int gethandle()const;
	void sethandle(int handle);

public:
	vector<CPoint> U_0_discrete_point;
	int U0_size=0,U1_size=0,V0_size=0,V1_size=0,W0_size=0,W1_size=0;
	vector<CPoint> U_1_discrete_point;
	vector<CPoint> V_0_discrete_point;
	vector<CPoint> V_1_discrete_point;
	vector<CPoint> W_0_discrete_point;
	vector<CPoint> W_1_discrete_point;
	vector<vector<CPoint>> U_0_point_to_calcuNormal;
	vector<vector<CPoint>> U_1_point_to_calcuNormal;
	vector<vector<CPoint>> V_0_point_to_calcuNormal;
	vector<vector<CPoint>> V_1_point_to_calcuNormal;
	vector<vector<CPoint>> W_0_point_to_calcuNormal;
	vector<vector<CPoint>> W_1_point_to_calcuNormal;
	vector<vector<Eigen::Vector3f>> U_0_Normal;
	vector<vector<Eigen::Vector3f>> U_1_Normal;
	vector<vector<Eigen::Vector3f>> V_0_Normal;
	vector<vector<Eigen::Vector3f>> V_1_Normal;
	vector<vector<Eigen::Vector3f>> W_0_Normal;
	vector<vector<Eigen::Vector3f>> W_1_Normal;
	vector<vector<Eigen::Vector3f>> Face_Normal;
	void calcuNormal() {

	}
	void Clearmyvector() {
		U_0_discrete_point.clear();
	    U_1_discrete_point.clear();
	    V_0_discrete_point.clear();
		V_1_discrete_point.clear();
		W_0_discrete_point.clear();
		W_1_discrete_point.clear();
	}
	/// \brief 求曲线上uv值对应的坐标点.
	///
	CPoint getPointAtParam(double u, double v,double w);

	/// \brief 求曲线上uv值对应的坐标点.
	///
	bool getPointAtParam(double u, double v,double w, CPoint& point);
	/// \brief 根据坐标点求相应参数
	///
	/// \return 若点不在曲面上，返回false.
	bool getParamAtPoint(const CPoint& point, double& u, double& v,double &w);
	void eval(double u, double v, double w, Eigen::VectorXd& basisValue);
	/// \brief 求曲面umin umax vmin vmax wmin wmax
	///
	bool getUVMinMax(double& umin, double& umax, double& vmin, double& vmax, double& wmin, double& wmax);
	void getLocation(double u, double v, double w, Eigen::MatrixXd& P);
	void updateUKnots();
	void updateVKnots();
	void updateWKnots();


	vector<double> m_UKnots_nonredundant;   ///< u节点矢量	
	vector<double> m_VKnots_nonredundant;   ///< v节点矢量
	vector<double> m_WKnots_nonredundant;   ///< w节点矢量
	vector<double> m_UKnots;   ///< u节点矢量	
	vector<double> m_VKnots;   ///< v节点矢量
	vector<double> m_WKnots;   ///< w节点矢量
public:
	// useful function
	bool getJacobi(double u, double v, double w, std::vector<int>& localIndex, Eigen::MatrixXd& jac);

	bool getJacobi(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, Eigen::MatrixXd& jac);

	bool getJacobiU(double u, double v, double w, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement,
		Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU);

	bool getJacobiU(const Eigen::MatrixXd& basisValue, std::vector<int>& localIndex, const Eigen::MatrixXd& displacement,
		Eigen::MatrixXd& jac, Eigen::MatrixXd& jacU);

	bool updateCoefs(const Eigen::MatrixXd& coefs);

	void eval_AllDers(std::vector<double> u, std::vector<double> v, std::vector<double> w, vector<Eigen::MatrixXd>& basisValue);

	void eval_into(std::vector<double> u, std::vector<double> v, std::vector<double> w, Eigen::MatrixXd& basisValue);

	void eval_2Ders(double u, double v, double w, vector<Eigen::VectorXd>& basisValue);

	void active(int u_interval, int v_interval, int w_interval, std::vector<int>& localIndex) const;
	void activeFromParam(double u, double v, double w, std::vector<int>& localIndex) const;
	CPoint getPointFromParam(double u, double v, double w);

	bool isBelongToTheFace(int dir, int id) { return facePointIndex[dir - 1].count(id) > 0; }
private:
	bool createBSplineSolid();
	//void computeArea();
	


private:
	int SolidHandle;
	vector<vector<vector<CPoint>>> controlPoints;
	Eigen::VectorXd basisValue;
	Eigen::MatrixXd P;
	Eigen::MatrixXd Location;
	int m_iTypeCoord;  ///< 当前使用的坐标系 : 1 - 笛卡尔坐标(默认) / 0 - 极坐标/柱坐标系
	//double m_dArea;  ///< 曲面面积
	bool m_bRational;  ///< true: 有理 
	bool m_bUPeriodic;  ///< false:非周期 
	bool m_bVPeriodic;  ///< false:非周期 
	bool m_bWPeriodic;  ///< false:非周期 
	vector<vector<vector<CPoint>>> m_ControlPoints;	///< 记录控制顶点
	vector<CPoint> m_Coefs;
	
	vector<vector<vector<double>>> m_Weights;	///< 权重
	vector<int>    m_UMultis;   ///< 重复度
	vector<int>    m_VMultis;   //< 重复度
	vector<int>    m_WMultis;   //< 重复度

	int m_iUDegree; ///< u次数
	int m_iVDegree; ///< v次数
	int m_iWDegree; ///< w次数
	//Handle_Geom_BSplineSurface m_pSurface;///< 曲面
	double* u_knots;        // 完整节点
	double* v_knots;
	double* w_knots;
	int u_size = 0;
	int v_size = 0;
	int w_size = 0;
	std::vector<std::unordered_set<int>> facePointIndex;
};

