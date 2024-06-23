#pragma once
#ifndef _CBSPLINECURVE_H_
#define _CBSPLINECURVE_H_

#include <Geom_BSplineCurve.hxx>
#include <vector>
#include "CGeomCurve.h"
#include "CPoint.h"

/// \brief b样条曲线类.
///
/// 提供一些基本b样条曲线操作.
/// 当前点的处理面向XYZ坐标系.
///
class CBSplineCurve : public CGeomCurve
{
	friend class CAidFunction;
	friend class CBSplineSurface;
	friend class SurfaceOfLinearExtrusion;
	friend class SurfaceOfRevolution;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	CBSplineCurve();
	virtual ~CBSplineCurve();
	
	/// \brief 构造.
	/// 
	/// \param splineCurve 拷贝构造.
	CBSplineCurve(const CBSplineCurve &splineCurve);

	/// \brief 构造一条无理b样条曲线.
	/// 
	/// \param degree 次数.
	/// \param controlPoints 控制点.
	/// \param knots 节点序列（不重复的非递减序列）.
	/// \param multis 记录每个knot的重复度.
	/// \param periodic 周期.
	/// \note size(knots) = size(multis); sum(multis) = 1 + degree + size(controlPoints)
	CBSplineCurve(const int degree, const vector<CPoint> &controlPoints, const vector<double> &knots, const vector<int> &multis, bool periodic = false);

	/// \brief 构造一条有理b样条曲线.
	/// 
	/// \param degree 次数.
	/// \param controlPoints 控制点.
	/// \param knots 节点序列（不重复的非递减序列）.
	/// \param multis 记录每个knot的重复度.
	/// \param weights 每个控制点的权重，必须大于0.
	/// \param periodic 周期.
	/// \note size(knots) = size(multis); size(controlPoints) = size(weights); sum(multis) = 1 + degree + number(controlPoints)
	CBSplineCurve(const int degree, const vector<CPoint> &controlPoints, const vector<double> &knots, const vector<int> &multis, const vector<double> &weights, bool periodic = false/*, int typeCoord = 1*/);
	
	/// \brief 插值构造一条无理b样条曲线.
	/// 
	/// \param curvePoints 插值点.
	/// \param tol 插值精度，默认为1e-10.
	CBSplineCurve(const vector<CPoint> &curvePoints, double tol = 1.e-10);

	/// \brief 插值构造一条有理b样条曲线.
	/// 
	/// \param curvePoints 插值点.
	/// \param weights 插值点权重，必须大于0
	/// \param tol 插值精度，默认为1e-10.
	CBSplineCurve(const vector<CPoint> &curvePoints, const vector<double> &weights, double tol = 1.e-10);

	/// \brief 曲线拼接，要求每段曲线前后点连续.
	/// 
	/// \param curves 待合并曲线.
	/// \param tol 插值精度，默认为1e-10.
	CBSplineCurve(const vector<CBSplineCurve> &curves, double tol = 1.e-10);

public:

	/// \brief 获取b样条曲线类型.
	/// 
	/// \return  BSPLINECURVE = 2.
	int getCurveType();

public:
	
	/// \brief 设置曲线次数.
	/// 
	/// \param degree 次数.
	/// \return true if the spline curve is not constructed, then the spline curve's degree will be set.
	/// \note increase degree is not allowed.
	bool setDegree(int degree);  

	/// \brief 设置/修改曲线控制点.
	/// 
	/// \param controlPoints 控制点.
	/// \return true if the spline curve's control points are set.
	/// \note if the curve is not null, the number of new control points should be the same as the old one .
	bool setControlPoints(const vector<CPoint> &controlPoints);

	/// \brief 设置/修改曲线控制点+权重.
	/// 
	/// \param controlPoints 控制点.
	/// \param weights 权重.
	/// \return true if the spline curve's control points and weights are set.
	/// \note if the curve is not null, the number of new control points should be the same as the old one, weights is the same .
	bool setControlPoints(const vector<CPoint> &controlPoints, const vector<double> &weights);

	/// \brief 修改曲线单个控制点.
	/// 
	/// \param index controlPoints 中的索引，from 1 - numOfControlPoints.
	/// \param controlPoint 控制点.
	bool setSingleControlPoint(int index, const CPoint &controlPoint);

	/// \brief 修改曲线单个控制点+权重.
	/// 
	/// \param index controlPoints 中的索引，from 1 - numOfControlPoints.
	/// \param controlPoint 控制点
	/// \param weight 权重.
	bool setSingleControlPoint(int index, const CPoint &controlPoint, double weight); 

	/// \brief 设置/修改曲线权重.
	/// 
	/// \param weights 权重, 可以给无理b样条添加权重.
	/// \return true if the spline curve's weights are set.
	bool setWeights(const vector<double> &weights);

	/// \brief 修改曲线权重.
	/// 
	/// \param index from 1 - numOfControlPoints.
	/// \param weight 权重.
	bool setSingleWeight(int index, double weight); 

	/// \brief 设置/修改曲线节点序列.
	/// 
	/// \param knots 非递减不重复序列.
	/// \return true if the spline curve's knots are set.
	bool setKnots(const vector<double> &knots);

	/// \brief 设置/修改曲线节点序列.
	/// 
	/// \param index from 1 - numOfKnots.
	/// \param knot 非递减不重复序列.
	bool setSingleKnot(int index, double knot);

	/// \brief 设置/修改曲线节点序列重复度.
	/// 
	/// \param multis 序列重复度.
	/// \return true if the spline curve's multis are set.
	bool setMultis(const vector<int> &multis);

	/// \brief 设置/修改曲线节点序列.
	/// 
	/// \param index from 1 - numOfKnots.
	/// \param multi 必须大于原值，小于等于degree.
	bool setSingleMultis(int index, int multi);

	/// \brief 设置曲线周期.
	///
	bool setPeriodic(bool periodic);

	/// \brief 插入节点
	/// 
	/// \param u u∈[umin, umax].
	/// \param multi 插入节点重复度
	bool insertKnot(double u, int multi); 

	/// \brief 设置/修改曲线节点序列.
	/// 
	/// \param index from 1 - numOfKnots.
	/// \param multi 删除节点重复度，0表示全删.
	bool removeKnot(int index, int multi); 

	/// \brief 构造曲线，当b样条为空，后期配置好各个参数，可调用这个方法构造.
	/// \return true when the degree, control points, knots, multis (weights) are set and obey the role : sum (multis) = degree + number(control points) + 1. 
	bool buildCurve();

public:
	/// \brief 获取曲线控制点.
	///
	vector<CPoint> getControlPoints() const;

	/// \brief 获取曲线控制点数目.
	///
	int getSizeControlPoints() const;

	/// \brief 获取曲线knots.
	///
	vector<double> getKnots() const;

	/// \brief 获取曲线knots数目.
	///
	int getSizeKnots() const;

	/// \brief 获取曲线权重.
	///
	vector<double> getWeights() const;

	/// \brief 获取曲线权重数目.
	///
	int getSizeWeight() const;

	/// \brief 获取曲线multis.
	///
	vector<int> getMultis() const;

	/// \brief 获取曲线multis数目.
	///
	int getSizeMultis() const;

	/// \brief 获取曲线单个控制点
	///
	/// \param index from 1 - numOfControlPoints.
	CPoint getSingleControlPoint(int index); 

	/// \brief 获取曲线单个权重
	///
	/// \param index from 1 - numOfControlPoints.
	double getSingleWeight(int index);

	/// \brief 获取曲线knot
	///
	/// \param index // from 1 - numOfKnots.
	double getSingleKnot(int index);

	/// \brief 获取曲线knot的multi
	///
	/// \param index // from 1 - numOfKnots.
	int getSingleMulti(int index);

	/// \brief 获取曲线长度.
	/// 
	double getLength() const;  

	/// \brief 获取曲线有理状态.
	///
	bool   getRational() const;

	/// \brief 获取曲线次数
	///
	int    getDegree() const;
	/// \brief 获取曲线闭合状态.
	///
	/// \return true closed   false  not-closed
	bool   getClosed() const;

	/// \brief 获取曲线周期状态.
	///
	bool   getPeriodic() const;

public:

	/// \brief 获取曲线上对应u值点的切向量 | 一阶导.
	/// 
	/// \param u .
	/// \param vector 返回一阶导.
	bool getDerivation1(double u, CVector &vector);

	/// \brief 获取曲线上对应u值点的二阶导，曲线至少C2连续.
	/// 
	/// \param u .
	/// \param vector 返回二阶导.
	bool getDerivation2(double u, CVector &vector);

	/// \brief 获取曲线上对应u值点的三阶导，曲线至少C3连续.
	/// 
	/// \param u .
	/// \param vector 返回三阶导.
	bool getDerivation3(double u, CVector &vector);

	/// \brief 求曲线上u值对应的坐标点.
	///
	CPoint getPointAtParam(double u);

	/// \brief 求曲线上u值对应的坐标点.
	///
	bool getPointAtParam(double u, CPoint &point);
	
	/// \brief 根据坐标点求相应参数
	///
	/// \return 若点不在曲线上，返回false.
	bool getParamAtPoint(const CPoint &point, double &u); 

	/// \brief 求曲线umin.
	///
	double getUmin();

	/// \brief 求曲线umax.
	///
	double getUmax();

	/// \brief umax - umin.
	///
	double getDeltaU();

	/// \brief 获取曲线起点.
	///
	CPoint getStartPoint();

	/// \brief 获取曲线起点.
	///
	bool getStartPoint(CPoint &point);

	/// \brief 获取曲线终点.
	///
	CPoint getEndPoint();

	/// \brief 获取曲线终点.
    ///
	bool getEndPoint(CPoint &point);



public:

	/// \brief 获取平移曲线
	///
	/// \param vec 平移向量
	/// \param spline 返回平移后的曲线，当前曲线无变化.
	bool setTranslation(CGeomPoint & vec, CBSplineCurve & spline);

	/// \brief 获取旋转曲线
	///
	/// \param axisOrigin 旋转轴的起点
	/// \param axisDirection  旋转轴的方向
	/// \param angle 旋转弧度 逆时针旋转
	/// \param spline 返回旋转后的曲线，当前曲线无变化.
	bool setRotation(CGeomPoint &axisOrigin, CGeomPoint &axisDirection, double angle, CBSplineCurve &spline);

	/// \brief 裁剪曲线
	///
	/// \param u1
	/// \param u2   u1, u2 需要在[umin, umax], u2 > u1.
	bool setTrimmed(double u1, double u2);

	/// \brief 曲线对称.
	/// \param plane 平面.
	/// \param spline 返回对称后的曲线.
	bool setMirror(CGeom &plane, CBSplineCurve &spline);

	/// \brief 曲线反向修正
	void reverse ();

	// \brief 曲线延伸
	// \param length 曲线尾点延伸长度控制.
	// \param spline 返回延伸后的曲线，当前曲线无变化.
	// \return 延伸是否成功.
	// \note 曲线尾端切方向延伸,曲线要求为准均匀样条曲线或者插值实现的曲线.
	bool setExtendLength(double length, CBSplineCurve &spline);

protected:
	Handle_Geom_BSplineCurve getGeomCurve() const;
	void setGeomCurve(Handle_Geom_BSplineCurve &c) ;

private:
	bool createBSplineCurve();
	void computeLength();
	void updateCurveDate();
	void updateControlPoints();
	void updateWeights();
	void updateKnots();
	void updateMultis();

private:
	int m_iTypeCoord;  ///< 当前使用的坐标系 : 1 - 笛卡尔坐标(默认) / 0 - 极坐标/柱坐标系.
	double m_dLength;  ///< 曲线长.
	bool m_bRational;  ///< true: 有理 .
	bool m_bPeriodic;  ///< false:非周期 .
	vector<CPoint> m_ControlPoints;	///< 记录控制顶点.
	vector<double> m_Knots;   ///< 节点矢量	.
	vector<double> m_Weights;	///< 权重.
	vector<int>    m_Multis;   ///< 重复度.

	int m_iDegree; ///< 次数.
	Handle_Geom_BSplineCurve m_pCurve;///< 曲线.
};

#endif // !_CBSPLINECURVE_H_