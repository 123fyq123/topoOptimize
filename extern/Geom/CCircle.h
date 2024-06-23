#pragma once
#ifndef _CCIRCLE_H_
#define _CCIRCLE_H_
#include <Geom_TrimmedCurve.hxx>
#include <Geom_Circle.hxx>
#include "CGeomCurve.h"
#include "CPoint.h"

class CCircle : public CGeomCurve
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;

public:
	CCircle();
	virtual ~CCircle();

	/// \brief 三点构造圆.
	/// 
	/// \param p1 点 | 圆心.
	/// \param p2 点 | 点.
	/// \param p3 点 | 点.
	/// \param type 0 = 圆上三点构造圆，1=圆心+两点构造圆
	CCircle(CPoint &p1, CPoint &p2, CPoint &p3, int type);

	/// \brief 三点构造圆弧.
	/// 
	/// \param p1 点 | 圆心.
	/// \param p2 点 | 点.
	/// \param p3 点 | 点.
	/// \param type 0 = 三点构造圆弧，1=圆心+两端点构造圆弧
	void buildArc(CPoint &p1, CPoint &p2, CPoint &p3, int type);

public:

	/// \brief 获取曲线类型.
	/// 
	/// \return  CIRCLE = 4 ARC = 5.
	int getCurveType();

public:

	/// \brief 获取圆心.
	/// 
	/// \param center 接收圆心.
	void getCenter(CPoint &center);
	
	/// \brief 获取圆心.
	/// 
	/// \return 返回圆心.
	CPoint getCenter();

	/// \brief 获取半径.
	/// 
	/// \param radius 接收半径.
	void getRadius(double &radius);
	
	/// \brief 获取半径.
	/// 
	/// \return 返回半径.
	double getRadius();

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
	/// \note 仅对圆弧有效.
	CPoint getStartPoint();

	/// \brief 获取曲线起点.
	///
	/// \note 仅对圆弧有效.
	bool getStartPoint(CPoint &point);

	/// \brief 获取曲线终点.
	///
	/// \note 仅对圆弧有效.
	CPoint getEndPoint();

	/// \brief 获取曲线终点.
	///
	/// \note 仅对圆弧有效.
	bool getEndPoint(CPoint &point);

protected:
	Handle_Geom_Curve getGeomCurve();
	//void setGeomCurve(Handle_Geom_Curve &c);

private:
	/// \brief 三点计算圆心，半径，轴.
	/// 
	/// \param p1 点.
	/// \param p2 点.
	/// \param p3 点.
	/// \note 要求三点共面
	bool compute(CPoint &p1, CPoint &p2, CPoint &p3);
	/// \brief 检查三点是否共线
	bool check(CPoint &p1, CPoint &p2, CPoint &p3);

private:
	int m_iTypeCoord;  ///< 当前使用的坐标系 : 1 - 笛卡尔坐标(默认) / 0 - 极坐标/柱坐标系.
	double m_dRadius;  ///<所属圆的半径.
	CPoint m_cCenter;  ///<圆心.
	CVector m_cAxis;   ///<轴.
	Handle_Geom_Circle m_pCircle;///所属圆.
	Handle_Geom_TrimmedCurve m_pArc;///圆弧.
};
#endif//_CCIRCLE_H_
