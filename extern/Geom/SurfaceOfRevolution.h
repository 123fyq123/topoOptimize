#pragma once
#ifndef _SURFACEOFREVOLUTION_H_
#define _SURFACEOFREVOLUTION_H_

//创建回转曲面
/*
如果基准线是一条线或一条被修剪的线，则为圆锥面；
如果基准线是一条与旋转轴平行的直线，则为圆柱面；
如果基准线是垂直于旋转轴的线，则为一个二次曲面。
*/
#include "CGeomSurface.h"
#include"CBSplineCurve.h"
#include<Geom_SurfaceOfRevolution.hxx>
#include <vector>
class SurfaceOfRevolution :
	public CGeomSurface{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	/// \brief 默认构造函数
	SurfaceOfRevolution ();
	/// \brief 析构函数
	virtual~SurfaceOfRevolution ();
public:
	/// \brief 创建回转曲面
	/// \param curve 是基准曲线
	/// \param point 是基准点
	/// \param axl 旋转方向
	SurfaceOfRevolution (CBSplineCurve &curve,CPoint &point,CVector &axl);

	/// \brief 设置旋转轴
	/// \param point 是基准点
	/// \param axl 旋转方向
	void setDirection (CPoint&point,CVector &axl);

	/// \brief 设置基准曲线
	void setBasisBSplineCurve (CBSplineCurve &curve);

	/// \brief 获取旋转轴的位置点
	CPoint getLocation ();

	/// \brief 获取旋转轴方向
	CVector getAxis ();

	/// \brief 获取曲面边界
	/// \param U1 u方向最小参数
	/// \param U2 u方向最大参数
	/// \param V1 v方向最小参数
	/// \param V2 v方向最大参数
	void getBounds (double &U1,double &U2,double &V1,double &V2);

	/// \brief 判断基准曲线是否封闭
	/// \return true 封闭 false 非封闭
	bool isClosed ();

	/// \brief 判断曲线是否满足N连续
	/// \param N 连续性
	/// \note N应该满足N>=0
	bool isCNc (int N);

	/// \brief 计算曲面参数处坐标
	/// \param u U方向参数值
	/// \param v V方向参数值
	/// \param point 返回参数处坐标点
	/// \note u v参数必须在参数边界范围内
	void getPointAtParam (double &u,double &v, CPoint &point);

	/// \brief 根据坐标点求相应参数
	///
	/// \return 若点不在曲面上，返回false.
	bool getParamAtPoint (const CPoint &point,double &u,double &v);

	/// \brief 计算曲面参数处的一阶导
	/// \param u U方向参数值
	/// \param v V方向参数值
	/// \param vectors 返回一阶导 0索引为u向一阶导，1索引为v向一阶导
	/// \note u v参数必须在参数边界范围内
	void getDerivation1 (CPoint &point,double &u,double &v,vector<CVector> &vectors);

	/// \brief 计算曲面参数处的二阶导
	/// \param u U方向参数值
	/// \param v V方向参数值
	/// \param vectors 返回二阶导 0索引为u向二阶导，1索引为v向二阶导，2索引为uv二阶导，先偏导u，后偏导v
	/// \note u v参数必须在参数边界范围内
	void getDerivation2 (CPoint &point,double &u,double &v,vector<CVector> &vectors);

	/// \brief 计算曲面参数处的三阶导
    /// \param u U方向参数值
    /// \param v V方向参数值
	/// \param vectors 返回三阶导 0索引为u向三阶导，1索引为v向三阶导，2索引为uuv三阶导，3索引为uvv三阶导
    /// \note u v参数必须在参数边界范围内
	void getDerivation3 (CPoint &point,double &u,double &v,vector<CVector> &vectors);

public:
	/// \brief 平移曲面.
	/// \param vec  平移向量.
	/// \note 返回旋转后曲面的新的控制点和权值.
	void setTranslation (CVector vec);

	/// \brief 旋转曲面.
	/// \param p 定位点.
	/// \param ax  旋转轴方向向量.
	/// \param ang  旋转角度（逆时针）.
	/// \note 返回旋转后曲面的新的控制点和权值.
	void setRotation (CPoint p,CVector ax,double ang);
public:
	/// \brief 获取曲面类型.
	int getSurfaceType ();
protected:
	Handle_Geom_SurfaceOfRevolution getGeomSurface ();
private:
	Handle_Geom_SurfaceOfRevolution mSurface;
	CPoint mPoint;
	CBSplineCurve mCurve;

};
#endif // !_SURFACEOFREVOLUTION_H_
