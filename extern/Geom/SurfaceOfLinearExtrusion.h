#pragma once
#ifndef _SURFACEOFLINEAREXTRUSION_H_
#define _SURFACEOFLINEAREXTRUSION_H_

#include "CGeomSurface.h"
#include"CBSplineCurve.h"
#include<Geom_SurfaceOfLinearExtrusion.hxx>
#include <Geom_RectangularTrimmedSurface.hxx>
#include <vector>
/// \brief 线性拉伸曲面——沿直线放样
///
///如果挤出的曲线是圆,则曲面为圆柱面；
///如果挤出曲线是一条直线，则为平面。
///
class SurfaceOfLinearExtrusion :
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
	SurfaceOfLinearExtrusion ();

	/// \brief 析构函数
	virtual~SurfaceOfLinearExtrusion ();

	/// \brief 构造拉伸曲面。给定曲线和拉伸方向
	/// \param curve 基准曲线 
	/// \param axl 拉伸方向
	/// \note 如果拉伸方向与基准曲线在同一平面，拉伸出的曲面为曲边平面
	SurfaceOfLinearExtrusion (CBSplineCurve &curve, CVector &axl);

	/// \brief 设置方向
	/// \param axl 新的拉伸方向
	void setDirection ( CVector &axl);

	/// \brief 设置基准曲线
	/// \param curve 新的基准曲线
	/// \note 支持B样条和NURBS曲线
	void setBasisBSplineCurve ( CBSplineCurve &curve);

	/// \brief 获取曲面参数边界
	/// \param U1 u方向最小参数
	/// \param U2 u方向最大参数
	/// \param V1 v方向最小参数
	/// \param V2 v方向最大参数
	void getBounds (double &U1,double &U2,double &V1,double &V2);

	/// \brief u方向参数修正
	void uReverse ();

	/// \brief v方向参数修正
	void vReverse ();

	/// \brief 判断基准曲线是否封闭
	/// \return true 封闭  false 非封闭
	bool isClosed ();

	/// \brief 判断曲线是否满足N连续
	/// \param N 指定连续性
	/// \note N需要是>=0的整数
	bool isCNc (int N);

	/// \brief 计算拉伸曲面参数处坐标
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
	void getDerivation1 (double &u,double &v,vector<CVector> &vectors);

	/// \brief 计算曲面参数处的二阶导
	/// \param u U方向参数值
	/// \param v V方向参数值
	/// \param vectors 返回二阶导 0索引为u向二阶导，1索引为v向二阶导，2索引为uv二阶导，先偏导u，后偏导v
	/// \note u v参数必须在参数边界范围内
	void getDerivation2 (double &u,double &v,vector<CVector> &vectors);

	/// \brief 计算曲面参数处的三阶导
	/// \param u U方向参数值
	/// \param v V方向参数值
	/// \param vectors 返回三阶导 0索引为u向三阶导，1索引为v向三阶导，2索引为uuv三阶导，3索引为uvv三阶导
	/// \note u v参数必须在参数边界范围内
	void getDerivation3(double &u,double &v,vector<CVector> &vectors);
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
	//获取曲面对象
	Handle_Geom_RectangularTrimmedSurface getGeomSurface ();
private:
	bool getParamAtPointBasis (const CPoint &point,double &u,double &v);
private:
	CBSplineCurve mcurve;
	CVector maxl;
	Handle_Geom_SurfaceOfLinearExtrusion mSurface;
	Handle_Geom_RectangularTrimmedSurface mRSurface;
};

#endif // !_SURFACEOFLINEAREXTRUSION_H_