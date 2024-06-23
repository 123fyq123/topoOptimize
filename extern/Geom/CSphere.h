#pragma once
#ifndef _CSPHERE_H_
#define _CSPHERE_H_

/// \brief 创建基本图形：球体曲面
#include<Geom_SphericalSurface.hxx>
#include "CGeomSurface.h"
#include"CPoint.h"
#include"DefineType.h"
#include <vector>
class CSphere :
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
	CSphere ();

	virtual~CSphere ();

	/// \brief 构造函数：创建球体
	/// \param point 球体位置坐标点
	/// \param radius 球体半径
	/// \note 半径必须是>0的参数
	CSphere (CPoint point,double radius);

public:
	/// \brief 设置位置
    /// \param point 新的位置
	void setLocation (CPoint point);

	/// \brief 设置半径
	/// \param radius 新的半径
	void setRadius (double radius);

	/// \brief 计算球表面积
	/// \return 球表面积
	double computeArea ();

	/// \brief 计算球体体积
	/// \return 球体积
	double computeVolume ();

	/// \brief 计算边界参数
	/// \return umin u方向最小参数边界
	/// \return vmin v方向最小参数边界
	void getUVmin (double&umin,double&vmin);


	/// \brief 计算边界参数
	/// \return umax u方向最大参数边界
	/// \return vmax v方向最大参数边界
	void getUVmax (double &umax,double &vmax);

	/// \brief 计算绝对笛卡尔坐标系下二次曲面隐式方程的系数
	/// 
    /// A1.X**2 + A2.Y**2 + A3.Z**2 + 2.(B1.X.Y + B2.X.Z + B3.Y.Z) + 
	/// 2.(C1.X + C2.Y + C3.Z) + D = 0.0
	void computeCoefficients (double &A1,double &A2,double &A3,double &B1,double &B2,double &B3,double &C1,double &C2,double &C3,double &D);

	/// \brief 判断U方向封闭性
	bool IsUClosed ();

	/// \brief 判断V方向封闭性
	bool IsVClosed ();

	/// \brief 判断U方向周期性
	bool IsUPeriodic ();

	/// \brief 判断V方向周期性
	bool IsVPeriodic ();

	/// \brief 计算参数处的坐标点
	/// \param u 参数值 
	/// \param v 参数值
	/// \param point 返回坐标点
	/// \note 参数值必须在边界范围内
	void getPointAtParam (double u,double v, CPoint &point);

	/// \brief 根据坐标点求相应参数
	///
	/// \return 若点不在曲面上，返回false.
	bool getParamAtPoint (const CPoint &point,double &u,double &v);

	/// \brief 计算曲面参数处的一阶导
	void getDerivation1 (double &u,double &v, vector<CVector> &vectors);

	/// \brief 计算曲面参数处的二阶导
	void getDerivation2 (double &u,double &v, vector<CVector> &vectors);

	/// \brief 计算曲面参数处的三阶导
	void getDerivation3 (double &u,double &v, vector<CVector> &vectors);

public:
	/// \brief 获取平移球体
	///
	/// \param vec 平移向量
	/// \param sphere 返回平移后的球体，当前球体无变化
	bool setTranslation (CVector &vec,CSphere &sphere);

	/// \brief 获取曲面类型
	int getSurfaceType ();

protected:
	//获取球面对象
	Handle_Geom_SphericalSurface getGeomSurface();

	//更新球体数据
	void upCSpheredata ();

private:
	CPoint mPoint;
	double mRadius;
	Handle_Geom_SphericalSurface mSurface;

};

#endif // !_CSPHERE_H_