#pragma once
#ifndef _CCYLINDRICALSURFACE_H_
#define _CCYLINDRICALSURFACE_H_

#include <Geom_CylindricalSurface.hxx>
#include <vector>
#include <Geom_RectangularTrimmedSurface.hxx>
#include "CGeomSurface.h"
#include "CPoint.h"
/// \brief 圆柱曲面类.
///
/// 提供一些基本圆柱曲面操作.
/// 当前点的处理面向XYZ坐标系.
///
class CCylindricalSurface :public CGeomSurface
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	CCylindricalSurface();

	/// \brief 构造.
	/// 
	/// \param surface 拷贝构造.
	CCylindricalSurface(const CCylindricalSurface &surface);

	/// \brief 构造圆柱曲面.
	/// 
	/// \param point 起点.
	/// \param vec 轴向,设定轴长为圆柱高度.
	/// \param radius 半径.
	CCylindricalSurface(const CPoint &point, const CVector &vec, const double &radius);
	
	virtual ~CCylindricalSurface();

public:
	/// \brief 获取圆柱曲面类型.
	/// 
	/// \return  CYLINDRICALSURFACE = 4.
	int getSurfaceType();

public:
	/// \brief 设置圆柱面半径
	///
	/// \param radius 半径.
	bool setRadius(double radius);

	/// \brief 设置圆柱面起点.
	///
	bool setLocation(CPoint &location);

	/// \brief 设置圆柱面轴向,设定轴长为圆柱高度.
	///
	bool setDirection(CVector &direction);

	/// \brief 构造曲面
	///
	/// \note 半径 起点 轴向要求已经合理设置.
	bool buildSurface();
public:
	/// \brief 获取曲面半径.
	///
	bool getRadius(double &radius);

	/// \brief 获取曲面半径.
	///
	double getRadius();

	/// \brief 获取曲面起点.
	///
	bool getLocation(CPoint &point);

	/// \brief 获取曲面轴向.
	///
	bool getDirection(CVector &vec);

	/// \brief 获取曲面上对应uv值点的切向量 | 一阶导.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回一阶导 0索引为u向一阶导，1索引为v向一阶导.
	bool getDerivation1(double u, double v, vector<CVector> &vectors);

	/// \brief 获取曲面上对应uv值点的二阶导，曲面至少C2连续.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回二阶导 0索引为u向二阶导，1索引为v向二阶导，2索引为uv二阶导，先偏导u，后偏导v.
	bool getDerivation2(double u, double v, vector<CVector> &vectors);

	/// \brief 获取曲面上对应uv值点的三阶导，曲面至少C3连续.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回三阶导 0索引为u向三阶导，1索引为v向三阶导，2索引为uuv三阶导，3索引为uvv三阶导.
	bool getDerivation3(double u, double v, vector<CVector> &vectors);

	/// \brief 求曲线上uv值对应的坐标点.
    ///
	bool getPointAtParam(double u, double v, CPoint &point);
	CPoint getPointAtParam(double u,double v);
	/// \brief 根据坐标点求相应参数
	///
	/// \return 若点不在曲面上，返回false.
	bool getParamAtPoint(const CPoint &point, double &u, double &v); 

	/// \brief 求曲面umin umax vmin vmax
	///
	bool getUVMinMax(double &umin, double &umax, double &vmin, double &vmax);


public:
	/// \brief 获取平移曲面
	///
	/// \param vec 平移向量
	/// \param surface 返回平移后的曲面，当前曲面无变化
	/// \note 目前仅对圆柱起点操作
	bool setTranslation(CGeomPoint & vec, CCylindricalSurface &surface);

	/// \brief 获取旋转曲面
	///
	/// \param axisOrigin 旋转轴的起点
	/// \param axisDirection  旋转轴的方向
	/// \param angle 旋转弧度 逆时针旋转
	/// \param surface 返回旋转后的曲面，当前曲面无变化
	/// \note 目前仅对圆柱起点操作
	bool setRotation(CGeomPoint &axisOrigin, CGeomPoint &axisDirection, double angle, CCylindricalSurface &surface);

protected:
	Handle_Geom_RectangularTrimmedSurface getGeomSurface() const;
	

private:
	bool createCylindricalSurface();
	bool getParamAtPointBasis(const CPoint &point, double &u, double &v);

private:
	CPoint m_Location;   ///< 圆柱起点，typecoord = -1 为无效.
	CVector m_Direction; ///< 圆柱轴向，typecoord = -1 为无效，轴长为0无效.
	double m_dRadius;    ///< 圆柱半径.
	Handle_Geom_RectangularTrimmedSurface  m_pSurface;  ///< 有界曲面.
	Handle_Geom_CylindricalSurface m_pSurfaceBasis; ///< 无界曲面.
};

#endif // !_CCYLINDRICALSURFACE_H_