#pragma once
#ifndef _CPLANE_H_
#define _CPLANE_H_
#include "CGeomSurface.h"
#include "CPoint.h"
#include<Geom_Plane.hxx>
#include<vector>
#include <Geom_RectangularTrimmedSurface.hxx>
/// \brief 平面类.
///
/// 仅支持XYZ坐标系.
/// 
class CPlane : public CGeomSurface
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	CPlane();

	/// \brief 点法式构造平面.
	///
	/// \param point 位于平面上的点.
	/// \param vec 平面法向.
	CPlane( CPoint &point,  CVector &vec);

	/// \brief 三点构成平面.
	/// 
	CPlane( CPoint &point1,  CPoint &point2,  CPoint &point3);

	virtual ~CPlane();

public:

	/// \brief 获取法向.
	/// 
	/// \param normal    返回平面法向.
	void getNormal(CGeomPoint &normal);

	/// \brief 获取点.
	/// 
	/// \param point    返回平面的点.
	void getPoint(CGeomPoint &point);

	//获取边界
	void getBounds (double &u1,double &u2,double &v1,double&v2);
	//获取参数处点坐标
	void getPointAtParam (CPoint &point,double u,double v);
	//获取边界点
	void getBoundPoint (vector<CPoint> &boundpoint);

	/// \brief 获取平面类型.
	/// 
	/// \return  PLANE = 1.
	int getSurfaceType();

	bool getParamAtPointBasis(const CPoint &point, double &u, double &v);

protected:
	//获取平面对象
	Handle_Geom_Plane getGeomPlane ();
	Handle_Geom_RectangularTrimmedSurface getRGeomPlane ();

private:
	CPoint m_cPoint;      ///< 记录点法式中的点.
	CVector m_cNormal;	  ///< 记录点法式中的法向.
	CPoint m_cPoint1;     ///< 记录三点式中的点.
	CPoint m_cPoint2;     ///< 记录三点式中的点.
	CPoint m_cPoint3;     ///< 记录三点式中的点.
	Handle_Geom_Plane Plane;//occ平面对象
	Handle_Geom_RectangularTrimmedSurface mRPlane;//有界平面
};

#endif // !_CPLANE_H_