#pragma once
#ifndef _CLINE_H_
#define _CLINE_H_
#include <Geom_Line.hxx>
#include <Geom_TrimmedCurve.hxx>
#include "CPoint.h"
#include "CGeomCurve.h"
#include "CGeomSurface.h"

/// \brief 直线类.
///
/// 仅支持XYZ坐标系.
/// 
class CLine : public CGeomCurve
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	CLine();

	/// \brief 点向式构造直线.
	///
	/// \param point 位于线上的点.
	/// \param vec 线矢量.
	CLine(CGeomPoint &point, CGeomPoint &vec);

	/// \brief 平面相交构造直线.
	///
	/// \param plane1 平面1.
	/// \param plane2 平面2.
	CLine(CGeomSurface &plane1, CGeomSurface &plane2);  //3D中可以两个平面构造一条直线
	virtual ~CLine();

public:

	/// \brief 获取直线类型.
	/// 
	/// \return  LINE = 1.
	int getCurveType() const;

public:
	/// \brief 设置线上的点.
	///
	/// \param point .
	void setStartPoint(CPoint &point);

	/// \brief 设置线矢量.
	///
	/// \param vec .
	void setDirection(CVector &vec);

	/// \brief 获取线上的点.
	///
	/// \param u .
	/// \return 返回点.
	CPoint getPointAtParam(double u);

	/// \brief 获取线的起点， u = 0.
	///
	/// \return 返回点.
	CPoint getStartPoint(void);

	/// \brief 获取线的终点， 起点+线矢量.
	///
	/// \return 返回点.
	CPoint getEndPoint(void);

	/// \brief 获取线起点参数， u = 0.
	///
	double getUmin(void);

	/// \brief 获取线终点参数， 线矢量长度.
	///
	double getUmax(void);

	/// \brief 获取线矢量.
	///
	/// \param vec 矢量.
	void getDirection(CGeomPoint &vec);

	/// \brief 获取类型.
	///
	/// \return LINE = 1.
	int getCurveType(); 

public:
	/// \brief 直线旋转.
	/// \param line 线.
	/// \param angle .
	/// \param typeAngle angle类型.
	/// \param result 返回旋转后的直线.
	void rotate(CGeom &line, double angle, int typeAngle, CLine &result);

	/// \brief 直线对称.
	/// \param plane 平面.
	/// \param result 返回对称后的直线.
	void mirror(CGeom &plane, CLine &result);

	/// \brief 直线平移.
	/// \param vec 平移方向及大小.
	/// \param result 返回平移后的直线.
	void translate(CGeomPoint &vec, CLine &result);

	/// \brief 线段长度计算.
	/// \param u1 线段参数1.
	/// \param u2 线段参数2.
	/// \return 线段长度.
	double computeLength(double u1, double u2);

	/// \brief 获取线段.
	/// \param u1 线段起始点.
	/// \param u2 线段终点.
	/// \param line 返回裁剪直线.
	void setTrimmed(double u1, double u2, CLine &line);

protected:
	Handle_Geom_Curve getGeomCurve();
	void setGeomCurve(Handle_Geom_Line &c);

private:
	double m_dStart;///< 线段起始点，默认0
	double m_dEnd;///< 线段终点，默认vector的长度
	CPoint m_cPoint;    ///< 起始点
	CVector m_cVector;  ///< 线矢量,包含线段长度信息

	Handle_Geom_Line m_pCurve; ///< 直线
	//Handle_Geom_TrimmedCurve m_pTrimmedCurve; ///< 线段
};

#endif // !_CLINE_H_