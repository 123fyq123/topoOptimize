#pragma once
#ifndef _CPOINT_H_
#define _CPOINT_H_
#include <gp_Pnt.hxx>
#include "CVector.h"
/// \brief 点类.
///
/// 点类操作.
/// 当前点的处理面向XYZ坐标系，柱坐标（其他类不支持）.
///
class CPoint : public CGeomPoint
{
public:
	/// \brief 构造.
	/// 
	/// 默认构造一个点，
	/// （0, 0, 0）.
	///
	CPoint();

	/// \brief 构造.
	/// 
	/// \param p1         x坐标   | r.
	/// \param p2         y坐标   | theta.
	/// \param p3         z坐标.
	/// \param typeCoord 1 XYZ坐标| 0 柱坐标.
	/// \param typeAngle 1 degree | 0 radian.
	CPoint(double p1, double p2, double p3, int typeCoord = 1, int typeAngle = 0);

	/// \brief 构造.
	/// 
	/// \param point 拷贝构造.
	CPoint(const CPoint &point);

	/// \brief 构造.
	/// 
	/// \param pnt 待转换的occ点类.
	CPoint(const gp_Pnt &pnt);

	/// \brief 构造.
	/// 
	/// \param vec 待转换的occ向量类.
	CPoint(const gp_Vec &vec);

	virtual ~CPoint();

	/// \brief 设置坐标类型
	/// \param type 1 XYZ坐标 | 0 柱坐标 | -1 无效.
	bool setTypeCoord(int type);

	/// \brief 获取坐标类型.
	/// 
	/// \return 返回坐标系类型.
	int getTypeCoord() const;
public: 

	/// \brief 设置/修改x坐标.
	/// 
	/// \param x x坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setX(double x);

	/// \brief 获取x坐标.
	/// 
	/// \return 返回x坐标值.
	double getX() const;

	/// \brief 设置/修改y坐标.
	/// 
	/// \param y y坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setY(double y);

	/// \brief 获取y坐标.
	/// 
	/// \return  返回y坐标值.
	double getY() const;

	/// \brief 设置/修改z坐标.
	/// 
	/// \param z z坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setZ(double z);

	/// \brief 获取z坐标.
	/// 
	/// \return  返回z坐标值.
	double getZ() const;

	/// \brief 设置/修改xyz坐标.
	/// 
	/// \param x x坐标.
	/// \param y y坐标.
	/// \param z z坐标.
	/// \param iTypeCoord 坐标系类型，目前仅支持XYZ坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setXYZ(double x, double y, double z, int iTypeCoord = 1);

	/// \brief 获取xyz坐标.
	/// 
	/// \param x 返回x坐标.
	/// \param y 返回y坐标.
	/// \param z 返回z坐标.
	/// \return   true 获取成功，false 获取失败.
	bool getXYZ(double &x, double &y, double &z);

public:
	/// \brief 设置/修改r.
	/// 
	/// \param r 径长.
	/// \return  true 设置成功，false 设置失败.
	bool setR(double r);

	/// \brief 获取r.
	/// 
	/// \return  返回r.
	double getR(void);

	/// \brief 获取angle类型.
	/// 
	/// \return 返回angle类型.
	int getTypeAngle(void);

	/// \brief 设置/修改r.
	/// 
	/// \param theta degree | redian.
	/// \param typeAngle 1  |   0.
	/// \return  true 设置成功，false 设置失败.
	bool setTheta(double theta, int typeAngle);

	/// \brief 获取theta.
	/// 
	/// \return 返回theta.
	double getTheta(void);

	/// \brief 设置/修改xyz坐标.
	/// 
	/// \param r .
	/// \param theta .
	/// \param z z坐标.
	/// \param typeAngle angle类型.
	/// \return  true 设置成功，false 设置失败.
	bool setPolar(double r, double theta, double z, int typeAngle);

	/// \brief 获取xyz坐标.
	/// 
	/// \param r         返回r.
	/// \param theta     返回theta.
	/// \param z         返回z坐标.
	/// \param typeAngle 返回angle类型
	/// \return   true 获取成功，false 获取失败.
	bool getPolar(double &r, double & theta, double &z, int &typeAngle);

public:
	/// \brief +重载.
	/// 
	CPoint operator+(const CPoint &point);

	/// \brief -重载.
	/// 
	CPoint operator-(const CPoint &point);
	
	/// \brief *重载.
	/// 
	CPoint operator*(const CPoint &point);

	/// \brief /重载.
	/// 
	CPoint operator/(const CPoint &point);

	/// \brief +重载.
	///
	CPoint operator+(const CVector &vec);

	/// \brief *重载.
	/// 
	CPoint operator*(const double d);
	bool operator ==(const CPoint& point);

public:
	/// \brief 返回距原点的距离.
	/// 
	double norm(void) const; 

	/// \brief 点取反.
	/// 
	bool neg(void);

	/// \brief 单位化.
	/// 
	void normaliz(); 

	/// \brief 点旋转.
	/// \param line 线.
	/// \param angle .
	/// \param typeAngle angle类型.
	/// \param point 返回旋转后的点.
	void rotate(CGeom &line, double angle, int typeAngle, CGeomPoint &point); 

	/// \brief 点对称.
	/// \param plane 平面.
	/// \param point 返回对称后的点.
	void mirror(CGeom &plane, CGeomPoint &point); //对称

	/// \brief 获取对应垂直向量.
	/// 
	void getPerpendicular(CGeomPoint &point);

	/// \brief 给定两点以及中间参数位置，求中间点
	///
	/// \param p1 第一个坐标点
	/// \param p2 第二个坐标点
	/// \param u  参数位置
	/// \return point 新的坐标点
	/// \note 参数位置需要满足<0,1>
	void getPointAtParam (CPoint &p1,CPoint &p2,CPoint &point,double &u);

	/// \brief 中点公式求坐标点
	/// \param pa 第一个坐标点
	/// \param pb 第二个坐标点
	/// \return mpoint 中点坐标
	void getMidPoint (CPoint &pa,CPoint &pb,CPoint &mpoint);

	/// \brief 求原坐标关于指定坐标位置对称点
	/// \param p 原坐标点
	/// \param sp 指定坐标位置
	/// \return newpoint 对称点
	void getSymmetryPoint (CPoint &p,CPoint &sp,CPoint &newpoint);

    double getDistance(const CPoint& p);
public:
	/// \brief 转换occ.
	/// 
	void getPnt(gp_Pnt &pnt) const; 

	/// \brief 转换occ.
	/// 
	gp_Pnt getPnt() const;
	double m_dX;	///< x坐标.
	double m_dY;	///< y坐标.
	double m_dZ;	 ///< z坐标.
private:
	inline void computePolar();
	inline void computeCartesian();

private:
	
	double m_dR;	 ///< r 极坐标.
	double m_dTheta;	///< θ极坐标.
	int m_iTypeCoord;  ///< 坐标系类别 : 1 - 笛卡尔坐标 / 0 - 极坐标/柱坐标系.
	int m_iTypeAngle;  ///< 当前角度表示 ： 1 - 度数/ 0 - 弧度.
};

#endif // !_CPOINT_H_