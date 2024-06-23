#pragma once
#ifndef _CVECTOR_H_
#define _CVECTOR_H_
#include <gp_Vec.hxx>
#include "CGeomPoint.h"
/// \brief 向量类.
///
/// 仅支持XYZ坐标系.
/// 
class CVector : public CGeomPoint
{
public:
	/// \brief 默认构造.
	/// 
	/// 默认构造一个向量.
	///
	CVector();

	/// \brief 构造.
	/// 
	/// 根据参入的参数，
	/// 构造一个向量.
	///
	/// \param x x坐标.
	/// \param y y坐标.
	/// \param z z坐标.
	/// \param typeCoord 坐标系，目前仅支持XYZ坐标系.
	CVector(double x, double y, double z, int typeCoord = 1);

	/// \brief 拷贝构造.
	/// 
	/// 构造一个和传入参数一样的
	/// 向量.
	///
	/// \param vec 被拷贝的向量.
	CVector(const CVector &vec);

	/// \brief 构造.
	/// 
	/// 将occ下的vector向量转成
	/// 当前向量.
	/// \param vec 待转换的occ向量类.
	CVector(const gp_Vec &vec);

	/// \brief 析构.
	/// 
	/// 析构向量.
	/// 
	virtual ~CVector();

public:
	/// \brief 设置/修改x坐标.
	/// 
	/// \param x x坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setX(double x);

	/// \brief 获取x坐标.
	/// 
	/// \return 返回x坐标值.
	double getX() const override;

	/// \brief 设置/修改y坐标.
	/// 
	/// \param y y坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setY(double y);

	/// \brief 获取y坐标.
	/// 
	/// \return  返回y坐标值.
	double getY() const override;

	/// \brief 设置/修改z坐标.
	/// 
	/// \param z z坐标.
	/// \return  true 设置成功，false 设置失败.
	bool setZ(double z);

	/// \brief 设置坐标系类型
	///
	/// \param type 1 xyz坐标系 | -1 无效.
	bool setTypeCoord(int type);

	/// \brief 获取z坐标.
	/// 
	/// \return  返回z坐标值.
	double getZ() const override;

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
	bool getXYZ(double &x, double &y, double &z) ;

	/// \brief 获取坐标系类型.
	/// 
	/// \return  返回向量当前使用的坐标系.
	int getTypeCoord();

public:
	/// \brief %重载.
	///
	CVector operator%( const CVector &point );
	CVector operator*( double t );
	CVector operator/( int num );
	CVector operator+( const CVector &point );
public:
	/// \brief 向量单位化.
	/// 
	/// 
	void normalize();

	/// \brief 获取对应垂直向量.
	/// 
	/// 
	void getPerpendicular(CGeomPoint &);

public:
	/*/// \brief 向量缩放.
	/// 
	/// \param t 缩放倍数.
	/// \return  返回缩放后的新向量，原向量无变化.*/


public:

	/// \brief 获取occ下的vector向量.
	/// 
	/// \return  返回occ下的vector向量.
	gp_Vec getVec();

private:
	double m_dX;          ///< x坐标.
	double m_dY;          ///< y坐标.
	double m_dZ;          ///< z坐标.
	int m_iTypeCoord;     ///< 坐标系类别，目前仅为XYZ.
};
#endif // !_CVECTOR_H_
