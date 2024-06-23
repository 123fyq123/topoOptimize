#pragma once
#include<Bnd_Box.hxx>
#include<gp_Lin.hxx>
#include<gp_Pln.hxx>
#include"CPoint.h"
#include"CLine.h"
#include"CPlane.h"
/// \brief 描述三维空间的图形的包围盒
class CBndBox {
	friend class CAidFunction;
public:
	/// \brief 生成空的包围盒
	CBndBox ();

	/// \brief 更新包围盒的范围
	void updateBox (const double aXmin,const double aYmin,const double aZmin,const double aXmax,const double aYmax,const double aZmax);

	/// \brief 获取包围盒的边界
	void getBound (double& theXmin,double& theYmin,double& theZmin,double& theXmax,double& theYmax,double& theZmax) const;

	/// \brief 获取包围盒的下边角
	CPoint getCornerMin ()const;

	/// \brief 获取包围盒的上边角
	CPoint getCornerMax ()const;

public:
	/// \brief 判断点是否在包围盒外
	bool isOut (const CPoint& point)const;

	/// \brief 判断线是否在包围盒外
	bool isOut (CLine line)const;

	/// \brief 判断平面是否在包围盒外
	bool isOut (CPlane plane)const;

	/// \brief 判断其他包围盒是否在包围盒外
	bool isOut (const CBndBox& other)const;

	/// \brief 判断两个包围盒之间的距离
	double getDistance (const CBndBox& other)const;

	/// \brief 计算包围盒对角线的长度
	double getDiagonal ()const;

protected:
	Bnd_Box mBox;
	Bnd_Box getmBndBox ();
	void setmBndBox (Bnd_Box Box);
	void updateBoxParam ();
private:
	double Xmin;
	double Xmax;
	double Ymin;
	double Ymax;
	double Zmin;
	double Zmax;
};

