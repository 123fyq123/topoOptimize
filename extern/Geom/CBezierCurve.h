#pragma once
#ifndef _CBezierCurve_H_
#define _CBezierCurve_H_
#include<Geom_BezierCurve.hxx>
#include<GeomAPI_Interpolate.hxx>
#include<GeomConvert_BSplineCurveToBezierCurve.hxx>
#include<GeomAdaptor_Curve.hxx>
#include<GCPnts_AbscissaPoint.hxx>
#include<gp_Trsf.hxx>
#include<gp_Ax1.hxx>
#include<vector>
#include<memory>
#include"CPoint.h"
#include"CGeomCurve.h"
using namespace std;

/// \brief  Bezier曲线类.
class CBezierCurve: public CGeomCurve
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:

	/// \brief 默认构造函数.
	CBezierCurve();
	/// \brief 默认析构函数.
	virtual~CBezierCurve();

	/// \brief  构造函数,创建无理Bezier曲线.
	/// \param controlPoints     传入控制顶点集合，集合长度>=2.
	CBezierCurve(vector<CPoint>&controlPoints);

	/// \brief  构造函数，创建有理Bezier曲线 .
	/// \param controlPoints     传入控制顶点集合，集合长度>=2.
	/// \param weight   传入权值集合，权值大小在【0.0-1.0】之间，长度和控制顶点长度相同.
	CBezierCurve(vector<CPoint>&controlPoints,vector<double>&weight);

	/// \brief 插值构造多段Bezier曲线.
	/// \param points  传入用于插值的数据点.
	/// \param periodic  生成曲线的周期性，true为周期性,false为非周期性.
	/// \param tolerance  插值精度，一般设置为1.0e-7~1.0e-10.
vector<CBezierCurve> interPolate(vector<CPoint>points, bool periodic, double tolerance);

	/// \brief  获取所有控制点.
	vector<CPoint> getControlPoints() const;

	/// \brief  获取权值集合.
	vector<double> getWeights() const;

	/// \brief  获取控制点个数.
	int getPointNumber();

	/// \brief  判断有理性.
	bool getRational() const;

	/// \brief返回曲线次数.
	int getDegree() const;

	/// \brief 判断曲线封闭性.
	bool getClosed() const;

	/// \brief 增加曲线次数.
	/// \param Degree   增加到的次数，应大于原曲线次数.
    void increase(const int Degree);

	/// \brief 在指定位置插入控制点.
	/// \param Index    插入控制点的位置，应满足0<Index<控制点个数.
	/// \param P  插入的控制点.
	/// \param type   插入位置   type==0时在index处的控制点前插入，type==1时在index处的控制点后插入.
	void insertControlPoint(const int Index,CPoint& P,const int type);

	/// \brief 在指定位置插入控制点以及它的权值.
	/// \param Index    插入控制点的位置，应满足0<Index<=控制点个数.
	/// \param P  插入的控制点.
	/// \param Weight  设置的权值.
	/// \param type   插入位置   type==0时在index处的控制点前插入，type==1时在index处的控制点后插入.
	void insertControlPoint(const int Index,CPoint& P, const double Weight,const int type);
	
	/// \brief 移除索引处的控制点.
	/// \param Index   删除的控制点位置，应满足0<Index<=控制点个数.
	void removeControlPoint(const int Index);

	/// \brief 分割曲线.
	/// \param U1   分割起点参数.
	/// \param U2   分割终点参数，应满足0.0<=U1<U2<=1.0.
	void setTrimmed(const double U1, const double U2);

	/// \brief  替换索引位置的控制点.
	/// \param Index   要替换的控制点位置，应满足0<Index<=控制点个数.
	/// \param P   要替换的新的控制点.
	void setControlPoint(const int Index,  CPoint& P);

	/// \brief  替换索引位置的控制点和权值.
	/// \param Index   要替换的控制点位置，应满足0<Index<=控制点个数.
	/// \param P   要替换的新的控制点.
	/// \param weight   要替换的新的控制点的权值.
	void setControlPoint(const int Index, CPoint& P,double weight);

	/// \brief  替换索引位置控制点处的权值.
    /// \param Index   要替换的控制点位置，应满足0<Index<=控制点个数.
	/// \param Weight   要设置的新的权值.
	void setWeight(const int Index, const double Weight);

	/// \brief  判断周期性.
	bool isPeriodic();

	/// \brief  获取指定参数处的顶点.
	/// \param U  参数,应满足 0.0<=U<=1.0.
	/// \param P   返回参数U处的顶点.
	void getPointAtParam(const double U, CPoint& P) const;

	/// \brief  获取指定顶点处的参数.
	/// \param P   指定的顶点位置.
	/// \param U  返回当前顶点处的参数.
	void getParamAtPoint(CPoint P, double&U);

	/// \brief  求指定参数处曲线一阶导.
	/// \param U   指定参数，应满足 0.0<=U<=1.0.
	/// \param V1  返回参数处曲线一阶导.
    void getDerivation1(const double U, CVector& V1) const;
	
	/// \brief  求指定参数处曲线二阶导.
	/// \param U   指定参数，应满足 0.0<=U<=1.0.
	/// \param V2  返回参数处曲线二阶导.
	void getDerivation2(const double U, CVector& V2) const;

	/// \brief  求指定参数处曲线三阶导.
	/// \param U   指定参数，应满足 0.0<=U<=1.0.
	/// \param V3  返回参数处曲线三阶导.
	void getDerivation3(const double U, CVector& V3) const;
	
	/// \brief  平移曲线.
	/// \param vec  平移的向量(方向和距离).
	/// \param beziercurve   返回平移后的曲线.
	void setTranslation( CVector vec,CBezierCurve&beziercurve);

	/// \brief  旋转曲线.
	/// \param p   旋转点.
	/// \param ax   旋转轴方向向量.
	/// \param ang   旋转角度.
	/// \param  beziercurve  返回旋转后的曲线.
	void setRotation(CPoint p, CVector ax, double ang,CBezierCurve&beziercurve);

	/// \brief  获取曲线起点.
	CPoint getStartPoint();

	/// \brief  获取曲线终点.
	CPoint getEndPoint();

	/// \brief  获取曲线控制点个数.
	int getNumPoint();

	/// \brief  获取曲线长度.
	double getCurveLength();

	//获取曲线类型.
	int getCurveType();
	//曲线反向
	void reverse ();
protected:
	//转换数据结构到自定义数据类型，更新成员变量
	void TransStruct();
	//转换vector<CPoint>到TColgp_HArray1OfPnt
	Handle_TColgp_HArray1OfPnt transStruct(vector<CPoint>P);
	//转换vector<double>到OCC权值类型
	TColStd_HArray1OfReal transStruct(vector<double>W);
	//设置曲线对象
	void setBezierCurve(Handle_Geom_BezierCurve beziercurve);
	//获取曲线对象
	Handle_Geom_BezierCurve getBezierCurve ()const;
private:
	void createBezierCurve(vector<CPoint>&controlPoints, vector<double>&weight);//创建曲线
	vector<CPoint> m_ControlPoints;//控制顶点集合
	bool m_Rational;//有理性
	vector<double> m_Weights;//权值集合
	bool m_Closed;//封闭性
	double maxderivinv;
	bool maxderivinvok;
	Handle_Geom_BezierCurve BezierCurve;//Bezier曲线对象
};
#endif