#pragma once
#ifndef _CBezierSurface_H_
#define _CBezierSurface_H_
#include<vector>
#include<Geom_SurfaceOfLinearExtrusion.hxx>
#include<GeomAPI_PointsToBSplineSurface.hxx>
#include<GeomConvert_BSplineSurfaceToBezierSurface.hxx>
#include<Geom_BezierSurface.hxx>
#include"CBezierCurve.h"
#include "CGeomSurface.h"
#include"CPoint.h"
#include"CVector.h"
#include"DefineType.h"
using namespace std;
/// \brief  Bezier曲面类.
class CBezierSurface :
	public CGeomSurface
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:

	CBezierSurface();

	virtual ~CBezierSurface();

	/// \brief 构造无理Bezier曲面.
	/// \param controlPoints     传入控制顶点集合，集合每个方向上满足row>=2&&line>=2.
	CBezierSurface(vector<vector<CPoint>>&controlPoints);

	/// \brief 构造有理Bezier曲面.
	/// \param controlPoints     传入控制顶点集合，集合每个方向上满足row>=2&&line>=2.
	/// \param weights   传入权值集合，权值大小在【0.0-1.0】之间，长度和控制顶点长度相同.
	CBezierSurface(vector<vector<CPoint>>&controlPoints,  vector<vector<double>>&weights);

	/// \brief 插值构造多片Bezier曲面.
	/// \param points  传入用于插值的数据点.
	vector<vector<CBezierSurface>>  interPolate(vector<vector<CPoint>>&points);

	/// \brief 获取控制点集合.
	vector<vector<CPoint>> getControlPoints()const;

	/// \brief 获取权值集合.
	vector<vector<double>> getWeights()const;

	/// \brief 获取U方向控制点个数.
	int getSizeUControlPoints()const;

    /// \brief 获取V方向控制点个数.
	int getSizeVControlPoints()const;

	/// \brief 获取指定索引处的控制顶点.
	/// \param Uindex   row方向索引，应满足0<=Uindex<row方向控制点个数.
	/// \param Vindex   line方向索引，应满足0<=Vindex<line方向控制点个数.
	CPoint getPoint(const int Uindex, const int Vindex)const;

	/// \brief 获取曲面在U方向的次数.
	int getUDegree()const;

	/// \brief 获取曲面在V方向的次数.
	int getVDegree()const;

	/// \brief 获取指定索引处的权值.
	double getSingleWeight(const int Uindex, const int Vindex)const;

	/// \brief 判断U方向的有理性.
	bool getURational()const;

	/// \brief 判断V方向的有理性.
	bool getVRational()const;

	/// \brief 改变Bezier曲面上的U V方向.
	void exchangeUV();

	/// \brief 在两个参数方向上提高Bezier曲面的阶数.
	/// \param Udegree   传入U方向要提高到的次数，不小于原U方向次数.
	/// \param Vdegree   传入V方向要提高到的次数，不小于原V方向次数.
	void increaseDgree(const int Udegree, const int Vdegree);

	/// \brief 在指定索引处插入纵列控制顶点.
	/// \param Vindex  索引位置，应满足0<=Vindex<line方向控制点个数.
	/// \param point  要插入的新的控制点集合.
	///  \param type   插入位置，type==0时在Vindex处的控制点前插入，type==1时在Vindex处的控制点后插入.
	void insertVControlPoint(const int Vindex, vector<CPoint>point,int type);

	/// \brief 插入纵列指定索引处控制点以及它的权值.
    /// \param Vindex  索引位置，应满足0<=Vindex<line方向控制点个数.
	/// \param P  要插入的新的控制点集合.
	///  \param Weight 插入的新控制点的权值集合.
	/// \param type  插入位置，type==0时在Vindex处的控制点前插入，type==1时在Vindex处的控制点后插入.
	void insertVControlPoint(const int Vindex, vector<CPoint>& P, const vector<double> Weight, const int type);

	/// \brief 插入横向控制顶点.
	/// \param Uindex  索引位置，应满足0<=Uindex<row方向控制点个数.
	/// \param P  要插入的新的控制点集合.
	///  \param type   插入位置，type==0时在Vindex处的控制点前插入，type==1时在Vindex处的控制点后插入.
	void insertUControlPoint(const int Uindex, vector<CPoint>P, int type);

	/// \brief 插入横向控制点以及它的权值.
	/// \param Uindex  索引位置，应满足0<=Uindex<row方向控制点个数.
	/// \param P  要插入的新的控制点集合.
	///  \param Weight 插入的新控制点的权值集合.
	/// \param type  插入位置，type==0时在Vindex处的控制点前插入，type==1时在Vindex处的控制点后插入.
    void insertUControlPoint(const int Uindex,  vector<CPoint>P, const vector<double> Weight, const int type);

	/// \brief 移除第Vindex列控制顶点.
	/// \param Vindex  索引位置，应满足1<=Vindex<=line方向控制点个数.
	void removeVControlPoint(const int Vindex);

	/// \brief 移除第Uindex行控制顶点.
	/// \param Uindex  索引位置，应满足1<=Vindex<=row方向控制点个数.
	void removeUControlPoint(const int Uindex);

	/// \brief 分割裁剪曲面.
	/// \note  U1,U2,V1,U2应满足：0.0<U1<U2<1.0,0.0<V1<V2<1.0.
	void setTrimmed(const double U1, const double U2, const double V1, const double V2);

	/// \brief 修改单个控制点.
	/// \note  索引位置，应满足1<=Vindex(Uindex)<=row方向控制点个数.
	void setSingleControlPoint(const int Uindex, const int Vindex, const CPoint&point);

	/// \brief 修改单个控制点以及权值..
	/// \note  索引位置，应满足1<=Vindex(Uindex)<=row方向控制点个数.
	void setSingleControlPoint(const int Uindex, const int Vindex, const CPoint&point,const double weight);

	/// \brief 修改纵向控制点.
	/// \note  索引位置，应满足1<=Vindex<=row方向控制点个数.
	void setVControlPoint(const int Vindex, const vector<CPoint>P);

	/// \brief 修改纵向控制点以及权值.
	/// \note  索引位置，应满足1<=Vindex<=row方向控制点个数.
	void setVControlPoint(const int Vindex, const vector<CPoint>P,const vector<double>weights);

	/// \brief 修改横向控制点.
	/// \note  索引位置，应满足1<=(Uindex)<=row方向控制点个数.
	void setUControlPoint(const int Uindex, const vector<CPoint>P);

	/// \brief 修改横向控制点以及权值.
	/// \note  索引位置，应满足1<=(Uindex)<=row方向控制点个数.
	void setUControlPoint(const int Uindex, const vector<CPoint>P, const vector<double>weights);

	/// \brief 修改单个权值.
	/// \note  索引位置，应满足1<=Vindex(Uindex)<=row方向控制点个数.
	void setSingleWeight(const int Uindex, const int Vindex, const double weight);

	/// \brief 修改纵向权值.
	/// \note  索引位置，应满足1<=Vindex<=row方向控制点个数.
	void setVWeight(const int Vindex, vector<double>weights);

	/// \brief 修改横向权值.
	/// \note 索引位置，应满足1<=(Uindex)<=row方向控制点个数.
	void setUWeight(const int Uindex, vector<double>weights);

	/// \brief 获取U,V处的数据顶点.
	/// \note 应满足：0.0 <= U< =1.0, 0.0 <= V< =1.0.
	void getPointAtParam(const double U, const double V, CPoint&P)const;

	/// \brief 获取U,V处的顶点与一阶导.
	/// 
	/// \param U .
	/// \param V .
	/// \param vectors 返回一阶导 0索引为u向一阶导，1索引为v向一阶导.
	/// \note 应满足：0.0 <= U< =1.0, 0.0 <= V< =1.0.
	void getDerivation1(const double U, const double V, vector<CVector> &vectors)const;

	/// \brief 获取U,V处的顶点与一阶导和二阶导.
	/// 
	/// \param U .
	/// \param V .
	/// \param vectors 返回二阶导 0索引为u向二阶导，1索引为v向二阶导，2索引为uv二阶导，先偏导u，后偏导v.
	/// \note 应满足：0.0 <= U< =1.0, 0.0 <= V< =1.0.
	void getDerivation2(const double U, const double V, vector<CVector> &vectors)const;

	/// \brief 获取U,V处的顶点与一阶导、二阶导和三阶导.
	/// 
	/// \param U .
	/// \param V .
	/// \param vectors 返回三阶导 0索引为u向三阶导，1索引为v向三阶导，2索引为uuv三阶导，3索引为uvv三阶导.
	/// \note 应满足：0.0 <= U< =1.0, 0.0 <= V< =1.0.
	void getDerivation3(const double U, const double V, vector<CVector> &vectors) const;
	
	/// \brief 平移曲面.
	/// \param vec  平移向量.
	/// \note 返回旋转后曲面的新的控制点和权值.
	void setTranslation(CVector vec);

	/// \brief 旋转曲面.
	/// \param p 定位点.
	/// \param ax  旋转轴方向向量.
	/// \param ang  旋转角度（逆时针）.
	/// \note 返回旋转后曲面的新的控制点和权值.
	void setRotation(CPoint p, CVector ax, double ang);
	// 获取曲面类型.
	int getSurfaceType();
protected:
	//更新成员变量
	void TransStruct();
	//转换vector<CPoint>到TColgp_Array1OfPnt
	TColgp_Array1OfPnt transStruct(vector<CPoint>P);
	//转换vector<double>到OCC权值类型
	TColStd_Array1OfReal transStruct(vector<double>W);
	TColgp_Array2OfPnt transStruct(vector<vector<CPoint>>P);
	//传入曲面
	void setBezierSurface(Handle_Geom_BezierSurface BS);
	// 获取曲线对象
	Handle_Geom_BezierSurface getGeomSurface ();

private:

	void createBezierSurface(vector<vector<CPoint>>&controlPoints, vector<vector<double>>&weights);//创建曲面
	bool m_Urational;//u方向有理性
	bool m_Vrational;//v方向有理性
	vector<vector<CPoint>>m_ControlPoints;//控制点集合
	vector<vector<double>> m_Weights;//权值集合
	double umaxderivinv;
	double vmaxderivinv;
	bool maxderivinvok;
	Handle_Geom_BezierSurface BezierSurface;
};

#endif // !_