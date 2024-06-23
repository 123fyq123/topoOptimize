#pragma once
#include <Geom_BSplineSurface.hxx>
#include "CGeomSurface.h"
#include "CBSplineCurve.h"
#include "nrbNurbs.h"
/// \brief b样条曲面类.
///
/// 提供一些基本b样条曲面操作.
/// 当前点的处理面向XYZ坐标系.
///
class CBSplineSurface : public CGeomSurface
{
	friend class CAidFunction;
	friend class BladeExterior;
	friend class BladeInterior;
	friend class Film_hole;
	friend class Partition;
	friend class TrailingCrack;
	friend class Turb_column;
public:
	CBSplineSurface();
	virtual ~CBSplineSurface();

	/// \brief 构造.
	/// 
	/// \param splineSurface 拷贝构造.
	CBSplineSurface(const CBSplineSurface &splineSurface);

	/// \brief 构造一条无理b样条曲面.
	/// 
	/// \param udegree u次数.
	/// \param vdegree v次数.
	/// \param controlPoints 控制点.
	/// \param uknots 节点序列（不重复的非递减序列）.
	/// \param vknots 节点序列（不重复的非递减序列）.
	/// \param umultis 记录每个uknot的重复度.
	/// \param vmultis 记录每个vknot的重复度.
	/// \param uperiodic u向周期.
	/// \param vperiodic v向周期.
	CBSplineSurface(const int &udegree, const int &vdegree, const vector<vector<CPoint>> &controlPoints, const vector<double> &uknots, const vector<double> &vknots,
		const vector<int> &umultis, const vector<int> &vmultis, bool uperiodic = false, bool vperiodic = false);
	
	/// \brief 构造一条有理b样条曲面.
	/// 
	/// \param udegree u次数.
	/// \param vdegree v次数.
	/// \param controlPoints 控制点.
	/// \param uknots 节点序列（不重复的非递减序列）.
	/// \param vknots 节点序列（不重复的非递减序列）.
	/// \param umultis 记录每个uknot的重复度.
	/// \param vmultis 记录每个vknot的重复度.
	/// \param weights 权重.
	/// \param uperiodic u向周期.
	/// \param vperiodic v向周期.
	CBSplineSurface(const int &udegree, const int &vdegree, const vector<vector<CPoint>> &controlPoints, const vector<double> &uknots, const vector<double> &vknots, 
		const vector<int> &umultis,  const vector<int> &vmultis, const vector<vector<double>> &weights, bool uperiodic = false, bool vperiodic = false);
	
	/// \brief 插值构造一条无理b样条曲面.
	/// 
	/// \param surfacePoints 插值点.
	CBSplineSurface(const vector<vector<CPoint>> &surfacePoints);

	/// \brief 插值构造一条有理b样条曲面.
	/// 
	/// \param surfacePoints 插值点.
	/// \param weights 权重.
	CBSplineSurface(const vector<vector<CPoint>> &surfacePoints, const vector<vector<double>> &weights);

	/// \brief 两条曲线构成一个曲面.
	/// 
	/// \param curve1 曲线1.
	/// \param curve2 曲线2.
	CBSplineSurface(const CBSplineCurve &curve1, const CBSplineCurve &curve2);

	/// \brief 边界曲线构成曲面
	///
	/// \param curve1 边界曲线1
	/// \param curve2 边界曲线2
	/// \param curve3 边界曲线3
	/// \param curve4 边界曲线4
	/// \note 必须保证4条边界曲线能够围成严格封闭的曲线区域（边界点相同）
	CBSplineSurface (const CBSplineCurve &curve1,const CBSplineCurve &curve2,const CBSplineCurve &curve3,const CBSplineCurve &curve4);

	/// \brief 构造固定截面的扫掠曲面。通过沿另一条曲线（路径）扫掠固定半径（截面）来构建曲面
	///
	/// \param pathcurve 扫掠过的路径曲线
	/// \param radius    固定截面半径
	CBSplineSurface ( CBSplineCurve &pathcurve,  double radius);

	/// \brief 构造指定截面曲线的扫掠曲面。通过沿另一条曲线（路径）扫掠一条曲线（截面）来构建曲面
	/// \param pathcurve 扫掠过的路径曲线
	/// \param sectioncurve  截面曲线
	/// \param sweep    ==true时才能成功放样
	/// \note 路径曲线为直线时，是沿直线进行扫掠
	CBSplineSurface ( CBSplineCurve &pathcurve, CBSplineCurve &sectioncurve,bool sweep);

public:
	/// \brief 获取b样条曲面类型.
	/// 
	/// \return  BSPLINESURFACE = 2.
	int getSurfaceType();

public:
	/// \brief 插值构造一条无理b样条曲面.
	/// 
	/// \param surfacePoints 插值点.
	/// \note 不用调用buildSurface方法,直接生成曲面，成功返回true，失败返回false.
	bool setInterpolatePoints(const vector<vector<CPoint>> &surfacePoints);

	/// \brief 获取曲面边界
	///
	/// \param curves [0]为vmin对应的(umin, umax);[0]为vmax对应的(umin, umax);[0]为umin对应的(vmin, vmax);[0]为umax对应的(vmin, vmax).
	void getBoundCurves (vector<CBSplineCurve> &curves);

	//对有洞的曲面进行曲面填充
	//bool PlateSurface (CBSplineSurface *bsplinesurface,int degree = 3);
	/// \brief 设置曲面u向次数.
	/// 
	/// \param udegree u向次数.
	/// \return true if the spline surface is not constructed, then the spline surface's udegree will be set.
	/// \note increase degree is not allowed.
	bool setUDegree(int udegree);

	/// \brief 设置曲面v向次数.
	/// 
	/// \param vdegree v向次数.
	/// \return true if the spline surface is not constructed, then the spline surface's vdegree will be set.
	/// \note increase degree is not allowed.
	bool setVDegree(int vdegree);

	/// \brief 设置/修改曲面控制点.
	/// 
	/// \param controlPoints 控制点.
	bool setControlPoints(const vector<vector<CPoint>> &controlPoints);

	/// \brief 设置/修改曲面控制点+权重.
	/// 
	/// \param controlPoints 控制点.
	/// \param weights 权重.
	bool setControlPoints(const vector<vector<CPoint>> &controlPoints, const vector<vector<double>> &weights);

	/// \brief 修改曲面单个控制点 controlPoints[u][v].
	/// \param uindex controlPoints 中的u索引,从1开始.
	/// \param vindex controlPoints 中的v索引，从1开始.
	/// \param controlPoint 新控制点值.
	bool setSingleControlPoint(int &uindex, int &vindex, const CPoint &controlPoint);

	/// \brief 修改曲面单个控制点 controlPoints[u][v].
	/// \param uindex controlPoints 中的u索引,从1开始.
	/// \param vindex controlPoints 中的v索引，从1开始.
	/// \param controlPoint 新控制点值.
	/// \param weight 权重.
	bool setSingleControlPoint(int &uindex, int &vindex, const CPoint &controlPoint, double &weight);

	/// \brief 设置/修改曲面权重.
	/// 
	/// \param weights 权重, 可以给无理b样条添加权重.
	/// \return true if the spline surface's weights are set.
	bool setWeights(const vector<vector<double>> &weights);

	/// \brief 修改曲线权重.
	/// 
	bool setSingleWeight(int uindex, int vindex, double weight); 

	/// \brief 设置/修改曲面u节点序列.
	/// 
	/// \param uknots 非递减不重复序列.
	/// \return true if the spline surface's uknots are set.
	bool setUKnots(const vector<double> &uknots);

	/// \brief 设置/修改曲面v节点序列.
	/// 
	/// \param vknots 非递减不重复序列.
	/// \return true if the spline surface's vknots are set.
	bool setVKnots(const vector<double> &vknots);

	/// \brief 修改曲面uknot.
	/// 
	/// \param uindex (1, size(uknots)).
	/// \param uknot 新的uknot值.
	bool setSingleUKnot(int uindex, double uknot);

	/// \brief 修改曲面vknot.
	/// 
	/// \param vindex (1, size(vknots)).
	/// \param vknot 新的vknot值.
	bool setSingleVKnot(int vindex, double vknot);

	/// \brief 设置/修改曲面u节点序列重复度.
	/// 
	bool setUMultis(const vector<int> &umultis);

	/// \brief 设置/修改曲面v节点序列重复度.
	/// 
	bool setVMultis(const vector<int> &vmultis);

	/// \brief 修改曲面u节点序列重复度.
	/// 
	/// \param uindex .
	/// \param multi 大于原值，最大可设置到udegree.
	bool setSingleUMultis(int uindex,int multi);

	/// \brief 修改曲面v节点序列重复度.
	/// 
	/// \param vindex .
	/// \param multi 大于原值，最大可设置到vdegree.
	bool setSingleVMultis(int vindex, int multi);

	/// \brief 设置曲面u向周期
	///
	bool setUPeriodic(bool uperiodic);

	/// \brief 设置曲面v向周期
	///
	bool setVPeriodic(bool vperiodic);

	/// \brief 插入u节点
	/// 
	/// \param u u∈[umin, umax].
	/// \param umulti .
	bool insertUKnot(double u, int umulti); 
	
	/// \brief 设置/修改u节点序列.
	/// 
	/// \param uindex from 1 - numOfUKnots.
	/// \param umulti 删除节点重复度，0表示全删.
	bool removeUKnot(int uindex, int umulti); 

	/// \brief 插入v节点
	/// 
	/// \param v v∈[vmin, vmax].
	/// \param vmulti .
	bool insertVKnot(double v, int vmulti); 

	/// \brief 设置/修改v节点序列.
	/// 
	/// \param vindex from 1 - numOfVKnots.
	/// \param vmulti 删除节点重复度，0表示全删.
	bool removeVKnot(int vindex, int vmulti);  

	/// \brief 构造曲面，当b样条为空，后期配置好各个参数，可调用这个方法构造.
	///
	bool buildSurface();

public:
	/// \brief 获取曲面控制点.
	///
	vector<vector<CPoint>> getControlPoints() const;

	/// \brief 获取曲面u控制点数目.
	///
	int getSizeUControlPoints() const;
	
	/// \brief 获取曲面v控制点数目.
	///
	int getSizeVControlPoints() const;

	/// \brief 获取曲面uknots.
	///
	vector<double> getUKnots() const;

	/// \brief 获取曲面uknots数目.
	///
	int getSizeUKnots() const;

	/// \brief 获取曲面vknots数目.
	///
	int getSizeVKnots() const;

	/// \brief 获取曲面vknots.
	///
	vector<double> getVKnots() const;

	/// \brief 获取曲面权重.
	///
	vector<vector<double>> getWeights() const;

	/// \brief 获取曲面u权重数目.
	///
	int getSizeUWeights() const;

	/// \brief 获取曲面v权重数目.
	///
	int getSizeVWeights() const;

	/// \brief 获取曲面u节点序列重复度.
	///
	vector<int>    getUMultis() const;

	/// \brief 获取曲面u节点序列重复度数目.
	///
	int getSizeUMultis() const;

	/// \brief 获取曲面v节点序列重复度.
	///
	vector<int>    getVMultis() const;

	/// \brief 获取曲面v节点序列重复度数目.
	///
	int getSizeVMultis() const;

	/// \brief 获取单个控制点 索引从1开始.
	///
	CPoint getSingleControlPoint(int uindex, int vindex); 

	/// \brief 获取单个权重 索引从1开始.
	/// 
	double getSingleWeight(int uindex, int vindex);

	/// \brief 获取单个ukont 索引从1开始.
	/// 
	double getSingleUKnot(int uindex);

	/// \brief 获取单个vkont 索引从1开始.
	/// 
	double getSingleVKnot(int vindex);

	/// \brief 获取单个umulti 索引从1开始.
	/// 
	int getSingleUMulti(int uindex);

	/// \brief 获取单个vmulti 索引从1开始.
	/// 
	int getSingleVMulti(int vindex);

	/// \brief 获取曲面面积, 当前无效.
	///
	double getArea() const; 

	/// \brief 返回曲面有理状态.
	///
	bool   getRational() const;
	
	/// \brief 获取udegree.
	///
	int    getUDegree() const;
	
	/// \brief 获取vdegree.
	///
	int    getVDegree() const;
	
	/// \brief 获取u向闭合状态 true 表示闭合.
	///
	bool   getUClosed() const; 
	
	/// \brief 获取v向闭合状态
	///
	bool   getVClosed() const;  

	/// \brief 获取u向周期.
	///
	bool   getUPeriodic() const;
	
	/// \brief 获取v向周期.
	///
	bool   getVPeriodic() const;

public:

	/// \brief 获取曲面上对应uv值点的切向量 | 一阶导.
	/// 
	vector<CVector> getDerivation1(double u, double v);

	/// \brief 获取曲面上对应uv值点的切向量 | 一阶导.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回一阶导 0索引为u向一阶导，1索引为v向一阶导.
	bool getDerivation1(double u, double v, vector<CVector> &vectors);

	/// \brief 获取曲面上对应uv值点的二阶导，曲面至少C2连续.
	/// 
	vector<CVector> getDerivation2(double u, double v);

	/// \brief 获取曲面上对应uv值点的二阶导，曲面至少C2连续.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回二阶导 0索引为u向二阶导，1索引为v向二阶导，2索引为uv二阶导，先偏导u，后偏导v.
	bool getDerivation2(double u, double v, vector<CVector> &vectors);

	/// \brief 获取曲面上对应uv值点的三阶导，曲面至少C3连续.
	/// 
	vector<CVector> getDerivation3(double u, double v);

	/// \brief 获取曲面上对应uv值点的三阶导，曲面至少C3连续.
	/// 
	/// \param u .
	/// \param v .
	/// \param vectors 返回三阶导 0索引为u向三阶导，1索引为v向三阶导，2索引为uuv三阶导，3索引为uvv三阶导.
	bool getDerivation3(double u, double v, vector<CVector> &vectors);

	/// \brief 求曲线上uv值对应的坐标点.
	///
	CPoint getPointAtParam(double u, double v);

	/// \brief 求曲线上uv值对应的坐标点.
	///
	bool getPointAtParam(double u, double v, CPoint &point);

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
	bool setTranslation(CGeomPoint & vec, CBSplineSurface &surface);
	
	/// \brief 获取旋转曲面
	///
	/// \param axisOrigin 旋转轴的起点
	/// \param axisDirection  旋转轴的方向
	/// \param angle 旋转弧度 逆时针旋转
	/// \param surface 返回旋转后的曲面，当前曲面无变化
	bool setRotation(CGeomPoint &axisOrigin, CGeomPoint &axisDirection, double angle, CBSplineSurface &surface); 

	/// \brief 裁剪曲面
	///
	/// \param u1 .
	/// \param u2   u1, u2 需要在[umin, umax], u2 > u1.
	/// \param v1 .
	/// \param v2 类似u.
	bool setTrimmed(double u1, double u2, double v1, double v2);//  u1, u2 需要在[umin, umax]内，且 u2 > u1, v类似

	/// \brief 曲面对称.
	/// \param plane 平面.
	/// \param surface 返回对称后的曲面.
	bool setMirror(CGeom &plane, CBSplineSurface &surface);

protected:
	Handle_Geom_BSplineSurface getGeomSurface() const;
	void setGeomSurface (Handle_Geom_BSplineSurface surface);

private:
	bool createBSplineSurface();
	//void computeArea();
	void updateSurfaceDate();
	void updateControlPoints();
	void updateWeights();
	void updateUKnots();
	void updateVKnots();
	void updateUMultis();
	void updateVMultis();


private:
	int m_iTypeCoord;  ///< 当前使用的坐标系 : 1 - 笛卡尔坐标(默认) / 0 - 极坐标/柱坐标系
	double m_dArea;  ///< 曲面面积
	bool m_bRational;  ///< true: 有理 
	bool m_bUPeriodic;  ///< false:非周期 
	bool m_bVPeriodic;  ///< false:非周期 
	vector<vector<CPoint>> m_ControlPoints;	///< 记录控制顶点
	vector<double> m_UKnots;   ///< u节点矢量	
	vector<double> m_VKnots;   ///< v节点矢量
	vector<vector<double>> m_Weights;	///< 权重
	vector<int>    m_UMultis;   ///< 重复度
	vector<int>    m_VMultis;   //< 重复度

	int m_iUDegree; ///< u次数
	int m_iVDegree; ///< v次数
	Handle_Geom_BSplineSurface m_pSurface;///< 曲面
};

