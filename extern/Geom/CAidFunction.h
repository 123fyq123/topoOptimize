#pragma once
#ifndef _CAIDFUNCTION_H_
#define _CAIDFUNCTION_H_
#include "CGeomCurve.h"
#include "CGeomSurface.h"
#include "CPoint.h"
#include <vector>
//#include <Geom_Curve.hxx>
//#include <Geom_Surface.hxx>
#include "CBSplineCurve.h"
#include "CBezierSurface.h"
#include "CBSplineSurface.h"
#include "CBezierCurve.h"
#include "CCylindricalSurface.h"
#include"CCircle.h"
#include"CSphere.h"
#include"SurfaceOfLinearExtrusion.h"
#include"SurfaceOfRevolution.h"                          
#include"CPlane.h"
#include"CBndBox.h"
/// \brief 辅助类.
///
/// 提供一些基本几何计算，相交运算，几何文件导出.
/// 当前点的处理面向XYZ坐标系.
///
class CAidFunction
{
public:
	CAidFunction();
	virtual ~CAidFunction();
	/// \brief 计算2D下的两点点积.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \return  返回两点点积值.
	static double dot2D(CPoint &point1, CPoint &point2);

	/// \brief 计算2D下两点之间的角度.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \param typeAngle   指定返回角度类型 1 : degree  0 : radian.
	/// \return  返回两点间角度.
	static double angle2D (CPoint &point1,CPoint &point2,int typeAngle);

	/// \brief 计算2D下两点距离.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \return  返回两点距离.
	static double distance2D (CPoint &point1,CPoint &point2);

	/// \brief 计算2D下两点叉积.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \return  返回两点叉积.
	static double crossProduct2D (CPoint &point1,CPoint &point2);

	/// \brief 计算3D下的两点点积.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \return  返回两点点积值.
	static double dot3D (CGeomPoint &point1,CGeomPoint &point2);

	/// \brief 计算3D下两点之间的角度.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \param typeAngle   指定返回角度类型 1 : degree  0 : radian.
	/// \return  返回两点间角度.
	static double angle3D (CPoint &point1,CPoint &point2,int typeAngle);

	/// \brief 计算3D下两点距离.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \return  返回两点距离.
	static double distance3D (CPoint &point1,CPoint &point2);

	/// \brief 计算3D下两点叉积.
	/// 
	/// \param point1 .
	/// \param point2 .
	/// \param vec 返回两点叉积后的值.
	static void crossProduct3D (CGeomPoint &point1,CGeomPoint &point2,CGeomPoint &vec);

	/// \brief 直线平面求交.
	/// 
	/// \param line .
	/// \param plane .
	/// \return 返回t参数，用于line上求点.
	static double intersectLinePlane (CGeomCurve &line,CGeomSurface &plane);

	/// \brief 线线求交.
	/// 
	/// \param curve1 .
	/// \param curve2 .
	/// \param points 返回交点.
	/// \return true 求交成功 false 求交失败.
	/// \note line相交默认裁剪线段.
	static bool getCCIntersection (CGeomCurve &curve1,CGeomCurve &curve2,vector<CPoint> &points);

	/// \brief 线线求交.
	/// 
	/// \param curve1 .
	/// \param curve2 .
	/// \param points 返回交点.
	/// \param umin1 >= curve1.umin.
	/// \param umax1 <= curve1.umax.
	/// \param umin2 >= curve2.umin.
	/// \param umax2 <= curve2.umax.
	/// \return true 求交成功 false 求交失败.
	/// \note line直线通过参数配置，不通过裁剪线段.
	static bool getCCIntersection (CGeomCurve &curve1,CGeomCurve &curve2,vector<CPoint> &points,double umin1,double umax1,double umin2,double umax2);

	/// \brief 线面求交.
	/// 
	/// \param curve .
	/// \param surface .
	/// \param points 返回交点.
	/// \return true 求交成功 false 求交失败.
	/// \note line相交默认裁剪线段，平面默认全平面.
	static bool getCSIntersection (CGeomCurve &curve,CGeomSurface &surface,vector<CPoint> &points);

	/// \brief 面面求交.
	/// 
	/// \param surface1 .
	/// \param surface2 .
	/// \param curves 返回交线.
	/// \param tol 容差.
	/// \return true 求交成功 false 求交失败.
	/// \note 平面默认全平面|当面与面有重叠面时，底层默认两平面平行.
	static bool getSSIntersection(CGeomSurface &surface1, CGeomSurface &surface2, vector<CBSplineCurve> &curves, double tol = 1.e-10);

	/// \brief 线面投影.
	///
	/// \param surface .
	/// \param curve 待投影曲线 .
	/// \param projectCurve 返回投影线.
	/// \return true 投影成功 false 投影失败.
	/// \note 直线支持裁剪线段，平面支持全屏面.
	static bool getPorjectCurve (CGeomSurface &surface,CGeomCurve &curve,CBSplineCurve &projectCurve);

	/// \brief 合并两条曲线.
	/// \param curve1 待合并曲线1.
	/// \param curve2 待合并曲线2.
	/// \param finalCurve 最终曲线.
	/// \param tol 容差.
	static void mergeCurves (CGeomCurve &curve1,CGeomCurve &curve2,CBSplineCurve &finalCurve,double tol = 1.e-10);

	/// \brief 曲面内部抽线
	/// \param surface 输入曲面
	/// \param ORIENTATION  ==0 指定U方向，ORIENTATION==1指定V方向
	/// \param Param 指定参数位置
	/// \param curve 输出
	/// \note Param应该在参数范围内
	static bool extractCurve (CGeomSurface &surface,int ORIENTATIONC,double Param,CBSplineCurve &curve);

	/// \brief 抽取中线
	///
	/// \param curve1 边曲线1
	/// \param curve2 边曲线2
	/// \param midcurve 输出curve1和curve2物理空间中间位置的曲线
	/// \note 给定的两张曲线的位置关系要合理，否则无法抽出中间线
	static void extractmidCurve (CGeomCurve&curve1,CGeomCurve&curve2,CBSplineCurve&midcurve);
	/// \brief 抽取中面
	///
	/// \param surface1 边曲面1
	/// \param surface2 边曲面2
	/// \param midsurface 输出surface1和surface2物理空间中间位置的曲面
	/// \note 给定的两张曲面的位置关系要合理，否则无法抽出中间面

	static void extractmidSurface (CGeomSurface&surface1,CGeomSurface&surface2,CBSplineSurface&midsurface);

	/// \brief 合并两个曲面
	///
	/// \param surface1 边曲面1
	/// \param surface2 边曲面2
	/// \param mergesurface 输出surface1和surface2组成的一张曲面
	/// \note 给定的两张曲面要有共边，并且参数方向要一致

	static void mergeSurfaces (CGeomSurface &surface1,CGeomSurface &surface2,CBSplineSurface& mergesurface);

	/// \brief 曲线转换器——将几何曲线转换成B样条曲线
	///
	/// \param curve 要转换的曲线
	/// \param bsplinecurve 输出
	/// \note 支持Bezier曲线、圆、圆弧转成B样条曲线
	static void adaptorCurve(CGeomCurve &curve,CBSplineCurve&bsplinecurve);

	/// \brief 曲面转换器——将几何曲面转换成B样条曲面
	///
	/// \param surface 要转换的曲面
	/// \param bsplinesurface 输出
	/// \note 支持Bezier曲面、圆柱曲面、拉伸曲面、回转面、球面转成B样条曲面
	static void adaptorSurface (CGeomSurface &surface,CBSplineSurface &bsplinesurface);

	/// \brief 计算参数曲面的参数点的法线
	/// 
	/// \param surface 参数曲面，仅支持参数表示的曲面
	/// \param u 参数点
	/// \param v 参数点
	/// \param normal 返回的法线
	static void geomNormal (CGeomSurface &surface,double u,double v,CVector&normal);

    /// \brief 构造曲线的包围盒.
	/// 
	/// \param curve 向包围盒中添加的曲线.
	/// \param tol 精度.
	/// \param box 保存曲线集合的包围盒. 
	static void addCurveToBox (CGeomCurve& curve,double tol,CBndBox& box);

	/// \brief 构造曲面的包围盒.
	/// 
	/// \param surface 向包围盒中添加的曲面.
	/// \param tol 精度.
	/// \param box 保存曲面集合的包围盒. 
	static void addSurfaceToBox (CGeomSurface& surface,double tol,CBndBox& box);

	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param curve 待输出曲线.
	/// \note python转换未添加 . 
	static void getIGESCurve (string fileName,CGeomCurve &curve);
	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param curves 待输出曲线集. 
	/// \note python转换未添加.
	static void getIGESCurve (string fileName,vector<CBSplineCurve> curves);

	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param surface 待输出曲面. 
	/// \note python转换未添加.
	static void getIGESSurface (string fileName,CGeomSurface &surface);
	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param surfaces 待输出曲面集. 
	/// \note python转换未添加.
	static void getIGESSurface (string fileName,vector<CBSplineSurface > &surfaces);

	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param point 待输出点. 
	/// \note  不支持CVector类，python转换未添加.
	static void getIGESPoint (string fileName,CPoint &point);
	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param points 待输出点集. 
	/// \note  不支持CVector类，python转换未添加.
	static void getIGESPoint (string fileName,vector<CPoint > &points);

	/// \brief IEGS文件输出.
	/// 
	/// \param fileName 包含后缀名.
	/// \param points 待输出点集. 
	/// \param curves 待输出曲线集. 
	/// \param surfaces 待输出曲面集. 
	/// \note 不支持CVector类，python转换未添加.
	static void getIGESPntCurSur (string fileName,vector<CPoint > &points,vector<CBSplineCurve > &curves,vector<CBSplineSurface > &surfaces);
	/// \brief 读取IGES文件
	///
	/// \param fileName  文件路径/文件名，支持iges
	/// \param points 读入点集.
	/// \return true - 文件读取成功
	static bool readIGESGeomPoint(string &fileName, vector<CPoint > & points);

	/// \brief 读取IGES文件
	///
	/// \param fileName  文件路径/文件名，支持iges
	/// \param curves 读入曲线集. 
	/// \return true - 文件读取成功
	static bool readIGESGeomCurve(string &fileName, vector<CBSplineCurve > & curves);

	/// \brief 读取IGES文件
	///
	/// \param fileName  文件路径/文件名，支持iges
	/// \param surfaces 读入曲面集. 
	/// \return true - 文件读取成功
	static bool readIGESGeomSurface(string &fileName, vector<CBSplineSurface > & surfaces);

	/// \brief Dump控制台输出.
	/// 
	/// \param curve 待输出曲线. 
	static void getDumpCurve (CGeomCurve &curve);

	/// \brief Dump控制台输出.
	/// 
	/// \param surface 待输出曲面. 
	static void getDumpSurface (CGeomSurface &surface);

	/*/// \brief IEGS文件输出，返回对应的几何类.
	/// 
	/// \param filePath 待读入文件路径. 
	/// \param fileName 输出的文件名，包含后缀名.
	/// \param typeSurface 指定输出的曲面类型.
	/// \param surfaces 输出的几何曲面集合.
	/// \note 目前python暂不可用.
	static void readGeomTurbo (string filePath,string fileName,int typeSurface,vector<CBSplineSurface *> &surfaces);*/

	/*/// \brief IEGS文件输出，返回对应的几何类.
	/// 
	/// \param filePath 待读入文件路径. 
	/// \param fileName 输出的文件名，包含后缀名.
	/// \param surface1 输出的几何子午面1.
	/// \param surface2 输出的几何子午面2.
	/// \param surface3 输出的几何b样条曲面1.
	/// \param surface4 输出的几何b样条曲面2.
	static void readGeomTurbo (string filePath,string fileName,CCylindricalSurface &surface1,CCylindricalSurface &surface2,CBSplineSurface &surface3,CBSplineSurface &surface4);*/

private:
	//static void getCCIntersect(Handle_Geom_Curve &curve1, Handle_Geom_Curve &curve2, vector<CPoint> &points);
	static bool getCCIntersect (Handle_Geom_Curve &curve1,Handle_Geom_Curve &curve2,vector<CPoint> &points,double umin1,double umax1,double umin2,double umax2);
	static bool getCSIntersect (Handle_Geom_Curve &curve,Handle_Geom_Surface &surface,vector<CPoint> &points,double u1,double u2,double umin1,double umax1,double umin2,double umax2);
	//static void getCPIntersect (Handle_Geom_Curve &curve,Handle_Geom_Plane &plane,vector<CPoint> &points,double u1,double u2,double umin1,double umax1,double umin2,double umax2);

	static bool getSSIntersect (Handle_Geom_Surface &surface1,Handle_Geom_Surface &surface2,vector<CBSplineCurve> &curves,double tol = 1.e-10);
	static void getPSProject (Handle_Geom_Surface &surface,vector<CPoint> &points,vector<CPoint> &pnts);

	static bool  getIGESGeom(int type, string &fileName, vector<CPoint > points, vector<CBSplineCurve > & curves, vector<CBSplineSurface > &surfaces);
};


#endif // !_CAIDFUNCTION_H_
