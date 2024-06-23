#pragma once
#ifndef _DEFINETYPE_H_
#define _DEFINETYPE_H_

enum CURVETYPE
{
	LINE = 1,			///< 枚举 线.
	BSPLINECURVE,       ///< 枚举 b样条曲线.
	BEZIERCURVE,         ///< 枚举 bezier曲线.
	CIRCLE,             ///< 枚举 Circle曲线.
	ARC                ///< 枚举 Arc曲线.
	
};

enum SURFACETYPE
{
	PLANE = 1,          ///< 枚举 平面.
	BSPLINESURFACE,     ///< 枚举 b样条曲面.
	BEZIERSURFACE,      ///< 枚举 bezier曲面.
	CYLINDRICALSURFACE,  ///< 枚举 圆柱曲面.
	LINESURFACE,        ///< 枚举 拉伸曲面.
	REVOSURFACE,        ///< 枚举 回转曲面.
	SPHERE              ///< 枚举 球体
};


#endif // !_DEFINETYPE_H_
