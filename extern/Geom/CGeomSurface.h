#pragma once
#ifndef _CGEOMSURFACE_H_
#define _CGEOMSURFACE_H_
#include "CGeom.h"
#include "CGeomPoint.h"

/// \brief ÇúÃæ»ùÀà.
///
/// Ô¤Áô.
///
class CGeomSurface : public CGeom
{
public:
	CGeomSurface();
	virtual ~CGeomSurface();

	virtual void getNormal(CGeomPoint &) {}
	virtual void getPoint(CGeomPoint &) {}

	virtual int getSurfaceType () {
		return 0;
	}
};

#endif // !_CGEOMSURFACE_H_