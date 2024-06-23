#include "BezierSolid.h"
V3f BezierSolid::p(double r, std::vector<V3f> data)
{
	return V3f((1 - r)*(1 - r)*(1 - r)*data[0].x + 3 * r*(1 - r)*(1 - r)*data[1].x + 3 * r*r*(1 - r)*data[2].x + r*r*r*data[3].x,
		(1 - r)*(1 - r)*(1 - r)*data[0].y + 3 * r*(1 - r)*(1 - r)*data[1].y + 3 * r*r*(1 - r)*data[2].y + r*r*r*data[3].y,
		(1 - r)*(1 - r)*(1 - r)*data[0].z + 3 * r*(1 - r)*(1 - r)*data[1].z + 3 * r*r*(1 - r)*data[2].z + r*r*r*data[3].z);
}

V3f BezierSolid::p1(double r, std::vector<V3f> data)
{
	return V3f((r*r - 2 * r + 1)*data[0].x + (-2 * r*r + 2 * r)*data[1].x + r * r*data[2].x, 
		(r*r - 2 * r + 1)*data[0].y + (-2 * r*r + 2 * r)*data[1].y + r * r*data[2].y
		, (r*r - 2 * r + 1)*data[0].z + (-2 * r*r + 2 * r)*data[1].z + r * r*data[2].z);
}

V3f BezierSolid::P(double u, double v, double w)
{
	std::vector<V3f> fdata;
	for (std::vector<std::vector<std::vector<V3f>>>::iterator sdata_it = SolidData.begin();sdata_it != SolidData.end();sdata_it++)
	{
		std::vector<std::vector<V3f>> temp = *sdata_it;
		std::vector<V3f> ldata;
		for (std::vector<std::vector<V3f>>::iterator temp_it = temp.begin();temp_it != temp.end();temp_it++)
		{
			ldata.push_back(p(u, *temp_it));
		}
		fdata.push_back(p(v, ldata));
	}
	return p(w, fdata);
}

V3f BezierSolid::P1(double u, double v, double w)
{
	std::vector<V3f> fdata;
	for (std::vector<std::vector<std::vector<V3f>>>::iterator sdata_it = SolidData.begin(); sdata_it != SolidData.end(); sdata_it++)
	{
		std::vector<std::vector<V3f>> temp = *sdata_it;
		std::vector<V3f> ldata;
		for (std::vector<std::vector<V3f>>::iterator temp_it = temp.begin(); temp_it != temp.end(); temp_it++)
		{
			ldata.push_back(p1(u, *temp_it));
		}
		fdata.push_back(p1(v, ldata));
	}
	return p1(w, fdata);
}