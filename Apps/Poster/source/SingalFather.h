#ifndef __SINGALFATHER_H__
#define __SINGALFATHER_H__
//#include "node.hpp"
#include <QObject>
#include "Common/GeometryKernel.h"
class SingalFather : public QObject {
	Q_OBJECT
public:
	static SingalFather *getInstance();
signals:
	
	/* */
	void scheckProcess(int pre_status, int now_status, bool& checked);
	void supdateGL();
	void sdrawNormalSkel();
	void genCoarseMeshFinish();
	void getData(double x, double y, double z, double angleX, double angleY, double angleZ, double length, double width, double height);

private:
	static SingalFather *instance;
	SingalFather();
};

#endif
