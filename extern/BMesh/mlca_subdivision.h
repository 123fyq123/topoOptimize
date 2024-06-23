#ifndef __MLCA_SUBDIVISION_H__
#define __MLCA_SUBDIVISION_H__
#include "Common/GeometryKernel.h"
namespace BMesh {

class MLCASubdivision
{
public:
	MLCASubdivision() {};
	HexV3fMesh& Spliting(HexV3fMesh &mesh);
	HexV3fMesh& Averaging(HexV3fMesh &mesh);
protected:
	std::vector<V3f> new_vv;
	std::vector<V3f> new_ev;
};

}
#endif