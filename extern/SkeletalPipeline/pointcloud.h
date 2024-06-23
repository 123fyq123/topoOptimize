#ifndef __POINT_CLOUD_H__
#define __POINT_CLOUD_H__
#include "common.h"

#define MAX_LEN 999999
class OBBBoundingBox_Extend
{
public:
	OBBBoundingBox_Extend();
	~OBBBoundingBox_Extend();

	//2d w-h 3d i-j-k
	std::vector<double> half_axi_len;
	//center
	Eigen::VectorXd center;

	//eigen values and eigen vectors
	Eigen::MatrixXd eigen_vectors;
	Eigen::VectorXd eigen_value;
	  
	//mesh nodes -1 0 1
	Eigen::VectorXi flag_nodes;

	//inner mesh nodes 
	Eigen::MatrixXd inner_nodes;

	//OBB coordinates
	Eigen::MatrixXd OBB_coords;


	Eigen::MatrixXd Compute_covariance_matrix(Eigen::MatrixXd& points);
	

	void Compute_eigen_value(Eigen::MatrixXd& cov, Eigen::MatrixXd& points);
	
	void Orthonormalize(Eigen::MatrixXd& ColVecs);
	

	void Compute_OBB();

	inline double Alingned_axi_val(Eigen::VectorXd& axi, Eigen::VectorXd& point) {
		return (axi.dot(point)/axi.norm());
	}
};



#endif // !__POINT_CLOUD_H__

