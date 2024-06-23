#include "pointcloud.h"


OBBBoundingBox_Extend::OBBBoundingBox_Extend()
{
}

OBBBoundingBox_Extend::~OBBBoundingBox_Extend()
{
}
Eigen::MatrixXd OBBBoundingBox_Extend::Compute_covariance_matrix(Eigen::MatrixXd & points)
{
	Eigen::MatrixXd covariance;
	Eigen::MatrixXd temp_pv;
	temp_pv.resize(points.rows(),points.cols());
	covariance.resize(2,2);
	covariance.setZero();
	
	//E(x)
	Eigen::VectorXd ave;
	ave.resize(2);
	ave.setZero();

	for (auto i = 0; i < points.cols(); ++i) {
		ave += points.col(i) / points.cols();
	}

	//x-E(x)
	for (auto i = 0; i < points.cols(); ++i) {
		temp_pv.col(i) = points.col(i) - ave;
	}

	//compute covariance
	for (auto xy_i = 0; xy_i < points.rows(); ++xy_i) {
		for (auto xy_j = 0; xy_j < points.rows(); ++xy_j) {
			double temp = 0.00;
			for (auto i = 0; i < points.cols(); ++i) {
				temp += temp_pv.col(i)[xy_i] * temp_pv.col(i)[xy_j];
			}
			temp /= temp_pv.cols();
			covariance(xy_i,xy_j) = temp;
		}
	}

	return covariance;
}

void OBBBoundingBox_Extend::Compute_eigen_value(Eigen::MatrixXd & cov ,Eigen::MatrixXd& points)
{
	center.resize(cov.cols());
	center.setZero();
	half_axi_len.resize(cov.cols());
	half_axi_len.assign(cov.cols(),0.00);
	//eigen_vectors.resize(cov.cols(),cov.cols());
	eigen_value.resize(cov.cols());
	eigen_value.setZero();

	Eigen::EigenSolver<Eigen::MatrixXd> es(cov);
	eigen_vectors = (es.eigenvectors().real());
	eigen_value = es.eigenvalues().real();

	for (uint32_t i = 0; i < cov.cols(); ++i) {
		
		eigen_vectors.col(i).normalize();
	}

	//Schmidt orthogonal
	//Orthonormalize(eigen_vectors);
	eigen_vectors.col(0).normalize();
	eigen_vectors.col(1).normalize();


	Eigen::Vector2d center_init(0.00,0.00);
	for (uint32_t i = 0; i < points.cols(); ++i) {
		center_init += points.col(i);
	}
	center_init = (1.00/points.cols())*center_init;

	std::vector<double> min_val(cov.cols(),MAX_LEN);
	std::vector<double> max_val(cov.cols(),-MAX_LEN);
	
	double temp_len;
	//std::cout << "\naxi:"<< eigen_vectors <<"\n";
	for (uint32_t i = 0; i < points.cols(); ++i) {
		for (uint32_t j = 0; j < points.rows(); ++j) {
			Eigen::VectorXd temp_axi = eigen_vectors.col(j);
			Eigen::VectorXd temp_p = points.col(i);
			if ((temp_len=Alingned_axi_val(temp_axi,temp_p))>max_val[j]) {
				max_val[j] = temp_len;
			}
			if ((temp_len = Alingned_axi_val(temp_axi, temp_p))<min_val[j]) {
				min_val[j] = temp_len;
			}

		}
		
	}
	//adjust the oritention of max and min value


	for (auto i = 0; i < cov.cols(); ++i) {
		
		half_axi_len[i] = (max_val[i] - min_val[i]) / 2;
		
	}


	double axi0_dis = (min_val[0] + max_val[0]) / 2;
	double axi1_dis = (min_val[1] + max_val[1]) / 2;

	center = axi0_dis * eigen_vectors.col(0)+axi1_dis*eigen_vectors.col(1);
	//if (center_init[0] * center[0] < 0) center[0] = -center[0];
	//if (center_init[1] * center[1] < 0) center[1] = -center[1];

	
	//center = center_init;
	//center = center_init;
	//std::cout << "\nhalf_axi:" << half_axi_len[0] << " "<<half_axi_len[1]<<" center:"<<center<<"\n";
}

void OBBBoundingBox_Extend::Orthonormalize(Eigen::MatrixXd & ColVecs)
{
	
	ColVecs.col(0).normalize();
	double temp;
	for (std::size_t k = 0; k != ColVecs.cols() - 1; ++k)
	{
		for (std::size_t j = 0; j != k + 1; ++j)
		{
			temp = ColVecs.col(j).transpose() * ColVecs.col(k + 1);
			ColVecs.col(k + 1) -= ColVecs.col(j) * temp;
		}
		ColVecs.col(k + 1).normalize();
	}
	
}

void OBBBoundingBox_Extend::Compute_OBB()
{
	OBB_coords.resize(eigen_vectors.rows(),4);
	std::vector<Eigen::VectorXd> tsr_coor;
	if (eigen_vectors.cols() == 2) {
		tsr_coor = {
			Eigen::Vector2d(-1.00, -1.00),
			Eigen::Vector2d(-1.00, 1.00),
			Eigen::Vector2d(1.00, 1.00),
			Eigen::Vector2d(1.00, -1.00) };

	}
	else {
		tsr_coor = {
		Eigen::Vector3d(-1.00, -1.00, -1.00),
		Eigen::Vector3d(1.00, -1.00, -1.00),
		Eigen::Vector3d(1.00, 1.00, -1.00),
		Eigen::Vector3d(-1.00, 1.00, -1.00),
		Eigen::Vector3d(-1.00, -1.00, 1.00),
		Eigen::Vector3d(1.00, -1.00, 1.00),
		Eigen::Vector3d(1.00, 1.00, 1.00),
		Eigen::Vector3d(-1.00, 1.00, 1.00) };
	}


	std::cout << "\ncenter:"<<center<<"\n";
	/*
	for (uint32_t i = 0; i < tsr_coor.size(); ++i) {
		Eigen::Vector2d temp;
		//temp.resize(tsr_coor[0].size());
		for (uint32_t j = 0; j < tsr_coor[0].size(); ++j) {
			temp(j) = center[j]+tsr_coor[i][j]*half_axi_len[j];
		}
		OBB_coords.col(i) = temp;
	}
	*/


	for (uint32_t i = 0; i < tsr_coor.size(); ++i) {
		OBB_coords.col(i) = center + tsr_coor[i][0] * half_axi_len[0] * eigen_vectors.col(0) + tsr_coor[i][1] * half_axi_len[1] * eigen_vectors.col(1);
	}


}

