#ifndef __POINTCLUSTERCAL_H__
#define __POINTCLUSTERCAL_H__

#include"common.h"
#include"DistanceMap.h"
#include"SkeletonEx.h"
#include <nlopt/src/api/nlopt.h>

#define DIM_PRO 2
#define MIN_PHI 10e-2
#define INF (10000000)

class PointClusters
{
public:
	
	PointClusters(SkeletonEx& skeleton, double dmax):se(skeleton),delta_d(dmax) {
		val_num = 0 ;
		for(uint32_t i = 0 ; i<skeleton.ske_points.size(); ++i){
			if(skeleton.ske_points[i].type == SkPType::Regular){
				++val_num;
			}
		}
	};
	~PointClusters() {};
	void ExtractOrderedVList();
	
	void CalSparseCluster();
	void GetLeastSquareLine();
	void CalSymmetryAxis();
	
	inline double GetSymmetryError(uint32_t from, uint32_t to);
	inline double IntergrateSError(uint32_t index, uint32_t from, uint32_t to);
	
	double EnergyF(uint32_t branchIndex);
	double ErrorApproE(uint32_t a);
	//void DisplayField(FieldType type_field);

	//init optimizer
	void init_sparse_optimizer() ;
	//approximation error term + sparse term
    double appro_term_fun(unsigned n, const double *x, double *grad);
    //constraints 
    double boudary_constraints(unsigned n, const double *x , double *grad);

	//multi-constraints
	void multi_boudary_constraints(unsigned m, double *result, unsigned n, const double* x, double* grad);

	// lower bounds and upper bounds of (x_ori,y_ori)
    void set_lower_upper_bounds();
	//set stop threshold
	void set_stop_it(double xtol, double ftol, double force);
	
	void set_initial_val() ;
	std::vector<std::vector<uint32_t>> pos_split ;
	SkeletonEx& se;
	std::vector<std::pair<uint32_t,std::vector<Eigen::Vector2d>>> branchClusteres;
	
private:
	nlopt_opt global_opter,local_opter;
	uint32_t iteration_times = 0;
	int val_num; 
	std::vector<double> weight_w;
	std::vector<double> result_show ;
	
	std::vector<SkeletonPoint> copy_ves;
	double des_c = 0.6;//854;
	//分隔点筛选阀值
	double del_pp = 0.06;//0.08
	//点对称误差
	double delta_d = 1.00;//1.20;
	std::vector<double> delta_error;
	//inline double get_

};

#endif // !__POINTCLUSTERCAL_H__
