#ifndef __SHAPE_APPROXIMATE_H__
#define __SHAPE_APPROXIMATE_H__
#define Power3(c) (c*c*c) 
#include "common.h"
#include "mmc.h"
#include "PointClustersCal.h"
#include "ComponentSymmetryAxi.h"
#include <nlopt/src/api/nlopt.h>

class ShapeApproximate
{
private:
    std::string paraMMC;
    std::string meshMMC;
    uint32_t multiples;
    //对称轴
    std::vector<std::vector<SymAxi>>& axis ;
    //骨架线起、终点
    std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>>& sym_axi_ends;
    
    SkeletonEx& se;
    std::vector<std::vector<uint32_t>>& pos_split ;
    uint32_t count_x;
    nlopt_opt global_opter,local_opter;
    uint32_t iteration_times = 0 ;
    std::vector<Eigen::Vector2d> inner_p ;
    std::vector<Eigen::Vector2d> outter_p ;
    std::vector<Eigen::Vector2d> bdy_p ;

    //查找所在单元
    inline void FindCell(Eigen::Vector2d& coor, Eigen::Vector2i& p_cloest, uint32_t& p_cloest_i){
        
        p_cloest(0) = std::abs( int(coor(0)));
        p_cloest(1) = std::abs( int(coor(1)));
        p_cloest_i = se.grid_index(p_cloest(0),p_cloest(1));
        
    }
    
    double Heaviside(double x)
	{
        double ph = 0.1;
        double alpha = 1e-6;
		if (x < -ph)
		{
			return alpha;
		}
		else if (x > ph)
		{
			return 1.00;
		}
		else
		{
			return (3.00 / 4.00) * (1.00 - alpha) * ((x / ph) - (Power3(x) / (3.00 * Power3(ph)))) + ((1.00 + alpha) / 2.00);
		}
	};
    
public:
    MMC_list mmc_list ;
    ShapeApproximate(SkeletonEx& se, std::vector<std::vector<uint32_t>>& pos_split ,std::vector<std::vector<SymAxi>>& axis,std::vector<std::pair<Eigen::Vector2d,Eigen::Vector2d>>& sym_axi_ends,std::string paraMMC,std::string meshMMC,uint32_t multiples):
    se(se),pos_split(pos_split),axis(axis),sym_axi_ends(sym_axi_ends),paraMMC(paraMMC),meshMMC(meshMMC),multiples(multiples){};
    ~ShapeApproximate(){};
    
    void searchNLayers(uint32_t N) ;
    
    void InitMMCList();
    //init optimizer
	void init_sparse_optimizer() ;
	//approximation error term + sparse term
    double appro_term_fun(unsigned n, const double *x, double *grad);
    //constraints 
    double boudary_constraints(unsigned n, const double *x , double *grad);

	//multi-constraints
	void multi_boudary_constraints(unsigned m, double *result, unsigned n, const double* x, double* grad);

	// lower bounds and upper bounds of (x_ori,y_ori)
    void set_lower_upper_bounds(std::pair<double,double> x_ori, std::pair<double,double> y_ori);
	//set stop threshold
	void set_stop_it(double xtol, double ftol, double force);
	uint32_t max_iterations;
	void set_initial_val(uint32_t max_iterations) ;

};



#endif