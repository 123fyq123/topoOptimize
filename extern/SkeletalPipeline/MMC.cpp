#include "MMC.h"


double getRotateAngle(double x1, double y1, double x2, double y2)
{
	const double epsilon = 1.0e-6;

	double dist, dot, degree, angle;

	// normalize
	dist = sqrt(x1 * x1 + y1 * y1);
	x1 /= dist;
	y1 /= dist;
	dist = sqrt(x2 * x2 + y2 * y2);
	x2 /= dist;
	y2 /= dist;
	// dot product
	dot = x1 * x2 + y1 * y2;
	if (std::fabs(dot - 1.0) <= epsilon)
		angle = 0.0;
	else if (std::fabs(dot + 1.0) <= epsilon)
		angle = M_PI;
	else {
		double cross;

		angle = std::acos(dot);
		//cross product
		cross = x1 * y2 - x2 * y1;
		// vector p2 is clockwise from vector p1
		// with respect to the origin (0.0)
		if (cross < 0) {
			angle = 2 * M_PI - angle;
		}
	}
	//std::cout << angle;
	return angle;
}


void MMC_list::initialize_one_component(uint32_t i, uint32_t i_bran, uint32_t i_part, std::vector<double> ti, Eigen::Vector2d vo, double hl, Eigen::Vector2d an){
	elem_list[i].index = i ;

	elem_list[i].index_branch = i_bran ;
	
	elem_list[i].index_part = i_part ;
	//TODO: 不确定是否要assign
	elem_list[i].ti = ti ;
	elem_list[i].x0 = vo(0) ;
	elem_list[i].y0 = vo(1) ;
	elem_list[i].hl = hl;
	elem_list[i].an = getRotateAngle(1,0,an(0), an(1));
	//elem_list[i].nei_mmcs_l = left_mmc;
	//elem_list[i].nei_mmcs_r = right_mmc;
	
}

void MMC_list::initialize_components(const std::vector<uint32_t>& i_bran, const std::vector<uint32_t>& i_part, const std::vector<std::vector<double>>& ti, const std::vector<Eigen::Vector2d>& vo, const std::vector<double>& hl, const std::vector<Eigen::Vector2d>& an)
{
	for(uint32_t i = 0 ; i<num_of_mmc ; ++i){
		elem_list[i].index_branch = i_bran[i] ;
		elem_list[i].index_part = i_part[i] ;
		//TODO: 不确定是否要assign
		elem_list[i].ti = ti[i] ;
		elem_list[i].x0 = vo[i](0) ;
		elem_list[i].y0 = vo[i](1) ;
		elem_list[i].hl = hl[i];
		elem_list[i].an = getRotateAngle(1,0,an[i](0), an[i](1));
	}

}

