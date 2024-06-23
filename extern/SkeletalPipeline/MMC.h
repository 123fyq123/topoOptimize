#ifndef __MMC__H__
#define __MMC__H__

#include "common.h"
//#include "DataVisualization.h"
//const static uint32_t SO = 6;
const static uint32_t SO = 18;// default: 10
const static double N_EM = 20.00;




struct MMC_2delem
{
	uint32_t index_branch;
	int32_t index_part ;
	uint32_t index;
	uint32_t is_inner;
	//ti(half length of axi), t2(...), half length of L and angle
	std::vector<double> ti;
	double hl, an;
	double x0, y0;
	
	//左右邻接组件(pair-<组件号,邻接对称轴的首点或尾点（非本对称轴）>)
    std::vector<std::pair<uint32_t,uint32_t>> nei_mmcs_l;
	std::vector<std::pair<uint32_t,uint32_t>> nei_mmcs_r;

	//有约束邻接组件
	std::set<uint32_t> cons_nei;
};

class MMC_list
{
public:
	MMC_list(){};
	MMC_list(uint32_t num_of_mmc) : num_of_mmc(num_of_mmc){
		elem_list.resize(num_of_mmc) ;
		
	};
	~MMC_list(){};
	
	double num_of_mmc;
	//bounding box - num of MMC component, compoonents, coordinates of corners, coordinates of origin point, axis, length of half axi
	
	std::vector<MMC_2delem> elem_list;

	void initialize_one_component(uint32_t i, uint32_t i_bran, uint32_t i_part, std::vector<double> ti, Eigen::Vector2d vo, double hl, Eigen::Vector2d an);

	void initialize_components(const std::vector<uint32_t>& i_bran, const std::vector<uint32_t>& i_part, const std::vector<std::vector<double>>& ti, const std::vector<Eigen::Vector2d>& vo, const std::vector<double>& hl, const std::vector<Eigen::Vector2d>& an);

	//get the width of phi(x,y) in sample point
	inline double get_f_x(uint32_t i_mmc, double x_){
		
		double return_val = 0.00;
		double temp = 1.00;
		double temp_x = 0.00;
		uint32_t n = elem_list[i_mmc].ti.size()-1;
		
		double l2 = 2*elem_list[i_mmc].hl;
		for(int i = 0 ; i<=n ; ++i){
			temp = elem_list[i_mmc].ti[i] ;
			for(int j = 0; j<=n ; ++j){
				if(i==j) continue ;
				temp_x = ((n*x_/l2)+(n/2.00)-j)/(i-j);
				
				temp=temp*temp_x ;
			}
			
			return_val+=temp ;
		}
		return return_val;
	}
	inline double get_f(uint32_t i_mmc, const Eigen::Vector2d &xy)
	{
		
		double x_ = std::cos(elem_list[i_mmc].an) * (xy(0) - elem_list[i_mmc].x0) + std::sin(elem_list[i_mmc].an) * (xy(1) - elem_list[i_mmc].y0);
		
		double return_val = 0.00;
		double temp = 1.00;
		double temp_x = 0.00;
		uint32_t n = elem_list[i_mmc].ti.size()-1;
		
		double l2 = 2*elem_list[i_mmc].hl;
		for(int i = 0 ; i<=n ; ++i){
			temp = elem_list[i_mmc].ti[i] ;
			for(int j = 0; j<=n ; ++j){
				if(i==j) continue ;
				temp_x = ((n*x_/l2)+(n/2.00)-j)/(i-j);
				
				temp=temp*temp_x ;
			}
			
			return_val+=temp ;
		}
		return return_val;
	}

	//非递归快速幂
	inline double qpow(double a, long long n)
	{
		double ans = 1; // 赋值为乘法单位元，可能要根据构造函数修改
		while (n)
		{
			if (n & 1)
				ans = ans * a; 
			n >>= 1;
			a = a * a;
		}
		return ans;
	}
	//get the value of MMC function phi(x,y) in sample point
	inline double get_nord_phi_xy(uint32_t i, const Eigen::Vector2d &xy)
	{
		double x_ = std::cos(elem_list[i].an) * (xy[0] - elem_list[i].x0) + std::sin(elem_list[i].an) * (xy[1] - elem_list[i].y0);
		double y_ = -std::sin(elem_list[i].an) * (xy[0] - elem_list[i].x0) + std::cos(elem_list[i].an) * (xy[1] - elem_list[i].y0);
		double fx = get_f(i, xy);
		//std::cout<<"i:"<<i<<"x_:"<<x_<<"y_:"<<y_<<"fx:"<<fx <<"\n";
		return 1.00 - qpow(x_ / elem_list[i].hl, SO) - qpow(y_ / fx, SO);
	}

	//快速exp
	inline double fa_exp2(double x) {
		x = 1.0 + x / 1024;
		x *= x; x *= x; x *= x; x *= x;
		x *= x; x *= x; x *= x; x *= x;
		x *= x; x *= x;
		return x;
	}
	

	//get the value of EM(F) in sample point2
	inline double get_EM_xy(const Eigen::Vector2d &xy)
	{
		double max_val = -9999999;
		for (uint32_t i = 0; i < elem_list.size(); ++i)
		{
			double maxt = (*(std::max_element(elem_list[i].ti.begin(), elem_list[i].ti.end()))) + 3.00;
			if ((xy - Eigen::Vector2d(elem_list[i].x0, elem_list[i].y0)).squaredNorm() > maxt * maxt + elem_list[i].hl * elem_list[i].hl + 20) continue;

			double cck = get_nord_phi_xy(i, xy);
			if (cck > max_val) max_val = cck;
		}
		if (max_val < 0) return 0;
		return max_val;
	}
	inline void DisplayBoundary(double sampleXrot,std::vector<std::vector<Eigen::Vector2d>>& boundaryLines){
		boundaryLines.resize(elem_list.size()) ;
		for(uint32_t i = 0; i<elem_list.size(); ++i){
			uint32_t samples = (int)(2*elem_list[i].hl/sampleXrot);
			double xrot = 2*elem_list[i].hl/samples;
			boundaryLines[i].resize(2*(samples+1)+1) ;
			for(uint32_t xi = 0; xi<samples+1; ++xi){
				double xc =  -elem_list[i].hl + xi*xrot;
				double yc = -get_f_x(i,xc);
				Eigen::Vector2d xy_convert;
				xy_convert(0) = elem_list[i].x0+std::cos(elem_list[i].an)*xc-std::sin(elem_list[i].an)*yc;
				xy_convert(1) = elem_list[i].y0+std::sin(elem_list[i].an)*xc+std::cos(elem_list[i].an)*yc;
				boundaryLines[i][xi] = xy_convert;
			}
			for(uint32_t xi = 0; xi<samples+1; ++xi){
				double xc =  elem_list[i].hl - xi*xrot;
				double yc = get_f_x(i,xc);
				Eigen::Vector2d xy_convert;
				xy_convert(0) = elem_list[i].x0+std::cos(elem_list[i].an)*xc-std::sin(elem_list[i].an)*yc;
				xy_convert(1) = elem_list[i].y0+std::sin(elem_list[i].an)*xc+std::cos(elem_list[i].an)*yc;
				boundaryLines[i][xi+samples+1] = xy_convert;
			}
			boundaryLines[i][2*(samples+1)] = boundaryLines[i][0];		
		}
	}
	/*
	inline void display_all(){
		DataVisualization d;
		std::vector<std::vector<double>> x,y,z ;
		
		for(double xi = -802.00; xi<0.00; xi+=0.5){
			std::vector<double> x_row,y_row,z_row;
			for(double yi = -272.00; yi<0.00; yi+=0.5){
				x_row.push_back(xi) ;
				y_row.push_back(yi) ;
				//get_EM_xy(Eigen::Vector2d(x,y))
				z_row.push_back(get_EM_xy(Eigen::Vector2d(xi,yi))) ;
				
			}
			x.push_back(x_row);
			y.push_back(y_row);
			z.push_back(z_row);
		}
		d.ShowSurface("mmc",x,y,z) ;
	}*/
	//写出mmc文件
	inline void write_to_file(std::string filename){
		std::fstream f(filename,std::ios::out) ;
		int size_of_til = 0;
		int size_of_tir = 0;
		f<<num_of_mmc<<"\n";
		for(uint32_t i = 0; i<num_of_mmc; ++i){
			size_of_til = elem_list[i].ti.size();
			size_of_tir = elem_list[i].ti.size();
			f<<elem_list[i].index<<" "<<elem_list[i].is_inner<<" "<<size_of_til<<" "<<size_of_tir<<" "<<elem_list[i].nei_mmcs_l.size()<<" "<<elem_list[i].nei_mmcs_r.size()<<" "<<elem_list[i].cons_nei.size()<<"\n";
			f<<elem_list[i].x0<<" "<<elem_list[i].y0<<" "<<elem_list[i].hl<<" "<<elem_list[i].an<<"\n";
			for(uint32_t lr = 0 ; lr<2; ++lr){
				for(uint32_t tii = 0; tii<elem_list[i].ti.size(); ++tii){
					f<<elem_list[i].ti[tii]<<" ";
				}
				f<<"\n";
			}
			for(uint32_t ci = 0; ci<elem_list[i].nei_mmcs_l.size(); ++ci){
				f<<elem_list[i].nei_mmcs_l[ci].first<<" "<<elem_list[i].nei_mmcs_l[ci].second;
				if(ci!=elem_list[i].nei_mmcs_l.size()-1){
					f<<" ";
				}else{
					f<<"\n";
				}
			}
			for(uint32_t ci = 0; ci<elem_list[i].nei_mmcs_r.size(); ++ci){
				f<<elem_list[i].nei_mmcs_r[ci].first<<" "<<elem_list[i].nei_mmcs_r[ci].second;
				if(ci!=elem_list[i].nei_mmcs_r.size()-1){
					f<<" ";
				}else{
					f<<"\n";
				}
			}
			for(auto ci = elem_list[i].cons_nei.begin(); ci!=elem_list[i].cons_nei.end(); ++ci){
				f<<*ci<<" ";
			}
			if(elem_list[i].cons_nei.size() != 0)
				f<<"\n";
		}
		f.close() ;
	}
	inline void write_vtk(std::string filename,Eigen::MatrixXd& V){
		std::fstream f(filename, std::ios::out);

		//f << "# vtk DataFile Version 2.0" << std::endl << "hex mesh vtk data" << std::endl;
		//f << "ASCII" << std::endl;
		//f << "DATASET UNSTRUCTURED_GRID" << std::endl;

		//f << "POINTS " << V.cols() << " double" << std::endl;
		f << "x,y,p" << std::endl;
		for (uint32_t i = 0; i<V.cols(); i++)
			f << V(0, i) << "," << V(1, i) <<","<<((V(2,i)>0)?1:0) << std::endl;

		//f << "POINT_DATA " << V.cols() << std::endl;
		//f << "SCALARS fixed int" << std::endl;
		//f << "LOOKUP_TABLE default" << std::endl;
		//for (uint32_t i = 0; i<V.cols(); i++)
		//	f << V(2,i) << std::endl;

		f.close();
	}
};

#endif // !__MMC__H__
