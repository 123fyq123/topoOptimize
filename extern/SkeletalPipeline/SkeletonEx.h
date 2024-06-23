#ifndef __SKELETON_H__
#define __SKELETON_H__
#include "common.h"
#include "DistanceMap.h"
//#include "DataVisualization.h"
//#define NUMBER_SAMPLES 60
//60
enum SkPType {
	Regular = 0,
	PerJunction = 1,
	Junction = 2,
	End=3
};
enum FieldType {
	Flag = 0,
	Distance = 1,
	AOF = 2,
	AOF_BDY = 3,
	Ske_AOF = 4,
	Ske_End = 5,
	Ske_AOF_D = 6,
	Ske_BRANCH_P = 7,
	Ske_BRANCH = 8,
	Ske_SEGLINE = 9,
	Ske_NOR_VECTOR = 10,
	Sym_POINTS = 11,
	Sym_OPT = 12,
	Sym_BROLINES = 13,
	sym_SEGLINE_CO = 14,
	Com_SYM_PS = 15,
	Com_SYM_LINES = 16,
	Com_SYM_P_VEC = 17
};



struct SkeletonPoint {
	SkPType type;
	Eigen::Vector2i coor;
	Eigen::Vector2d coord;
	uint32_t index;
	Eigen::Vector2d cloest_p;
	//std::pair<uint32_t, uint32_t> n2pIndex;
	//std::pair<Eigen::Vector2i, Eigen::Vector2i> n2p;
	//std::vector<Eigen::Vector2i> n2p_junction;
	double dis;
	std::vector<double> alphas;
	std::vector<std::pair<uint32_t,Eigen::Vector2d>> normal_brans;
	std::vector<double> rotation_ans;
};
//PreJunction-PerJunction PerJunction-Junction Junction-Junction PerJunction-end Junction-end
struct SkeletonBranch {
	uint32_t id;
	std::vector<uint32_t> vlist;
	//branch index-Junction branch index-PreJunction
	std::vector<std::pair<uint32_t,uint32_t>> neiBranch;
	
};

class SkeletonEx
{

public:
	bool readed; 
	int NUMBER_SAMPLES = 60;
	double fluxThreshold = -0.41;
	double epsilon = 1e-6;
	uint32_t width, height;
	//OpenMesh::PropertyManager<OpenMesh::FPropHandleT<int>> propsBdy;
	
	std::vector<Eigen::Vector2d> grad;
	std::vector<Eigen::Vector2i> grid_xy;
	Eigen::MatrixXi grid_index;
	//-1 bdy 0 bdyline(bdy) 1 innner
	Eigen::MatrixXi grid_value;
	Grid distance_map;
	Eigen::MatrixXd flux;
	//拓扑细化相关
	Eigen::MatrixXi newm;
	Eigen::MatrixXi del_elem;

	
	std::map<Eigen::Vector2i, uint32_t, Cmp_map_br> ske_p_index;
	//std::map<uint32_t, std::vector<uint32_t>> p_ske_index;

	//顺序：自上而下 自左向右 y-x
	inline void indexToXYZ(uint32_t index, Eigen::Vector2i& xy) {
		xy(0) = index % height;
		xy(1) = width - 1 - (index / height);
	}


	polyMesh mesh;
	std::string bdy_flag = "bdy_flag";
	std::string v_bdy_flag = "v_bdy_flag" ;
	std::vector<SkeletonBranch> branches;
	std::vector<SkeletonPoint> ske_points;
	SkeletonEx(){};
	SkeletonEx(uint32_t samplesCir, double fluxThreshold, double epsilon, uint32_t width, uint32_t height) :NUMBER_SAMPLES(samplesCir), fluxThreshold(fluxThreshold), epsilon(epsilon), width(width), height(height) {
		grid_xy.resize(width*height);
		grid_index.resize(height,width);
		grid_index.setZero();
		grid_value.resize(height,width);
		grid_value.setZero();
	};

	//read mesh with face_color flag
	//already generate the scalar field with height*width size
	void ReadTopMesh(std::string mesh_path,bool need_c);

	//init the boundary of design domain
	void InitBoundary(bool is_8_nei);

	//compute distance field in design domain
	void ComDistanceField();

	void SampleSphere2D(Eigen::MatrixXd& samples);

	//compute AOF in each grid
	void ComGradField();

	//compute R for each skeleton point
	void BindCloestPoints();
	
	//extract thick skeleton following 4 and 5 steps
	void ExThickSkeleton(bool is_8_nei);

	//homotopy thining applied in pipeline 
	void ExThinSkeleton(double pis);

	//clean and shift skeleton
	void ShiftSkeleton();

	//build each branch to SkeletonBranch
	void BuildBranches();

	//display skeleton
	void Display(uint32_t type_sk);

	//display field
	//void DisplayField(FieldType type_field);
};


#endif // !__SKELETON_H__
