#pragma once
#include <set>
#include <omp.h>
#include "MeshKernel/Mesh.h"
#include "Tools/meshmath.h"

struct MorchNode {
	FH pre_fh;// 上条交线所在的面
	EH eh;// 交点所在的边
	Vec dir;// 到达交点时的方向向量
	Vex pos;// 交点的坐标
};

struct SingularPoint {
	Vec pos;
	int handle = -1, mode = 0;// mode: 0 -> vh, 1 -> eh, 2 -> fh
};

struct Streamline {
	int beg_sid = -1, end_sid = -1;
	bool end_on_boundary = false;
	int end_edgehandle = -1;
	bool end_on_fault = false;
	std::vector<Vex> points;
	std::vector<double> length;
	std::vector<double> ratio;
	void init_length_ratio() {
		length = std::vector<double>(points.size(), 0);
		ratio = std::vector<double>(points.size(), 0);
		for (int i = 1; i < points.size(); ++i) {
			length[i] = length[i-1] + (points[i] - points[i-1]).norm();
		}
		for (int i = 1; i < points.size(); ++i) {
			ratio[i] = length[i] / length.back();
		}
	}
};

class QuadMesher {
public:
	QuadMesher(MeshKernel::SurfaceMesh _trimesh, std::vector<EH> _feature_ehs): trimesh(_trimesh), feature_ehs(_feature_ehs){};
	
	bool meshing(MeshKernel::SurfaceMesh& quadmesh);

	void get_streamlines(std::vector<Streamline>& res, std::vector<SingularPoint>& points) {
		res = streamlines;
		points = singular_points;
	}

private:
	MeshKernel::SurfaceMesh trimesh;// input parameter
	std::vector<EH> feature_ehs;

	std::vector<SingularPoint> singular_points;
	std::vector<std::vector<Vec>> crossfield;// 所有顶点的标架信息
	std::vector<Streamline> streamlines;// 所有的流线
	std::vector<EH> boundary_ehs;
	double perimeter = 0;

	std::unordered_map<VH, int> vh_to_spid;
	std::unordered_map<EH, int> eh_to_spid;
	std::unordered_map<FH, int> fh_to_spid;

	void calculate_cross_field();// 计算标架场 by Su
	void generate_initial_streamlines();// 生成流线
	void simplify_streamlines();// 流线的简化
	void quad_meshing();

	bool is_on_streamline(Vec v, EH eh, FH fh, Vec& march_dir, double& weight_vh1);// 检测是否有从奇异点 v 出发射向 eh 的流线
	bool is_on_streamline(Vec v, VH vh, EH eh, Vec& march_dir);
	void morch_on_streamline(const MorchNode&, Streamline& streamline, std::set<EH>& visited);// 生成流线
	std::vector<Vec> get_cross_dirs(const Vec&, EH eh, FH fh);// 返回该点在边上插值得到的标架场
	std::vector<Vec> get_cross_dirs(const Vec&);// 通过输入的一个标架场分量生成四个完整的标架场
	bool get_intersection_line2(const Vex& source, const Vec& dir, EH eh, Vex& intersection);// 返回二维射线与线段 eh 的交点, true 表示相交, false 表示未相交
	Vec get_most_similar_dir(const Vec& ref_dir, const std::vector<Vec>& dirs);// 返回 dirs 中与 ref_dir 方向最接近的方向向量
	std::unordered_set<EH> get_smallest_bounding_ehs(const EH& eh);// 返回这条边的最小多边形的边
	void get_streamline_itrp(const Streamline& sl0, const Streamline& sl1, Streamline& sl_new);// 由 sl0 和 sl1 插值得到新流线, 注意: sl0 须与 sl1 同向 
	void get_among_ehs(const EH& beg, const EH& end, std::set<EH>& res);// 返回连接 边beg 和 边end 之间的最小个数边界边集合(含 beg 和 end)
	void generate_initial_streamlines_v(const SingularPoint& point);// 生成流线
	void generate_initial_streamlines_e(const SingularPoint& point);// 生成流线
	void generate_initial_streamlines_f(const SingularPoint& point);// 生成流线
	bool get_instersection_polygon(const Vex& src, const Vex& dir, const EH& src_eh , const FH& src_fh, EH& intersect_eh, FH& intersect_fh, Vex& intersection);

};