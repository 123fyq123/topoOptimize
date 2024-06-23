#ifndef SUBDIVISION_H
#define SUBDIVISION_H

#include "Common/ComHeader.h"
#include "BSpline3D/BSplineSolid.h"
class Split //该类完成分割操作
{
public:
	HexV3fMesh* mesh;
	HexV3fMesh* newmesh;
	std::vector<V3f> v_crease;
	std::vector<double> material_distribution;
public:
	Split() {};
	~Split() {};
	Split(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :mesh(_mesh), newmesh(_newmesh) {};
	Split(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<V3f> _v_crease) :mesh(_mesh), newmesh(_newmesh), v_crease(_v_crease) {};
	Split(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<double> _material_distribution) :mesh(_mesh), newmesh(_newmesh)
		, material_distribution(_material_distribution) {};
	HexV3fMesh& Spliting();
	HexV3fMesh& Spliting_crease();
};

class Average  //该类完成平均操作
{
public:
	HexV3fMesh* mesh;
	std::vector<V3f> v_crease;
	std::vector<V3f> v_c_crease;
	std::vector<VertexHandle> v_singular;
	std::vector<V3f> v3f_singular;
	Average() {};
	~Average() {};
	Average(HexV3fMesh* _mesh) :mesh(_mesh) {};
	Average(HexV3fMesh* _mesh, std::vector<VertexHandle> _v_singular) :mesh(_mesh), v_singular(_v_singular) {};
	Average(HexV3fMesh* _mesh, std::vector<V3f> _v3f_singular) :mesh(_mesh), v3f_singular(_v3f_singular) {};
	Average(HexV3fMesh* _mesh, std::vector<V3f> _v_crease, std::vector<V3f> _v_c_crease) :mesh(_mesh), v_crease(_v_crease), v_c_crease(_v_c_crease) {};
public:
	HexV3fMesh& Averaging();//普通细分
	HexV3fMesh& Averaging_creases();//加入折痕（有错误）
	HexV3fMesh& Averaging_singular();//加入奇异点
	HexV3fMesh& Averaging_test();//用于测试
};

class MLCASubdivision
{
public:
	HexV3fMesh* mesh;
	HexV3fMesh* newmesh;
public:
	MLCASubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :mesh(_mesh), newmesh(_newmesh) {};
	HexV3fMesh& Spliting();
	HexV3fMesh& Averaging();
protected:
	std::vector<V3f> new_vv;
	std::vector<V3f> new_ev;
public:
	HexV3fMesh& Subdivising();
};


class Subdivision //该类用插值细分直接进行细分操作
{
public:
	HexV3fMesh* mesh;
	HexV3fMesh* newmesh;
	const double w = 1.0 / 16;
	Subdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :mesh(_mesh), newmesh(_newmesh) {};
private:
	bool check_face(FaceHandle _f);//判断是否用线型细分
	bool check_edge(EdgeHandle _e);
	V3f cal_cell_mid(const int cell_num);//局部变量不能返回引用
	V3f cal_face_mid(const int cell_num, const int face_num);
	V3f cal_edge_mid(const int cell_num, const int ver1, const int ver2);
public:
	HexV3fMesh& Subdivising();
};

class CCSubdivision//该类用CC细分
{
public:
	HexV3fMesh* mesh;
	HexV3fMesh* newmesh;
	//std::vector<bool> e_cs;//折痕边
	std::vector<bool> v_cs;//折痕点
	CCSubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :mesh(_mesh), newmesh(_newmesh) {};
	CCSubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<bool> _v_cs, std::vector<double> _material_distribution) :
		mesh(_mesh), newmesh(_newmesh), v_cs(_v_cs), material_distribution(_material_distribution) {};
	CCSubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<double> _material_distribution)
		:mesh(_mesh), newmesh(_newmesh), material_distribution(_material_distribution) {};
public:
	V3f cal_cell_mid(const int cell_num);
	V3f cal_cell_midv(CellHandle cell_id);

	V3f cal_face_mid(const int cell_num, const int face_num);
	V3f cal_face_mid(FaceHandle face_id);
	V3f cal_face_midv(FaceHandle face_id);

	V3f cal_edge_mid(const int cell_num, const int ver1, const int ver2);
	V3f cal_edge_mid_1(const int cell_num, const int ver1, const int ver2);
	V3f cal_edge_mid(EdgeHandle edge_id);
	V3f cal_edge_midv(EdgeHandle edge_id);

	V3f cal_new_vertex(const int cell_num, const int ver_num);
	V3f cal_new_vertex_1(const int cell_num, const int ver_num);
	V3f cal_new_vertex(VertexHandle ver_id);
public:
	V3f cal_new_csv(const int cell_num, const int ver_num);
	V3f cal_es_mid(const int cell_num, const int ver1, const int ver2);
	EdgeHandle check_edge(const int cell_num, const int ver1, const int ver2);
public:
	void check_cse()
	{
		//int tmp = 0;
		//e_cs.resize(mesh->edges_.size(), false);
		//for (int i = 0;i < mesh->edges_.size();i++)
		//{
		//	if (v_cs[mesh->edges_[i].from_h] == true && v_cs[mesh->edges_[i].to_h] == true)
		//	{
		//		tmp++;
		//		e_cs[i] = true;
		//	}
		//		//e_cs[i] = true;
		//}
		//mesh->e_cs = e_cs;
		mesh->v_cs = v_cs;
	}
	bool is_cs(VertexHandle _v)
	{
		return v_cs[_v];
	}
	bool is_cs(EdgeHandle _e)
	{
		if (v_cs[mesh->edges_[_e].from_h] == true && v_cs[mesh->edges_[_e].to_h] == true)
			return true;
		else
			return false;
		//return e_cs[_e];
	}
public:
	HexV3fMesh& CCSubdivising();
	HexV3fMesh& CCSubdivising_1();
	HexV3fMesh& CCSubWithCrease();//在CC细分中引入折痕
	HexV3fMesh& FirstCCSubdivising();//用于第一次细分 标记点 用于与Bezier体的比较
	HexV3fMesh& SecondCCSubdivising();//用于随后的细分 标记点 用于与Bezier体的比较
public:
	std::vector<double> material_distribution;
};

class InterpolateCC : public CCSubdivision //基于CC细分的插值细分
{
public:
	/*HexV3fMesh *ICmesh;
	HexV3fMesh *ICnewmesh;*/
	InterpolateCC(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :CCSubdivision(_mesh, _newmesh) {};
	InterpolateCC(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, HexV3fMesh* _tempmesh) :CCSubdivision(_mesh, _newmesh), tempmesh(_tempmesh) {};
	void PBfirstep();//利用CC细分得到要用的点存入firstver
	void NBfirstep();//求每个点的向量
	void MLCAfirstep();
private:
	/*-----------Push-Back法-----------*/
	std::vector<V3f> firstver;//第一步得到的点
	std::vector<V3f> detver;//第一步得到的点与原来点的差值
	/*std::vector<V3f> firstedge;//第一步得到的边点
	std::vector<V3f> detedge;//第一步得到的面点与原来边点的差值
	std::vector<V3f> firstface;//第一步得到的面点
	std::vector<V3f> detface;//第一步得到的面点与原来面点的差值
	std::vector<V3f> oriv;//原始模型的点*/
	std::vector<V3f> ev;
	std::vector<V3f> fv;
	std::vector<V3f> cv;
	std::vector<V3f> vv;
	/*const double le = 0.05;
	const double uf = 0.05;
	const double dc = 0.05;*/
	const double le = 0.25;
	const double uf = 0.125;
	const double dc = 0.0625;
	const double sv = 0.5;
	/*const double le = 0.75;
	const double uf = 0.625;
	const double dc = 0.5;*/
	/*const double le = 0.5;
	const double uf = 0.5;
	const double dc = 0.25;*/
	V3f cal_CIPB_ev(const int cell_num, const int ver1, const int ver2);
	V3f cal_CIPB_ev(EdgeHandle edge_id);
	V3f cal_CIPB_fv(const int cell_num, const int face_num);
	V3f cal_CIPB_fv(FaceHandle face_id);
	V3f cal_CIPB_cv(const int cell_num);
	V3f cal_CIPB_cv(CellHandle cell_id);
	V3f cal_CIPB_nv(const int cell_num, const int ver_num);
	V3f cal_CIPB_nv(VertexHandle ver_id);
	V3f cal_CIPB_fnv(VertexHandle ver_id);
	/*-----------Normal-Based法-----------*/
	std::vector<V3f> normalver;//存入每个点的向量
	const double we = 0.25;
	const double vf = 0.125;
	V3f cal_CINB_ev(const int cell_num, const int ver1, const int ver2);
	V3f cal_CINB_ev(EdgeHandle edge_id);
	V3f cal_CINB_fv(const int cell_num, const int face_num);
	V3f cal_CINB_fv(FaceHandle face_id);
	V3f cal_CINB_nv(const int cell_num, const int ver_num);
	/*-----------MLCA部分-----------*/
	HexV3fMesh* tempmesh;
	void Spliting();
	//void Averaging();
	/*std::vector<V3f> new_vv;
	std::vector<V3f> new_ev;
	std::vector<V3f> new_fv;
	std::vector<V3f> new_cv;*/
	std::vector<V3f> mlcadet;
	std::vector<VertexHandle> vid;
	std::vector<EdgeHandle> eid;
	std::vector<FaceHandle> fid;
	std::vector<CellHandle> cid;
	V3f cal_MLCA_ev(const int cell_num, const int ver1, const int ver2);
	V3f cal_MLCA_ev(EdgeHandle edge_id);
	V3f cal_MLCA_fv(const int cell_num, const int face_num);
	V3f cal_MLCA_fv(FaceHandle face_id);
	V3f cal_MLCA_cv(const int cell_num);
	V3f cal_MLCA_cv(CellHandle cell_id);
	V3f cal_MLCA_nv(const int cell_num, const int ver_num);
	V3f cal_mlca_vertex(VertexHandle ver_id);
public:
	HexV3fMesh& CCInterpolatingPB();
	HexV3fMesh& CCInterpolatingPB1();
	HexV3fMesh& CCInterpolatingNB();
	HexV3fMesh& CCInterpolatingMLCA();
};

class InterpolateMLCA : public MLCASubdivision
{
public:
	InterpolateMLCA(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :MLCASubdivision(_mesh, _newmesh) {};
	void firstep();//利用MLCA细分得到要用的点存入firstver
private:

};

class ApproximationCC : public CCSubdivision
{
public:
	ApproximationCC(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :CCSubdivision(_mesh, _newmesh) {};
	ApproximationCC(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<std::vector<std::vector<std::vector<V3f>>>> _SolidData)
		:CCSubdivision(_mesh, _newmesh), SolidData(_SolidData) {};
	ApproximationCC(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, std::vector<EdgeHandle> _sharpedge)
		:CCSubdivision(_mesh, _newmesh), sharpedge(_sharpedge) {};
	std::vector<std::vector<V3f>> surv;
	std::vector<V3f> nor_surv;//法向量
	std::vector<std::vector<std::vector<V3f>>> beiv;
	std::vector<std::vector<V3f>> bei_surv;
	std::vector<V3f> nor_bei_surv;
	void approximationing();
	void cal_nor_sur();
	void constructBezier();
	void cal_nor_beisur();
	std::vector<std::vector<std::vector<V3f>>> isoline;//等参线
	std::vector<EdgeHandle> sharpedge;
public:
	V3f cal_corner_point(VertexHandle _v);
	//V3f cal_corner_point_1(VertexHandle _v,CellHandle _c);
	V3f cal_edge_point(EdgeHandle _e, VertexHandle _v);
	//V3f cal_edge_point_1(EdgeHandle _e, VertexHandle _v, CellHandle _c);
	V3f cal_inter_point(FaceHandle _f, VertexHandle _v);
	V3f cal_inter_point(CellHandle _c, VertexHandle _v);
	V3f solid_edge_point(EdgeHandle _e, VertexHandle _v);
	V3f solid_edge_point_1(EdgeHandle _e, VertexHandle _v, CellHandle _c);
	V3f solid_face_point(FaceHandle _f, VertexHandle _v);
	V3f solid_corner_point(VertexHandle _v);
	V3f solid_corner_point_1(VertexHandle _v, CellHandle _c);
	V3f solid_edge_point(int cell_num, int edge_num, int vertex_num);
	V3f solid_edge_point_1(int cell_num, int edge_num, int vertex_num);
	V3f solid_face_point(int cell_num, int face_num, int vertex_num);
	/*不考虑奇异情况，进行比较*/
	V3f cal_inter_point1(FaceHandle _f, VertexHandle _v);
	V3f cal_inter_point1(CellHandle _c, VertexHandle _v);
	V3f solid_edge_point1(EdgeHandle _e, VertexHandle _v);
	V3f solid_face_point1(FaceHandle _f, VertexHandle _v);
	V3f solid_corner_point1(VertexHandle _v);
	V3f solid_edge_point1(int cell_num, int edge_num, int vertex_num);
	V3f solid_face_point1(int cell_num, int face_num, int vertex_num);
public:
	V3f cal_corner_pointsharp(VertexHandle _v);
	V3f cal_edge_pointsharp(EdgeHandle _e, VertexHandle _v);
	V3f cal_inter_pointsharp(FaceHandle _f, VertexHandle _v);
public:
	std::vector<std::vector<std::vector<std::vector<V3f>>>> SolidData;
	HexV3fMesh& CCapproximation();
	HexV3fMesh& CCapproximation1();
	HexV3fMesh& CCapproximation_1();//用于取消非流形
	HexV3fMesh constructBezierSolid();
	HexV3fMesh constructBezierSolid1();
	std::vector<HexV3fMesh>& CCapproximationDAT();
	std::vector<HexV3fMesh>& constructBezierSolidDAT();
	std::vector<HexV3fMesh> vmesh;
	std::vector<HexV3fMesh> vBmesh;
	std::vector<V3f> ccver;//存模型中每个点对应CC细分的极限点，用于比较
	//计算一个点用该方法得到的点
	void cal_ccver();
	//计算一个点对应的CC细分的极限点
	void cal_limitccver();
	/*不考虑奇异情况，进行比较*/
	HexV3fMesh& CCapproximation2();
	HexV3fMesh& CCapproximation3();
	std::vector<HexV3fMesh>& CCapproximationDAT1();
	V3f selectapoint();
	BSplineSolid bss;
	void construct_bss();//4*4*4个控制点
	void construct_bss1();//3*3*3个控制点
	void construct_bss_new();
	void construct_bss_new_new();
	void construct_bss1_new();//无需重复查找
	void construct_bss1_new_new();
};

class TwoNSubdivision : public CCSubdivision
{
public:
	TwoNSubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh) :CCSubdivision(_mesh, _newmesh) {};
	TwoNSubdivision(HexV3fMesh* _mesh, HexV3fMesh* _newmesh, int _iteration_num, std::vector<V3f> _oriv)
		:CCSubdivision(_mesh, _newmesh), iteration_num(_iteration_num), oriv(_oriv) {};
private:
	int iteration_num;
	std::vector<std::vector<V3f>> edgepoint;
	void cal_edgepoint_withn();
	V3f cal_facepoint_withn(FaceHandle _f);
	V3f cal_facepoint_withn(FaceHandle _f, int n);
	void cal_facepoint_withn_solid();
	std::vector<std::vector<V3f>> facepointsolid;
	std::vector<V3f> select_face_line(FaceHandle _f, EdgeHandle _e);//将一个方向的所有边点全部加入容器
	std::vector<V3f> oriv;//用于判断那些点是初始点
	V3f cal_cellpoint_withn(CellHandle _c);
	std::vector<V3f> select_cell_line(CellHandle _c, FaceHandle _f);//将一个方向的所有面点全部加入容器
private:
	V3f cal_twon_edge_point(const int cell_num, const int ver1, const int ver2);
	V3f cal_twon_face_point(const int cell_num, const int face_num);
	V3f cal_twon_new_vertex(const int cell_num, const int ver_num);
public:
	HexV3fMesh& TwoNSubdiving();
};

#endif // !SUBDIVISION_H

