#include "GeometryKernel.h"
#include "math.h"
#if !defined(M_PI)
#  define M_PI 3.1415926535897932
#endif

void triangleMesh::ReorderMesh() {
	if (faces.size() == 0) return;
	flag.resize(faces.size(), 0);
	flag[0] = 1;
	for (int i = 0; i < fe[0].size(); ++i) {
		for (int j = 0; j < ef[fe[0][i]].size(); ++j) {
			if (ef[fe[0][i]][j] != 0) {
				DFS(ef[fe[0][i]][j], faces[0][i], faces[0][(i + 1) % 3]);
				break;
			}
		}
	}
}
void triangleMesh::DFS(int f, int v1, int v2) {

	if (flag[f]) return;
	flag[f] = true;
	if ((v1 == faces[f][0] && v2 == faces[f][1])) {
		faces[f][0] = v2;
		faces[f][1] = v1;
		int t = fe[f][1];
		fe[f][1] = fe[f][2];
		fe[f][2] = t;
	}
	else if ((v1 == faces[f][1] && v2 == faces[f][2])) {
		faces[f][1] = v2;
		faces[f][2] = v1;
		int t = fe[f][0];
		fe[f][0] = fe[f][2];
		fe[f][2] = t;
	}
	else if ((v1 == faces[f][2] && v2 == faces[f][0])) {
		faces[f][2] = v2;
		faces[f][0] = v1;
		int t = fe[f][0];
		fe[f][0] = fe[f][1];
		fe[f][1] = t;
	}
	for (int i = 0; i < fe[f].size(); ++i) {
		for (int j = 0; j < ef[fe[f][i]].size(); ++j) {
			if (ef[fe[f][i]][j] != f) {
				DFS(ef[fe[f][i]][j], faces[f][i], faces[f][(i + 1) % 3]);
			}

		}
	}
	return;
}

VertexHandle HexV3fMesh::add_vertex(const V3f& _p) {
	//省得麻烦不校验了
	//从1开始
	std::vector<V3f>::iterator v_id;
	//if(quick_v.find(_p) == quick_v.end())
	//if (quick_v.count(_p) == 0)
	if ((v_id = std::find(vertices_.begin(), vertices_.end(), _p)) == vertices_.end())
	{
		//quick_v[_p] = dis_v++;
		//quick_v.insert(std::pair<V3f, int>(_p, dis_v++));
		vertices_.push_back(_p);
		neighbor_v.push_back(std::set<EdgeHandle, compare_OVM>());
		neighbor_v_f.push_back(std::set<FaceHandle, compare_OVM>());
		neighbor_v_c.push_back(std::set<CellHandle, compare_OVM>());
		return VertexHandle(vertices_.size() - 1);
		//return VertexHandle(dis_v - 1);
	}
	else
	{
		return VertexHandle(abs(std::distance(vertices_.begin(), v_id)));
		//return VertexHandle(quick_v[_p]);
	}
}

FaceHandle HexV3fMesh::add_face(const std::vector<VertexHandle>& v_ve) {
	//加入面和对应边
	Face e=Face();
	e.vs_f = v_ve;
	faces_.push_back(e);
	//给每个面上边的邻接面加入引用
	return FaceHandle(faces_.size()-1);
}

CellHandle HexV3fMesh::add_cell(const std::vector<VertexHandle>& v_ve) {
	//加入边 和邻接关系 v-e
	std::vector<EdgeHandle> e_handle(12);
	for (uint16_t i = 0; i < 12; ++i) {
		Edge temp_e(v_ve[EDGE_LINK_V_INDEX[i * 2]],
			v_ve[EDGE_LINK_V_INDEX[i * 2 + 1]]);
		std::vector<Edge>::iterator it;
		if ((it = std::find(edges_.begin(), edges_.end(), temp_e)) != edges_.end()) {
			//已存在
			e_handle[i] = EdgeHandle(it - edges_.begin());
		}
		else {
			//不存在时
			edges_.push_back(temp_e);
			neighbor_e.push_back(std::set<FaceHandle, compare_OVM>());
			e_handle[i] = EdgeHandle(edges_.size()-1);
			neighbor_v[temp_e.from_h].insert(e_handle[i]);
			neighbor_v[temp_e.to_h].insert(e_handle[i]);
		}
		
	}
	//加入面 和邻接关系 f-e
	std::vector<FaceHandle> f_handle(6);
	for (uint16_t i = 0; i < 6; ++i) {
		std::vector<EdgeHandle> edge_vec = {
			e_handle[FACE_LINK_E_INDEX[i * 4]],
			e_handle[FACE_LINK_E_INDEX[i * 4 + 1]],
			e_handle[FACE_LINK_E_INDEX[i * 4 + 2]],
			e_handle[FACE_LINK_E_INDEX[i * 4 + 3]]
		};
		std::vector<VertexHandle> v_vec = {
			v_ve[FACE_LINK_V_INDEX[i * 4]],
			v_ve[FACE_LINK_V_INDEX[i * 4 + 1]],
			v_ve[FACE_LINK_V_INDEX[i * 4 + 2]],
			v_ve[FACE_LINK_V_INDEX[i * 4 + 3]]
		};
		Face face(edge_vec, v_vec);
		std::vector<Face>::iterator it;
		if ((it = std::find(faces_.begin(), faces_.end(), face)) != faces_.end()) {
			f_handle[i] = FaceHandle(it - faces_.begin());
		}
		else {
			//不存在时
			faces_.push_back(face);
			neighbor_f.push_back(std::set<CellHandle, compare_OVM>());
			f_handle[i] = FaceHandle(faces_.size()-1);
			for (uint16_t arr_i = 0; arr_i < 4; ++arr_i) {
				neighbor_e[edge_vec[arr_i]].insert(f_handle[i]);
				neighbor_v_f[v_vec[arr_i]].insert(f_handle[i]);
			}
		}
		
		
	}
	//加入单元
	cells_.push_back(Cell(f_handle, v_ve, e_handle));
	//f-c邻接关系
	for (uint16_t arr_i = 0; arr_i < 6; ++arr_i) {
		neighbor_f[f_handle[arr_i]].insert(CellHandle(cells_.size()-1));
	}
	//v_c邻接关系
	for (uint16_t arr_i = 0; arr_i < 8; ++arr_i) {
		this->neighbor_v_c[v_ve[arr_i]].insert(CellHandle(cells_.size() - 1));
	}
	return CellHandle(cells_.size() - 1);
}

CellHandle HexV3fMesh::add_cell(const std::vector<FaceHandle>& f_ve) {
	//加入单元

	//加入邻接关系 面-单元
	return CellHandle();
}

void HexV3fMesh::examine_crease()
{
	uint16_t vec_num = vertices_.size();
	for (uint16_t i = 0; i < vec_num; i++)
	{
		if (is_bdy(VertexHandle(i)))
		{
			if (neighbor_v[VertexHandle(i)].size() == 3)
			{
				v_c_crease.push_back(vertices_[VertexHandle(i)]);
			}
			if (neighbor_v[VertexHandle(i)].size() >= 4)
			{
				v_crease.push_back(vertices_[VertexHandle(i)]);
			}
		}
	}
}

bool HexV3fMesh::is_singular(EdgeHandle _e)
{
	int count = neighbor_e[_e].size();
	if (is_bdy(_e))
	{
		if (count == REGULAR_EDGE_BDY_N)
			return false;
		else
			return true;
	}
	else
	{
		if (count == REGULAR_EDGE_INNER_N)
			return false;
		else
			return true;
	}
}

bool HexV3fMesh::is_singular(VertexHandle _v)
{
	int count = neighbor_v[_v].size();
	if (is_bdy(_v))
	{
		if (count == REGULAR_VERTEX_BDY_N)
			return false;
		else
			return true;
	}
	else
	{
		if (count == REGULAR_VERTEX_INNER_N)
			return false;
		else
			return true;
	}
}

bool HexV3fMesh::is_f_in_c(FaceHandle _f, CellHandle _c)
{
	std::vector<FaceHandle> c_f = cells_[_c].faces_;
	for (std::vector<FaceHandle>::iterator c_f_it = c_f.begin();c_f_it != c_f.end();c_f_it++)
	{
		if (_f == *c_f_it)
			return true;
	}
	return false;
}

bool HexV3fMesh::is_e_in_f(EdgeHandle _e, FaceHandle _f)
{
	std::vector<EdgeHandle> f_e = faces_[_f].edges_;
	for (std::vector<EdgeHandle>::iterator f_e_it = f_e.begin();f_e_it != f_e.end();f_e_it++)
	{
		if (_e == *f_e_it)
			return true;
	}
	return false;
}

bool HexV3fMesh::is_e_in_c(EdgeHandle _e, CellHandle _c)
{
	std::vector<FaceHandle> c_f = cells_[_c].faces_;
	for (std::vector<FaceHandle>::iterator c_f_it = c_f.begin();c_f_it != c_f.end();c_f_it++)
	{
		if (is_e_in_f(_e, *c_f_it))
			return true;
	}
	return false;
}

void HexV3fMesh::examine_singular()
{
	v3f_singular.clear();
	uint16_t vec_num = vertices_.size();
	for (uint16_t i = 0; i < vec_num; i++)
	{
		std::set<EdgeHandle, compare_OVM> e_v = neighbor_v[VertexHandle(i)];
		std::vector<EdgeHandle> bd_e_v;
		if (is_bdy(VertexHandle(i)))
		{
			if (neighbor_v_c[VertexHandle(i)].size() >= 4)//边角点不算入奇异点
			{
				for (std::set<EdgeHandle, compare_OVM>::iterator e_v_it = e_v.begin();e_v_it != e_v.end();e_v_it++)
				{
					if (is_bdy(*e_v_it))
						bd_e_v.push_back(*e_v_it);
				}
				if (bd_e_v.size() != 4)
				{
					//V3f v = vertices_[VertexHandle(i)];
					v3f_singular.push_back(vertices_[VertexHandle(i)]);
				}
			}
			else
				continue;
		/*	if (neighbor_v[VertexHandle(i)].size() >= 6)
			{
				V3f v = vertices_[VertexHandle(i)];
				v3f_singular.push_back(vertices_[VertexHandle(i)]);
			}*/
		}
		else
		{
			if (neighbor_v[VertexHandle(i)].size() != 6)
			{
				//V3f v = vertices_[VertexHandle(i)];
				v3f_singular.push_back(vertices_[VertexHandle(i)]);
			}
		}
	}
}

int HexV3fMesh::val(VertexHandle v1, VertexHandle v2, bool _is_bdy)
{
	if (_is_bdy == true)
	{
		std::set<FaceHandle, compare_OVM> v1_c = neighbor_v_f[v1];
		std::set<FaceHandle, compare_OVM> v1_c_bdy;
		for (auto v1_c_it = v1_c.begin();v1_c_it != v1_c.end();v1_c_it++)
		{
			if (is_bdy(*v1_c_it))
				v1_c_bdy.insert(*v1_c_it);
		}
		std::set<FaceHandle, compare_OVM> v2_c = neighbor_v_f[v2];
		std::set<FaceHandle, compare_OVM> v2_c_bdy;
		for (auto v2_c_it = v2_c.begin();v2_c_it != v2_c.end();v2_c_it++)
		{
			if (is_bdy(*v2_c_it))
				v2_c_bdy.insert(*v2_c_it);
		}
		std::vector<FaceHandle> vc;
		std::set_intersection(v1_c_bdy.begin(), v1_c_bdy.end(), v2_c_bdy.begin(), v2_c_bdy.end(), std::inserter(vc, vc.begin()));
		//std::cout << "vc.size():" << vc.size() << std::endl;
		return vc.size();
	}
	else
	{
		std::set<CellHandle, compare_OVM> v1_c = neighbor_v_c[v1];
		std::set<CellHandle, compare_OVM> v2_c = neighbor_v_c[v2];
		std::vector<CellHandle> vc;
		std::set_intersection(v1_c.begin(), v1_c.end(), v2_c.begin(), v2_c.end(), std::inserter(vc, vc.begin()));
		//std::cout << "vc.size():" << vc.size() << std::endl;
		return vc.size();
	}
	
}

int HexV3fMesh::val(VertexHandle v , bool _is_bdy)
{
	if (_is_bdy == true)
	{
		std::set<FaceHandle, compare_OVM> v_c = neighbor_v_f[v];
		std::set<FaceHandle, compare_OVM> v_c_bdy;
		for (auto v_c_it = v_c.begin();v_c_it != v_c.end();v_c_it++)
		{
			if (is_bdy(*v_c_it))
				v_c_bdy.insert(*v_c_it);
		}
		return v_c_bdy.size();
	}
	else
	{
		std::set<CellHandle, compare_OVM> v_c = neighbor_v_c[v];
		return v_c.size();
	}
}

int HexV3fMesh::dim(VertexHandle v1, VertexHandle v2)
{
	std::set<EdgeHandle, compare_OVM> v1_e = neighbor_v[v1];
	//std::set<EdgeHandle, compare_OVM> v2_e = neighbor_v[v2];
	std::set<VertexHandle, compare_OVM> v1_e_v;
	std::set<FaceHandle, compare_OVM> v1_f = neighbor_v_f[v1];
	//std::set<FaceHandle, compare_OVM> v2_f = neighbor_v_f[v2];
	std::set<VertexHandle, compare_OVM> v1_f_v;
	std::set<CellHandle, compare_OVM> v1_c = neighbor_v_c[v1];
	//std::set<CellHandle, compare_OVM> v2_c = neighbor_v_c[v2];
	std::set<VertexHandle, compare_OVM> v1_c_v;
	for (std::set<EdgeHandle, compare_OVM>::iterator it1 = v1_e.begin(); it1 != v1_e.end();it1++)
	{
		v1_e_v.insert(edges_[*it1].from_h);
		v1_e_v.insert(edges_[*it1].to_h);
	}
	for (std::set<FaceHandle, compare_OVM>::iterator it2 = v1_f.begin(); it2 != v1_f.end();it2++)
	{
		std::vector<VertexHandle> vf = faces_[*it2].vs_f;
		for (std::vector<VertexHandle>::iterator it = vf.begin();it != vf.end();it++)
		{
			v1_f_v.insert(*it);
		}
	}
	for (std::set<CellHandle, compare_OVM>::iterator it3 = v1_c.begin();it3 != v1_c.end();it3++)
	{
		std::vector<VertexHandle> vc = cells_[*it3].vertex_;
		for (std::vector<VertexHandle>::iterator it = vc.begin();it != vc.end();it++)
		{
			v1_c_v.insert(*it);
		}
	}
	if (v1 == v2)
		return 0;
	else if (std::find(v1_e_v.begin(), v1_e_v.end(), v2) != v1_e_v.end())
		return 1;
	else if (std::find(v1_f_v.begin(), v1_f_v.end(), v2) != v1_f_v.end())
		return 2;
	else if (std::find(v1_c_v.begin(), v1_c_v.end(), v2) != v1_c_v.end())
		return 3;
	else
	{
		std::cout << "error" << std::endl;
		return 9;
	}
}

void HexV3fMesh::examine_bdfv()
{
	boundary_face_v.clear();
	boundary_face.clear();
	boundary_face_VID.clear();
	//boundary_face.clear();
	//boundary_face_v.clear();//每一次计算前先清空上一次保留的值，不然画图出错
	for (int bdfv_num = 0;bdfv_num < faces_.size();bdfv_num++)
	{
		if (is_bdy(FaceHandle(bdfv_num)))
		{
			std::vector<V3f> bdfv_1;
			std::vector<int> bdfv_id;
			bdfv_1.push_back(vertices_[faces_[bdfv_num].vs_f[0]]);
			bdfv_1.push_back(vertices_[faces_[bdfv_num].vs_f[1]]);
			bdfv_1.push_back(vertices_[faces_[bdfv_num].vs_f[2]]);
			bdfv_1.push_back(vertices_[faces_[bdfv_num].vs_f[3]]);
			bdfv_id.push_back(faces_[bdfv_num].vs_f[0]);
			bdfv_id.push_back(faces_[bdfv_num].vs_f[1]);
			bdfv_id.push_back(faces_[bdfv_num].vs_f[2]);
			bdfv_id.push_back(faces_[bdfv_num].vs_f[3]);
			boundary_face_v.push_back(bdfv_1);
			boundary_face.push_back(FaceHandle(bdfv_num));
			boundary_face_VID.push_back(bdfv_id);
		 /* V3f a;V3f b;V3f n;V3f temp;
			a.x = vertices_[faces_[bdfv_num].vs_f[1]].x - vertices_[faces_[bdfv_num].vs_f[0]].x;
			a.y = vertices_[faces_[bdfv_num].vs_f[1]].y - vertices_[faces_[bdfv_num].vs_f[0]].y;
			a.z = vertices_[faces_[bdfv_num].vs_f[1]].z - vertices_[faces_[bdfv_num].vs_f[0]].z;
			b.x = vertices_[faces_[bdfv_num].vs_f[3]].x - vertices_[faces_[bdfv_num].vs_f[0]].x;
			b.y = vertices_[faces_[bdfv_num].vs_f[3]].y - vertices_[faces_[bdfv_num].vs_f[0]].y;
			b.z = vertices_[faces_[bdfv_num].vs_f[3]].z - vertices_[faces_[bdfv_num].vs_f[0]].z;
			n.x = a.y*b.z - a.z*b.y;
			n.y = -(a.x*b.z - a.z*b.x);
			n.z = a.x*b.y - a.y*b.x;
			temp.x =  vertices_[faces_[bdfv_num].vs_f[0]].x - center.x;
			temp.y =  vertices_[faces_[bdfv_num].vs_f[0]].y - center.y;
			temp.z =  vertices_[faces_[bdfv_num].vs_f[0]].z - center.z;
			double cos;
			cos = (n.x*temp.x + n.y*temp.y + n.z*temp.z) / (n.norm()*temp.norm());
			if (cos >= 0)
			{
				normal_boundary_face_v.push_back(n);
			}
			else
			{
				normal_boundary_face_v.push_back(V3f((n.x),(n.y),(n.z)));
			}*/
		}
	}
}

void HexV3fMesh::cal_bdnorm()
{
	normal_boundary_face_v.clear();
	for (int bdn_num = 0;bdn_num < boundary_face.size();bdn_num++)
	{
		V3f n;
		n = cal_norm(boundary_face[bdn_num]);
		normal_boundary_face_v.push_back(n);
	}
}

V3f HexV3fMesh::cal_norm(FaceHandle f)
{
	V3f a; V3f b; V3f n;
	std::vector<VertexHandle> fv = faces_[f].vs_f;
	a = vertices_[fv[1]] - vertices_[fv[0]];
	b = vertices_[fv[3]] - vertices_[fv[0]];
	n = a % b;
	n = n.normalize();
	std::set<CellHandle, compare_OVM> fc = neighbor_f[f];
	CellHandle fcc = *(fc.begin());
	//使法向量向外
	std::vector<VertexHandle> fcv = cells_[fcc].vertex_;
	V3f fc_center;
	for (std::vector<VertexHandle>::iterator fcv_it = fcv.begin();fcv_it != fcv.end();fcv_it++)
	{
		fc_center = fc_center + vertices_[*fcv_it];
	}
	fc_center = fc_center / 8;
	V3f c = fc_center - vertices_[fv[0]];
	double temp = n.x*c.x + n.y*c.y + n.z*c.z;
	if (temp > 0)
		n = n * -1;
	return n;


	/*std::vector<FaceHandle> fcf = cells_[fcc].faces_;
	int i = 0;
	for (std::vector<FaceHandle>::iterator fcf_it = fcf.begin();fcf_it != fcf.end();fcf_it++)
	{
		if (*fcf_it == f)
		{
			if (i == 0 || i == 1 || i == 2 || i == 3 || i == 4)
			{
				i++;
				return n * -1;
			}
			else
			{
				i++;
				return n;
			}
		}
		else
			i++;
	}*/

	//a = vertices_[fv[1]] - vertices_[fv[0]];
	//b = vertices_[fv[3]] - vertices_[fv[0]];
	//n = a % b;
	///*V3f c = center - n;
	//double temp = n.x*c.x + n.y*c.y + n.z*c.z;
	//if (temp < 0)
	//	n = n * -1;*/
	//n = n.normalize();
	////n = n / n.norm();
	////n.x /= n.norm();n.y /= n.norm();n.z /= n.norm();//MD 原来一直错在这里 fuck
	//return n;

}

void HexV3fMesh::cal_cen()
{
	double cen_x = 0, cen_y = 0, cen_z = 0;
	double max = 0;
	double max_x = -999, max_y = -999, max_z = -999;
	double min_x = 999, min_y = 999, min_z = 999;
	int num = vertices_.size();
	for (std::vector<V3f>::iterator it = vertices_.begin();it != vertices_.end();it++)
	{
		cen_x += (*it).x;
		cen_y += (*it).y;
		cen_z += (*it).z;
	}
	center.x = cen_x / num;
	center.y = cen_y / num;
	center.z = cen_z / num;
	for (std::vector<V3f>::iterator it = vertices_.begin();it != vertices_.end();it++)
	{
		//double maxr = V3f((*it).x - center.x, (*it).y - center.y, (*it).z - center.z).norm();
		double maxr = ((*it) - center).norm();
		average_r += maxr;
		if (maxr > max)
			max = maxr;
		if ((*it).x > max_x) max_x = (*it).x;
		if ((*it).x < min_x) min_x = (*it).x;
		if ((*it).y > max_y) max_y = (*it).y;
		if ((*it).y < min_y) min_y = (*it).y;
		if ((*it).z > max_z) max_z = (*it).z;
		if ((*it).z < min_z) min_z = (*it).z;
	}
	average_r /= num;
	r = max;
	bbMax = V3f(max_x, max_y, max_z);
	bbMin = V3f(min_x, min_y, min_z);
	//double tempt = sqrt(r / 3);
	//bbMax = V3f(center.x + tempt, center.y + tempt, center.z + tempt);
	//bbMin = V3f(center.x - tempt, center.y - tempt, center.z - tempt);
}

void HexV3fMesh::cal_cell_cen()
{
	cell_center.clear();
	cell_r.clear();
	//cell_center.clear();//每一次计算前先清空上一次保留的值，不然画图出错
	for (std::vector<Cell>::iterator c_it = cells_.begin();c_it != cells_.end();c_it++)
	{
		V3f center;
		int num = (*c_it).vertex_.size();
		for (std::vector<VertexHandle>::iterator v_it = (*c_it).vertex_.begin();v_it != (*c_it).vertex_.end();v_it++)
		{
			center.x += vertices_[*v_it].x;
			center.y += vertices_[*v_it].y;
			center.z += vertices_[*v_it].z;
		}
		center.x /= num;
		center.y /= num;
		center.z /= num;
		cell_center.push_back(center);
	}
	for (int i = 0; i < cells_.size(); i++)
	{
		double temp_r = 0;
		for (int j = 0; j < cells_[i].vertex_.size(); j++)
		{
			V3f temp_v = vertices_[cells_[i].vertex_[j]];
			temp_r += sqrt(pow(cell_center[i].x - temp_v.x, 2) + pow(cell_center[i].y - temp_v.y, 2) +
				pow(cell_center[i].z - temp_v.z, 2));
		}
		temp_r /= cells_[i].vertex_.size();
		cell_r.push_back(temp_r);
	}
}

void HexV3fMesh::cal_mami_ed()
{
	double max = 0, min = 999;
	for (std::vector<Edge>::iterator e_it = edges_.begin();e_it != edges_.end();e_it++)
	{
		/*double dis = V3f(vertices_[(*e_it).from_h].x - vertices_[(*e_it).to_h].x,
			vertices_[(*e_it).from_h].y - vertices_[(*e_it).to_h].y,
			vertices_[(*e_it).from_h].z - vertices_[(*e_it).to_h].z).norm();*/
		double dis = (vertices_[(*e_it).from_h] - vertices_[(*e_it).to_h]).norm();
		if (dis > max)
			max = dis;
		if (dis < min)
			min = dis;
	}
	minEdgeLen = min;
	maxEdgeLen = max;
}

//std::vector<CellHandle> HexV3fMesh::incident_cell(FaceHandle _f)
//{
//	std::vector<CellHandle> cc;
//	int c_num = 0;
//	for (std::vector<Cell>::iterator c_it = cells_.begin();c_it != cells_.end(); c_it++)
//	{
//		std::vector<FaceHandle> tempf = (*c_it).faces_;
//		if (find(tempf.begin(), tempf.end(), _f) != tempf.end())
//		{
//			cc.push_back(CellHandle(c_num));
//		}
//		c_num++;
//	}
//	return cc;
//}

VertexHandle HexV3fMesh::find_op_vertex(VertexHandle _v, FaceHandle _f)
{
	std::vector<EdgeHandle> fe = faces_[_f].edges_;
	std::set<EdgeHandle, compare_OVM> ve = neighbor_v[_v];
	std::vector<VertexHandle> temp;
	temp.push_back(_v);
	for (std::vector<EdgeHandle>::iterator fe_it = fe.begin();fe_it != fe.end();fe_it++)
	{
		if (find(ve.begin(), ve.end(), *fe_it) != ve.end())
		{
			if (edges_[*fe_it].from_h == _v)
				temp.push_back(edges_[*fe_it].to_h);
			else
				temp.push_back(edges_[*fe_it].from_h);
		}
	}
	std::vector<VertexHandle> fv = faces_[_f].vs_f;
	for (std::vector<VertexHandle>::iterator fv_it = fv.begin();fv_it != fv.end();fv_it++)
	{
		if (find(temp.begin(), temp.end(), *fv_it) == temp.end())
			return (*fv_it);
	}
}

bool HexV3fMesh::is_two_e_in_same_face(EdgeHandle e1, EdgeHandle e2)
{
	std::set<FaceHandle> e1_f = neighbor_e_bdy(e1);
	std::set<FaceHandle> e2_f = neighbor_e_bdy(e2);
	for (std::set<FaceHandle, compare_OVM>::iterator e1_f_it = e1_f.begin();e1_f_it != e1_f.end();e1_f_it++)
	{
		for (std::set<FaceHandle, compare_OVM>::iterator e2_f_it = e2_f.begin();e2_f_it != e2_f.end();e2_f_it++)
		{
			if (*e1_f_it == *e2_f_it)
				return true;
		}
	}
	return false;
}

bool HexV3fMesh::is_two_e_in_same_cell(EdgeHandle e1, EdgeHandle e2)
{
	std::set<CellHandle> e1_c = neighbor_e_c(e1);
	std::set<CellHandle> e2_c = neighbor_e_c(e2);
	for (std::set<CellHandle>::iterator e1_c_it = e1_c.begin();e1_c_it != e1_c.end();e1_c_it++)
	{
		for (std::set<CellHandle>::iterator e2_c_it = e2_c.begin();e2_c_it != e2_c.end();e2_c_it++)
		{
			if (*e1_c_it == *e2_c_it)
				return true;
		}
	}
	return false;
}

bool HexV3fMesh::is_two_f_in_same_cell(FaceHandle f1, FaceHandle f2)
{
	std::set<CellHandle, compare_OVM> f_c1 = neighbor_f[f1];
	std::set<CellHandle, compare_OVM> f_c2 = neighbor_f[f2];
	for (std::set<CellHandle, compare_OVM>::iterator f_c1_it = f_c1.begin();f_c1_it != f_c1.end();f_c1_it++)
	{
		for (std::set<CellHandle, compare_OVM>::iterator f_c2_it = f_c2.begin();f_c2_it != f_c2.end();f_c2_it++)
		{
			if (*f_c1_it == *f_c2_it)
				return true;
		}
	}
	return false;
}

CellHandle HexV3fMesh::delete_cell(CellHandle cellHandle)
{
	//拓扑逆序删除元素

	//先删除邻接关系 再删除元素
	std::vector<FaceHandle> f_vec = cells_[cellHandle].faces_;
	std::set<EdgeHandle, compare_OVM> e_set;
	std::set<FaceHandle, compare_OVM> f_delete_set;
	std::set<EdgeHandle, compare_OVM> e_delete_set;
	std::set<VertexHandle, compare_OVM> v_delete_set;
	//面邻接
	for (std::vector<FaceHandle>::iterator f_vec_it = f_vec.begin(); f_vec_it != f_vec.end(); f_vec_it++)
	{
		FaceHandle f = *f_vec_it;
		neighbor_f[f].erase(cellHandle);
		std::vector<EdgeHandle> f_e = faces_[f].edges_;
		for (std::vector<EdgeHandle>::iterator f_e_it = f_e.begin(); f_e_it != f_e.end(); f_e_it++)
		{
			EdgeHandle e = *f_e_it;
			e_set.insert(e);
		}
		if (neighbor_f[f].size() == 0)
		{
			f_delete_set.insert(f);
			std::vector<VertexHandle> f_v = faces_[f].vs_f;
			for (std::vector<VertexHandle>::iterator f_v_it = f_v.begin(); f_v_it != f_v.end(); f_v_it++)
			{
				VertexHandle v = *f_v_it;
				neighbor_v_f[v].erase(f);
			}
		}
	}
	//边邻接
	for (std::set<EdgeHandle, compare_OVM>::iterator e_set_it = e_set.begin(); e_set_it != e_set.end(); e_set_it++)
	{
		EdgeHandle e = *e_set_it;
		if (neighbor_e[e].size() == 2)
		{
			e_delete_set.insert(e);
		}
		for (std::set<FaceHandle, compare_OVM>::iterator f_delete_set_it = f_delete_set.begin(); f_delete_set_it != f_delete_set.end();
			f_delete_set_it++)
		{
			FaceHandle f = *f_delete_set_it;
			neighbor_e[e].erase(f);
		}
	}
	//点邻接删除
	std::vector<VertexHandle> c_v = cells_[cellHandle].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		VertexHandle v = *c_v_it;
		neighbor_v_c[v].erase(cellHandle);
		std::set<EdgeHandle, compare_OVM>& edge_set = neighbor_v[v];
		//只有角点可以删掉
		if (edge_set.size() == 3) 
		{
			v_delete_set.insert(v);
		}
		for (std::set<EdgeHandle, compare_OVM>::iterator edge_iter = edge_set.begin(); edge_iter != edge_set.end();)
		{
			if (e_delete_set.find(*edge_iter) != e_delete_set.end()) 
			{
				edge_set.erase(edge_iter++);
			}
			else
			{
				++edge_iter;
			}
		}
	}
	//删除元素 c->f->e->v
	cells_[cellHandle].is_valid = false;
	for (std::set<FaceHandle, compare_OVM>::iterator f_delete_set_it = f_delete_set.begin(); f_delete_set_it != f_delete_set.end();
		f_delete_set_it++)
	{
		faces_[*f_delete_set_it].is_valid = false;
	}
	for (std::set<EdgeHandle, compare_OVM>::iterator e_delete_set_it = e_delete_set.begin(); e_delete_set_it != e_delete_set.end();
		e_delete_set_it++)
	{
		edges_[*e_delete_set_it].is_valid = false;
	}
	for (std::set<VertexHandle, compare_OVM>::iterator v_delete_set_it = v_delete_set.begin(); v_delete_set_it != v_delete_set.end();
		v_delete_set_it++)
	{
		vertices_[*v_delete_set_it].is_valid = false;
	}
	//TODO:可能换个地方比较好 太费时 放这
	examine_bdy();
	return cellHandle;
}

void HexV3fMesh::smooth_mesh()
{
	
	std::vector<V3f> newV = vertices_;
	for (int vi = 0; vi < vertices_.size(); ++vi)
	{
		V3f tempV(0, 0, 0);
		int tempn = 0;
		if (is_bdy(VertexHandle(vi)))
		{
			std::set<EdgeHandle, compare_OVM> ve = neighbor_v[vi];
			if (ve.size() == 0)
			{
				continue;
			}
			for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
			{
				if (is_bdy(*ve_it))
				{
					if (edges_[*ve_it].from_h == VertexHandle(vi))
					{
						tempV = tempV + vertices_[edges_[*ve_it].to_h];
						tempn++;
					}
					else
					{
						tempV = tempV + vertices_[edges_[*ve_it].from_h];
						tempn++;
					}
				}
			}
			if (!is_cs(VertexHandle(vi)))
			{
				tempV = tempV / tempn;
				newV[vi] = tempV;
				/*std::cout << "is_cs " <<vi << ": " << is_cs(VertexHandle(vi)) << std::endl;*/
			}
		}
		/*else
		{
			std::set<EdgeHandle, compare_OVM> ve = neighbor_v[vi];
			if (ve.size() == 0)
			{
				continue;
			}
			for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
			{
				if (edges_[*ve_it].from_h == VertexHandle(vi))
				{
					tempV = tempV + vertices_[edges_[*ve_it].to_h];
					tempn++;
				}
				else
				{
					tempV = tempV + vertices_[edges_[*ve_it].from_h];
					tempn++;
				}
			}
			tempV = tempV / tempn;
			newV[vi] = tempV;
		}*/
	}
	//vertices_ = newV;
	////内部点需要重新来一次
	//for (int vi = 0; vi < vertices_.size(); ++vi)
	//{
	//	V3f tempV(0, 0, 0);
	//	int tempn = 0;
	//	if (!is_bdy(VertexHandle(vi)))
	//	{
	//		std::set<EdgeHandle, compare_OVM> ve = neighbor_v[vi];
	//		if (ve.size() == 0)
	//		{
	//			continue;
	//		}
	//		for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	//		{
	//			if (edges_[*ve_it].from_h == VertexHandle(vi))
	//			{
	//				tempV = tempV + vertices_[edges_[*ve_it].to_h];
	//				tempn++;
	//			}
	//			else
	//			{
	//				tempV = tempV + vertices_[edges_[*ve_it].from_h];
	//				tempn++;
	//			}
	//		}
	//		tempV = tempV / tempn;
	//		newV[vi] = tempV;
	//	}
	//}
	vertices_ = newV;

}