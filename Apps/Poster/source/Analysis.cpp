#include "TopOpt.h"
#include "Bdy.h"


void Analysis(ModelDraw* widget, std::string _filename, std::vector<std::vector<int>>& _pickface, std::vector<vector<int>>& _pickedge, 
	std::vector<std::vector<double>>& x, std::vector<std::vector<double>>& y, std::vector<std::vector<double>>& z)
{

	if (_pickedge.empty() || _pickface.empty()) {
		std::cerr << "Pick data is empty!" << std::endl;
		return;
	}

	/*-----------------读取文件-----------------*/
	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(_filename);

	widget->loadHexMesh(_filename);


	file.read2HexMeshObj(mesh, true, true, true);       //读取初始文件

	mesh->examine_bdy();
	mesh->cal_cell_cen();
	mesh->cal_cen();
	std::cout << "INIT MESH OK" << std::endl;
	EdgeHandle e(1);

	std::vector<VertexHandle> pickvertex;
	std::vector<EdgeHandle> pickedge;
	std::vector<FaceHandle> pickface;
	std::vector<CellHandle> pickcell;
	for (int i = 0; i < _pickface.size(); ++i) {
		std::sort(_pickface[i].begin(), _pickface[i].end());
		for (int j = 0; j < mesh->faces_.size(); ++j) {
			auto fv = mesh->faces_[j].vs_f;
			std::sort(fv.begin(), fv.end());
			bool flag = true;
			for (int k = 0; k < fv.size(); ++k) {
				if (fv[k].idx() != _pickface[i][k]) {
					flag = false;
					break;
				}
			}
			if (flag) { pickface.push_back((FaceHandle)j); break; }
		}
	}
	for (int i = 0; i < _pickedge.size(); ++i) {
		std::sort(_pickedge[i].begin(), _pickedge[i].end());
		for (int j = 0; j < mesh->edges_.size(); ++j) {
			std::vector<VertexHandle> ev{ mesh->edges_[j].from_h, mesh->edges_[j].to_h };
			std::sort(ev.begin(), ev.end());
			bool flag = true;
			for (int k = 0; k < ev.size(); ++k) {
				if (ev[k].idx() != _pickedge[i][k]) {
					flag = false;
					break;
				}
			}
			if (flag) { pickedge.push_back((EdgeHandle)j); break; }
		}
	}
	
	/*-----------------获取边界条件csv-----------------*/
	mesh->v_cs.resize(mesh->vertices_.size(), false);                    //折痕
	/*-----------------获取边界条件sharp-----------------*/
	std::vector<V3f> edgemid;
	std::vector<EdgeHandle> edgesharp;
	std::cout << edgesharp.size() << std::endl;
	std::cout << mesh->v_cs.size() << std::endl;
	


	HexV3fMesh* newmesh = new HexV3fMesh();
	mesh->cal_cell_cen();
	mesh->cal_cen();

	/*-----------------下一步的细分网格-----------------*/
	
	std::cout << "Start building BSplineSolid data structure" << std::endl;
	ApproximationCC AppCCSIMP(mesh, newmesh, edgesharp);          //构造bezier体逼近
	std::cout << "Start building BSplineSolid data structure" << std::endl;
	AppCCSIMP.construct_bss_new();
	BSplineSolid m_solids = AppCCSIMP.bss;
	std::cout << "Building BSplineSolid data structure completely" << std::endl;


	std::vector<std::vector<double>> local_bigrhs;
	std::vector<Matrix> local_matrix;
	std::vector<std::vector<double>> m_bdycoef;
	std::vector<BdyBase> bdybaseid;
	std::vector<std::vector<int>> m_colids;
	std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
	m_bdycoef.resize(m_solids.Local_V3f.size());
	std::vector<double> _material_distribution(m_solids.Local_V3f.size(), 1);

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);    //父网格中的子网格乘3
	}
	
	/*-----liner------*/
	HexV3fMesh* newmesh1 = new HexV3fMesh();
	int PickedVecNum = 0;
	ApproximationCC AppCC1(mesh, newmesh1);
	std::vector<V3f> chosev;
	std::vector<VertexHandle> tempcheckv;
	std::vector<EdgeHandle> tempchecke;
	std::vector<FaceHandle> tempcheckf;
	std::vector<CellHandle> tempcheckc;
	int chosevsiz = 0;
	std::cout << pickedge.size() << std::endl;
	for (auto v : pickedge) std::cout << v << std::endl;
	for (int i = 0; i < pickedge.size(); i++)
	{
		EdgeHandle e = pickedge[i];
		chosev.push_back(AppCC1.solid_edge_point(e, mesh->edges_[e].from_h));
		chosev.push_back(AppCC1.solid_edge_point(e, mesh->edges_[e].to_h));
		if (std::find(tempcheckv.begin(), tempcheckv.end(), mesh->edges_[e].from_h) == tempcheckv.end())
		{
			chosev.push_back(AppCC1.solid_corner_point(mesh->edges_[e].from_h));
			tempcheckv.push_back(mesh->edges_[e].from_h);
		}
		if (std::find(tempcheckv.begin(), tempcheckv.end(), mesh->edges_[e].to_h) == tempcheckv.end())
		{
			chosev.push_back(AppCC1.solid_corner_point(mesh->edges_[e].to_h));
			tempcheckv.push_back(mesh->edges_[e].to_h);
		}
	}
	chosevsiz = chosev.size();
	std::cout << "chosevsiz " << chosevsiz << std::endl;
	for (int i = 0; i < pickface.size(); i++)
	{
		FaceHandle f = pickface[i];
		for (int v_num = 0; v_num < 4; v_num++)
		{
			chosev.push_back(AppCC1.solid_face_point(f, mesh->faces_[f].vs_f[v_num]));
			if (std::find(tempcheckv.begin(), tempcheckv.end(), mesh->faces_[f].vs_f[v_num]) == tempcheckv.end())
			{
				chosev.push_back(AppCC1.solid_corner_point(mesh->faces_[f].vs_f[v_num]));
				tempcheckv.push_back(mesh->faces_[f].vs_f[v_num]);
			}
		}
		for (int e_num = 0; e_num < 4; e_num++)
		{
			EdgeHandle e = mesh->faces_[f].edges_[e_num];
			if (std::find(tempchecke.begin(), tempchecke.end(), e) == tempchecke.end())
			{
				chosev.push_back(AppCC1.solid_edge_point(e, mesh->edges_[e].from_h));
				chosev.push_back(AppCC1.solid_edge_point(e, mesh->edges_[e].to_h));
				tempchecke.push_back(e);
			}
		}
	}

	for (int cen_i = 0; cen_i < chosev.size(); cen_i++)
	{
		V3f m_glo_cen1 = chosev[cen_i];
		int mindis_id = 0; double min_dis = 999999999; int mindis_cellid = 0;
		for (int i = 0; i < m_solids.Local_V3f.size(); i++)
		{
			for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
			{
				V3f temp_v = m_solids.Local_V3f[i][j];
				double dis = distance(temp_v, m_glo_cen1);
				if (dis < min_dis)
				{
					min_dis = dis;
					mindis_id = j;
					mindis_cellid = i;
				}
			}
		}
		if (cen_i < chosevsiz)
			m_bdycoef[mindis_cellid][3 * mindis_id + 1] = -1000;//Y方向
		else
		{

			m_bdycoef[mindis_cellid][3 * mindis_id] = 0;//X方向
			m_bdycoef[mindis_cellid][3 * mindis_id + 1] = 0;//Y方向
			m_bdycoef[mindis_cellid][3 * mindis_id + 2] = 0;//Z方向
		}
		bdybaseid.push_back(BdyBase(mindis_id, m_solids.Local_to_Global[mindis_cellid][mindis_id], mindis_cellid));
		if (cen_i == chosevsiz - 1)
			bdybaseid = bdybaseid;
	}
	/*-----开始拓扑优化------*/

	FunctionBSpline fb(m_solids, _material_distribution);
	fb.subnum = chosevsiz;

	// 第一次计算刚度矩阵并建立对应关系
	std::vector<double> m_bigrhs;
	std::vector<int> rows, cols;          //稀疏存储的元素行列信息
	std::vector<double> elems;             //稀疏存储的元素信息
	fb.Global_MatrixAndRhs_Sparse_liner(m_bdycoef, bdybaseid, local_matrix,
		m_bigrhs, local_bigrhs, rows, cols, elems, m_colids, m_bdyquick);

	///*-----------------求解线型方程组-----------------*/
	std::vector<double> m_coef;
	fb.SolveLinearEquation(-2, m_solids.Global_size * 3, rows, cols, elems, m_bigrhs, m_coef);

	for (int Local_V3f_id = 0; Local_V3f_id < m_solids.Local_V3f.size(); Local_V3f_id++)
		for (int v_id = 0; v_id < m_solids.Local_V3f[Local_V3f_id].size(); v_id++)
			for (int k = 0; k < 3; k++)
				m_bdycoef[Local_V3f_id][3 * v_id + k] = m_coef[3 * m_solids.Local_to_Global[Local_V3f_id][v_id] + k];

	std::vector<std::vector<double>> m_bdycoef_x, m_bdycoef_y, m_bdycoef_z;
	m_bdycoef_x.resize(m_solids.Local_V3f.size());
	m_bdycoef_y.resize(m_solids.Local_V3f.size());
	m_bdycoef_z.resize(m_solids.Local_V3f.size());
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		m_bdycoef_x[i].resize(8, 0);
		m_bdycoef_y[i].resize(8, 0);
		m_bdycoef_z[i].resize(8, 0);
	}

	for (int i = 0; i < m_solids.Local_V3f.size(); ++i)
	{
		for (int j = 0; j < 8; ++j)
		{
			m_bdycoef_x[i][j] = m_bdycoef[i][3 * j];
			m_bdycoef_y[i][j] = m_bdycoef[i][3 * j + 1];
			m_bdycoef_z[i][j] = m_bdycoef[i][3 * j + 2];
		}
	}
	x = m_bdycoef_x;
	y = m_bdycoef_y;
	z = m_bdycoef_z;
	
	/*std::cout << std::endl;
	for (int i = 0; i < mesh->cells_.size(); ++i) {
		CellHandle ch(i);
		std::cout << "_hex ch = " << ch << ", vhs: ";
		auto vhs = mesh->cells_[ch].vertex_;
		for (auto& vh : vhs) {
			std::cout << vh << " ";
		}
		std::cout << std::endl;
	}*/

	return;
}
