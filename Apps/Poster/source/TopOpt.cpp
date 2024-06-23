#include"TopOpt.h"
#include "IGA_CC_Helper.h"
#include "CatmullClarkBezierSolver.h"

using namespace CCBezier;

const double TargetRate = 0.5;       //优化的最终比例
const double MaxLoopTime = 20;
const double RemoveStep = 0.02;
double MinD(double a, double b) {
	if (a > b) return b;
	else return a;
}
double MaxD(double a, double b) {
	if (a > b) return a;
	else return b;
}

void readCe(std::vector<double>& Ce) {
	int n = Ce.size();

	std::string filename = "C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\Ce";
	std::ifstream in(filename);

	for (int i = 0; i < n; ++i) {
		in >> Ce[i];
	}
}
void write2density(std::vector<double>& density) {
	std::string _filename = "C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\density";
	std::ofstream off(_filename.c_str(), std::ios::out | std::ios::trunc);
	
	for (auto& d : density) {
		off << d << " ";
	}
	off.close();
}

void OC(std::vector<double>& Density, const double TargetRate, const std::vector<double>& dc) {
	double a = 0.0;
	double b = 100000.0;
	double MoveStep = 0.2;
	std::vector<double> TempDensity = Density;
	while (b - a > 0.0001) {
		double Lmid = 0.5 * (a + b);
		double SumDensity = 0.0;
		for (int i = 0; i < Density.size(); ++i) {
			TempDensity[i] = MaxD(0.001, MaxD(Density[i] - MoveStep, MinD(1.0, MinD(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
			SumDensity += TempDensity[i];
		}
		if (SumDensity - TargetRate * Density.size() > 0) {
			a = Lmid;
		}
		else {
			b = Lmid;
		}
	}
	Density = TempDensity;
}

void RemoveOneStep(std::vector<double>& Density,  const std::vector<double>& dc, const std::vector<int>& isPicked,
	std::vector<int>& isRemoved) {

	std::vector<int> ids;
	for (int i = 0; i < Density.size(); ++i) {
		ids.push_back(i);
	}
	// remove的排在后面
	// picked的排在后面
	sort(ids.begin(), ids.end(), [&](int a, int b) {
		return  (isPicked[a] < isPicked[b]
			|| (isPicked[a] == isPicked[b] && isRemoved[a] < isRemoved[b])
			|| (isPicked[a] == isPicked[b] && isRemoved[a] == isRemoved[b] && -dc[a] < -dc[b]));
		});

	int k = Density.size() * RemoveStep;
	if (k % 2) k--;
	for (int i = 0; i < k; ++i) {
		int ii = ids[i];
		if (isPicked[ii] || isRemoved[ii]) break;
		Density[ii] = 0.000001;
		isRemoved[ii] = 1;
	}
}
void OC(std::vector<double>& Density, const double TargetRate, const std::vector<double>& dc, const std::vector<int>& isPicked) {
	double a = 0.0;
	double b = 100000.0;
	double MoveStep = 0.2;

	std::vector<double> TempDensity = Density;
	while (b - a > 1e-10) {
		//std::cout << "[OC]: " << "(a,b) = " << a << "," << b << "\n";
		double Lmid = 0.5 * (a + b);
		double SumDensity = 0.0;
		for (int i = 0; i < Density.size(); ++i) {
			/*if (isPicked[i] == 1) TempDensity[i] = 1;
			else TempDensity[i] = MaxD(0.0001, MaxD(Density[i] - MoveStep, MinD(1.0, MinD(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));*/
			TempDensity[i] = MaxD(0.001, MaxD(Density[i] - MoveStep, MinD(1.0, MinD(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
			SumDensity += TempDensity[i];
		}
		if (SumDensity - TargetRate * Density.size() > 0) {
			a = Lmid;
		}
		else {
			b = Lmid;
		}
	}
	double SumDensityOutput = 0.0;
	for (int i = 0; i < Density.size(); ++i) {
		SumDensityOutput += TempDensity[i];
	}
	std::cout << "SumDensityOutput = " << SumDensityOutput / Density.size() << std::endl;

	Density = TempDensity;
}

struct PassPara {
public:
	PassPara(BSplineSolid _m_solids, std::vector<double> _material_distribution) :fb(_m_solids, _material_distribution) {}    //传到求解器的参数
	FunctionBSpline fb;       //解等几何分析的类
	std::vector<std::vector<double>> m_bdycoef;      //收到的外力
	std::vector<double> m_bigrhs;
	std::vector<std::vector<double>> u;              //计算出来的位移u
	std::vector<double> u_vec;                 
	std::vector<BdyBase> bdybaseid;
	std::vector<std::vector<Matrix>> local_matrixs;                       // 用于多分辨率
	std::vector<Matrix> local_matrixs_single_resolution;                  // 用于单分辨率
	std::vector<std::vector<double>> local_bigrhs_single_resolution;
	std::vector<std::vector<std::vector<double>>> local_bigrhs;
	std::vector<std::vector<int>> m_colids;
	std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
	std::vector<std::vector<int>> celllocal2global;
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	HexV3fMesh nextmesh;
};

struct LoopData {
public:
	LoopData(CatmullClarkBezier _model, IGA::CCBezierBdy _bdy) :solver(_model, _bdy) {}    //传到求解器的参数
	
	CatmullClarkBezierLinerElasticlSolver solver;
	std::vector<std::vector<double>> u;              //计算出来的位移u
	BSplineSolid m_solids;
	std::vector<Matrix> local_matrixs_single_resolution;                  // 用于单分辨率
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;

};

void FilterSensity(std::vector<double>& dc, const std::vector<double>& x, const std::vector<std::vector<CellHandle>>& aroundcell, const std::vector<std::vector<double>>& aroundcell_weight) {
	//对灵敏度进行过滤操作
	//std::cout << "Filter" << std::endl;
	auto Tempdc = dc;                         //定义一个暂时存储dc的变量
	int n = dc.size();                        //总的自由度数目
	for (int i = 0; i < n; ++i) {             //循环每一个小块
		double weight = 0;
		double temp_dc = 0;
		int aroundnum = aroundcell[i].size();        //当前块周围的cell数目
		for (int j = 0; j < aroundnum; ++j) {        //遍历这些周围的cell，取出weight，并与对应的块相乘，加到当前块专属的temp_dc中
			weight += aroundcell_weight[i][j];
			temp_dc += x[aroundcell[i][j]] * dc[aroundcell[i][j]] * aroundcell_weight[i][j];
		}
		temp_dc /= weight * x[i];
		Tempdc[i] = temp_dc;          //temp_dc赋值
	}
	dc = Tempdc;       //将计算的结果返回
}

void FilterSensity(std::vector<double>& x, const std::vector<std::vector<CellHandle>>& aroundcell, const std::vector<std::vector<double>>& aroundcell_weight) {
	//对灵敏度进行过滤操作
	//std::cout << "Filter" << std::endl;
	auto Tempdx = x;                         //定义一个暂时存储dc的变量
	int n = x.size();                        //总的自由度数目
	for (int i = 0; i < n; ++i) {             //循环每一个小块
		double weight = 0;
		double temp_dx = 0;
		int aroundnum = aroundcell[i].size();        //当前块周围的cell数目
		for (int j = 0; j < aroundnum; ++j) {        //遍历这些周围的cell，取出weight，并与对应的块相乘，加到当前块专属的temp_dc中
			weight += aroundcell_weight[i][j];
			temp_dx += x[aroundcell[i][j]] * aroundcell_weight[i][j];
		}
		temp_dx /= weight;
		Tempdx[i] = temp_dx;          //temp_dc赋值
	}
	x = Tempdx;       //将计算的结果返回
}

double Diff(const std::vector<double>& v1, const std::vector<double>& v2) {
	std::cout << "Diff" << std::endl;
	double res = 0.0;
	int n = v1.size();
	for (int i = 0; i < n; ++i) {
		res += fabs(v1[i] - v2[i]);
	}
	return res;
}

double myFunc(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data) {
	//在myFunc中应该根据x的值，修改刚度矩阵，并进行u的求解
	//std::cout << "myfunc" << std::endl;
	PassPara* PP = reinterpret_cast<PassPara*>(my_func_data);

	//修改刚度矩阵
	std::vector<double> m_bigrhs;
	std::vector<int> rows, cols;          //稀疏存储的元素行列信息
	std::vector<double> elems;            //稀疏存储的元素信息
	PP->fb.Global_MatrixAndRhs_Sparse_liner_mutires_quick(PP->m_bdycoef, PP->bdybaseid, PP->local_matrixs,
		m_bigrhs, PP->local_bigrhs, rows, cols, elems, PP->m_colids, PP->m_bdyquick);

	//求解PDE
	std::vector<double> m_coef;       //求解结果

	PP->fb.SolveLinearEquation(-2, PP->fb.m_solids.Global_size * 3, rows, cols, elems, m_bigrhs, m_coef);


	for (int Local_V3f_id = 0; Local_V3f_id < PP->fb.m_solids.Local_V3f.size(); Local_V3f_id++) {

		for (int v_id = 0; v_id < PP->fb.m_solids.Local_V3f[Local_V3f_id].size(); v_id++) {
			for (int k = 0; k < 3; k++)
			{
				PP->u[Local_V3f_id][3 * v_id + k] = m_coef[3 * PP->fb.m_solids.Local_to_Global[Local_V3f_id][v_id] + k];
			}
		}

	}
	//计算utku
	int n = x.size();
	std::vector<double> Ce(n);
	std::cout << std::endl;
	for (int i = 0; i < PP->fb.m_solids.Local_V3f.size(); i++)
	{
		for (int temsen = 0; temsen < 8; ++temsen)
		{
			int cell_cenid = PP->celllocal2global[i][temsen];
			Ce[cell_cenid] = PP->fb.uTKu(PP->u[i], PP->local_matrixs[i][temsen]);
			
		}

	}
	double f = 0.0;
	std::vector<double> dc(n);
	
	for (int i = 0; i < n; ++i) {
		dc[i] = (-3 * x[i] * x[i] * (1 - 0.0001)) * Ce[i];
		f += (0.0001 + x[i] * x[i] * x[i] * (1 - 0.0001)) * Ce[i];

	}
	
	

	if (grad.size() > 0) {
		FilterSensity(dc, x, PP->aroundcell, PP->aroundcell_weight);
		for (int i = 0; i < n; ++i) {
			grad[i] = dc[i];
		}
	}
	std::cout << std::endl;
	std::cout << "function:" << f << std::endl;
	return f;
}

double uTKu(std::vector<double> u, Matrix K0)
{
	int size = u.size();
	std::vector<double> cal_uTK(size, 0);
	for (int i = 0; i < size; i++)
	{
		for (int j = 0; j < size; j++)
		{
			cal_uTK[i] += u[j] * K0[i][j];
		}
	}
	double cal_uTKu = 0;
	for (int i = 0; i < size; i++)
	{
		cal_uTKu += cal_uTK[i] * u[i];
	}
	return cal_uTKu;
}

double single_resolution_topopt_func(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data) {
	//在myFunc中应该根据x的值，修改刚度矩阵，并进行u的求解

	PassPara* PP = reinterpret_cast<PassPara*>(my_func_data);

	std::vector<int> rows, cols;          //稀疏存储的元素行列信息
	std::vector<double> elems;             //稀疏存储的元素信息
	std::vector<double> m_coef;

	//修改刚度矩阵
	//PP->fb.Global_MatrixAndRhs_Sparse_liner_quick(PP->m_bdycoef, PP->bdybaseid, PP->local_matrixs_single_resolution,
	//	PP->m_bigrhs, PP->local_bigrhs_single_resolution, rows, cols, elems, PP->m_colids, PP->m_bdyquick, PP->u_vec);
	//std::cout << "组装完毕\n";

	/*PP->fb.Single_Resolution_TopOpt_Init_Step(PP->m_bdycoef, PP->bdybaseid, PP->local_matrixs_single_resolution,
		PP->m_bigrhs, PP->m_colids, PP->u_vec);*/
	PP->fb.Single_Resolution_TopOpt_Loop_Step(PP->local_matrixs_single_resolution, PP->m_bigrhs, PP->m_colids, PP->bdybaseid, PP->u_vec);

	std::cout << "Single_Resolution_TopOpt_Loop_Step end" << std::endl;

	double avg_u = 0.0;
	for (int Local_V3f_id = 0; Local_V3f_id < PP->fb.m_solids.Local_V3f.size(); Local_V3f_id++) {
		for (int v_id = 0; v_id < PP->fb.m_solids.Local_V3f[Local_V3f_id].size(); v_id++) {
			for (int k = 0; k < 3; k++)
			{
				PP->u[Local_V3f_id][3 * v_id + k] = PP->u_vec[3 * PP->fb.m_solids.Local_to_Global[Local_V3f_id][v_id] + k];
				avg_u += fabs(PP->u[Local_V3f_id][3 * v_id + k]);
			}
		}
	}
	avg_u /= PP->fb.m_solids.Local_V3f.size() * 64 * 3;
	std::cout << "avg_u = " << avg_u << std::endl;
	//计算utku
	int n = x.size();
	std::vector<double> Ce(n);
	for (int i = 0; i < PP->fb.m_solids.Local_V3f.size(); i++)
	{
		Ce[i] = PP->fb.uTKu(PP->u[i], PP->local_matrixs_single_resolution[i]);
	}

	double f = 0.0;
	std::vector<double> dc(n);

	for (int i = 0; i < n; ++i) {
		
		dc[i] = (-3 * x[i] * x[i] * (1 - 0.0001)) * Ce[i];
		f += (0.0001 + x[i] * x[i] * x[i] * (1 - 0.0001)) * Ce[i];

	}
	
	if (grad.size() > 0) {
		FilterSensity(dc, x, PP->aroundcell, PP->aroundcell_weight);
		for (int i = 0; i < n; ++i) {
			grad[i] = dc[i];
		}
	}
	std::cout << std::endl;
	std::cout << "function:" << f << std::endl;
	return f;
}


double single_resolution_topopt_func_CatmullClark(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data) {
	//在myFunc中应该根据x的值，修改刚度矩阵，并进行u的求解
	LoopData* PP = reinterpret_cast<LoopData*>(my_func_data);

	PP->solver.SolvingWithDensity(x, PP->local_matrixs_single_resolution);

	PP->u.resize(PP->m_solids.Local_V3f.size(), std::vector<double>(3 * 64, 0));

	for (int patch_id = 0; patch_id < PP->m_solids.Local_V3f.size(); patch_id++)
		for(int ww = 0; ww < 4; ++ww)
			for(int vv = 0; vv < 4; ++vv)
				for (int uu = 0; uu < 4; ++uu)
				{
					int cnt = uu + vv * 4 + ww * 16;

					int cur_local_id = PP->m_solids.Control_Point[patch_id][ww][vv][uu];
					int cur_global_id = PP->m_solids.Local_to_Global[patch_id][cur_local_id];

					for (int k = 0; k < 3; k++) {
						PP->u[patch_id][cnt + k * 64] 
							= PP->solver.m_u[cur_global_id + k * PP->m_solids.Global_size];
					}
						
				}
			

	//计算utku
	int n = x.size();
	std::vector<double> Ce(n);
	for (int i = 0; i < PP->m_solids.Local_V3f.size(); i++)
	{
		Ce[i] = uTKu(PP->u[i], PP->local_matrixs_single_resolution[i]);
	}

	double f = 0.0;
	std::vector<double> dc(n);

	for (int i = 0; i < n; ++i) {

		dc[i] = (-3 * x[i] * x[i]) * Ce[i];
		f += (x[i] * x[i] * x[i]) * Ce[i];
		std::cout << dc[i] << " ";
	}
	std::cout << "\n";
	std::cout << "grad.size() = " << grad.size() << "/n";
	if (grad.size() > 0) {
		FilterSensity(dc, x, PP->aroundcell, PP->aroundcell_weight);
		for (int i = 0; i < n; ++i) {
			grad[i] = dc[i];
		}
	}
	std::cout << std::endl;
	std::cout << "function:" << f << std::endl;
	return f;
}
double single_resolution_topopt_func_CatmullClark_GPUSolver(const std::vector<double>& x, std::vector<double>& grad, void* my_func_data) {
	//在myFunc中应该根据x的值，修改刚度矩阵，并进行u的求解
	LoopData* PP = reinterpret_cast<LoopData*>(my_func_data);

	system("C:\\Users\\ASUS\\Desktop\\Dev\\Dependence\\GPMeshLib\\cmake-build-release-visual-studio\\tests\\core\\core.exe");

	
	int n = x.size();
	std::vector<double> Ce(n);
	readCe(Ce);

	double f = 0.0;
	std::vector<double> dc(n);

	for (int i = 0; i < n; ++i) {

		dc[i] = (-3 * x[i] * x[i]) * Ce[i];
		f += (x[i] * x[i] * x[i]) * Ce[i];
		std::cout << dc[i] << " ";
	}
	std::cout << "\n";
	std::cout << "grad.size() = " << grad.size() << "/n";
	if (grad.size() > 0) {
		FilterSensity(dc, x, PP->aroundcell, PP->aroundcell_weight);
		for (int i = 0; i < n; ++i) {
			grad[i] = dc[i];
		}
	}
	std::cout << std::endl;
	std::cout << "function:" << f << std::endl;
	return f;
}

void SaveCurMesh(const std::vector<double>& Density, HexV3fMesh* NextMesh, std::string& FileName, std::atomic_bool& is_need_update_topopt_model) {
	//std::cout << "Save Mesh" << std::endl;
	//if (is_need_update_topopt_model == false) {
		HexV3fMesh* PrintMesh = new HexV3fMesh();
		*PrintMesh = *NextMesh;
		for (int i = 0; i < Density.size(); ++i) {
			//std::cout << Density[i] << std::endl;
			if (Density[i] < 0.5) {
				PrintMesh->delete_cell(CellHandle(i));
			}
		}
		FileHandler filenew = FileHandler();
		string str = "topopt_temp.mesh";
		filenew.write2VTKFile(FileName, PrintMesh);
		if (is_need_update_topopt_model == false) {
			filenew.writeFile(str, PrintMesh);
			is_need_update_topopt_model = true;
		}

		delete PrintMesh;
		
	//}
	
}

void SaveCurMesh(const std::vector<double>& Density, HexV3fMesh* NextMesh, std::string& FileName) {
	//std::cout << "Save Mesh" << std::endl;
	
		HexV3fMesh* PrintMesh = new HexV3fMesh();
		*PrintMesh = *NextMesh;
		for (int i = 0; i < Density.size(); ++i) {
			//std::cout << Density[i] << std::endl;
			if (Density[i] < 0.5) {
				PrintMesh->delete_cell(CellHandle(i));
			}
		}
		FileHandler filenew = FileHandler();
		string str = "topopt_temp.mesh";
		filenew.write2VTKFile(FileName, PrintMesh);
		filenew.writeFile(str, PrintMesh);

		delete PrintMesh;
}




//std::string ModelName = "SimpleTestBeam";
std::string ModelName = "TestModel4";
//std::string ModelName = "SimpleTestBeam_1";
//std::string ModelFile = "C:\\Users\\74680\\Desktop\\TopOpt\\TopOpt\\Data\\TestModel4.mesh";
//std::string PickedDomain = "C:\\Users\\74680\\Desktop\\TopOpt\\TopOpt\\Pick\\" + ModelName + ".ppa";
//std::string SharpEdgeFile = "C:\\Users\\74680\\Desktop\\TopOpt\\TopOpt\\Pick\\" + ModelName + "_sharp.ppe";
std::string OutPutFile = "./Result/TopOptResult/" +  ModelName + "/";
//std::string OutPutFile = "C:\\Users\\74680\\Desktop\\TopOpt\\TopOpt\\" + ModelName + "\\";

void TopOpt(ModelDraw* widget, std::string _filename, std::vector<std::vector<int>>& _pickface, std::vector<vector<int>>& _pickedge, 
	std::vector<std::vector<int>>& protectcell)
{

	int filenum = 0;
	/*-----------------读取文件-----------------*/
	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(_filename);

	widget->loadHexMesh(_filename);
	

	file.read2HexMeshObj(mesh, true, true, true);       //读取初始文件
	std::cout << "read2HexMeshObj OK" << std::endl;
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
			if (flag) {pickface.push_back((FaceHandle)j); break;}
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
	for (int i = 0; i < protectcell.size(); ++i) {
		std::sort(protectcell[i].begin(), protectcell[i].end());
		for (int j = 0; j < mesh->cells_.size(); ++j) {
			auto cv = mesh->cells_[j].vertex_;
			std::sort(cv.begin(), cv.end());
			bool flag = true;
			for (int k = 0; k < cv.size(); ++k) {
				if (cv[k].idx() != protectcell[i][k]) {
					flag = false;
					break;
				}
			}
			if (flag) { pickcell.push_back((CellHandle)j); break; }
		}
	}
	/*FileHandler file1(PickedDomain);
	file1.readpickedall(pickvertex, pickedge, pickface, pickcell);
	std::cout << "READ BC OK" << std::endl;*/
	/*-----------------brige-----------------*/
	/*int f1n = 0, f2n = 0;
	for (int i = 0; i < mesh->faces_.size(); ++i)
	{
		if (mesh->is_bdy(FaceHandle(i)))
		{
			V3f fv = mesh->cal_face_center(FaceHandle(i));
			if (fv.z == 12)
			{
				f1n++;
				pickface.push_back(FaceHandle(i));
			}
		}
	}
	std::cout << "f1n: " << f1n << std::endl;
	for (int i = 0; i < mesh->faces_.size(); ++i)
	{
		if (mesh->is_bdy(FaceHandle(i)))
		{
			V3f fv = mesh->cal_face_center(FaceHandle(i));
			if (fv.z < 0)
			{
				f2n++;
				pickface.push_back(FaceHandle(i));
			}
		}
	}
	std::cout << "f2n: " << f2n << std::endl;
	std::cout << "pickface size :" << pickface.size() << std::endl;*/
	/*-----------------获取边界条件csv-----------------*/
	mesh->v_cs.resize(mesh->vertices_.size(), false);                    //折痕
	/*-----------------获取边界条件sharp-----------------*/
	std::vector<V3f> edgemid;
	std::vector<EdgeHandle> edgesharp;
	/*FileHandler filesharp(SharpEdgeFile);
	filesharp.readpickededge(edgesharp); */                     //读取边界边（保留）
	std::cout << "READ sharp OK" << std::endl;
	std::cout << edgesharp.size() << std::endl;
	std::cout << mesh->v_cs.size() << std::endl;
	//for (int i = 0; i < edgesharp.size(); ++i)
	//{
	//	mesh->v_cs[mesh->edges_[edgesharp[i]].from_h] = true;                 //边界边上的点保留
	//	mesh->v_cs[mesh->edges_[edgesharp[i]].to_h] = true;
	//}
	//std::cout << "Sharp edge size :" << edgesharp.size() << std::endl;


	HexV3fMesh* newmesh = new HexV3fMesh();
	mesh->cal_cell_cen();
	mesh->cal_cen();
	std::vector<double> material_distribution;
	material_distribution.resize(mesh->cells_.size(), 1);	           //所有材料属性设置为1
	/*-----------------下一步的细分网格-----------------*/
	std::cout << "Start nest subdivid mesh" << std::endl;
	std::vector<double> material_distribution_mutires;
	int tosubtime = 1;
	/*--------------一次-------------------*/
	std::cout << "Start first subdivid mesh" << std::endl;
	HexV3fMesh* nextsubmesh = new HexV3fMesh();
	CCSubdivision subcsv(mesh, nextsubmesh, mesh->v_cs, material_distribution);
	//*nextsubmesh = subcsv.CCSubWithCrease();                    //带有折痕的细分一次
	*nextsubmesh = subcsv.CCSubdivising();

	subcsv.check_cse();
	nextsubmesh->examine_bdy();
	nextsubmesh->cal_cell_cen();
	material_distribution_mutires = subcsv.material_distribution;

	std::cout << "Start building BSplineSolid data structure" << std::endl;
	ApproximationCC AppCCSIMP(mesh, newmesh, edgesharp);          //构造bezier体逼近
	std::cout << "Start building BSplineSolid data structure" << std::endl;
	AppCCSIMP.construct_bss_new();
	BSplineSolid m_solids = AppCCSIMP.bss;

	std::cout << "Building BSplineSolid data structure completely" << std::endl;

	Matrix m_bigmatrix;
	std::vector<std::vector<std::vector<double>>> local_bigrhs;
	std::vector<std::vector<Matrix>> local_matrixs;
	std::vector<std::vector<double>> m_bdycoef;
	std::vector<std::vector<double>> m_bdycoef_1;
	std::vector<BdyBase> bdybaseid;
	std::vector<BdyBase> bdybaseid_1;
	std::vector<std::vector<std::pair<int, int>>> m_bdyquick;
	std::vector<std::vector<std::vector<std::pair<std::pair<int, int>, std::pair<int, int>>>>> matrixs_quick;
	std::vector<std::vector<std::vector<int>>> matrixs_quick_check;
	std::vector<std::vector<V3f>> centers;
	std::vector<std::vector<int>> m_colids;

	m_bdycoef.resize(m_solids.Local_V3f.size());                 //父网格的受力

	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);    //父网格中的子网格乘3
	}


	m_bdycoef_1 = m_bdycoef;
	std::vector<std::vector<int>> celllocal2global;              //局部到全局
	std::vector<std::vector<CellHandle>> aroundcell;             //子网格一个cell周围的cell
	std::vector<std::vector<double>> aroundcell_weight;          //子网格周围网格的权重
	std::cout << "Init OK" << std::endl;
	/*-----liner------*/
	HexV3fMesh* newmesh1 = new HexV3fMesh();
	std::vector<int> isPickedVec(nextsubmesh->cells_.size(), 0);     //是否被保护
	int PickedVecNum = 0;
	ApproximationCC AppCC1(mesh, newmesh1);
	std::vector<V3f> chosev;
	std::vector<VertexHandle> tempcheckv;
	std::vector<EdgeHandle> tempchecke;
	std::vector<FaceHandle> tempcheckf;
	std::vector<CellHandle> tempcheckc;
	std::set<CellHandle> pick_coarse_cells;
	int chosevsiz = 0;
	std::cout << pickedge.size() << std::endl;
	for (auto v : pickedge) std::cout << v << std::endl;

	for (int i = 0; i < pickedge.size(); i++)
	{
		EdgeHandle e = pickedge[i];
		auto ec = mesh->neighbor_e_c(e);


		for (auto& cc : ec) pick_coarse_cells.insert(cc);


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
		auto fc = mesh->neighbor_f[f];
		for (auto& cc : fc) pick_coarse_cells.insert(cc);
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
			m_bdycoef_1[mindis_cellid][3 * mindis_id + 1] = -1000;//Y方向
		else
		{
			m_bdycoef_1[mindis_cellid][3 * mindis_id] = 0;//X方向
			m_bdycoef_1[mindis_cellid][3 * mindis_id + 1] = 0;//Y方向
			m_bdycoef_1[mindis_cellid][3 * mindis_id + 2] = 0;//Z方向
		}

		bdybaseid.push_back(BdyBase(mindis_id, m_solids.Local_to_Global[mindis_cellid][mindis_id], mindis_cellid));


		if (cen_i == chosevsiz - 1)
			bdybaseid_1 = bdybaseid;
	}
	/*-----开始拓扑优化------*/
	

	PassPara PP(m_solids, material_distribution);
	PassPara* PP_point = &PP;
	PP.fb.subnum = chosevsiz;
	//PP.fb.subnum = fvn;
	PP.m_bdycoef = m_bdycoef_1;
	PP.bdybaseid = bdybaseid;
	PP.local_matrixs = local_matrixs;
	PP.local_bigrhs = local_bigrhs;
	PP.m_colids = m_colids;
	PP.m_bdyquick = m_bdyquick;
	PP.u = m_bdycoef;
	std::cout << "Start iteration" << std::endl;

	// 第一次计算刚度矩阵并建立对应关系
	std::vector<double> m_bigrhs;
	std::vector<int> rows, cols;          //稀疏存储的元素行列信息
	std::vector<double> elems;            //稀疏存储的元素信息
	PP.fb.Global_MatrixAndRhs_Sparse_liner_mutires(PP.m_bdycoef, PP.bdybaseid, PP.local_matrixs,
		m_bigrhs, PP.local_bigrhs, rows, cols, elems, PP.m_colids, PP.m_bdyquick, centers);
	celllocal2global.resize(m_solids.Local_V3f.size(), std::vector<int>(pow(8, tosubtime), 0));

	PP.fb.material_distribution_mutires.resize(m_solids.Local_V3f.size(), std::vector<double>(pow(8, tosubtime), 1));
	std::cout << "第一次刚度矩阵组装完毕" << std::endl;
	double InitDensity = (nextsubmesh->cells_.size() * TargetRate - PickedVecNum) / nextsubmesh->cells_.size();
	std::unordered_map<int, int> cnt;
	int errmap = 0;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		for (int temsen = 0; temsen < pow(8, tosubtime); ++temsen)
		{
			//找到相应块的下标
			V3f tempcenter = centers[i][temsen];
			int cell_cenid; double mindis = 99999;
			for (int cell_it = 0; cell_it < nextsubmesh->cells_.size(); ++cell_it)
			{
				V3f temp_cell_center = nextsubmesh->cell_center[cell_it];
				double dis = distance(tempcenter, temp_cell_center);
				if (dis < mindis)
				{
					mindis = dis;
					cell_cenid = cell_it;
				}
			}
			//std::cout << mindis << std::endl;
			celllocal2global[i][temsen] = cell_cenid;
			cnt[cell_cenid]++;
			if (cnt[cell_cenid] > 1) errmap++;
			if (pick_coarse_cells.count((CellHandle)i) || (std::find(pickcell.begin(), pickcell.end(), i) != pickcell.end())) isPickedVec[cell_cenid] = 1;
			if (isPickedVec[cell_cenid] == 1) PP.fb.material_distribution_mutires[i][temsen] = 1;
			//else PP.fb.material_distribution_mutires[i][temsen] = InitDensity;     //更新子网格的密度值
			else PP.fb.material_distribution_mutires[i][temsen] = 1;
		}
	}
	std::cout << "errmap:" << errmap << std::endl;
	PP.celllocal2global = celllocal2global;
	std::cout << "全局到局部关系建立" << std::endl;
	for (int i = 0; i < nextsubmesh->cells_.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = 0;
		std::set<CellHandle> c_c = nextsubmesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			rmin = rmin + nextsubmesh->cell_r[i] + nextsubmesh->cell_r[*c_c_it];
		}
		rmin = rmin / c_c.size();
		rmin = rmin * 1.5;
		for (int k = 0; k < nextsubmesh->cells_.size(); k++)
		{
			double rdis = distance(nextsubmesh->cell_center[i], nextsubmesh->cell_center[k]);
			if (rdis < rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
	}
	PP.aroundcell = aroundcell;
	PP.aroundcell_weight = aroundcell_weight;
	std::cout << "邻接权重关系建立" << std::endl;
	double DensityChange = 100.0;
	int LoopTime = 1;
	std::vector<double> Density(nextsubmesh->cells_.size(), InitDensity);
	std::vector<double> DensityPre(nextsubmesh->cells_.size(), InitDensity);
	//设置选定的块的密度为1
	for (int i = 0; i < nextsubmesh->cells_.size(); ++i) {
		//if (isPickedVec[i] == 1) {
		Density[i] = DensityPre[i] = 1.0;
		//}
	}
	std::vector<double> dc(nextsubmesh->cells_.size(), 0.0);
	int NumOfVar = nextsubmesh->cells_.size();
	double f = 0.0;
	std::cout << "开始迭代" << std::endl;
	
	while (DensityChange > 0.1 && LoopTime < 1000) {
		//计算梯度信息
		f = myFunc(Density, dc, PP_point);
		//OC法计算密度
		OC(Density, TargetRate, dc, isPickedVec);
		//更新fb中的密度信息
		std::cout << "更新并输出密度" << std::endl;
		for (int i = 0; i < m_solids.Local_V3f.size(); i++)
		{
			for (int temsen = 0; temsen < pow(8, tosubtime); ++temsen)
			{
				int cell_cenid = celllocal2global[i][temsen];
				PP.fb.material_distribution_mutires[i][temsen] = Density[cell_cenid];
				if (cell_cenid % 500 == 0) {
					std::cout << Density[cell_cenid] << " ";
				}
			}
		}
		std::cout << std::endl;

		//输出信息
		DensityChange = Diff(Density, DensityPre);
		DensityPre = Density;
		std::cout << "第" << LoopTime << "次循环后DensityChange为：" << DensityChange << std::endl;
		std::cout << "第" << LoopTime << "次循环后柔度为：" << f << std::endl;
		std::string FileName(OutPutFile + std::to_string(LoopTime) + ".vtk");
		//SaveCurMesh(Density, nextsubmesh, FileName, widget);
		SaveCurMesh(Density, nextsubmesh, FileName);
		string str = "temp.mesh";
		//widget->loadHexMesh(str);
		widget->loadHexMesh(str, false);
		string qstr = string("第") + to_string(LoopTime) + string("/30次迭代");
		//info->printString(QString(QString::fromLocal8Bit(qstr.data())));
		std::cout << "第" << LoopTime << "次模型保存完毕" << std::endl;
		if (LoopTime == 100) return;
		LoopTime++;
		filenum++;
	}
}



void TopOpt_CC_Bezier_V1(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model) {

	int filenum = 0;

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	// 读文件
	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(_filename);
	widget->loadHexMesh(_filename);


	// 逼近bezier体
	// Todo: sharpedge

	std::cout << "read2HexMeshObj" << std::endl;
	file.read2HexMeshObj(mesh, true, true, true);       //读取初始文件
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	BSplineSolid m_solids = AppCCSIMP.bss;

	//// 得到CC边界条件
	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert(_filename, hexbdy, ccbdy);

	// 初始化仿真模块
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
		m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);  
	}
	std::cout << "ccbdy.size() = " << ccbdy.size() << std::endl;

	std::vector<int> isProtectCell(m_solids.Local_V3f.size(), 0);
	for (auto& bdy : ccbdy)
	{
		if (bdy.cell_id == 0) {
			std::cout << "bdy.cell_id = 0:" << bdy.global_base_id << " " << bdy.local_base_id << " " << bdy.bdy_value << std::endl;
		}
		if (bdy.bdy_type == 2) {
			for (int i = 0; i < 3; ++i)	m_bdycoef[bdy.cell_id][3 * bdy.local_base_id + i] = bdy.bdy_value[i];
			//m_bdycoef[bdy.cell_id][3 * bdy.local_base_id + 1] = -1;
		}
		isProtectCell[bdy.cell_id] = 1;
		bdybaseid.push_back(BdyBase(bdy.local_base_id, bdy.global_base_id, bdy.cell_id, static_cast<int>(bdy.bdy_type)));
	}

	std::vector<double> m_bigrhs;
	std::vector<int> rows, cols;          //稀疏存储的元素行列信息
	std::vector<double> elems;             //稀疏存储的元素信息
	std::vector<double> m_coef;

	PassPara PP(m_solids, _material_distribution);

	PassPara* PP_pointer = &PP;
	PP.bdybaseid = bdybaseid;
	PP.local_matrixs_single_resolution = local_matrix;
	PP.m_colids = m_colids;
	PP.m_bdyquick = m_bdyquick;
	PP.u_vec = m_coef;
	PP.m_bigrhs = m_bigrhs;
	PP.u = m_bdycoef;
	PP.m_bdycoef = m_bdycoef;

	/*PP.fb.Global_MatrixAndRhs_Sparse_liner(PP.m_bdycoef, PP.bdybaseid, PP.local_matrixs_single_resolution,
		PP.m_bigrhs, PP.local_bigrhs_single_resolution, rows,
		cols, elems, PP.m_colids, PP.m_bdyquick);*/
	PP.fb.Single_Resolution_TopOpt_Init_Step(PP.m_bdycoef, PP.bdybaseid, PP.local_matrixs_single_resolution,
		PP.m_bigrhs, PP.m_colids, PP.u_vec);

	std::unordered_map<int, int> cnt;

	// 施加保护体
	
	int errmap = 0;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		auto it = find(IGA::protected_cell.begin(), IGA::protected_cell.end(), i);
		if (it != IGA::protected_cell.end()) isProtectCell[i] = 1;
	}


	// 施加filter需要的信息
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = 999999999;
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}
		rmin = rmin * 1.5;

		std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis < rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}
	
	PP.aroundcell = aroundcell;
	PP.aroundcell_weight = aroundcell_weight;


	// 初始化密度 和 梯度
	double DensityChange = 100.0;
	int LoopTime = 1;
	std::vector<double> Density(mesh->cells_.size(), 1.0);
	std::vector<double> DensityPre(mesh->cells_.size(), 1.0);
	std::vector<double> dc(mesh->cells_.size(), 0.0);
	int NumOfVar = mesh->cells_.size();

	double f = 0.0;
	std::cout << "开始迭代" << std::endl;
	std::vector<int> isRemoved(mesh->cells_.size(), 0);
	while (DensityChange > 0.1 && LoopTime < MaxLoopTime) {
		//计算梯度信息
		f = single_resolution_topopt_func(Density, dc, PP_pointer);
		std::cout << "single_resolution_topopt_func end" << std::endl;
		//迭代法
		RemoveOneStep(Density, dc, isProtectCell, isRemoved);
		//OC法计算密度
		//OC(Density, TargetRate, dc, isProtectCell);
		//更新fb中的密度信息
		std::cout << "更新并输出密度" << std::endl;
		for (int i = 0; i < m_solids.Local_V3f.size(); i++)
		{
			PP.fb.material_distribution[i] = Density[i];
			//std::cout << Density[i] << " " << std::endl;
		}
		std::cout << std::endl;

		//输出信息
		DensityChange = Diff(Density, DensityPre);
		DensityPre = Density;
		std::cout << "第" << LoopTime << "次循环后DensityChange为：" << DensityChange << std::endl;
		std::cout << "第" << LoopTime << "次循环后柔度为：" << f << std::endl;
		std::string FileName(OutPutFile + std::to_string(LoopTime) + ".vtk");

		SaveCurMesh(Density, mesh, FileName, is_need_update_topopt_model);



		string qstr = string("第") + to_string(LoopTime) + string("/") + to_string(MaxLoopTime) + "次迭代";
		
		std::cout << "第" << LoopTime << "次模型保存完毕" << std::endl;
		LoopTime++;
		filenum++;
	}
}


void TopOpt_CC_Bezier_V2(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model) {

	int filenum = 0;

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	// 读文件
	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(_filename);
	widget->loadHexMesh(_filename);


	// 逼近bezier体
	// Todo: sharpedge

	std::cout << "read2HexMeshObj" << std::endl;
	file.read2HexMeshObj(mesh, true, true, true);       //读取初始文件
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	BSplineSolid m_solids = AppCCSIMP.bss;

	//// 得到CC边界条件
	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert(_filename, hexbdy, ccbdy);
	std::vector<int> isProtectCell(m_solids.Local_V3f.size(), 0);

	for (auto& bdy : ccbdy)
	{
		
		isProtectCell[bdy.cell_id] = 1;
	}
	// 初始化仿真模块
	
	std::vector<std::vector<double>> m_bdycoef;
	m_bdycoef.resize(m_solids.Local_V3f.size());
	std::vector<double> _material_distribution(m_solids.Local_V3f.size(), 1);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);
	}
	

	std::unordered_map<int, int> cnt;

	// 施加保护体

	int errmap = 0;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		auto it = find(IGA::protected_cell.begin(), IGA::protected_cell.end(), i);
		if (it != IGA::protected_cell.end()) isProtectCell[i] = 1;
	}

	
	// 施加filter需要的信息
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = 999999999;
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}

		rmin = rmin * 1.5;
		std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}

	CatmullClarkBezier model;
	model.control_vertices.resize(m_solids.Global_size);
	for (int patch_id = 0; patch_id < m_solids.Control_Point.size(); ++patch_id) {
		std::vector<int> patch_global_ids;
		for (int w_id = 0; w_id < 4; ++w_id) {
			for (int v_id = 0; v_id < 4; ++v_id) {
				for (int u_id = 0; u_id < 4; ++u_id) {
					int cur_local_id = m_solids.Control_Point[patch_id][w_id][v_id][u_id];
					int cur_global_id = m_solids.Local_to_Global[patch_id][cur_local_id];
					patch_global_ids.push_back(cur_global_id);
					assert(cur_global_id >= 0);
					double cx = m_solids.Local_V3f[patch_id][cur_local_id].x;
					double cy = m_solids.Local_V3f[patch_id][cur_local_id].y;
					double cz = m_solids.Local_V3f[patch_id][cur_local_id].z;
					model.control_vertices[cur_global_id] = Eigen::Vector3d(cx, cy, cz);
				}
			}
		}
		model.global_id.push_back(patch_global_ids);
		assert(model.global_id.size() == 64);
	}
	std::cout << "model generate" << std::endl;
	LoopData PP(model, ccbdy);
	PP.aroundcell = aroundcell;
	PP.aroundcell_weight = aroundcell_weight;
	PP.m_solids = m_solids;
	LoopData* PP_pointer = &PP;

	// 初始化密度 和 梯度
	double DensityChange = 100.0;
	int LoopTime = 1;
	std::vector<double> Density(mesh->cells_.size(), 1.0);
	std::vector<double> DensityPre(mesh->cells_.size(), 1.0);
	std::vector<double> dc(mesh->cells_.size(), 0.0);
	int NumOfVar = mesh->cells_.size();

	double f = 0.0;
	std::cout << "开始迭代" << std::endl;
	std::vector<int> isRemoved(mesh->cells_.size(), 0);
	while (DensityChange > 0.1 && LoopTime < MaxLoopTime) {
		//计算梯度信息
		f = single_resolution_topopt_func_CatmullClark(Density, dc, PP_pointer);
		std::cout << "single_resolution_topopt_func end" << std::endl;
		//迭代法
		//RemoveOneStep(Density, dc, isProtectCell, isRemoved);
		OC(Density, TargetRate, dc, isProtectCell);
		//OC(Density, TargetRate, dc);
		DensityChange = Diff(Density, DensityPre);
		DensityPre = Density;
		std::cout << "第" << LoopTime << "次循环后DensityChange为：" << DensityChange << std::endl;
		std::cout << "第" << LoopTime << "次循环后柔度为：" << f << std::endl;
		std::string FileName(OutPutFile + std::to_string(LoopTime) + ".vtk");
		SaveCurMesh(Density, mesh, FileName, is_need_update_topopt_model);
		string qstr = string("第") + to_string(LoopTime) + string("/") + to_string(MaxLoopTime) + "次迭代";
		std::cout << "第" << LoopTime << "次模型保存完毕" << std::endl;
		LoopTime++;
		filenum++;
	}
	widget->hex_density = Density;
	is_need_update_topopt_model = false;
	
	std::cout << "拓扑优化结束" << std::endl;

}

void TopOpt_CC_Bezier_GPUSolver(ModelDraw* widget,
	std::string _filename,
	IGA::HexBdy& hexbdy,
	std::atomic_bool& is_need_update_topopt_model) {

	int filenum = 0;

	if (hexbdy.empty()) {
		std::cerr << "No Boundary Condition!" << std::endl;
		return;
	}

	// 读文件
	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(_filename);
	widget->loadHexMesh(_filename);


	// 逼近bezier体
	// Todo: sharpedge

	std::cout << "read2HexMeshObj" << std::endl;
	file.read2HexMeshObj(mesh, true, true, true);       //读取初始文件
	ApproximationCC AppCCSIMP(mesh, nullptr);
	AppCCSIMP.construct_bss_new();
	BSplineSolid m_solids = AppCCSIMP.bss;
	IGA_CC_Helper::write2Bezier("",m_solids);

	//// 得到CC边界条件
	IGA::CCBezierBdy ccbdy;
	IGA::BdyConvert(_filename, hexbdy, ccbdy);
    IGA_CC_Helper::write2bdy("", ccbdy, m_solids);

	std::vector<int> isProtectCell(m_solids.Local_V3f.size(), 0);

	for (auto& bdy : ccbdy)
	{
		isProtectCell[bdy.cell_id] = 1;
	}
	// 初始化仿真模块

	std::vector<std::vector<double>> m_bdycoef;
	m_bdycoef.resize(m_solids.Local_V3f.size());
	std::vector<double> _material_distribution(m_solids.Local_V3f.size(), 1);
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		m_bdycoef[i].resize(m_solids.Local_V3f[i].size() * 3, 0);
	}


	std::unordered_map<int, int> cnt;

	// 施加保护体

	int errmap = 0;
	for (int i = 0; i < m_solids.Local_V3f.size(); i++)
	{
		auto it = find(IGA::protected_cell.begin(), IGA::protected_cell.end(), i);
		if (it != IGA::protected_cell.end()) isProtectCell[i] = 1;
	}


	// 施加filter需要的信息
	std::vector<std::vector<CellHandle>> aroundcell;
	std::vector<std::vector<double>> aroundcell_weight;
	for (int i = 0; i < m_solids.Local_V3f.size(); ++i) {
		std::vector<CellHandle> temp_aroundcell;
		std::vector<double> temp_aroundcell_weight;
		double rmin = 999999999;
		std::set<CellHandle> c_c = mesh->neighbor_c_c(CellHandle(i));
		for (std::set<CellHandle>::iterator c_c_it = c_c.begin(); c_c_it != c_c.end(); c_c_it++)
		{
			// rmin = rmin + mesh->cell_r[i] + mesh->cell_r[*c_c_it];
			//rmin = rmin + distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]);
			rmin = fmin(rmin, distance(mesh->cell_center[i], mesh->cell_center[*c_c_it]));
		}

		rmin = rmin * 1.5;
		std::cout << "rmin = " << rmin << std::endl;
		for (int k = 0; k < mesh->cells_.size(); k++)
		{
			//double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			double rdis = distance(mesh->cell_center[i], mesh->cell_center[k]);
			if (rdis <= rmin)
			{
				temp_aroundcell.push_back(CellHandle(k));
				temp_aroundcell_weight.push_back(rmin - rdis);
			}
		}
		aroundcell.push_back(temp_aroundcell);
		aroundcell_weight.push_back(temp_aroundcell_weight);
		std::cout << "aroundcell.size() = " << aroundcell[i].size() << std::endl;
	}

	CatmullClarkBezier model;
	model.control_vertices.resize(m_solids.Global_size);
	for (int patch_id = 0; patch_id < m_solids.Control_Point.size(); ++patch_id) {
		std::vector<int> patch_global_ids;
		for (int w_id = 0; w_id < 4; ++w_id) {
			for (int v_id = 0; v_id < 4; ++v_id) {
				for (int u_id = 0; u_id < 4; ++u_id) {
					int cur_local_id = m_solids.Control_Point[patch_id][w_id][v_id][u_id];
					int cur_global_id = m_solids.Local_to_Global[patch_id][cur_local_id];
					patch_global_ids.push_back(cur_global_id);
					assert(cur_global_id >= 0);
					double cx = m_solids.Local_V3f[patch_id][cur_local_id].x;
					double cy = m_solids.Local_V3f[patch_id][cur_local_id].y;
					double cz = m_solids.Local_V3f[patch_id][cur_local_id].z;
					model.control_vertices[cur_global_id] = Eigen::Vector3d(cx, cy, cz);
				}
			}
		}
		model.global_id.push_back(patch_global_ids);
		assert(model.global_id.size() == 64);
	}
	std::cout << "model generate" << std::endl;
	LoopData PP(model, ccbdy);
	PP.aroundcell = aroundcell;
	PP.aroundcell_weight = aroundcell_weight;
	PP.m_solids = m_solids;
	LoopData* PP_pointer = &PP;

	// 初始化密度 和 梯度
	double DensityChange = 100.0;
	int LoopTime = 1;
	std::vector<double> Density(mesh->cells_.size(), 1.0);
	std::vector<double> DensityPre(mesh->cells_.size(), 1.0);
	std::vector<double> dc(mesh->cells_.size(), 0.0);
	int NumOfVar = mesh->cells_.size();

	double f = 0.0;
	std::cout << "开始迭代" << std::endl;
	std::vector<int> isRemoved(mesh->cells_.size(), 0);
	while (DensityChange > 0.1 && LoopTime < MaxLoopTime) {
		//计算梯度信息
		f = single_resolution_topopt_func_CatmullClark_GPUSolver(Density, dc, PP_pointer);
		std::cout << "single_resolution_topopt_func end" << std::endl;
		//迭代法
		//RemoveOneStep(Density, dc, isProtectCell, isRemoved);
		OC(Density, TargetRate, dc, isProtectCell);
		//OC(Density, TargetRate, dc);
		DensityChange = Diff(Density, DensityPre);
		DensityPre = Density;
		std::cout << "第" << LoopTime << "次循环后DensityChange为：" << DensityChange << std::endl;
		std::cout << "第" << LoopTime << "次循环后柔度为：" << f << std::endl;
		std::string FileName(OutPutFile + std::to_string(LoopTime) + ".vtk");
		SaveCurMesh(Density, mesh, FileName, is_need_update_topopt_model);
		string qstr = string("第") + to_string(LoopTime) + string("/") + to_string(MaxLoopTime) + "次迭代";
		std::cout << "第" << LoopTime << "次模型保存完毕" << std::endl;
		LoopTime++;
		filenum++;
	}
	widget->hex_density = Density;
	is_need_update_topopt_model = false;

	std::cout << "拓扑优化结束" << std::endl;

}