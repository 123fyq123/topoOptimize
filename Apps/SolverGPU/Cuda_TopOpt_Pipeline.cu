//
// Created by XieJin on 10/27/2023.
//
#include "Cuda_TopOpt_Pipeline.h"
#include<Hyperelasticity/utils.h>
#include<Hyperelasticity/TimeIntegrator.h>
#include<chrono>
#include<Hyperelasticity/collision_detect.h>
#include <iostream>
#include "Util/GP_Log.h"
#include "Util/GP_Macros.h"
#include "Util/GP_Util.h"
#include "GP_Local.h"
#include "GP_Patch_Info.h"
#include "GP_Context.h"
#include "Patcher/GP_Patcher.h"
#include "GP_Mesh_Static.h"
#include "GP_vertex_normal_kernel.cuh"
#include "GP_Spline.h"
#include "GP_Spline_Static.h"
#include "GP_Spline_Test.h"
#include "GP_Spline_v0.h"
#include "GP_Invoke.h"
#include "GP_Spline_t.h"
#include "GP_Param.h"
#include "GP_BoundaryCPoint.h"
#include "Poster/source/Timer.h"
#include "omp.h"
#include "Poster/source/nanoflann.h"
#include <algorithm>
#include <vector>
#include <map>
#include <Eigen/Dense>
#include "MMASolver.h"
#include <fstream>
namespace Cuda_TopOpt_Pipeline_nsp {
	using namespace nanoflann;
	struct Point {
		double x, y, z;
	};
	template <typename T>
	struct PointCloud
	{

		using coord_t = T;  //!< The type of each coordinate

		std::vector<Point> pts;

		// Must return the number of data points
		inline size_t kdtree_get_point_count() const { return pts.size(); }

		// Returns the dim'th component of the idx'th point in the class:
		// Since this is inlined and the "dim" argument is typically an immediate
		// value, the
		//  "if/else's" are actually solved at compile time.
		inline T kdtree_get_pt(const size_t idx, const size_t dim) const
		{
			if (dim == 0)
				return pts[idx].x;
			else if (dim == 1)
				return pts[idx].y;
			else
				return pts[idx].z;
		}

		// Optional bounding-box computation: return false to default to a standard
		// bbox computation loop.
		//   Return true if the BBOX was already computed by the class and returned
		//   in "bb" so it can be avoided to redo it again. Look at bb.size() to
		//   find out the expected dimensionality (e.g. 2 or 3 for point clouds)
		template <class BBOX>
		bool kdtree_get_bbox(BBOX& /* bb */) const
		{
			return false;
		}
	};

	typedef KDTreeSingleIndexAdaptor<
		L2_Simple_Adaptor<double, PointCloud<double>>,
		PointCloud<double>,
		3 /* dim */
	> KDTree;

	void GetNearestIds(std::vector<Point>& p, std::vector<Point>& q, std::vector<int>& ans) {
		ans.resize(p.size());
		PointCloud<double> A, B;
		A.pts = p;
		B.pts = q;
		KDTree index(3, B, KDTreeSingleIndexAdaptorParams(10 /* max leaf */));
		index.buildIndex();
		int a_index = 0;
		for (const auto& a : A.pts) {
			size_t ret_index;
			double out_dist_sqr;
			nanoflann::KNNResultSet<double> resultSet(1);
			resultSet.init(&ret_index, &out_dist_sqr);
			index.findNeighbors(resultSet, &a.x);
			ans[a_index++] = ret_index;
		}
	}

}

void GPU_MMA(
	vector<double>& designVariables,        // 设计变量数组
	const vector<double>& df0dx,            // 目标函数灵敏度数组
	const vector<double>& volfaceDelta,           // 目标体积分数
	const vector<double>& dfdx,             // 约束灵敏度数组
	const vector<double>& xmin,             // 下限数组
	const vector<double>& xmax)            // 上限数组         
{
	int nDesign = designVariables.size();
	MMASolver* mmaSolver = new MMASolver(nDesign, 1); // 1个约束（体积）
	mmaSolver->Update(designVariables.data(), df0dx.data(), volfaceDelta.data(), dfdx.data(), xmin.data(), xmax.data());
	delete mmaSolver;
}

void GPU_OC(std::vector<double>& Density, const double TargetRate, const std::vector<double>& dc, const std::vector<int>& isPicked) {
	double a = 0.0;
	double b = 10000000000;
	double MoveStep = 0.2;

	std::vector<double> TempDensity = Density;

	int max_iter = 50;
	while (b - a > 1e-10 && max_iter-- > 0) {
		std::cout << "[OC]: " << "(a,b) = " << a << "," << b << "\n";
		double Lmid = 0.5 * (a + b);
		double SumDensity = 0.0;
		for (int i = 0; i < Density.size(); ++i) {
			if (isPicked[i] == 1) TempDensity[i] = 1;
			else TempDensity[i] = fmax(0.0001, fmax(Density[i] - MoveStep, fmin(1.0, fmin(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
			//TempDensity[i] = fmax(0.001, fmax(Density[i] - MoveStep, fmin(1.0, fmin(Density[i] + MoveStep, Density[i] * sqrt(-dc[i] / Lmid)))));
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
	std::cout << "TargetRate = " << TargetRate << std::endl;
	Density = TempDensity;
}

void Cuda_TopOpt_Pipeline(const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model) {
	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;

	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml");

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;

	// add bdy info
	{
		int control_size = 0;
		std::ifstream in(bdy_file);
		in >> control_size;
		std::cout << "control_size = " << control_size << "\n";

		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			CPoint cp(x, y, z);
			CVector cv(vx, vy, vz);
			bdy_points.push_back(cp);
			bdy_value.push_back(cv);
			bdy_type.push_back(cur_type);
		}
	}
	timer.out("add bdy info");
	gpmesh::GPLog::init(spdlog::level::trace);

	CUDA_ERROR(cudaDeviceReset());
	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;

	timer.out("deal with bdy info");
	std::vector<double> dc_pre(patch.size(), 0);
	auto exist_num = (double)patch.size();
	auto tot_num = (double)patch.size();
	//while(1){
	// data init
	gpmesh::GPSplineSolidV0 m_spline;
	//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
	//    timer.out("TopOpt_GPU_Init");
	m_spline.init_patch(patch);
	timer.out("TopOpt_GPU_Geometry_Init");

	timer.out("before sys");
	const SparseSystem& temp_sys = m_spline.m_sparse_system;
	timer.out("after sys");
	{
		std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
		std::vector<int> global_ids;

		for (auto& p : bdy_points) {
			A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}

		GetNearestIds(A, B, global_ids);
		timer.out("after GetNearestIds");
		std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
		auto& global_local_map = temp_sys.get_l2g_map();

		for (int i = 0; i < bdy_points.size(); ++i) {
			auto& p = bdy_points[i];
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id];

			if (global_id == -1) {
				std::cerr << "cant find control points\n";
				//std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
				std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
			}
			else {
				if (bdy_type[i] == 2) {
					neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
				}
				else if (bdy_type[i] == 1) {
					dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
				}
			}
		}
	}

	m_spline.init_dirichlet_points(dirichlet_points);

	std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
	// first assembly
	m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
	timer.out("Get_Local_Matrix");
	m_spline.m_sparse_system.rhs().setZero();

	m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
	m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


	int topopt_iter = 0;

	std::vector<int> isPicked(patch.size(), 0);
	for (auto& xx : protect_cell_ids) {
		isPicked[xx] = 1;
	}

	while (topopt_iter < max_iter) {

		std::cout << topopt_iter++ << std::endl;
		int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
		if (reduce_num < 1) break;
		Eigen::MatrixXd U;
		/// todo: need to check


		m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

		Eigen::MatrixXd UU(patch.size(), 64 * 3);

		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

				UU(p_id, local_id) = U(global_id, 0);
				UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
				UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
			}
		}

		std::vector<double> Ce(patch.size());
		std::vector<double> dc(patch.size());
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
			if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
			dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
		}


		// fliter
		{
			auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
			int n = patch.size();                        //�ܵ����ɶ���Ŀ
			for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
				double weight = 0;
				double temp_dx = 0;
				int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
				for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
					weight += neighbor_info[i][j].second;
					temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
				}
				temp_dx /= weight;
				Tempdc[i] = temp_dx;          //temp_dc��ֵ
			}
			dc = Tempdc;       //������Ľ������
		}

		//        {
		//            for(int i = 0; i < tot_num; ++i){
		//                dc[i] = (dc[i] + dc_pre[i]) / 2;
		//                dc_pre[i] = dc[i];
		//            }
		//        }

		GPU_OC(density, target_rate, dc, isPicked);

		//        std::vector<std::pair<double, int>> dc_ids;
		//        for(int i = 0; i < patch.size(); ++i){
		//            dc_ids.emplace_back(-dc[i], i);
		//            //std::cout <<density[i] << " " <<dc[i] <<" ";
		//        }
		//
		//        std::sort(dc_ids.begin(), dc_ids.end());
		//
		//        int remove_num = 0;
		//
		//        int kk = -1;
		//        while(remove_num < reduce_num && kk < (int)dc_ids.size()){
		//            kk++;
		//            int id = dc_ids[kk].second;
		//           //std::cout <<"dc_ids_id = " << id << std::endl;
		//            if(density[id] < 0.8) continue;
		//            //std::cout <<"NOOOOOOOOOOOOOOO\n";
		//            if(std::find(protect_cell_ids.begin(), protect_cell_ids.end(), id) != protect_cell_ids.end()) continue;
		//           //std::cout <<"YESSSSSSSSSSSSSSSSSSS\n";
		//            density[id] = 0.0001;
		//            remove_num++;
		//            exist_num--;
		//        }
		//        std::cout << "reduce_num = " << reduce_num << std::endl;
		//        std::cout << exist_num << std::endl;
		is_need_update_model = true;
	}
}


void Cuda_TopOpt_Discrete(const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model) {
	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;

	timer.out("read_xml start");
	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml end");

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;

	// add bdy info
	{
		int control_size = 0;
		std::ifstream in(bdy_file);
		in >> control_size;
		std::cout << "control_size = " << control_size << "\n";

		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			CPoint cp(x, y, z);
			CVector cv(vx, vy, vz);
			bdy_points.push_back(cp);
			bdy_value.push_back(cv);
			bdy_type.push_back(cur_type);
		}
	}
	gpmesh::GPLog::init(spdlog::level::trace);

	CUDA_ERROR(cudaDeviceReset());
	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;

	timer.out("deal with bdy info");
	std::vector<double> dc_pre(patch.size(), 0);
	auto exist_num = (double)patch.size();
	auto tot_num = (double)patch.size();
	//while(1){
	// data init
	gpmesh::GPSplineSolidV0 m_spline;
	//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
	//    timer.out("TopOpt_GPU_Init");
	timer.out("Init_patch start");
	m_spline.init_patch(patch);
	timer.out("Init_patch end");

	const SparseSystem& temp_sys = m_spline.m_sparse_system;
	{
		std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
		std::vector<int> global_ids;

		for (auto& p : bdy_points) {
			A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}

		GetNearestIds(A, B, global_ids);
		std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
		auto& global_local_map = temp_sys.get_l2g_map();

		for (int i = 0; i < bdy_points.size(); ++i) {
			auto& p = bdy_points[i];
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id];

			if (global_id == -1) {
				std::cerr << "cant find control points\n";
				//std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
				std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
			}
			else {
				if (bdy_type[i] == 2) {
					neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
				}
				else if (bdy_type[i] == 1) {
					dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
				}
			}
		}
	}

	m_spline.init_dirichlet_points(dirichlet_points);

	std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
	// first assembly
	m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
	timer.out("Get_Local_Matrix");
	m_spline.m_sparse_system.rhs().setZero();

	m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
	m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


	int topopt_iter = 0;

	std::vector<int> isPicked(patch.size(), 0);
	for (auto& xx : protect_cell_ids) {
		isPicked[xx] = 1;
	}

	Eigen::MatrixXd U;
	m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

	Eigen::MatrixXd UU(patch.size(), 64 * 3);

	for (int p_id = 0; p_id < patch.size(); ++p_id) {
		for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
			int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

			UU(p_id, local_id) = U(global_id, 0);
			UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
			UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
		}
	}

	std::vector<double> Ce(patch.size());
	std::vector<double> dc(patch.size());
	for (int p_id = 0; p_id < patch.size(); ++p_id) {
		Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
		if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
		dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
	}


	// fliter
	{
		auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
		int n = patch.size();                        //�ܵ����ɶ���Ŀ
		for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
			double weight = 0;
			double temp_dx = 0;
			int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
			for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
				weight += neighbor_info[i][j].second;
				temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
			}
			temp_dx /= weight;
			Tempdc[i] = temp_dx;          //temp_dc��ֵ
		}
		dc = Tempdc;       //������Ľ������
	}

	double maxdc = fabs(dc[0]);
	for (int i = 1; i < dc.size(); i++) {
		if (fabs(dc[i]) > maxdc) {
			maxdc = fabs(dc[i]);
		}
	}
	for (int i = 0; i < dc.size(); i++) {
		dc[i] /= maxdc;
	}

	double c_d;
	double c_old = 0, c_new = 0, sumnew = 0;
	c_d = density.size();

	double beta, v;
	std::vector<double> sigma;
	std::vector<double> olddensity = density;
	std::vector<double> newdensity = density;

	sigma.resize(density.size(), 0);
	beta = 2000;
	v = 1;
	double lambda = 1e-3;

	double w1 = 1e-8, Vr;
	Vr = density.size() * 0.98;

	bool flag = false;
	int pareto = 0;

	timer.out("Iteration start");

	Timer itertimer;
	while (topopt_iter < max_iter) {

		std::cout << topopt_iter++ << std::endl;
		int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
		if (reduce_num < 1) break;

		/// todo: need to check

		olddensity = density;



		while (true) {
			double suma = 0, sumb = 0;
			std::vector<double> s(density.size());
			for (int i = 0; i < density.size(); ++i) {
				s[i] = lambda + dc[i];
			}

			for (int i = 0; i < s.size(); ++i) {
				if (fabs(s[i]) <= 1E-8) {
					s[i] = 1E-8 * (s[i] >= 0 ? 1 : -1);
				}
			}
			// 计算 delta
			std::complex<double> delta = pow(beta, 2) / 27;
			// 计算 r
			std::vector<std::complex<double>> r(s.size());
			for (int i = 0; i < s.size(); ++i) {
				std::complex<double> sqrt_term = sqrt(std::complex<double>(pow(s[i], 2) * (pow(s[i], 2) - delta)));
				r[i] = pow(delta, -1.0 / 3) * pow(std::complex<double>(2 * pow(s[i], 2) - delta + sqrt_term + sqrt_term), 1.0 / 3);
			}
			// 计算 sigma
			for (int i = 0; i < sigma.size(); ++i) {
				sigma[i] = (1.0 / 6) * beta * (-1 + r[i].real() * 2);
				suma += v - dc[i] / sigma[i];
				sumb += 1 / sigma[i];
			}

			lambda = (suma - 2 * Vr) / sumb;

			c_old = 0;
			c_new = 0;
			sumnew = 0;
			for (int i = 0; i < density.size(); i++) {
				if (isPicked[i] == 1) {
					newdensity[i] = 1;
				}
				else {
					newdensity[i] = 0.5 * (1.0 - (lambda + dc[i]) / sigma[i]);
				}

				sumnew += newdensity[i];
				c_old += olddensity[i] * dc[i];
				c_new += newdensity[i] * dc[i];
			}
			if (fabs(c_old - c_new) / fabs(c_old) <= w1) {
				break;
			}
			else {
				olddensity = newdensity;
			}
		}

		for (int i = 0; i < newdensity.size(); i++) {
			if (newdensity[i] <= 0) {
				newdensity[i] = 0.001;
			}
			else if (newdensity[i] >= 1) {
				newdensity[i] = 1;
			}

			newdensity[i] = round(newdensity[i]);
		}

		//density = newdensity;
		//Vr *= 0.9;
		m_spline.Solver_With_Density(neumann_points, dirichlet_points, newdensity, U);

		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

				UU(p_id, local_id) = U(global_id, 0);
				UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
				UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
			}
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
			if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
			dc[p_id] = (-3 * newdensity[p_id] * newdensity[p_id] * (1 - 0.0001)) * Ce[p_id];
		}

		// fliter
		{
			auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
			int n = patch.size();                        //�ܵ����ɶ���Ŀ
			for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
				double weight = 0;
				double temp_dx = 0;
				int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
				for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
					weight += neighbor_info[i][j].second;
					temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
				}
				temp_dx /= weight;
				Tempdc[i] = temp_dx;          //temp_dc��ֵ
			}
			dc = Tempdc;       //������Ľ������
		}

		double maxdc = fabs(dc[0]);
		for (int i = 1; i < dc.size(); i++) {
			if (fabs(dc[i]) > maxdc) {
				maxdc = fabs(dc[i]);
			}
		}
		for (int i = 0; i < dc.size(); i++) {
			dc[i] /= maxdc;
		}

		if (fabs(c_new - c_d) / fabs(c_d) > 5e-3) {
			density = newdensity;
			c_d = c_new;
		}
		else {
			if (flag) {
				break;
			}
			density = newdensity;
			if (Vr / density.size() <= target_rate) {
				Vr = density.size() * max(0.3, Vr / density.size());
				flag = true;
				continue;
			}
			else {
				Vr *= 0.98;
			}

		}

		GPU_OC(density, target_rate, dc, isPicked);

		is_need_update_model = true;
	}
	std::cout << "Finsh" << std::endl;
	itertimer.out("Iteration time");
}


void Cuda_TopOpt_Multi_Pipeline(const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	int sub_time,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model) {
	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;

	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml");

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;

	// add bdy info
	{
		int control_size = 0;
		std::ifstream in(bdy_file);
		in >> control_size;
		std::cout << "control_size = " << control_size << "\n";

		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			CPoint cp(x, y, z);
			CVector cv(vx, vy, vz);
			bdy_points.push_back(cp);
			bdy_value.push_back(cv);
			bdy_type.push_back(cur_type);
		}
	}
	timer.out("add bdy info");
	gpmesh::GPLog::init(spdlog::level::trace);

	CUDA_ERROR(cudaDeviceReset());
	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;

	timer.out("deal with bdy info");
	std::vector<double> dc_pre(patch.size(), 0);
	auto exist_num = (double)patch.size();
	auto tot_num = (double)patch.size();

	{
		assert(sub_time < 5);
		int sub_size = 1 << sub_time;
		std::vector<std::vector<Eigen::MatrixXd>> sub_local_matrix(sub_size, std::vector<Eigen::MatrixXd>(patch.size(), Eigen::MatrixXd::Zero(192, 192)));

	}





	//while(1){
	// data init
	gpmesh::GPSplineSolidV0 m_spline;
	//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
	//    timer.out("TopOpt_GPU_Init");
	m_spline.init_patch(patch);
	timer.out("TopOpt_GPU_Geometry_Init");

	timer.out("before sys");
	const SparseSystem& temp_sys = m_spline.m_sparse_system;
	timer.out("after sys");
	{
		std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
		std::vector<int> global_ids;

		for (auto& p : bdy_points) {
			A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}

		GetNearestIds(A, B, global_ids);
		timer.out("after GetNearestIds");
		std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
		auto& global_local_map = temp_sys.get_l2g_map();

		for (int i = 0; i < bdy_points.size(); ++i) {
			auto& p = bdy_points[i];
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id];

			if (global_id == -1) {
				std::cerr << "cant find control points\n";
				//std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
				std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
			}
			else {
				if (bdy_type[i] == 2) {
					neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
				}
				else if (bdy_type[i] == 1) {
					dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
				}
			}
		}
	}

	m_spline.init_dirichlet_points(dirichlet_points);

	std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192));
	// first assembly
	m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
	timer.out("Get_Local_Matrix");
	m_spline.m_sparse_system.rhs().setZero();

	m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
	m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);


	int topopt_iter = 0;

	std::vector<int> isPicked(patch.size(), 0);
	for (auto& xx : protect_cell_ids) {
		isPicked[xx] = 1;
	}

	while (topopt_iter < max_iter) {

		std::cout << topopt_iter++ << std::endl;
		int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
		if (reduce_num < 1) break;
		Eigen::MatrixXd U;
		/// todo: need to check


		m_spline.Solver_With_Density(neumann_points, dirichlet_points, density, U);

		Eigen::MatrixXd UU(patch.size(), 64 * 3);

		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

				UU(p_id, local_id) = U(global_id, 0);
				UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
				UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
			}
		}

		std::vector<double> Ce(patch.size());
		std::vector<double> dc(patch.size());
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
			if (Ce[p_id] < 0) std::cout << "Ce is wrong!\n";
			dc[p_id] = (-3 * density[p_id] * density[p_id] * (1 - 0.0001)) * Ce[p_id];
		}


		// fliter
		{
			auto Tempdc = dc;                         //����һ����ʱ�洢dc�ı���
			int n = patch.size();                        //�ܵ����ɶ���Ŀ
			for (int i = 0; i < n; ++i) {             //ѭ��ÿһ��С��
				double weight = 0;
				double temp_dx = 0;
				int aroundnum = neighbor_info[i].size();        //��ǰ����Χ��cell��Ŀ
				for (int j = 0; j < aroundnum; ++j) {        //������Щ��Χ��cell��ȡ��weight�������Ӧ�Ŀ���ˣ��ӵ���ǰ��ר����temp_dc��
					weight += neighbor_info[i][j].second;
					temp_dx += dc[neighbor_info[i][j].first] * neighbor_info[i][j].second;
				}
				temp_dx /= weight;
				Tempdc[i] = temp_dx;          //temp_dc��ֵ
			}
			dc = Tempdc;       //������Ľ������
		}

		//        {
		//            for(int i = 0; i < tot_num; ++i){
		//                dc[i] = (dc[i] + dc_pre[i]) / 2;
		//                dc_pre[i] = dc[i];
		//            }
		//        }

		GPU_OC(density, target_rate, dc, isPicked);

		//        std::vector<std::pair<double, int>> dc_ids;
		//        for(int i = 0; i < patch.size(); ++i){
		//            dc_ids.emplace_back(-dc[i], i);
		//            //std::cout <<density[i] << " " <<dc[i] <<" ";
		//        }
		//
		//        std::sort(dc_ids.begin(), dc_ids.end());
		//
		//        int remove_num = 0;
		//
		//        int kk = -1;
		//        while(remove_num < reduce_num && kk < (int)dc_ids.size()){
		//            kk++;
		//            int id = dc_ids[kk].second;
		//           //std::cout <<"dc_ids_id = " << id << std::endl;
		//            if(density[id] < 0.8) continue;
		//            //std::cout <<"NOOOOOOOOOOOOOOO\n";
		//            if(std::find(protect_cell_ids.begin(), protect_cell_ids.end(), id) != protect_cell_ids.end()) continue;
		//           //std::cout <<"YESSSSSSSSSSSSSSSSSSS\n";
		//            density[id] = 0.0001;
		//            remove_num++;
		//            exist_num--;
		//        }
		//        std::cout << "reduce_num = " << reduce_num << std::endl;
		//        std::cout << exist_num << std::endl;
		is_need_update_model = true;
	}
}


void Cuda_Linear_Elastic_Solver(const std::string& input_model_file, const std::string& boundary_file,
	const std::string& output_file, double Yangs, double Poissons) {
	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;
	bool bH = true;


	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml");

	bool bUseBoundaryPoint = true;

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;

	// add bdy info
	{
		int control_size = 0;
		std::ifstream in(boundary_file);
		in >> control_size;
		std::cout << "control_size = " << control_size << "\n";

		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			CPoint cp(x, y, z);
			CVector cv(vx, vy, vz);
			bdy_points.push_back(cp);
			bdy_value.push_back(cv);
			bdy_type.push_back(cur_type);
		}
	}
	timer.out("add bdy info");
	gpmesh::GPLog::init(spdlog::level::trace);

	CUDA_ERROR(cudaDeviceReset());
	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;
	gpmesh::GPSplineSolidV0 m_spline;
	//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
	//    timer.out("TopOpt_GPU_Init");
	m_spline.init_patch(patch);
	timer.out("Linear Elasticity GPU Geometry Init");

	timer.out("before sys");
	const SparseSystem& temp_sys = m_spline.m_sparse_system;
	timer.out("after sys");
	{
		std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
		std::vector<int> global_ids;

		for (auto& p : bdy_points) {
			A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}

		GetNearestIds(A, B, global_ids);
		timer.out("after GetNearestIds");
		std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
		auto& global_local_map = temp_sys.get_l2g_map();

		for (int i = 0; i < bdy_points.size(); ++i) {
			auto& p = bdy_points[i];
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id];

			if (global_id == -1) {
				std::cerr << "cant find control points\n";
				//std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
				std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
			}
			else {
				if (bdy_type[i] == 2) {
					neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
				}
				else if (bdy_type[i] == 1) {
					dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
				}
			}
		}
	}

	m_spline.init_dirichlet_points(dirichlet_points);

	Eigen::MatrixXd UU;
	m_spline.debug_IGA_linear_solver(patch, neumann_points, dirichlet_points, Yangs, Poissons, UU);
	timer.out("solver");

	stiffAssembler stiff(patch, dirichlet_points, neumann_points);

	Eigen::MatrixXd originUU = Eigen::MatrixXd::Zero(stiff.system().getControlNums() * stiff.system().getDim(), 1);

	std::string origin_output_path = std::string(IGA_OUTPUT_DIR) + "Origin/Origin_";

	stiff.WriteParaview(patch, origin_output_path, UU,
		1000, false);

	std::string output_path = std::string(IGA_OUTPUT_DIR) + "Test/test_linear" + to_string(0) + "_";

	Eigen::MatrixXd displacement1;

	int begin = 0;
	stiff.constructCoefs(UU, displacement1);

	for (size_t i = 0; i < patch.size(); ++i) {
		int num = patch[i].getSizeControlPoints();
		patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
		begin += num;
	}

	stiff.WriteParaview(patch, output_path, UU,
		1000, false);

	begin = 0;
	displacement1 *= -1;
	for (size_t i = 0; i < patch.size(); ++i) {
		int num = patch[i].getSizeControlPoints();
		patch[i].updateCoefs(displacement1.block(begin, 0, num, 3));
		begin += num;
	}

	// 输出
	{
		std::ofstream out(output_file, std::ofstream::out | std::ofstream::trunc);
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];

				out << UU(global_id, 0) << " " << UU(global_id + m_spline.m_sparse_system.getControlNums(), 0)
					<< " " << UU(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0) << "\n";
			}
		}
		out.close();
	}
}



// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// -----------------------------------------------3DMMC------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
// ----------------------------------------------------------------------------------------------------------
/// 3DMMC组件结构体 - 表示一个椭球体杆件
struct MMC_Component {
	// 中心位置坐标
	double x0, y0, z0;

	// 三个半轴长度（长轴、半长轴、短轴）
	double l1, l2, l3;

	// 三个欧拉旋转角（绕x、y、z轴的旋转角度，单位：弧度）
	double alpha, beta, gamma;

	/// 默认构造函数
	MMC_Component() : x0(0), y0(0), z0(0),
		l1(0), l2(0), l3(0),
		alpha(0), beta(0), gamma(0) {
	}

	/// 带参数的构造函数
	MMC_Component(double x, double y, double z,
		double length1, double length2, double length3,
		double rot_alpha, double rot_beta, double rot_gamma)
		: x0(x), y0(y), z0(z),
		l1(length1), l2(length2), l3(length3),
		alpha(rot_alpha), beta(rot_beta), gamma(rot_gamma) {
	}
};

const int p = 6;
const double EPS = 1e-9;
// 1 * 1 单元
double Heaviside(double phi, double alpha, double epsilon) {
	if (phi > epsilon) return 1.0;
	if (phi < -epsilon) return alpha;
	return 3 * (1 - alpha) / 4 * (phi / epsilon - pow(phi, 3) / (3 * pow(epsilon, 3))) + (1 + alpha) / 2;
}
// xInt: x杆件放置间隔
std::vector<MMC_Component> init(const std::vector<double> xyzRange,
	double xInt, double yInt, double zInt,
	const std::vector<double>& vInt, int type = -1)
{
	std::vector<MMC_Component> comps;
	if (vInt.size() < 6) return comps;

	double DL = xyzRange[1] - xyzRange[0];
	double DW = xyzRange[3] - xyzRange[2];
	double DH = xyzRange[5] - xyzRange[4];

	// 生成 x,y,z 采样位置：从 initial 到 <= D，步长 2*Int
	std::vector<double> xs, ys, zs;
	int lx = xyzRange[0] + xInt, rx = xyzRange[1];
	int ly = xyzRange[2] + yInt, ry = xyzRange[3];
	int lz = xyzRange[4] + zInt, rz = xyzRange[5];
	for (double x = lx; x <= rx; x += 2.0 * xInt) xs.push_back(x);
	for (double y = ly; y <= ry; y += 2.0 * yInt) ys.push_back(y);
	for (double z = lz; z <= rz; z += 2.0 * zInt) zs.push_back(z);

	// 每个格点放置 8 个杆件，beta/gamma 的符号按 Matlab 模式：
	// bet signs:  [ 1 -1  1 -1  1 -1  1 -1 ] * vInt[4]
	// gam signs:  [ 1  1 -1 -1  1  1 -1 -1 ] * vInt[5]
	const int patternN = 8;
	const double alphaSigns[patternN] = { 1.0,  1.0, 1.0, 1.0,  -1.0,  -1.0, -1.0, -1.0 };
	const double betSigns[patternN] = { 1.0, -1.0,  1.0, -1.0,  1.0, -1.0,  1.0, -1.0 };
	const double gamSigns[patternN] = { 1.0,  1.0, -1.0, -1.0,  1.0,  1.0, -1.0, -1.0 };
	// 读取 vInt：l1,l2,l3, alpha, beta_amp, gamma_amp
	double l1 = vInt[0];
	double l2 = vInt[1];
	double l3 = vInt[2];
	double alp_tmp = vInt[3];
	double bet_tmp = vInt[4];
	double gam_tmp = vInt[5];

	for (size_t iz = 0; iz < zs.size(); ++iz) {
		for (size_t iy = 0; iy < ys.size(); ++iy) {
			for (size_t ix = 0; ix < xs.size(); ++ix) {
				if (type != -1) {
					double x0 = xs[ix];
					double y0 = ys[iy];
					double z0 = zs[iz];
					//for (int t = 0; t < 2; ++t) {
					MMC_Component c1(x0, y0, z0,
						l1, l2, l3,
						alp_tmp, bet_tmp, gam_tmp);
					comps.push_back(c1);
					MMC_Component c2(x0, y0, z0,
						l1, l2, l3,
						-alp_tmp, bet_tmp, gam_tmp);
					comps.push_back(c2);
					MMC_Component c3(x0, y0, z0,
						l1, l2, l3,
						-alp_tmp, bet_tmp, gam_tmp);
					comps.push_back(c3);
					MMC_Component c4(x0, y0, z0,
						l1, l2, l3,
						-alp_tmp, bet_tmp, gam_tmp);
					comps.push_back(c4);
					//}
				}
				else
				{
					for (int i = 0; i < patternN; ++i) {
						double x0 = xs[ix];
						double y0 = ys[iy];
						double z0 = zs[iz];
						auto alp = alp_tmp * alphaSigns[i];
						double bet = bet_tmp * betSigns[i];
						double gam = gam_tmp * gamSigns[i];
						MMC_Component c(x0, y0, z0,
							l1, l2, l3,
							alp, bet, gam);
						comps.push_back(c);
					}
				}
			}
		}
	}

	return comps;
}

void calc_Phi_Node(
	const std::map<int, Eigen::Vector3d>& allVertex, // 网格中心
	const double minL,
	std::vector<MMC_Component> components, // 所有组件
	vector<double>& componentPhi, // 第i个组件的phi值  节点
	std::vector<std::vector<double>>& componentPhiDrv, // 第i个组件的phi导数  节点
	int i, // 第i个组件
	int p, // TDF的p值
	double epsilon)
{
	int nVertex = allVertex.size();
	componentPhi.resize(nVertex);
	componentPhiDrv.resize(nVertex, std::vector<double>(9, 0.0));
	double x0 = components[i].x0;
	double y0 = components[i].y0;
	double z0 = components[i].z0;
	double l1 = components[i].l1 + std::numeric_limits<double>::epsilon();
	double l2 = components[i].l2 + std::numeric_limits<double>::epsilon();
	double l3 = components[i].l3 + std::numeric_limits<double>::epsilon();
	double alp = components[i].alpha;
	double bet = components[i].beta;
	double gam = components[i].gamma;
	//if (l1 <= minL || l2 <= minL || l3 <= minL) return;
	// 旋转矩阵
	double sa = std::sin(alp), ca = std::cos(alp);
	double sb = std::sin(bet), cb = std::cos(bet);
	double sg = std::sin(gam), cg = std::cos(gam);

	Eigen::Matrix3d R;
	R << cb * cg, cb* sg, -sb,
		sa* sb* cg - ca * sg, sa* sb* sg + ca * cg, sa* cb,
		ca* sb* cg + sa * sg, ca* sb* sg - sa * cg, ca* cb;

	// 为数值稳定加一个小常数
	const double tiny = 1e-12;

	// 存储所有网格点的坐标差值和旋转后坐标
	std::vector<Eigen::Vector3d> lc_vec(nVertex);
	std::vector<Eigen::Vector3d> xyz_vec(nVertex);
	// 存储所有网格点的旋转后坐标
	std::vector<double> x1_vec(nVertex), y1_vec(nVertex), z1_vec(nVertex);
	// x1/l1+x2/l2+x3/l3
	std::vector<double> temp_vec(nVertex);

	for (auto vertexPair : allVertex)
	{
		int vertexTag = vertexPair.first;
		Eigen::Vector3d lc = vertexPair.second - Eigen::Vector3d(x0, y0, z0); // 局部向量
		Eigen::Vector3d xyz = R * lc;
		lc_vec[vertexTag] = lc;
		double x1 = xyz[0] + tiny;
		double y1 = xyz[1] + tiny;
		double z1 = xyz[2] + tiny;
		xyz_vec[vertexTag] = xyz;
		x1_vec[vertexTag] = x1;
		y1_vec[vertexTag] = y1;
		z1_vec[vertexTag] = z1;

		double t1 = std::pow(x1 / l1, p);
		double t2 = std::pow(y1 / l2, p);
		double t3 = std::pow(z1 / l3, p);

		double temp = t1 + t2 + t3;
		temp_vec[vertexTag] = temp;
		// 公式2
		double phi = 1.0 - std::pow(temp, 1.0 / p);

		componentPhi[vertexTag] = phi;
	}


	// 计算旋转矩阵对 α, β, γ 的导数
	Eigen::Matrix3d Ra, Rb, Rg;
	// 公式23
	Ra << 0, 0, 0,
		R(2, 0), R(2, 1), R(2, 2),
		-R(1, 0), -R(1, 1), -R(1, 2);

	// 公式24
	Rb << -sb * cg, -sb * sg, -cb,
		sa* cb* cg, sa* cb* sg, -sa * sb,
		ca* cb* cg, ca* cb* sg, -ca * sb;

	// 公式25
	Rg << -cb * sg, cb* cg, 0,
		-sa * sb * sg - ca * cg, sa* sb* cg - ca * sg, 0,
		-ca * sb * sg + sa * cg, ca* sb* cg + sa * sg, 0;

	for (auto vertexPair : allVertex) {
		int vertexTag = vertexPair.first;
		const double x1 = x1_vec[vertexTag];
		const double y1 = y1_vec[vertexTag];
		const double z1 = z1_vec[vertexTag];
		const double temp = temp_vec[vertexTag];
		const Eigen::Vector3d& lc = lc_vec[vertexTag];  // 使用引用避免拷贝

		const double temp1 = std::pow(temp, 1.0 / p - 1.0) * std::pow(x1 / l1, p - 1);
		const double temp2 = std::pow(temp, 1.0 / p - 1.0) * std::pow(y1 / l2, p - 1);
		const double temp3 = std::pow(temp, 1.0 / p - 1.0) * std::pow(z1 / l3, p - 1);

		// 式26 R11/R12...前面的系数
		const double dpdx1 = temp1 / l1;
		const double dpdy1 = temp2 / l2;
		const double dpdz1 = temp3 / l3;

		// 对半轴长度的导数
		const double dpdl1 = temp1 * (x1 / (l1 * l1));
		const double dpdl2 = temp2 * (y1 / (l2 * l2));
		const double dpdl3 = temp3 * (z1 / (l3 * l3));

		// 计算对旋转角度的导数分量
		const double dx1_dalpha = Ra.row(0) * lc;
		const double dy1_dalpha = Ra.row(1) * lc;
		const double dz1_dalpha = Ra.row(2) * lc;

		const double dx1_dbeta = Rb.row(0) * lc;
		const double dy1_dbeta = Rb.row(1) * lc;
		const double dz1_dbeta = Rb.row(2) * lc;

		const double dx1_dgamma = Rg.row(0) * lc;
		const double dy1_dgamma = Rg.row(1) * lc;
		const double dz1_dgamma = Rg.row(2) * lc;

		// 组装导数向量
		std::vector<double>& drv = componentPhiDrv[vertexTag];

		// 对中心坐标的导数
		drv[0] = dpdx1 * R(0, 0) + dpdy1 * R(1, 0) + dpdz1 * R(2, 0);
		drv[1] = dpdx1 * R(0, 1) + dpdy1 * R(1, 1) + dpdz1 * R(2, 1);
		drv[2] = dpdx1 * R(0, 2) + dpdy1 * R(1, 2) + dpdz1 * R(2, 2);

		// 对半轴长度的导数
		drv[3] = dpdl1;
		drv[4] = dpdl2;
		drv[5] = dpdl3;

		// 对旋转角度的导数
		drv[6] = -(dpdx1 * dx1_dalpha + dpdy1 * dy1_dalpha + dpdz1 * dz1_dalpha);
		drv[7] = -(dpdx1 * dx1_dbeta + dpdy1 * dy1_dbeta + dpdz1 * dz1_dbeta);
		drv[8] = -(dpdx1 * dx1_dgamma + dpdy1 * dy1_dgamma + dpdz1 * dz1_dgamma);
	}
}


void ComputeGlobalTDF(
	const std::vector<std::vector<double>>& allPhis, // 所有组件的Phi
	double lmd, // KS 参数
	std::vector<double>& globalPhi,
	std::vector<std::vector<double>>& phimaxDphi) // 全局Phi
{
	const int nComp = (int)allPhis.size();
	const int nGrid = (int)allPhis[0].size();

	globalPhi.assign(nGrid, 0.0);
	phimaxDphi.assign(nGrid, std::vector<double>(nComp, 0.0));
	const double tiny = 1e-300;

	for (int g = 0; g < nGrid; ++g) {
		std::vector<double> expValues(nComp);
		double sumExp = 0.0;
		for (int c = 0; c < nComp; ++c) {
			expValues[c] = std::exp(lmd * allPhis[c][g]);
			sumExp += expValues[c];
		}
		globalPhi[g] = std::max(-1e3, std::log(sumExp) / lmd);   // 公式4  K-S函数

		for (int c = 0; c < nComp; ++c) {
			phimaxDphi[g][c] = expValues[c] / (sumExp + tiny);   // 公式10  
		}
	}
}

void ComputeSensitivitiesByAnaytic(
	int iter,
	const double minL,
	const std::vector<double>& minx,
	const std::vector<double>& maxx,
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::vector<double>& originGlobalPhi,           // 全局Phi值
	const std::vector<MMC_Component>& mmcComponents,
	const Eigen::MatrixXd& U,                       // 位移向量
	const std::vector<double>& density,
	const std::vector<Eigen::MatrixXd>& local_matrix, // 局部刚度矩阵
	double alpha,                                   // Heaviside参数
	double epsilon,                                 // Heaviside参数
	double DL, double DW, double DH,                // 设计域尺寸
	std::vector<double>& df0dx,                     // 目标函数对设计变量的灵敏度
	std::vector<double>& dfdx)                      // 体积对设计变量的灵敏度
{
	const int nVertex = allVertex.size();
	const int nComp = mmcComponents.size();
	const int nVarsPerComp = 9; // 每个组件9个设计变量
	const int nGrid = cell2V.size();

	// 初始化灵敏度向量
	df0dx.assign(nComp * nVarsPerComp, 0.0);
	dfdx.assign(nComp * nVarsPerComp, 0.0);

	// 计算原始Heaviside值和应变能
	std::vector<double> originalVertexH(nVertex);
	std::vector<double> originalGridH(nGrid);
	for (int i = 0; i < nVertex; ++i) {
		originalVertexH[i] = Heaviside(originGlobalPhi[i], alpha, epsilon);
	}

	// 节点分配密度给单元
	for (int i = 0; i < nGrid; ++i)
	{
		auto pairVec = cell2V.at(i);
		originalGridH[i] = 0;
		for (int j = 0; j < 8; ++j)
		{
			int tag = pairVec[j].first;
			originalGridH[i] += originalVertexH[tag];
		}
		originalGridH[i] /= 8.0;
	}

	std::vector<double> Ce(nGrid);
	for (int p_id = 0; p_id < nGrid; ++p_id) {
		Ce[p_id] = U.row(p_id) * local_matrix[p_id] * U.row(p_id).transpose();
	}

	// 保存每个组件每个设计变量在每个vertex的dH值
	std::vector<std::vector<std::vector<double>>> dH_values(nComp,
		std::vector<std::vector<double>>(nVarsPerComp, std::vector<double>(nVertex, 0.0)));

	// 对每个组件的每个设计变量进行有限差分
	for (int compIdx = 0; compIdx < nComp; ++compIdx) {
		for (int varIdx = 0; varIdx < nVarsPerComp; ++varIdx) {
			// 创建正向扰动后的组件集合
			std::vector<MMC_Component> perturbedComponentsPos = mmcComponents;
			double range = maxx[varIdx] - minx[varIdx];
			double delta = range * 0.05;
			// 正向扰动 1.02倍
			switch (varIdx) {
			case 0: perturbedComponentsPos[compIdx].x0 = mmcComponents[compIdx].x0 + delta; break;
			case 1: perturbedComponentsPos[compIdx].y0 = mmcComponents[compIdx].y0 + delta; break;
			case 2: perturbedComponentsPos[compIdx].z0 = mmcComponents[compIdx].z0 + delta; break;
			case 3: perturbedComponentsPos[compIdx].l1 = mmcComponents[compIdx].l1 + delta; break;
			case 4: perturbedComponentsPos[compIdx].l2 = mmcComponents[compIdx].l2 + delta; break;
			case 5: perturbedComponentsPos[compIdx].l3 = mmcComponents[compIdx].l3 + delta; break;
			case 6: perturbedComponentsPos[compIdx].alpha = mmcComponents[compIdx].alpha + delta; break;
			case 7: perturbedComponentsPos[compIdx].beta = mmcComponents[compIdx].beta + delta; break;
			case 8: perturbedComponentsPos[compIdx].gamma = mmcComponents[compIdx].gamma + delta; break;
			}

			// 计算正向扰动后的局部TDF
			std::vector<std::vector<double>> perturbedPhisPos(nComp);
			for (int i = 0; i < nComp; ++i) {
				std::vector<double> componentPhi;
				std::vector<std::vector<double>> componentPhiDrv;
				calc_Phi_Node(allVertex, minL, perturbedComponentsPos, componentPhi, componentPhiDrv, i, p, epsilon);
				perturbedPhisPos[i] = componentPhi;
			}

			// 计算正向扰动后的全局TDF
			std::vector<double> perturbedGlobalPhiPos;
			std::vector<std::vector<double>> dummyPhimaxDphiPos;
			ComputeGlobalTDF(perturbedPhisPos, 100, perturbedGlobalPhiPos, dummyPhimaxDphiPos);

			// 计算正向扰动后的Heaviside值
			std::vector<double> perturbedHPos(nVertex);
			for (int i = 0; i < nVertex; ++i) {
				perturbedHPos[i] = Heaviside(perturbedGlobalPhiPos[i], alpha, epsilon);
			}

			for (int vIndex = 0; vIndex < nVertex; ++vIndex) {
				// 保存每个网格点的dH值
				dH_values[compIdx][varIdx][vIndex] = (perturbedHPos[vIndex] - originalVertexH[vIndex]) / delta;
			}
		}
	}

	std::vector<double> CeNode(nVertex, 0);
	for (int p_id = 0; p_id < nGrid; ++p_id)
	{
		for (int j = 0; j < 8; ++j)
		{
			CeNode[cell2V.at(p_id)[j].first] += Ce[p_id] / 8.0;
		}
	}

	// 计算最终的灵敏度值
	for (int compIdx = 0; compIdx < nComp; ++compIdx) {
		for (int varIdx = 0; varIdx < nVarsPerComp; ++varIdx) {
			double df0_approx = 0.0;
			double dfd_approx = 0.0;

			for (int gridIdx = 0; gridIdx < nVertex; ++gridIdx) {
				// 目标函数灵敏度（柔顺度）
				df0_approx += -CeNode[gridIdx] * dH_values[compIdx][varIdx][gridIdx];

				// 体积约束灵敏度
				dfd_approx += dH_values[compIdx][varIdx][gridIdx] * 1.0 / 8.0;
			}

			// 存储灵敏度值
			df0dx[compIdx * nVarsPerComp + varIdx] = df0_approx;
			dfdx[compIdx * nVarsPerComp + varIdx] = dfd_approx;
		}
	}

	std::ofstream outFile1("df0dx anaytic.txt", std::ios::out | std::ios::app);
	outFile1 << iter << "th " << "df0dx: " << std::endl;
	int cnt = 0;
	for (auto df0 : df0dx)
	{
		cnt++;
		outFile1 << cnt << ":" << df0 << " ";
	}
	outFile1 << std::endl;

	std::ofstream outFile2("dfdx anaytic.txt", std::ios::out | std::ios::app);
	cnt = 0;
	outFile2 << iter << "th " << "dfdx: " << std::endl;
	for (auto df : dfdx)
	{
		cnt++;
		outFile2 << cnt << ":" << df << " ";
	}
	outFile2 << std::endl;
}

void FlattenDesignVariables(
	const std::vector<MMC_Component>& components,
	std::vector<double>& vars)  // 输出：预分配的数组
{
	for (const auto& comp : components) {
		// 位置变量
		vars.emplace_back(comp.x0);
		vars.emplace_back(comp.y0);
		vars.emplace_back(comp.z0);

		// 尺寸变量
		vars.emplace_back(comp.l1);
		vars.emplace_back(comp.l2);
		vars.emplace_back(comp.l3);

		// 角度变量
		vars.emplace_back(comp.alpha);
		vars.emplace_back(comp.beta);
		vars.emplace_back(comp.gamma);
	}
}

void UpdateComponentsFromVector(
	std::vector<MMC_Component>& components,
	std::vector<double>& x)  // 输入：设计变量数组
{
	int idx = 0;
	for (auto& comp : components) {
		// 位置变量
		comp.x0 = x[idx++];
		comp.y0 = x[idx++];
		comp.z0 = x[idx++];

		// 尺寸变量
		comp.l1 = x[idx++];
		comp.l2 = x[idx++];
		comp.l3 = x[idx++];

		// 角度变量
		comp.alpha = x[idx++];
		comp.beta = x[idx++];
		comp.gamma = x[idx++];
	}
}

void Cuda_TopOpt_MMC_Node(
	const std::map<int, Eigen::Vector3d>& allVertex,
	const std::map<int, std::vector<std::pair<int, Eigen::Vector3d>>>& cell2V,
	const std::string& input_model_file,
	const std::string& bdy_file,
	const std::vector<int>& protect_cell_ids,
	const std::vector<std::vector<std::pair<int, double>>>& neighbor_info,
	double target_rate,
	double max_iter,
	double Yangs,
	double Poissons,
	std::vector<double>& density,
	std::atomic_bool& is_need_update_model,
	int mmcType) {
	Timer timer;
	Eigen::setNbThreads(24);
	Eigen::initParallel();

	std::vector<CBSplineSolid> patch;

	ToolFunc::readFileXml(input_model_file, patch);
	timer.out("read_xml");

	std::vector<CPoint> bdy_points;
	std::vector<int> bdy_type;
	std::vector<CVector> bdy_value;

	// add bdy info
	{
		int control_size = 0;
		std::ifstream in(bdy_file);
		in >> control_size;
		std::cout << "control_size = " << control_size << "\n";

		for (int i = 0; i < control_size; ++i) {
			double x, y, z;
			int cur_type;
			double vx, vy, vz;
			in >> x >> y >> z >> cur_type >> vx >> vy >> vz;
			CPoint cp(x, y, z);
			CVector cv(vx, vy, vz);
			bdy_points.push_back(cp);
			bdy_value.push_back(cv);   // 边界值
			bdy_type.push_back(cur_type); // 边界类型
		}
	}
	timer.out("add bdy info");
	gpmesh::GPLog::init(spdlog::level::trace);

	CUDA_ERROR(cudaDeviceReset());
	std::vector<boundary_point> neumann_points;
	std::vector<boundary_point> dirichlet_points;

	timer.out("deal with bdy info");
	std::vector<double> dc_pre(patch.size(), 0);
	auto exist_num = (double)patch.size();
	auto tot_num = (double)patch.size();
	//while(1){
	// data init
	gpmesh::GPSplineSolidV0 m_spline;
	//    m_spline.TopOpt_GPU_Init(patch, neumann_points, dirichlet_points);
	//    timer.out("TopOpt_GPU_Init");
	m_spline.init_patch(patch);
	timer.out("TopOpt_GPU_Geometry_Init");

	double xmax = -1e9, xmin = 1e9, ymax = -1e9, ymin = 1e9, zmax = -1e9, zmin = 1e9;
	for (auto v : allVertex) {
		xmax = std::max(xmax, v.second.x());
		xmin = std::min(xmin, v.second.x());
		ymax = std::max(ymax, v.second.y());
		ymin = std::min(ymin, v.second.y());
		zmax = std::max(zmax, v.second.z());
		zmin = std::min(zmin, v.second.z());
	}
	vector<double> xyzRange{ xmin, xmax, ymin, ymax, zmin, zmax };
	const double DL = xmax - xmin;
	const double DW = ymax - ymin;
	const double DH = zmax - zmin;
	//std::cout << xmax << " " << xmin << " " << ymax << " " << ymin << " " << zmax << " " << zmin << std::endl;


	timer.out("before sys");
	const SparseSystem& temp_sys = m_spline.m_sparse_system;
	// l1/2,l2/2,l3/2,alp,bet,gam
	// 1: z / y 2: x / z 3: y / x
	// 2.5维 参数1
	std::vector<double> BaseComponent;
	// DL,DW,DH,xInt,yInt,zInt
	std::vector<MMC_Component> mmcComponents;
	// type：1/2——2.5维
	if (mmcType == 1) {
		// 1 * 30 * 20
		BaseComponent = { 0.6, 7, 0.8, atan(1), 0, 0 };
		mmcComponents = init(xyzRange, 0.5, 5, 5, BaseComponent, mmcType);
	}
	else if (mmcType == 2) {
		// 1 * 60 * 20
		BaseComponent = { 0.6, 15, 1, atan(0.5), 0, 0 };
		mmcComponents = init(xyzRange, 0.5, 10, 5, BaseComponent, mmcType);
	}
	else if (mmcType == 3) {
		// 1 * 60 * 60
		BaseComponent = { 0.6, 15, 1, atan(1.0), 0, 0 };
		mmcComponents = init(xyzRange, 0.5, 10, 10, BaseComponent, mmcType);
	}
	else if (mmcType == 4) {
		// 2.5维 参数2
		// 2021教育论文组件参数
		// 1 * 128 * 64
		BaseComponent = { 24, 5, 0.6, 0, atan(1.0), 0 };
		mmcComponents = init(xyzRange, 16, 0.5, 16, BaseComponent);
		mmcComponents = init(xyzRange, 16, 0.5, 16, BaseComponent);
	}
	else {
		BaseComponent = { 0.6, 8, 1, atan(1.0), 0, 0 };
		mmcComponents = init(xyzRange, 0.5, 5, 5, BaseComponent, -1);
	}

	// 2.5维 参数2
	// 2021教育论文组件参数
	//std::vector<double> BaseComponent = { 24, 5, 0.6, 0, atan(1.0), 0 };
	//std::vector<MMC_Component> mmcComponents = init(xyzRange, 16, 0.5, 16, BaseComponent);


	// 3维 参数1
	/*std::vector<double> BaseComponent = { 1, 9, 1, atan(0.5), atan(0.8), atan(2.5)};
	std::vector<MMC_Component> mmcComponents = init(xyzRange, 2, 5, 2.5, BaseComponent);*/


	timer.out("after sys");
	{
		// B: 所有全局控制点  A: 边界点
		std::vector<Cuda_TopOpt_Pipeline_nsp::Point> A, B;
		std::vector<int> global_ids;

		for (auto& p : bdy_points) {
			A.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ p.getX(), p.getY(), p.getZ() });
		}
		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				auto& pp = patch[p_id].coefs()[local_id];
				B.emplace_back(Cuda_TopOpt_Pipeline_nsp::Point{ pp.getX(), pp.getY(), pp.getZ() });
			}
		}

		// 找到距离A[i]最近的点并返回索引
		GetNearestIds(A, B, global_ids);
		timer.out("after GetNearestIds");
		std::cout << "bdy_points.size() = " << bdy_points.size() << std::endl;
		auto& global_local_map = temp_sys.get_l2g_map();

		for (int i = 0; i < bdy_points.size(); ++i) {
			auto& p = bdy_points[i];
			int patch_id = global_ids[i] / 64;
			int local_id = global_ids[i] % 64;
			int global_id = global_local_map[patch_id][local_id]; // 控制点全局索引

			if (global_id == -1) {
				std::cerr << "cant find control points\n";
				//std::cout << cc.getX() <<" " <<cc.getY() <<" " << cc.getZ() <<"\n";
				std::cout << p.getX() << " " << p.getY() << " " << p.getZ() << "\n";
			}
			else {
				if (bdy_type[i] == 2) {
					neumann_points.emplace_back(condition_type::neumann, bdy_value[i], global_id);
				}
				else if (bdy_type[i] == 1) {
					dirichlet_points.emplace_back(condition_type::dirichlet, CVector(0, 0, 0), global_id);
				}
			}
		}
	}

	m_spline.init_dirichlet_points(dirichlet_points);

	std::vector<Eigen::MatrixXd> local_matrix(patch.size(), Eigen::MatrixXd::Zero(192, 192)); // 64 * 3
	// first assembly
	m_spline.Get_Local_Matrix(local_matrix, Yangs, Poissons);
	timer.out("Get_Local_Matrix");
	m_spline.m_sparse_system.rhs().setZero();

	m_spline.helper_set_neumann_point_in_sparse_system_rhs(neumann_points);
	m_spline.helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);

	int topopt_iter = 0;

	std::vector<int> isPicked(patch.size(), 0);
	for (auto& xx : protect_cell_ids) {
		isPicked[xx] = 1;
	}

	// mma参数 
	const int nComp = mmcComponents.size(); // 组件数量

	const double minL = 0.;
	// 单个组件的边界定义
	std::vector<double> single_xmin = {
		xmin, ymin, zmin,  // x0,y0,z0 下限
		minL, minL, minL,              // l1,l2,l3 下限
		-M_PI, -M_PI, -M_PI         // alp,bet,gam 下限
	};

	std::vector<double> single_xmax = {
		xmax, ymax, zmax,     // x0,y0,z0 上限
		sqrt(DL * DL + DW * DW + DH * DH) / 2.0, // l1 上限
		sqrt(DL * DL + DW * DW + DH * DH) / 2.0, // l2 上限
		sqrt(DL * DL + DW * DW + DH * DH) / 2.0, // l3 上限
		M_PI, M_PI, M_PI            // alp,bet,gam 上限
	};


	// 扩展到所有组件
	std::vector<double> xmin_all(nComp * 9);
	std::vector<double> xmax_all(nComp * 9);

	for (int i = 0; i < nComp; ++i) {
		std::copy(single_xmin.begin(), single_xmin.end(), xmin_all.begin() + i * 9);
		std::copy(single_xmax.begin(), single_xmax.end(), xmax_all.begin() + i * 9);
	}

	std::vector<std::vector<double>> allPhis;   // 所有组件的TDF nVertex * nComponent
	std::vector<std::vector<std::vector<double>>> allPhiDrv; // 所有组件的TDF对设计变量的导数 nVertex * nComponent * 9
	std::vector<double> globalPhi; // 全局TDF nGrid * 1
	const double alhpa = 1e-3, epsilon = 0.20; // Heaviside 函数参数
	MMASolver mmaSolver(9 * mmcComponents.size(), 1); // 1个约束（体积）
	const int nVertex = allVertex.size();
	const int nGrid = density.size();
	const int nVarsPerComp = 9;
	allPhiDrv.resize(nComp);
	allPhis.resize(nComp);
	globalPhi.resize(nVertex);
	while (topopt_iter < max_iter) {
		std::cout << topopt_iter++ << std::endl;
		int reduce_num = (int)(0.02 * (exist_num - tot_num * target_rate));
		if (reduce_num < 1) break;
		Eigen::MatrixXd U;

		// 计算所有杆件局部TDF
		for (int i = 0; i < mmcComponents.size(); ++i)
		{
			std::vector<double> componentPhi;
			std::vector<std::vector<double>> componentPhiDrv;
			calc_Phi_Node(allVertex, minL, mmcComponents, componentPhi, componentPhiDrv, i, p, EPS);
			allPhiDrv[i] = componentPhiDrv;
			allPhis[i] = componentPhi;
		}

		// 计算全局TDF
		const int lmd = 100;
		std::vector<std::vector<double>> phimaxDphi; // 全局TDF对局部TDF的灵敏度  nVertex * nComponent
		{
			globalPhi.assign(nVertex, 0.0);
			phimaxDphi.assign(nVertex, std::vector<double>(nComp, 0.0));
			const double tiny = 1e-9;

			for (int g = 0; g < nVertex; ++g) {
				std::vector<double> expValues(nComp);
				double sumExp = 0.0;
				for (int c = 0; c < nComp; ++c) {
					expValues[c] = std::exp(lmd * allPhis[c][g]);
					sumExp += expValues[c];
				}
				globalPhi[g] = std::max(-1e3, std::log(sumExp) / lmd);   // 公式4  K-S函数

				for (int c = 0; c < nComp; ++c) {
					phimaxDphi[g][c] = expValues[c] / (sumExp + tiny);   // 公式10  
				}
			}
		}

		std::vector<double> vertexH;
		vertexH.resize(nVertex);
		// 计算heaviside值
		for (int i = 0; i < globalPhi.size(); ++i)
		{
			vertexH[i] = Heaviside(globalPhi[i], alhpa, epsilon);
		}

		// 节点分配密度给单元
		for (int i = 0; i < nGrid; ++i)
		{
			auto pairVec = cell2V.at(i);
			density[i] = 0;
			for (int j = 0; j < 8; ++j)
			{
				int tag = pairVec[j].first;
				density[i] += vertexH[tag];
			}
			density[i] /= 8.0;
		}

		m_spline.Solver_With_Density_MMC(neumann_points, dirichlet_points, density, U);

		Eigen::MatrixXd UU(patch.size(), 64 * 3);

		for (int p_id = 0; p_id < patch.size(); ++p_id) {
			for (int local_id = 0; local_id < patch[p_id].coefs().size(); ++local_id) {
				int global_id = m_spline.m_sparse_system.get_l2g_map()[p_id][local_id];
				UU(p_id, local_id) = U(global_id, 0);
				UU(p_id, local_id + 64) = U(global_id + m_spline.m_sparse_system.getControlNums(), 0);
				UU(p_id, local_id + 128) = U(global_id + 2 * m_spline.m_sparse_system.getControlNums(), 0);
			}
		}

		std::vector<double> df0dx; // 目标函数灵敏度
		std::vector<double> dfdx; // 体积约束灵敏度

		double vol = 0;
		for (int i = 0; i < density.size(); ++i)
		{
			vol += density[i];
		}
		vol = vol / density.size();
		std::vector<double> delta(1, vol - target_rate);

		// 计算灵敏度
		{
			// 公式 10 Heaviside函数的导数
			std::vector<double> delta_H(nVertex);
			for (int i = 0; i < nVertex; ++i) {
				if (std::abs(globalPhi[i]) > epsilon) {
					delta_H[i] = 0.0;
				}
				else {
					delta_H[i] = 3.0 * (1.0 - alhpa) / (4.0 * epsilon) *
						(1.0 - globalPhi[i] * globalPhi[i] / (epsilon * epsilon));
				}
			}

			double allCe = 0.;
			// 单元应变能Ce
			std::vector<double> Ce(nGrid);
			for (int p_id = 0; p_id < nGrid; ++p_id) {
				Ce[p_id] = UU.row(p_id) * local_matrix[p_id] * UU.row(p_id).transpose();
				allCe += Ce[p_id];
			}

			std::cout << "iter: " << topopt_iter << " C: " << allCe << std::endl;

			// 节点应变能CeNode
			std::vector<double> CeNode(nVertex, 0);
			for (int p_id = 0; p_id < nGrid; ++p_id)
			{
				for (int j = 0; j < 8; ++j)
				{
					CeNode[cell2V.at(p_id)[j].first] += Ce[p_id] / 8.0;
				}
			}


			// 定义 PhimaxDrvAct，维度为 nComponent * nVertex * 9
			std::vector<std::vector<std::vector<double>>> PhimaxDrvAct(nComp, std::vector<std::vector<double>>(nVertex, std::vector<double>(9, 0.0)));

			// 遍历每个组件
			for (int c = 0; c < nComp; ++c) {
				// 遍历每个顶点
				for (int v = 0; v < nVertex; ++v) {
					// 遍历每个设计变量（9个）
					for (int d = 0; d < 9; ++d) {
						// 计算 PhimaxDrvAct[c][v][d]
						PhimaxDrvAct[c][v][d] = phimaxDphi[v][c] * allPhiDrv[c][v][d];
					}
				}
			}

			df0dx.assign(nComp * nVarsPerComp, 0.0); // 调整为实际设计变量数
			dfdx.assign(nComp * nVarsPerComp, 0.0);

			std::vector<double> volNod(nVertex, 0.);
			for (auto cellPair : cell2V)
			{
				for (int i = 0; i < 8; ++i)
				{
					int tag = cellPair.second[i].first;
					volNod[tag] += 1.0 / 8.0;
				}
			}

			// 计算df0dx （目标函数灵敏度）
			for (int i = 0; i < nComp; ++i)
			{
				for (int j = 0; j < nVarsPerComp; ++j)
				{
					for (int k = 0; k < nVertex; ++k)
					{
						df0dx[i * nVarsPerComp + j] += -CeNode[k] * delta_H[k] * PhimaxDrvAct[i][k][j];
					}
				}
			}

			// 计算dfdx (体积灵敏度)
			for (int i = 0; i < nComp; ++i) {
				for (int j = 0; j < nVarsPerComp; ++j) {
					for (int k = 0; k < nVertex; ++k) {
						dfdx[i * nVarsPerComp + j] += volNod[k] * delta_H[k] * PhimaxDrvAct[i][k][j] * 1.0 / (DL * DW * DH); // 公式9
					}
				}
			}
		}

		vector<double> designVariables;
		FlattenDesignVariables(mmcComponents, designVariables);
		mmaSolver.Update(designVariables.data(), df0dx.data(), delta.data(), dfdx.data(), xmin_all.data(), xmax_all.data());
		UpdateComponentsFromVector(mmcComponents, designVariables);
		is_need_update_model = true;
	}
}
