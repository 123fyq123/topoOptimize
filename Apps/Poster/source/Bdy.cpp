#include "Bdy.h"
#include <algorithm>
#include "Subdivision.h"
#include "FileHandler.h"
#include "Timer.h"

namespace IGA {

	void BdyConvert(const std::string& filename, HexBdy& hexbdy, CCBezierBdy& ccbdy) {
		// 输入mesh信息  辅助转换
		HexV3fMesh* mesh = new HexV3fMesh();
		std::cout << filename << std::endl;
		FileHandler file(filename);
		file.read2HexMeshObj(mesh, true, true, true);
		// 得到样条信息 m_solids
		ApproximationCC AppCC(mesh, nullptr);
		AppCC.construct_bss_new();
		BSplineSolid m_solids = AppCC.bss;
		// 已经插入的控制点 不需要再插入
		std::set<V3f> have_insert;         
		// 遍历每一个六面体边界条件
		
		for (auto& hbdy : hexbdy) {
			// 求出所在元素对应的控制顶点
			std::set<V3f> select_control_points;
			//std::cout << "hbdy.element_type = " << hbdy.element_type << std::endl;
			// 分类别分别计算影响的控制顶点
			if (hbdy.element_type == HexBdyBase::ElementType::Point) {
				auto pos = hbdy.element_points[0];	
				VertexHandle vh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->vertices_.size(); ++i) {
					double dx = (pos[0] - (mesh->vertices_[i]).x);
					double dy = (pos[1] - (mesh->vertices_[i]).y);
					double dz = (pos[2] - (mesh->vertices_[i]).z);

					double dist = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (dist < min_point_dis) {
						min_point_dis = dist;
						vh = (VertexHandle)i;
					}
				}
				auto control_point = AppCC.cal_corner_point(vh);
				if (have_insert.count(control_point) == 0) {
					select_control_points.insert(control_point);
					have_insert.insert(control_point);
				}
				
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Line) {
				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);	
				}

				sort(pos.begin(), pos.end());
				EdgeHandle eh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->edges_.size(); ++i) {
					std::vector<std::vector<double>>  e_pos;
					auto vh1 = mesh->edges_[i].from_h;
					auto vh2 = mesh->edges_[i].to_h;

					e_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh1].x,mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });
					e_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });
					

					sort(e_pos.begin(), e_pos.end());

					double dist = 0;

					for (int j = 0; j < e_pos.size(); ++j) {
						for (int kk = 0; kk < 3; ++kk) {
							double dd = (e_pos[j][kk] - pos[j][kk]);
							dist += dd * dd;
						}
						
					}
					dist = std::sqrt(dist);

					if (dist < min_point_dis) {
						min_point_dis = dist;
						eh = (EdgeHandle)i;
					}
				}
				auto vh1 = mesh->edges_[eh].from_h;
				auto vh2 = mesh->edges_[eh].to_h;

				std::vector<V3f> control_points = std::vector<V3f>{ AppCC.cal_edge_point(eh, vh1), AppCC.cal_edge_point(eh, vh2), AppCC.cal_corner_point(vh1), AppCC.cal_corner_point(vh2) };

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Face) {
				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);
				}
				
				sort(pos.begin(), pos.end());

				FaceHandle fh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->faces_.size(); ++i) {
					std::vector<std::vector<double>> f_pos;
					auto vh1 = mesh->faces_[i].vs_f[0];
					auto vh2 = mesh->faces_[i].vs_f[1];
					auto vh3 = mesh->faces_[i].vs_f[2];
					auto vh4 = mesh->faces_[i].vs_f[3];

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh1].x, mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh3].x, mesh->vertices_[vh3].y, mesh->vertices_[vh3].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh4].x, mesh->vertices_[vh4].y, mesh->vertices_[vh4].z });

					sort(f_pos.begin(), f_pos.end());

					double dist = 0;

					for (int j = 0; j < f_pos.size(); ++j) {
						for (int kk = 0; kk < 3; ++kk) {
							double dd = (f_pos[j][kk] - pos[j][kk]);
							dist += dd * dd;
						}

					}
					dist = std::sqrt(dist);

					if (dist < min_point_dis) {
						min_point_dis = dist;
						fh = (FaceHandle)i;
					}
				}
				//std::cout << "匹配对应误差 = " << min_point_dis << std::endl;
				std::vector<V3f> control_points;
				for (int v_num = 0; v_num < 4; v_num++)
				{
					control_points.push_back(AppCC.solid_face_point(fh, mesh->faces_[fh].vs_f[v_num]));
					control_points.push_back(AppCC.solid_corner_point(mesh->faces_[fh].vs_f[v_num]));
				}
				for (int e_num = 0; e_num < 4; e_num++)
				{
					EdgeHandle e = mesh->faces_[fh].edges_[e_num];
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].from_h));
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].to_h));	
				}

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Cell) {
				// 暂时没用， 有需要再加吧

			}


			//std::cout <<"select_control_points.size() = "  << select_control_points.size() << std::endl;
			// 对于每个控制顶点 加一个ccBezier的边界条件
			for (auto& control_point : select_control_points) {

				CCBezierBdyBase ccbase;
				ccbase.bdy_value = hbdy.bdy_value;
				ccbase.bdy_type = hbdy.bdy_type;


				int mindis_id = -1; double min_dis = std::numeric_limits<float>::max(); int mindis_cellid = -1;
				for (int i = 0; i < m_solids.Local_V3f.size(); i++)
				{
					for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
					{
						V3f temp_v = m_solids.Local_V3f[i][j];
						double dis = distance(temp_v, control_point);
						if (dis < min_dis)
						{
							min_dis = dis;
							mindis_id = j;
							mindis_cellid = i;
						}
					}
				}
				//std::cout << "控制点对应误差 = " << min_dis << std::endl;
				ccbase.cell_id = mindis_cellid;
				ccbase.local_base_id = mindis_id;
				ccbase.global_base_id = m_solids.Local_to_Global[mindis_cellid][mindis_id];
				ccbdy.push_back(ccbase);
			}
		}

	}

	void BdyConvert_V2(const std::string& filename, HexBdy& hexbdy, CCBezierBdy& ccbdy, ApproximationCC& AppCC) {
		// 输入mesh信息  辅助转换
		HexV3fMesh* mesh = new HexV3fMesh();
		std::cout << filename << std::endl;
		FileHandler file(filename);
		file.read2HexMeshObj(mesh, true, true, true);
		// 得到样条信息 m_solids
		//ApproximationCC AppCC(mesh, nullptr);
		//AppCC.construct_bss_new();
		BSplineSolid m_solids = AppCC.bss;
		// 已经插入的控制点 不需要再插入
		std::set<V3f> have_insert;
		// 遍历每一个六面体边界条件

		for (auto& hbdy : hexbdy) {
			// 求出所在元素对应的控制顶点
			std::set<V3f> select_control_points;
			//std::cout << "hbdy.element_type = " << hbdy.element_type << std::endl;
			// 分类别分别计算影响的控制顶点
			if (hbdy.element_type == HexBdyBase::ElementType::Point) {
				auto pos = hbdy.element_points[0];
				VertexHandle vh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->vertices_.size(); ++i) {
					double dx = (pos[0] - (mesh->vertices_[i]).x);
					double dy = (pos[1] - (mesh->vertices_[i]).y);
					double dz = (pos[2] - (mesh->vertices_[i]).z);

					double dist = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (dist < min_point_dis) {
						min_point_dis = dist;
						vh = (VertexHandle)i;
					}
				}
				auto control_point = AppCC.cal_corner_point(vh);
				if (have_insert.count(control_point) == 0) {
					select_control_points.insert(control_point);
					have_insert.insert(control_point);
				}

			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Line) {
				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);
				}

				sort(pos.begin(), pos.end());
				EdgeHandle eh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->edges_.size(); ++i) {
					std::vector<std::vector<double>>  e_pos;
					auto vh1 = mesh->edges_[i].from_h;
					auto vh2 = mesh->edges_[i].to_h;

					e_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh1].x, mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });
					e_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });


					sort(e_pos.begin(), e_pos.end());

					double dist = 0;

					for (int j = 0; j < e_pos.size(); ++j) {
						for (int kk = 0; kk < 3; ++kk) {
							double dd = (e_pos[j][kk] - pos[j][kk]);
							dist += dd * dd;
						}

					}
					dist = std::sqrt(dist);

					if (dist < min_point_dis) {
						min_point_dis = dist;
						eh = (EdgeHandle)i;
					}
				}
				auto vh1 = mesh->edges_[eh].from_h;
				auto vh2 = mesh->edges_[eh].to_h;

				std::vector<V3f> control_points = std::vector<V3f>{ AppCC.cal_edge_point(eh, vh1), AppCC.cal_edge_point(eh, vh2), AppCC.cal_corner_point(vh1), AppCC.cal_corner_point(vh2) };

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Face) {
				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);
				}

				sort(pos.begin(), pos.end());

				FaceHandle fh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->faces_.size(); ++i) {
					std::vector<std::vector<double>> f_pos;
					auto vh1 = mesh->faces_[i].vs_f[0];
					auto vh2 = mesh->faces_[i].vs_f[1];
					auto vh3 = mesh->faces_[i].vs_f[2];
					auto vh4 = mesh->faces_[i].vs_f[3];

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh1].x, mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh3].x, mesh->vertices_[vh3].y, mesh->vertices_[vh3].z });

					f_pos.emplace_back(std::vector<double>{ mesh->vertices_[vh4].x, mesh->vertices_[vh4].y, mesh->vertices_[vh4].z });

					sort(f_pos.begin(), f_pos.end());

					double dist = 0;

					for (int j = 0; j < f_pos.size(); ++j) {
						for (int kk = 0; kk < 3; ++kk) {
							double dd = (f_pos[j][kk] - pos[j][kk]);
							dist += dd * dd;
						}

					}
					dist = std::sqrt(dist);

					if (dist < min_point_dis) {
						min_point_dis = dist;
						fh = (FaceHandle)i;
					}
				}
				//std::cout << "匹配对应误差 = " << min_point_dis << std::endl;
				std::vector<V3f> control_points;
				for (int v_num = 0; v_num < 4; v_num++)
				{
					control_points.push_back(AppCC.solid_face_point(fh, mesh->faces_[fh].vs_f[v_num]));
					control_points.push_back(AppCC.solid_corner_point(mesh->faces_[fh].vs_f[v_num]));
				}
				for (int e_num = 0; e_num < 4; e_num++)
				{
					EdgeHandle e = mesh->faces_[fh].edges_[e_num];
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].from_h));
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].to_h));
				}

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Cell) {
				// 暂时没用， 有需要再加吧

			}


			//std::cout <<"select_control_points.size() = "  << select_control_points.size() << std::endl;
			// 对于每个控制顶点 加一个ccBezier的边界条件
			for (auto& control_point : select_control_points) {

				CCBezierBdyBase ccbase;
				ccbase.bdy_value = hbdy.bdy_value;
				ccbase.bdy_type = hbdy.bdy_type;


				int mindis_id = -1; double min_dis = std::numeric_limits<float>::max(); int mindis_cellid = -1;
				for (int i = 0; i < m_solids.Local_V3f.size(); i++)
				{
					for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
					{
						V3f temp_v = m_solids.Local_V3f[i][j];
						double dis = distance(temp_v, control_point);
						if (dis < min_dis)
						{
							min_dis = dis;
							mindis_id = j;
							mindis_cellid = i;
						}
					}
				}
				//std::cout << "控制点对应误差 = " << min_dis << std::endl;
				ccbase.cell_id = mindis_cellid;
				ccbase.local_base_id = mindis_id;
				ccbase.global_base_id = m_solids.Local_to_Global[mindis_cellid][mindis_id];
				ccbdy.push_back(ccbase);
			}
		}

	}

	void BdyConvert_V3(const HexV3fMesh* mesh, HexBdy& hexbdy, CCBezierBdy& ccbdy, ApproximationCC& AppCC) {
		BSplineSolid m_solids = AppCC.bss;
		// 已经插入的控制点 不需要再插入
		std::set<V3f> have_insert;

		// 遍历每一个六面体边界条件
		std::vector<HexBdyBase> Point_hbdy;
		std::vector<HexBdyBase> Line_hbdy;
		std::vector<HexBdyBase> Face_hbdy;
		std::vector<HexBdyBase> Cell_hbdy;
		for (auto& hbdy : hexbdy) {
			if (hbdy.element_type == HexBdyBase::ElementType::Point) {
				Point_hbdy.push_back(hbdy);
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Line) {
				Line_hbdy.push_back(hbdy);
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Face) {
				Face_hbdy.push_back(hbdy);
			}
			else if (hbdy.element_type == HexBdyBase::ElementType::Cell) {
				Cell_hbdy.push_back(hbdy);
			}
		}

		if (!Point_hbdy.empty()) {
			for (auto& hbdy : Point_hbdy) {
				std::set<V3f> select_control_points;
				auto pos = hbdy.element_points[0];
				VertexHandle vh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->vertices_.size(); ++i) {
					double dx = (pos[0] - (mesh->vertices_[i]).x);
					double dy = (pos[1] - (mesh->vertices_[i]).y);
					double dz = (pos[2] - (mesh->vertices_[i]).z);

					double dist = std::sqrt((dx * dx) + (dy * dy) + (dz * dz));
					if (dist < min_point_dis) {
						min_point_dis = dist;
						vh = (VertexHandle)i;
					}
				}
				auto control_point = AppCC.cal_corner_point(vh);
				if (have_insert.count(control_point) == 0) {
					select_control_points.insert(control_point);
					have_insert.insert(control_point);
				}

				for (auto& control_point : select_control_points) {

					CCBezierBdyBase ccbase;
					ccbase.bdy_value = hbdy.bdy_value;
					ccbase.bdy_type = hbdy.bdy_type;


					int mindis_id = -1; double min_dis = std::numeric_limits<float>::max(); int mindis_cellid = -1;
					for (int i = 0; i < m_solids.Local_V3f.size(); i++)
					{
						for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
						{
							V3f temp_v = m_solids.Local_V3f[i][j];
							double dis = distance(temp_v, control_point);
							if (dis < min_dis)
							{
								min_dis = dis;
								mindis_id = j;
								mindis_cellid = i;
							}
						}
					}
					//std::cout << "控制点对应误差 = " << min_dis << std::endl;
					ccbase.cell_id = mindis_cellid;
					ccbase.local_base_id = mindis_id;
					ccbase.global_base_id = m_solids.Local_to_Global[mindis_cellid][mindis_id];
					ccbdy.push_back(ccbase);
				}
			}
		}

		if (!Line_hbdy.empty()) {
			std::vector<std::vector<std::vector<double>>>  e_pos(mesh->edges_.size());
			for (int i = 0; i < mesh->edges_.size(); ++i) {
				auto vh1 = mesh->edges_[i].from_h;
				auto vh2 = mesh->edges_[i].to_h;

				e_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh1].x, mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });
				e_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });

				sort(e_pos[i].begin(), e_pos[i].end());
			}

			for (auto& hbdy : Line_hbdy) {
				std::set<V3f> select_control_points;

				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);
				}

				sort(pos.begin(), pos.end());
				EdgeHandle eh;
				double min_point_dis = std::numeric_limits<float>::max();
				for (int i = 0; i < mesh->edges_.size(); ++i) {
					double dist = 0;

					for (int j = 0; j < e_pos.size(); ++j) {
						for (int kk = 0; kk < 3; ++kk) {
							double dd = (e_pos[i][j][kk] - pos[j][kk]);
							dist += dd * dd;
						}

					}
					//dist = std::sqrt(dist);
					if (dist < min_point_dis) {
						min_point_dis = dist;
						eh = (EdgeHandle)i;
					}
				}
				auto vh1 = mesh->edges_[eh].from_h;
				auto vh2 = mesh->edges_[eh].to_h;

				std::vector<V3f> control_points = std::vector<V3f>{ AppCC.cal_edge_point(eh, vh1), AppCC.cal_edge_point(eh, vh2), AppCC.cal_corner_point(vh1), AppCC.cal_corner_point(vh2) };

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}

				for (auto& cp : select_control_points) {

					CCBezierBdyBase ccbase;
					ccbase.bdy_value = hbdy.bdy_value;
					ccbase.bdy_type = hbdy.bdy_type;

					int mindis_id = -1; double min_dis = std::numeric_limits<float>::max(); int mindis_cellid = -1;
					for (int i = 0; i < m_solids.Local_V3f.size(); i++)
					{
						for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
						{
							V3f temp_v = m_solids.Local_V3f[i][j];
							double dis = distance(temp_v, cp);
							if (dis < min_dis)
							{
								min_dis = dis;
								mindis_id = j;
								mindis_cellid = i;
							}
						}
					}
					//std::cout << "控制点对应误差 = " << min_dis << std::endl;
					ccbase.cell_id = mindis_cellid;
					ccbase.local_base_id = mindis_id;
					ccbase.global_base_id = m_solids.Local_to_Global[mindis_cellid][mindis_id];
					ccbdy.push_back(ccbase);
				}
			}
		}

		if (!Face_hbdy.empty()) {
			std::vector < std::vector<std::vector<double>>> f_pos(mesh->faces_.size());
			std::vector<int> valid_i;
			for (int i = 0; i < mesh->faces_.size(); ++i) {
				if (mesh->neighbor_f[i].size() >= 2) {
					continue;
				}
				valid_i.emplace_back(i);
				auto vh1 = mesh->faces_[i].vs_f[0];
				auto vh2 = mesh->faces_[i].vs_f[1];
				auto vh3 = mesh->faces_[i].vs_f[2];
				auto vh4 = mesh->faces_[i].vs_f[3];

				f_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh1].x, mesh->vertices_[vh1].y, mesh->vertices_[vh1].z });
				f_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh2].x, mesh->vertices_[vh2].y, mesh->vertices_[vh2].z });
				f_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh3].x, mesh->vertices_[vh3].y, mesh->vertices_[vh3].z });
				f_pos[i].emplace_back(std::vector<double>{ mesh->vertices_[vh4].x, mesh->vertices_[vh4].y, mesh->vertices_[vh4].z });

				sort(f_pos[i].begin(), f_pos[i].end());
			}

			for (auto& hbdy : Face_hbdy) {
				std::set<V3f> select_control_points;

				std::vector<std::vector<double>> pos;
				for (auto& p : hbdy.element_points) {
					pos.push_back(p);
				}

				sort(pos.begin(), pos.end());

				FaceHandle fh;
				bool find_i = false;
				for (const int& i : valid_i) {

					find_i = true;
					for (int j = 0; j < f_pos[i].size(); ++j) {
						if (f_pos[i][j][0] != pos[j][0] || f_pos[i][j][1] != pos[j][1] || f_pos[i][j][2] != pos[j][2]) {
							find_i = false;
							break;
						}
					}
					if (find_i) {
						fh = (FaceHandle)i;
						break;
					}
				}
				std::vector<V3f> control_points;
				for (int v_num = 0; v_num < 4; v_num++)
				{
					control_points.push_back(AppCC.solid_face_point(fh, mesh->faces_[fh].vs_f[v_num]));
					control_points.push_back(AppCC.solid_corner_point(mesh->faces_[fh].vs_f[v_num]));
				}
				for (int e_num = 0; e_num < 4; e_num++)
				{
					EdgeHandle e = mesh->faces_[fh].edges_[e_num];
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].from_h));
					control_points.push_back(AppCC.solid_edge_point(e, mesh->edges_[e].to_h));
				}

				for (auto& cp : control_points) {
					if (have_insert.count(cp) == 0) {
						select_control_points.insert(cp);
						have_insert.insert(cp);
					}
				}

				for (auto& cp : select_control_points) {

					CCBezierBdyBase ccbase;
					ccbase.bdy_value = hbdy.bdy_value;
					ccbase.bdy_type = hbdy.bdy_type;

					bool find_id = false;
					int mindis_id = -1;
					int mindis_cellid = -1;
					for (int i = 0; i < m_solids.Local_V3f.size(); i++)
					{
						for (int j = 0; j < m_solids.Local_V3f[i].size(); j++)
						{
							V3f temp_v = m_solids.Local_V3f[i][j];
							if (cp == temp_v) {
								mindis_id = j;
								mindis_cellid = i;
								find_id = true;
								break;
							}
						}
						if (find_id) {
							break;
						}
					}
					ccbase.cell_id = mindis_cellid;
					ccbase.local_base_id = mindis_id;
					ccbase.global_base_id = m_solids.Local_to_Global[mindis_cellid][mindis_id];
					ccbdy.push_back(ccbase);
				}
			}
		}

		if (!Cell_hbdy.empty()) {

		}


	}

	void WriteBdy(const std::string& filename, HexBdy& hexbdy, std::vector<int>& pro)
	{
		///保存所有点选产生的边界信息
		///名字一直固定？

		std::ofstream file(filename);
		double size = hexbdy.size();
		file << size << "\n";
		for (auto& thisBdy : hexbdy)
		{
			file << thisBdy.element_type << "\n";
			for (auto& vec2 : thisBdy.element_points)
			{
				file << vec2[0] << " " << vec2[1] << " " << vec2[2] << " " << "\n";
			}

			file << thisBdy.bdy_type << "\n";
			for (int i = 0; i < 3; i++) {
				file << thisBdy.bdy_value[i] << " ";
			}
			file << "\n";
		}
        std::cout <<pro.size() << "\n";
		file << pro.size() << "\n";
		for (int i = 0; i < pro.size(); i++)
		{
			file << pro[i] << " ";
		}
		file << "\n";


	}

	void ReadBdy(const std::string& filename, HexBdy& hexbdy, std::vector<int>& pro)
	{
		std::ifstream file(filename, std::ios::in);
		std::string line;


		std::getline(file, line);
		std::stringstream ss(line);
		int num;
		ss >> num;

		for (int i = 0; i < num; ++i) {
			std::getline(file, line);
			ss.clear();
			ss.str(line);
			HexBdyBase bdy;
			// element type
			int cnt[5] = { 0, 1, 2, 4, 8 };

			bdy.element_type = (IGA::HexBdyBase::ElementType)(line[0] - '0');

			for (int i = 0; i < cnt[bdy.element_type]; ++i) {
				std::getline(file, line);
				ss.clear();
				ss.str(line);
				double x, y, z;
				ss >> x >> y >> z;
				bdy.element_points.push_back({ x, y, z });
				//std::cout << x << " " << y << " " << z << std::endl;
			}
			std::getline(file, line);
			ss.clear();
			ss.str(line);
			bdy.bdy_type = (IGA::BdyBase::BdyType)(line[0] - '0');

			std::getline(file, line);
			ss.clear();
			ss.str(line);
			double x, y, z;
			ss >> x >> y >> z;
			bdy.bdy_value = std::vector<double>{ x, y, z };

			hexbdy.push_back(bdy);
		}

		std::getline(file, line);
		ss.clear();
		int protect_cnt = stoi(line);

		std::getline(file, line);
		ss.clear();
		ss.str(line);

		pro.clear();
		for (int i = 0; i < protect_cnt; ++i) {
			int cell_id; ss >> cell_id;
			pro.push_back(cell_id);
		}
	}

}