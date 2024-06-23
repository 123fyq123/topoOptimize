#include"Subdivision.h"
#include "BezierSolid.h"
#include "FileHandler.h"
#include <string>

void ApproximationCC::construct_bss1_new()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<bool> check_vertex(mesh->vertices_.size(), false);
	std::vector<bool> check_edge(mesh->edges_.size(), false);
	std::vector<bool> check_face(mesh->faces_.size(), false);
	std::map<V3f, int> quick_v3f;
	int global_v = 0;
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		if (i % 1000 == 0)
			std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.local_bdy[i].resize(27, false);
		bss.Local_V3f[i].resize(27);
		bss.Local_to_Global[i].resize(27);
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = mesh->vertices_[mesh->cells_[i].vertex_[vn]];
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].vertex_[vn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_vertex[mesh->cells_[i].vertex_[vn]] == false)
			{
				quick_v3f[v] = global_v++;
				bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
				check_vertex[mesh->cells_[i].vertex_[vn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = quick_v3f[v];
		}
		for (int en = 0; en < 12; en++)
		{
			V3f v = cal_edge_midv(mesh->cells_[i].edges_[en]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].edges_[en]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_edge[mesh->cells_[i].edges_[en]] == false)
			{
				quick_v3f[v] = global_v++;
				bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
				check_edge[mesh->cells_[i].edges_[en]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = quick_v3f[v];
		}
		for (int fn = 0; fn < 6; fn++)
		{
			V3f v = cal_face_midv(mesh->cells_[i].faces_[fn]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].faces_[fn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_face[mesh->cells_[i].faces_[fn]] == false)
			{
				quick_v3f[v] = global_v++;
				bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
				check_face[mesh->cells_[i].faces_[fn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = quick_v3f[v];
		}
		V3f cell_midv = cal_cell_midv(CellHandle(i));
		bss.Local_V3f[i][tempsize++] = cell_midv;
		bss.Local_to_Global[i][tempsize - 1] = global_v++;

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<int>>> sc(3);
		std::vector<std::vector<int>> sf(3);
		std::vector<int> sl(3);
		//第一层
		sl[0] = 5; sl[1] = 9; sl[2] = 1;
		sf[0] = sl;
		sl[0] = 17; sl[1] = 22; sl[2] = 16;
		sf[1] = sl;
		sl[0] = 6; sl[1] = 15; sl[2] = 2;
		sf[2] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = 10; sl[1] = 20; sl[2] = 8;
		sf[0] = sl;
		sl[0] = 24; sl[1] = 26; sl[2] = 25;
		sf[1] = sl;
		sl[0] = 14; sl[1] = 21; sl[2] = 12;
		sf[2] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = 4; sl[1] = 11; sl[2] = 0;
		sf[0] = sl;
		sl[0] = 19; sl[1] = 23; sl[2] = 18;
		sf[1] = sl;
		sl[0] = 7; sl[1] = 13; sl[2] = 3;
		sf[2] = sl;
		sc[2] = sf;
		bss.Control_Point[i] = sc;

		//构建局部的每个小单元的u,v,w三个方向的次数
		std::vector<int> order(3); order[0] = 2; order[1] = 2; order[2] = 2;
		bss.Local_Order[i] = order;

		//构建局部的每个小单元的u,v,w三个方向的节点向量
		std::vector<std::vector<int>> node(3);
		node[0].push_back(0); node[0].push_back(0); node[0].push_back(0);
		node[0].push_back(1); node[0].push_back(1); node[0].push_back(1);
		node[1].push_back(0); node[1].push_back(0); node[1].push_back(0);
		node[1].push_back(1); node[1].push_back(1); node[1].push_back(1);
		node[2].push_back(0); node[2].push_back(0); node[2].push_back(0);
		node[2].push_back(1); node[2].push_back(1); node[2].push_back(1);
		bss.Local_Node[i] = node;
	}
	int bdy_sum = 0;
	for (int i = 0; i < mesh->vertices_.size(); i++)
	{
		if (mesh->is_bdy(VertexHandle(i)))
			bdy_sum++;
	}
	for (int i = 0; i < mesh->edges_.size(); i++)
	{
		if (mesh->is_bdy(EdgeHandle(i)))
			bdy_sum++;
	}
	for (int i = 0; i < mesh->faces_.size(); i++)
	{
		if (mesh->is_bdy(FaceHandle(i)))
			bdy_sum++;
	}
	bss.Global_size = global_v;
	bss.Global_size_without_bdy = bss.Global_size - bdy_sum;
}

void ApproximationCC::construct_bss_new()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<bool> check_vertex(mesh->vertices_.size(), false);
	std::vector<int> check_edge(mesh->edges_.size(), 0);
	//std::vector<bool> check_vertex_and_edge(mesh->vertices_.size(), false);
	std::vector<int> check_face(mesh->faces_.size(), 0);
	//std::vector<bool> check_vertex_and_face(mesh->vertices_.size(), false);
	std::map<V3f, int> quick_v3f;
	int global_v = 0;
	bss.init(cell_num);
	//for(int i = cell_num - 1; i >= 0; i--)
	for (uint64_t i = 0; i < cell_num; i++)
	{
//		if (i % 100 == 0)
//			std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.local_bdy[i].resize(64, false);
		bss.Local_V3f[i].resize(64);
		bss.Local_to_Global[i].resize(64);
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = solid_corner_point(mesh->cells_[i].vertex_[vn]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].vertex_[vn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_vertex[mesh->cells_[i].vertex_[vn]] == false)
			{
				quick_v3f[v] = global_v++;
				bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
				check_vertex[mesh->cells_[i].vertex_[vn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = quick_v3f[v];
		}
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]);
			bss.Local_V3f[i][tempsize++] = v;
			quick_v3f[v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
		}
		V3f temp_v;
		//第一层
		temp_v = solid_edge_point(i, 9, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[9]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[9]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[5]] == false)
		if (check_edge[mesh->cells_[i].edges_[9]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[9]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 9, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[9]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[9]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[6]] == false)
		if (check_edge[mesh->cells_[i].edges_[9]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[9]]++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 1, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[1]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[5]] == false)
		if (check_edge[mesh->cells_[i].edges_[1]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[1]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 2, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[2]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[5]] == false)
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[2]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 2, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[2]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[6]] == false)
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[2]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 7, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[7]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[7]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[6]] == false)
		if (check_edge[mesh->cells_[i].edges_[7]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[7]]++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 1, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[1]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[1]] == false)
		if (check_edge[mesh->cells_[i].edges_[1]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[1]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 2, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[2]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[1]] == false)
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[2]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 2, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[2]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[2]] == false)
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[2]]++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 7, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[7]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[7]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[2]] == false)
		if (check_edge[mesh->cells_[i].edges_[7]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[7]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 8, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[8]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[8]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[1]] == false)
		if (check_edge[mesh->cells_[i].edges_[8]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[8]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 8, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[8]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[8]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[2]] == false)
		if (check_edge[mesh->cells_[i].edges_[8]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[8]]++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];
		//第二层
		temp_v = solid_edge_point(i, 2, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[2]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[5]] == false)
		if (check_edge[mesh->cells_[i].edges_[2]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[2]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 4, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[4]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[5]] == false)
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[4]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 4, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[4]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[6]] == false)
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[4]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 6, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[6]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[6]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[6]] == false)
		if (check_edge[mesh->cells_[i].edges_[6]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[6]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 0, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[0]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[5]] == false)
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[0]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[5]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 1, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[1]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[6]] == false)
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[1]]++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[6]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 0, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[0]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[1]] == false)
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[0]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 1, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[1]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[2]] == false)
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[1]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 0, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[0]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[1]] == false)
		if (check_edge[mesh->cells_[i].edges_[0]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[0]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 5, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[5]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[1]] == false)
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[5]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[1]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 5, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[5]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[2]] == false)
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[5]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 4, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[4]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[2]] == false)
		if (check_edge[mesh->cells_[i].edges_[4]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[4]]++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[2]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];
		//第三层
		temp_v = solid_edge_point(i, 2, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[2]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[4]] == false)
		if (check_edge[mesh->cells_[i].edges_[2]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[2]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 4, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[4]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[4]] == false)
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[4]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 4, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[4]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[7]] == false)
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[4]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 6, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[6]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[6]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[7]] == false)
		if (check_edge[mesh->cells_[i].edges_[6]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[6]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 0, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[0]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[4]] == false)
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[0]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 1, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[1]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[7]] == false)
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[1]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 0, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[0]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[0]] == false)
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[0]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 1, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[1]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[3]] == false)
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[1]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 0, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[0]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[0]] == false)
		if (check_edge[mesh->cells_[i].edges_[0]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[0]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 5, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[5]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[0]] == false)
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[5]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 5, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[5]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[3]] == false)
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[5]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 4, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[4]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[3]] == false)
		if (check_edge[mesh->cells_[i].edges_[4]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[4]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];
		//第4层
		temp_v = solid_edge_point(i, 11, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[11]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[11]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[4]] == false)
		if (check_edge[mesh->cells_[i].edges_[11]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[11]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 11, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[11]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[11]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[7]] == false)
		if (check_edge[mesh->cells_[i].edges_[11]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[11]]++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 3, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[3]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[4]] == false)
		if (check_edge[mesh->cells_[i].edges_[3]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[3]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 3, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[3]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[4]] == false)
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[3]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[4]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 3, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[3]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[7]] == false)
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[3]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 5, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[5]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[7]] == false)
		if (check_edge[mesh->cells_[i].edges_[5]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[5]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[7]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 3, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[3]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[0]] == false)
		if (check_edge[mesh->cells_[i].edges_[3]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[3]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 3, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[3]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[0]] == false)
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[3]] ++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_face_point(i, 3, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_face[mesh->cells_[i].faces_[3]] == false || check_vertex_and_face[mesh->cells_[i].vertex_[3]] == false)
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_face[mesh->cells_[i].faces_[3]]++;
			//check_vertex_and_face[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 5, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[5]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[3]] == false)
		if (check_edge[mesh->cells_[i].edges_[5]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[5]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 10, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[10]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[10]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[0]] == false)
		if (check_edge[mesh->cells_[i].edges_[10]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[10]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[0]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		temp_v = solid_edge_point(i, 10, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[10]))
			bss.local_bdy[i][tempsize - 1] = true;
		//if (check_edge[mesh->cells_[i].edges_[10]] == false || check_vertex_and_edge[mesh->cells_[i].vertex_[3]] == false)
		if (check_edge[mesh->cells_[i].edges_[10]] != 2)
		{
			quick_v3f[temp_v] = global_v++;
			bss.Local_to_Global[i][tempsize - 1] = global_v - 1;
			check_edge[mesh->cells_[i].edges_[10]] ++;
			//check_vertex_and_edge[mesh->cells_[i].vertex_[3]] = true;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = quick_v3f[temp_v];

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<int>>> sc(4);
		std::vector<std::vector<int>> sf(4);
		std::vector<int> sl(4);
		//第一层
		sl[0] = 5; sl[1] = 18; sl[2] = 22; sl[3] = 1;
		sf[0] = sl;
		sl[0] = 16; sl[1] = 19; sl[2] = 23; sl[3] = 26;
		sf[1] = sl;
		sl[0] = 17; sl[1] = 20; sl[2] = 24; sl[3] = 27;
		sf[2] = sl;
		sl[0] = 6; sl[1] = 21; sl[2] = 25; sl[3] = 2;
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = 28; sl[1] = 32; sl[2] = 34; sl[3] = 36;
		sf[0] = sl;
		sl[0] = 29; sl[1] = 13; sl[2] = 9; sl[3] = 37;
		sf[1] = sl;
		sl[0] = 30; sl[1] = 14; sl[2] = 10; sl[3] = 38;
		sf[2] = sl;
		sl[0] = 31; sl[1] = 33; sl[2] = 35; sl[3] = 39;
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = 40; sl[1] = 44; sl[2] = 46; sl[3] = 48;
		sf[0] = sl;
		sl[0] = 41; sl[1] = 12; sl[2] = 8; sl[3] = 49;
		sf[1] = sl;
		sl[0] = 42; sl[1] = 15; sl[2] = 11; sl[3] = 50;
		sf[2] = sl;
		sl[0] = 43; sl[1] = 45; sl[2] = 47; sl[3] = 51;
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = 4; sl[1] = 54; sl[2] = 58; sl[3] = 0;
		sf[0] = sl;
		sl[0] = 52; sl[1] = 55; sl[2] = 59; sl[3] = 62;
		sf[1] = sl;
		sl[0] = 53; sl[1] = 56; sl[2] = 60; sl[3] = 63;
		sf[2] = sl;
		sl[0] = 7; sl[1] = 57; sl[2] = 61; sl[3] = 3;
		sf[3] = sl;
		sc[3] = sf;
		bss.Control_Point[i] = sc;

		//构建局部的每个小单元的u,v,w三个方向的次数
		std::vector<int> order(3); order[0] = 3; order[1] = 3; order[2] = 3;
		bss.Local_Order[i] = order;

		//构建局部的每个小单元的u,v,w三个方向的节点向量
		std::vector<std::vector<int>> node(3);
		node[0].push_back(0); node[0].push_back(0); node[0].push_back(0); node[0].push_back(0);
		node[0].push_back(1); node[0].push_back(1); node[0].push_back(1); node[0].push_back(1);
		node[1].push_back(0); node[1].push_back(0); node[1].push_back(0); node[1].push_back(0);
		node[1].push_back(1); node[1].push_back(1); node[1].push_back(1); node[1].push_back(1);
		node[2].push_back(0); node[2].push_back(0); node[2].push_back(0); node[2].push_back(0);
		node[2].push_back(1); node[2].push_back(1); node[2].push_back(1); node[2].push_back(1);
		bss.Local_Node[i] = node;
	}
	int bdy_sum = 0;
	int v_bdy = 0;
	for (int i = 0; i < mesh->vertices_.size(); i++)
	{
		if (mesh->is_bdy(VertexHandle(i)))
		{
			bdy_sum++;
			v_bdy++;
		}
	}
	int e_bdy = 0;
	for (int i = 0; i < mesh->edges_.size(); i++)
	{
		if (mesh->is_bdy(EdgeHandle(i)))
		{
			bdy_sum += 2;
			e_bdy++;
		}
	}
	int f_bdy = 0;
	for (int i = 0; i < mesh->faces_.size(); i++)
	{
		if (mesh->is_bdy(FaceHandle(i)))
		{
			f_bdy++;
			bdy_sum += 4;
		}
	}
	bss.Global_size = global_v;
	bss.Global_size_without_bdy = bss.Global_size - bdy_sum;
}

void ApproximationCC::construct_bss_new_new()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<bool> check_vertex(mesh->vertices_.size(), false);
	std::vector<int> check_edge(mesh->edges_.size(), 0);
	std::vector<int> check_face(mesh->faces_.size(), 0);
	std::vector<V3f> slow_v3f;
	int global_v = 0;
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		if (i % 100 == 0)
			std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.local_bdy[i].resize(64, false);
		bss.Local_V3f[i].resize(64);
		bss.Local_to_Global[i].resize(64);
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = solid_corner_point(mesh->cells_[i].vertex_[vn]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].vertex_[vn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_vertex[mesh->cells_[i].vertex_[vn]] == false)
			{
				slow_v3f.push_back(v);
				bss.Local_to_Global[i][tempsize - 1] = global_v++;
				check_vertex[mesh->cells_[i].vertex_[vn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), v) - slow_v3f.begin();
		}
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]);
			bss.Local_V3f[i][tempsize++] = v;
			slow_v3f.push_back(v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
		}
		V3f temp_v;
		//第一层
		temp_v = solid_edge_point(i, 9, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[9]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[9]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[9]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 9, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[9]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[9]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[9]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 1, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[1]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[1]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 2, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[2]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 2, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[2]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 7, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[7]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[7]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[7]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 1, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[1]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[1]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 2, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[2]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 2, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[2]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[2]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 7, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[7]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[7]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[7]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 8, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[8]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[8]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[8]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 8, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[8]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[8]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[8]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();
		//第二层
		temp_v = solid_edge_point(i, 2, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[2]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[2]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 4, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[4]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 4, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[4]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 6, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[6]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[6]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[6]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 0, 5);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[0]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 1, 6);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[1]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 0, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[0]] ++;;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 1, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[1]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 0, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[0]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[0]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 5, 1);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 5, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 4, 2);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[4]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[4]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();
		//第三层
		temp_v = solid_edge_point(i, 2, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[2]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[2]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[2]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 4, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[4]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 4, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[4]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[4]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 6, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[6]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[6]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[6]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 0, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[0]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 1, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[1]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 0, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[0]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[0]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 1, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[1]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[1]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[1]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 0, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[0]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[0]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[0]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 5, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 5, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[5]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 4, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[4]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[4]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[4]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();
		//第4层
		temp_v = solid_edge_point(i, 11, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[11]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[11]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[11]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 11, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[11]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[11]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[11]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 3, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[3]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[3]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 3, 4);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[3]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 3, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[3]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 5, 7);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[5]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 3, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[3]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[3]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 3, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[3]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_face_point(i, 3, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].faces_[3]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_face[mesh->cells_[i].faces_[3]] != 4)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_face[mesh->cells_[i].faces_[3]]++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 5, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[5]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[5]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[5]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 10, 0);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[10]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[10]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[10]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		temp_v = solid_edge_point(i, 10, 3);
		bss.Local_V3f[i][tempsize++] = temp_v;
		if (mesh->is_bdy(mesh->cells_[i].edges_[10]))
			bss.local_bdy[i][tempsize - 1] = true;
		if (check_edge[mesh->cells_[i].edges_[10]] != 2)
		{
			slow_v3f.push_back(temp_v);
			bss.Local_to_Global[i][tempsize - 1] = global_v++;
			check_edge[mesh->cells_[i].edges_[10]] ++;
		}
		else
			bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), temp_v) - slow_v3f.begin();

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<int>>> sc(4);
		std::vector<std::vector<int>> sf(4);
		std::vector<int> sl(4);
		//第一层
		sl[0] = 5; sl[1] = 18; sl[2] = 22; sl[3] = 1;
		sf[0] = sl;
		sl[0] = 16; sl[1] = 19; sl[2] = 23; sl[3] = 26;
		sf[1] = sl;
		sl[0] = 17; sl[1] = 20; sl[2] = 24; sl[3] = 27;
		sf[2] = sl;
		sl[0] = 6; sl[1] = 21; sl[2] = 25; sl[3] = 2;
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = 28; sl[1] = 32; sl[2] = 34; sl[3] = 36;
		sf[0] = sl;
		sl[0] = 29; sl[1] = 13; sl[2] = 9; sl[3] = 37;
		sf[1] = sl;
		sl[0] = 30; sl[1] = 14; sl[2] = 10; sl[3] = 38;
		sf[2] = sl;
		sl[0] = 31; sl[1] = 33; sl[2] = 35; sl[3] = 39;
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = 40; sl[1] = 44; sl[2] = 46; sl[3] = 48;
		sf[0] = sl;
		sl[0] = 41; sl[1] = 12; sl[2] = 8; sl[3] = 49;
		sf[1] = sl;
		sl[0] = 42; sl[1] = 15; sl[2] = 11; sl[3] = 50;
		sf[2] = sl;
		sl[0] = 43; sl[1] = 45; sl[2] = 47; sl[3] = 51;
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = 4; sl[1] = 54; sl[2] = 58; sl[3] = 0;
		sf[0] = sl;
		sl[0] = 52; sl[1] = 55; sl[2] = 59; sl[3] = 62;
		sf[1] = sl;
		sl[0] = 53; sl[1] = 56; sl[2] = 60; sl[3] = 63;
		sf[2] = sl;
		sl[0] = 7; sl[1] = 57; sl[2] = 61; sl[3] = 3;
		sf[3] = sl;
		sc[3] = sf;
		bss.Control_Point[i] = sc;

		//构建局部的每个小单元的u,v,w三个方向的次数
		std::vector<int> order(3); order[0] = 3; order[1] = 3; order[2] = 3;
		bss.Local_Order[i] = order;

		//构建局部的每个小单元的u,v,w三个方向的节点向量
		std::vector<std::vector<int>> node(3);
		node[0].push_back(0); node[0].push_back(0); node[0].push_back(0); node[0].push_back(0);
		node[0].push_back(1); node[0].push_back(1); node[0].push_back(1); node[0].push_back(1);
		node[1].push_back(0); node[1].push_back(0); node[1].push_back(0); node[1].push_back(0);
		node[1].push_back(1); node[1].push_back(1); node[1].push_back(1); node[1].push_back(1);
		node[2].push_back(0); node[2].push_back(0); node[2].push_back(0); node[2].push_back(0);
		node[2].push_back(1); node[2].push_back(1); node[2].push_back(1); node[2].push_back(1);
		bss.Local_Node[i] = node;
	}
	int bdy_sum = 0;
	int v_bdy = 0;
	for (int i = 0; i < mesh->vertices_.size(); i++)
	{
		if (mesh->is_bdy(VertexHandle(i)))
		{
			bdy_sum++;
			v_bdy++;
		}
	}
	int e_bdy = 0;
	for (int i = 0; i < mesh->edges_.size(); i++)
	{
		if (mesh->is_bdy(EdgeHandle(i)))
		{
			bdy_sum += 2;
			e_bdy++;
		}
	}
	int f_bdy = 0;
	for (int i = 0; i < mesh->faces_.size(); i++)
	{
		if (mesh->is_bdy(FaceHandle(i)))
		{
			f_bdy++;
			bdy_sum += 4;
		}
	}
	bss.Global_size = global_v;
	bss.Global_size_without_bdy = bss.Global_size - bdy_sum;
}

void ApproximationCC::construct_bss1_new_new()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<bool> check_vertex(mesh->vertices_.size(), false);
	std::vector<bool> check_edge(mesh->edges_.size(), false);
	std::vector<bool> check_face(mesh->faces_.size(), false);
	std::vector<V3f> slow_v3f;
	int global_v = 0;
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		if (i % 1000 == 0)
			std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.local_bdy[i].resize(27, false);
		bss.Local_V3f[i].resize(27);
		bss.Local_to_Global[i].resize(27);
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			V3f v = mesh->vertices_[mesh->cells_[i].vertex_[vn]];
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].vertex_[vn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_vertex[mesh->cells_[i].vertex_[vn]] == false)
			{
				slow_v3f.push_back(v);
				bss.Local_to_Global[i][tempsize - 1] = global_v++;
				check_vertex[mesh->cells_[i].vertex_[vn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), v) - slow_v3f.begin();
		}
		for (int en = 0; en < 12; en++)
		{
			V3f v = cal_edge_midv(mesh->cells_[i].edges_[en]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].edges_[en]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_edge[mesh->cells_[i].edges_[en]] == false)
			{
				slow_v3f.push_back(v);
				bss.Local_to_Global[i][tempsize - 1] = global_v++;
				check_edge[mesh->cells_[i].edges_[en]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), v) - slow_v3f.begin();
		}
		for (int fn = 0; fn < 6; fn++)
		{
			V3f v = cal_face_midv(mesh->cells_[i].faces_[fn]);
			bss.Local_V3f[i][tempsize++] = v;
			if (mesh->is_bdy(mesh->cells_[i].faces_[fn]))
				bss.local_bdy[i][tempsize - 1] = true;
			if (check_face[mesh->cells_[i].faces_[fn]] == false)
			{
				slow_v3f.push_back(v);
				bss.Local_to_Global[i][tempsize - 1] = global_v++;
				check_face[mesh->cells_[i].faces_[fn]] = true;
			}
			else
				bss.Local_to_Global[i][tempsize - 1] = std::find(slow_v3f.begin(), slow_v3f.end(), v) - slow_v3f.begin();
		}
		V3f cell_midv = cal_cell_midv(CellHandle(i));
		bss.Local_V3f[i][tempsize++] = cell_midv;
		bss.Local_to_Global[i][tempsize - 1] = global_v++;

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<int>>> sc(3);
		std::vector<std::vector<int>> sf(3);
		std::vector<int> sl(3);
		//第一层
		sl[0] = 5; sl[1] = 9; sl[2] = 1;
		sf[0] = sl;
		sl[0] = 17; sl[1] = 22; sl[2] = 16;
		sf[1] = sl;
		sl[0] = 6; sl[1] = 15; sl[2] = 2;
		sf[2] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = 10; sl[1] = 20; sl[2] = 8;
		sf[0] = sl;
		sl[0] = 24; sl[1] = 26; sl[2] = 25;
		sf[1] = sl;
		sl[0] = 14; sl[1] = 21; sl[2] = 12;
		sf[2] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = 4; sl[1] = 11; sl[2] = 0;
		sf[0] = sl;
		sl[0] = 19; sl[1] = 23; sl[2] = 18;
		sf[1] = sl;
		sl[0] = 7; sl[1] = 13; sl[2] = 3;
		sf[2] = sl;
		sc[2] = sf;
		bss.Control_Point[i] = sc;

		//构建局部的每个小单元的u,v,w三个方向的次数
		std::vector<int> order(3); order[0] = 2; order[1] = 2; order[2] = 2;
		bss.Local_Order[i] = order;

		//构建局部的每个小单元的u,v,w三个方向的节点向量
		std::vector<std::vector<int>> node(3);
		node[0].push_back(0); node[0].push_back(0); node[0].push_back(0);
		node[0].push_back(1); node[0].push_back(1); node[0].push_back(1);
		node[1].push_back(0); node[1].push_back(0); node[1].push_back(0);
		node[1].push_back(1); node[1].push_back(1); node[1].push_back(1);
		node[2].push_back(0); node[2].push_back(0); node[2].push_back(0);
		node[2].push_back(1); node[2].push_back(1); node[2].push_back(1);
		bss.Local_Node[i] = node;
	}
	int bdy_sum = 0;
	for (int i = 0; i < mesh->vertices_.size(); i++)
	{
		if (mesh->is_bdy(VertexHandle(i)))
			bdy_sum++;
	}
	for (int i = 0; i < mesh->edges_.size(); i++)
	{
		if (mesh->is_bdy(EdgeHandle(i)))
			bdy_sum++;
	}
	for (int i = 0; i < mesh->faces_.size(); i++)
	{
		if (mesh->is_bdy(FaceHandle(i)))
			bdy_sum++;
	}
	bss.Global_size = global_v;
	bss.Global_size_without_bdy = bss.Global_size - bdy_sum;
}

V3f ApproximationCC::cal_corner_pointsharp(VertexHandle _v)
{
	bool check = false;//是否为尖锐边
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	for (auto ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (std::find(sharpedge.begin(), sharpedge.end(), *ve_it) != sharpedge.end())
		{
			check = true;
			break;
		}
	}
	V3f newv(0, 0, 0);
	if (check == true)
	{
		newv = newv + mesh->vertices_[_v] * 4;
		for (auto ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
		{
			if (std::find(sharpedge.begin(), sharpedge.end(), *ve_it) != sharpedge.end())
			{
				EdgeHandle e = *ve_it;
				if (mesh->edges_[e].from_h == _v)
					newv = newv + mesh->vertices_[mesh->edges_[e].to_h];
				else
					newv = newv + mesh->vertices_[mesh->edges_[e].from_h];
			}
		}
		newv = newv / 6;
	}
	else
		newv = cal_corner_point(_v);
	return newv;
}

V3f ApproximationCC::cal_edge_pointsharp(EdgeHandle _e, VertexHandle _v)
{
	bool check = false;//是否为尖锐边
	if (std::find(sharpedge.begin(), sharpedge.end(), _e) != sharpedge.end())
		check = true;
	V3f newv(0, 0, 0);
	if (check == true)
	{
		newv = newv + mesh->vertices_[_v] * 2;
		if (mesh->edges_[_e].from_h == _v)
			newv = newv + mesh->vertices_[mesh->edges_[_e].to_h];
		else
			newv = newv + mesh->vertices_[mesh->edges_[_e].from_h];
		newv = newv / 3;
	}
	else
		newv = cal_edge_point(_e, _v);
	return newv;
}

V3f ApproximationCC::cal_inter_pointsharp(FaceHandle _f, VertexHandle _v)
{
	bool check = false;//是否为尖锐边
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	for (auto ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (std::find(sharpedge.begin(), sharpedge.end(), *ve_it) != sharpedge.end())
		{
			check = true;
			break;
		}
	}
	V3f newv(0, 0, 0);
	if (check == true)
	{
		newv = newv + mesh->vertices_[_v] * 4;
		std::vector<EdgeHandle> fe = mesh->faces_[_f].edges_;
		for (std::vector<EdgeHandle>::iterator fe_it = fe.begin(); fe_it != fe.end(); fe_it++)
		{
			if (mesh->edges_[*fe_it].from_h == _v)
				newv = newv + mesh->vertices_[mesh->edges_[*fe_it].to_h] * 2;
			else if (mesh->edges_[*fe_it].to_h == _v)
				newv = newv + mesh->vertices_[mesh->edges_[*fe_it].from_h] * 2;
		}
		VertexHandle opv = mesh->find_op_vertex(_v, _f);
		newv = newv + mesh->vertices_[opv];
		newv = newv / 9;
	}
	else
		newv = cal_inter_point(_f, _v);
	return newv;
}