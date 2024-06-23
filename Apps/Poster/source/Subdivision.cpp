#include "Subdivision.h"
#include "BezierSolid.h"
#include "FileHandler.h"
#include <string>
#ifndef PI
#define PI  3.1425926
#endif // !PI


//-----Split-----//
HexV3fMesh& Split::Spliting()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<double> material_distribution1(8 * cell_num, 1);
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		std::vector<V3f>::iterator v_id;

		// 傻逼了 直接for不就好了

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//再插入每条边的中点 暴力插入 以免顺序错误
		//第一条边的中点
		V3f edge_mid_1;
		edge_mid_1.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[1]].x) * 1.0 / 2 * 1.0;
		edge_mid_1.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[1]].y) * 1.0 / 2 * 1.0;
		edge_mid_1.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[1]].z) * 1.0 / 2 * 1.0;
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2;
		edge_mid_2.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[5]].x) * 1.0 / 2 * 1.0;
		edge_mid_2.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[5]].y) * 1.0 / 2 * 1.0;
		edge_mid_2.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[5]].z) * 1.0 / 2 * 1.0;
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3;
		edge_mid_3.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[4]].x) * 1.0 / 2 * 1.0;
		edge_mid_3.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[4]].y) * 1.0 / 2 * 1.0;
		edge_mid_3.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[4]].z) * 1.0 / 2 * 1.0;
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4;
		edge_mid_4.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[0]].x) * 1.0 / 2 * 1.0;
		edge_mid_4.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[0]].y) * 1.0 / 2 * 1.0;
		edge_mid_4.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[0]].z) * 1.0 / 2 * 1.0;
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5;
		edge_mid_5.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_5.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_5.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6;
		edge_mid_6.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_6.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_6.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7;
		edge_mid_7.x = (mesh->vertices_[mesh->cells_[i].vertex_[6]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_7.y = (mesh->vertices_[mesh->cells_[i].vertex_[6]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_7.z = (mesh->vertices_[mesh->cells_[i].vertex_[6]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8;
		edge_mid_8.x = (mesh->vertices_[mesh->cells_[i].vertex_[2]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_8.y = (mesh->vertices_[mesh->cells_[i].vertex_[2]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_8.z = (mesh->vertices_[mesh->cells_[i].vertex_[2]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9;
		edge_mid_9.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_9.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_9.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10;
		edge_mid_10.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_10.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_10.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11;
		edge_mid_11.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[3]].x) * 1.0 / 2 * 1.0;
		edge_mid_11.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[3]].y) * 1.0 / 2 * 1.0;
		edge_mid_11.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[3]].z) * 1.0 / 2 * 1.0;
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12;
		edge_mid_12.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_12.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_12.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);

		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid;
			face_mid.x = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].x) * 1.0 / 4 * 1.0;
			face_mid.y = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].y) * 1.0 / 4 * 1.0;
			face_mid.z = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].z) * 1.0 / 4 * 1.0;
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}

		//插入单元的内心
		V3f cell_mid;
		double x_sum = 0;  double y_sum = 0;  double z_sum = 0;
		for (int j = 0; j < 8; j++)
		{
			x_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].x;
			y_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].y;
			z_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].z;
		}
		cell_mid.x = x_sum * 1.0 / 8 * 1.0;
		cell_mid.y = y_sum * 1.0 / 8 * 1.0;
		cell_mid.z = z_sum * 1.0 / 8 * 1.0;
		temp_vex[26] = newmesh->add_vertex(cell_mid);

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		if (material_distribution[i] != 1)
		{
			material_distribution1[newmesh->add_cell(sub_cell[0])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[1])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[2])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[3])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[4])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[5])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[6])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[7])] = 0.0001;
		}
		else
		{
			newmesh->add_cell(sub_cell[0]);
			newmesh->add_cell(sub_cell[1]);
			newmesh->add_cell(sub_cell[2]);
			newmesh->add_cell(sub_cell[3]);
			newmesh->add_cell(sub_cell[4]);
			newmesh->add_cell(sub_cell[5]);
			newmesh->add_cell(sub_cell[6]);
			newmesh->add_cell(sub_cell[7]);
		}
	}
	material_distribution = material_distribution1;
	return *newmesh;
}

HexV3fMesh& Split::Spliting_crease()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		std::vector<V3f>::iterator v_id;

		// 傻逼了 直接for不就好了

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]])) == newmesh->vertices_.end())
			{
				temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
			}
			else
			{
				temp_vex[mesh_vex_num] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
			}
		}
		//再插入每条边的中点 暴力插入 以免顺序错误
		//第一条边的中点
		V3f edge_mid_1;
		edge_mid_1.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[1]].x) * 1.0 / 2 * 1.0;
		edge_mid_1.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[1]].y) * 1.0 / 2 * 1.0;
		edge_mid_1.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[1]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_1)) == newmesh->vertices_.end())
		{
			temp_vex[8] = newmesh->add_vertex(edge_mid_1);
		}
		else
		{
			temp_vex[8] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		//判断该点是否在折痕上
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[0]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[1]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[1];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_1);
		}
		//第二条边的中点
		V3f edge_mid_2;
		edge_mid_2.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[5]].x) * 1.0 / 2 * 1.0;
		edge_mid_2.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[5]].y) * 1.0 / 2 * 1.0;
		edge_mid_2.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[5]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_2)) == newmesh->vertices_.end())
		{
			temp_vex[9] = newmesh->add_vertex(edge_mid_2);
		}
		else
		{
			temp_vex[9] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[1]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[5]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[2];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_2);
		}
		//第三条边的中点
		V3f edge_mid_3;
		edge_mid_3.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[4]].x) * 1.0 / 2 * 1.0;
		edge_mid_3.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[4]].y) * 1.0 / 2 * 1.0;
		edge_mid_3.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[4]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_3)) == newmesh->vertices_.end())
		{
			temp_vex[10] = newmesh->add_vertex(edge_mid_3);
		}
		else
		{
			temp_vex[10] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[5]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[4]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[3];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_3);
		}
		//第四条边的中点
		V3f edge_mid_4;
		edge_mid_4.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[0]].x) * 1.0 / 2 * 1.0;
		edge_mid_4.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[0]].y) * 1.0 / 2 * 1.0;
		edge_mid_4.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[0]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_4)) == newmesh->vertices_.end())
		{
			temp_vex[11] = newmesh->add_vertex(edge_mid_4);
		}
		else
		{
			temp_vex[11] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[4]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[0]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[4];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_4);
		}
		//第五条边的中点
		V3f edge_mid_5;
		edge_mid_5.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_5.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_5.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_5)) == newmesh->vertices_.end())
		{
			temp_vex[12] = newmesh->add_vertex(edge_mid_5);
		}
		else
		{
			temp_vex[12] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[3]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[2]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[5];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_5);
		}
		//第六条边的中点
		V3f edge_mid_6;
		edge_mid_6.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_6.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_6.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_6)) == newmesh->vertices_.end())
		{
			temp_vex[13] = newmesh->add_vertex(edge_mid_6);
		}
		else
		{
			temp_vex[13] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[3]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[7]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[6];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_6);
		}
		//第七条边的中点
		V3f edge_mid_7;
		edge_mid_7.x = (mesh->vertices_[mesh->cells_[i].vertex_[6]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_7.y = (mesh->vertices_[mesh->cells_[i].vertex_[6]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_7.z = (mesh->vertices_[mesh->cells_[i].vertex_[6]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_7)) == newmesh->vertices_.end())
		{
			temp_vex[14] = newmesh->add_vertex(edge_mid_7);
		}
		else
		{
			temp_vex[14] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[6]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[7]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[7];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_7);
		}
		//第八条边的中点
		V3f edge_mid_8;
		edge_mid_8.x = (mesh->vertices_[mesh->cells_[i].vertex_[2]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_8.y = (mesh->vertices_[mesh->cells_[i].vertex_[2]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_8.z = (mesh->vertices_[mesh->cells_[i].vertex_[2]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_8)) == newmesh->vertices_.end())
		{
			temp_vex[15] = newmesh->add_vertex(edge_mid_8);
		}
		else
		{
			temp_vex[15] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[2]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[6]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[8];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_8);
		}
		//第九条边的中点
		V3f edge_mid_9;
		edge_mid_9.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_9.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_9.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_9)) == newmesh->vertices_.end())
		{
			temp_vex[16] = newmesh->add_vertex(edge_mid_9);
		}
		else
		{
			temp_vex[16] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[1]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[2]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[9];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_9);
		}
		//第十条边的中点
		V3f edge_mid_10;
		edge_mid_10.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_10.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_10.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_10)) == newmesh->vertices_.end())
		{
			temp_vex[17] = newmesh->add_vertex(edge_mid_10);
		}
		else
		{
			temp_vex[17] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[5]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[6]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[10];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_10);
		}
		//第十一条边的中点
		V3f edge_mid_11;
		edge_mid_11.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[3]].x) * 1.0 / 2 * 1.0;
		edge_mid_11.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[3]].y) * 1.0 / 2 * 1.0;
		edge_mid_11.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[3]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_11)) == newmesh->vertices_.end())
		{
			temp_vex[18] = newmesh->add_vertex(edge_mid_11);
		}
		else
		{
			temp_vex[18] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[0]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[3]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[11];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_11);
		}
		//第十二条边的中点
		V3f edge_mid_12;
		edge_mid_12.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_12.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_12.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), edge_mid_12)) == newmesh->vertices_.end())
		{
			temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		}
		else
		{
			temp_vex[19] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		if ((std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[4]]) != v_crease.end()) &&
			(std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[mesh->cells_[i].vertex_[7]]) != v_crease.end()))
		{
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[12];
			if (ef.size() > 2)
				v_crease.push_back(edge_mid_12);
		}
		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid;
			face_mid.x = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].x) * 1.0 / 4 * 1.0;
			face_mid.y = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].y) * 1.0 / 4 * 1.0;
			face_mid.z = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].z) * 1.0 / 4 * 1.0;
			if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), face_mid)) == newmesh->vertices_.end())
			{
				temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
			}
			else
			{
				temp_vex[20 + mesh_face_num] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
			}
		}
		//插入单元的内心
		V3f cell_mid;
		double x_sum = 0;  double y_sum = 0;  double z_sum = 0;
		for (int j = 0; j < 8; j++)
		{
			x_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].x;
			y_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].y;
			z_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].z;
		}
		cell_mid.x = x_sum * 1.0 / 8 * 1.0;
		cell_mid.y = y_sum * 1.0 / 8 * 1.0;
		cell_mid.z = z_sum * 1.0 / 8 * 1.0;
		if ((v_id = std::find(newmesh->vertices_.begin(), newmesh->vertices_.end(), cell_mid)) == newmesh->vertices_.end())
		{
			temp_vex[26] = newmesh->add_vertex(cell_mid);
		}
		else
		{
			temp_vex[26] = VertexHandle(abs(std::distance(newmesh->vertices_.begin(), v_id)));
		}
		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

//-----Avgerage-----//
HexV3fMesh& Average::Averaging()
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh->vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		if (mesh->is_bdy(VertexHandle(vex_num)))
		{
			std::set<FaceHandle, compare_OVM> neighbor_vf = mesh->neighbor_v_f[vex_num];
			std::vector<FaceHandle> bdy_face;
			std::set<FaceHandle, compare_OVM>::iterator vf_id;
			for (vf_id = neighbor_vf.begin(); vf_id != neighbor_vf.end(); vf_id++)
			{
				if (mesh->is_bdy(*vf_id))
				{
					bdy_face.push_back(*vf_id);
				}
			}
			int vf_num = bdy_face.size();
			std::vector<V3f> vf_cen(vf_num);
			std::vector<FaceHandle>::iterator bdy_face_id;
			int u = 0;
			for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++)
			{
				std::vector<VertexHandle> vf = mesh->faces_[*bdy_face_id].vs_f;
				vf_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x) / 4 * 1.0;
				vf_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y) / 4 * 1.0;
				vf_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z) / 4 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex.x += vf_cen[vf_cen_num].x;
				new_vex.y += vf_cen[vf_cen_num].y;
				new_vex.z += vf_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num * 1.0;
			new_vex.y /= vf_num * 1.0;
			new_vex.z /= vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
		else
		{
			std::set<CellHandle, compare_OVM> neighbor_vc = mesh->neighbor_v_c[vex_num];
			std::set<CellHandle, compare_OVM>::iterator vc_id;
			int vf_num = neighbor_vc.size();
			std::vector<V3f> vc_cen(vf_num);
			int u = 0;
			for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++)
			{
				std::vector<VertexHandle> vf = mesh->cells_[*vc_id].vertex_;
				vc_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x +
					mesh->vertices_[vf[4]].x + mesh->vertices_[vf[5]].x +
					mesh->vertices_[vf[6]].x + mesh->vertices_[vf[7]].x) / 8 * 1.0;
				vc_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y +
					mesh->vertices_[vf[4]].y + mesh->vertices_[vf[5]].y +
					mesh->vertices_[vf[6]].y + mesh->vertices_[vf[7]].y) / 8 * 1.0;
				vc_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z +
					mesh->vertices_[vf[4]].z + mesh->vertices_[vf[5]].z +
					mesh->vertices_[vf[6]].z + mesh->vertices_[vf[7]].z) / 8 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex.x += vc_cen[vf_cen_num].x;
				new_vex.y += vc_cen[vf_cen_num].y;
				new_vex.z += vc_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num * 1.0;
			new_vex.y /= vf_num * 1.0;
			new_vex.z /= vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
	}
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		mesh->vertices_[vex_num] = changed_vex[vex_num];
	}
	return *mesh;
}

HexV3fMesh& Average::Averaging_creases()
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh->vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		if (mesh->is_bdy(VertexHandle(vex_num)))
		{
			if (std::find(v_crease.begin(), v_crease.end(), mesh->vertices_[VertexHandle(vex_num)]) != v_crease.end())
			{
				std::vector<EdgeHandle> ve;
				std::set<EdgeHandle, compare_OVM> crease_ve = mesh->neighbor_v[VertexHandle(vex_num)];
				//只要边界的
				for (std::set<EdgeHandle, compare_OVM>::iterator it = crease_ve.begin(); it != crease_ve.end(); it++)
				{
					if (mesh->is_bdy(*it))
					{
						ve.push_back(*it);
					}
				}
				std::vector<V3f> temp_vex;
				for (std::vector<EdgeHandle>::iterator it = ve.begin(); it != ve.end(); it++)
				{
					V3f v;
					v.x = (mesh->vertices_[mesh->edges_[*it].from_h].x + mesh->vertices_[mesh->edges_[*it].to_h].x) / 2;
					v.y = (mesh->vertices_[mesh->edges_[*it].from_h].y + mesh->vertices_[mesh->edges_[*it].to_h].y) / 2;
					v.z = (mesh->vertices_[mesh->edges_[*it].from_h].z + mesh->vertices_[mesh->edges_[*it].to_h].z) / 2;
					temp_vex.push_back(v);
				}
				V3f result_v;
				int num = temp_vex.size();
				for (int temp_vex_id = 0; temp_vex_id < num; temp_vex_id++)
				{
					result_v.x += temp_vex[temp_vex_id].x;
					result_v.y += temp_vex[temp_vex_id].y;
					result_v.z += temp_vex[temp_vex_id].z;
				}
				result_v.x /= num;
				result_v.y /= num;
				result_v.z /= num;
				changed_vex[vex_num] = result_v;
				v_crease.push_back(result_v);
			}
			else if (std::find(v_c_crease.begin(), v_c_crease.end(), mesh->vertices_[VertexHandle(vex_num)]) != v_c_crease.end())
			{
				changed_vex[vex_num] = mesh->vertices_[vex_num];
			}
			else
			{
				std::set<FaceHandle, compare_OVM> neighbor_vf = mesh->neighbor_v_f[vex_num];
				std::vector<FaceHandle> bdy_face;
				std::set<FaceHandle, compare_OVM>::iterator vf_id;
				for (vf_id = neighbor_vf.begin(); vf_id != neighbor_vf.end(); vf_id++)
				{
					if (mesh->is_bdy(*vf_id))
					{
						bdy_face.push_back(*vf_id);
					}
				}
				int vf_num = bdy_face.size();
				std::vector<V3f> vf_cen(vf_num);
				std::vector<FaceHandle>::iterator bdy_face_id;
				int u = 0;
				for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++)
				{
					std::vector<VertexHandle> vf = mesh->faces_[*bdy_face_id].vs_f;
					vf_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
						+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x) / 4;
					vf_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
						+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y) / 4;
					vf_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
						+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z) / 4;
					u++;
				}
				V3f new_vex;
				for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
				{
					new_vex.x += vf_cen[vf_cen_num].x;
					new_vex.y += vf_cen[vf_cen_num].y;
					new_vex.z += vf_cen[vf_cen_num].z;
				}
				new_vex.x /= vf_num;
				new_vex.y /= vf_num;
				new_vex.z /= vf_num;
				changed_vex[vex_num] = new_vex;
			}
		}
		else
		{
			std::set<CellHandle, compare_OVM> neighbor_vc = mesh->neighbor_v_c[vex_num];
			std::set<CellHandle, compare_OVM>::iterator vc_id;
			int vf_num = neighbor_vc.size();
			std::vector<V3f> vc_cen(vf_num);
			int u = 0;
			for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++)
			{
				std::vector<VertexHandle> vf = mesh->cells_[*vc_id].vertex_;
				vc_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x +
					mesh->vertices_[vf[4]].x + mesh->vertices_[vf[5]].x +
					mesh->vertices_[vf[6]].x + mesh->vertices_[vf[7]].x) / 8 * 1.0;
				vc_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y +
					mesh->vertices_[vf[4]].y + mesh->vertices_[vf[5]].y +
					mesh->vertices_[vf[6]].y + mesh->vertices_[vf[7]].y) / 8 * 1.0;
				vc_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z +
					mesh->vertices_[vf[4]].z + mesh->vertices_[vf[5]].z +
					mesh->vertices_[vf[6]].z + mesh->vertices_[vf[7]].z) / 8 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex.x += vc_cen[vf_cen_num].x;
				new_vex.y += vc_cen[vf_cen_num].y;
				new_vex.z += vc_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num * 1.0;
			new_vex.y /= vf_num * 1.0;
			new_vex.z /= vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
	}
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		mesh->vertices_[vex_num] = changed_vex[vex_num];
	}
	return *mesh;
}

HexV3fMesh& Average::Averaging_singular()
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh->vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		if (mesh->is_bdy(VertexHandle(vex_num)))
		{
			std::set<FaceHandle, compare_OVM> neighbor_vf = mesh->neighbor_v_f[vex_num];
			std::vector<FaceHandle> bdy_face;
			std::set<FaceHandle, compare_OVM>::iterator vf_id;
			for (vf_id = neighbor_vf.begin(); vf_id != neighbor_vf.end(); vf_id++)
			{
				if (mesh->is_bdy(*vf_id))
				{
					bdy_face.push_back(*vf_id);
				}
			}
			std::vector<FaceHandle>::iterator bdy_face_id;
			if (std::find(v3f_singular.begin(), v3f_singular.end(), mesh->vertices_[VertexHandle(vex_num)]) != v3f_singular.end())
			{
				V3f v = mesh->vertices_[VertexHandle(vex_num)];
				std::cout << "奇异点" << mesh->vertices_[VertexHandle(vex_num)].x << " " << mesh->vertices_[VertexHandle(vex_num)].y << " " << mesh->vertices_[VertexHandle(vex_num)].z << std::endl;
				std::set<VertexHandle, compare_OVM> v_bdy;
				for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++)
				{
					std::vector<VertexHandle> temp;
					temp = mesh->faces_[*bdy_face_id].vs_f;
					for (std::vector<VertexHandle>::iterator it = temp.begin(); it != temp.end(); it++)
					{
						v_bdy.insert(*it);
					}
				}
				int v_bdy_size = v_bdy.size();
				std::cout << v_bdy_size << std::endl;
				V3f new_vertex;
				for (std::set<VertexHandle, compare_OVM>::iterator tt = v_bdy.begin(); tt != v_bdy.end(); tt++)
				{
					double coe1 = 1 * 1.0 / ((mesh->val(VertexHandle(vex_num), true)) * (pow(4, 2)));
					//std::cout << "coe1:" << coe1 << std::endl;
					double coe2 = pow(3, (2 - mesh->dim(VertexHandle(vex_num), *tt)));
					//std::cout << "coe2:" << coe2 << std::endl;
					double coe3 = mesh->val(VertexHandle(vex_num), *tt, true);
					//std::cout << "coe3:" << coe3 << std::endl;
					double weight = coe1 * coe2 * coe3;
					//std::cout <<"weight:"<< weight << std::endl;
					new_vertex.x += weight * mesh->vertices_[*tt].x;
					new_vertex.y += weight * mesh->vertices_[*tt].y;
					new_vertex.z += weight * mesh->vertices_[*tt].z;
				}
				std::cout << "处理后的奇异点" << new_vertex.x << " " << new_vertex.y << " " << new_vertex.z << std::endl;
				changed_vex[vex_num] = new_vertex;
				v3f_singular.push_back(new_vertex);
				//经过测试 表面奇异点也根据全局计算
				//std::cout << "奇异点" << mesh->vertices_[VertexHandle(vex_num)].x << " " << mesh->vertices_[VertexHandle(vex_num)].y << " " << mesh->vertices_[VertexHandle(vex_num)].z << std::endl;
				//std::set<VertexHandle, compare_OVM> v_in;
				//for (auto v_in_id = mesh->neighbor_v_c[VertexHandle(vex_num)].begin();v_in_id != mesh->neighbor_v_c[VertexHandle(vex_num)].end();v_in_id++)
				//{
				//	std::vector<VertexHandle> temp;
				//	temp = mesh->cells_[*v_in_id].vertex_;
				//	for (std::vector<VertexHandle>::iterator it = temp.begin();it != temp.end();it++)
				//	{
				//		v_in.insert(*it);
				//	}
				//}
				//int v_in_size = v_in.size();
				//std::cout << v_in_size << std::endl;
				//V3f new_vertex;
				//for (std::set<VertexHandle, compare_OVM>::iterator tt = v_in.begin();tt != v_in.end();tt++)
				//{
				//	double coe1 = 1 * 1.0 / ((mesh->val(VertexHandle(vex_num), false)) * (pow(4, 3)));
				//	//std::cout << "coe1:" << coe1 << std::endl;
				//	double coe2 = pow(3, (3 - mesh->dim(VertexHandle(vex_num), *tt)));
				//	//std::cout << "coe2:" << coe2 << std::endl;
				//	double coe3 = mesh->val(VertexHandle(vex_num), *tt, false);
				//	//std::cout << "coe3:" << coe3 << std::endl;
				//	double weight = coe1 * coe2 * coe3;
				//	//std::cout << "weight:" << weight << std::endl;
				//	new_vertex.x += weight * mesh->vertices_[*tt].x;
				//	new_vertex.y += weight * mesh->vertices_[*tt].y;
				//	new_vertex.z += weight * mesh->vertices_[*tt].z;
				//}
				//std::cout << "处理后的奇异点" << new_vertex.x << " " << new_vertex.y << " " << new_vertex.z << std::endl;
				//changed_vex[vex_num] = new_vertex;
				//v3f_singular.push_back(new_vertex);
			}
			else
			{
				int vf_num = bdy_face.size();
				std::vector<V3f> vf_cen(vf_num);
				int u = 0;
				for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++)
				{
					std::vector<VertexHandle> vf = mesh->faces_[*bdy_face_id].vs_f;
					vf_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
						+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x) / 4 * 1.0;
					vf_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
						+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y) / 4 * 1.0;
					vf_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
						+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z) / 4 * 1.0;
					u++;
				}
				V3f new_vex;
				for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
				{
					new_vex.x += vf_cen[vf_cen_num].x;
					new_vex.y += vf_cen[vf_cen_num].y;
					new_vex.z += vf_cen[vf_cen_num].z;
				}
				new_vex.x /= vf_num * 1.0;
				new_vex.y /= vf_num * 1.0;
				new_vex.z /= vf_num * 1.0;
				changed_vex[vex_num] = new_vex;
			}
		}
		else
		{
			if (std::find(v3f_singular.begin(), v3f_singular.end(), mesh->vertices_[VertexHandle(vex_num)]) != v3f_singular.end())
			{
				std::cout << "奇异点" << mesh->vertices_[VertexHandle(vex_num)].x << " " << mesh->vertices_[VertexHandle(vex_num)].y << " " << mesh->vertices_[VertexHandle(vex_num)].z << std::endl;
				std::set<VertexHandle, compare_OVM> v_in;
				for (auto v_in_id = mesh->neighbor_v_c[VertexHandle(vex_num)].begin(); v_in_id != mesh->neighbor_v_c[VertexHandle(vex_num)].end(); v_in_id++)
				{
					std::vector<VertexHandle> temp;
					temp = mesh->cells_[*v_in_id].vertex_;
					for (std::vector<VertexHandle>::iterator it = temp.begin(); it != temp.end(); it++)
					{
						v_in.insert(*it);
					}
				}
				int v_in_size = v_in.size();
				std::cout << v_in_size << std::endl;
				V3f new_vertex;
				for (std::set<VertexHandle, compare_OVM>::iterator tt = v_in.begin(); tt != v_in.end(); tt++)
				{
					double coe1 = 1 * 1.0 / ((mesh->val(VertexHandle(vex_num), false)) * (pow(4, 3)));
					//std::cout << "coe1:" << coe1 << std::endl;
					double coe2 = pow(3, (3 - mesh->dim(VertexHandle(vex_num), *tt)));
					//std::cout << "coe2:" << coe2 << std::endl;
					double coe3 = mesh->val(VertexHandle(vex_num), *tt, false);
					//std::cout << "coe3:" << coe3 << std::endl;
					double weight = coe1 * coe2 * coe3;
					//std::cout << "weight:" << weight << std::endl;
					new_vertex.x += weight * mesh->vertices_[*tt].x;
					new_vertex.y += weight * mesh->vertices_[*tt].y;
					new_vertex.z += weight * mesh->vertices_[*tt].z;
				}
				std::cout << "处理后的奇异点" << new_vertex.x << " " << new_vertex.y << " " << new_vertex.z << std::endl;
				changed_vex[vex_num] = new_vertex;
				v3f_singular.push_back(new_vertex);
			}
			else
			{
				std::set<CellHandle, compare_OVM> neighbor_vc = mesh->neighbor_v_c[vex_num];
				std::set<CellHandle, compare_OVM>::iterator vc_id;
				int vf_num = neighbor_vc.size();
				std::vector<V3f> vc_cen(vf_num);
				int u = 0;
				for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++)
				{
					std::vector<VertexHandle> vf = mesh->cells_[*vc_id].vertex_;
					vc_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
						+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x +
						mesh->vertices_[vf[4]].x + mesh->vertices_[vf[5]].x +
						mesh->vertices_[vf[6]].x + mesh->vertices_[vf[7]].x) / 8 * 1.0;
					vc_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
						+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y +
						mesh->vertices_[vf[4]].y + mesh->vertices_[vf[5]].y +
						mesh->vertices_[vf[6]].y + mesh->vertices_[vf[7]].y) / 8 * 1.0;
					vc_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
						+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z +
						mesh->vertices_[vf[4]].z + mesh->vertices_[vf[5]].z +
						mesh->vertices_[vf[6]].z + mesh->vertices_[vf[7]].z) / 8 * 1.0;
					u++;
				}
				V3f new_vex;
				for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
				{
					new_vex.x += vc_cen[vf_cen_num].x;
					new_vex.y += vc_cen[vf_cen_num].y;
					new_vex.z += vc_cen[vf_cen_num].z;
				}
				new_vex.x /= vf_num * 1.0;
				new_vex.y /= vf_num * 1.0;
				new_vex.z /= vf_num * 1.0;
				changed_vex[vex_num] = new_vex;
			}
		}
	}
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		mesh->vertices_[vex_num] = changed_vex[vex_num];
	}
	return *mesh;
}

HexV3fMesh& Average::Averaging_test()
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh->vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		/*if (mesh->is_bdy(VertexHandle(vex_num)))
		{
			std::set<FaceHandle, compare_OVM> neighbor_vf = mesh->neighbor_v_f[vex_num];
			std::vector<FaceHandle> bdy_face;
			std::set<FaceHandle, compare_OVM>::iterator vf_id;
			for (vf_id = neighbor_vf.begin();vf_id != neighbor_vf.end();vf_id++)
			{
				if (mesh->is_bdy(*vf_id))
				{
					bdy_face.push_back(*vf_id);
				}
			}
			int vf_num = bdy_face.size();
			std::vector<V3f> vf_cen(vf_num);
			std::vector<FaceHandle>::iterator bdy_face_id;
			int u = 0;
			for (bdy_face_id = bdy_face.begin();bdy_face_id != bdy_face.end();bdy_face_id++)
			{
				std::vector<VertexHandle> vf = mesh->faces_[*bdy_face_id].vs_f;
				vf_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x) / 4 * 1.0;
				vf_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y) / 4 * 1.0;
				vf_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z) / 4 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0;vf_cen_num < vf_num;vf_cen_num++)
			{
				new_vex.x += vf_cen[vf_cen_num].x;
				new_vex.y += vf_cen[vf_cen_num].y;
				new_vex.z += vf_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num*1.0;
			new_vex.y /= vf_num*1.0;
			new_vex.z /= vf_num*1.0;
			changed_vex[vex_num] = new_vex;
		}
		else
		{*/
		std::set<CellHandle, compare_OVM> neighbor_vc = mesh->neighbor_v_c[vex_num];
		std::set<CellHandle, compare_OVM>::iterator vc_id;
		int vf_num = neighbor_vc.size();
		std::vector<V3f> vc_cen(vf_num);
		int u = 0;
		for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++)
		{
			std::vector<VertexHandle> vf = mesh->cells_[*vc_id].vertex_;
			vc_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
				+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x +
				mesh->vertices_[vf[4]].x + mesh->vertices_[vf[5]].x +
				mesh->vertices_[vf[6]].x + mesh->vertices_[vf[7]].x) / 8 * 1.0;
			vc_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
				+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y +
				mesh->vertices_[vf[4]].y + mesh->vertices_[vf[5]].y +
				mesh->vertices_[vf[6]].y + mesh->vertices_[vf[7]].y) / 8 * 1.0;
			vc_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
				+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z +
				mesh->vertices_[vf[4]].z + mesh->vertices_[vf[5]].z +
				mesh->vertices_[vf[6]].z + mesh->vertices_[vf[7]].z) / 8 * 1.0;
			u++;
		}
		V3f new_vex;
		for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
		{
			new_vex.x += vc_cen[vf_cen_num].x;
			new_vex.y += vc_cen[vf_cen_num].y;
			new_vex.z += vc_cen[vf_cen_num].z;
		}
		new_vex.x /= vf_num * 1.0;
		new_vex.y /= vf_num * 1.0;
		new_vex.z /= vf_num * 1.0;
		changed_vex[vex_num] = new_vex;
	}
	//}
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		mesh->vertices_[vex_num] = changed_vex[vex_num];
	}
	return *mesh;
}

//-----Subdivision-----//
HexV3fMesh& Subdivision::Subdivising()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		std::cout << "Subdivising:" << i << "/" << cell_num << std::endl;

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_edge_mid(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_edge_mid(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_edge_mid(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_edge_mid(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_edge_mid(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_edge_mid(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_edge_mid(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_edge_mid(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_edge_mid(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_edge_mid(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_edge_mid(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_edge_mid(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

V3f Subdivision::cal_cell_mid(const int cell_num)
{
	V3f cell_mid;
	std::vector<V3f> cmp_vertex;
	std::vector<V3f> cmq_vertex;
	if (mesh->is_bdy(CellHandle(cell_num)))
	{
		////加入p点集合
		//std::vector<VertexHandle> cm_ver = mesh->cells_[cell_num].vertex_;
		//for (std::vector<VertexHandle>::iterator cmv_it = cm_ver.begin();cmv_it != cm_ver.end();cmv_it++)
		//{
		//	cmp_vertex.push_back(mesh->vertices_[*cmv_it]);
		//}
		////加入q点集合
		//std::vector<FaceHandle> cm_f = mesh->cells_[cell_num].faces_;
		//std::vector<FaceHandle> bdy_and_op_f;//只算同一层的
		////第一遍遍历寻找不加入q点的面
		//for (std::vector<FaceHandle>::iterator cmf_it = cm_f.begin();cmf_it != cm_f.end();cmf_it++)
		//{
		//	if (mesh->is_bdy(*cmf_it))
		//	{
		//		bdy_and_op_f.push_back(*cmf_it);
		//		bdy_and_op_f.push_back(mesh->get_opposite_face_in_cell(*cmf_it, CellHandle(cell_num)));
		//	}
		//}
		////第二遍遍历进行加入q点集合
		//for (std::vector<FaceHandle>::iterator cmf_it = cm_f.begin();cmf_it != cm_f.end();cmf_it++)
		//{
		//	if (std::find(bdy_and_op_f.begin(), bdy_and_op_f.end(), *cmf_it) == bdy_and_op_f.end())
		//	{
		//		std::set<CellHandle, compare_OVM> cmf_c = mesh->neighbor_f[*cmf_it];
		//		CellHandle cmfc;
		//		for (std::set<CellHandle, compare_OVM>::iterator cmf_c_it = cmf_c.begin();cmf_c_it != cmf_c.end(); cmf_c_it++)
		//		{
		//			if ((*cmf_c_it) != CellHandle(cell_num))
		//			{
		//				cmfc = *cmf_c_it;
		//				break;
		//			}
		//		}
		//		FaceHandle cmf_op = mesh->get_opposite_face_in_cell(*cmf_it, cmfc);
		//		std::vector<VertexHandle> cmf_op_v = mesh->faces_[cmf_op].vs_f;
		//		for (std::vector<VertexHandle>::iterator cmf_op_v_it = cmf_op_v.begin();cmf_op_v_it != cmf_op_v.end();cmf_op_v_it++)
		//		{
		//			cmq_vertex.push_back(mesh->vertices_[*cmf_op_v_it]);
		//		}
		//	}
		//}
		////计算块的中心点
		//int cmp_size = cmp_vertex.size();
		//int cmq_size = cmq_vertex.size();
		//*double cmp_weight = (6 * w + 1) / 8;
		//double cmq_weight = w / 4;*/
	 // /*double cmp_weight = ((cmq_size*1.0 / (cmq_size + cmp_size)) * 8 * w + 1) / 8;
		//double cmq_weight = ((cmp_size*1.0 / (cmq_size + cmp_size)) * 8 * w) / 8;*/
		//double cmp_weight = (6 * w + 1) / 8;
		//double cmq_weight;
		//if (cmq_size == 0) cmq_weight = 0;
		//else cmq_weight = w / (cmq_size *1.0 / 6);
		//V3f cmp_sum;V3f cmq_sum;
		//for (int cmp_i = 0;cmp_i < cmp_size;cmp_i++)
		//{
		//	cmp_sum = cmp_sum + cmp_vertex[cmp_i];
		//}
		//for (int cmq_i = 0;cmq_i < cmq_size;cmq_i++)
		//{
		//	cmq_sum = cmq_sum + cmq_vertex[cmq_i];
		//}
		//cell_mid = cmp_sum*cmp_weight - cmq_sum*cmq_weight;
		std::vector<VertexHandle> cm_ver = mesh->cells_[cell_num].vertex_;
		for (std::vector<VertexHandle>::iterator cmv_it = cm_ver.begin(); cmv_it != cm_ver.end(); cmv_it++)
		{
			cmp_vertex.push_back(mesh->vertices_[*cmv_it]);
		}
		int cmp_size = cmp_vertex.size();
		V3f cmp_sum;
		for (int cmp_i = 0; cmp_i < cmp_size; cmp_i++)
		{
			cmp_sum = cmp_sum + cmp_vertex[cmp_i];
		}
		cell_mid = cmp_sum * (1.0 / 8);
	}
	else
	{
		//加入p点集合
		std::vector<VertexHandle> cm_ver = mesh->cells_[cell_num].vertex_;
		for (std::vector<VertexHandle>::iterator cmv_it = cm_ver.begin(); cmv_it != cm_ver.end(); cmv_it++)
		{
			cmp_vertex.push_back(mesh->vertices_[*cmv_it]);
		}
		//加入q点集合
		std::vector<FaceHandle> cm_f = mesh->cells_[cell_num].faces_;
		for (std::vector<FaceHandle>::iterator cmf_it = cm_f.begin(); cmf_it != cm_f.end(); cmf_it++)
		{
			std::set<CellHandle, compare_OVM> cmf_c = mesh->neighbor_f[*cmf_it];
			CellHandle cmfc;
			for (std::set<CellHandle, compare_OVM>::iterator cmf_c_it = cmf_c.begin(); cmf_c_it != cmf_c.end(); cmf_c_it++)
			{
				if ((*cmf_c_it) != CellHandle(cell_num))
				{
					cmfc = *cmf_c_it;
					break;
				}
			}
			FaceHandle cmf_op = mesh->get_opposite_face_in_cell(*cmf_it, cmfc);
			std::vector<VertexHandle> cmf_op_v = mesh->faces_[cmf_op].vs_f;
			for (std::vector<VertexHandle>::iterator cmf_op_v_it = cmf_op_v.begin(); cmf_op_v_it != cmf_op_v.end(); cmf_op_v_it++)
			{
				cmq_vertex.push_back(mesh->vertices_[*cmf_op_v_it]);
			}
		}
		//计算块的中心点
		int cmp_size = cmp_vertex.size();
		int cmq_size = cmq_vertex.size();
		double cmp_weight = (6 * w + 1) / 8;
		double cmq_weight = w / 4;
		V3f cmp_sum; V3f cmq_sum;
		for (int cmp_i = 0; cmp_i < cmp_size; cmp_i++)
		{
			cmp_sum = cmp_sum + cmp_vertex[cmp_i];
		}
		for (int cmq_i = 0; cmq_i < cmq_size; cmq_i++)
		{
			cmq_sum = cmq_sum + cmq_vertex[cmq_i];
		}
		cell_mid = cmp_sum * cmp_weight - cmq_sum * cmq_weight;
	}
	return cell_mid;
}

V3f Subdivision::cal_face_mid(const int cell_num, const int face_num)
{
	V3f face_mid;
	std::vector<V3f> fmp_vertex;
	std::vector<V3f> fmq_vertex;
	std::vector<V3f> fmr_vertex;
	std::vector<V3f> fms_vertex;
	std::vector<V3f> fmt_vertex;
	if (mesh->is_bdy(mesh->cells_[cell_num].faces_[face_num]))
	{
		//加入p点集合
		std::vector<VertexHandle> fm_ver = mesh->faces_[mesh->cells_[cell_num].faces_[face_num]].vs_f;
		for (std::vector<VertexHandle>::iterator fmv_it = fm_ver.begin(); fmv_it != fm_ver.end(); fmv_it++)
		{
			fmp_vertex.push_back(mesh->vertices_[*fmv_it]);
		}
		//加入r点集合
		std::vector<EdgeHandle> fm_f_e = mesh->faces_[mesh->cells_[cell_num].faces_[face_num]].edges_;
		for (std::vector<EdgeHandle>::iterator fm_f_e_it = fm_f_e.begin(); fm_f_e_it != fm_f_e.end(); fm_f_e_it++)
		{
			std::set<FaceHandle, compare_OVM> fm_f_e_f = mesh->neighbor_e[*fm_f_e_it];
			for (std::set<FaceHandle, compare_OVM>::iterator fm_f_e_f_it = fm_f_e_f.begin(); fm_f_e_f_it != fm_f_e_f.end(); fm_f_e_f_it++)
			{
				if (mesh->is_bdy(*fm_f_e_f_it) && (*fm_f_e_f_it != mesh->cells_[cell_num].faces_[face_num]))
				{
					EdgeHandle fm_f_e_f_ope = mesh->get_opposite_edge_in_face(*fm_f_e_it, *fm_f_e_f_it);
					fmr_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_ope].from_h]);
					fmr_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_ope].to_h]);
					break;
				}
			}
		}
		//计算面的中心点
		int fmp_size = fmp_vertex.size();
		int fmr_size = fmr_vertex.size();
		/*double fmp_weight = ((fmr_size*1.0 / (fmr_size + fmp_size)) * 4 * w + 1) / 4;
		double fmr_weight = (fmp_size*1.0 / (fmr_size + fmp_size)) * 4 * w / 4;*/
		double fmp_weight = (3 * w + 1) / 4;
		double fmr_weight = (3 * w) / 8;
		V3f fmp_sum; V3f fmr_sum;
		for (int fmp_i = 0; fmp_i < fmp_size; fmp_i++)
		{
			fmp_sum = fmp_sum + fmp_vertex[fmp_i];
		}
		for (int fmr_i = 0; fmr_i < fmr_size; fmr_i++)
		{
			fmr_sum = fmr_sum + fmr_vertex[fmr_i];
		}
		face_mid = fmp_sum * fmp_weight - fmr_sum * fmr_weight;
	}
	else if (check_face(mesh->cells_[cell_num].faces_[face_num]))//当该面在边界块上时，用线型细分
	{
		std::vector<VertexHandle> fm_ver = mesh->faces_[mesh->cells_[cell_num].faces_[face_num]].vs_f;
		for (std::vector<VertexHandle>::iterator fmv_it = fm_ver.begin(); fmv_it != fm_ver.end(); fmv_it++)
		{
			fmp_vertex.push_back(mesh->vertices_[*fmv_it]);
		}
		int fmp_size = fmp_vertex.size();
		V3f fmp_sum;
		for (int fmp_i = 0; fmp_i < fmp_size; fmp_i++)
		{
			fmp_sum = fmp_sum + fmp_vertex[fmp_i];
		}
		face_mid = fmp_sum * (1.0 / 4);
	}
	else
	{
		//加入p点集合
		std::vector<VertexHandle> fm_ver = mesh->faces_[mesh->cells_[cell_num].faces_[face_num]].vs_f;
		for (std::vector<VertexHandle>::iterator fmv_it = fm_ver.begin(); fmv_it != fm_ver.end(); fmv_it++)
		{
			fmp_vertex.push_back(mesh->vertices_[*fmv_it]);
		}
		//加入q点集合
		std::set<CellHandle, compare_OVM> fm_c = mesh->neighbor_f[mesh->cells_[cell_num].faces_[face_num]];
		for (std::set<CellHandle, compare_OVM>::iterator fm_c_it = fm_c.begin(); fm_c_it != fm_c.end(); fm_c_it++)
		{
			FaceHandle fm_c_op = mesh->get_opposite_face_in_cell(mesh->cells_[cell_num].faces_[face_num], *fm_c_it);
			std::vector<VertexHandle> fm_c_op_v = mesh->faces_[fm_c_op].vs_f;
			for (std::vector<VertexHandle>::iterator fm_c_op_v_it = fm_c_op_v.begin(); fm_c_op_v_it != fm_c_op_v.end(); fm_c_op_v_it++)
			{
				fmq_vertex.push_back(mesh->vertices_[*fm_c_op_v_it]);
			}
		}
		//加入r、s和t点集合
		std::set<FaceHandle, compare_OVM> fm_c_fall;//用于插入面相邻块的所有面
		for (std::set<CellHandle, compare_OVM>::iterator fm_c_it = fm_c.begin(); fm_c_it != fm_c.end(); fm_c_it++)
		{
			std::vector<FaceHandle> fm_c_f = mesh->cells_[*fm_c_it].faces_;
			for (std::vector<FaceHandle>::iterator fm_c_f_it = fm_c_f.begin(); fm_c_f_it != fm_c_f.end(); fm_c_f_it++)
			{
				fm_c_fall.insert(*fm_c_f_it);
			}
		}
		std::vector<EdgeHandle> fm_f_e = mesh->faces_[mesh->cells_[cell_num].faces_[face_num]].edges_;
		int N = 0;//记录奇异边个数
		for (std::vector<EdgeHandle>::iterator fm_f_e_it = fm_f_e.begin(); fm_f_e_it != fm_f_e.end(); fm_f_e_it++)
		{
			if (!mesh->is_singular(*fm_f_e_it))
			{
				//加入点r
				std::set<FaceHandle, compare_OVM> fm_f_e_f = mesh->neighbor_e[*fm_f_e_it];
				for (std::set<FaceHandle, compare_OVM>::iterator fm_f_e_f_it = fm_f_e_f.begin(); fm_f_e_f_it != fm_f_e_f.end(); fm_f_e_f_it++)
				{
					if (std::find(fm_c_fall.begin(), fm_c_fall.end(), *fm_f_e_f_it) == fm_c_fall.end())
					{
						EdgeHandle fm_f_ope = mesh->get_opposite_edge_in_face(*fm_f_e_it, *fm_f_e_f_it);
						fmr_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_ope].from_h]);
						fmr_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_ope].to_h]);
						break;
					}
				}
				//加入点s
				//CellHandle fm_ci[2];//两个与该面相邻的块
				//int fm_c_size = fm_c.size();int fm_c_n = 0;
				//for (std::set<CellHandle, compare_OVM>::iterator fm_c_it = fm_c.begin();fm_c_it != fm_c.end();fm_c_it++)
				//{
				//	fm_ci[fm_c_n] = *fm_c_it;
				//	fm_c_n++;
				//}
				//std::set<FaceHandle, compare_OVM> fm_f_e_f = mesh->neighbor_e[*fm_f_e_it];
				for (std::set<FaceHandle, compare_OVM>::iterator fm_f_e_f_it = fm_f_e_f.begin(); fm_f_e_f_it != fm_f_e_f.end(); fm_f_e_f_it++)
				{
					if (*fm_f_e_f_it == mesh->cells_[cell_num].faces_[face_num])
						continue;
					//else if (mesh->is_f_in_c(*fm_f_e_f_it, fm_ci[0]) || mesh->is_f_in_c(*fm_f_e_f_it, fm_ci[1]))
					else if (std::find(fm_c_fall.begin(), fm_c_fall.end(), *fm_f_e_f_it) != fm_c_fall.end())
						continue;
					else
					{
						std::set<CellHandle, compare_OVM> fm_f_e_f_c = mesh->neighbor_f[*fm_f_e_f_it];//这两个块包含点s
						for (std::set<CellHandle, compare_OVM>::iterator fm_f_e_f_c_it = fm_f_e_f_c.begin(); fm_f_e_f_c_it != fm_f_e_f_c.end(); fm_f_e_f_c_it++)
						{
							FaceHandle fm_f_e_f_c_opf = mesh->get_opposite_face_in_cell(*fm_f_e_f_it, *fm_f_e_f_c_it);
							std::vector<VertexHandle> fm_f_e_f_c_opf_v = mesh->faces_[fm_f_e_f_c_opf].vs_f;
							for (std::vector<VertexHandle>::iterator fm_f_e_f_c_opf_v_it = fm_f_e_f_c_opf_v.begin();
								fm_f_e_f_c_opf_v_it != fm_f_e_f_c_opf_v.end(); fm_f_e_f_c_opf_v_it++)
							{
								if (std::find(fmr_vertex.begin(), fmr_vertex.end(), mesh->vertices_[*fm_f_e_f_c_opf_v_it]) == fmr_vertex.end())
									fms_vertex.push_back(mesh->vertices_[*fm_f_e_f_c_opf_v_it]);
							}
						}
					}
				}
			}
			else//若该奇异边度大于4.则r点分裂，s点不变，若该奇异边度数小于4，则s点合并，r点消失。
			{
				N += 1;
				CellHandle fm_ci[2];//两个与该 面相邻的块
				int fm_c_size = fm_c.size(); int fm_c_n = 0;
				for (std::set<CellHandle, compare_OVM>::iterator fm_c_it = fm_c.begin(); fm_c_it != fm_c.end(); fm_c_it++)
				{
					fm_ci[fm_c_n] = *fm_c_it;
					fm_c_n++;
				}
				std::set<FaceHandle, compare_OVM> fm_f_e_f = mesh->neighbor_e[*fm_f_e_it];
				int sin_size = fm_f_e_f.size();
				for (std::set<FaceHandle, compare_OVM>::iterator fm_f_e_f_it = fm_f_e_f.begin(); fm_f_e_f_it != fm_f_e_f.end(); fm_f_e_f_it++)
				{
					if (*fm_f_e_f_it == mesh->cells_[cell_num].faces_[face_num])
						continue;
					else if (mesh->is_f_in_c(*fm_f_e_f_it, fm_ci[0]))
					{
						std::set<CellHandle, compare_OVM> fm_f_e_f_c = mesh->neighbor_f[*fm_f_e_f_it];
						for (std::set<CellHandle, compare_OVM>::iterator fm_f_e_f_c_it = fm_f_e_f_c.begin(); fm_f_e_f_c_it != fm_f_e_f_c.end(); fm_f_e_f_c_it++)
						{
							if (*fm_f_e_f_c_it == fm_ci[0])
								continue;
							else
							{
								//该块为t、s点所在块
								std::vector<FaceHandle> fm_f_e_f_c_f = mesh->cells_[*fm_f_e_f_c_it].faces_;
								for (std::vector<FaceHandle>::iterator fm_f_e_f_c_f_it = fm_f_e_f_c_f.begin(); fm_f_e_f_c_f_it != fm_f_e_f_c_f.end(); fm_f_e_f_c_f_it++)
								{
									if ((std::find(fm_f_e_f.begin(), fm_f_e_f.end(), *fm_f_e_f_c_f_it) != fm_f_e_f.end()) &&
										(*fm_f_e_f_c_f_it != *fm_f_e_f_it))//若边的度数小于4，则不存在r点的面
									{
										//该面为t点所在面
										//加入点t
										if (sin_size > 4)
										{
											EdgeHandle fm_f_e_f_c_f_e = mesh->get_opposite_edge_in_face(*fm_f_e_it, *fm_f_e_f_c_f_it);
											fmt_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_c_f_e].from_h]);
											fmt_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_c_f_e].to_h]);
										}
										//该面所在块的对面为点s所在面
										//加入点s
										FaceHandle fm_f_e_f_c_opf = mesh->get_opposite_face_in_cell(*fm_f_e_f_c_f_it, *fm_f_e_f_c_it);
										std::vector<VertexHandle> fm_f_e_f_c_opf_v = mesh->faces_[fm_f_e_f_c_opf].vs_f;
										for (std::vector<VertexHandle>::iterator fm_f_e_f_c_opf_v_it = fm_f_e_f_c_opf_v.begin();
											fm_f_e_f_c_opf_v_it != fm_f_e_f_c_opf_v.end(); fm_f_e_f_c_opf_v_it++)
										{
											if ((std::find(fmq_vertex.begin(), fmq_vertex.end(), mesh->vertices_[*fm_f_e_f_c_opf_v_it]) == fmq_vertex.end())
												&& (std::find(fms_vertex.begin(), fms_vertex.end(), mesh->vertices_[*fm_f_e_f_c_opf_v_it]) == fms_vertex.end()))
												fms_vertex.push_back(mesh->vertices_[*fm_f_e_f_c_opf_v_it]);
										}
										break;
									}
								}
							}
						}
					}
					else if (mesh->is_f_in_c(*fm_f_e_f_it, fm_ci[1]))
					{
						std::set<CellHandle, compare_OVM> fm_f_e_f_c = mesh->neighbor_f[*fm_f_e_f_it];
						for (std::set<CellHandle, compare_OVM>::iterator fm_f_e_f_c_it = fm_f_e_f_c.begin(); fm_f_e_f_c_it != fm_f_e_f_c.end(); fm_f_e_f_c_it++)
						{
							if (*fm_f_e_f_c_it == fm_ci[1])
								continue;
							else
							{
								//该块为t、s点所在块
								std::vector<FaceHandle> fm_f_e_f_c_f = mesh->cells_[*fm_f_e_f_c_it].faces_;
								for (std::vector<FaceHandle>::iterator fm_f_e_f_c_f_it = fm_f_e_f_c_f.begin(); fm_f_e_f_c_f_it != fm_f_e_f_c_f.end(); fm_f_e_f_c_f_it++)
								{
									if ((std::find(fm_f_e_f.begin(), fm_f_e_f.end(), *fm_f_e_f_c_f_it) != fm_f_e_f.end()) &&
										(*fm_f_e_f_c_f_it != *fm_f_e_f_it))//若边的度数小于4，则不存在r点的面
									{
										//该面为t点所在面
										//加入点t
										if (sin_size > 4)
										{
											EdgeHandle fm_f_e_f_c_f_e = mesh->get_opposite_edge_in_face(*fm_f_e_it, *fm_f_e_f_c_f_it);
											fmt_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_c_f_e].from_h]);
											fmt_vertex.push_back(mesh->vertices_[mesh->edges_[fm_f_e_f_c_f_e].to_h]);
										}
										//该面所在块的对面为点s所在面
										//加入点s
										FaceHandle fm_f_e_f_c_opf = mesh->get_opposite_face_in_cell(*fm_f_e_f_c_f_it, *fm_f_e_f_c_it);
										std::vector<VertexHandle> fm_f_e_f_c_opf_v = mesh->faces_[fm_f_e_f_c_opf].vs_f;
										for (std::vector<VertexHandle>::iterator fm_f_e_f_c_opf_v_it = fm_f_e_f_c_opf_v.begin();
											fm_f_e_f_c_opf_v_it != fm_f_e_f_c_opf_v.end(); fm_f_e_f_c_opf_v_it++)
										{
											if ((std::find(fmq_vertex.begin(), fmq_vertex.end(), mesh->vertices_[*fm_f_e_f_c_opf_v_it]) == fmq_vertex.end())
												&& (std::find(fms_vertex.begin(), fms_vertex.end(), mesh->vertices_[*fm_f_e_f_c_opf_v_it]) == fms_vertex.end()))
												fms_vertex.push_back(mesh->vertices_[*fm_f_e_f_c_opf_v_it]);
										}
										break;
									}
								}
							}
						}
					}
					else
						continue;
				}
			}
		}
		//计算面的中心点
		int fmp_size = fmp_vertex.size();
		int fmq_size = fmq_vertex.size();
		int fmr_size = fmr_vertex.size();
		int fms_size = fms_vertex.size();
		int fmt_size = fmt_vertex.size();

		double fmp_weight = (2 * w + 1) / 4;
		double fmq_weight = w / 4;
		double fmr_weight;
		double fms_weight;
		if (fmr_size == 0 && fmt_size == 0) fmr_weight = 0;
		else fmr_weight = w / ((fmt_size / 2 + fmr_size) / 2);
		if (fms_size == 0) fms_weight = 0;
		else fms_weight = w / (fms_size / 2);
		double fmt_weight = fmr_weight / 2;

		/*double fmp_weight = (2 * w + 1) / 4;
		double fmq_weight = w / 4;
		double fmr_weight = w / 4;
		double fms_weight = w / 8;
		double fmt_weight = w / 8;*/

		/*double fmp_weight, fmq_weight, fmr_weight, fms_weight, fmt_weight;
		fmp_weight = (2 * w + 1) / 4;
		fmq_weight = w / 4;
		if((fmt_size / 2 + fmr_size) !=  8)
		fmr_weight = w / 4;
		else
		{
		if (fmr_size == 0 && fmt_size == 0) fmr_weight = 0;
		else fmr_weight = w / ((fmt_size / 2 + fmr_size) / 2);
		}
		if(fms_size != 16)
		fms_weight = w / 8;
		else
		{
		if (fms_size == 0) fms_weight = 0;
		else fms_weight = w / (fms_size / 2);
		}
		fmt_weight = fmr_weight / 2;*/

		V3f fmp_sum; V3f fmq_sum; V3f fmr_sum; V3f fms_sum; V3f fmt_sum;
		for (int fmp_i = 0; fmp_i < fmp_size; fmp_i++)
		{
			fmp_sum = fmp_sum + fmp_vertex[fmp_i];
		}
		for (int fmq_i = 0; fmq_i < fmq_size; fmq_i++)
		{
			fmq_sum = fmq_sum + fmq_vertex[fmq_i];
		}
		for (int fmr_i = 0; fmr_i < fmr_size; fmr_i++)
		{
			fmr_sum = fmr_sum + fmr_vertex[fmr_i];
		}
		for (int fms_i = 0; fms_i < fms_size; fms_i++)
		{
			fms_sum = fms_sum + fms_vertex[fms_i];
		}
		for (int fmt_i = 0; fmt_i < fmt_size; fmt_i++)
		{
			fmt_sum = fmt_sum + fmt_vertex[fmt_i];
		}
		face_mid = fmp_sum * fmp_weight + fmq_sum * fmq_weight - fmr_sum * fmr_weight - fms_sum * fms_weight - fmt_sum * fmt_weight;
	}
	return face_mid;
}

V3f Subdivision::cal_edge_mid(const int cell_num, const int ver1, const int ver2)
{
	V3f edge_mid;
	std::vector<V3f> emp_vertex;
	std::vector<V3f> emq_vertex;
	std::vector<V3f> emr_vertex;
	std::vector<V3f> ems_vertex;
	//首先确定该边
	EdgeHandle em_e;
	bool is_find = false;
	std::vector<FaceHandle> em_c_f = mesh->cells_[cell_num].faces_;
	for (std::vector<FaceHandle>::iterator em_c_f_it = em_c_f.begin(); em_c_f_it != em_c_f.end(); em_c_f_it++)
	{
		if (is_find == true)
			break;
		std::vector<EdgeHandle> em_c_f_e = mesh->faces_[*em_c_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_c_f_e_it = em_c_f_e.begin(); em_c_f_e_it != em_c_f_e.end(); em_c_f_e_it++)
		{
			if (((mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver2]))
				|| ((mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver2])))
			{
				em_e = *em_c_f_e_it;
				is_find = true;
				break;
			}
		}
	}
	if (mesh->is_bdy(em_e))
	{
		const int N = mesh->neighbor_e[em_e].size();
		//表面的不算入奇异边
		//加入p点集合
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver1]]);
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver2]]);
		//加入q、r点集合
		std::set<FaceHandle, compare_OVM> em_e_f = mesh->neighbor_e[em_e];
		for (std::set<FaceHandle, compare_OVM>::iterator em_e_f_it = em_e_f.begin(); em_e_f_it != em_e_f.end(); em_e_f_it++)
		{
			if (mesh->is_bdy(*em_e_f_it))
			{
				//加入q点集合
				EdgeHandle em_e_f_ope = mesh->get_opposite_edge_in_face(em_e, *em_e_f_it);
				emq_vertex.push_back(mesh->vertices_[mesh->edges_[em_e_f_ope].from_h]);
				emq_vertex.push_back(mesh->vertices_[mesh->edges_[em_e_f_ope].to_h]);
				/*std::vector<VertexHandle> em_e_f_v = mesh->faces_[*em_e_f_it].vs_f;
				for (std::vector<VertexHandle>::iterator em_e_f_v_it = em_e_f_v.begin();em_e_f_v_it != em_e_f_v.end();em_e_f_v_it++)
				{
					if ((*em_e_f_v_it != mesh->cells_[cell_num].vertex_[ver2]) && (*em_e_f_v_it != mesh->cells_[cell_num].vertex_[ver1]))
						emq_vertex.push_back(mesh->vertices_[*em_e_f_v_it]);
				}*/
				//加入r点集合
				std::vector<EdgeHandle> em_e_f_e = mesh->faces_[*em_e_f_it].edges_;
				//EdgeHandle em_ope = mesh->get_opposite_edge_in_face(em_e,*em_e_f_it);
				//确定r点
				for (std::vector<EdgeHandle>::iterator em_e_f_e_it = em_e_f_e.begin(); em_e_f_e_it != em_e_f_e.end(); em_e_f_e_it++)
				{
					if ((*em_e_f_e_it != em_e) && (*em_e_f_e_it != em_e_f_ope))
					{
						//确定这条边上的q点
						VertexHandle em_e_f_e_v;
						if (std::find(emq_vertex.begin(), emq_vertex.end(), mesh->vertices_[mesh->edges_[*em_e_f_e_it].from_h]) != emq_vertex.end())
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].from_h;
						else
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].to_h;
						std::set<FaceHandle, compare_OVM> em_e_f_e_f = mesh->neighbor_e[*em_e_f_e_it];
						for (std::set<FaceHandle, compare_OVM>::iterator em_e_f_e_f_it = em_e_f_e_f.begin(); em_e_f_e_f_it != em_e_f_e_f.end(); em_e_f_e_f_it++)
						{
							if (mesh->is_bdy(*em_e_f_e_f_it) && (*em_e_f_e_f_it != *em_e_f_it))
							{
								//该面为r点所在面
								//确定r点所在边
								std::vector<EdgeHandle> em_e_f_e_f_e = mesh->faces_[*em_e_f_e_f_it].edges_;
								for (std::vector<EdgeHandle>::iterator em_e_f_e_f_e_it = em_e_f_e_f_e.begin(); em_e_f_e_f_e_it != em_e_f_e_f_e.end(); em_e_f_e_f_e_it++)
								{
									if (*em_e_f_e_f_e_it == *em_e_f_e_it)
										continue;
									else if (mesh->edges_[*em_e_f_e_f_e_it].from_h == em_e_f_e_v)
										emr_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_e_f_e_it].to_h]);
									else if (mesh->edges_[*em_e_f_e_f_e_it].to_h == em_e_f_e_v)
										emr_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_e_f_e_it].from_h]);
									else
										continue;
								}
							}
						}
					}
				}
			}
		}
		//计算边的中心点
		int emp_size = emp_vertex.size();
		int emq_size = emq_vertex.size();
		int emr_size = emr_vertex.size();
		double emp_weight = 1.0 / 2;
		double emq_weight = w / 2;
		double emr_weight = w / 2;
		V3f emp_sum; V3f emq_sum; V3f emr_sum;
		for (int emp_i = 0; emp_i < emp_size; emp_i++)
		{
			emp_sum = emp_sum + emp_vertex[emp_i];
		}
		for (int emq_i = 0; emq_i < emq_size; emq_i++)
		{
			emq_sum = emq_sum + emq_vertex[emq_i];
		}
		for (int emr_i = 0; emr_i < emr_size; emr_i++)
		{
			emr_sum = emr_sum + emr_vertex[emr_i];
		}
		edge_mid = emp_sum * emp_weight + emq_sum * emq_weight - emr_sum * emr_weight;
	}
	else if (check_edge(em_e))//当该边在边界块上时，用线型细分
	{
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver1]]);
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver2]]);
		int emp_size = emp_vertex.size();
		V3f emp_sum;
		for (int emp_i = 0; emp_i < emp_size; emp_i++)
		{
			emp_sum = emp_sum + emp_vertex[emp_i];
		}
		edge_mid = emp_sum * (1.0 / 2);
	}
	else
	{
		//加入p点集合
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver1]]);
		emp_vertex.push_back(mesh->vertices_[mesh->cells_[cell_num].vertex_[ver2]]);
		//加入q、r、s点集合
		std::set<FaceHandle, compare_OVM> em_e_f = mesh->neighbor_e[em_e];
		const int N = em_e_f.size();
		for (std::set<FaceHandle, compare_OVM>::iterator em_e_f_it = em_e_f.begin(); em_e_f_it != em_e_f.end(); em_e_f_it++)
		{
			//为判断将面的邻接块放入数组
			CellHandle em_e_f_cche[2];
			std::set<CellHandle, compare_OVM> em_e_f_c = mesh->neighbor_f[*em_e_f_it];
			int em_e_f_c_num = 0;
			for (std::set<CellHandle, compare_OVM>::iterator em_e_f_c_it = em_e_f_c.begin(); em_e_f_c_it != em_e_f_c.end(); em_e_f_c_it++)
			{
				em_e_f_cche[em_e_f_c_num++] = *em_e_f_c_it;
			}
			//加入q点集合
			EdgeHandle em_e_f_ope = mesh->get_opposite_edge_in_face(em_e, *em_e_f_it);
			emq_vertex.push_back(mesh->vertices_[mesh->edges_[em_e_f_ope].from_h]);
			emq_vertex.push_back(mesh->vertices_[mesh->edges_[em_e_f_ope].to_h]);
			/*std::vector<VertexHandle> em_e_f_v = mesh->faces_[*em_e_f_it].vs_f;
			for (std::vector<VertexHandle>::iterator em_e_f_v_it = em_e_f_v.begin();em_e_f_v_it != em_e_f_v.end();em_e_f_v_it++)
			{
				if ((*em_e_f_v_it != mesh->cells_[cell_num].vertex_[ver2]) && (*em_e_f_v_it != mesh->cells_[cell_num].vertex_[ver1]))
					emq_vertex.push_back(mesh->vertices_[*em_e_f_v_it]);
			}*/
			//加入r、s点集合
			std::vector<EdgeHandle> em_e_f_e = mesh->faces_[*em_e_f_it].edges_;
			//EdgeHandle em_ope = mesh->get_opposite_edge_in_face(em_e, *em_e_f_it);
			for (std::vector<EdgeHandle>::iterator em_e_f_e_it = em_e_f_e.begin(); em_e_f_e_it != em_e_f_e.end(); em_e_f_e_it++)
			{
				if ((*em_e_f_e_it != em_e) && (*em_e_f_e_it != em_e_f_ope))
				{
					int sin_size = mesh->neighbor_e[*em_e_f_e_it].size();
					//如果该邻接边不是奇异边，则加入r点
					if (!mesh->is_singular(*em_e_f_e_it))
					{
						//确定r点
						//确定这条边上的q点
						VertexHandle em_e_f_e_v;
						if (std::find(emq_vertex.begin(), emq_vertex.end(), mesh->vertices_[mesh->edges_[*em_e_f_e_it].from_h]) != emq_vertex.end())
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].from_h;
						else
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].to_h;
						std::set<FaceHandle, compare_OVM> em_e_f_e_f = mesh->neighbor_e[*em_e_f_e_it];
						for (std::set<FaceHandle, compare_OVM>::iterator em_e_f_e_f_it = em_e_f_e_f.begin(); em_e_f_e_f_it != em_e_f_e_f.end(); em_e_f_e_f_it++)
						{
							if (!mesh->is_f_in_c(*em_e_f_e_f_it, em_e_f_cche[0]) && !mesh->is_f_in_c(*em_e_f_e_f_it, em_e_f_cche[1]))
							{
								//该面为r点所在面
								//确定r点所在边
								std::vector<EdgeHandle> em_e_f_e_f_e = mesh->faces_[*em_e_f_e_f_it].edges_;
								for (std::vector<EdgeHandle>::iterator em_e_f_e_f_e_it = em_e_f_e_f_e.begin(); em_e_f_e_f_e_it != em_e_f_e_f_e.end(); em_e_f_e_f_e_it++)
								{
									if (*em_e_f_e_f_e_it == *em_e_f_e_it)
										continue;
									else if (mesh->edges_[*em_e_f_e_f_e_it].from_h == em_e_f_e_v)
										emr_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_e_f_e_it].to_h]);
									else if (mesh->edges_[*em_e_f_e_f_e_it].to_h == em_e_f_e_v)
										emr_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_e_f_e_it].from_h]);
									else
										continue;
								}
							}
						}
					}
					//如果为奇异边，则需加入s点
					if (mesh->is_singular(*em_e_f_e_it) && sin_size > 4)//度数小于4的奇异边没有s点
					{
						//确定s点
						//确定这条边上的q点
						VertexHandle em_e_f_e_v;
						if (std::find(emq_vertex.begin(), emq_vertex.end(), mesh->vertices_[mesh->edges_[*em_e_f_e_it].from_h]) != emq_vertex.end())
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].from_h;
						else
							em_e_f_e_v = mesh->edges_[*em_e_f_e_it].to_h;
						//首先取第一个位于em_e_f_cche[0]一侧的s点
						std::vector<FaceHandle> em_e_f_c_f0 = mesh->cells_[em_e_f_cche[0]].faces_;
						for (std::vector<FaceHandle>::iterator em_e_f_c_f_it = em_e_f_c_f0.begin(); em_e_f_c_f_it != em_e_f_c_f0.end(); em_e_f_c_f_it++)
						{
							if (mesh->is_e_in_f(*em_e_f_e_it, *em_e_f_c_f_it) && (*em_e_f_it != *em_e_f_c_f_it))
							{
								//该面为s点所在块的面，现在取这个块
								std::set<CellHandle, compare_OVM> em_e_f_c_f_c = mesh->neighbor_f[*em_e_f_c_f_it];
								for (std::set<CellHandle, compare_OVM>::iterator em_e_f_c_f_c_it = em_e_f_c_f_c.begin();
									em_e_f_c_f_c_it != em_e_f_c_f_c.end(); em_e_f_c_f_c_it++)
								{
									if (*em_e_f_c_f_c_it != em_e_f_cche[0])
									{
										//这个块为s点所在块，现在取s点所在面
										std::vector<FaceHandle> em_e_f_c_f_c_f = mesh->cells_[*em_e_f_c_f_c_it].faces_;
										for (std::vector<FaceHandle>::iterator em_e_f_c_f_c_f_it = em_e_f_c_f_c_f.begin();
											em_e_f_c_f_c_f_it != em_e_f_c_f_c_f.end(); em_e_f_c_f_c_f_it++)
										{
											if (mesh->is_e_in_f(*em_e_f_e_it, *em_e_f_c_f_c_f_it) && (*em_e_f_c_f_c_f_it != *em_e_f_c_f_it))
											{
												//该面为s点所在面
												//确定s点所在边
												std::vector<EdgeHandle> em_e_f_c_f_c_f_e = mesh->faces_[*em_e_f_c_f_c_f_it].edges_;
												for (std::vector<EdgeHandle>::iterator em_e_f_c_f_c_f_e_it = em_e_f_c_f_c_f_e.begin();
													em_e_f_c_f_c_f_e_it != em_e_f_c_f_c_f_e.end(); em_e_f_c_f_c_f_e_it++)
												{
													if (*em_e_f_c_f_c_f_e_it == *em_e_f_e_it)
														continue;
													else if (mesh->edges_[*em_e_f_c_f_c_f_e_it].from_h == em_e_f_e_v)
														ems_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_c_f_c_f_e_it].to_h]);
													else if (mesh->edges_[*em_e_f_c_f_c_f_e_it].to_h == em_e_f_e_v)
														ems_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_c_f_c_f_e_it].from_h]);
													else
														continue;
												}
											}
										}
									}
								}
							}
						}
						//再取第二个位于em_e_f_cche[1]一侧的s点
						std::vector<FaceHandle> em_e_f_c_f1 = mesh->cells_[em_e_f_cche[1]].faces_;
						for (std::vector<FaceHandle>::iterator em_e_f_c_f_it = em_e_f_c_f1.begin(); em_e_f_c_f_it != em_e_f_c_f1.end(); em_e_f_c_f_it++)
						{
							if (mesh->is_e_in_f(*em_e_f_e_it, *em_e_f_c_f_it) && (*em_e_f_it != *em_e_f_c_f_it))
							{
								//该面为s点所在块的面，现在取这个块
								std::set<CellHandle, compare_OVM> em_e_f_c_f_c = mesh->neighbor_f[*em_e_f_c_f_it];
								for (std::set<CellHandle, compare_OVM>::iterator em_e_f_c_f_c_it = em_e_f_c_f_c.begin();
									em_e_f_c_f_c_it != em_e_f_c_f_c.end(); em_e_f_c_f_c_it++)
								{
									if (*em_e_f_c_f_c_it != em_e_f_cche[1])
									{
										//这个块为s点所在块，现在取s点所在面
										std::vector<FaceHandle> em_e_f_c_f_c_f = mesh->cells_[*em_e_f_c_f_c_it].faces_;
										for (std::vector<FaceHandle>::iterator em_e_f_c_f_c_f_it = em_e_f_c_f_c_f.begin();
											em_e_f_c_f_c_f_it != em_e_f_c_f_c_f.end(); em_e_f_c_f_c_f_it++)
										{
											if (mesh->is_e_in_f(*em_e_f_e_it, *em_e_f_c_f_c_f_it) && (*em_e_f_c_f_c_f_it != *em_e_f_c_f_it))
											{
												//该面为s点所在面
												//确定s点所在边
												std::vector<EdgeHandle> em_e_f_c_f_c_f_e = mesh->faces_[*em_e_f_c_f_c_f_it].edges_;
												for (std::vector<EdgeHandle>::iterator em_e_f_c_f_c_f_e_it = em_e_f_c_f_c_f_e.begin();
													em_e_f_c_f_c_f_e_it != em_e_f_c_f_c_f_e.end(); em_e_f_c_f_c_f_e_it++)
												{
													if (*em_e_f_c_f_c_f_e_it == *em_e_f_e_it)
														continue;
													else if (mesh->edges_[*em_e_f_c_f_c_f_e_it].from_h == em_e_f_e_v)
														ems_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_c_f_c_f_e_it].to_h]);
													else if (mesh->edges_[*em_e_f_c_f_c_f_e_it].to_h == em_e_f_e_v)
														ems_vertex.push_back(mesh->vertices_[mesh->edges_[*em_e_f_c_f_c_f_e_it].from_h]);
													else
														continue;
												}
											}
										}
									}
								}
							}
						}
					}
				}
			}
		}
		//计算边的中心点
		int emp_size = emp_vertex.size();
		int emq_size = emq_vertex.size();
		int emr_size = emr_vertex.size();
		int ems_size = ems_vertex.size();

		double emp_weight = 1.0 / 2;
		double emq_weight = w / N;
		/*double emr_weight = w / N;
		double ems_weight = w / (2 * N);*/
		double emr_weight = w / ((emr_size + ems_size / 2) / 2);
		double ems_weight = emr_weight / 2;

		/*double emp_weight, emq_weight, emr_weight, ems_weight;
		emp_weight = 1.0 / 2;
		emq_weight = w / N;
		if((emr_size + ems_size / 2) != 2*N)
			emr_weight = w / N;
		else
			emr_weight = w / ((emr_size + ems_size / 2) / 2);
		ems_weight = emr_weight / 2;*/

		V3f emp_sum; V3f emq_sum; V3f emr_sum; V3f ems_sum;
		for (int emp_i = 0; emp_i < emp_size; emp_i++)
		{
			emp_sum = emp_sum + emp_vertex[emp_i];
		}
		for (int emq_i = 0; emq_i < emq_size; emq_i++)
		{
			emq_sum = emq_sum + emq_vertex[emq_i];
		}
		for (int emr_i = 0; emr_i < emr_size; emr_i++)
		{
			emr_sum = emr_sum + emr_vertex[emr_i];
		}
		for (int ems_i = 0; ems_i < ems_size; ems_i++)
		{
			ems_sum = ems_sum + ems_vertex[ems_i];
		}
		edge_mid = emp_sum * emp_weight + emq_sum * emq_weight - emr_sum * emr_weight - ems_sum * ems_weight;
	}
	return edge_mid;
}

bool Subdivision::check_face(FaceHandle _f)
{
	std::set<CellHandle, compare_OVM> fm_c = mesh->neighbor_f[_f];
	bool check = false;
	for (std::set<CellHandle, compare_OVM>::iterator fm_c_it = fm_c.begin(); fm_c_it != fm_c.end(); fm_c_it++)
	{
		if (mesh->is_bdy(*fm_c_it))
			check = true;
		break;
	}
	return check;
}

bool Subdivision::check_edge(EdgeHandle _e)
{
	bool check = false;
	std::set<FaceHandle, compare_OVM> em_f = mesh->neighbor_e[_e];
	for (std::set<FaceHandle, compare_OVM>::iterator em_f_it = em_f.begin(); em_f_it != em_f.end(); em_f_it++)
	{
		if (check == true) break;
		std::vector<EdgeHandle> em_f_e = mesh->faces_[*em_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_f_e_it = em_f_e.begin(); em_f_e_it != em_f_e.end(); em_f_e_it++)
		{
			if (mesh->is_bdy(*em_f_e_it))
			{
				check = true;
				break;
			}
		}
	}
	return check;
}

//-----CCSubdivision-----//
HexV3fMesh& CCSubdivision::CCSubdivising()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<double> material_distribution1(8 * cell_num, 1);
	for (uint64_t i = 0; i < cell_num; i++)
	{

		std::cout << "Subdivising: " << i << "/" << cell_num << std::endl;
		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_new_vertex(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_edge_mid(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_edge_mid(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_edge_mid(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_edge_mid(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_edge_mid(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_edge_mid(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_edge_mid(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_edge_mid(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_edge_mid(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_edge_mid(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_edge_mid(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_edge_mid(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		if (material_distribution[i] != 1)
		{
			material_distribution1[newmesh->add_cell(sub_cell[0])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[1])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[2])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[3])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[4])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[5])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[6])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[7])] = 0.0001;
		}
		else
		{
			newmesh->add_cell(sub_cell[0]);
			newmesh->add_cell(sub_cell[1]);
			newmesh->add_cell(sub_cell[2]);
			newmesh->add_cell(sub_cell[3]);
			newmesh->add_cell(sub_cell[4]);
			newmesh->add_cell(sub_cell[5]);
			newmesh->add_cell(sub_cell[6]);
			newmesh->add_cell(sub_cell[7]);
		}
	}
	material_distribution = material_distribution1;
	return *newmesh;
}

HexV3fMesh& CCSubdivision::CCSubdivising_1()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_new_vertex_1(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_edge_mid_1(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_edge_mid_1(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_edge_mid_1(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_edge_mid_1(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_edge_mid_1(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_edge_mid_1(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_edge_mid_1(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_edge_mid_1(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_edge_mid_1(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_edge_mid_1(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_edge_mid_1(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_edge_mid_1(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

V3f CCSubdivision::cal_cell_mid(const int cell_num)
{
	V3f cell_mid;
	std::vector<VertexHandle> c_v = mesh->cells_[cell_num].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it];
	}
	cell_mid = cell_mid / 8;
	return cell_mid;
}

V3f CCSubdivision::cal_cell_midv(CellHandle cell_id)
{
	V3f cell_mid;
	std::vector<VertexHandle> c_v = mesh->cells_[cell_id].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it];
	}
	cell_mid = cell_mid / 8;
	return cell_mid;
}

V3f CCSubdivision::cal_face_midv(FaceHandle face_id)
{
	V3f face_mid;
	std::vector<VertexHandle> f_v = mesh->faces_[face_id].vs_f;
	for (std::vector<VertexHandle>::iterator f_v_it = f_v.begin(); f_v_it != f_v.end(); f_v_it++)
	{
		face_mid = face_mid + mesh->vertices_[*f_v_it];
	}
	face_mid = face_mid / 4;
	return face_mid;
}

V3f CCSubdivision::cal_face_mid(const int cell_num, const int face_num)
{
	V3f face_mid;
	V3f face_sum;
	FaceHandle face = mesh->cells_[cell_num].faces_[face_num];//确定该面
	if (mesh->is_bdy(face))
	{
		face_mid = cal_face_midv(face);
	}
	else
	{
		V3f face_ver = cal_face_midv(face);//面点
		face_sum = face_sum + face_ver * 2;
		//计算邻接块点
		std::set<CellHandle, compare_OVM> f_c = mesh->neighbor_f[face];
		for (std::set<CellHandle, compare_OVM>::iterator f_c_it = f_c.begin(); f_c_it != f_c.end(); f_c_it++)
		{
			face_sum = face_sum + cal_cell_midv(*f_c_it);
		}
		face_mid = face_sum / 4;
	}
	return face_mid;
}

V3f CCSubdivision::cal_face_mid(FaceHandle face_id)
{
	V3f face_mid;
	V3f face_sum;
	if (mesh->is_bdy(face_id))
	{
		face_mid = cal_face_midv(face_id);
	}
	else
	{
		V3f face_ver = cal_face_midv(face_id);//面点
		face_sum = face_sum + face_ver * 2;
		//计算邻接块点
		std::set<CellHandle, compare_OVM> f_c = mesh->neighbor_f[face_id];
		for (std::set<CellHandle, compare_OVM>::iterator f_c_it = f_c.begin(); f_c_it != f_c.end(); f_c_it++)
		{
			face_sum = face_sum + cal_cell_midv(*f_c_it);
		}
		face_mid = face_sum / 4;
	}
	return face_mid;
}

V3f CCSubdivision::cal_edge_midv(EdgeHandle edge_id)
{
	V3f edge_mid;
	edge_mid = edge_mid + mesh->vertices_[mesh->edges_[edge_id].from_h];
	edge_mid = edge_mid + mesh->vertices_[mesh->edges_[edge_id].to_h];
	edge_mid = edge_mid / 2;
	return edge_mid;
}

V3f CCSubdivision::cal_edge_mid(const int cell_num, const int ver1, const int ver2)
{
	V3f edge_mid;
	V3f edge_sum;
	//首先确定该边
	EdgeHandle edge;
	bool is_find = false;
	std::vector<FaceHandle> em_c_f = mesh->cells_[cell_num].faces_;
	for (std::vector<FaceHandle>::iterator em_c_f_it = em_c_f.begin(); em_c_f_it != em_c_f.end(); em_c_f_it++)
	{
		if (is_find == true)
			break;
		std::vector<EdgeHandle> em_c_f_e = mesh->faces_[*em_c_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_c_f_e_it = em_c_f_e.begin(); em_c_f_e_it != em_c_f_e.end(); em_c_f_e_it++)
		{
			if (((mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver2]))
				|| ((mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver2])))
			{
				edge = *em_c_f_e_it;
				is_find = true;
				break;
			}
		}
	}
	if (mesh->is_bdy(edge))
	{
		edge_sum = edge_sum + cal_edge_midv(edge) * 2;
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge];
		int e_f_n = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			if (mesh->is_bdy(*e_f_it))
			{
				e_f_n++;
				edge_sum = edge_sum + cal_face_midv(*e_f_it);
			}
		}
		edge_mid = edge_sum / (2 + e_f_n);
	}
	else
	{
		std::set<CellHandle, compare_OVM> e_c;//存边的邻接块
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge];
		const int n = e_f.size();
		edge_sum = edge_sum + cal_edge_midv(edge) * (n - 3);
		//edge_sum = edge_sum + cal_edge_midv(edge);
		V3f e_f_mid;//计算邻接平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			e_f_mid = e_f_mid + cal_face_midv(*e_f_it);
			std::set<CellHandle, compare_OVM> e_f_c = mesh->neighbor_f[*e_f_it];//根据边的邻接面的邻接块加入边的邻接块
			for (std::set<CellHandle, compare_OVM>::iterator e_f_c_it = e_f_c.begin(); e_f_c_it != e_f_c.end(); e_f_c_it++)
			{
				e_c.insert(*e_f_c_it);
			}
		}
		e_f_mid = e_f_mid / n;
		V3f e_c_mid;//计算边的邻接平均块点
		const int e_c_n = e_c.size();
		for (std::set<CellHandle, compare_OVM>::iterator e_c_it = e_c.begin(); e_c_it != e_c.end(); e_c_it++)
		{
			e_c_mid = e_c_mid + cal_cell_midv(*e_c_it);
		}
		e_c_mid = e_c_mid / e_c_n;
		edge_sum = edge_sum + e_f_mid * 2 + e_c_mid;
		edge_mid = edge_sum / n;
		//edge_mid = edge_sum / 4;
	}
	return edge_mid;
}

V3f CCSubdivision::cal_edge_mid_1(const int cell_num, const int ver1, const int ver2)
{
	V3f edge_mid;
	V3f edge_sum;
	//首先确定该边
	EdgeHandle edge;
	bool is_find = false;
	std::vector<FaceHandle> em_c_f = mesh->cells_[cell_num].faces_;
	for (std::vector<FaceHandle>::iterator em_c_f_it = em_c_f.begin(); em_c_f_it != em_c_f.end(); em_c_f_it++)
	{
		if (is_find == true)
			break;
		std::vector<EdgeHandle> em_c_f_e = mesh->faces_[*em_c_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_c_f_e_it = em_c_f_e.begin(); em_c_f_e_it != em_c_f_e.end(); em_c_f_e_it++)
		{
			if (((mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver2]))
				|| ((mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver2])))
			{
				edge = *em_c_f_e_it;
				is_find = true;
				break;
			}
		}
	}
	if (mesh->is_bdy(edge))
	{
		std::set<FaceHandle, compare_OVM> e_f_temp = mesh->neighbor_e[edge];
		int ef_temp = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f_temp.begin(); e_f_it != e_f_temp.end(); e_f_it++)
		{
			if (mesh->is_bdy(*e_f_it))
				ef_temp++;
		}
		if (ef_temp == 4 /*|| ef_temp > 1*/)
		{
			edge_sum = edge_sum + cal_edge_midv(edge) * 2;
			std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge];
			int e_f_n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
			{
				if (mesh->is_bdy(*e_f_it) && mesh->is_f_in_c(*e_f_it, CellHandle(cell_num)))
				{
					e_f_n++;
					edge_sum = edge_sum + cal_face_midv(*e_f_it);
				}
			}
			edge_mid = edge_sum / (2 + e_f_n);
		}
		else
		{
			edge_sum = edge_sum + cal_edge_midv(edge) * 2;
			std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge];
			int e_f_n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
			{
				if (mesh->is_bdy(*e_f_it))
				{
					e_f_n++;
					edge_sum = edge_sum + cal_face_midv(*e_f_it);
				}
			}
			edge_mid = edge_sum / (2 + e_f_n);
		}

	}
	else
	{
		std::set<CellHandle, compare_OVM> e_c;//存边的邻接块
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge];
		const int n = e_f.size();
		edge_sum = edge_sum + cal_edge_midv(edge) * (n - 3);
		//edge_sum = edge_sum + cal_edge_midv(edge);
		V3f e_f_mid;//计算邻接平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			e_f_mid = e_f_mid + cal_face_midv(*e_f_it);
			std::set<CellHandle, compare_OVM> e_f_c = mesh->neighbor_f[*e_f_it];//根据边的邻接面的邻接块加入边的邻接块
			for (std::set<CellHandle, compare_OVM>::iterator e_f_c_it = e_f_c.begin(); e_f_c_it != e_f_c.end(); e_f_c_it++)
			{
				e_c.insert(*e_f_c_it);
			}
		}
		e_f_mid = e_f_mid / n;
		V3f e_c_mid;//计算边的邻接平均块点
		const int e_c_n = e_c.size();
		for (std::set<CellHandle, compare_OVM>::iterator e_c_it = e_c.begin(); e_c_it != e_c.end(); e_c_it++)
		{
			e_c_mid = e_c_mid + cal_cell_midv(*e_c_it);
		}
		e_c_mid = e_c_mid / e_c_n;
		edge_sum = edge_sum + e_f_mid * 2 + e_c_mid;
		edge_mid = edge_sum / n;
		//edge_mid = edge_sum / 4;
	}
	return edge_mid;
}

V3f CCSubdivision::cal_edge_mid(EdgeHandle edge_id)
{
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge_id))
	{
		edge_sum = edge_sum + cal_edge_midv(edge_id) * 2;
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge_id];
		int e_f_n = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			if (mesh->is_bdy(*e_f_it))
			{
				e_f_n++;
				edge_sum = edge_sum + cal_face_midv(*e_f_it);
			}
		}
		edge_mid = edge_sum / (2 + e_f_n);
	}
	else
	{
		std::set<CellHandle, compare_OVM> e_c;//存边的邻接块
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[edge_id];
		const int n = e_f.size();
		edge_sum = edge_sum + cal_edge_midv(edge_id) * (n - 3);
		//edge_sum = edge_sum + cal_edge_midv(edge_id);
		V3f e_f_mid;//计算邻接平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			e_f_mid = e_f_mid + cal_face_midv(*e_f_it);
			std::set<CellHandle, compare_OVM> e_f_c = mesh->neighbor_f[*e_f_it];//根据边的邻接面的邻接块加入边的邻接块
			for (std::set<CellHandle, compare_OVM>::iterator e_f_c_it = e_f_c.begin(); e_f_c_it != e_f_c.end(); e_f_c_it++)
			{
				e_c.insert(*e_f_c_it);
			}
		}
		e_f_mid = e_f_mid / n;
		V3f e_c_mid;//计算边的邻接平均块点
		const int e_c_n = e_c.size();
		for (std::set<CellHandle, compare_OVM>::iterator e_c_it = e_c.begin(); e_c_it != e_c.end(); e_c_it++)
		{
			e_c_mid = e_c_mid + cal_cell_midv(*e_c_it);
		}
		e_c_mid = e_c_mid / e_c_n;
		edge_sum = edge_sum + e_f_mid * 2 + e_c_mid;
		edge_mid = edge_sum / n;
		//edge_mid = edge_sum / 4;
	}
	return edge_mid;
}

V3f CCSubdivision::cal_new_vertex(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int v_e_n = 0;
		V3f v_e_mid;//计算点邻接边的平均边点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				v_e_n++;
				v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
			}
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = vertex_sum + v_e_mid * 2;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		int v_f_n = 0;
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				v_f_n++;
				v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
			}
		}
		v_f_mid = v_f_mid / v_f_n;
		vertex_sum = vertex_sum + v_f_mid;
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (v_e_n - 3);
		new_vertex = vertex_sum / v_e_n;
	}
	else
	{
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		const int v_c_n = v_c.size();
		V3f v_c_mid;//计算点邻接块的平均块点
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		}
		v_c_mid = v_c_mid / v_c_n;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		const int v_f_n = v_f.size();
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		}
		v_f_mid = v_f_mid / v_f_n;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		const int v_e_n = v_e.size();
		V3f v_e_mid;//计算点邻接面的平均面点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = v_c_mid + v_f_mid * 3 + v_e_mid * 3 + mesh->vertices_[vertex];
		new_vertex = vertex_sum / 8;
	}
	return new_vertex;
}

V3f CCSubdivision::cal_new_vertex_1(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		std::set<EdgeHandle, compare_OVM> v_e_temp = mesh->neighbor_v[vertex];
		int bdy_edge = 0;
		int not_bdy_edge = 0;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e_temp.begin(); v_e_it != v_e_temp.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
				bdy_edge++;
			else
				not_bdy_edge++;
		}
		if (bdy_edge > 5 || not_bdy_edge < 1)
		{
			std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
			int v_e_n = 0;
			V3f v_e_mid;//计算点邻接边的平均边点
			std::set<CellHandle> c_c_temp = mesh->neighbor_c_c(CellHandle(cell_num));
			c_c_temp.insert(CellHandle(cell_num));
			for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
			{
				for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin(); c_c_temp_it != c_c_temp.end(); c_c_temp_it++)
				{
					if (mesh->is_bdy(*v_e_it) && mesh->is_e_in_c(*v_e_it, *c_c_temp_it))
					{
						v_e_n++;
						v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
						break;
					}
				}

			}
			v_e_mid = v_e_mid / v_e_n;
			vertex_sum = vertex_sum + v_e_mid * 2;
			std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
			int v_f_n = 0;
			V3f v_f_mid;//计算点邻接面的平均面点
			for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
			{
				for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin(); c_c_temp_it != c_c_temp.end(); c_c_temp_it++)
				{
					if (mesh->is_bdy(*v_f_it) && mesh->is_f_in_c(*v_f_it, *c_c_temp_it))
					{
						v_f_n++;
						v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
						break;
					}
				}

			}
			v_f_mid = v_f_mid / v_f_n;
			vertex_sum = vertex_sum + v_f_mid;
			vertex_sum = vertex_sum + mesh->vertices_[vertex] * (v_e_n - 3);
			new_vertex = vertex_sum / v_e_n;
		}
		else
		{
			std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
			int v_e_n = 0;
			V3f v_e_mid;//计算点邻接边的平均边点
			for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
			{
				if (mesh->is_bdy(*v_e_it))
				{
					v_e_n++;
					v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
				}
			}
			v_e_mid = v_e_mid / v_e_n;
			vertex_sum = vertex_sum + v_e_mid * 2;
			std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
			int v_f_n = 0;
			V3f v_f_mid;//计算点邻接面的平均面点
			for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
			{
				if (mesh->is_bdy(*v_f_it))
				{
					v_f_n++;
					v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
				}
			}
			v_f_mid = v_f_mid / v_f_n;
			vertex_sum = vertex_sum + v_f_mid;
			vertex_sum = vertex_sum + mesh->vertices_[vertex] * (v_e_n - 3);
			new_vertex = vertex_sum / v_e_n;
		}
	}
	else
	{
		//if (mesh->is_singular(vertex))
		//{
		//	int n = mesh->neighbor_v[vertex].size();//度数
		//	std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		//	const int v_c_n = v_c.size();
		//	V3f v_c_mid;//计算点邻接块的平均块点
		//	for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin();v_c_it != v_c.end();v_c_it++)
		//	{
		//		v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		//	}
		//	v_c_mid = v_c_mid / ((n - 2) * 16);
		//	std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		//	V3f v_f_mid;//计算点邻接面的平均面点
		//	for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin();v_f_it != v_f.end();v_f_it++)
		//	{
		//		v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		//	}
		//	v_f_mid = v_f_mid / ((n - 2) * 8);
		//	std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		//	V3f v_e_mid;//计算点邻接面的平均面点
		//	for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin();v_e_it != v_e.end();v_e_it++)
		//	{
		//		v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		//	}
		//	v_e_mid = v_e_mid / ((n - 2) * 4);
		//	new_vertex = new_vertex + v_c_mid + v_f_mid + v_e_mid + mesh->vertices_[vertex] / 8;
		//}
		//else
		//{
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		const int v_c_n = v_c.size();
		V3f v_c_mid;//计算点邻接块的平均块点
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		}
		v_c_mid = v_c_mid / v_c_n;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		const int v_f_n = v_f.size();
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		}
		v_f_mid = v_f_mid / v_f_n;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		const int v_e_n = v_e.size();
		V3f v_e_mid;//计算点邻接面的平均面点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = v_c_mid + v_f_mid * 3 + v_e_mid * 3 + mesh->vertices_[vertex];
		new_vertex = vertex_sum / 8;
		//}
	}
	return new_vertex;
}

V3f CCSubdivision::cal_new_vertex(VertexHandle ver_id)
{
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(ver_id))
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[ver_id];
		int v_e_n = 0;
		V3f v_e_mid;//计算点邻接边的平均边点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				v_e_n++;
				v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
			}
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = vertex_sum + v_e_mid * 2;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[ver_id];
		int v_f_n = 0;
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				v_f_n++;
				v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
			}
		}
		v_f_mid = v_f_mid / v_f_n;
		vertex_sum = vertex_sum + v_f_mid;
		vertex_sum = vertex_sum + mesh->vertices_[ver_id] * (v_e_n - 3);
		new_vertex = vertex_sum / v_e_n;
	}
	else
	{
		//if (mesh->is_singular(ver_id))
		//{
		//	int n = mesh->neighbor_v[ver_id].size();//度数
		//	std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[ver_id];
		//	const int v_c_n = v_c.size();
		//	V3f v_c_mid;//计算点邻接块的平均块点
		//	for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin();v_c_it != v_c.end();v_c_it++)
		//	{
		//		v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		//	}
		//	v_c_mid = v_c_mid / ((n - 2) * 16);
		//	std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[ver_id];
		//	V3f v_f_mid;//计算点邻接面的平均面点
		//	for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin();v_f_it != v_f.end();v_f_it++)
		//	{
		//		v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		//	}
		//	v_f_mid = v_f_mid / ((n - 2) * 8);
		//	std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[ver_id];
		//	V3f v_e_mid;//计算点邻接面的平均面点
		//	for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin();v_e_it != v_e.end();v_e_it++)
		//	{
		//		v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		//	}
		//	v_e_mid = v_e_mid / ((n - 2) * 4);
		//	new_vertex = new_vertex + v_c_mid + v_f_mid + v_e_mid + mesh->vertices_[ver_id] / 8;
		//}
		//else
		//{
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[ver_id];
		const int v_c_n = v_c.size();
		V3f v_c_mid;//计算点邻接块的平均块点
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		}
		v_c_mid = v_c_mid / v_c_n;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[ver_id];
		const int v_f_n = v_f.size();
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		}
		v_f_mid = v_f_mid / v_f_n;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[ver_id];
		const int v_e_n = v_e.size();
		V3f v_e_mid;//计算点邻接面的平均面点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = v_c_mid + v_f_mid * 3 + v_e_mid * 3 + mesh->vertices_[ver_id];
		new_vertex = vertex_sum / 8;
		//}
	}
	return new_vertex;
}

HexV3fMesh& CCSubdivision::CCSubWithCrease()
{
	uint64_t cell_num = mesh->cells_.size();
	std::vector<bool> new_v_cs;//折痕点
	new_v_cs.resize(18 * cell_num + 9, false);
	std::vector<double> material_distribution1(8 * cell_num, 1);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		if (i % 10 == 0) {
			std::cout << "Subdivising: " << i << "/" << cell_num << std::endl;
		}

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_new_csv(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
			if (v_cs[mesh->cells_[i].vertex_[mesh_vex_num]])
			{
				new_v_cs[temp_vex[mesh_vex_num]] = true;
			}
		}

		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_es_mid(i, 0, 1);
		EdgeHandle edm1 = check_edge(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);
		if (is_cs(edm1))
		{
			new_v_cs[temp_vex[8]] = true;
		}

		//第二条边的中点
		V3f edge_mid_2 = cal_es_mid(i, 1, 5);
		EdgeHandle edm2 = check_edge(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);
		if (is_cs(edm2))
		{
			new_v_cs[temp_vex[9]] = true;
		}

		//第三条边的中点
		V3f edge_mid_3 = cal_es_mid(i, 5, 4);
		EdgeHandle edm3 = check_edge(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);
		if (is_cs(edm3))
		{
			new_v_cs[temp_vex[10]] = true;
		}

		//第四条边的中点
		V3f edge_mid_4 = cal_es_mid(i, 4, 0);
		EdgeHandle edm4 = check_edge(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);
		if (is_cs(edm4))
		{
			new_v_cs[temp_vex[11]] = true;
		}

		//第五条边的中点
		V3f edge_mid_5 = cal_es_mid(i, 3, 2);
		EdgeHandle edm5 = check_edge(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);
		if (is_cs(edm5))
		{
			new_v_cs[temp_vex[12]] = true;
		}

		//第六条边的中点
		V3f edge_mid_6 = cal_es_mid(i, 3, 7);
		EdgeHandle edm6 = check_edge(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);
		if (is_cs(edm6))
		{
			new_v_cs[temp_vex[13]] = true;
		}

		//第七条边的中点
		V3f edge_mid_7 = cal_es_mid(i, 6, 7);
		EdgeHandle edm7 = check_edge(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);
		if (is_cs(edm7))
		{
			new_v_cs[temp_vex[14]] = true;
		}

		//第八条边的中点
		V3f edge_mid_8 = cal_es_mid(i, 2, 6);
		EdgeHandle edm8 = check_edge(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);
		if (is_cs(edm8))
		{
			new_v_cs[temp_vex[15]] = true;
		}

		//第九条边的中点
		V3f edge_mid_9 = cal_es_mid(i, 1, 2);
		EdgeHandle edm9 = check_edge(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);
		if (is_cs(edm9))
		{
			new_v_cs[temp_vex[16]] = true;
		}

		//第十条边的中点
		V3f edge_mid_10 = cal_es_mid(i, 5, 6);
		EdgeHandle edm10 = check_edge(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);
		if (is_cs(edm10))
		{
			new_v_cs[temp_vex[17]] = true;
		}

		//第十一条边的中点
		V3f edge_mid_11 = cal_es_mid(i, 0, 3);
		EdgeHandle edm11 = check_edge(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);
		if (is_cs(edm11))
		{
			new_v_cs[temp_vex[18]] = true;
		}

		//第十二条边的中点
		V3f edge_mid_12 = cal_es_mid(i, 4, 7);
		EdgeHandle edm12 = check_edge(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		if (is_cs(edm12))
		{
			new_v_cs[temp_vex[19]] = true;
		}
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		if (material_distribution[i] != 1)
		{
			material_distribution1[newmesh->add_cell(sub_cell[0])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[1])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[2])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[3])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[4])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[5])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[6])] = 0.0001;
			material_distribution1[newmesh->add_cell(sub_cell[7])] = 0.0001;
		}
		else
		{
			newmesh->add_cell(sub_cell[0]);
			newmesh->add_cell(sub_cell[1]);
			newmesh->add_cell(sub_cell[2]);
			newmesh->add_cell(sub_cell[3]);
			newmesh->add_cell(sub_cell[4]);
			newmesh->add_cell(sub_cell[5]);
			newmesh->add_cell(sub_cell[6]);
			newmesh->add_cell(sub_cell[7]);
		}
	}
	v_cs = new_v_cs;
	return *newmesh;
}

V3f CCSubdivision::cal_new_csv(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	if (is_cs(vertex))
	{
		new_vertex = new_vertex + mesh->vertices_[vertex] * (3.0 / 4);
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (is_cs(*v_e_it))
			{
				if (mesh->vertices_[mesh->edges_[*v_e_it].from_h] == mesh->vertices_[vertex])
				{
					//if (is_cs(mesh->edges_[*v_e_it].to_h))
					//{
					new_vertex = new_vertex + mesh->vertices_[mesh->edges_[*v_e_it].to_h] * (1.0 / 8);
					//}
				}
				else
				{
					//if (is_cs(mesh->edges_[*v_e_it].from_h))
					//{
					new_vertex = new_vertex + mesh->vertices_[mesh->edges_[*v_e_it].from_h] * (1.0 / 8);
					//}
				}
			}
		}
	}
	else
	{
		new_vertex = cal_new_vertex(cell_num, ver_num);
	}
	return new_vertex;
}

V3f CCSubdivision::cal_es_mid(const int cell_num, const int ver1, const int ver2)
{
	V3f edge_mid;
	EdgeHandle edge;
	bool is_find = false;
	std::vector<FaceHandle> em_c_f = mesh->cells_[cell_num].faces_;
	for (std::vector<FaceHandle>::iterator em_c_f_it = em_c_f.begin(); em_c_f_it != em_c_f.end(); em_c_f_it++)
	{
		if (is_find == true)
			break;
		std::vector<EdgeHandle> em_c_f_e = mesh->faces_[*em_c_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_c_f_e_it = em_c_f_e.begin(); em_c_f_e_it != em_c_f_e.end(); em_c_f_e_it++)
		{
			if (((mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver2]))
				|| ((mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver2])))
			{
				edge = *em_c_f_e_it;
				is_find = true;
				break;
			}
		}
	}
	if (is_cs(edge))
		edge_mid = cal_edge_midv(edge);
	else
		edge_mid = cal_edge_mid(cell_num, ver1, ver2);
	return edge_mid;
}

EdgeHandle CCSubdivision::check_edge(const int cell_num, const int ver1, const int ver2)
{
	EdgeHandle edge;
	bool is_find = false;
	std::vector<FaceHandle> em_c_f = mesh->cells_[cell_num].faces_;
	for (std::vector<FaceHandle>::iterator em_c_f_it = em_c_f.begin(); em_c_f_it != em_c_f.end(); em_c_f_it++)
	{
		if (is_find == true)
			break;
		std::vector<EdgeHandle> em_c_f_e = mesh->faces_[*em_c_f_it].edges_;
		for (std::vector<EdgeHandle>::iterator em_c_f_e_it = em_c_f_e.begin(); em_c_f_e_it != em_c_f_e.end(); em_c_f_e_it++)
		{
			if (((mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver2]))
				|| ((mesh->edges_[*em_c_f_e_it].to_h == mesh->cells_[cell_num].vertex_[ver1]) && (mesh->edges_[*em_c_f_e_it].from_h == mesh->cells_[cell_num].vertex_[ver2])))
			{
				edge = *em_c_f_e_it;
				is_find = true;
				break;
			}
		}
	}
	return edge;
}

HexV3fMesh& CCSubdivision::FirstCCSubdivising()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_new_vertex(i, mesh_vex_num);
			if (mesh_vex_num == 5)
				new_v.bezierid = 123;
			else if (mesh_vex_num == 1)
				new_v.bezierid = 23;
			else if (mesh_vex_num == 6)
				new_v.bezierid = 13;
			else if (mesh_vex_num == 4)
				new_v.bezierid = 12;
			else
				new_v.bezierid = -1;
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_edge_mid(i, 0, 1);
		edge_mid_1.bezierid = 2;
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_edge_mid(i, 1, 5);
		edge_mid_2.bezierid = 23;
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_edge_mid(i, 5, 4);
		edge_mid_3.bezierid = 12;
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_edge_mid(i, 4, 0);
		edge_mid_4.bezierid = 2;
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_edge_mid(i, 3, 2);
		edge_mid_5.bezierid = -1;
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_edge_mid(i, 3, 7);
		edge_mid_6.bezierid = -1;
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_edge_mid(i, 6, 7);
		edge_mid_7.bezierid = 1;
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_edge_mid(i, 2, 6);
		edge_mid_8.bezierid = 3;
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_edge_mid(i, 1, 2);
		edge_mid_9.bezierid = 3;
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_edge_mid(i, 5, 6);
		edge_mid_10.bezierid = 13;
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_edge_mid(i, 0, 3);
		edge_mid_11.bezierid = -1;
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_edge_mid(i, 4, 7);
		edge_mid_12.bezierid = 1;
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			if (mesh_face_num == 0)
				face_mid.bezierid = 2;
			else if (mesh_face_num == 2)
				face_mid.bezierid = 3;
			else if (mesh_face_num == 4)
				face_mid.bezierid = 1;
			else
				face_mid.bezierid = -1;
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

HexV3fMesh& CCSubdivision::SecondCCSubdivising()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_new_vertex(i, mesh_vex_num);
			new_v.bezierid = mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]].bezierid;
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_edge_mid(i, 0, 1);
		int id1 = mesh->vertices_[mesh->edges_[mesh->cells_[i].edges_[1]].from_h].bezierid;
		int id2 = mesh->vertices_[mesh->edges_[mesh->cells_[i].edges_[1]].to_h].bezierid;
		if (id1 == -1 || id2 == -1)
			edge_mid_1.bezierid = -1;
		//else if()
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_edge_mid(i, 1, 5);
		edge_mid_2.bezierid = 23;
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_edge_mid(i, 5, 4);
		edge_mid_3.bezierid = 12;
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_edge_mid(i, 4, 0);
		edge_mid_4.bezierid = 2;
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_edge_mid(i, 3, 2);
		edge_mid_5.bezierid = -1;
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_edge_mid(i, 3, 7);
		edge_mid_6.bezierid = -1;
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_edge_mid(i, 6, 7);
		edge_mid_7.bezierid = 1;
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_edge_mid(i, 2, 6);
		edge_mid_8.bezierid = 3;
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_edge_mid(i, 1, 2);
		edge_mid_9.bezierid = 3;
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_edge_mid(i, 5, 6);
		edge_mid_10.bezierid = 13;
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_edge_mid(i, 0, 3);
		edge_mid_11.bezierid = -1;
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_edge_mid(i, 4, 7);
		edge_mid_12.bezierid = 1;
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_mid(i, mesh_face_num);
			if (mesh_face_num == 0)
				face_mid.bezierid = 2;
			else if (mesh_face_num == 2)
				face_mid.bezierid = 3;
			else if (mesh_face_num == 4)
				face_mid.bezierid = 1;
			else
				face_mid.bezierid = -1;
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

//-----InterpolateCC-----//
/*-----------Push-Back法-----------*/
HexV3fMesh& InterpolateCC::CCInterpolatingPB()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_CIPB_nv(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(firstver[mesh->cells_[i].vertex_[mesh_vex_num]]);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_CIPB_ev(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_CIPB_ev(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_CIPB_ev(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_CIPB_ev(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_CIPB_ev(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_CIPB_ev(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_CIPB_ev(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_CIPB_ev(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_CIPB_ev(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_CIPB_ev(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_CIPB_ev(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_CIPB_ev(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_CIPB_fv(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_CIPB_cv(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

HexV3fMesh& InterpolateCC::CCInterpolatingPB1()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			VertexHandle ver = mesh->cells_[i].vertex_[mesh_vex_num];
			V3f new_v = cal_CIPB_fnv(ver);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		EdgeHandle e1 = check_edge(i, 0, 1);
		V3f edge_mid_1 = ev[e1];
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		EdgeHandle e2 = check_edge(i, 1, 5);
		V3f edge_mid_2 = ev[e2];
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		EdgeHandle e3 = check_edge(i, 5, 4);
		V3f edge_mid_3 = ev[e3];
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		EdgeHandle e4 = check_edge(i, 4, 0);
		V3f edge_mid_4 = ev[e4];
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		EdgeHandle e5 = check_edge(i, 3, 2);
		V3f edge_mid_5 = ev[e5];
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		EdgeHandle e6 = check_edge(i, 3, 7);
		V3f edge_mid_6 = ev[e6];
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		EdgeHandle e7 = check_edge(i, 6, 7);
		V3f edge_mid_7 = ev[e7];
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		EdgeHandle e8 = check_edge(i, 2, 6);
		V3f edge_mid_8 = ev[e8];
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		EdgeHandle e9 = check_edge(i, 1, 2);
		V3f edge_mid_9 = ev[e9];
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		EdgeHandle e10 = check_edge(i, 5, 6);
		V3f edge_mid_10 = ev[e10];
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		EdgeHandle e11 = check_edge(i, 0, 3);
		V3f edge_mid_11 = ev[e11];
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		EdgeHandle e12 = check_edge(i, 4, 7);
		V3f edge_mid_12 = ev[e12];
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			FaceHandle new_fv = mesh->cells_[i].faces_[mesh_face_num];
			V3f face_mid = fv[new_fv];
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cv[i];
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

void InterpolateCC::PBfirstep()
{
	firstver.clear();
	detver.clear();
	/*firstedge.swap(std::vector<V3f>());
	detedge.swap(std::vector<V3f>());
	firstface.swap(std::vector<V3f>());
	detface.swap(std::vector<V3f>());*/
	//oriv.swap(std::vector<V3f>());
	int ver_size = mesh->vertices_.size();
	for (int i = 0; i < ver_size; i++)
	{
		V3f new_vertex = cal_new_vertex(VertexHandle(i));
		V3f det = mesh->vertices_[i] - new_vertex;
		firstver.push_back(new_vertex);
		detver.push_back(det);
		//oriv.push_back(mesh->vertices_[i]);
	}
	/*int edge_size = mesh->edges_.size();
	for (int i = 0;i < edge_size;i++)
	{
		V3f new_vertex = cal_edge_mid(EdgeHandle(i));
		V3f det = cal_edge_midv(EdgeHandle(i)) - new_vertex;
		firstedge.push_back(new_vertex);
		detedge.push_back(det);
	}
	int face_size = mesh->faces_.size();
	for (int i = 0;i < face_size;i++)
	{
		V3f new_vertex = cal_face_mid(FaceHandle(i));
		V3f det = cal_face_midv(FaceHandle(i)) - new_vertex;
		firstface.push_back(new_vertex);
		detface.push_back(det);
	}*/
	/*vv.resize(ver_size);
	for (int i = 0; i < ver_size; i++)
	{
		V3f temp = cal_CIPB_nv(VertexHandle(i));
		vv[i] = temp;
	}
	int edge_size = mesh->edges_.size();
	ev.resize(edge_size);
	for (int i = 0;i < edge_size;i++)
	{
		V3f temp = cal_CIPB_ev(EdgeHandle(i));
		ev[i] = temp;
	}
	int face_size = mesh->faces_.size();
	fv.resize(face_size);
	for (int i = 0;i < face_size;i++)
	{
		V3f temp = cal_CIPB_fv(FaceHandle(i));
		fv[i] = temp;
	}
	int cell_size = mesh->cells_.size();
	cv.resize(cell_size);
	for (int i = 0;i < cell_size;i++)
	{
		V3f temp = cal_CIPB_cv(CellHandle(i));
		cv[i] = temp;
	}*/
}

V3f InterpolateCC::cal_CIPB_ev(const int cell_num, const int ver1, const int ver2)
{
	EdgeHandle edge = check_edge(cell_num, ver1, ver2);
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge))
	{
		edge_mid = cal_edge_midv(edge) + (detver[mesh->edges_[edge].from_h] + detver[mesh->edges_[edge].to_h]) * le;
	}
	else
	{
		//edge_mid = cal_edge_mid(cell_num, ver1, ver2);
		//edge_mid = cal_edge_midv(edge) + detedge[edge] * le;
		edge_mid = cal_edge_midv(edge) + (detver[mesh->edges_[edge].from_h] + detver[mesh->edges_[edge].to_h]) * le;
	}
	return edge_mid;
}

V3f InterpolateCC::cal_CIPB_ev(EdgeHandle edge_id)
{
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge_id))
	{
		edge_mid = cal_edge_midv(edge_id) + (detver[mesh->edges_[edge_id].from_h] + detver[mesh->edges_[edge_id].to_h]) * le;
	}
	else
	{
		//edge_mid = cal_edge_mid(edge_id);
		//edge_mid = cal_edge_midv(edge_id) + detedge[edge_id] * le;
		edge_mid = cal_edge_midv(edge_id) + (detver[mesh->edges_[edge_id].from_h] + detver[mesh->edges_[edge_id].to_h]) * le;
	}
	return edge_mid;
}

V3f InterpolateCC::cal_CIPB_fv(const int cell_num, const int face_num)
{
	V3f face_mid;
	V3f face_sum;
	FaceHandle face = mesh->cells_[cell_num].faces_[face_num];//确定该面
	if (mesh->is_bdy(face))
	{
		std::vector<VertexHandle> fv = mesh->faces_[face].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mesh->vertices_[*fv_it] / 4.0;
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
		face_mid = face_sum;
	}
	else
	{
		//face_mid = cal_face_mid(cell_num, face_num);
		//face_mid = cal_face_midv(face) + detface[face] * 2 * uf;
		face_mid = face_mid + cal_face_midv(face);
		std::vector<VertexHandle> fv = mesh->faces_[face].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
	}
	return face_mid;
}

V3f InterpolateCC::cal_CIPB_fv(FaceHandle face_id)
{
	V3f face_mid;
	V3f face_sum;
	if (mesh->is_bdy(face_id))
	{
		std::vector<VertexHandle> fv = mesh->faces_[face_id].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mesh->vertices_[*fv_it] / 4.0;
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
		face_mid = face_sum;
	}
	else
	{
		//face_mid = cal_face_mid(face_id);
		//face_mid = cal_face_midv(face_id) + detface[face_id] * 2 * uf;
		face_mid = face_mid + cal_face_midv(face_id);
		std::vector<VertexHandle> fv = mesh->faces_[face_id].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
	}
	return face_mid;
}

V3f InterpolateCC::cal_CIPB_cv(const int cell_num)
{
	V3f cell_mid;
	//cell_mid = cal_cell_mid(cell_num);
	std::vector<VertexHandle> c_v = mesh->cells_[cell_num].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it] / 8.0;
		cell_mid = cell_mid + detver[*c_v_it] * 2 * dc;
	}
	//std::set<EdgeHandle> c_e;
	////将块的所含边放入容器
	//std::vector<FaceHandle> c_f = mesh->cells_[cell_num].faces_;
	//for (std::vector<FaceHandle>::iterator c_f_it = c_f.begin();c_f_it != c_f.end();c_f_it++)
	//{
	//	std::vector<EdgeHandle> c_f_e = mesh->faces_[*c_f_it].edges_;
	//	for (std::vector<EdgeHandle>::iterator c_f_e_it = c_f_e.begin();c_f_e_it != c_f_e.end();c_f_e_it++)
	//	{
	//		c_e.insert(*c_f_e_it);
	//	}
	//}
	//for (std::set<EdgeHandle>::iterator c_e_it = c_e.begin();c_e_it != c_e.end();c_e_it++)
	//{
	//	cell_mid = cell_mid + cal_edge_midv(*c_e_it) / 8.0;
	//	cell_mid = cell_mid + detedge[*c_e_it] * 2 * dc;
	//}
	return cell_mid;
}

V3f InterpolateCC::cal_CIPB_cv(CellHandle cell_id)
{
	V3f cell_mid;
	//cell_mid = cal_cell_mid(cell_id);
	std::vector<VertexHandle> c_v = mesh->cells_[cell_id].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it] / 8.0;
		cell_mid = cell_mid + detver[*c_v_it] * 2 * dc;
	}
	//std::set<EdgeHandle> c_e;
	////将块的所含边放入容器
	//std::vector<FaceHandle> c_f = mesh->cells_[cell_id].faces_;
	//for (std::vector<FaceHandle>::iterator c_f_it = c_f.begin();c_f_it != c_f.end();c_f_it++)
	//{
	//	std::vector<EdgeHandle> c_f_e = mesh->faces_[*c_f_it].edges_;
	//	for (std::vector<EdgeHandle>::iterator c_f_e_it = c_f_e.begin();c_f_e_it != c_f_e.end();c_f_e_it++)
	//	{
	//		c_e.insert(*c_f_e_it);
	//	}
	//}
	//for (std::set<EdgeHandle>::iterator c_e_it = c_e.begin();c_e_it != c_e.end();c_e_it++)
	//{
	//	cell_mid = cell_mid + cal_edge_midv(*c_e_it) / 8.0;
	//	cell_mid = cell_mid + detedge[*c_e_it] * 2 * dc;
	//}
	return cell_mid;
}

V3f InterpolateCC::cal_CIPB_nv(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		//const int n = mesh->neighbor_v[vertex].size();
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * n * (n + 5);
		int en = 0;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int v_e_n = v_e.size();
		V3f e_av;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				en++;
				V3f temp_ev = cal_CIPB_ev(*v_e_it);
				e_av = e_av + temp_ev;
				//vertex_sum = vertex_sum - temp_ev * 4;
			}
		}
		e_av = e_av / en;
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (en + 5);
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		V3f f_av;
		//int fn = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				//fn++;
				V3f temp_fv = cal_CIPB_fv(*v_f_it);
				f_av = f_av + temp_fv;
				//vertex_sum = vertex_sum - temp_fv;
			}
		}
		f_av = f_av / en;
		new_vertex = vertex_sum - e_av * 4 - f_av;
		new_vertex = new_vertex / en;
		//new_vertex = mesh->vertices_[vertex];
	}
	else
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int en = v_e.size();
		int sum = 0;
		sum += 16 * (en - 2);
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * (4 * 4 * (en-2)+ 38 * (en - 2));
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			V3f temp_ev = cal_CIPB_ev(*v_e_it);
			int m = mesh->neighbor_e[*v_e_it].size();
			sum += 4 * m;
			vertex_sum = vertex_sum - temp_ev * 4 * m;
			/*if (mesh->is_singular(*v_e_it))
			{
				vertex_sum = vertex_sum - temp_ev * 4 * (en - 2);
			}
			else
			{
				vertex_sum = vertex_sum - temp_ev * 4 * 4;
			}*/
		}
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			V3f temp_fv = cal_CIPB_fv(*v_f_it);
			vertex_sum = vertex_sum - temp_fv * 4;
			sum += 4;
		}
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			V3f temp_cv = cal_CIPB_cv(*v_c_it);
			vertex_sum = vertex_sum - temp_cv;
			sum += 1;
		}
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * sum;
		new_vertex = vertex_sum / ((en - 2) * 4 * 4);
	}
	return new_vertex;
}

V3f InterpolateCC::cal_CIPB_nv(VertexHandle ver_id)
{
	VertexHandle vertex = ver_id;
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		//const int n = mesh->neighbor_v[vertex].size();
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * n * (n + 5);
		int en = 0;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int v_e_n = v_e.size();
		V3f e_av;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				en++;
				V3f temp_ev = cal_CIPB_ev(*v_e_it);
				e_av = e_av + temp_ev;
				//vertex_sum = vertex_sum - temp_ev * 4;
			}
		}
		e_av = e_av / en;
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (en + 5);
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		V3f f_av;
		//int fn = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				//fn++;
				V3f temp_fv = cal_CIPB_fv(*v_f_it);
				f_av = f_av + temp_fv;
				//vertex_sum = vertex_sum - temp_fv;
			}
		}
		f_av = f_av / en;
		new_vertex = vertex_sum - e_av * 4 - f_av;
		new_vertex = new_vertex / en;
		//new_vertex = mesh->vertices_[vertex];
	}
	else
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int en = v_e.size();
		int sum = 0;
		sum += 16 * (en - 2);
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * (4 * 4 * (en-2)+ 38 * (en - 2));
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			V3f temp_ev = cal_CIPB_ev(*v_e_it);
			int m = mesh->neighbor_e[*v_e_it].size();
			sum += 4 * m;
			vertex_sum = vertex_sum - temp_ev * 4 * m;
			/*if (mesh->is_singular(*v_e_it))
			{
			vertex_sum = vertex_sum - temp_ev * 4 * (en - 2);
			}
			else
			{
			vertex_sum = vertex_sum - temp_ev * 4 * 4;
			}*/
		}
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			V3f temp_fv = cal_CIPB_fv(*v_f_it);
			vertex_sum = vertex_sum - temp_fv * 4;
			sum += 4;
		}
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			V3f temp_cv = cal_CIPB_cv(*v_c_it);
			vertex_sum = vertex_sum - temp_cv;
			sum += 1;
		}
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * sum;
		new_vertex = vertex_sum / ((en - 2) * 4 * 4);
	}
	return new_vertex;
}

V3f InterpolateCC::cal_CIPB_fnv(VertexHandle ver_id)
{
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(ver_id))
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[ver_id];
		int v_e_n = 0;
		V3f v_e_mid;//计算点邻接边的平均边点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				v_e_n++;
				v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
			}
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = vertex_sum + v_e_mid * 2;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[ver_id];
		int v_f_n = 0;
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				v_f_n++;
				v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
			}
		}
		v_f_mid = v_f_mid / v_f_n;
		vertex_sum = vertex_sum + v_f_mid;
		vertex_sum = vertex_sum + vv[ver_id] * (v_e_n - 3);
		new_vertex = vertex_sum / v_e_n;
	}
	else
	{
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[ver_id];
		const int v_c_n = v_c.size();
		V3f v_c_mid;//计算点邻接块的平均块点
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			v_c_mid = v_c_mid + cal_cell_midv(*v_c_it);
		}
		v_c_mid = v_c_mid / v_c_n;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[ver_id];
		const int v_f_n = v_f.size();
		V3f v_f_mid;//计算点邻接面的平均面点
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			v_f_mid = v_f_mid + cal_face_midv(*v_f_it);
		}
		v_f_mid = v_f_mid / v_f_n;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[ver_id];
		const int v_e_n = v_e.size();
		V3f v_e_mid;//计算点邻接面的平均面点
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			v_e_mid = v_e_mid + cal_edge_midv(*v_e_it);
		}
		v_e_mid = v_e_mid / v_e_n;
		vertex_sum = v_c_mid + v_f_mid * 3 + v_e_mid * 3 + vv[ver_id];
		new_vertex = vertex_sum / 8;
		//}
	}
	return new_vertex;
}

/*-----------Normal-Based法-----------*/
HexV3fMesh& InterpolateCC::CCInterpolatingNB()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_CINB_nv(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(firstver[mesh->cells_[i].vertex_[mesh_vex_num]]);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_CINB_ev(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_CINB_ev(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_CINB_ev(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_CINB_ev(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_CINB_ev(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_CINB_ev(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_CINB_ev(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_CINB_ev(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_CINB_ev(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_CINB_ev(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_CINB_ev(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_CINB_ev(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_CINB_fv(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

void InterpolateCC::NBfirstep()
{
	normalver.clear();
	int ver_size = mesh->vertices_.size();
	for (int i = 0; i < ver_size; i++)
	{
		VertexHandle vertex = VertexHandle(i);
		double fi; V3f nor; V3f nx;
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			std::vector<EdgeHandle> ne;
			if (mesh->is_bdy(*v_f_it))
			{
				std::vector<EdgeHandle> v_f_e = mesh->faces_[*v_f_it].edges_;
				for (std::vector<EdgeHandle>::iterator v_f_e_it = v_f_e.begin(); v_f_e_it != v_f_e.end(); v_f_e_it++)
				{
					if (mesh->edges_[*v_f_e_it].from_h == vertex || mesh->edges_[*v_f_e_it].to_h == vertex)
						ne.push_back(*v_f_e_it);
				}
				fi = mesh->cal_angle(ne[0], ne[1], vertex);
				nor = mesh->cal_norm(*v_f_it);
				//修正向量
				std::set<CellHandle, compare_OVM> v_f_c = mesh->neighbor_f[*v_f_it];
				CellHandle vfc = *(v_f_c.begin());
				V3f vfc_mid = cal_cell_midv(vfc);
				V3f checkn = vfc_mid - mesh->vertices_[vertex];
				double temp = nor.x * checkn.x + nor.y * checkn.y + nor.z * checkn.z;
				if (temp > 0)
					nor = nor * -1;
				nx = nx + nor * fi;
			}
		}
		nx = nx.normalize();
		normalver.push_back(nx);
	}
}

V3f InterpolateCC::cal_CINB_ev(const int cell_num, const int ver1, const int ver2)
{
	EdgeHandle edge = check_edge(cell_num, ver1, ver2);
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge))
	{
		V3f edmid = cal_edge_midv(edge);
		edge_sum = edge_sum + edmid;
		V3f di = (mesh->vertices_[mesh->edges_[edge].from_h] - mesh->vertices_[mesh->edges_[edge].to_h]) * 0.5;
		V3f nori = normalver[mesh->edges_[edge].from_h];
		V3f dj = (mesh->vertices_[mesh->edges_[edge].to_h] - mesh->vertices_[mesh->edges_[edge].from_h]) * 0.5;
		V3f norj = normalver[mesh->edges_[edge].to_h];
		V3f temp = nori * (di.x * nori.x + di.y * nori.y + di.z * nori.z) + norj * (dj.x * norj.x + dj.y * norj.y + dj.z * norj.z);
		edge_sum = edge_sum + temp * we * 0.5;
		edge_mid = edge_sum;
	}
	else
		edge_mid = cal_edge_mid(edge);
	return edge_mid;
}

V3f InterpolateCC::cal_CINB_ev(EdgeHandle edge_id)
{
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge_id))
	{
		V3f edmid = cal_edge_midv(edge_id);
		edge_sum = edge_sum + edmid;
		V3f di = (mesh->vertices_[mesh->edges_[edge_id].from_h] - mesh->vertices_[mesh->edges_[edge_id].to_h]) * 0.5;
		V3f nori = normalver[mesh->edges_[edge_id].from_h];
		V3f dj = (mesh->vertices_[mesh->edges_[edge_id].to_h] - mesh->vertices_[mesh->edges_[edge_id].from_h]) * 0.5;
		V3f norj = normalver[mesh->edges_[edge_id].to_h];
		V3f temp = nori * (di.x * nori.x + di.y * nori.y + di.z * nori.z) + norj * (dj.x * norj.x + dj.y * norj.y + dj.z * norj.z);
		edge_sum = edge_sum + temp * we * 0.5;
		edge_mid = edge_sum;
	}
	else
		edge_mid = cal_edge_mid(edge_id);
	return edge_mid;
}

V3f InterpolateCC::cal_CINB_fv(const int cell_num, const int face_num)
{
	V3f face_mid;
	V3f face_sum;
	FaceHandle face = mesh->cells_[cell_num].faces_[face_num];//确定该面
	if (mesh->is_bdy(face))
	{
		V3f face_mid = cal_face_midv(face);
		face_sum = face_sum + face_mid;
		std::vector<VertexHandle> fv = mesh->faces_[face].vs_f;
		V3f temp;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			V3f temp1 = mesh->vertices_[*fv_it] - face_mid;
			V3f nori = normalver[*fv_it];
			temp = temp + nori * (temp1.x * nori.x + temp1.y * nori.y + temp1.z * nori.z);
		}
		face_sum = face_sum + temp * vf / 4;
		face_mid = face_sum;
	}
	else
		face_mid = cal_face_mid(face);
	return face_mid;
}

V3f InterpolateCC::cal_CINB_fv(FaceHandle face_id)
{
	V3f face_mid;
	V3f face_sum;
	if (mesh->is_bdy(face_id))
	{
		V3f face_mid = cal_face_midv(face_id);
		face_sum = face_sum + face_mid;
		std::vector<VertexHandle> fv = mesh->faces_[face_id].vs_f;
		V3f temp;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			V3f temp1 = mesh->vertices_[*fv_it] - face_mid;
			V3f nori = normalver[*fv_it];
			temp = temp + nori * (temp1.x * nori.x + temp1.y * nori.y + temp1.z * nori.z);
		}
		face_sum = face_sum + temp * vf / 4;
		face_mid = face_sum;
	}
	else
		face_mid = cal_face_mid(face_id);
	return face_mid;
}

V3f InterpolateCC::cal_CINB_nv(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		//const int n = mesh->neighbor_v[vertex].size();
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * n * (n + 5);
		int en = 0;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int v_e_n = v_e.size();
		V3f e_av;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				en++;
				V3f temp_ev = cal_CINB_ev(*v_e_it);
				e_av = e_av + temp_ev;
				//vertex_sum = vertex_sum - temp_ev * 4;
			}
		}
		e_av = e_av / en;
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (en + 5);
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		V3f f_av;
		//int fn = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				//fn++;
				V3f temp_fv = cal_CINB_fv(*v_f_it);
				f_av = f_av + temp_fv;
				//vertex_sum = vertex_sum - temp_fv;
			}
		}
		f_av = f_av / en;
		new_vertex = vertex_sum - e_av * 4 - f_av;
		new_vertex = new_vertex / en;
		//new_vertex = mesh->vertices_[vertex];
	}
	else
		new_vertex = cal_new_vertex(cell_num, ver_num);
	return new_vertex;
}

/*-----------MLCA部分-----------*/
void InterpolateCC::MLCAfirstep()
{
	Spliting();
	mlcadet.clear();
	int ver_size = mesh->vertices_.size();
	for (int i = 0; i < ver_size; i++)
	{
		V3f new_vertex = cal_mlca_vertex(VertexHandle(i));
		V3f det = mesh->vertices_[i] - new_vertex;
		mlcadet.push_back(det);
	}
}

void InterpolateCC::Spliting()
{
	uint64_t cell_num = mesh->cells_.size();
	vid.resize(cell_num * 27, VertexHandle(-1));
	/*eid.resize(cell_num * 69);
	fid.resize(cell_num * 36);
	cid.resize(cell_num * 8);*/
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		std::vector<V3f>::iterator v_id;

		// 傻逼了 直接for不就好了

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			temp_vex[mesh_vex_num] = tempmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
			//tempmesh->vertices_[temp_vex[mesh_vex_num]].lastid = mesh->cells_[i].vertex_[mesh_vex_num];
			vid[temp_vex[mesh_vex_num]] = mesh->cells_[i].vertex_[mesh_vex_num];
			/*eid[temp_vex[mesh_vex_num]] = EdgeHandle(-1);
			fid[temp_vex[mesh_vex_num]] = FaceHandle(-1);
			cid[temp_vex[mesh_vex_num]] = CellHandle(-1);*/
		}
		//再插入每条边的中点 暴力插入 以免顺序错误
		//第一条边的中点
		V3f edge_mid_1;
		edge_mid_1.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[1]].x) * 1.0 / 2 * 1.0;
		edge_mid_1.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[1]].y) * 1.0 / 2 * 1.0;
		edge_mid_1.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[1]].z) * 1.0 / 2 * 1.0;
		temp_vex[8] = tempmesh->add_vertex(edge_mid_1);
		//eid[temp_vex[8]] = check_edge(i, 0, 1);

		//第二条边的中点
		V3f edge_mid_2;
		edge_mid_2.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[5]].x) * 1.0 / 2 * 1.0;
		edge_mid_2.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[5]].y) * 1.0 / 2 * 1.0;
		edge_mid_2.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[5]].z) * 1.0 / 2 * 1.0;
		temp_vex[9] = tempmesh->add_vertex(edge_mid_2);
		//eid[temp_vex[9]] = check_edge(i, 1, 5);

		//第三条边的中点
		V3f edge_mid_3;
		edge_mid_3.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[4]].x) * 1.0 / 2 * 1.0;
		edge_mid_3.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[4]].y) * 1.0 / 2 * 1.0;
		edge_mid_3.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[4]].z) * 1.0 / 2 * 1.0;
		temp_vex[10] = tempmesh->add_vertex(edge_mid_3);
		//eid[temp_vex[10]] = check_edge(i, 5, 4);

		//第四条边的中点
		V3f edge_mid_4;
		edge_mid_4.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[0]].x) * 1.0 / 2 * 1.0;
		edge_mid_4.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[0]].y) * 1.0 / 2 * 1.0;
		edge_mid_4.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[0]].z) * 1.0 / 2 * 1.0;
		temp_vex[11] = tempmesh->add_vertex(edge_mid_4);
		//eid[temp_vex[11]] = check_edge(i, 4, 0);

		//第五条边的中点
		V3f edge_mid_5;
		edge_mid_5.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_5.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_5.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[12] = tempmesh->add_vertex(edge_mid_5);
		///eid[temp_vex[12]] = check_edge(i, 3, 2);

		//第六条边的中点
		V3f edge_mid_6;
		edge_mid_6.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_6.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_6.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[13] = tempmesh->add_vertex(edge_mid_6);
		//eid[temp_vex[13]] = check_edge(i, 3, 7);

		//第七条边的中点
		V3f edge_mid_7;
		edge_mid_7.x = (mesh->vertices_[mesh->cells_[i].vertex_[6]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_7.y = (mesh->vertices_[mesh->cells_[i].vertex_[6]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_7.z = (mesh->vertices_[mesh->cells_[i].vertex_[6]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[14] = tempmesh->add_vertex(edge_mid_7);
		//eid[temp_vex[14]] = check_edge(i, 6, 7);

		//第八条边的中点
		V3f edge_mid_8;
		edge_mid_8.x = (mesh->vertices_[mesh->cells_[i].vertex_[2]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_8.y = (mesh->vertices_[mesh->cells_[i].vertex_[2]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_8.z = (mesh->vertices_[mesh->cells_[i].vertex_[2]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[15] = tempmesh->add_vertex(edge_mid_8);
		//eid[temp_vex[15]] = check_edge(i, 2, 6);

		//第九条边的中点
		V3f edge_mid_9;
		edge_mid_9.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_9.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_9.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[16] = tempmesh->add_vertex(edge_mid_9);
		//eid[temp_vex[16]] = check_edge(i, 1, 2);

		//第十条边的中点
		V3f edge_mid_10;
		edge_mid_10.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_10.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_10.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[17] = tempmesh->add_vertex(edge_mid_10);
		//eid[temp_vex[17]] = check_edge(i, 5, 6);

		//第十一条边的中点
		V3f edge_mid_11;
		edge_mid_11.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[3]].x) * 1.0 / 2 * 1.0;
		edge_mid_11.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[3]].y) * 1.0 / 2 * 1.0;
		edge_mid_11.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[3]].z) * 1.0 / 2 * 1.0;
		temp_vex[18] = tempmesh->add_vertex(edge_mid_11);
		//eid[temp_vex[18]] = check_edge(i, 0, 3);

		//第十二条边的中点
		V3f edge_mid_12;
		edge_mid_12.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_12.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_12.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[19] = tempmesh->add_vertex(edge_mid_12);
		//eid[temp_vex[19]] = check_edge(i, 4, 7);

		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid;
			face_mid.x = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].x) * 1.0 / 4 * 1.0;
			face_mid.y = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].y) * 1.0 / 4 * 1.0;
			face_mid.z = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].z) * 1.0 / 4 * 1.0;
			temp_vex[20 + mesh_face_num] = tempmesh->add_vertex(face_mid);
			//fid[temp_vex[20 + mesh_face_num]] = mesh->cells_[i].faces_[mesh_face_num];
		}

		//插入单元的内心
		V3f cell_mid;
		double x_sum = 0;  double y_sum = 0;  double z_sum = 0;
		for (int j = 0; j < 8; j++)
		{
			x_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].x;
			y_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].y;
			z_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].z;
		}
		cell_mid.x = x_sum * 1.0 / 8 * 1.0;
		cell_mid.y = y_sum * 1.0 / 8 * 1.0;
		cell_mid.z = z_sum * 1.0 / 8 * 1.0;
		temp_vex[26] = tempmesh->add_vertex(cell_mid);
		//cid[temp_vex[26]] = CellHandle(i);

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		tempmesh->add_cell(sub_cell[0]);
		tempmesh->add_cell(sub_cell[1]);
		tempmesh->add_cell(sub_cell[2]);
		tempmesh->add_cell(sub_cell[3]);
		tempmesh->add_cell(sub_cell[4]);
		tempmesh->add_cell(sub_cell[5]);
		tempmesh->add_cell(sub_cell[6]);
		tempmesh->add_cell(sub_cell[7]);
	}
	//return *tempmesh;
}

V3f InterpolateCC::cal_mlca_vertex(VertexHandle ver_id)
{
	V3f new_vertex;
	if (mesh->is_bdy(ver_id))
		new_vertex = cal_new_vertex(ver_id);
	else
	{
		int size = vid.size();
		for (int i = 0; i < size; i++)
		{
			if (vid[i] == VertexHandle(ver_id))
			{
				std::set<CellHandle, compare_OVM> vc = tempmesh->neighbor_v_c[i];
				int vc_num = vc.size();
				std::vector<V3f> vc_cen(vc_num);
				int u = 0;
				for (std::set<CellHandle, compare_OVM>::iterator vc_it = vc.begin(); vc_it != vc.end(); vc_it++)
				{
					vc_cen[u] = cal_cell_midv(*vc_it);
					u++;
				}
				for (int vf_cen_num = 0; vf_cen_num < vc_num; vf_cen_num++)
				{
					new_vertex = new_vertex + vc_cen[vf_cen_num];
				}
				new_vertex = new_vertex / (vc_num * 1.0);
				break;
			}
		}
	}
	return new_vertex;
}

V3f InterpolateCC::cal_MLCA_ev(const int cell_num, const int ver1, const int ver2)
{
	EdgeHandle edge = check_edge(cell_num, ver1, ver2);
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge))
	{
		edge_mid = cal_edge_midv(edge) + (detver[mesh->edges_[edge].from_h] + detver[mesh->edges_[edge].to_h]) * le;
	}
	else
	{
		//edge_mid = cal_edge_mid(cell_num, ver1, ver2);
		//edge_mid = cal_edge_midv(edge) + detedge[edge] * le;
		edge_mid = cal_edge_midv(edge) + (mlcadet[mesh->edges_[edge].from_h] + mlcadet[mesh->edges_[edge].to_h]) * le;
	}
	return edge_mid;
}

V3f InterpolateCC::cal_MLCA_ev(EdgeHandle edge_id)
{
	V3f edge_mid;
	V3f edge_sum;
	if (mesh->is_bdy(edge_id))
	{
		edge_mid = cal_edge_midv(edge_id) + (detver[mesh->edges_[edge_id].from_h] + detver[mesh->edges_[edge_id].to_h]) * le;
	}
	else
	{
		//edge_mid = cal_edge_mid(cell_num, ver1, ver2);
		//edge_mid = cal_edge_midv(edge) + detedge[edge] * le;
		edge_mid = cal_edge_midv(edge_id) + (mlcadet[mesh->edges_[edge_id].from_h] + mlcadet[mesh->edges_[edge_id].to_h]) * le;
	}
	return edge_mid;
}

V3f InterpolateCC::cal_MLCA_fv(const int cell_num, const int face_num)
{
	V3f face_mid;
	V3f face_sum;
	FaceHandle face = mesh->cells_[cell_num].faces_[face_num];//确定该面
	if (mesh->is_bdy(face))
	{
		std::vector<VertexHandle> fv = mesh->faces_[face].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mesh->vertices_[*fv_it] / 4.0;
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
		face_mid = face_sum;
	}
	else
	{
		//face_mid = cal_face_mid(cell_num, face_num);
		//face_mid = cal_face_midv(face) + detface[face] * 2 * uf;
		face_mid = face_mid + cal_face_midv(face);
		std::vector<VertexHandle> fv = mesh->faces_[face].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mlcadet[*fv_it] * 2 * uf;
		}
	}
	return face_mid;
}

V3f InterpolateCC::cal_MLCA_fv(FaceHandle face_id)
{
	V3f face_mid;
	V3f face_sum;
	if (mesh->is_bdy(face_id))
	{
		std::vector<VertexHandle> fv = mesh->faces_[face_id].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mesh->vertices_[*fv_it] / 4.0;
			face_sum = face_sum + detver[*fv_it] * 2 * uf;
		}
		face_mid = face_sum;
	}
	else
	{
		//face_mid = cal_face_mid(face_id);
		//face_mid = cal_face_midv(face_id) + detface[face_id] * 2 * uf;
		face_mid = face_mid + cal_face_midv(face_id);
		std::vector<VertexHandle> fv = mesh->faces_[face_id].vs_f;
		for (std::vector<VertexHandle>::iterator fv_it = fv.begin(); fv_it != fv.end(); fv_it++)
		{
			face_sum = face_sum + mlcadet[*fv_it] * 2 * uf;
		}
	}
	return face_mid;
}

V3f InterpolateCC::cal_MLCA_cv(const int cell_num)
{
	V3f cell_mid;
	std::vector<VertexHandle> c_v = mesh->cells_[cell_num].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it] / 8.0;
		cell_mid = cell_mid + mlcadet[*c_v_it] * 2 * dc;
	}
	return cell_mid;
}

V3f InterpolateCC::cal_MLCA_cv(CellHandle cell_id)
{
	V3f cell_mid;
	std::vector<VertexHandle> c_v = mesh->cells_[cell_id].vertex_;
	for (std::vector<VertexHandle>::iterator c_v_it = c_v.begin(); c_v_it != c_v.end(); c_v_it++)
	{
		cell_mid = cell_mid + mesh->vertices_[*c_v_it] / 8.0;
		cell_mid = cell_mid + mlcadet[*c_v_it] * 2 * dc;
	}
	return cell_mid;
}

V3f InterpolateCC::cal_MLCA_nv(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	V3f new_vertex;
	V3f vertex_sum;
	if (mesh->is_bdy(vertex))
	{
		//const int n = mesh->neighbor_v[vertex].size();
		//vertex_sum = vertex_sum + mesh->vertices_[vertex] * n * (n + 5);
		int en = 0;
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int v_e_n = v_e.size();
		V3f e_av;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
			{
				en++;
				V3f temp_ev = cal_CIPB_ev(*v_e_it);
				e_av = e_av + temp_ev;
				//vertex_sum = vertex_sum - temp_ev * 4;
			}
		}
		e_av = e_av / en;
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (en + 5);
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		V3f f_av;
		//int fn = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			if (mesh->is_bdy(*v_f_it))
			{
				//fn++;
				V3f temp_fv = cal_CIPB_fv(*v_f_it);
				f_av = f_av + temp_fv;
				//vertex_sum = vertex_sum - temp_fv;
			}
		}
		f_av = f_av / en;
		new_vertex = vertex_sum - e_av * 4 - f_av;
		new_vertex = new_vertex / en;
		//new_vertex = mesh->vertices_[vertex];
	}
	else
	{
		std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
		int en = v_e.size();
		vertex_sum = vertex_sum + mesh->vertices_[vertex] * (27 + 20 * (en - 2));
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
		{
			V3f temp_ev = cal_CIPB_ev(*v_e_it);
			if (mesh->is_singular(*v_e_it))
			{
				vertex_sum = vertex_sum - temp_ev * (en - 2) * 4;
			}
			else
			{
				vertex_sum = vertex_sum - temp_ev * 3 * 3;
			}
		}
		std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
		for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
		{
			V3f temp_fv = cal_CIPB_fv(*v_f_it);
			vertex_sum = vertex_sum - temp_fv * 3;
		}
		std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
		for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
		{
			V3f temp_cv = cal_CIPB_cv(*v_c_it);
			vertex_sum = vertex_sum - temp_cv;
		}
		new_vertex = vertex_sum / 27;
		if (mesh->is_singular(vertex))
		{
			//修正奇异点
			//new_vertex = new_vertex - detver[vertex] * 4 * sv;
			/*int size1 = mesh->neighbor_v_f[vertex].size();
			int size2 = mesh->neighbor_v_c[vertex].size();
			int size3 = mesh->neighbor_v[vertex].size();*/
		}
	}
	return new_vertex;
}

HexV3fMesh& InterpolateCC::CCInterpolatingMLCA()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		//std::vector<V3f>::iterator v_id;

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_v = cal_MLCA_nv(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_v);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(firstver[mesh->cells_[i].vertex_[mesh_vex_num]]);
			//temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_MLCA_ev(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_MLCA_ev(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_MLCA_ev(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_MLCA_ev(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_MLCA_ev(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_MLCA_ev(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_MLCA_ev(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_MLCA_ev(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_MLCA_ev(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_MLCA_ev(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_MLCA_ev(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_MLCA_ev(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_MLCA_fv(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cell_mid(i);
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

//-----MLCASubdivision-----//
HexV3fMesh& MLCASubdivision::Spliting()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		std::vector<V3f>::iterator v_id;

		// 傻逼了 直接for不就好了

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			temp_vex[mesh_vex_num] = newmesh->add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[mesh_vex_num]]);
		}
		//再插入每条边的中点 暴力插入 以免顺序错误
		//第一条边的中点
		V3f edge_mid_1;
		edge_mid_1.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[1]].x) * 1.0 / 2 * 1.0;
		edge_mid_1.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[1]].y) * 1.0 / 2 * 1.0;
		edge_mid_1.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[1]].z) * 1.0 / 2 * 1.0;
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2;
		edge_mid_2.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[5]].x) * 1.0 / 2 * 1.0;
		edge_mid_2.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[5]].y) * 1.0 / 2 * 1.0;
		edge_mid_2.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[5]].z) * 1.0 / 2 * 1.0;
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3;
		edge_mid_3.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[4]].x) * 1.0 / 2 * 1.0;
		edge_mid_3.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[4]].y) * 1.0 / 2 * 1.0;
		edge_mid_3.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[4]].z) * 1.0 / 2 * 1.0;
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4;
		edge_mid_4.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[0]].x) * 1.0 / 2 * 1.0;
		edge_mid_4.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[0]].y) * 1.0 / 2 * 1.0;
		edge_mid_4.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[0]].z) * 1.0 / 2 * 1.0;
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5;
		edge_mid_5.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_5.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_5.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6;
		edge_mid_6.x = (mesh->vertices_[mesh->cells_[i].vertex_[3]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_6.y = (mesh->vertices_[mesh->cells_[i].vertex_[3]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_6.z = (mesh->vertices_[mesh->cells_[i].vertex_[3]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7;
		edge_mid_7.x = (mesh->vertices_[mesh->cells_[i].vertex_[6]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_7.y = (mesh->vertices_[mesh->cells_[i].vertex_[6]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_7.z = (mesh->vertices_[mesh->cells_[i].vertex_[6]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8;
		edge_mid_8.x = (mesh->vertices_[mesh->cells_[i].vertex_[2]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_8.y = (mesh->vertices_[mesh->cells_[i].vertex_[2]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_8.z = (mesh->vertices_[mesh->cells_[i].vertex_[2]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9;
		edge_mid_9.x = (mesh->vertices_[mesh->cells_[i].vertex_[1]].x + mesh->vertices_[mesh->cells_[i].vertex_[2]].x) * 1.0 / 2 * 1.0;
		edge_mid_9.y = (mesh->vertices_[mesh->cells_[i].vertex_[1]].y + mesh->vertices_[mesh->cells_[i].vertex_[2]].y) * 1.0 / 2 * 1.0;
		edge_mid_9.z = (mesh->vertices_[mesh->cells_[i].vertex_[1]].z + mesh->vertices_[mesh->cells_[i].vertex_[2]].z) * 1.0 / 2 * 1.0;
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10;
		edge_mid_10.x = (mesh->vertices_[mesh->cells_[i].vertex_[5]].x + mesh->vertices_[mesh->cells_[i].vertex_[6]].x) * 1.0 / 2 * 1.0;
		edge_mid_10.y = (mesh->vertices_[mesh->cells_[i].vertex_[5]].y + mesh->vertices_[mesh->cells_[i].vertex_[6]].y) * 1.0 / 2 * 1.0;
		edge_mid_10.z = (mesh->vertices_[mesh->cells_[i].vertex_[5]].z + mesh->vertices_[mesh->cells_[i].vertex_[6]].z) * 1.0 / 2 * 1.0;
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11;
		edge_mid_11.x = (mesh->vertices_[mesh->cells_[i].vertex_[0]].x + mesh->vertices_[mesh->cells_[i].vertex_[3]].x) * 1.0 / 2 * 1.0;
		edge_mid_11.y = (mesh->vertices_[mesh->cells_[i].vertex_[0]].y + mesh->vertices_[mesh->cells_[i].vertex_[3]].y) * 1.0 / 2 * 1.0;
		edge_mid_11.z = (mesh->vertices_[mesh->cells_[i].vertex_[0]].z + mesh->vertices_[mesh->cells_[i].vertex_[3]].z) * 1.0 / 2 * 1.0;
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12;
		edge_mid_12.x = (mesh->vertices_[mesh->cells_[i].vertex_[4]].x + mesh->vertices_[mesh->cells_[i].vertex_[7]].x) * 1.0 / 2 * 1.0;
		edge_mid_12.y = (mesh->vertices_[mesh->cells_[i].vertex_[4]].y + mesh->vertices_[mesh->cells_[i].vertex_[7]].y) * 1.0 / 2 * 1.0;
		edge_mid_12.z = (mesh->vertices_[mesh->cells_[i].vertex_[4]].z + mesh->vertices_[mesh->cells_[i].vertex_[7]].z) * 1.0 / 2 * 1.0;
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);

		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid;
			face_mid.x = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].x
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].x) * 1.0 / 4 * 1.0;
			face_mid.y = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].y
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].y) * 1.0 / 4 * 1.0;
			face_mid.z = (mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[0]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[1]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[2]].z
				+ mesh->vertices_[mesh->faces_[mesh->cells_[i].faces_[mesh_face_num]].vs_f[3]].z) * 1.0 / 4 * 1.0;
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}

		//插入单元的内心
		V3f cell_mid;
		double x_sum = 0;  double y_sum = 0;  double z_sum = 0;
		for (int j = 0; j < 8; j++)
		{
			x_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].x;
			y_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].y;
			z_sum += mesh->vertices_[mesh->cells_[i].vertex_[j]].z;
		}
		cell_mid.x = x_sum * 1.0 / 8 * 1.0;
		cell_mid.y = y_sum * 1.0 / 8 * 1.0;
		cell_mid.z = z_sum * 1.0 / 8 * 1.0;
		temp_vex[26] = newmesh->add_vertex(cell_mid);

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}

HexV3fMesh& MLCASubdivision::Averaging()
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh->vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		if (mesh->is_bdy(VertexHandle(vex_num)))
		{
			std::set<FaceHandle, compare_OVM> neighbor_vf = mesh->neighbor_v_f[vex_num];
			std::vector<FaceHandle> bdy_face;
			std::set<FaceHandle, compare_OVM>::iterator vf_id;
			for (vf_id = neighbor_vf.begin(); vf_id != neighbor_vf.end(); vf_id++)
			{
				if (mesh->is_bdy(*vf_id))
				{
					bdy_face.push_back(*vf_id);
				}
			}
			int vf_num = bdy_face.size();
			std::vector<V3f> vf_cen(vf_num);
			std::vector<FaceHandle>::iterator bdy_face_id;
			int u = 0;
			for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++)
			{
				std::vector<VertexHandle> vf = mesh->faces_[*bdy_face_id].vs_f;
				vf_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x) / 4 * 1.0;
				vf_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y) / 4 * 1.0;
				vf_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z) / 4 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex.x += vf_cen[vf_cen_num].x;
				new_vex.y += vf_cen[vf_cen_num].y;
				new_vex.z += vf_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num * 1.0;
			new_vex.y /= vf_num * 1.0;
			new_vex.z /= vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
		else
		{
			std::set<CellHandle, compare_OVM> neighbor_vc = mesh->neighbor_v_c[vex_num];
			std::set<CellHandle, compare_OVM>::iterator vc_id;
			int vf_num = neighbor_vc.size();
			std::vector<V3f> vc_cen(vf_num);
			int u = 0;
			for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++)
			{
				std::vector<VertexHandle> vf = mesh->cells_[*vc_id].vertex_;
				vc_cen[u].x = (mesh->vertices_[vf[0]].x + mesh->vertices_[vf[1]].x
					+ mesh->vertices_[vf[2]].x + mesh->vertices_[vf[3]].x +
					mesh->vertices_[vf[4]].x + mesh->vertices_[vf[5]].x +
					mesh->vertices_[vf[6]].x + mesh->vertices_[vf[7]].x) / 8 * 1.0;
				vc_cen[u].y = (mesh->vertices_[vf[0]].y + mesh->vertices_[vf[1]].y
					+ mesh->vertices_[vf[2]].y + mesh->vertices_[vf[3]].y +
					mesh->vertices_[vf[4]].y + mesh->vertices_[vf[5]].y +
					mesh->vertices_[vf[6]].y + mesh->vertices_[vf[7]].y) / 8 * 1.0;
				vc_cen[u].z = (mesh->vertices_[vf[0]].z + mesh->vertices_[vf[1]].z
					+ mesh->vertices_[vf[2]].z + mesh->vertices_[vf[3]].z +
					mesh->vertices_[vf[4]].z + mesh->vertices_[vf[5]].z +
					mesh->vertices_[vf[6]].z + mesh->vertices_[vf[7]].z) / 8 * 1.0;
				u++;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex.x += vc_cen[vf_cen_num].x;
				new_vex.y += vc_cen[vf_cen_num].y;
				new_vex.z += vc_cen[vf_cen_num].z;
			}
			new_vex.x /= vf_num * 1.0;
			new_vex.y /= vf_num * 1.0;
			new_vex.z /= vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
	}
	for (uint64_t vex_num = 0; vex_num < mesh->vertices_.size(); vex_num++)
	{
		mesh->vertices_[vex_num] = changed_vex[vex_num];
	}
	return *mesh;
}

//-----ApproximationCC-----//
V3f ApproximationCC::cal_corner_point(VertexHandle _v)
{
	V3f newv;
	int n = 0;
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_bdy(*ve_it))
			n++;
	}
	newv = newv + mesh->vertices_[_v] * (n * n);
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_bdy(*ve_it))
		{
			if (_v == mesh->edges_[*ve_it].from_h)
				newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 4;
			else
				newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 4;
		}
	}
	std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
	for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin(); vf_it != vf.end(); vf_it++)
	{
		if (mesh->is_bdy(*vf_it))
		{
			V3f opv = mesh->vertices_[mesh->find_op_vertex(_v, *vf_it)];
			newv = newv + opv;
		}
	}
	newv = newv / (n * n + 5 * n);
	return newv;
}

//V3f ApproximationCC::cal_corner_point_1(VertexHandle _v,CellHandle _c)
//{
//	V3f newv;
//	int n = 0;
//	int bdy_edge = 0;
//	int not_bdy_edge = 0;
//	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
//	std::set<CellHandle> c_c_temp = mesh->neighbor_c_c(_c);
//	c_c_temp.insert(_c);
//	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//	{
//		if (mesh->is_bdy(*ve_it))
//		{
//			n++;
//			for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin();c_c_temp_it != c_c_temp.end();c_c_temp_it++)
//			{
//				if (mesh->is_e_in_c(*ve_it, *c_c_temp_it))
//				{
//					bdy_edge++;
//					break;
//				}
//			}	
//		}
//		else
//			not_bdy_edge++;	
//	}
//	if (n > 5 || not_bdy_edge < 1)
//	{
//		newv = newv + mesh->vertices_[_v] * (bdy_edge * bdy_edge);
//		for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//		{
//			for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin();c_c_temp_it != c_c_temp.end();c_c_temp_it++)
//			{
//				if (mesh->is_bdy(*ve_it) && mesh->is_e_in_c(*ve_it, *c_c_temp_it))
//				{
//					if (_v == mesh->edges_[*ve_it].from_h)
//					{
//						newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 4;
//						break;
//					}
//					else
//					{
//						newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 4;
//						break;
//					}		
//				}
//			}	
//		}
//		std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
//		for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin();vf_it != vf.end();vf_it++)
//		{
//			for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin();c_c_temp_it != c_c_temp.end();c_c_temp_it++)
//			{
//				if (mesh->is_bdy(*vf_it) && mesh->is_f_in_c(*vf_it, *c_c_temp_it))
//				{
//					V3f opv = mesh->vertices_[mesh->find_op_vertex(_v, *vf_it)];
//					newv = newv + opv;
//					break;
//				}
//			}		
//		}
//		newv = newv / (bdy_edge * bdy_edge + 5 * bdy_edge);
//		return newv;
//	}
//	else
//	{
//		newv = newv + mesh->vertices_[_v] * (n * n);
//		for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//		{
//			if (mesh->is_bdy(*ve_it))
//			{
//				if (_v == mesh->edges_[*ve_it].from_h)
//					newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 4;
//				else
//					newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 4;
//			}
//		}
//		std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
//		for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin();vf_it != vf.end();vf_it++)
//		{
//			if (mesh->is_bdy(*vf_it))
//			{
//				V3f opv = mesh->vertices_[mesh->find_op_vertex(_v, *vf_it)];
//				newv = newv + opv;
//			}
//		}
//		newv = newv / (n * n + 5 * n);
//		return newv;
//	}	
//}

V3f ApproximationCC::cal_edge_point(EdgeHandle _e, VertexHandle _v)
{
	V3f newv;
	int n = 0;
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_bdy(*ve_it))
			n++;
	}
	newv = newv + mesh->vertices_[_v] * (2 * n);
	std::set<EdgeHandle, compare_OVM> opve;
	VertexHandle opv;
	if (mesh->edges_[_e].from_h == _v)
	{
		opv = mesh->edges_[_e].to_h;
		newv = newv + mesh->vertices_[opv] * 4;
		opve = mesh->neighbor_v[opv];
	}
	else
	{
		opv = mesh->edges_[_e].from_h;
		newv = newv + mesh->vertices_[opv] * 4;
		opve = mesh->neighbor_v[opv];
	}
	std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
	int m = 0;
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_bdy(*ve_it) && (*ve_it) != _e)
		{
			for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin(); ef_it != ef.end(); ef_it++)
			{
				if (mesh->is_bdy(*ef_it))
				{
					if (mesh->is_e_in_f(*ve_it, *ef_it))
					{
						m++;
						if (mesh->edges_[*ve_it].from_h == _v)
						{
							newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 2;
							break;
						}
						else
						{
							newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 2;
							break;
						}
					}
				}
			}
		}
	}
	for (std::set<EdgeHandle, compare_OVM>::iterator opve_it = opve.begin(); opve_it != opve.end(); opve_it++)
	{
		if (mesh->is_bdy(*opve_it) && (*opve_it) != _e)
		{
			for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin(); ef_it != ef.end(); ef_it++)
			{
				if (mesh->is_bdy(*ef_it))
				{
					if (mesh->is_e_in_f(*opve_it, *ef_it))
					{
						if (mesh->edges_[*opve_it].from_h == opv)
						{
							newv = newv + mesh->vertices_[mesh->edges_[*opve_it].to_h];
							break;
						}
						else
						{
							newv = newv + mesh->vertices_[mesh->edges_[*opve_it].from_h];
							break;
						}
					}
				}
			}
		}
	}
	newv = newv / (2 * n + 4 + 3 * m);
	return newv;
}

//V3f ApproximationCC::cal_edge_point_1(EdgeHandle _e, VertexHandle _v, CellHandle _c)
//{
//	V3f newv;
//	int n = 0;
//	int bdy_edge = 0;
//	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
//	std::set<CellHandle> c_c_temp = mesh->neighbor_c_c(_c);
//	c_c_temp.insert(_c);
//	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//	{
//		if (mesh->is_bdy(*ve_it))
//		{
//			n++;
//			for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin();c_c_temp_it != c_c_temp.end();c_c_temp_it++)
//			{
//				if (mesh->is_e_in_c(*ve_it, *c_c_temp_it))
//				{
//					bdy_edge++;
//					break;
//				}
//			}
//		}
//	}
//	std::set<FaceHandle, compare_OVM> e_f_temp = mesh->neighbor_e[_e];
//	int ef_temp = 0;
//	for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f_temp.begin();e_f_it != e_f_temp.end();e_f_it++)
//	{
//		if (mesh->is_bdy(*e_f_it))
//			ef_temp++;
//	}
//	if (ef_temp == 4)
//	{
//		newv = newv + mesh->vertices_[_v] * (2 * bdy_edge);
//		std::set<EdgeHandle, compare_OVM> opve;
//		VertexHandle opv;
//		if (mesh->edges_[_e].from_h == _v)
//		{
//			opv = mesh->edges_[_e].to_h;
//			newv = newv + mesh->vertices_[opv] * 4;
//			opve = mesh->neighbor_v[opv];
//		}
//		else
//		{
//			opv = mesh->edges_[_e].from_h;
//			newv = newv + mesh->vertices_[opv] * 4;
//			opve = mesh->neighbor_v[opv];
//		}
//		std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
//		int m = 0;
//		for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//		{
//			if (mesh->is_bdy(*ve_it) && (*ve_it) != _e && mesh->is_e_in_c(*ve_it,_c))
//			{
//				for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin();ef_it != ef.end();ef_it++)
//				{
//					if (mesh->is_bdy(*ef_it))
//					{
//						if (mesh->is_e_in_f(*ve_it, *ef_it))
//						{
//							m++;
//							if (mesh->edges_[*ve_it].from_h == _v)
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 2;
//								break;
//							}
//							else
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 2;
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//		for (std::set<EdgeHandle, compare_OVM>::iterator opve_it = opve.begin();opve_it != opve.end();opve_it++)
//		{
//			if (mesh->is_bdy(*opve_it) && (*opve_it) != _e&& mesh->is_e_in_c(*opve_it, _c))
//			{
//				for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin();ef_it != ef.end();ef_it++)
//				{
//					if (mesh->is_bdy(*ef_it))
//					{
//						if (mesh->is_e_in_f(*opve_it, *ef_it))
//						{
//							if (mesh->edges_[*opve_it].from_h == opv)
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*opve_it].to_h];
//								break;
//							}
//							else
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*opve_it].from_h];
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//		newv = newv / (2 * bdy_edge + 4 + 3 * m);
//		return newv;
//	}
//	else
//	{
//		newv = newv + mesh->vertices_[_v] * (2 * n);
//		std::set<EdgeHandle, compare_OVM> opve;
//		VertexHandle opv;
//		if (mesh->edges_[_e].from_h == _v)
//		{
//			opv = mesh->edges_[_e].to_h;
//			newv = newv + mesh->vertices_[opv] * 4;
//			opve = mesh->neighbor_v[opv];
//		}
//		else
//		{
//			opv = mesh->edges_[_e].from_h;
//			newv = newv + mesh->vertices_[opv] * 4;
//			opve = mesh->neighbor_v[opv];
//		}
//		std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
//		int m = 0;
//		for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin();ve_it != ve.end();ve_it++)
//		{
//			if (mesh->is_bdy(*ve_it) && (*ve_it) != _e)
//			{
//				for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin();ef_it != ef.end();ef_it++)
//				{
//					if (mesh->is_bdy(*ef_it))
//					{
//						if (mesh->is_e_in_f(*ve_it, *ef_it))
//						{
//							m++;
//							if (mesh->edges_[*ve_it].from_h == _v)
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 2;
//								break;
//							}
//							else
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 2;
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//		for (std::set<EdgeHandle, compare_OVM>::iterator opve_it = opve.begin();opve_it != opve.end();opve_it++)
//		{
//			if (mesh->is_bdy(*opve_it) && (*opve_it) != _e)
//			{
//				for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin();ef_it != ef.end();ef_it++)
//				{
//					if (mesh->is_bdy(*ef_it))
//					{
//						if (mesh->is_e_in_f(*opve_it, *ef_it))
//						{
//							if (mesh->edges_[*opve_it].from_h == opv)
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*opve_it].to_h];
//								break;
//							}
//							else
//							{
//								newv = newv + mesh->vertices_[mesh->edges_[*opve_it].from_h];
//								break;
//							}
//						}
//					}
//				}
//			}
//		}
//		newv = newv / (2 * n + 4 + 3 * m);
//		return newv;
//	}	
//}

V3f ApproximationCC::cal_inter_point(FaceHandle _f, VertexHandle _v)
{
	V3f newv;
	int n = 0;
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_bdy(*ve_it))
			n++;
	}
	newv = newv + mesh->vertices_[_v] * n;
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
	newv = newv / (n + 5);
	return newv;
}

V3f ApproximationCC::cal_inter_point(CellHandle _c, VertexHandle _v)
{
	return cal_inter_point1(_c, _v);
}

void ApproximationCC::approximationing()
{
	int face_num = mesh->faces_.size();
	/*std::vector<V3f> check;
	std::vector<V3f> check1;*/
	for (int i = 0; i < face_num; i++)//以面为单位遍历
	{
		if (mesh->is_bdy(FaceHandle(i)))
		{
			std::vector<V3f> temp;
			std::vector<EdgeHandle> fe = mesh->faces_[i].edges_;
			//VertexHandle markv = mesh->edges_[*(fe.begin())].from_h;//确定方向
			VertexHandle fev, feopv;
			V3f cornerv, interv, epv, epopv;
			//边排序
			/*第1条边*/
			fev = mesh->edges_[fe[0]].from_h; feopv = mesh->edges_[fe[0]].to_h;
			cornerv = cal_corner_point(fev);
			interv = cal_inter_point(FaceHandle(i), fev);
			epv = cal_edge_point(fe[0], fev);
			epopv = cal_edge_point(fe[0], feopv);
			temp.push_back(cornerv);
			temp.push_back(interv);
			temp.push_back(epv);
			temp.push_back(epopv);
			std::vector<EdgeHandle> tempev;
			EdgeHandle tempe = fe[0];
			tempev.push_back(tempe);
			/*后3条边*/
			for (int laternum = 0; laternum < 3; laternum++)
			{
				for (std::vector<EdgeHandle>::iterator fe_it = fe.begin(); fe_it != fe.end(); fe_it++)
				{
					if (find(tempev.begin(), tempev.end(), *fe_it) == tempev.end())
					{
						if (mesh->edges_[*fe_it].from_h == feopv)
						{
							tempe = *fe_it;
							tempev.push_back(tempe);
							fev = mesh->edges_[*fe_it].from_h;
							feopv = mesh->edges_[*fe_it].to_h;
							break;
						}
						else if (mesh->edges_[*fe_it].to_h == feopv)
						{
							tempe = *fe_it;
							tempev.push_back(tempe);
							fev = mesh->edges_[*fe_it].to_h;
							feopv = mesh->edges_[*fe_it].from_h;
							break;
						}
						else
							continue;
					}
				}
				cornerv = cal_corner_point(fev);
				interv = cal_inter_point(FaceHandle(i), fev);
				epv = cal_edge_point(tempe, fev);
				epopv = cal_edge_point(tempe, feopv);
				temp.push_back(cornerv);
				temp.push_back(interv);
				temp.push_back(epv);
				temp.push_back(epopv);
			}

			//for (std::vector<EdgeHandle>::iterator fe_it = fe.begin(); fe_it != fe.end();fe_it++)//以边来遍历这16个控制点
			//{
			//	VertexHandle frv = mesh->edges_[*fe_it].from_h;
			//	VertexHandle tov = mesh->edges_[*fe_it].to_h;
			//	if (markv == frv)
			//	{
			//		fev = frv;
			//		feopv = tov;
			//		markv = tov;
			//	}
			//	else
			//	{
			//		fev = tov;
			//		feopv = frv;
			//		markv = frv;
			//	}
			//	V3f cornerv = cal_corner_point(fev);
			//	V3f interv = cal_inter_point(FaceHandle(i), fev);
			//	V3f epv = cal_edge_point(*fe_it, fev);
			//	V3f epopv = cal_edge_point(*fe_it, feopv);
			//	temp.push_back(cornerv);
			//	temp.push_back(interv);
			//	temp.push_back(epv);
			//	temp.push_back(epopv);
			//	if (find(check.begin(), check.end(), epv) == check.end())
			//		check.push_back(epv);
			//	check1.push_back(epv);
			//}
			/*int tempsize = temp.size();
			for (int i = 0; i < tempsize;i++)
			{
				if (find(check.begin(), check.end(), temp[i]) == check.end())
					check.push_back(temp[i]);
				check1.push_back(temp[i]);
			}*/
			std::vector<V3f> f1, f2, f3, f4, f5, f6, f7, f8, f9;
			f1.push_back(temp[0]); f1.push_back(temp[2]); f1.push_back(temp[1]); f1.push_back(temp[15]);
			f2.push_back(temp[2]); f2.push_back(temp[3]); f2.push_back(temp[5]); f2.push_back(temp[1]);
			f3.push_back(temp[3]); f3.push_back(temp[4]); f3.push_back(temp[6]); f3.push_back(temp[5]);
			f4.push_back(temp[15]); f4.push_back(temp[1]); f4.push_back(temp[13]); f4.push_back(temp[14]);
			f5.push_back(temp[1]); f5.push_back(temp[5]); f5.push_back(temp[9]); f5.push_back(temp[13]);
			f6.push_back(temp[5]); f6.push_back(temp[6]); f6.push_back(temp[7]); f6.push_back(temp[9]);
			f7.push_back(temp[14]); f7.push_back(temp[13]); f7.push_back(temp[11]); f7.push_back(temp[12]);
			f8.push_back(temp[13]); f8.push_back(temp[9]); f8.push_back(temp[10]); f8.push_back(temp[11]);
			f9.push_back(temp[9]); f9.push_back(temp[7]); f9.push_back(temp[8]); f9.push_back(temp[10]);
			surv.push_back(f1); surv.push_back(f2); surv.push_back(f3); surv.push_back(f4); surv.push_back(f5);
			surv.push_back(f6); surv.push_back(f7); surv.push_back(f8); surv.push_back(f9);
			std::vector<V3f> beiv1, beiv2, beiv3, beiv4;
			std::vector<std::vector<V3f>> beif;
			beiv1.push_back(temp[0]); beiv1.push_back(temp[15]); beiv1.push_back(temp[14]); beiv1.push_back(temp[12]);
			beiv2.push_back(temp[2]); beiv2.push_back(temp[1]); beiv2.push_back(temp[13]); beiv2.push_back(temp[11]);
			beiv3.push_back(temp[3]); beiv3.push_back(temp[5]); beiv3.push_back(temp[9]); beiv3.push_back(temp[10]);
			beiv4.push_back(temp[4]); beiv4.push_back(temp[6]); beiv4.push_back(temp[7]); beiv4.push_back(temp[8]);
			beif.push_back(beiv1); beif.push_back(beiv2); beif.push_back(beiv3); beif.push_back(beiv4);
			beiv.push_back(beif);
		}
	}
	/*std::vector<V3f> check;
	for (int i = 0; i < surv.size(); ++i)
	{
		for (int j = 0;j < surv[i].size();++j)
		{
			if (find(check.begin(), check.end(), surv[i][j]) == check.end())
				check.push_back(surv[i][j]);
		}
	}
	int checknum = check.size();*/
}

void ApproximationCC::cal_nor_sur()
{
	int sursize = surv.size();
	for (int i = 0; i < sursize; i++)
	{
		V3f a; V3f b; V3f n;
		a = surv[i][1] - surv[i][0];
		b = surv[i][3] - surv[i][0];
		n = a % b;
		n = n.normalize();
		nor_surv.push_back(n);
	}
}

void ApproximationCC::cal_nor_beisur()
{
	int beisursize = bei_surv.size();
	for (int i = 0; i < beisursize; i++)
	{
		V3f a; V3f b; V3f n;
		a = bei_surv[i][1] - bei_surv[i][0];
		b = bei_surv[i][3] - bei_surv[i][0];
		n = a % b;
		n = n.normalize();
		nor_bei_surv.push_back(n);
	}
}

//void ApproximationCC::constructBezier()
//{
//	for (std::vector<std::vector<std::vector<V3f>>>::iterator beiv_it = beiv.begin();beiv_it != beiv.end();beiv_it++)    
//	{
//		std::vector<std::vector<V3f>> f_beiv = *beiv_it;
//		Wm4::Vector3<double> **bv;
//		bv = WM4_NEW Wm4::Vector3<double>*[4];
//		bv[0] = WM4_NEW Wm4::Vector3<double>[16];
//		for (int i = 1; i < 4; i++)
//		{
//			bv[i] = &bv[0][i*4];
//		}
//		/*int i = 0, j = 0;
//		for (std::vector<std::vector<V3f>>::iterator f_beiv_it = f_beiv.begin();f_beiv_it != f_beiv.end();f_beiv_it++)
//		{
//			std::vector<V3f> f_beiv_v = *f_beiv_it;
//			for (std::vector<V3f>::iterator f_beiv_v_it = f_beiv_v.begin();f_beiv_v_it != f_beiv_v.end();f_beiv_v_it++)
//			{
//				bv[i][j] = Wm4::Vector3<double>((*f_beiv_v_it).x, (*f_beiv_v_it).y, (*f_beiv_v_it).x);
//				j++;
//			}
//			i++;
//		}*/
//		bv[0][0] = Wm4::Vector3<double>(f_beiv[0][0].x, f_beiv[0][0].y, f_beiv[0][0].z);
//		bv[0][1] = Wm4::Vector3<double>(f_beiv[0][1].x, f_beiv[0][1].y, f_beiv[0][1].z);
//		bv[0][2] = Wm4::Vector3<double>(f_beiv[0][2].x, f_beiv[0][2].y, f_beiv[0][2].z);
//		bv[0][3] = Wm4::Vector3<double>(f_beiv[0][3].x, f_beiv[0][3].y, f_beiv[0][3].z);
//		bv[1][0] = Wm4::Vector3<double>(f_beiv[1][0].x, f_beiv[1][0].y, f_beiv[1][0].z);
//		bv[1][1] = Wm4::Vector3<double>(f_beiv[1][1].x, f_beiv[1][1].y, f_beiv[1][1].z);
//		bv[1][2] = Wm4::Vector3<double>(f_beiv[1][2].x, f_beiv[1][2].y, f_beiv[1][2].z);
//		bv[1][3] = Wm4::Vector3<double>(f_beiv[1][3].x, f_beiv[1][3].y, f_beiv[1][3].z);
//		bv[2][0] = Wm4::Vector3<double>(f_beiv[2][0].x, f_beiv[2][0].y, f_beiv[2][0].z);
//		bv[2][1] = Wm4::Vector3<double>(f_beiv[2][1].x, f_beiv[2][1].y, f_beiv[2][1].z);
//		bv[2][2] = Wm4::Vector3<double>(f_beiv[2][2].x, f_beiv[2][2].y, f_beiv[2][2].z);
//		bv[2][3] = Wm4::Vector3<double>(f_beiv[2][3].x, f_beiv[2][3].y, f_beiv[2][3].z);
//		bv[3][0] = Wm4::Vector3<double>(f_beiv[3][0].x, f_beiv[3][0].y, f_beiv[3][0].z);
//		bv[3][1] = Wm4::Vector3<double>(f_beiv[3][1].x, f_beiv[3][1].y, f_beiv[3][1].z);
//		bv[3][2] = Wm4::Vector3<double>(f_beiv[3][2].x, f_beiv[3][2].y, f_beiv[3][2].z);
//		bv[3][3] = Wm4::Vector3<double>(f_beiv[3][3].x, f_beiv[3][3].y, f_beiv[3][3].z);
//		Wm4::BSplineRectangle<double> tempbv(4, 4, bv, 3, 3, false, false, true, true);
//		int size = 20;//要生成size * size个面的控制网格 即(size +1) * (size + 1)个控制顶点
//		std::vector<std::vector<V3f>> f(size * size);
//		std::vector<V3f> bezierv((size +1) * (size + 1));
//		for (int i = 0;i < bezierv.size();i++)
//		{
//			bezierv[i] = tempbv.P(i % (size + 1) * (1.0 / size), i / (size + 1) * (1.0 / size));
//		}
//		for (int i = 0;i < f.size();i++)
//		{
//			int j = i + i / size;
//			f[i].push_back(bezierv[j]);f[i].push_back(bezierv[j+1]);
//			f[i].push_back(bezierv[j+size+2]);f[i].push_back(bezierv[j+size+1]);
//			bei_surv.push_back(f[i]);
//		}
//		/*bezierv[0] = tempbv.P(0, 0);bezierv[1] = tempbv.P(0.25, 0);bezierv[2] = tempbv.P(0.5, 0);bezierv[3] = tempbv.P(0.75, 0);bezierv[4] = tempbv.P(1, 0);
//		bezierv[5] = tempbv.P(0, 0.25);bezierv[6] = tempbv.P(0.25, 0.25);bezierv[7] = tempbv.P(0.5, 0.25);bezierv[8] = tempbv.P(0.75, 0.25);bezierv[9] = tempbv.P(1, 0.25);
//		bezierv[10] = tempbv.P(0, 0.5);bezierv[11] = tempbv.P(0.25, 0.5);bezierv[12] = tempbv.P(0.5, 0.5);bezierv[13] = tempbv.P(0.75, 0.5);bezierv[14] = tempbv.P(1, 0.5);
//		bezierv[15] = tempbv.P(0, 0.75);bezierv[16] = tempbv.P(0.25, 0.75);bezierv[17] = tempbv.P(0.5, 0.75);bezierv[18] = tempbv.P(0.75, 0.75);bezierv[19] = tempbv.P(1, 0.75);
//		bezierv[20] = tempbv.P(0, 1);bezierv[21] = tempbv.P(0.25, 1);bezierv[22] = tempbv.P(0.5, 1);bezierv[23] = tempbv.P(0.75, 1);bezierv[24] = tempbv.P(1, 1);*/
//		/*f[0].push_back(bezierv[0]);f[0].push_back(bezierv[1]);f[0].push_back(bezierv[6]);f[0].push_back(bezierv[5]);
//		f[1].push_back(bezierv[1]);f[1].push_back(bezierv[2]);f[1].push_back(bezierv[7]);f[1].push_back(bezierv[6]);
//		f[2].push_back(bezierv[2]);f[2].push_back(bezierv[3]);f[2].push_back(bezierv[8]);f[2].push_back(bezierv[7]);
//		f[3].push_back(bezierv[3]);f[3].push_back(bezierv[4]);f[3].push_back(bezierv[9]);f[3].push_back(bezierv[8]);
//		f[4].push_back(bezierv[5]);f[4].push_back(bezierv[6]);f[4].push_back(bezierv[11]);f[4].push_back(bezierv[10]);
//		f[5].push_back(bezierv[6]);f[5].push_back(bezierv[7]);f[5].push_back(bezierv[12]);f[5].push_back(bezierv[11]);
//		f[6].push_back(bezierv[7]);f[6].push_back(bezierv[8]);f[6].push_back(bezierv[13]);f[6].push_back(bezierv[12]);
//		f[7].push_back(bezierv[8]);f[7].push_back(bezierv[9]);f[7].push_back(bezierv[14]);f[7].push_back(bezierv[13]);
//		f[8].push_back(bezierv[10]);f[8].push_back(bezierv[11]);f[8].push_back(bezierv[16]);f[8].push_back(bezierv[15]);
//		f[9].push_back(bezierv[11]);f[9].push_back(bezierv[12]);f[9].push_back(bezierv[17]);f[9].push_back(bezierv[16]);
//		f[10].push_back(bezierv[12]);f[10].push_back(bezierv[13]);f[10].push_back(bezierv[18]);f[10].push_back(bezierv[17]);
//		f[11].push_back(bezierv[13]);f[11].push_back(bezierv[14]);f[11].push_back(bezierv[19]);f[11].push_back(bezierv[18]);
//		f[12].push_back(bezierv[15]);f[12].push_back(bezierv[16]);f[12].push_back(bezierv[21]);f[12].push_back(bezierv[20]);
//		f[13].push_back(bezierv[16]);f[13].push_back(bezierv[17]);f[13].push_back(bezierv[22]);f[13].push_back(bezierv[21]);
//		f[14].push_back(bezierv[17]);f[14].push_back(bezierv[18]);f[14].push_back(bezierv[23]);f[14].push_back(bezierv[22]);
//		f[15].push_back(bezierv[18]);f[15].push_back(bezierv[19]);f[15].push_back(bezierv[24]);f[15].push_back(bezierv[23]);*/
//		/*bei_surv.push_back(f[0]);bei_surv.push_back(f[4]);bei_surv.push_back(f[8]);bei_surv.push_back(f[12]);
//		bei_surv.push_back(f[1]);bei_surv.push_back(f[5]);bei_surv.push_back(f[9]);bei_surv.push_back(f[13]);
//		bei_surv.push_back(f[2]);bei_surv.push_back(f[6]);bei_surv.push_back(f[10]);bei_surv.push_back(f[14]);
//		bei_surv.push_back(f[3]);bei_surv.push_back(f[7]);bei_surv.push_back(f[11]);bei_surv.push_back(f[15]);*/
//	}
//}

V3f ApproximationCC::solid_edge_point(EdgeHandle _e, VertexHandle _v)
{
	if (mesh->is_bdy(_e))
	{
		// V3f edp;
		// std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
		// int n = 0;
		// for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin();ef_it != ef.end();ef_it++)
		// {
			// if (mesh->is_bdy(*ef_it))
			// {
				// edp = edp + cal_inter_point(*ef_it, _v);
				// n++;
			// }
		// }
		// edp = edp / n;
		// return edp;
		return cal_edge_pointsharp(_e, _v);
	}
	else
	{
		V3f edp;
		std::set<CellHandle> ec = mesh->neighbor_e_c(_e);
		int n = ec.size();
		for (std::set<CellHandle>::iterator ec_it = ec.begin(); ec_it != ec.end(); ec_it++)
		{
			edp = edp + cal_inter_point(*ec_it, _v);
		}
		edp = edp / n;
		return edp;
	}
}

V3f ApproximationCC::solid_edge_point_1(EdgeHandle _e, VertexHandle _v, CellHandle _c)
{
	if (mesh->is_bdy(_e))
	{
		std::set<FaceHandle, compare_OVM> e_f_temp = mesh->neighbor_e[_e];
		int ef_temp = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f_temp.begin(); e_f_it != e_f_temp.end(); e_f_it++)
		{
			if (mesh->is_bdy(*e_f_it))
				ef_temp++;
		}
		if (ef_temp == 4)
		{
			V3f edp;
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
			int n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin(); ef_it != ef.end(); ef_it++)
			{
				if (mesh->is_bdy(*ef_it) && mesh->is_f_in_c(*ef_it, _c))
				{
					edp = edp + cal_inter_point(*ef_it, _v);
					n++;
				}
			}
			edp = edp / n;
			return edp;
		}
		else
		{
			V3f edp;
			std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
			int n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin(); ef_it != ef.end(); ef_it++)
			{
				if (mesh->is_bdy(*ef_it))
				{
					edp = edp + cal_inter_point(*ef_it, _v);
					n++;
				}
			}
			edp = edp / n;
			return edp;
		}
	}
	else
	{
		V3f edp;
		std::set<CellHandle> ec = mesh->neighbor_e_c(_e);
		int n = ec.size();
		for (std::set<CellHandle>::iterator ec_it = ec.begin(); ec_it != ec.end(); ec_it++)
		{
			edp = edp + cal_inter_point(*ec_it, _v);
		}
		edp = edp / n;
		return edp;
	}
}

V3f ApproximationCC::solid_face_point(FaceHandle _f, VertexHandle _v)
{
	if (mesh->is_bdy(_f))
	{
		// V3f fp;
		// fp = cal_inter_point(_f, _v);
		// return fp;
		return cal_inter_pointsharp(_f, _v);
	}
	else
	{
		V3f fp;
		std::set<CellHandle, compare_OVM> fc = mesh->neighbor_f[_f];
		for (std::set<CellHandle, compare_OVM>::iterator fc_it = fc.begin(); fc_it != fc.end(); fc_it++)
		{
			fp = fp + cal_inter_point(*fc_it, _v);
		}
		fp = fp / 2;
		return fp;
	}
}

V3f ApproximationCC::solid_corner_point(VertexHandle _v)
{
	if (mesh->is_bdy(_v))
	{
		// V3f cp;
		// std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
		// int n = 0;
		// for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin();vf_it != vf.end();vf_it++)
		// {
			// if (mesh->is_bdy(*vf_it))
			// {
				// cp = cp + cal_inter_point(*vf_it, _v);
				// n++;
			// }
		// }
		// cp = cp / n;
		// return cp;
		return cal_corner_pointsharp(_v);
	}
	else
	{
		V3f cp;
		std::set<CellHandle, compare_OVM> vc = mesh->neighbor_v_c[_v];
		int n = vc.size();
		for (std::set<CellHandle, compare_OVM>::iterator vc_it = vc.begin(); vc_it != vc.end(); vc_it++)
		{
			cp = cp + cal_inter_point(*vc_it, _v);
		}
		cp = cp / n;
		return cp;
	}
}

V3f ApproximationCC::solid_corner_point_1(VertexHandle _v, CellHandle _c)
{
	if (mesh->is_bdy(_v))
	{
		std::set<EdgeHandle, compare_OVM> v_e_temp = mesh->neighbor_v[_v];
		int bdy_edge = 0;
		int not_bdy_edge = 0;
		for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e_temp.begin(); v_e_it != v_e_temp.end(); v_e_it++)
		{
			if (mesh->is_bdy(*v_e_it))
				bdy_edge++;
			else
				not_bdy_edge++;
		}
		std::set<CellHandle> c_c_temp = mesh->neighbor_c_c(_c);
		c_c_temp.insert(_c);
		if (bdy_edge > 5 || not_bdy_edge < 1)
		{
			V3f cp;
			std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
			int n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin(); vf_it != vf.end(); vf_it++)
			{
				if (mesh->is_bdy(*vf_it))
				{
					for (std::set<CellHandle>::iterator c_c_temp_it = c_c_temp.begin(); c_c_temp_it != c_c_temp.end(); c_c_temp_it++)
					{
						if (mesh->is_f_in_c(*vf_it, *c_c_temp_it))
						{
							cp = cp + cal_inter_point(*vf_it, _v);
							n++;
							break;
						}
					}
				}
			}
			cp = cp / n;
			return cp;
		}
		else
		{
			V3f cp;
			std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
			int n = 0;
			for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin(); vf_it != vf.end(); vf_it++)
			{
				if (mesh->is_bdy(*vf_it))
				{
					cp = cp + cal_inter_point(*vf_it, _v);
					n++;
				}
			}
			cp = cp / n;
			return cp;
		}


	}
	else
	{
		V3f cp;
		std::set<CellHandle, compare_OVM> vc = mesh->neighbor_v_c[_v];
		int n = vc.size();
		for (std::set<CellHandle, compare_OVM>::iterator vc_it = vc.begin(); vc_it != vc.end(); vc_it++)
		{
			cp = cp + cal_inter_point(*vc_it, _v);
		}
		cp = cp / n;
		return cp;
	}
}

V3f ApproximationCC::solid_edge_point(int cell_num, int edge_num, int vertex_num)
{
	return solid_edge_point(mesh->cells_[cell_num].edges_[edge_num], mesh->cells_[cell_num].vertex_[vertex_num]);
}

V3f ApproximationCC::solid_edge_point_1(int cell_num, int edge_num, int vertex_num)
{
	return solid_edge_point_1(mesh->cells_[cell_num].edges_[edge_num], mesh->cells_[cell_num].vertex_[vertex_num], CellHandle(cell_num));
}

V3f ApproximationCC::solid_face_point(int cell_num, int face_num, int vertex_num)
{
	return solid_face_point(mesh->cells_[cell_num].faces_[face_num], mesh->cells_[cell_num].vertex_[vertex_num]);
}

HexV3fMesh& ApproximationCC::CCapproximation()
{
	uint64_t cell_num = mesh->cells_.size();
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(solid_corner_point(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 9, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 9, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 1, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 7, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 1, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 7, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 8, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 8, 2));
		//第二层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 6, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 4, 2));
		//第三层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 2, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 6, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 4, 3));
		//第4层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 11, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 11, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 5, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 5, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 10, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 10, 3));

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(27);
		/*第一层*/
		//第1个单元
		sub_cell[0].push_back(temp[5]); sub_cell[0].push_back(temp[16]); sub_cell[0].push_back(temp[19]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[28]); sub_cell[0].push_back(temp[29]); sub_cell[0].push_back(temp[13]); sub_cell[0].push_back(temp[32]);
		//第2个单元
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[17]); sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[19]);
		sub_cell[1].push_back(temp[29]); sub_cell[1].push_back(temp[30]); sub_cell[1].push_back(temp[14]); sub_cell[1].push_back(temp[13]);
		//第3个单元
		sub_cell[2].push_back(temp[17]); sub_cell[2].push_back(temp[6]); sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[20]);
		sub_cell[2].push_back(temp[30]); sub_cell[2].push_back(temp[31]); sub_cell[2].push_back(temp[33]); sub_cell[2].push_back(temp[14]);
		//第4个单元
		sub_cell[3].push_back(temp[18]); sub_cell[3].push_back(temp[19]); sub_cell[3].push_back(temp[23]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[32]); sub_cell[3].push_back(temp[13]); sub_cell[3].push_back(temp[9]); sub_cell[3].push_back(temp[34]);
		//第5个单元
		sub_cell[4].push_back(temp[19]); sub_cell[4].push_back(temp[20]); sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[13]); sub_cell[4].push_back(temp[14]); sub_cell[4].push_back(temp[10]); sub_cell[4].push_back(temp[9]);
		//第6个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[21]); sub_cell[5].push_back(temp[25]); sub_cell[5].push_back(temp[24]);
		sub_cell[5].push_back(temp[14]); sub_cell[5].push_back(temp[33]); sub_cell[5].push_back(temp[35]); sub_cell[5].push_back(temp[10]);
		//第7个单元
		sub_cell[6].push_back(temp[22]); sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]); sub_cell[6].push_back(temp[1]);
		sub_cell[6].push_back(temp[34]); sub_cell[6].push_back(temp[9]); sub_cell[6].push_back(temp[37]); sub_cell[6].push_back(temp[36]);
		//第8个单元
		sub_cell[7].push_back(temp[23]); sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[27]); sub_cell[7].push_back(temp[26]);
		sub_cell[7].push_back(temp[9]); sub_cell[7].push_back(temp[10]); sub_cell[7].push_back(temp[38]); sub_cell[7].push_back(temp[37]);
		//第9个单元
		sub_cell[8].push_back(temp[24]); sub_cell[8].push_back(temp[25]); sub_cell[8].push_back(temp[2]); sub_cell[8].push_back(temp[27]);
		sub_cell[8].push_back(temp[10]); sub_cell[8].push_back(temp[35]); sub_cell[8].push_back(temp[39]); sub_cell[8].push_back(temp[38]);
		/*第二层*/
		//第10个单元
		sub_cell[9].push_back(temp[28]); sub_cell[9].push_back(temp[29]); sub_cell[9].push_back(temp[13]); sub_cell[9].push_back(temp[32]);
		sub_cell[9].push_back(temp[40]); sub_cell[9].push_back(temp[41]); sub_cell[9].push_back(temp[12]); sub_cell[9].push_back(temp[44]);
		//第11个单元
		sub_cell[10].push_back(temp[29]); sub_cell[10].push_back(temp[30]); sub_cell[10].push_back(temp[14]); sub_cell[10].push_back(temp[13]);
		sub_cell[10].push_back(temp[41]); sub_cell[10].push_back(temp[42]); sub_cell[10].push_back(temp[15]); sub_cell[10].push_back(temp[12]);
		//第12个单元
		sub_cell[11].push_back(temp[30]); sub_cell[11].push_back(temp[31]); sub_cell[11].push_back(temp[33]); sub_cell[11].push_back(temp[14]);
		sub_cell[11].push_back(temp[42]); sub_cell[11].push_back(temp[43]); sub_cell[11].push_back(temp[45]); sub_cell[11].push_back(temp[15]);
		//第13个单元
		sub_cell[12].push_back(temp[32]); sub_cell[12].push_back(temp[13]); sub_cell[12].push_back(temp[9]); sub_cell[12].push_back(temp[34]);
		sub_cell[12].push_back(temp[44]); sub_cell[12].push_back(temp[12]); sub_cell[12].push_back(temp[8]); sub_cell[12].push_back(temp[46]);
		//第14个单元
		sub_cell[13].push_back(temp[13]); sub_cell[13].push_back(temp[14]); sub_cell[13].push_back(temp[10]); sub_cell[13].push_back(temp[9]);
		sub_cell[13].push_back(temp[12]); sub_cell[13].push_back(temp[15]); sub_cell[13].push_back(temp[11]); sub_cell[13].push_back(temp[8]);
		//第15个单元
		sub_cell[14].push_back(temp[14]); sub_cell[14].push_back(temp[33]); sub_cell[14].push_back(temp[35]); sub_cell[14].push_back(temp[10]);
		sub_cell[14].push_back(temp[15]); sub_cell[14].push_back(temp[45]); sub_cell[14].push_back(temp[47]); sub_cell[14].push_back(temp[11]);
		//第16个单元
		sub_cell[15].push_back(temp[34]); sub_cell[15].push_back(temp[9]); sub_cell[15].push_back(temp[37]); sub_cell[15].push_back(temp[36]);
		sub_cell[15].push_back(temp[46]); sub_cell[15].push_back(temp[8]); sub_cell[15].push_back(temp[49]); sub_cell[15].push_back(temp[48]);
		//第17个单元
		sub_cell[16].push_back(temp[9]); sub_cell[16].push_back(temp[10]); sub_cell[16].push_back(temp[38]); sub_cell[16].push_back(temp[37]);
		sub_cell[16].push_back(temp[8]); sub_cell[16].push_back(temp[11]); sub_cell[16].push_back(temp[50]); sub_cell[16].push_back(temp[49]);
		//第18个单元
		sub_cell[17].push_back(temp[10]); sub_cell[17].push_back(temp[35]); sub_cell[17].push_back(temp[39]); sub_cell[17].push_back(temp[38]);
		sub_cell[17].push_back(temp[11]); sub_cell[17].push_back(temp[47]); sub_cell[17].push_back(temp[51]); sub_cell[17].push_back(temp[50]);
		/*第三层*/
		//第19个单元
		sub_cell[18].push_back(temp[40]); sub_cell[18].push_back(temp[41]); sub_cell[18].push_back(temp[12]); sub_cell[18].push_back(temp[44]);
		sub_cell[18].push_back(temp[4]); sub_cell[18].push_back(temp[52]); sub_cell[18].push_back(temp[55]); sub_cell[18].push_back(temp[54]);
		//第20个单元
		sub_cell[19].push_back(temp[41]); sub_cell[19].push_back(temp[42]); sub_cell[19].push_back(temp[15]); sub_cell[19].push_back(temp[12]);
		sub_cell[19].push_back(temp[52]); sub_cell[19].push_back(temp[53]); sub_cell[19].push_back(temp[56]); sub_cell[19].push_back(temp[55]);
		//第21个单元
		sub_cell[20].push_back(temp[42]); sub_cell[20].push_back(temp[43]); sub_cell[20].push_back(temp[45]); sub_cell[20].push_back(temp[15]);
		sub_cell[20].push_back(temp[53]); sub_cell[20].push_back(temp[7]); sub_cell[20].push_back(temp[57]); sub_cell[20].push_back(temp[56]);
		//第22个单元
		sub_cell[21].push_back(temp[44]); sub_cell[21].push_back(temp[12]); sub_cell[21].push_back(temp[8]); sub_cell[21].push_back(temp[46]);
		sub_cell[21].push_back(temp[54]); sub_cell[21].push_back(temp[55]); sub_cell[21].push_back(temp[59]); sub_cell[21].push_back(temp[58]);
		//第23个单元
		sub_cell[22].push_back(temp[12]); sub_cell[22].push_back(temp[15]); sub_cell[22].push_back(temp[11]); sub_cell[22].push_back(temp[8]);
		sub_cell[22].push_back(temp[55]); sub_cell[22].push_back(temp[56]); sub_cell[22].push_back(temp[60]); sub_cell[22].push_back(temp[59]);
		//第24个单元
		sub_cell[23].push_back(temp[15]); sub_cell[23].push_back(temp[45]); sub_cell[23].push_back(temp[47]); sub_cell[23].push_back(temp[11]);
		sub_cell[23].push_back(temp[56]); sub_cell[23].push_back(temp[57]); sub_cell[23].push_back(temp[61]); sub_cell[23].push_back(temp[60]);
		//第25个单元
		sub_cell[24].push_back(temp[46]); sub_cell[24].push_back(temp[8]); sub_cell[24].push_back(temp[49]); sub_cell[24].push_back(temp[48]);
		sub_cell[24].push_back(temp[58]); sub_cell[24].push_back(temp[59]); sub_cell[24].push_back(temp[62]); sub_cell[24].push_back(temp[0]);
		//第26个单元
		sub_cell[25].push_back(temp[8]); sub_cell[25].push_back(temp[11]); sub_cell[25].push_back(temp[50]); sub_cell[25].push_back(temp[49]);
		sub_cell[25].push_back(temp[59]); sub_cell[25].push_back(temp[60]); sub_cell[25].push_back(temp[63]); sub_cell[25].push_back(temp[62]);
		//第27个单元
		sub_cell[26].push_back(temp[11]); sub_cell[26].push_back(temp[47]); sub_cell[26].push_back(temp[51]); sub_cell[26].push_back(temp[50]);
		sub_cell[26].push_back(temp[60]); sub_cell[26].push_back(temp[61]); sub_cell[26].push_back(temp[3]); sub_cell[26].push_back(temp[63]);

		//插入单元
		for (int j = 0; j < 27; j++)
		{
			newmesh->add_cell(sub_cell[j]);
		}
	}
	return *newmesh;
}

HexV3fMesh& ApproximationCC::CCapproximation1()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(solid_corner_point(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 9, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 9, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 1, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 7, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 1, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 7, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 8, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 8, 2));
		//第二层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 6, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 4, 2));
		//第三层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 2, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 6, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 4, 3));
		//第4层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 11, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 11, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 5, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 5, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 10, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point(i, 10, 3));

		//构建Bezier体的数据
		std::vector<std::vector<std::vector<V3f>>> sc(4);
		std::vector<std::vector<V3f>> sf(4);
		std::vector<V3f> sl(4);
		//第一层
		sl[0] = newmesh->vertices_[temp[5]]; sl[1] = newmesh->vertices_[temp[18]];
		sl[2] = newmesh->vertices_[temp[22]]; sl[3] = newmesh->vertices_[temp[1]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[16]]; sl[1] = newmesh->vertices_[temp[19]];
		sl[2] = newmesh->vertices_[temp[23]]; sl[3] = newmesh->vertices_[temp[26]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[17]]; sl[1] = newmesh->vertices_[temp[20]];
		sl[2] = newmesh->vertices_[temp[24]]; sl[3] = newmesh->vertices_[temp[27]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[6]]; sl[1] = newmesh->vertices_[temp[21]];
		sl[2] = newmesh->vertices_[temp[25]]; sl[3] = newmesh->vertices_[temp[2]];
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = newmesh->vertices_[temp[28]]; sl[1] = newmesh->vertices_[temp[32]];
		sl[2] = newmesh->vertices_[temp[34]]; sl[3] = newmesh->vertices_[temp[36]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[29]]; sl[1] = newmesh->vertices_[temp[13]];
		sl[2] = newmesh->vertices_[temp[9]]; sl[3] = newmesh->vertices_[temp[37]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[30]]; sl[1] = newmesh->vertices_[temp[14]];
		sl[2] = newmesh->vertices_[temp[10]]; sl[3] = newmesh->vertices_[temp[38]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[31]]; sl[1] = newmesh->vertices_[temp[33]];
		sl[2] = newmesh->vertices_[temp[35]]; sl[3] = newmesh->vertices_[temp[39]];
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = newmesh->vertices_[temp[40]]; sl[1] = newmesh->vertices_[temp[44]];
		sl[2] = newmesh->vertices_[temp[46]]; sl[3] = newmesh->vertices_[temp[48]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[41]]; sl[1] = newmesh->vertices_[temp[12]];
		sl[2] = newmesh->vertices_[temp[8]]; sl[3] = newmesh->vertices_[temp[49]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[42]]; sl[1] = newmesh->vertices_[temp[15]];
		sl[2] = newmesh->vertices_[temp[11]]; sl[3] = newmesh->vertices_[temp[50]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[43]]; sl[1] = newmesh->vertices_[temp[45]];
		sl[2] = newmesh->vertices_[temp[47]]; sl[3] = newmesh->vertices_[temp[51]];
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = newmesh->vertices_[temp[4]]; sl[1] = newmesh->vertices_[temp[54]];
		sl[2] = newmesh->vertices_[temp[58]]; sl[3] = newmesh->vertices_[temp[0]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[52]]; sl[1] = newmesh->vertices_[temp[55]];
		sl[2] = newmesh->vertices_[temp[59]]; sl[3] = newmesh->vertices_[temp[62]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[53]]; sl[1] = newmesh->vertices_[temp[56]];
		sl[2] = newmesh->vertices_[temp[60]]; sl[3] = newmesh->vertices_[temp[63]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[7]]; sl[1] = newmesh->vertices_[temp[57]];
		sl[2] = newmesh->vertices_[temp[61]]; sl[3] = newmesh->vertices_[temp[3]];
		sf[3] = sl;
		sc[3] = sf;

		SolidData.push_back(sc);
	}
	/*FileHandler filenew = FileHandler();
	std::string name("E://bezier.axl");
	filenew.write2AXLFile(name, SolidData);*/
	return *newmesh;
}

HexV3fMesh& ApproximationCC::CCapproximation_1()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(solid_corner_point_1(mesh->cells_[i].vertex_[vn], CellHandle(i)));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 9, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 9, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 1, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 7, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 1, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 2, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 7, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 8, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 8, 2));
		//第二层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 6, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 4, 2));
		//第三层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 2, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 4, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 6, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 1, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 5, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 4, 3));
		//第4层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 11, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 11, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 5, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point(i, 3, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 5, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 10, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point_1(i, 10, 3));

		//构建Bezier体的数据
		std::vector<std::vector<std::vector<V3f>>> sc(4);
		std::vector<std::vector<V3f>> sf(4);
		std::vector<V3f> sl(4);
		//第一层
		sl[0] = newmesh->vertices_[temp[5]]; sl[1] = newmesh->vertices_[temp[18]];
		sl[2] = newmesh->vertices_[temp[22]]; sl[3] = newmesh->vertices_[temp[1]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[16]]; sl[1] = newmesh->vertices_[temp[19]];
		sl[2] = newmesh->vertices_[temp[23]]; sl[3] = newmesh->vertices_[temp[26]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[17]]; sl[1] = newmesh->vertices_[temp[20]];
		sl[2] = newmesh->vertices_[temp[24]]; sl[3] = newmesh->vertices_[temp[27]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[6]]; sl[1] = newmesh->vertices_[temp[21]];
		sl[2] = newmesh->vertices_[temp[25]]; sl[3] = newmesh->vertices_[temp[2]];
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = newmesh->vertices_[temp[28]]; sl[1] = newmesh->vertices_[temp[32]];
		sl[2] = newmesh->vertices_[temp[34]]; sl[3] = newmesh->vertices_[temp[36]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[29]]; sl[1] = newmesh->vertices_[temp[13]];
		sl[2] = newmesh->vertices_[temp[9]]; sl[3] = newmesh->vertices_[temp[37]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[30]]; sl[1] = newmesh->vertices_[temp[14]];
		sl[2] = newmesh->vertices_[temp[10]]; sl[3] = newmesh->vertices_[temp[38]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[31]]; sl[1] = newmesh->vertices_[temp[33]];
		sl[2] = newmesh->vertices_[temp[35]]; sl[3] = newmesh->vertices_[temp[39]];
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = newmesh->vertices_[temp[40]]; sl[1] = newmesh->vertices_[temp[44]];
		sl[2] = newmesh->vertices_[temp[46]]; sl[3] = newmesh->vertices_[temp[48]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[41]]; sl[1] = newmesh->vertices_[temp[12]];
		sl[2] = newmesh->vertices_[temp[8]]; sl[3] = newmesh->vertices_[temp[49]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[42]]; sl[1] = newmesh->vertices_[temp[15]];
		sl[2] = newmesh->vertices_[temp[11]]; sl[3] = newmesh->vertices_[temp[50]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[43]]; sl[1] = newmesh->vertices_[temp[45]];
		sl[2] = newmesh->vertices_[temp[47]]; sl[3] = newmesh->vertices_[temp[51]];
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = newmesh->vertices_[temp[4]]; sl[1] = newmesh->vertices_[temp[54]];
		sl[2] = newmesh->vertices_[temp[58]]; sl[3] = newmesh->vertices_[temp[0]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[52]]; sl[1] = newmesh->vertices_[temp[55]];
		sl[2] = newmesh->vertices_[temp[59]]; sl[3] = newmesh->vertices_[temp[62]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[53]]; sl[1] = newmesh->vertices_[temp[56]];
		sl[2] = newmesh->vertices_[temp[60]]; sl[3] = newmesh->vertices_[temp[63]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[7]]; sl[1] = newmesh->vertices_[temp[57]];
		sl[2] = newmesh->vertices_[temp[61]]; sl[3] = newmesh->vertices_[temp[3]];
		sf[3] = sl;
		sc[3] = sf;

		SolidData.push_back(sc);
	}
	/*FileHandler filenew = FileHandler();
	std::string name("E://bezier.axl");
	filenew.write2AXLFile(name, SolidData);*/
	return *newmesh;
}

HexV3fMesh ApproximationCC::constructBezierSolid()
{
	HexV3fMesh beziermesh;
	int socB = 0;
	std::cout << "Start construct BezierSolid" << std::endl;
	for (std::vector<std::vector<std::vector<std::vector<V3f>>>>::iterator SolidData_it = SolidData.begin();
		SolidData_it != SolidData.end(); SolidData_it++)
	{
		if (socB % 100 == 0)
			std::cout << "Schedule of construct BezierSolid " << socB + 1 << ":" << SolidData.size() << std::endl;
		socB++;
		std::vector<std::vector<std::vector<V3f>>> bsdata = *SolidData_it;
		BezierSolid BS(bsdata);
		int size = 2;//要生成size*size*size个体
		// std::cout << "Start add control point" << std::endl;
		std::vector<std::vector<V3f>> bezierv(size + 1);
		for (int i = 0; i < bezierv.size(); i++)
		{
			for (int j = 0; j < (size + 1) * (size + 1); j++)
			{
				// std::cout << "Schedule of add control point " << i + 1 << ":"
					// << bezierv.size() << " " << j + 1 << ":" << (size + 1)*(size + 1) << std::endl;
				bezierv[i].push_back(BS.P(j % (size + 1) * (1.0 / size), j / (size + 1) * (1.0 / size), i % (size + 1) * (1.0 / size)));
			}
		}
		// std::cout << "Add control point completely" << std::endl;
		// std::cout << "Start add mesh" << std::endl;
		for (int j = 0; j < size; j++)
		{
			for (int i = 0; i < size * size; i++)
			{
				// std::cout << "Schedule of add mesh " << j + 1 << ":"
					// << size << " " << i + 1 << ":" << (size + 1)*(size + 1) << std::endl;
				std::vector<VertexHandle> ctemp;
				int t = i + i / size;
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 1]));
				beziermesh.add_cell(ctemp);
			}
		}
		// std::cout << "Add mesh completely" << std::endl;
	}
	std::cout << "Construct BezierSolid completely" << std::endl;
	/*int bdyfn = 0;
	for (int i = 0;i < beziermesh.faces_.size();i++)
	{
		if (beziermesh.is_bdy(FaceHandle(i)))
			bdyfn++;
	}*/
	return beziermesh;
}

HexV3fMesh ApproximationCC::constructBezierSolid1()
{
	HexV3fMesh beziermesh;
	std::cout << "Start construct BezierSolid" << std::endl;
	for (std::vector<std::vector<std::vector<std::vector<V3f>>>>::iterator SolidData_it = SolidData.begin();
		SolidData_it != SolidData.end(); SolidData_it++)
	{
		int socB = 0;
		std::cout << "Schedule of construct BezierSolid " << ++socB << ":" << SolidData.size() << std::endl;
		std::vector<std::vector<std::vector<V3f>>> bsdata = *SolidData_it;
		BezierSolid BS(bsdata);
		int size = 30;//要生成size*size*size个体
		std::cout << "Start add control point" << std::endl;
		std::vector<std::vector<V3f>> bezierv(size + 1);
		for (int i = 0; i < bezierv.size(); i++)
		{
			for (int j = 0; j < (size + 1) * (size + 1); j++)
			{
				std::cout << "Schedule of add control point " << i + 1 << ":"
					<< bezierv.size() << " " << j + 1 << ":" << (size + 1) * (size + 1) << std::endl;
				bezierv[i].push_back(BS.P1(j % (size + 1) * (1.0 / size), j / (size + 1) * (1.0 / size), i % (size + 1) * (1.0 / size)));
			}
		}
		std::cout << "Add control point completely" << std::endl;
		std::cout << "Start add mesh" << std::endl;
		for (int j = 0; j < size; j++)
		{
			for (int i = 0; i < size * size; i++)
			{
				std::cout << "Schedule of add mesh " << j + 1 << ":"
					<< size << " " << i + 1 << ":" << size * size << std::endl;
				std::vector<VertexHandle> ctemp;
				int t = i + i / size;
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 1]));
				beziermesh.add_cell(ctemp);
			}
		}
		std::cout << "Add mesh completely" << std::endl;
	}
	std::cout << "Construct BezierSolid completely" << std::endl;
	/*int bdyfn = 0;
	for (int i = 0;i < beziermesh.faces_.size();i++)
	{
		if (beziermesh.is_bdy(FaceHandle(i)))
			bdyfn++;
	}*/
	return beziermesh;
}

std::vector<HexV3fMesh>& ApproximationCC::CCapproximationDAT()
{
	uint64_t cell_num = mesh->cells_.size();
	vmesh.resize(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = vmesh[i].add_vertex(solid_corner_point(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = vmesh[i].add_vertex(cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 9, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 9, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 1, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 2, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 2, 6));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 7, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 1, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 2, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 2, 2));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 7, 2));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 8, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 8, 2));
		//第二层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 2, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 4, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 4, 6));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 6, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 0, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 1, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 0, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 1, 2));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 0, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 5, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 5, 2));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 4, 2));
		//第三层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 2, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 4, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 4, 7));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 6, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 0, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 1, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 0, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 1, 3));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 0, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 5, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 5, 3));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 4, 3));
		//第4层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 11, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 11, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 3, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 3, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 3, 7));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 5, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 3, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 3, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point(i, 3, 3));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 5, 3));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 10, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point(i, 10, 3));
		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(27);
		/*第一层*/
		//第1个单元
		sub_cell[0].push_back(temp[5]); sub_cell[0].push_back(temp[16]); sub_cell[0].push_back(temp[19]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[28]); sub_cell[0].push_back(temp[29]); sub_cell[0].push_back(temp[13]); sub_cell[0].push_back(temp[32]);
		//第2个单元
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[17]); sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[19]);
		sub_cell[1].push_back(temp[29]); sub_cell[1].push_back(temp[30]); sub_cell[1].push_back(temp[14]); sub_cell[1].push_back(temp[13]);
		//第3个单元
		sub_cell[2].push_back(temp[17]); sub_cell[2].push_back(temp[6]); sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[20]);
		sub_cell[2].push_back(temp[30]); sub_cell[2].push_back(temp[31]); sub_cell[2].push_back(temp[33]); sub_cell[2].push_back(temp[14]);
		//第4个单元
		sub_cell[3].push_back(temp[18]); sub_cell[3].push_back(temp[19]); sub_cell[3].push_back(temp[23]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[32]); sub_cell[3].push_back(temp[13]); sub_cell[3].push_back(temp[9]); sub_cell[3].push_back(temp[34]);
		//第5个单元
		sub_cell[4].push_back(temp[19]); sub_cell[4].push_back(temp[20]); sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[13]); sub_cell[4].push_back(temp[14]); sub_cell[4].push_back(temp[10]); sub_cell[4].push_back(temp[9]);
		//第6个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[21]); sub_cell[5].push_back(temp[25]); sub_cell[5].push_back(temp[24]);
		sub_cell[5].push_back(temp[14]); sub_cell[5].push_back(temp[33]); sub_cell[5].push_back(temp[35]); sub_cell[5].push_back(temp[10]);
		//第7个单元
		sub_cell[6].push_back(temp[22]); sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]); sub_cell[6].push_back(temp[1]);
		sub_cell[6].push_back(temp[34]); sub_cell[6].push_back(temp[9]); sub_cell[6].push_back(temp[37]); sub_cell[6].push_back(temp[36]);
		//第8个单元
		sub_cell[7].push_back(temp[23]); sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[27]); sub_cell[7].push_back(temp[26]);
		sub_cell[7].push_back(temp[9]); sub_cell[7].push_back(temp[10]); sub_cell[7].push_back(temp[38]); sub_cell[7].push_back(temp[37]);
		//第9个单元
		sub_cell[8].push_back(temp[24]); sub_cell[8].push_back(temp[25]); sub_cell[8].push_back(temp[2]); sub_cell[8].push_back(temp[27]);
		sub_cell[8].push_back(temp[10]); sub_cell[8].push_back(temp[35]); sub_cell[8].push_back(temp[39]); sub_cell[8].push_back(temp[38]);
		/*第二层*/
		//第10个单元
		sub_cell[9].push_back(temp[28]); sub_cell[9].push_back(temp[29]); sub_cell[9].push_back(temp[13]); sub_cell[9].push_back(temp[32]);
		sub_cell[9].push_back(temp[40]); sub_cell[9].push_back(temp[41]); sub_cell[9].push_back(temp[12]); sub_cell[9].push_back(temp[44]);
		//第11个单元
		sub_cell[10].push_back(temp[29]); sub_cell[10].push_back(temp[30]); sub_cell[10].push_back(temp[14]); sub_cell[10].push_back(temp[13]);
		sub_cell[10].push_back(temp[41]); sub_cell[10].push_back(temp[42]); sub_cell[10].push_back(temp[15]); sub_cell[10].push_back(temp[12]);
		//第12个单元
		sub_cell[11].push_back(temp[30]); sub_cell[11].push_back(temp[31]); sub_cell[11].push_back(temp[33]); sub_cell[11].push_back(temp[14]);
		sub_cell[11].push_back(temp[42]); sub_cell[11].push_back(temp[43]); sub_cell[11].push_back(temp[45]); sub_cell[11].push_back(temp[15]);
		//第13个单元
		sub_cell[12].push_back(temp[32]); sub_cell[12].push_back(temp[13]); sub_cell[12].push_back(temp[9]); sub_cell[12].push_back(temp[34]);
		sub_cell[12].push_back(temp[44]); sub_cell[12].push_back(temp[12]); sub_cell[12].push_back(temp[8]); sub_cell[12].push_back(temp[46]);
		//第14个单元
		sub_cell[13].push_back(temp[13]); sub_cell[13].push_back(temp[14]); sub_cell[13].push_back(temp[10]); sub_cell[13].push_back(temp[9]);
		sub_cell[13].push_back(temp[12]); sub_cell[13].push_back(temp[15]); sub_cell[13].push_back(temp[11]); sub_cell[13].push_back(temp[8]);
		//第15个单元
		sub_cell[14].push_back(temp[14]); sub_cell[14].push_back(temp[33]); sub_cell[14].push_back(temp[35]); sub_cell[14].push_back(temp[10]);
		sub_cell[14].push_back(temp[15]); sub_cell[14].push_back(temp[45]); sub_cell[14].push_back(temp[47]); sub_cell[14].push_back(temp[11]);
		//第16个单元
		sub_cell[15].push_back(temp[34]); sub_cell[15].push_back(temp[9]); sub_cell[15].push_back(temp[37]); sub_cell[15].push_back(temp[36]);
		sub_cell[15].push_back(temp[46]); sub_cell[15].push_back(temp[8]); sub_cell[15].push_back(temp[49]); sub_cell[15].push_back(temp[48]);
		//第17个单元
		sub_cell[16].push_back(temp[9]); sub_cell[16].push_back(temp[10]); sub_cell[16].push_back(temp[38]); sub_cell[16].push_back(temp[37]);
		sub_cell[16].push_back(temp[8]); sub_cell[16].push_back(temp[11]); sub_cell[16].push_back(temp[50]); sub_cell[16].push_back(temp[49]);
		//第18个单元
		sub_cell[17].push_back(temp[10]); sub_cell[17].push_back(temp[35]); sub_cell[17].push_back(temp[39]); sub_cell[17].push_back(temp[38]);
		sub_cell[17].push_back(temp[11]); sub_cell[17].push_back(temp[47]); sub_cell[17].push_back(temp[51]); sub_cell[17].push_back(temp[50]);
		/*第三层*/
		//第19个单元
		sub_cell[18].push_back(temp[40]); sub_cell[18].push_back(temp[41]); sub_cell[18].push_back(temp[12]); sub_cell[18].push_back(temp[44]);
		sub_cell[18].push_back(temp[4]); sub_cell[18].push_back(temp[52]); sub_cell[18].push_back(temp[55]); sub_cell[18].push_back(temp[54]);
		//第20个单元
		sub_cell[19].push_back(temp[41]); sub_cell[19].push_back(temp[42]); sub_cell[19].push_back(temp[15]); sub_cell[19].push_back(temp[12]);
		sub_cell[19].push_back(temp[52]); sub_cell[19].push_back(temp[53]); sub_cell[19].push_back(temp[56]); sub_cell[19].push_back(temp[55]);
		//第21个单元
		sub_cell[20].push_back(temp[42]); sub_cell[20].push_back(temp[43]); sub_cell[20].push_back(temp[45]); sub_cell[20].push_back(temp[15]);
		sub_cell[20].push_back(temp[53]); sub_cell[20].push_back(temp[7]); sub_cell[20].push_back(temp[57]); sub_cell[20].push_back(temp[56]);
		//第22个单元
		sub_cell[21].push_back(temp[44]); sub_cell[21].push_back(temp[12]); sub_cell[21].push_back(temp[8]); sub_cell[21].push_back(temp[46]);
		sub_cell[21].push_back(temp[54]); sub_cell[21].push_back(temp[55]); sub_cell[21].push_back(temp[59]); sub_cell[21].push_back(temp[58]);
		//第23个单元
		sub_cell[22].push_back(temp[12]); sub_cell[22].push_back(temp[15]); sub_cell[22].push_back(temp[11]); sub_cell[22].push_back(temp[8]);
		sub_cell[22].push_back(temp[55]); sub_cell[22].push_back(temp[56]); sub_cell[22].push_back(temp[60]); sub_cell[22].push_back(temp[59]);
		//第24个单元
		sub_cell[23].push_back(temp[15]); sub_cell[23].push_back(temp[45]); sub_cell[23].push_back(temp[47]); sub_cell[23].push_back(temp[11]);
		sub_cell[23].push_back(temp[56]); sub_cell[23].push_back(temp[57]); sub_cell[23].push_back(temp[61]); sub_cell[23].push_back(temp[60]);
		//第25个单元
		sub_cell[24].push_back(temp[46]); sub_cell[24].push_back(temp[8]); sub_cell[24].push_back(temp[49]); sub_cell[24].push_back(temp[48]);
		sub_cell[24].push_back(temp[58]); sub_cell[24].push_back(temp[59]); sub_cell[24].push_back(temp[62]); sub_cell[24].push_back(temp[0]);
		//第26个单元
		sub_cell[25].push_back(temp[8]); sub_cell[25].push_back(temp[11]); sub_cell[25].push_back(temp[50]); sub_cell[25].push_back(temp[49]);
		sub_cell[25].push_back(temp[59]); sub_cell[25].push_back(temp[60]); sub_cell[25].push_back(temp[63]); sub_cell[25].push_back(temp[62]);
		//第27个单元
		sub_cell[26].push_back(temp[11]); sub_cell[26].push_back(temp[47]); sub_cell[26].push_back(temp[51]); sub_cell[26].push_back(temp[50]);
		sub_cell[26].push_back(temp[60]); sub_cell[26].push_back(temp[61]); sub_cell[26].push_back(temp[3]); sub_cell[26].push_back(temp[63]);

		//插入单元
		for (int j = 0; j < 27; j++)
		{
			vmesh[i].add_cell(sub_cell[j]);
		}
	}
	return vmesh;
}

std::vector<HexV3fMesh>& ApproximationCC::constructBezierSolidDAT()
{
	int datsize = SolidData.size();
	int datnum = 0;
	isoline.resize(datsize);
	for (std::vector<std::vector<std::vector<std::vector<V3f>>>>::iterator SolidData_it = SolidData.begin();
		SolidData_it != SolidData.end(); SolidData_it++)
	{
		HexV3fMesh beziermesh;
		std::vector<std::vector<std::vector<V3f>>> bsdata = *SolidData_it;
		BezierSolid BS(bsdata);
		int size = 8;//要生成size*size*size个体
		std::vector<std::vector<V3f>> bezierv(size + 1);
		for (int i = 0; i < bezierv.size(); i++)
		{
			for (int j = 0; j < (size + 1) * (size + 1); j++)
			{
				bezierv[i].push_back(BS.P(j % (size + 1) * (1.0 / size), j / (size + 1) * (1.0 / size), i % (size + 1) * (1.0 / size)));
			}
		}
		//添加等参线
		/*std::vector<V3f> face1;
		for (int i = 0;i < 25;i++)
		{
			face1.push_back(BS.P(i % 5 * (1.0 / 4), i / 5 * (1.0 / 4), 0));
		}
		std::vector<V3f> face2;
		for (int i = 0;i < 25;i++)
		{
			face2.push_back(BS.P(i % 5 * (1.0 / 4), i / 5 * (1.0 / 4), 1));
		}
		std::vector<V3f> face3;
		for (int i = 0;i < 25;i++)
		{
			face3.push_back(BS.P(0, i / 5 * (1.0 / 4), i % 5 * (1.0 / 4)));
		}
		std::vector<V3f> face4;
		for (int i = 0;i < 25;i++)
		{
			face4.push_back(BS.P(1, i / 5 * (1.0 / 4), i % 5 * (1.0 / 4)));
		}
		std::vector<V3f> face5;
		for (int i = 0;i < 25;i++)
		{
			face5.push_back(BS.P(i % 5 * (1.0 / 4), 0, i / 5 * (1.0 / 4)));
		}
		std::vector<V3f> face6;
		for (int i = 0;i < 25;i++)
		{
			face6.push_back(BS.P(i % 5 * (1.0 / 4), 1, i / 5 * (1.0 / 4)));
		}
		isoline[datnum].push_back(face1);isoline[datnum].push_back(face2);isoline[datnum].push_back(face3);
		isoline[datnum].push_back(face4);isoline[datnum].push_back(face5);isoline[datnum].push_back(face6);*/

		std::vector<V3f> line1;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line1.push_back(BS.P(0.25, j * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line1.push_back(BS.P(0.25, 1, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line1.push_back(BS.P(0.25, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line1.push_back(BS.P(0.25, 0, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line2;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line2.push_back(BS.P(0.5, j * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line2.push_back(BS.P(0.5, 1, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line2.push_back(BS.P(0.5, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line2.push_back(BS.P(0.5, 0, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line3;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line3.push_back(BS.P(0.75, j * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line3.push_back(BS.P(0.75, 1, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line3.push_back(BS.P(0.75, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line3.push_back(BS.P(0.75, 0, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line11;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line11.push_back(BS.P(0.0, j * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line11.push_back(BS.P(0.0, 1, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line11.push_back(BS.P(0.0, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line11.push_back(BS.P(0.0, 0, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line12;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line12.push_back(BS.P(1.0, j * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line12.push_back(BS.P(1.0, 1, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line12.push_back(BS.P(1.0, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line12.push_back(BS.P(1.0, 0, (size - j) * (1.0 / size)));
		}


		std::vector<V3f> line4;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line4.push_back(BS.P(j * (1.0 / size), 0.25, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line4.push_back(BS.P(1, 0.25, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line4.push_back(BS.P((size - j) * (1.0 / size), 0.25, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line4.push_back(BS.P(0, 0.25, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line5;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line5.push_back(BS.P(j * (1.0 / size), 0.5, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line5.push_back(BS.P(1, 0.5, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line5.push_back(BS.P((size - j) * (1.0 / size), 0.5, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line5.push_back(BS.P(0, 0.5, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line6;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line6.push_back(BS.P(j * (1.0 / size), 0.75, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line6.push_back(BS.P(1, 0.75, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line6.push_back(BS.P((size - j) * (1.0 / size), 0.75, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line6.push_back(BS.P(0, 0.75, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line21;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line21.push_back(BS.P(j * (1.0 / size), 0.0, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line21.push_back(BS.P(1, 0.0, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line21.push_back(BS.P((size - j) * (1.0 / size), 0.0, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line21.push_back(BS.P(0, 0.0, (size - j) * (1.0 / size)));
		}
		std::vector<V3f> line22;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line22.push_back(BS.P(j * (1.0 / size), 1.0, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line22.push_back(BS.P(1, 1.0, j * (1.0 / size)));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line22.push_back(BS.P((size - j) * (1.0 / size), 1.0, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line22.push_back(BS.P(0, 1.0, (size - j) * (1.0 / size)));
		}


		std::vector<V3f> line7;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line7.push_back(BS.P(j * (1.0 / size), 1, 0.5));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line7.push_back(BS.P(1, (size - j) * (1.0 / size), 0.5));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line7.push_back(BS.P((size - j) * (1.0 / size), 0, 0.5));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line7.push_back(BS.P(0, j * (1.0 / size), 0.5));
		}
		std::vector<V3f> line8;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line8.push_back(BS.P(j * (1.0 / size), 1, 0.75));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line8.push_back(BS.P(1, (size - j) * (1.0 / size), 0.75));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line8.push_back(BS.P((size - j) * (1.0 / size), 0, 0.75));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line8.push_back(BS.P(0, j * (1.0 / size), 0.75));
		}
		std::vector<V3f> line9;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line9.push_back(BS.P(j * (1.0 / size), 1, 0.25));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line9.push_back(BS.P(1, (size - j) * (1.0 / size), 0.25));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line9.push_back(BS.P((size - j) * (1.0 / size), 0, 0.25));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line9.push_back(BS.P(0, j * (1.0 / size), 0.25));
		}
		std::vector<V3f> line31;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line31.push_back(BS.P(j * (1.0 / size), 1, 0));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line31.push_back(BS.P(1, (size - j) * (1.0 / size), 0));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line31.push_back(BS.P((size - j) * (1.0 / size), 0, 0));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line31.push_back(BS.P(0, j * (1.0 / size), 0));
		}
		std::vector<V3f> line32;
		for (int j = 0; j < size + 1; j++)//第1个面
		{
			line32.push_back(BS.P(j * (1.0 / size), 1, 1));
		}
		for (int j = 0; j < size + 1; j++)//第2个面
		{
			line32.push_back(BS.P(1, (size - j) * (1.0 / size), 1));
		}
		for (int j = 0; j < size + 1; j++)//第3个面
		{
			line32.push_back(BS.P((size - j) * (1.0 / size), 0, 1));
		}
		for (int j = 0; j < size + 1; j++)//第4个面
		{
			line32.push_back(BS.P(0, j * (1.0 / size), 1));
		}

		isoline[datnum].push_back(line1); isoline[datnum].push_back(line2); isoline[datnum].push_back(line3);
		isoline[datnum].push_back(line4); isoline[datnum].push_back(line5); isoline[datnum].push_back(line6);
		isoline[datnum].push_back(line7); isoline[datnum].push_back(line8); isoline[datnum].push_back(line9);
		isoline[datnum].push_back(line11); isoline[datnum].push_back(line12); isoline[datnum].push_back(line21); isoline[datnum].push_back(line22);
		isoline[datnum].push_back(line31); isoline[datnum].push_back(line32);

		for (int j = 0; j < size; j++)
		{
			for (int i = 0; i < size * size; i++)
			{
				std::vector<VertexHandle> ctemp;
				int t = i + i / size;
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j][t + size + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + 1]));
				ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 2])); ctemp.push_back(beziermesh.add_vertex(bezierv[j + 1][t + size + 1]));
				beziermesh.add_cell(ctemp);
			}
		}
		vBmesh.push_back(beziermesh);
		datnum++;
	}
	return vBmesh;
}

void ApproximationCC::cal_ccver()
{
	int mv_size = mesh->vertices_.size();
	for (int i = 0; i < mv_size; i++)
	{
		ccver.push_back(solid_corner_point(VertexHandle(i)));
	}
}

void ApproximationCC::cal_limitccver()
{
	int mv_size = mesh->vertices_.size();
	for (int i = 0; i < mv_size; i++)
	{
		VertexHandle vertex = VertexHandle(i);
		V3f new_vertex;
		V3f vertex_sum;
		if (mesh->is_bdy(vertex))
		{
			int en = 0;
			std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
			int v_e_n = v_e.size();
			V3f e_av;
			for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
			{
				if (mesh->is_bdy(*v_e_it))
				{
					en++;
					V3f temp_ev = cal_edge_mid(*v_e_it);
					e_av = e_av + temp_ev;
				}
			}
			//e_av = e_av / en;
			//vertex_sum = vertex_sum + mesh->vertices_[vertex] * en * en;
			std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
			V3f f_av;
			for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
			{
				if (mesh->is_bdy(*v_f_it))
				{
					V3f temp_fv = cal_face_mid(*v_f_it);
					f_av = f_av + temp_fv;
				}
			}
			//f_av = f_av / en;
			vertex_sum = vertex_sum + cal_new_vertex(vertex) * en * en + e_av * 4 + f_av;
			new_vertex = vertex_sum / ((en + 5) * en);
		}
		else
		{
			std::set<EdgeHandle, compare_OVM> v_e = mesh->neighbor_v[vertex];
			int en = v_e.size();
			int sum = 0;
			sum += 16 * (en - 2);
			for (std::set<EdgeHandle, compare_OVM>::iterator v_e_it = v_e.begin(); v_e_it != v_e.end(); v_e_it++)
			{
				V3f temp_ev = cal_edge_mid(*v_e_it);
				int m = mesh->neighbor_e[*v_e_it].size();
				sum += 4 * m;
				vertex_sum = vertex_sum + temp_ev * 4 * m;
			}
			std::set<FaceHandle, compare_OVM> v_f = mesh->neighbor_v_f[vertex];
			for (std::set<FaceHandle, compare_OVM>::iterator v_f_it = v_f.begin(); v_f_it != v_f.end(); v_f_it++)
			{
				V3f temp_fv = cal_face_mid(*v_f_it);
				vertex_sum = vertex_sum + temp_fv * 4;
				sum += 4;
			}
			std::set<CellHandle, compare_OVM> v_c = mesh->neighbor_v_c[vertex];
			for (std::set<CellHandle, compare_OVM>::iterator v_c_it = v_c.begin(); v_c_it != v_c.end(); v_c_it++)
			{
				V3f temp_cv = cal_cell_mid(*v_c_it);
				vertex_sum = vertex_sum + temp_cv;
				sum += 1;
			}
			vertex_sum = vertex_sum + cal_new_vertex(vertex) * ((en - 2) * 4 * 4);
			new_vertex = vertex_sum / sum;
		}
		ccver.push_back(new_vertex);
	}
}

V3f ApproximationCC::selectapoint()
{
	int i = 0;
	for (std::vector<V3f>::iterator mv = mesh->vertices_.begin(); mv != mesh->vertices_.end(); mv++)
	{
		if (!mesh->is_bdy(VertexHandle(i)))
		{
			if (mesh->neighbor_v[i].size() == 5)
			{
				V3f v = solid_corner_point(VertexHandle(i));
				return v;
			}
			else
				i++;
		}
		else
			i++;

	}
}

void ApproximationCC::construct_bss()
{
	uint64_t cell_num = mesh->cells_.size();
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.CP_Local_to_Global[i].resize(64);
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_corner_point(mesh->cells_[i].vertex_[vn]));
			//构造局部的小单元的每个点对应于全局模型中的点
			bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
				bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_inter_point(CellHandle(i), mesh->cells_[i].vertex_[vn]));
			bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
				bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		}
		//第一层
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 9, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 9, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 1, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 2, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 2, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 7, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 1, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 2, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 2, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 7, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 8, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 8, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		//第二层
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 2, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 4, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 4, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 6, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 0, 5));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 1, 6));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 0, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 1, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 0, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 5, 1));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 5, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 4, 2));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		//第三层
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 2, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 4, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 4, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 6, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 0, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 1, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 0, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 1, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 0, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 5, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 5, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 4, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		//第4层
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 11, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 11, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 3, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 3, 4));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 3, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 5, 7));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 3, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 3, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_face_point(i, 3, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 5, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 10, 0));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(solid_edge_point(i, 10, 3));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<VertexHandle>>> sc(4);
		std::vector<std::vector<VertexHandle>> sf(4);
		std::vector<VertexHandle> sl(4);
		//第一层
		sl[0] = temp[5]; sl[1] = temp[18]; sl[2] = temp[22]; sl[3] = temp[1];
		sf[0] = sl;
		sl[0] = temp[16]; sl[1] = temp[19]; sl[2] = temp[23]; sl[3] = temp[26];
		sf[1] = sl;
		sl[0] = temp[17]; sl[1] = temp[20]; sl[2] = temp[24]; sl[3] = temp[27];
		sf[2] = sl;
		sl[0] = temp[6]; sl[1] = temp[21]; sl[2] = temp[25]; sl[3] = temp[2];
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = temp[28]; sl[1] = temp[32]; sl[2] = temp[34]; sl[3] = temp[36];
		sf[0] = sl;
		sl[0] = temp[29]; sl[1] = temp[13]; sl[2] = temp[9]; sl[3] = temp[37];
		sf[1] = sl;
		sl[0] = temp[30]; sl[1] = temp[14]; sl[2] = temp[10]; sl[3] = temp[38];
		sf[2] = sl;
		sl[0] = temp[31]; sl[1] = temp[33]; sl[2] = temp[35]; sl[3] = temp[39];
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = temp[40]; sl[1] = temp[44]; sl[2] = temp[46]; sl[3] = temp[48];
		sf[0] = sl;
		sl[0] = temp[41]; sl[1] = temp[12]; sl[2] = temp[8]; sl[3] = temp[49];
		sf[1] = sl;
		sl[0] = temp[42]; sl[1] = temp[15]; sl[2] = temp[11]; sl[3] = temp[50];
		sf[2] = sl;
		sl[0] = temp[43]; sl[1] = temp[45]; sl[2] = temp[47]; sl[3] = temp[51];
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = temp[4]; sl[1] = temp[54]; sl[2] = temp[58]; sl[3] = temp[0];
		sf[0] = sl;
		sl[0] = temp[52]; sl[1] = temp[55]; sl[2] = temp[59]; sl[3] = temp[62];
		sf[1] = sl;
		sl[0] = temp[53]; sl[1] = temp[56]; sl[2] = temp[60]; sl[3] = temp[63];
		sf[2] = sl;
		sl[0] = temp[7]; sl[1] = temp[57]; sl[2] = temp[61]; sl[3] = temp[3];
		sf[3] = sl;
		sc[3] = sf;
		bss.Control_Point_id[i] = sc;

		//构建Bezier体的数据
		std::vector<std::vector<std::vector<V3f>>> sc1(4);
		std::vector<std::vector<V3f>> sf1(4);
		std::vector<V3f> sl1(4);
		//第一层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[5]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[18]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[22]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[1]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[16]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[19]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[23]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[26]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[17]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[20]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[24]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[27]];
		sf1[2] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[6]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[21]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[25]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[2]];
		sf1[3] = sl1;
		sc1[0] = sf1;
		//第二层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[28]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[32]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[34]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[36]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[29]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[13]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[9]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[37]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[30]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[14]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[10]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[38]];
		sf1[2] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[31]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[33]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[35]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[39]];
		sf1[3] = sl1;
		sc1[1] = sf1;
		//第三层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[40]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[44]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[46]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[48]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[41]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[12]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[8]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[49]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[42]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[15]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[11]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[50]];
		sf1[2] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[43]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[45]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[47]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[51]];
		sf1[3] = sl1;
		sc1[2] = sf1;
		//第四层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[4]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[54]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[58]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[0]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[52]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[55]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[59]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[62]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[53]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[56]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[60]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[63]];
		sf1[2] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[7]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[57]];
		sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[61]]; sl1[3] = bss.Local_BSplineSolid[i].vertices_[temp[3]];
		sf1[3] = sl1;
		sc1[3] = sf1;
		SolidData.push_back(sc1);

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

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(27);
		std::vector<std::vector<VertexHandle>> sub_cell_global(27);
		/*第一层*/
		//第1个单元
		sub_cell[0].push_back(temp[5]); sub_cell[0].push_back(temp[16]); sub_cell[0].push_back(temp[19]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[28]); sub_cell[0].push_back(temp[29]); sub_cell[0].push_back(temp[13]); sub_cell[0].push_back(temp[32]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[5]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[16]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[19]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[18]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[28]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[29]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[13]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[32]]);
		//第2个单元
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[17]); sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[19]);
		sub_cell[1].push_back(temp[29]); sub_cell[1].push_back(temp[30]); sub_cell[1].push_back(temp[14]); sub_cell[1].push_back(temp[13]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[16]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[17]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[20]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[19]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[29]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[30]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[14]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		//第3个单元
		sub_cell[2].push_back(temp[17]); sub_cell[2].push_back(temp[6]); sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[20]);
		sub_cell[2].push_back(temp[30]); sub_cell[2].push_back(temp[31]); sub_cell[2].push_back(temp[33]); sub_cell[2].push_back(temp[14]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[17]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[6]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[21]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[20]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[30]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[31]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[33]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[14]]);
		//第4个单元
		sub_cell[3].push_back(temp[18]); sub_cell[3].push_back(temp[19]); sub_cell[3].push_back(temp[23]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[32]); sub_cell[3].push_back(temp[13]); sub_cell[3].push_back(temp[9]); sub_cell[3].push_back(temp[34]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[18]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[19]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[23]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[22]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[32]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[9]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[34]]);
		//第5个单元
		sub_cell[4].push_back(temp[19]); sub_cell[4].push_back(temp[20]); sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[13]); sub_cell[4].push_back(temp[14]); sub_cell[4].push_back(temp[10]); sub_cell[4].push_back(temp[9]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[19]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[20]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[24]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[23]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[13]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[14]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[10]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		//第6个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[21]); sub_cell[5].push_back(temp[25]); sub_cell[5].push_back(temp[24]);
		sub_cell[5].push_back(temp[14]); sub_cell[5].push_back(temp[33]); sub_cell[5].push_back(temp[35]); sub_cell[5].push_back(temp[10]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[20]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[21]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[25]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[24]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[14]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[33]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[35]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[10]]);
		//第7个单元
		sub_cell[6].push_back(temp[22]); sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]); sub_cell[6].push_back(temp[1]);
		sub_cell[6].push_back(temp[34]); sub_cell[6].push_back(temp[9]); sub_cell[6].push_back(temp[37]); sub_cell[6].push_back(temp[36]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[22]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[23]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[26]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[1]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[34]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[37]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[36]]);
		//第8个单元
		sub_cell[7].push_back(temp[23]); sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[27]); sub_cell[7].push_back(temp[26]);
		sub_cell[7].push_back(temp[9]); sub_cell[7].push_back(temp[10]); sub_cell[7].push_back(temp[38]); sub_cell[7].push_back(temp[37]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[23]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[24]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[27]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[26]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[9]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[10]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[38]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[37]]);
		//第9个单元
		sub_cell[8].push_back(temp[24]); sub_cell[8].push_back(temp[25]); sub_cell[8].push_back(temp[2]); sub_cell[8].push_back(temp[27]);
		sub_cell[8].push_back(temp[10]); sub_cell[8].push_back(temp[35]); sub_cell[8].push_back(temp[39]); sub_cell[8].push_back(temp[38]);
		sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[24]]); sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[25]]);
		sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[2]]); sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[27]]);
		sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[10]]); sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[35]]);
		sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[39]]); sub_cell_global[8].push_back(bss.CP_Local_to_Global[i][temp[38]]);
		/*第二层*/
		//第10个单元
		sub_cell[9].push_back(temp[28]); sub_cell[9].push_back(temp[29]); sub_cell[9].push_back(temp[13]); sub_cell[9].push_back(temp[32]);
		sub_cell[9].push_back(temp[40]); sub_cell[9].push_back(temp[41]); sub_cell[9].push_back(temp[12]); sub_cell[9].push_back(temp[44]);
		sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[28]]); sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[29]]);
		sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[13]]); sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[32]]);
		sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[40]]); sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[41]]);
		sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[12]]); sub_cell_global[9].push_back(bss.CP_Local_to_Global[i][temp[44]]);
		//第11个单元
		sub_cell[10].push_back(temp[29]); sub_cell[10].push_back(temp[30]); sub_cell[10].push_back(temp[14]); sub_cell[10].push_back(temp[13]);
		sub_cell[10].push_back(temp[41]); sub_cell[10].push_back(temp[42]); sub_cell[10].push_back(temp[15]); sub_cell[10].push_back(temp[12]);
		sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[29]]); sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[30]]);
		sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[14]]); sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[41]]); sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[42]]);
		sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[10].push_back(bss.CP_Local_to_Global[i][temp[12]]);
		//第12个单元
		sub_cell[11].push_back(temp[30]); sub_cell[11].push_back(temp[31]); sub_cell[11].push_back(temp[33]); sub_cell[11].push_back(temp[14]);
		sub_cell[11].push_back(temp[42]); sub_cell[11].push_back(temp[43]); sub_cell[11].push_back(temp[45]); sub_cell[11].push_back(temp[15]);
		sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[30]]); sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[31]]);
		sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[33]]); sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[14]]);
		sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[42]]); sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[43]]);
		sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[45]]); sub_cell_global[11].push_back(bss.CP_Local_to_Global[i][temp[15]]);
		//第13个单元
		sub_cell[12].push_back(temp[32]); sub_cell[12].push_back(temp[13]); sub_cell[12].push_back(temp[9]); sub_cell[12].push_back(temp[34]);
		sub_cell[12].push_back(temp[44]); sub_cell[12].push_back(temp[12]); sub_cell[12].push_back(temp[8]); sub_cell[12].push_back(temp[46]);
		sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[32]]); sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[9]]); sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[34]]);
		sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[44]]); sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[12]]);
		sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[8]]); sub_cell_global[12].push_back(bss.CP_Local_to_Global[i][temp[46]]);
		//第14个单元
		sub_cell[13].push_back(temp[13]); sub_cell[13].push_back(temp[14]); sub_cell[13].push_back(temp[10]); sub_cell[13].push_back(temp[9]);
		sub_cell[13].push_back(temp[12]); sub_cell[13].push_back(temp[15]); sub_cell[13].push_back(temp[11]); sub_cell[13].push_back(temp[8]);
		sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[13]]); sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[14]]);
		sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[10]]); sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[12]]); sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[15]]);
		sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[13].push_back(bss.CP_Local_to_Global[i][temp[8]]);
		//第15个单元
		sub_cell[14].push_back(temp[14]); sub_cell[14].push_back(temp[33]); sub_cell[14].push_back(temp[35]); sub_cell[14].push_back(temp[10]);
		sub_cell[14].push_back(temp[15]); sub_cell[14].push_back(temp[45]); sub_cell[14].push_back(temp[47]); sub_cell[14].push_back(temp[11]);
		sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[14]]); sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[33]]);
		sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[35]]); sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[10]]);
		sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[45]]);
		sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[47]]); sub_cell_global[14].push_back(bss.CP_Local_to_Global[i][temp[11]]);
		//第16个单元
		sub_cell[15].push_back(temp[34]); sub_cell[15].push_back(temp[9]); sub_cell[15].push_back(temp[37]); sub_cell[15].push_back(temp[36]);
		sub_cell[15].push_back(temp[46]); sub_cell[15].push_back(temp[8]); sub_cell[15].push_back(temp[49]); sub_cell[15].push_back(temp[48]);
		sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[34]]); sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[37]]); sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[36]]);
		sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[46]]); sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[8]]);
		sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[49]]); sub_cell_global[15].push_back(bss.CP_Local_to_Global[i][temp[48]]);
		//第17个单元
		sub_cell[16].push_back(temp[9]); sub_cell[16].push_back(temp[10]); sub_cell[16].push_back(temp[38]); sub_cell[16].push_back(temp[37]);
		sub_cell[16].push_back(temp[8]); sub_cell[16].push_back(temp[11]); sub_cell[16].push_back(temp[50]); sub_cell[16].push_back(temp[49]);
		sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[9]]); sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[10]]);
		sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[38]]); sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[37]]);
		sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[8]]); sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[11]]);
		sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[50]]); sub_cell_global[16].push_back(bss.CP_Local_to_Global[i][temp[49]]);
		//第18个单元
		sub_cell[17].push_back(temp[10]); sub_cell[17].push_back(temp[35]); sub_cell[17].push_back(temp[39]); sub_cell[17].push_back(temp[38]);
		sub_cell[17].push_back(temp[11]); sub_cell[17].push_back(temp[47]); sub_cell[17].push_back(temp[51]); sub_cell[17].push_back(temp[50]);
		sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[10]]); sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[35]]);
		sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[39]]); sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[38]]);
		sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[47]]);
		sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[51]]); sub_cell_global[17].push_back(bss.CP_Local_to_Global[i][temp[50]]);
		/*第三层*/
		//第19个单元
		sub_cell[18].push_back(temp[40]); sub_cell[18].push_back(temp[41]); sub_cell[18].push_back(temp[12]); sub_cell[18].push_back(temp[44]);
		sub_cell[18].push_back(temp[4]); sub_cell[18].push_back(temp[52]); sub_cell[18].push_back(temp[55]); sub_cell[18].push_back(temp[54]);
		sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[40]]); sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[41]]);
		sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[12]]); sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[44]]);
		sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[4]]); sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[52]]);
		sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[55]]); sub_cell_global[18].push_back(bss.CP_Local_to_Global[i][temp[54]]);
		//第20个单元
		sub_cell[19].push_back(temp[41]); sub_cell[19].push_back(temp[42]); sub_cell[19].push_back(temp[15]); sub_cell[19].push_back(temp[12]);
		sub_cell[19].push_back(temp[52]); sub_cell[19].push_back(temp[53]); sub_cell[19].push_back(temp[56]); sub_cell[19].push_back(temp[55]);
		sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[41]]); sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[42]]);
		sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[12]]);
		sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[52]]); sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[53]]);
		sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[56]]); sub_cell_global[19].push_back(bss.CP_Local_to_Global[i][temp[55]]);
		//第21个单元
		sub_cell[20].push_back(temp[42]); sub_cell[20].push_back(temp[43]); sub_cell[20].push_back(temp[45]); sub_cell[20].push_back(temp[15]);
		sub_cell[20].push_back(temp[53]); sub_cell[20].push_back(temp[7]); sub_cell[20].push_back(temp[57]); sub_cell[20].push_back(temp[56]);
		sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[42]]); sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[43]]);
		sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[45]]); sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[15]]);
		sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[53]]); sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[7]]);
		sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[57]]); sub_cell_global[20].push_back(bss.CP_Local_to_Global[i][temp[56]]);
		//第22个单元
		sub_cell[21].push_back(temp[44]); sub_cell[21].push_back(temp[12]); sub_cell[21].push_back(temp[8]); sub_cell[21].push_back(temp[46]);
		sub_cell[21].push_back(temp[54]); sub_cell[21].push_back(temp[55]); sub_cell[21].push_back(temp[59]); sub_cell[21].push_back(temp[58]);
		sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[44]]); sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[12]]);
		sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[8]]); sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[46]]);
		sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[54]]); sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[55]]);
		sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[59]]); sub_cell_global[21].push_back(bss.CP_Local_to_Global[i][temp[58]]);
		//第23个单元
		sub_cell[22].push_back(temp[12]); sub_cell[22].push_back(temp[15]); sub_cell[22].push_back(temp[11]); sub_cell[22].push_back(temp[8]);
		sub_cell[22].push_back(temp[55]); sub_cell[22].push_back(temp[56]); sub_cell[22].push_back(temp[60]); sub_cell[22].push_back(temp[59]);
		sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[12]]); sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[15]]);
		sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[8]]);
		sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[55]]); sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[56]]);
		sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[60]]); sub_cell_global[22].push_back(bss.CP_Local_to_Global[i][temp[59]]);
		//第24个单元
		sub_cell[23].push_back(temp[15]); sub_cell[23].push_back(temp[45]); sub_cell[23].push_back(temp[47]); sub_cell[23].push_back(temp[11]);
		sub_cell[23].push_back(temp[56]); sub_cell[23].push_back(temp[57]); sub_cell[23].push_back(temp[61]); sub_cell[23].push_back(temp[60]);
		sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[45]]);
		sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[47]]); sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[11]]);
		sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[56]]); sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[57]]);
		sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[61]]); sub_cell_global[23].push_back(bss.CP_Local_to_Global[i][temp[60]]);
		//第25个单元
		sub_cell[24].push_back(temp[46]); sub_cell[24].push_back(temp[8]); sub_cell[24].push_back(temp[49]); sub_cell[24].push_back(temp[48]);
		sub_cell[24].push_back(temp[58]); sub_cell[24].push_back(temp[59]); sub_cell[24].push_back(temp[62]); sub_cell[24].push_back(temp[0]);
		sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[46]]); sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[8]]);
		sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[49]]); sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[48]]);
		sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[58]]); sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[59]]);
		sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[62]]); sub_cell_global[24].push_back(bss.CP_Local_to_Global[i][temp[0]]);
		//第26个单元
		sub_cell[25].push_back(temp[8]); sub_cell[25].push_back(temp[11]); sub_cell[25].push_back(temp[50]); sub_cell[25].push_back(temp[49]);
		sub_cell[25].push_back(temp[59]); sub_cell[25].push_back(temp[60]); sub_cell[25].push_back(temp[63]); sub_cell[25].push_back(temp[62]);
		sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[8]]); sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[11]]);
		sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[50]]); sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[49]]);
		sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[59]]); sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[60]]);
		sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[63]]); sub_cell_global[25].push_back(bss.CP_Local_to_Global[i][temp[62]]);
		//第27个单元
		sub_cell[26].push_back(temp[11]); sub_cell[26].push_back(temp[47]); sub_cell[26].push_back(temp[51]); sub_cell[26].push_back(temp[50]);
		sub_cell[26].push_back(temp[60]); sub_cell[26].push_back(temp[61]); sub_cell[26].push_back(temp[3]); sub_cell[26].push_back(temp[63]);
		sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[47]]);
		sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[51]]); sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[50]]);
		sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[60]]); sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[61]]);
		sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[3]]); sub_cell_global[26].push_back(bss.CP_Local_to_Global[i][temp[63]]);

		//插入单元
		for (int j = 0; j < 27; j++)
		{
			bss.Local_BSplineSolid[i].add_cell(sub_cell[j]);
			bss.Global_BSplineSolid.add_cell(sub_cell_global[j]);
		}
	}
}

void ApproximationCC::construct_bss1()
{
	uint64_t cell_num = mesh->cells_.size();
	bss.init(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		std::cout << "Schedule of construct Bsplinesolid " << i + 1 << ":" << cell_num << std::endl;
		bss.CP_Local_to_Global[i].resize(64);
		VertexHandle temp[27];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(mesh->vertices_[mesh->cells_[i].vertex_[vn]]);
			//构造局部的小单元的每个点对应于全局模型中的点
			bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
				bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		}

		//第一条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[0]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第二条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[1]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第三条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[2]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第四条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[3]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第五条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[4]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第六条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[5]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第七条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[6]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第八条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[7]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第九条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[8]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第十条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[9]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第十一条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[10]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//第十二条边的中点
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cal_edge_midv(mesh->cells_[i].edges_[11]));
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_face_midv(mesh->cells_[i].faces_[mesh_face_num]);
			temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(face_mid);
			bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
				bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);
		}

		//插入单元的内心
		V3f cell_mid = cal_cell_midv(CellHandle(i));
		temp[tempsize++] = bss.Local_BSplineSolid[i].add_vertex(cell_mid);
		bss.CP_Local_to_Global[i][temp[tempsize - 1]] =
			bss.Global_BSplineSolid.add_vertex(bss.Local_BSplineSolid[i].vertices_[temp[tempsize - 1]]);

		//构建每一个块按顺序排序控制点的下标
		std::vector<std::vector<std::vector<VertexHandle>>> sc(3);
		std::vector<std::vector<VertexHandle>> sf(3);
		std::vector<VertexHandle> sl(3);
		//第一层
		sl[0] = temp[5]; sl[1] = temp[9]; sl[2] = temp[1];
		sf[0] = sl;
		sl[0] = temp[17]; sl[1] = temp[22]; sl[2] = temp[16];
		sf[1] = sl;
		sl[0] = temp[6]; sl[1] = temp[15]; sl[2] = temp[2];
		sf[2] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = temp[10]; sl[1] = temp[20]; sl[2] = temp[8];
		sf[0] = sl;
		sl[0] = temp[24]; sl[1] = temp[26]; sl[2] = temp[25];
		sf[1] = sl;
		sl[0] = temp[14]; sl[1] = temp[21]; sl[2] = temp[12];
		sf[2] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = temp[4]; sl[1] = temp[11]; sl[2] = temp[0];
		sf[0] = sl;
		sl[0] = temp[19]; sl[1] = temp[23]; sl[2] = temp[18];
		sf[1] = sl;
		sl[0] = temp[7]; sl[1] = temp[13]; sl[2] = temp[3];
		sf[2] = sl;
		sc[2] = sf;
		bss.Control_Point_id[i] = sc;

		//构建Bezier体的数据
		std::vector<std::vector<std::vector<V3f>>> sc1(3);
		std::vector<std::vector<V3f>> sf1(3);
		std::vector<V3f> sl1(3);
		//第一层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[5]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[9]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[1]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[17]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[22]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[16]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[6]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[15]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[2]];
		sf1[2] = sl1;
		sc1[0] = sf1;
		//第二层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[10]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[20]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[8]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[24]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[26]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[25]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[14]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[21]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[12]];
		sf1[2] = sl1;
		sc1[1] = sf1;
		//第三层
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[4]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[11]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[0]];
		sf1[0] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[19]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[23]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[18]];
		sf1[1] = sl1;
		sl1[0] = bss.Local_BSplineSolid[i].vertices_[temp[7]]; sl1[1] = bss.Local_BSplineSolid[i].vertices_[temp[13]]; sl1[2] = bss.Local_BSplineSolid[i].vertices_[temp[3]];
		sf1[2] = sl1;
		sc1[2] = sf1;
		SolidData.push_back(sc1);

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

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		std::vector<std::vector<VertexHandle>> sub_cell_global(8);
		//第一个单元
		sub_cell[0].push_back(temp[0]); sub_cell[0].push_back(temp[8]);
		sub_cell[0].push_back(temp[25]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[11]); sub_cell[0].push_back(temp[20]);
		sub_cell[0].push_back(temp[26]); sub_cell[0].push_back(temp[23]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[0]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[8]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[25]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[18]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[20]]);
		sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[26]]); sub_cell_global[0].push_back(bss.CP_Local_to_Global[i][temp[23]]);
		//第二个单元
		sub_cell[1].push_back(temp[8]); sub_cell[1].push_back(temp[1]);
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[25]);
		sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[9]);
		sub_cell[1].push_back(temp[22]); sub_cell[1].push_back(temp[26]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[8]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[1]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[16]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[25]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[20]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[22]]); sub_cell_global[1].push_back(bss.CP_Local_to_Global[i][temp[26]]);
		//第三个单元
		sub_cell[2].push_back(temp[18]); sub_cell[2].push_back(temp[25]);
		sub_cell[2].push_back(temp[12]); sub_cell[2].push_back(temp[3]);
		sub_cell[2].push_back(temp[23]); sub_cell[2].push_back(temp[26]);
		sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[13]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[18]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[25]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[12]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[3]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[23]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[26]]);
		sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[21]]); sub_cell_global[2].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		//第四个单元
		sub_cell[3].push_back(temp[25]); sub_cell[3].push_back(temp[16]);
		sub_cell[3].push_back(temp[2]); sub_cell[3].push_back(temp[12]);
		sub_cell[3].push_back(temp[26]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[15]); sub_cell[3].push_back(temp[21]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[25]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[16]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[2]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[12]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[26]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[22]]);
		sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[3].push_back(bss.CP_Local_to_Global[i][temp[21]]);
		//第五个单元
		sub_cell[4].push_back(temp[11]); sub_cell[4].push_back(temp[20]);
		sub_cell[4].push_back(temp[26]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[4]); sub_cell[4].push_back(temp[10]);
		sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[19]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[11]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[20]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[26]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[23]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[4]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[10]]);
		sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[24]]); sub_cell_global[4].push_back(bss.CP_Local_to_Global[i][temp[19]]);
		//第六个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[9]);
		sub_cell[5].push_back(temp[22]); sub_cell[5].push_back(temp[26]);
		sub_cell[5].push_back(temp[10]); sub_cell[5].push_back(temp[5]);
		sub_cell[5].push_back(temp[17]); sub_cell[5].push_back(temp[24]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[20]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[9]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[22]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[26]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[10]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[5]]);
		sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[17]]); sub_cell_global[5].push_back(bss.CP_Local_to_Global[i][temp[24]]);
		//第七个单元
		sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]);
		sub_cell[6].push_back(temp[21]); sub_cell[6].push_back(temp[13]);
		sub_cell[6].push_back(temp[19]); sub_cell[6].push_back(temp[24]);
		sub_cell[6].push_back(temp[14]); sub_cell[6].push_back(temp[7]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[23]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[26]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[21]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[13]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[19]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[24]]);
		sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[14]]); sub_cell_global[6].push_back(bss.CP_Local_to_Global[i][temp[7]]);
		//第八个单元
		sub_cell[7].push_back(temp[26]); sub_cell[7].push_back(temp[22]);
		sub_cell[7].push_back(temp[15]); sub_cell[7].push_back(temp[21]);
		sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[17]);
		sub_cell[7].push_back(temp[6]); sub_cell[7].push_back(temp[14]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[26]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[22]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[15]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[21]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[24]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[17]]);
		sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[6]]); sub_cell_global[7].push_back(bss.CP_Local_to_Global[i][temp[14]]);

		//插入单元
		for (int j = 0; j < 8; j++)
		{
			bss.Local_BSplineSolid[i].add_cell(sub_cell[j]);
			bss.Global_BSplineSolid.add_cell(sub_cell_global[j]);
		}
	}

}

/*不考虑奇异情况*/
V3f ApproximationCC::cal_inter_point1(FaceHandle _f, VertexHandle _v)
{
	V3f newv;
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
	return newv;
}

V3f ApproximationCC::cal_inter_point1(CellHandle _c, VertexHandle _v)
{
	V3f newv;
	std::set<EdgeHandle, compare_OVM> ve = mesh->neighbor_v[_v];
	std::vector<VertexHandle> temp;
	temp.push_back(_v);
	int sum = 0;
	newv = newv + mesh->vertices_[_v] * 8;
	sum = sum + 8;
	for (std::set<EdgeHandle, compare_OVM>::iterator ve_it = ve.begin(); ve_it != ve.end(); ve_it++)
	{
		if (mesh->is_e_in_c(*ve_it, _c))
		{
			if (mesh->edges_[*ve_it].from_h == _v)
			{
				newv = newv + mesh->vertices_[mesh->edges_[*ve_it].to_h] * 4;
				sum = sum + 4;
				temp.push_back(mesh->edges_[*ve_it].to_h);
			}
			else
			{
				newv = newv + mesh->vertices_[mesh->edges_[*ve_it].from_h] * 4;
				sum = sum + 4;
				temp.push_back(mesh->edges_[*ve_it].from_h);
			}
		}
	}
	std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
	for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin(); vf_it != vf.end(); vf_it++)
	{
		if (mesh->is_f_in_c(*vf_it, _c))
		{
			VertexHandle fopv = mesh->find_op_vertex(_v, *vf_it);
			newv = newv + mesh->vertices_[fopv] * 2;
			sum = sum + 2;
			temp.push_back(fopv);
		}
	}
	std::vector<VertexHandle> cv = mesh->cells_[_c].vertex_;
	for (std::vector<VertexHandle>::iterator cv_it = cv.begin(); cv_it != cv.end(); cv_it++)
	{
		if (std::find(temp.begin(), temp.end(), *cv_it) == temp.end())
		{
			newv = newv + mesh->vertices_[*cv_it];
			sum = sum + 1;
			break;
		}
	}
	newv = newv / sum;
	return newv;
}

V3f ApproximationCC::solid_edge_point1(EdgeHandle _e, VertexHandle _v)
{
	if (mesh->is_bdy(_e))
	{
		V3f edp;
		std::set<FaceHandle, compare_OVM> ef = mesh->neighbor_e[_e];
		int n = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator ef_it = ef.begin(); ef_it != ef.end(); ef_it++)
		{
			if (mesh->is_bdy(*ef_it))
			{
				edp = edp + cal_inter_point1(*ef_it, _v);
				n++;
			}
		}
		edp = edp / n;
		return edp;
	}
	else
	{
		V3f edp;
		std::set<CellHandle> ec = mesh->neighbor_e_c(_e);
		int n = ec.size();
		for (std::set<CellHandle>::iterator ec_it = ec.begin(); ec_it != ec.end(); ec_it++)
		{
			edp = edp + cal_inter_point1(*ec_it, _v);
		}
		edp = edp / n;
		return edp;
	}
}

V3f ApproximationCC::solid_face_point1(FaceHandle _f, VertexHandle _v)
{
	if (mesh->is_bdy(_f))
	{
		V3f fp;
		fp = cal_inter_point1(_f, _v);
		return fp;
	}
	else
	{
		V3f fp;
		std::set<CellHandle, compare_OVM> fc = mesh->neighbor_f[_f];
		for (std::set<CellHandle, compare_OVM>::iterator fc_it = fc.begin(); fc_it != fc.end(); fc_it++)
		{
			fp = fp + cal_inter_point1(*fc_it, _v);
		}
		fp = fp / 2;
		return fp;
	}
}

V3f ApproximationCC::solid_corner_point1(VertexHandle _v)
{
	if (mesh->is_bdy(_v))
	{
		V3f cp;
		std::set<FaceHandle, compare_OVM> vf = mesh->neighbor_v_f[_v];
		int n = 0;
		for (std::set<FaceHandle, compare_OVM>::iterator vf_it = vf.begin(); vf_it != vf.end(); vf_it++)
		{
			if (mesh->is_bdy(*vf_it))
			{
				cp = cp + cal_inter_point1(*vf_it, _v);
				n++;
			}
		}
		cp = cp / n;
		return cp;
	}
	else
	{
		V3f cp;
		std::set<CellHandle, compare_OVM> vc = mesh->neighbor_v_c[_v];
		int n = vc.size();
		for (std::set<CellHandle, compare_OVM>::iterator vc_it = vc.begin(); vc_it != vc.end(); vc_it++)
		{
			cp = cp + cal_inter_point1(*vc_it, _v);
		}
		cp = cp / n;
		return cp;
	}
}

V3f ApproximationCC::solid_edge_point1(int cell_num, int edge_num, int vertex_num)
{
	return solid_edge_point1(mesh->cells_[cell_num].edges_[edge_num], mesh->cells_[cell_num].vertex_[vertex_num]);
}

V3f ApproximationCC::solid_face_point1(int cell_num, int face_num, int vertex_num)
{
	return solid_face_point1(mesh->cells_[cell_num].faces_[face_num], mesh->cells_[cell_num].vertex_[vertex_num]);
}

HexV3fMesh& ApproximationCC::CCapproximation2()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(solid_corner_point1(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(cal_inter_point1(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 9, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 9, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 1, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 7, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 1, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 7, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 8, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 8, 2));
		//第二层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 6, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 4, 2));
		//第三层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 2, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 6, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 4, 3));
		//第4层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 11, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 11, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 5, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 5, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 10, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 10, 3));
		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(27);
		/*第一层*/
		//第1个单元
		sub_cell[0].push_back(temp[5]); sub_cell[0].push_back(temp[16]); sub_cell[0].push_back(temp[19]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[28]); sub_cell[0].push_back(temp[29]); sub_cell[0].push_back(temp[13]); sub_cell[0].push_back(temp[32]);
		//第2个单元
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[17]); sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[19]);
		sub_cell[1].push_back(temp[29]); sub_cell[1].push_back(temp[30]); sub_cell[1].push_back(temp[14]); sub_cell[1].push_back(temp[13]);
		//第3个单元
		sub_cell[2].push_back(temp[17]); sub_cell[2].push_back(temp[6]); sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[20]);
		sub_cell[2].push_back(temp[30]); sub_cell[2].push_back(temp[31]); sub_cell[2].push_back(temp[33]); sub_cell[2].push_back(temp[14]);
		//第4个单元
		sub_cell[3].push_back(temp[18]); sub_cell[3].push_back(temp[19]); sub_cell[3].push_back(temp[23]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[32]); sub_cell[3].push_back(temp[13]); sub_cell[3].push_back(temp[9]); sub_cell[3].push_back(temp[34]);
		//第5个单元
		sub_cell[4].push_back(temp[19]); sub_cell[4].push_back(temp[20]); sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[13]); sub_cell[4].push_back(temp[14]); sub_cell[4].push_back(temp[10]); sub_cell[4].push_back(temp[9]);
		//第6个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[21]); sub_cell[5].push_back(temp[25]); sub_cell[5].push_back(temp[24]);
		sub_cell[5].push_back(temp[14]); sub_cell[5].push_back(temp[33]); sub_cell[5].push_back(temp[35]); sub_cell[5].push_back(temp[10]);
		//第7个单元
		sub_cell[6].push_back(temp[22]); sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]); sub_cell[6].push_back(temp[1]);
		sub_cell[6].push_back(temp[34]); sub_cell[6].push_back(temp[9]); sub_cell[6].push_back(temp[37]); sub_cell[6].push_back(temp[36]);
		//第8个单元
		sub_cell[7].push_back(temp[23]); sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[27]); sub_cell[7].push_back(temp[26]);
		sub_cell[7].push_back(temp[9]); sub_cell[7].push_back(temp[10]); sub_cell[7].push_back(temp[38]); sub_cell[7].push_back(temp[37]);
		//第9个单元
		sub_cell[8].push_back(temp[24]); sub_cell[8].push_back(temp[25]); sub_cell[8].push_back(temp[2]); sub_cell[8].push_back(temp[27]);
		sub_cell[8].push_back(temp[10]); sub_cell[8].push_back(temp[35]); sub_cell[8].push_back(temp[39]); sub_cell[8].push_back(temp[38]);
		/*第二层*/
		//第10个单元
		sub_cell[9].push_back(temp[28]); sub_cell[9].push_back(temp[29]); sub_cell[9].push_back(temp[13]); sub_cell[9].push_back(temp[32]);
		sub_cell[9].push_back(temp[40]); sub_cell[9].push_back(temp[41]); sub_cell[9].push_back(temp[12]); sub_cell[9].push_back(temp[44]);
		//第11个单元
		sub_cell[10].push_back(temp[29]); sub_cell[10].push_back(temp[30]); sub_cell[10].push_back(temp[14]); sub_cell[10].push_back(temp[13]);
		sub_cell[10].push_back(temp[41]); sub_cell[10].push_back(temp[42]); sub_cell[10].push_back(temp[15]); sub_cell[10].push_back(temp[12]);
		//第12个单元
		sub_cell[11].push_back(temp[30]); sub_cell[11].push_back(temp[31]); sub_cell[11].push_back(temp[33]); sub_cell[11].push_back(temp[14]);
		sub_cell[11].push_back(temp[42]); sub_cell[11].push_back(temp[43]); sub_cell[11].push_back(temp[45]); sub_cell[11].push_back(temp[15]);
		//第13个单元
		sub_cell[12].push_back(temp[32]); sub_cell[12].push_back(temp[13]); sub_cell[12].push_back(temp[9]); sub_cell[12].push_back(temp[34]);
		sub_cell[12].push_back(temp[44]); sub_cell[12].push_back(temp[12]); sub_cell[12].push_back(temp[8]); sub_cell[12].push_back(temp[46]);
		//第14个单元
		sub_cell[13].push_back(temp[13]); sub_cell[13].push_back(temp[14]); sub_cell[13].push_back(temp[10]); sub_cell[13].push_back(temp[9]);
		sub_cell[13].push_back(temp[12]); sub_cell[13].push_back(temp[15]); sub_cell[13].push_back(temp[11]); sub_cell[13].push_back(temp[8]);
		//第15个单元
		sub_cell[14].push_back(temp[14]); sub_cell[14].push_back(temp[33]); sub_cell[14].push_back(temp[35]); sub_cell[14].push_back(temp[10]);
		sub_cell[14].push_back(temp[15]); sub_cell[14].push_back(temp[45]); sub_cell[14].push_back(temp[47]); sub_cell[14].push_back(temp[11]);
		//第16个单元
		sub_cell[15].push_back(temp[34]); sub_cell[15].push_back(temp[9]); sub_cell[15].push_back(temp[37]); sub_cell[15].push_back(temp[36]);
		sub_cell[15].push_back(temp[46]); sub_cell[15].push_back(temp[8]); sub_cell[15].push_back(temp[49]); sub_cell[15].push_back(temp[48]);
		//第17个单元
		sub_cell[16].push_back(temp[9]); sub_cell[16].push_back(temp[10]); sub_cell[16].push_back(temp[38]); sub_cell[16].push_back(temp[37]);
		sub_cell[16].push_back(temp[8]); sub_cell[16].push_back(temp[11]); sub_cell[16].push_back(temp[50]); sub_cell[16].push_back(temp[49]);
		//第18个单元
		sub_cell[17].push_back(temp[10]); sub_cell[17].push_back(temp[35]); sub_cell[17].push_back(temp[39]); sub_cell[17].push_back(temp[38]);
		sub_cell[17].push_back(temp[11]); sub_cell[17].push_back(temp[47]); sub_cell[17].push_back(temp[51]); sub_cell[17].push_back(temp[50]);
		/*第三层*/
		//第19个单元
		sub_cell[18].push_back(temp[40]); sub_cell[18].push_back(temp[41]); sub_cell[18].push_back(temp[12]); sub_cell[18].push_back(temp[44]);
		sub_cell[18].push_back(temp[4]); sub_cell[18].push_back(temp[52]); sub_cell[18].push_back(temp[55]); sub_cell[18].push_back(temp[54]);
		//第20个单元
		sub_cell[19].push_back(temp[41]); sub_cell[19].push_back(temp[42]); sub_cell[19].push_back(temp[15]); sub_cell[19].push_back(temp[12]);
		sub_cell[19].push_back(temp[52]); sub_cell[19].push_back(temp[53]); sub_cell[19].push_back(temp[56]); sub_cell[19].push_back(temp[55]);
		//第21个单元
		sub_cell[20].push_back(temp[42]); sub_cell[20].push_back(temp[43]); sub_cell[20].push_back(temp[45]); sub_cell[20].push_back(temp[15]);
		sub_cell[20].push_back(temp[53]); sub_cell[20].push_back(temp[7]); sub_cell[20].push_back(temp[57]); sub_cell[20].push_back(temp[56]);
		//第22个单元
		sub_cell[21].push_back(temp[44]); sub_cell[21].push_back(temp[12]); sub_cell[21].push_back(temp[8]); sub_cell[21].push_back(temp[46]);
		sub_cell[21].push_back(temp[54]); sub_cell[21].push_back(temp[55]); sub_cell[21].push_back(temp[59]); sub_cell[21].push_back(temp[58]);
		//第23个单元
		sub_cell[22].push_back(temp[12]); sub_cell[22].push_back(temp[15]); sub_cell[22].push_back(temp[11]); sub_cell[22].push_back(temp[8]);
		sub_cell[22].push_back(temp[55]); sub_cell[22].push_back(temp[56]); sub_cell[22].push_back(temp[60]); sub_cell[22].push_back(temp[59]);
		//第24个单元
		sub_cell[23].push_back(temp[15]); sub_cell[23].push_back(temp[45]); sub_cell[23].push_back(temp[47]); sub_cell[23].push_back(temp[11]);
		sub_cell[23].push_back(temp[56]); sub_cell[23].push_back(temp[57]); sub_cell[23].push_back(temp[61]); sub_cell[23].push_back(temp[60]);
		//第25个单元
		sub_cell[24].push_back(temp[46]); sub_cell[24].push_back(temp[8]); sub_cell[24].push_back(temp[49]); sub_cell[24].push_back(temp[48]);
		sub_cell[24].push_back(temp[58]); sub_cell[24].push_back(temp[59]); sub_cell[24].push_back(temp[62]); sub_cell[24].push_back(temp[0]);
		//第26个单元
		sub_cell[25].push_back(temp[8]); sub_cell[25].push_back(temp[11]); sub_cell[25].push_back(temp[50]); sub_cell[25].push_back(temp[49]);
		sub_cell[25].push_back(temp[59]); sub_cell[25].push_back(temp[60]); sub_cell[25].push_back(temp[63]); sub_cell[25].push_back(temp[62]);
		//第27个单元
		sub_cell[26].push_back(temp[11]); sub_cell[26].push_back(temp[47]); sub_cell[26].push_back(temp[51]); sub_cell[26].push_back(temp[50]);
		sub_cell[26].push_back(temp[60]); sub_cell[26].push_back(temp[61]); sub_cell[26].push_back(temp[3]); sub_cell[26].push_back(temp[63]);

		//插入单元
		for (int j = 0; j < 27; j++)
		{
			newmesh->add_cell(sub_cell[j]);
		}
	}
	return *newmesh;
}

HexV3fMesh& ApproximationCC::CCapproximation3()
{
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(solid_corner_point1(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = newmesh->add_vertex(cal_inter_point1(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 9, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 9, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 1, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 7, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 1, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 2, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 7, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 8, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 8, 2));
		//第二层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 2, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 6));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 6, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 5));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 6));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 2));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 0, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 1));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 2));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 4, 2));
		//第三层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 2, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 4, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 6, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 1, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 0, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 5, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 4, 3));
		//第4层
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 11, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 11, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 4));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 7));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 5, 7));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_face_point1(i, 3, 3));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 5, 3));

		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 10, 0));
		temp[tempsize++] = newmesh->add_vertex(solid_edge_point1(i, 10, 3));

		//构建Bezier体的数据
		std::vector<std::vector<std::vector<V3f>>> sc(4);
		std::vector<std::vector<V3f>> sf(4);
		std::vector<V3f> sl(4);
		//第一层
		sl[0] = newmesh->vertices_[temp[5]]; sl[1] = newmesh->vertices_[temp[18]];
		sl[2] = newmesh->vertices_[temp[22]]; sl[3] = newmesh->vertices_[temp[1]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[16]]; sl[1] = newmesh->vertices_[temp[19]];
		sl[2] = newmesh->vertices_[temp[23]]; sl[3] = newmesh->vertices_[temp[26]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[17]]; sl[1] = newmesh->vertices_[temp[20]];
		sl[2] = newmesh->vertices_[temp[24]]; sl[3] = newmesh->vertices_[temp[27]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[6]]; sl[1] = newmesh->vertices_[temp[21]];
		sl[2] = newmesh->vertices_[temp[25]]; sl[3] = newmesh->vertices_[temp[2]];
		sf[3] = sl;
		sc[0] = sf;
		//第二层
		sl[0] = newmesh->vertices_[temp[28]]; sl[1] = newmesh->vertices_[temp[32]];
		sl[2] = newmesh->vertices_[temp[34]]; sl[3] = newmesh->vertices_[temp[36]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[29]]; sl[1] = newmesh->vertices_[temp[13]];
		sl[2] = newmesh->vertices_[temp[9]]; sl[3] = newmesh->vertices_[temp[37]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[30]]; sl[1] = newmesh->vertices_[temp[14]];
		sl[2] = newmesh->vertices_[temp[10]]; sl[3] = newmesh->vertices_[temp[38]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[31]]; sl[1] = newmesh->vertices_[temp[33]];
		sl[2] = newmesh->vertices_[temp[35]]; sl[3] = newmesh->vertices_[temp[39]];
		sf[3] = sl;
		sc[1] = sf;
		//第三层
		sl[0] = newmesh->vertices_[temp[40]]; sl[1] = newmesh->vertices_[temp[44]];
		sl[2] = newmesh->vertices_[temp[46]]; sl[3] = newmesh->vertices_[temp[48]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[41]]; sl[1] = newmesh->vertices_[temp[12]];
		sl[2] = newmesh->vertices_[temp[8]]; sl[3] = newmesh->vertices_[temp[49]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[42]]; sl[1] = newmesh->vertices_[temp[15]];
		sl[2] = newmesh->vertices_[temp[11]]; sl[3] = newmesh->vertices_[temp[50]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[43]]; sl[1] = newmesh->vertices_[temp[45]];
		sl[2] = newmesh->vertices_[temp[47]]; sl[3] = newmesh->vertices_[temp[51]];
		sf[3] = sl;
		sc[2] = sf;
		//第四层
		sl[0] = newmesh->vertices_[temp[4]]; sl[1] = newmesh->vertices_[temp[54]];
		sl[2] = newmesh->vertices_[temp[58]]; sl[3] = newmesh->vertices_[temp[0]];
		sf[0] = sl;
		sl[0] = newmesh->vertices_[temp[52]]; sl[1] = newmesh->vertices_[temp[55]];
		sl[2] = newmesh->vertices_[temp[59]]; sl[3] = newmesh->vertices_[temp[62]];
		sf[1] = sl;
		sl[0] = newmesh->vertices_[temp[53]]; sl[1] = newmesh->vertices_[temp[56]];
		sl[2] = newmesh->vertices_[temp[60]]; sl[3] = newmesh->vertices_[temp[63]];
		sf[2] = sl;
		sl[0] = newmesh->vertices_[temp[7]]; sl[1] = newmesh->vertices_[temp[57]];
		sl[2] = newmesh->vertices_[temp[61]]; sl[3] = newmesh->vertices_[temp[3]];
		sf[3] = sl;
		sc[3] = sf;

		SolidData.push_back(sc);
	}
	return *newmesh;
}

std::vector<HexV3fMesh>& ApproximationCC::CCapproximationDAT1()
{
	uint64_t cell_num = mesh->cells_.size();
	vmesh.resize(cell_num);
	for (uint64_t i = 0; i < cell_num; i++)
	{
		VertexHandle temp[64];
		int tempsize = 0;
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = vmesh[i].add_vertex(solid_corner_point1(mesh->cells_[i].vertex_[vn]));
		}
		for (int vn = 0; vn < 8; vn++)
		{
			temp[tempsize++] = vmesh[i].add_vertex(cal_inter_point1(CellHandle(i), mesh->cells_[i].vertex_[vn]));
		}
		//第一层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 9, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 9, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 1, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 2, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 2, 6));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 7, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 1, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 2, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 2, 2));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 7, 2));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 8, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 8, 2));
		//第二层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 2, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 4, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 4, 6));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 6, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 0, 5));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 1, 6));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 0, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 1, 2));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 0, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 5, 1));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 5, 2));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 4, 2));
		//第三层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 2, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 4, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 4, 7));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 6, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 0, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 1, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 0, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 1, 3));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 0, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 5, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 5, 3));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 4, 3));
		//第4层
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 11, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 11, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 3, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 3, 4));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 3, 7));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 5, 7));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 3, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 3, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_face_point1(i, 3, 3));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 5, 3));

		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 10, 0));
		temp[tempsize++] = vmesh[i].add_vertex(solid_edge_point1(i, 10, 3));
		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(27);
		/*第一层*/
		//第1个单元
		sub_cell[0].push_back(temp[5]); sub_cell[0].push_back(temp[16]); sub_cell[0].push_back(temp[19]); sub_cell[0].push_back(temp[18]);
		sub_cell[0].push_back(temp[28]); sub_cell[0].push_back(temp[29]); sub_cell[0].push_back(temp[13]); sub_cell[0].push_back(temp[32]);
		//第2个单元
		sub_cell[1].push_back(temp[16]); sub_cell[1].push_back(temp[17]); sub_cell[1].push_back(temp[20]); sub_cell[1].push_back(temp[19]);
		sub_cell[1].push_back(temp[29]); sub_cell[1].push_back(temp[30]); sub_cell[1].push_back(temp[14]); sub_cell[1].push_back(temp[13]);
		//第3个单元
		sub_cell[2].push_back(temp[17]); sub_cell[2].push_back(temp[6]); sub_cell[2].push_back(temp[21]); sub_cell[2].push_back(temp[20]);
		sub_cell[2].push_back(temp[30]); sub_cell[2].push_back(temp[31]); sub_cell[2].push_back(temp[33]); sub_cell[2].push_back(temp[14]);
		//第4个单元
		sub_cell[3].push_back(temp[18]); sub_cell[3].push_back(temp[19]); sub_cell[3].push_back(temp[23]); sub_cell[3].push_back(temp[22]);
		sub_cell[3].push_back(temp[32]); sub_cell[3].push_back(temp[13]); sub_cell[3].push_back(temp[9]); sub_cell[3].push_back(temp[34]);
		//第5个单元
		sub_cell[4].push_back(temp[19]); sub_cell[4].push_back(temp[20]); sub_cell[4].push_back(temp[24]); sub_cell[4].push_back(temp[23]);
		sub_cell[4].push_back(temp[13]); sub_cell[4].push_back(temp[14]); sub_cell[4].push_back(temp[10]); sub_cell[4].push_back(temp[9]);
		//第6个单元
		sub_cell[5].push_back(temp[20]); sub_cell[5].push_back(temp[21]); sub_cell[5].push_back(temp[25]); sub_cell[5].push_back(temp[24]);
		sub_cell[5].push_back(temp[14]); sub_cell[5].push_back(temp[33]); sub_cell[5].push_back(temp[35]); sub_cell[5].push_back(temp[10]);
		//第7个单元
		sub_cell[6].push_back(temp[22]); sub_cell[6].push_back(temp[23]); sub_cell[6].push_back(temp[26]); sub_cell[6].push_back(temp[1]);
		sub_cell[6].push_back(temp[34]); sub_cell[6].push_back(temp[9]); sub_cell[6].push_back(temp[37]); sub_cell[6].push_back(temp[36]);
		//第8个单元
		sub_cell[7].push_back(temp[23]); sub_cell[7].push_back(temp[24]); sub_cell[7].push_back(temp[27]); sub_cell[7].push_back(temp[26]);
		sub_cell[7].push_back(temp[9]); sub_cell[7].push_back(temp[10]); sub_cell[7].push_back(temp[38]); sub_cell[7].push_back(temp[37]);
		//第9个单元
		sub_cell[8].push_back(temp[24]); sub_cell[8].push_back(temp[25]); sub_cell[8].push_back(temp[2]); sub_cell[8].push_back(temp[27]);
		sub_cell[8].push_back(temp[10]); sub_cell[8].push_back(temp[35]); sub_cell[8].push_back(temp[39]); sub_cell[8].push_back(temp[38]);
		/*第二层*/
		//第10个单元
		sub_cell[9].push_back(temp[28]); sub_cell[9].push_back(temp[29]); sub_cell[9].push_back(temp[13]); sub_cell[9].push_back(temp[32]);
		sub_cell[9].push_back(temp[40]); sub_cell[9].push_back(temp[41]); sub_cell[9].push_back(temp[12]); sub_cell[9].push_back(temp[44]);
		//第11个单元
		sub_cell[10].push_back(temp[29]); sub_cell[10].push_back(temp[30]); sub_cell[10].push_back(temp[14]); sub_cell[10].push_back(temp[13]);
		sub_cell[10].push_back(temp[41]); sub_cell[10].push_back(temp[42]); sub_cell[10].push_back(temp[15]); sub_cell[10].push_back(temp[12]);
		//第12个单元
		sub_cell[11].push_back(temp[30]); sub_cell[11].push_back(temp[31]); sub_cell[11].push_back(temp[33]); sub_cell[11].push_back(temp[14]);
		sub_cell[11].push_back(temp[42]); sub_cell[11].push_back(temp[43]); sub_cell[11].push_back(temp[45]); sub_cell[11].push_back(temp[15]);
		//第13个单元
		sub_cell[12].push_back(temp[32]); sub_cell[12].push_back(temp[13]); sub_cell[12].push_back(temp[9]); sub_cell[12].push_back(temp[34]);
		sub_cell[12].push_back(temp[44]); sub_cell[12].push_back(temp[12]); sub_cell[12].push_back(temp[8]); sub_cell[12].push_back(temp[46]);
		//第14个单元
		sub_cell[13].push_back(temp[13]); sub_cell[13].push_back(temp[14]); sub_cell[13].push_back(temp[10]); sub_cell[13].push_back(temp[9]);
		sub_cell[13].push_back(temp[12]); sub_cell[13].push_back(temp[15]); sub_cell[13].push_back(temp[11]); sub_cell[13].push_back(temp[8]);
		//第15个单元
		sub_cell[14].push_back(temp[14]); sub_cell[14].push_back(temp[33]); sub_cell[14].push_back(temp[35]); sub_cell[14].push_back(temp[10]);
		sub_cell[14].push_back(temp[15]); sub_cell[14].push_back(temp[45]); sub_cell[14].push_back(temp[47]); sub_cell[14].push_back(temp[11]);
		//第16个单元
		sub_cell[15].push_back(temp[34]); sub_cell[15].push_back(temp[9]); sub_cell[15].push_back(temp[37]); sub_cell[15].push_back(temp[36]);
		sub_cell[15].push_back(temp[46]); sub_cell[15].push_back(temp[8]); sub_cell[15].push_back(temp[49]); sub_cell[15].push_back(temp[48]);
		//第17个单元
		sub_cell[16].push_back(temp[9]); sub_cell[16].push_back(temp[10]); sub_cell[16].push_back(temp[38]); sub_cell[16].push_back(temp[37]);
		sub_cell[16].push_back(temp[8]); sub_cell[16].push_back(temp[11]); sub_cell[16].push_back(temp[50]); sub_cell[16].push_back(temp[49]);
		//第18个单元
		sub_cell[17].push_back(temp[10]); sub_cell[17].push_back(temp[35]); sub_cell[17].push_back(temp[39]); sub_cell[17].push_back(temp[38]);
		sub_cell[17].push_back(temp[11]); sub_cell[17].push_back(temp[47]); sub_cell[17].push_back(temp[51]); sub_cell[17].push_back(temp[50]);
		/*第三层*/
		//第19个单元
		sub_cell[18].push_back(temp[40]); sub_cell[18].push_back(temp[41]); sub_cell[18].push_back(temp[12]); sub_cell[18].push_back(temp[44]);
		sub_cell[18].push_back(temp[4]); sub_cell[18].push_back(temp[52]); sub_cell[18].push_back(temp[55]); sub_cell[18].push_back(temp[54]);
		//第20个单元
		sub_cell[19].push_back(temp[41]); sub_cell[19].push_back(temp[42]); sub_cell[19].push_back(temp[15]); sub_cell[19].push_back(temp[12]);
		sub_cell[19].push_back(temp[52]); sub_cell[19].push_back(temp[53]); sub_cell[19].push_back(temp[56]); sub_cell[19].push_back(temp[55]);
		//第21个单元
		sub_cell[20].push_back(temp[42]); sub_cell[20].push_back(temp[43]); sub_cell[20].push_back(temp[45]); sub_cell[20].push_back(temp[15]);
		sub_cell[20].push_back(temp[53]); sub_cell[20].push_back(temp[7]); sub_cell[20].push_back(temp[57]); sub_cell[20].push_back(temp[56]);
		//第22个单元
		sub_cell[21].push_back(temp[44]); sub_cell[21].push_back(temp[12]); sub_cell[21].push_back(temp[8]); sub_cell[21].push_back(temp[46]);
		sub_cell[21].push_back(temp[54]); sub_cell[21].push_back(temp[55]); sub_cell[21].push_back(temp[59]); sub_cell[21].push_back(temp[58]);
		//第23个单元
		sub_cell[22].push_back(temp[12]); sub_cell[22].push_back(temp[15]); sub_cell[22].push_back(temp[11]); sub_cell[22].push_back(temp[8]);
		sub_cell[22].push_back(temp[55]); sub_cell[22].push_back(temp[56]); sub_cell[22].push_back(temp[60]); sub_cell[22].push_back(temp[59]);
		//第24个单元
		sub_cell[23].push_back(temp[15]); sub_cell[23].push_back(temp[45]); sub_cell[23].push_back(temp[47]); sub_cell[23].push_back(temp[11]);
		sub_cell[23].push_back(temp[56]); sub_cell[23].push_back(temp[57]); sub_cell[23].push_back(temp[61]); sub_cell[23].push_back(temp[60]);
		//第25个单元
		sub_cell[24].push_back(temp[46]); sub_cell[24].push_back(temp[8]); sub_cell[24].push_back(temp[49]); sub_cell[24].push_back(temp[48]);
		sub_cell[24].push_back(temp[58]); sub_cell[24].push_back(temp[59]); sub_cell[24].push_back(temp[62]); sub_cell[24].push_back(temp[0]);
		//第26个单元
		sub_cell[25].push_back(temp[8]); sub_cell[25].push_back(temp[11]); sub_cell[25].push_back(temp[50]); sub_cell[25].push_back(temp[49]);
		sub_cell[25].push_back(temp[59]); sub_cell[25].push_back(temp[60]); sub_cell[25].push_back(temp[63]); sub_cell[25].push_back(temp[62]);
		//第27个单元
		sub_cell[26].push_back(temp[11]); sub_cell[26].push_back(temp[47]); sub_cell[26].push_back(temp[51]); sub_cell[26].push_back(temp[50]);
		sub_cell[26].push_back(temp[60]); sub_cell[26].push_back(temp[61]); sub_cell[26].push_back(temp[3]); sub_cell[26].push_back(temp[63]);

		//插入单元
		for (int j = 0; j < 27; j++)
		{
			vmesh[i].add_cell(sub_cell[j]);
		}
	}
	return vmesh;
}

//-----TwoNSubdivision-----//
void TwoNSubdivision::cal_edgepoint_withn()
{
	edgepoint.resize(iteration_num + 1, std::vector<V3f>(mesh->edges_.size(), V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D(iteration_num + 1, std::vector<V3f>(mesh->edges_.size(), V3f(0, 0, 0)));
	for (int i = 1; i <= iteration_num; i++)
	{
		if (i == 1)
		{
			for (int edge_num = 0; edge_num < mesh->edges_.size(); edge_num++)
			{
				/*if (!mesh->is_bdy(EdgeHandle(edge_num)))
					continue;*/
				edgepoint[i][edge_num] = cal_edge_midv(EdgeHandle(edge_num));
				D[i][edge_num] = edgepoint[i][edge_num] - edgepoint[i - 1][edge_num];
			}
		}
		else
		{
			for (int edge_num = 0; edge_num < mesh->edges_.size(); edge_num++)
			{
				if (!mesh->is_bdy(EdgeHandle(edge_num)))
				{
					VertexHandle vfrom = mesh->edges_[edge_num].from_h;
					VertexHandle vto = mesh->edges_[edge_num].to_h;
					std::set<EdgeHandle, compare_OVM> vfrom_e = mesh->neighbor_v[vfrom];
					std::set<EdgeHandle, compare_OVM> vto_e = mesh->neighbor_v[vto];
					int Nvfrom = vfrom_e.size();
					int Nvto = vto_e.size();
					std::vector<EdgeHandle> vfrom_e_sort;
					std::vector<EdgeHandle> vto_e_sort;//放置排序之后的边
					vfrom_e_sort.push_back(EdgeHandle(edge_num));
					vto_e_sort.push_back(EdgeHandle(edge_num));
					EdgeHandle from_e_first;//用于标记一个点的排序之后的邻接边的除了本边之外的第一条边
					V3f tempfrom, tempto;
					int temp_id = 0;
					EdgeHandle opfrom, opto;//用于记录与初始边不在一个块内的边
					for (std::set<EdgeHandle, compare_OVM>::iterator vfrom_e_it = vfrom_e.begin(); vfrom_e_it != vfrom_e.end(); vfrom_e_it++)
					{
						if (!mesh->is_two_e_in_same_cell(EdgeHandle(edge_num), *vfrom_e_it))
							opfrom = *vfrom_e_it;
					}
					for (std::set<EdgeHandle, compare_OVM>::iterator vto_e_it = vto_e.begin(); vto_e_it != vto_e.end(); vto_e_it++)
					{
						if (!mesh->is_two_e_in_same_cell(EdgeHandle(edge_num), *vto_e_it))
							opto = *vto_e_it;
					}
					//进行两个点的邻接边的排序
					bool check = false;//判断时候排序完成
					while (check == false)
					{
						for (std::set<EdgeHandle, compare_OVM>::iterator vfrom_e_it = vfrom_e.begin(); vfrom_e_it != vfrom_e.end(); vfrom_e_it++)
						{
							if (mesh->is_two_e_in_same_cell(*(vfrom_e_sort.end() - 1), *vfrom_e_it) && *vfrom_e_it != *(vfrom_e_sort.end() - 1)
								&& (*vfrom_e_it) != opfrom)
							{
								if (vfrom_e_sort.size() == Nvfrom - 1 && Nvfrom % 2 == 0)
								{
									vfrom_e_sort.push_back(opfrom);
									break;
								}
								else if (vfrom_e_sort.size() > 1)//要与上一条邻接边不同
								{
									if (std::find(vfrom_e_sort.begin(), vfrom_e_sort.end(), *vfrom_e_it) == vfrom_e_sort.end())
									{
										vfrom_e_sort.push_back(*vfrom_e_it);
										break;
									}
								}
								else
								{
									vfrom_e_sort.push_back(*vfrom_e_it);
									from_e_first = *vfrom_e_it;
									break;
								}
							}
						}
						if (vfrom_e_sort.size() == Nvfrom)
							check = true;
					}
					for (std::vector<EdgeHandle>::iterator vfrom_e_sort_it = vfrom_e_sort.begin(); vfrom_e_sort_it != vfrom_e_sort.end();
						vfrom_e_sort_it++)
					{
						double a_from = (4 / Nvfrom * 1.0) * cos(2 * PI * temp_id / Nvfrom * 1.0);
						temp_id++;
						tempfrom = tempfrom + D[i - 1][*vfrom_e_sort_it] * a_from;
					}
					//排序第二个点的邻接边
					check = false;
					while (check == false)
					{
						for (std::set<EdgeHandle, compare_OVM>::iterator vto_e_it = vto_e.begin(); vto_e_it != vto_e.end(); vto_e_it++)
						{
							if (mesh->is_two_e_in_same_cell(*(vto_e_sort.end() - 1), *vto_e_it) && *vto_e_it != *(vto_e_sort.end() - 1)
								&& (*vto_e_it) != opto)
							{
								if (vto_e_sort.size() == Nvto - 1 && Nvto % 2 == 0)
								{
									vto_e_sort.push_back(opto);
									break;
								}
								else if (vto_e_sort.size() > 1)//要与上一条邻接边不同
								{
									if (std::find(vto_e_sort.begin(), vto_e_sort.end(), *vto_e_it) == vto_e_sort.end())
									{
										vto_e_sort.push_back(*vto_e_it);
										break;
									}
								}
								else
								{
									if (mesh->is_two_e_in_same_cell(from_e_first, *vto_e_it))
									{
										vto_e_sort.push_back(*vto_e_it);
										break;
									}
								}
							}
						}
						if (vto_e_sort.size() == Nvto)
							check = true;
					}
					temp_id = 0;
					for (std::vector<EdgeHandle>::iterator vto_e_sort_it = vto_e_sort.begin(); vto_e_sort_it != vto_e_sort.end(); vto_e_sort_it++)
					{
						double a_to = (4 / Nvto * 1.0) * cos(2 * PI * temp_id / Nvto * 1.0);
						temp_id++;
						tempto = tempto + D[i - 1][*vto_e_sort_it] * a_to;
					}
					edgepoint[i][edge_num] = edgepoint[i - 1][edge_num] + (tempfrom + tempto) * ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
					D[i][edge_num] = edgepoint[i][edge_num] - edgepoint[i - 1][edge_num];
				}
				else
				{
					VertexHandle vfrom = mesh->edges_[edge_num].from_h;
					VertexHandle vto = mesh->edges_[edge_num].to_h;
					std::set<EdgeHandle> vfrom_e = mesh->neighbor_v_bdy(vfrom);
					std::set<EdgeHandle> vto_e = mesh->neighbor_v_bdy(vto);
					int Nvfrom = vfrom_e.size();
					int Nvto = vto_e.size();
					std::vector<EdgeHandle> vfrom_e_sort;
					std::vector<EdgeHandle> vto_e_sort;//放置排序之后的边
					vfrom_e_sort.push_back(EdgeHandle(edge_num));
					vto_e_sort.push_back(EdgeHandle(edge_num));
					EdgeHandle from_e_first;//用于标记一个点的排序之后的邻接边的除了本边之外的第一条边
					V3f tempfrom, tempto;
					int temp_id = 0;
					//进行两个点的邻接边的排序
					bool check = false;//判断时候排序完成
					while (check == false)
					{
						for (std::set<EdgeHandle>::iterator vfrom_e_it = vfrom_e.begin(); vfrom_e_it != vfrom_e.end(); vfrom_e_it++)
						{
							//EdgeHandle edgetemp = *(vfrom_e_sort.end() - 1);
							//EdgeHandle edgetemp1 = *vfrom_e_it;
							if (mesh->is_two_e_in_same_face(*(vfrom_e_sort.end() - 1), *vfrom_e_it) && *vfrom_e_it != *(vfrom_e_sort.end() - 1))
							{
								if (vfrom_e_sort.size() > 1)//要与上一条邻接边不同
								{
									if (*vfrom_e_it != *(vfrom_e_sort.end() - 2)  /*std::find(vfrom_e_sort.begin(), vfrom_e_sort.end(),*vfrom_e_it) == vfrom_e_sort.end()*/)
									{
										vfrom_e_sort.push_back(*vfrom_e_it);
										break;
									}
								}
								else
								{
									vfrom_e_sort.push_back(*vfrom_e_it);
									from_e_first = *vfrom_e_it;
									break;
								}
							}
						}
						if (vfrom_e_sort.size() == Nvfrom)
							check = true;
					}
					for (std::vector<EdgeHandle>::iterator vfrom_e_sort_it = vfrom_e_sort.begin(); vfrom_e_sort_it != vfrom_e_sort.end();
						vfrom_e_sort_it++)
					{
						double a_from = (4 / Nvfrom * 1.0) * cos(2 * PI * temp_id / Nvfrom * 1.0);
						temp_id++;
						tempfrom = tempfrom + D[i - 1][*vfrom_e_sort_it] * a_from;
					}
					//排序第二个点的邻接边
					check = false;
					while (check == false)
					{
						for (std::set<EdgeHandle>::iterator vto_e_it = vto_e.begin(); vto_e_it != vto_e.end(); vto_e_it++)
						{
							if (mesh->is_two_e_in_same_face(*(vto_e_sort.end() - 1), *vto_e_it) && *vto_e_it != *(vto_e_sort.end() - 1))
							{
								if (vto_e_sort.size() > 1)//要与上一条邻接边不同
								{
									if (*vto_e_it != *(vto_e_sort.end() - 2)  /*std::find(vto_e_sort.begin(), vto_e_sort.end(), *vto_e_it) == vto_e_sort.end()*/)
									{
										vto_e_sort.push_back(*vto_e_it);
										break;
									}
								}
								else
								{
									if (mesh->is_two_e_in_same_face(from_e_first, *vto_e_it))
									{
										vto_e_sort.push_back(*vto_e_it);
										break;
									}
								}
							}
						}
						if (vto_e_sort.size() == Nvto)
							check = true;
					}
					temp_id = 0;
					for (std::vector<EdgeHandle>::iterator vto_e_sort_it = vto_e_sort.begin(); vto_e_sort_it != vto_e_sort.end(); vto_e_sort_it++)
					{
						double a_to = (4 / Nvto * 1.0) * cos(2 * PI * temp_id / Nvto * 1.0);
						temp_id++;
						tempto = tempto + D[i - 1][*vto_e_sort_it] * a_to;
					}
					/*V3f t1 = tempfrom + tempto;
					double t2 = (2 * i - 1) / (8 * i * 1.0);
					V3f temp = t1*t2;*/
					edgepoint[i][edge_num] = edgepoint[i - 1][edge_num] + (tempfrom + tempto) * ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
					D[i][edge_num] = edgepoint[i][edge_num] - edgepoint[i - 1][edge_num];
				}
			}
		}
	}
}

V3f TwoNSubdivision::cal_facepoint_withn(FaceHandle _f)
{
	V3f face_point;
	std::vector<V3f> f_direction1 = select_face_line(_f, mesh->faces_[_f].edges_[0]);
	int f_direction1_num = f_direction1.size();
	std::vector<V3f> f_direction2 = select_face_line(_f, mesh->faces_[_f].edges_[1]);
	int f_direction2_num = f_direction2.size();
	//计算第一个方向
	std::vector<std::vector<V3f>> facepoint1(iteration_num + 1, std::vector<V3f>(f_direction1_num, V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D1(iteration_num + 1, std::vector<V3f>(f_direction1_num, V3f(0, 0, 0)));
	for (int i = 1; i <= iteration_num; i++)
	{
		if (i == 1)
		{
			for (int face_num = 0; face_num < f_direction1_num; face_num++)
			{
				facepoint1[i][face_num] = (f_direction1[face_num] + f_direction1[(face_num + 1) % f_direction1_num]) / 2;
				D1[i][face_num] = facepoint1[i][face_num] - facepoint1[i - 1][face_num];
			}
		}
		else
		{
			for (int face_num = 0; face_num < f_direction1_num; face_num++)
			{
				facepoint1[i][face_num] = facepoint1[i - 1][face_num] + (D1[i - 1][face_num] * 2
					- D1[i - 1][(face_num - 1 + f_direction1_num) % f_direction1_num] - D1[i - 1][(face_num + 1) % f_direction1_num])
					* ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
				D1[i][face_num] = facepoint1[i][face_num] - facepoint1[i - 1][face_num];
			}
		}
	}
	//计算第二个方向
	std::vector<std::vector<V3f>> facepoint2(iteration_num + 1, std::vector<V3f>(f_direction2_num, V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D2(iteration_num + 1, std::vector<V3f>(f_direction2_num, V3f(0, 0, 0)));
	for (int i = 1; i <= iteration_num; i++)
	{
		if (i == 1)
		{
			for (int face_num = 0; face_num < f_direction2_num; face_num++)
			{
				facepoint2[i][face_num] = (f_direction2[face_num] + f_direction2[(face_num + 1) % f_direction2_num]) / 2;
				D2[i][face_num] = facepoint2[i][face_num] - facepoint2[i - 1][face_num];
			}
		}
		else
		{
			for (int face_num = 0; face_num < f_direction2_num; face_num++)
			{
				facepoint2[i][face_num] = facepoint2[i - 1][face_num] + (D2[i - 1][face_num] * 2
					- D2[i - 1][(face_num - 1 + f_direction2_num) % f_direction2_num] - D2[i - 1][(face_num + 1) % f_direction2_num])
					* ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
				D2[i][face_num] = facepoint2[i][face_num] - facepoint2[i - 1][face_num];
			}
		}
	}
	face_point = (facepoint1[iteration_num][f_direction1_num - 1] + facepoint2[iteration_num][f_direction2_num - 1]) / 2;
	return face_point;
}

V3f TwoNSubdivision::cal_facepoint_withn(FaceHandle _f, int n)
{
	V3f face_point;
	std::vector<V3f> f_direction1 = select_face_line(_f, mesh->faces_[_f].edges_[0]);
	int f_direction1_num = f_direction1.size();
	std::vector<V3f> f_direction2 = select_face_line(_f, mesh->faces_[_f].edges_[1]);
	int f_direction2_num = f_direction2.size();
	//计算第一个方向
	std::vector<std::vector<V3f>> facepoint1(n + 1, std::vector<V3f>(f_direction1_num, V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D1(n + 1, std::vector<V3f>(f_direction1_num, V3f(0, 0, 0)));
	for (int i = 1; i <= n; i++)
	{
		if (i == 1)
		{
			for (int face_num = 0; face_num < f_direction1_num; face_num++)
			{
				facepoint1[i][face_num] = (f_direction1[face_num] + f_direction1[(face_num + 1) % f_direction1_num]) / 2;
				D1[i][face_num] = facepoint1[i][face_num] - facepoint1[i - 1][face_num];
			}
		}
		else
		{
			for (int face_num = 0; face_num < f_direction1_num; face_num++)
			{
				facepoint1[i][face_num] = facepoint1[i - 1][face_num] + (D1[i - 1][face_num] * 2
					- D1[i - 1][(face_num - 1 + f_direction1_num) % f_direction1_num] - D1[i - 1][(face_num + 1) % f_direction1_num])
					* ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
				D1[i][face_num] = facepoint1[i][face_num] - facepoint1[i - 1][face_num];
			}
		}
	}
	//计算第二个方向
	std::vector<std::vector<V3f>> facepoint2(n + 1, std::vector<V3f>(f_direction2_num, V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D2(n + 1, std::vector<V3f>(f_direction2_num, V3f(0, 0, 0)));
	for (int i = 1; i <= n; i++)
	{
		if (i == 1)
		{
			for (int face_num = 0; face_num < f_direction2_num; face_num++)
			{
				facepoint2[i][face_num] = (f_direction2[face_num] + f_direction2[(face_num + 1) % f_direction2_num]) / 2;
				D2[i][face_num] = facepoint2[i][face_num] - facepoint2[i - 1][face_num];
			}
		}
		else
		{
			for (int face_num = 0; face_num < f_direction2_num; face_num++)
			{
				facepoint2[i][face_num] = facepoint2[i - 1][face_num] + (D2[i - 1][face_num] * 2
					- D2[i - 1][(face_num - 1 + f_direction2_num) % f_direction2_num] - D2[i - 1][(face_num + 1) % f_direction2_num])
					* ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
				D2[i][face_num] = facepoint2[i][face_num] - facepoint2[i - 1][face_num];
			}
		}
	}
	face_point = (facepoint1[n][f_direction1_num - 1] + facepoint2[n][f_direction2_num - 1]) / 2;
	return face_point;
}

std::vector<V3f> TwoNSubdivision::select_face_line(FaceHandle _f, EdgeHandle _e)
{
	std::vector<V3f> face_line;
	EdgeHandle last_edge = mesh->get_opposite_edge_in_face(_e, _f);
	face_line.push_back(edgepoint[iteration_num][_e]);
	FaceHandle previous_f = _f;
	EdgeHandle previous_e = _e;
	while (previous_e != last_edge)
	{
		std::set<FaceHandle, compare_OVM> e_f = mesh->neighbor_e[previous_e];
		for (std::set<FaceHandle, compare_OVM>::iterator e_f_it = e_f.begin(); e_f_it != e_f.end(); e_f_it++)
		{
			if (mesh->is_bdy(*e_f_it) && *e_f_it != previous_f)
			{
				EdgeHandle op_e = mesh->get_opposite_edge_in_face(previous_e, *e_f_it);
				face_line.push_back(edgepoint[iteration_num][op_e]);
				previous_e = op_e;
				previous_f = *e_f_it;
				break;
			}
		}
	}
	return face_line;
}

void TwoNSubdivision::cal_facepoint_withn_solid()
{
	facepointsolid.resize(iteration_num + 1, std::vector<V3f>(mesh->faces_.size(), V3f(0, 0, 0)));
	std::vector<std::vector<V3f>> D(iteration_num + 1, std::vector<V3f>(mesh->faces_.size(), V3f(0, 0, 0)));
	for (int i = 1; i <= iteration_num; i++)
	{
		if (i == 1)
		{
			for (int face_num = 0; face_num < mesh->faces_.size(); face_num++)
			{
				facepointsolid[i][face_num] = cal_face_midv(FaceHandle(face_num));
				D[i][face_num] = facepointsolid[i][face_num] - facepointsolid[i - 1][face_num];
			}
		}
		else
		{
			for (int face_num = 0; face_num < mesh->faces_.size(); face_num++)
			{
				if (mesh->is_bdy(FaceHandle(face_num)))
				{
					V3f t_result;
					std::vector<EdgeHandle> f_e = mesh->faces_[face_num].edges_;
					for (std::vector<EdgeHandle>::iterator f_e_it = f_e.begin(); f_e_it != f_e.end(); f_e_it++)
					{
						std::set<FaceHandle, compare_OVM> f_e_f = mesh->neighbor_e[*f_e_it];
						for (std::set<FaceHandle, compare_OVM>::iterator f_e_f_it = f_e_f.begin(); f_e_f_it != f_e_f.end(); f_e_f_it++)
						{
							if (mesh->is_bdy(*f_e_f_it) && *f_e_f_it != FaceHandle(face_num))
							{
								t_result = t_result - D[i - 1][*f_e_f_it];
								break;
							}
						}
					}
					t_result = t_result + D[i - 1][face_num] * 4;
					facepointsolid[i][face_num] = facepointsolid[i - 1][face_num] + t_result * ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
					D[i][face_num] = facepointsolid[i][face_num] - facepointsolid[i - 1][face_num];
					continue;
				}
				std::vector<std::vector<FaceHandle>> fe_sort;
				V3f temp_result;
				FaceHandle first_face;
				std::vector<EdgeHandle> f_e = mesh->faces_[face_num].edges_;
				int j = 0;
				for (std::vector<EdgeHandle>::iterator f_e_it = f_e.begin(); f_e_it != f_e.end(); f_e_it++)
				{
					std::set<FaceHandle, compare_OVM> f_e_f = mesh->neighbor_e[*f_e_it];
					std::vector<FaceHandle> fe_sort_temp;
					fe_sort_temp.push_back(FaceHandle(face_num));
					bool check = false;
					while (check == false)
					{
						for (std::set<FaceHandle, compare_OVM>::iterator f_e_f_it = f_e_f.begin(); f_e_f_it != f_e_f.end(); f_e_f_it++)
						{
							if (mesh->is_two_f_in_same_cell(*(fe_sort_temp.end() - 1), *f_e_f_it) && *f_e_f_it != *(fe_sort_temp.end() - 1))
							{
								if (fe_sort_temp.size() > 1)//要与上一条邻接面不同
								{
									if (*f_e_f_it != *(fe_sort_temp.end() - 2))
									{
										fe_sort_temp.push_back(*f_e_f_it);
										break;
									}
								}
								else
								{
									if (j == 0)
									{
										first_face = *f_e_f_it;
										fe_sort_temp.push_back(*f_e_f_it);
										j++;
										break;
									}
									else
									{
										if (mesh->is_two_f_in_same_cell(first_face, *f_e_f_it))
										{
											fe_sort_temp.push_back(*f_e_f_it);
											break;
										}
									}
								}
							}
							else
							{
								if (fe_sort_temp.size() == (f_e_f.size() - 1)
									&& std::find(fe_sort_temp.begin(), fe_sort_temp.end(), *f_e_f_it) == fe_sort_temp.end())
								{
									fe_sort_temp.push_back(*f_e_f_it);
									break;
								}
							}
						}
						if (fe_sort_temp.size() == f_e_f.size())
						{
							check = true;
						}
					}
					fe_sort.push_back(fe_sort_temp);
				}
				for (std::vector<std::vector<FaceHandle>>::iterator fe_sort_it = fe_sort.begin(); fe_sort_it != fe_sort.end(); fe_sort_it++)
				{
					std::vector<FaceHandle> temp = *fe_sort_it;
					int temp_id = 0;
					int N = temp.size();
					V3f tempto;
					for (std::vector<FaceHandle>::iterator temp_it = temp.begin(); temp_it != temp.end(); temp_it++)
					{
						double a = (4 / N * 1.0) * cos(2 * PI * temp_id / N * 1.0);
						temp_id++;
						tempto = tempto + D[i - 1][*temp_it] * a;
					}
					temp_result = temp_result + tempto;
				}
				facepointsolid[i][face_num] = facepointsolid[i - 1][face_num] + temp_result * ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
				D[i][face_num] = facepointsolid[i][face_num] - facepointsolid[i - 1][face_num];
			}
		}
	}
}

V3f TwoNSubdivision::cal_twon_edge_point(const int cell_num, const int ver1, const int ver2)
{
	EdgeHandle edge = check_edge(cell_num, ver1, ver2);
	//V3f edge1 = edgepoint[iteration_num][edge];
	//V3f edge2 = cal_edge_mid(edge);
	return edgepoint[iteration_num][edge];
	/*if (mesh->is_bdy(edge))
		return edgepoint[iteration_num][edge];
	else
		return cal_edge_mid(edge);*/
}

V3f TwoNSubdivision::cal_twon_face_point(const int cell_num, const int face_num)
{
	FaceHandle face = mesh->cells_[cell_num].faces_[face_num];//确定该面
	/*return facepointsolid[iteration_num][face];*/
	if (mesh->is_bdy(face))
		return cal_facepoint_withn(face);
	else
		return facepointsolid[iteration_num][face];
}

V3f TwoNSubdivision::cal_twon_new_vertex(const int cell_num, const int ver_num)
{
	VertexHandle vertex = mesh->cells_[cell_num].vertex_[ver_num];
	return mesh->vertices_[vertex];
	/*if (mesh->is_bdy(vertex))
		return mesh->vertices_[vertex];
	else
		return cal_new_vertex(vertex);*/
}

V3f TwoNSubdivision::cal_cellpoint_withn(CellHandle _c)
{
	V3f cell_point;
	std::vector<std::vector<V3f>> cell_line(3);//一个块的3个方向
	cell_line[0] = select_cell_line(_c, mesh->cells_[_c].faces_[0]);
	cell_line[1] = select_cell_line(_c, mesh->cells_[_c].faces_[2]);
	cell_line[2] = select_cell_line(_c, mesh->cells_[_c].faces_[4]);
	for (std::vector<std::vector<V3f>>::iterator cell_line_it = cell_line.begin(); cell_line_it != cell_line.end(); cell_line_it++)
	{
		std::vector<V3f> cellline_temp = *cell_line_it;
		int temp_num = cellline_temp.size();
		std::vector<std::vector<V3f>> cellpoint(iteration_num + 1, std::vector<V3f>(temp_num, V3f(0, 0, 0)));
		std::vector<std::vector<V3f>> D(iteration_num + 1, std::vector<V3f>(temp_num, V3f(0, 0, 0)));
		for (int i = 1; i <= iteration_num; i++)
		{
			if (i == 1)
			{
				for (int cell_num = 0; cell_num < temp_num; cell_num++)
				{
					cellpoint[i][cell_num] = (cellline_temp[cell_num] + cellline_temp[(cell_num + 1) % temp_num]) / 2;
					D[i][cell_num] = cellpoint[i][cell_num] - cellpoint[i - 1][cell_num];
				}
			}
			else
			{
				for (int cell_num = 0; cell_num < temp_num; cell_num++)
				{
					cellpoint[i][cell_num] = cellpoint[i - 1][cell_num] + (D[i - 1][cell_num] * 2
						- D[i - 1][(cell_num - 1 + temp_num) % temp_num] - D[i - 1][(cell_num + 1) % temp_num])
						* ((2 * (i - 1) - 1) / (8 * (i - 1) * 1.0));
					D[i][cell_num] = cellpoint[i][cell_num] - cellpoint[i - 1][cell_num];
				}
			}
		}
		cell_point = cell_point + cellpoint[iteration_num][temp_num - 1];
	}
	cell_point = cell_point / 3;
	return cell_point;
}

std::vector<V3f> TwoNSubdivision::select_cell_line(CellHandle _c, FaceHandle _f)
{
	std::vector<V3f> cell_line;
	std::vector<FaceHandle> check_face;
	FaceHandle last_face = mesh->get_opposite_face_in_cell(_f, _c);
	if (mesh->is_bdy(_f))
		cell_line.push_back(cal_facepoint_withn(_f));
	else
		cell_line.push_back(facepointsolid[iteration_num][_f]);
	check_face.push_back(_f);
	//cell_line.push_back(facepointsolid[iteration_num][_f]);
	CellHandle previous_c = _c;
	FaceHandle previous_f = _f;
	while (previous_f != last_face)
	{
		std::set<CellHandle, compare_OVM> f_c = mesh->neighbor_f[previous_f];
		for (std::set<CellHandle, compare_OVM>::iterator f_c_it = f_c.begin(); f_c_it != f_c.end(); f_c_it++)
		{
			if (*f_c_it != previous_c)
			{
				FaceHandle op_f = mesh->get_opposite_face_in_cell(previous_f, *f_c_it);
				if (mesh->is_bdy(op_f))
					cell_line.push_back(cal_facepoint_withn(op_f));
				else
					cell_line.push_back(facepointsolid[iteration_num][op_f]);
				check_face.push_back(op_f);
				previous_f = op_f;
				previous_c = *f_c_it;
				break;
			}
		}
		if (mesh->is_bdy(previous_f))
		{
			std::vector<V3f> temp_cell_line;
			std::vector<FaceHandle> temp_check;
			if (mesh->is_bdy(last_face))
				temp_cell_line.insert(temp_cell_line.begin(), cal_facepoint_withn(last_face));
			else
				temp_cell_line.insert(temp_cell_line.begin(), facepointsolid[iteration_num][last_face]);
			temp_check.insert(temp_check.begin(), last_face);
			//temp_cell_line.insert(temp_cell_line.begin(),facepointsolid[iteration_num][last_face]);
			CellHandle previous_c_1 = _c;
			FaceHandle previous_f_1 = last_face;
			while (!mesh->is_bdy(previous_f_1))
			{
				std::set<CellHandle, compare_OVM> f_c_1 = mesh->neighbor_f[previous_f_1];
				for (std::set<CellHandle, compare_OVM>::iterator f_c_1_it = f_c_1.begin(); f_c_1_it != f_c_1.end(); f_c_1_it++)
				{
					if (*f_c_1_it != previous_c_1)
					{
						FaceHandle op_f = mesh->get_opposite_face_in_cell(previous_f_1, *f_c_1_it);
						if (mesh->is_bdy(op_f))
							temp_cell_line.insert(temp_cell_line.begin(), cal_facepoint_withn(op_f));
						else
							temp_cell_line.insert(temp_cell_line.begin(), facepointsolid[iteration_num][op_f]);
						temp_check.insert(temp_check.begin(), op_f);
						previous_f_1 = op_f;
						previous_c_1 = *f_c_1_it;
						break;
					}
				}
			}
			cell_line.insert(cell_line.end(), temp_cell_line.begin(), temp_cell_line.end());
			check_face.insert(check_face.end(), temp_check.begin(), temp_check.end());
			previous_f = last_face;
		}
	}
	return cell_line;
}

HexV3fMesh& TwoNSubdivision::TwoNSubdiving()
{
	cal_edgepoint_withn();
	cal_facepoint_withn_solid();
	uint64_t cell_num = mesh->cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{

		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			V3f new_vertex = cal_twon_new_vertex(i, mesh_vex_num);
			temp_vex[mesh_vex_num] = newmesh->add_vertex(new_vertex);
		}
		//-------插入边的中点-------//
		//第一条边的中点
		V3f edge_mid_1 = cal_twon_edge_point(i, 0, 1);
		temp_vex[8] = newmesh->add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2 = cal_twon_edge_point(i, 1, 5);
		temp_vex[9] = newmesh->add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3 = cal_twon_edge_point(i, 5, 4);
		temp_vex[10] = newmesh->add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4 = cal_twon_edge_point(i, 4, 0);
		temp_vex[11] = newmesh->add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5 = cal_twon_edge_point(i, 3, 2);
		temp_vex[12] = newmesh->add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6 = cal_twon_edge_point(i, 3, 7);
		temp_vex[13] = newmesh->add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7 = cal_twon_edge_point(i, 6, 7);
		temp_vex[14] = newmesh->add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8 = cal_twon_edge_point(i, 2, 6);
		temp_vex[15] = newmesh->add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9 = cal_twon_edge_point(i, 1, 2);
		temp_vex[16] = newmesh->add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10 = cal_twon_edge_point(i, 5, 6);
		temp_vex[17] = newmesh->add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11 = cal_twon_edge_point(i, 0, 3);
		temp_vex[18] = newmesh->add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12 = cal_twon_edge_point(i, 4, 7);
		temp_vex[19] = newmesh->add_vertex(edge_mid_12);
		//-------插入边的中点结束-------//

		//-------插入面的中心-------//
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid = cal_twon_face_point(i, mesh_face_num);
			//判断加入数组
			temp_vex[20 + mesh_face_num] = newmesh->add_vertex(face_mid);
		}
		//-------插入面的中心结束-------//

		//-------插入单元的内心-------//
		V3f cell_mid = cal_cellpoint_withn(CellHandle(i));
		//判断加入数组
		temp_vex[26] = newmesh->add_vertex(cell_mid);
		//-------插入单元的内心结束-------//

		//接下来 插入分割后的子单元
		std::vector<std::vector<VertexHandle>> sub_cell(8);
		//第一个单元
		sub_cell[0].push_back(temp_vex[0]); sub_cell[0].push_back(temp_vex[8]);
		sub_cell[0].push_back(temp_vex[25]); sub_cell[0].push_back(temp_vex[18]);
		sub_cell[0].push_back(temp_vex[11]); sub_cell[0].push_back(temp_vex[20]);
		sub_cell[0].push_back(temp_vex[26]); sub_cell[0].push_back(temp_vex[23]);
		//第二个单元
		sub_cell[1].push_back(temp_vex[8]); sub_cell[1].push_back(temp_vex[1]);
		sub_cell[1].push_back(temp_vex[16]); sub_cell[1].push_back(temp_vex[25]);
		sub_cell[1].push_back(temp_vex[20]); sub_cell[1].push_back(temp_vex[9]);
		sub_cell[1].push_back(temp_vex[22]); sub_cell[1].push_back(temp_vex[26]);
		//第三个单元
		sub_cell[2].push_back(temp_vex[18]); sub_cell[2].push_back(temp_vex[25]);
		sub_cell[2].push_back(temp_vex[12]); sub_cell[2].push_back(temp_vex[3]);
		sub_cell[2].push_back(temp_vex[23]); sub_cell[2].push_back(temp_vex[26]);
		sub_cell[2].push_back(temp_vex[21]); sub_cell[2].push_back(temp_vex[13]);
		//第四个单元
		sub_cell[3].push_back(temp_vex[25]); sub_cell[3].push_back(temp_vex[16]);
		sub_cell[3].push_back(temp_vex[2]); sub_cell[3].push_back(temp_vex[12]);
		sub_cell[3].push_back(temp_vex[26]); sub_cell[3].push_back(temp_vex[22]);
		sub_cell[3].push_back(temp_vex[15]); sub_cell[3].push_back(temp_vex[21]);
		//第五个单元
		sub_cell[4].push_back(temp_vex[11]); sub_cell[4].push_back(temp_vex[20]);
		sub_cell[4].push_back(temp_vex[26]); sub_cell[4].push_back(temp_vex[23]);
		sub_cell[4].push_back(temp_vex[4]); sub_cell[4].push_back(temp_vex[10]);
		sub_cell[4].push_back(temp_vex[24]); sub_cell[4].push_back(temp_vex[19]);
		//第六个单元
		sub_cell[5].push_back(temp_vex[20]); sub_cell[5].push_back(temp_vex[9]);
		sub_cell[5].push_back(temp_vex[22]); sub_cell[5].push_back(temp_vex[26]);
		sub_cell[5].push_back(temp_vex[10]); sub_cell[5].push_back(temp_vex[5]);
		sub_cell[5].push_back(temp_vex[17]); sub_cell[5].push_back(temp_vex[24]);
		//第七个单元
		sub_cell[6].push_back(temp_vex[23]); sub_cell[6].push_back(temp_vex[26]);
		sub_cell[6].push_back(temp_vex[21]); sub_cell[6].push_back(temp_vex[13]);
		sub_cell[6].push_back(temp_vex[19]); sub_cell[6].push_back(temp_vex[24]);
		sub_cell[6].push_back(temp_vex[14]); sub_cell[6].push_back(temp_vex[7]);
		//第八个单元
		sub_cell[7].push_back(temp_vex[26]); sub_cell[7].push_back(temp_vex[22]);
		sub_cell[7].push_back(temp_vex[15]); sub_cell[7].push_back(temp_vex[21]);
		sub_cell[7].push_back(temp_vex[24]); sub_cell[7].push_back(temp_vex[17]);
		sub_cell[7].push_back(temp_vex[6]); sub_cell[7].push_back(temp_vex[14]);
		//插入单元
		newmesh->add_cell(sub_cell[0]);
		newmesh->add_cell(sub_cell[1]);
		newmesh->add_cell(sub_cell[2]);
		newmesh->add_cell(sub_cell[3]);
		newmesh->add_cell(sub_cell[4]);
		newmesh->add_cell(sub_cell[5]);
		newmesh->add_cell(sub_cell[6]);
		newmesh->add_cell(sub_cell[7]);
	}
	return *newmesh;
}
