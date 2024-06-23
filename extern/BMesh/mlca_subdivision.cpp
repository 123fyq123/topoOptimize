#include "mlca_subdivision.h"

namespace BMesh{
/**************************************************
@brief   : 感觉插入重复的点了或者没有，因为六面体网格就是由一个一个cell构成的，cell是其单位
		   如果这么来看的话，那么点就应该重复。所以我在ovm程序中遍历一个顶点只有3条边而不是很多条边。
@author  : lee
@input   ：none
@output  ：none
@time    : none
**************************************************/
//-----MLCASubdivision-----//
HexV3fMesh& MLCASubdivision::Spliting(HexV3fMesh &mesh)
{
	HexV3fMesh newmesh;
	uint64_t cell_num = mesh.cells_.size();
	for (uint64_t i = 0; i < cell_num; i++)
	{
		//暴力方法 将每个网格的各点顺序固定死 插入27个点
		VertexHandle temp_vex[27];
		std::vector<V3f>::iterator v_id;

		// 傻逼了 直接for不就好了

		//首先插入原始的8个点
		for (int mesh_vex_num = 0; mesh_vex_num < 8; mesh_vex_num++)
		{
			temp_vex[mesh_vex_num] = newmesh.add_vertex(mesh.vertices_[mesh.cells_[i].vertex_[mesh_vex_num]]);
		}
		//再插入每条边的中点 暴力插入 以免顺序错误
		//第一条边的中点
		V3f edge_mid_1;
		edge_mid_1 = (mesh.vertices_[mesh.cells_[i].vertex_[0]] + mesh.vertices_[mesh.cells_[i].vertex_[1]]) * 1.0 / 2 * 1.0;
		temp_vex[8] = newmesh.add_vertex(edge_mid_1);

		//第二条边的中点
		V3f edge_mid_2;
		edge_mid_2 = (mesh.vertices_[mesh.cells_[i].vertex_[1]] + mesh.vertices_[mesh.cells_[i].vertex_[5]]) * 1.0 / 2 * 1.0;
		temp_vex[9] = newmesh.add_vertex(edge_mid_2);

		//第三条边的中点
		V3f edge_mid_3;
		edge_mid_3 = (mesh.vertices_[mesh.cells_[i].vertex_[5]] + mesh.vertices_[mesh.cells_[i].vertex_[4]]) * 1.0 / 2 * 1.0;
		temp_vex[10] = newmesh.add_vertex(edge_mid_3);

		//第四条边的中点
		V3f edge_mid_4;
		edge_mid_4 = (mesh.vertices_[mesh.cells_[i].vertex_[4]] + mesh.vertices_[mesh.cells_[i].vertex_[0]]) * 1.0 / 2 * 1.0;
		temp_vex[11] = newmesh.add_vertex(edge_mid_4);

		//第五条边的中点
		V3f edge_mid_5;
		edge_mid_5 = (mesh.vertices_[mesh.cells_[i].vertex_[3]] + mesh.vertices_[mesh.cells_[i].vertex_[2]]) * 1.0 / 2 * 1.0;
		temp_vex[12] = newmesh.add_vertex(edge_mid_5);

		//第六条边的中点
		V3f edge_mid_6;
		edge_mid_6 = (mesh.vertices_[mesh.cells_[i].vertex_[3]] + mesh.vertices_[mesh.cells_[i].vertex_[7]])*1.0 / 2 * 1.0;
		temp_vex[13] = newmesh.add_vertex(edge_mid_6);

		//第七条边的中点
		V3f edge_mid_7;
		edge_mid_7 = (mesh.vertices_[mesh.cells_[i].vertex_[6]] + mesh.vertices_[mesh.cells_[i].vertex_[7]])*1.0 / 2 * 1.0;
		temp_vex[14] = newmesh.add_vertex(edge_mid_7);

		//第八条边的中点
		V3f edge_mid_8;
		edge_mid_8 = (mesh.vertices_[mesh.cells_[i].vertex_[2]] + mesh.vertices_[mesh.cells_[i].vertex_[6]])*1.0 / 2 * 1.0;
		temp_vex[15] = newmesh.add_vertex(edge_mid_8);

		//第九条边的中点
		V3f edge_mid_9;
		edge_mid_9 = (mesh.vertices_[mesh.cells_[i].vertex_[1]] + mesh.vertices_[mesh.cells_[i].vertex_[2]])*1.0 / 2 * 1.0;
		temp_vex[16] = newmesh.add_vertex(edge_mid_9);

		//第十条边的中点
		V3f edge_mid_10;
		edge_mid_10 = (mesh.vertices_[mesh.cells_[i].vertex_[5]] + mesh.vertices_[mesh.cells_[i].vertex_[6]])*1.0 / 2 * 1.0;
		temp_vex[17] = newmesh.add_vertex(edge_mid_10);

		//第十一条边的中点
		V3f edge_mid_11;
		edge_mid_11 = (mesh.vertices_[mesh.cells_[i].vertex_[0]] + mesh.vertices_[mesh.cells_[i].vertex_[3]])*1.0 / 2 * 1.0;
		temp_vex[18] = newmesh.add_vertex(edge_mid_11);

		//第十二条边的中点
		V3f edge_mid_12;
		edge_mid_12 = (mesh.vertices_[mesh.cells_[i].vertex_[4]] + mesh.vertices_[mesh.cells_[i].vertex_[7]])*1.0 / 2 * 1.0;
		temp_vex[19] = newmesh.add_vertex(edge_mid_12);

		//再插入每个面的中点
		for (int mesh_face_num = 0; mesh_face_num < 6; mesh_face_num++)
		{
			V3f face_mid;
			face_mid = (mesh.vertices_[mesh.faces_[mesh.cells_[i].faces_[mesh_face_num]].vs_f[0]]
				+ mesh.vertices_[mesh.faces_[mesh.cells_[i].faces_[mesh_face_num]].vs_f[1]]
				+ mesh.vertices_[mesh.faces_[mesh.cells_[i].faces_[mesh_face_num]].vs_f[2]]
				+ mesh.vertices_[mesh.faces_[mesh.cells_[i].faces_[mesh_face_num]].vs_f[3]])*1.0 / 4 * 1.0;
			temp_vex[20 + mesh_face_num] = newmesh.add_vertex(face_mid);
		}

		//插入单元的内心
		V3f cell_mid;
		for (int j = 0; j < 8; j++)
		{
			cell_mid = cell_mid + mesh.vertices_[mesh.cells_[i].vertex_[j]];
		}
		cell_mid = cell_mid * 1.0 / 8 * 1.0;
		temp_vex[26] = newmesh.add_vertex(cell_mid);

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
		for (int j = 0; j < 8; j++) {
			newmesh.add_cell(sub_cell[j]);
		}
	}
	mesh = newmesh;
	return mesh;
}

/**************************************************
@brief   : mlca 核心部分  也就是简单的 A Subdivision Scheme for Hexahedral Meshes 2.2小节 公式1
		 : 感觉这个算法不太对呀，因为所有点都是0，这就很坑了
@author  : lee
@input   ：none
@output  ：none
@time    : none
**************************************************/
HexV3fMesh& MLCASubdivision::Averaging(HexV3fMesh &mesh)
{
	std::vector<V3f> changed_vex;
	changed_vex.resize(mesh.vertices_.size());
	for (uint64_t vex_num = 0; vex_num < mesh.vertices_.size(); vex_num++)
	{
		if (mesh.is_bdy(VertexHandle(vex_num)))// 暂时不知道isbdy 是什么鬼？QU  GU:猜测应该是否是边界点的意思  应该就是边界的意思
		{
			std::set<FaceHandle, compare_OVM> neighbor_vf = mesh.neighbor_v_f[vex_num]; // 点的邻接面
			std::vector<FaceHandle> bdy_face;
			std::set<FaceHandle, compare_OVM>::iterator vf_id;

			for (vf_id = neighbor_vf.begin(); vf_id != neighbor_vf.end(); vf_id++)
			{
				if (mesh.is_bdy(*vf_id))
				{
					bdy_face.push_back(*vf_id);
				}
			}
			int vf_num = bdy_face.size();
			std::vector<V3f> vf_cen(vf_num);
			std::vector<FaceHandle>::iterator bdy_face_id;
			int u = 0;
			// 边界面的中心点
			// 如果采用边界面的中心点，似乎确实好一些些
			for (bdy_face_id = bdy_face.begin(); bdy_face_id != bdy_face.end(); bdy_face_id++, u++)
			{
				std::vector<VertexHandle> vf = mesh.faces_[*bdy_face_id].vs_f;
				vf_cen[u] = vf_cen[u] + (mesh.vertices_[vf[0]] + mesh.vertices_[vf[1]]
					+ mesh.vertices_[vf[2]] + mesh.vertices_[vf[3]]) / 4 * 1.0;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex = new_vex + vf_cen[vf_cen_num];
			}
			new_vex = new_vex / vf_num * 1.0;
			changed_vex[vex_num] = new_vex;
		}
		else
		{
			std::set<CellHandle, compare_OVM> neighbor_vc = mesh.neighbor_v_c[vex_num];
			std::set<CellHandle, compare_OVM>::iterator vc_id;
			int vf_num = neighbor_vc.size();
			std::vector<V3f> vc_cen(vf_num);
			int u = 0;
			// 处在中心点一般有8个cell，然后计算每个cell的点的和 / 8 ,也就是 每个cell 的点的平均值
			for (vc_id = neighbor_vc.begin(); vc_id != neighbor_vc.end(); vc_id++, u++)
			{
				std::vector<VertexHandle> vf = mesh.cells_[*vc_id].vertex_;
				vc_cen[u] = (mesh.vertices_[vf[0]] + mesh.vertices_[vf[1]]
					+ mesh.vertices_[vf[2]] + mesh.vertices_[vf[3]] +
					mesh.vertices_[vf[4]] + mesh.vertices_[vf[5]] +
					mesh.vertices_[vf[6]] + mesh.vertices_[vf[7]]) / 8.0;
			}
			V3f new_vex;
			for (int vf_cen_num = 0; vf_cen_num < vf_num; vf_cen_num++)
			{
				new_vex = new_vex + vc_cen[vf_cen_num];
			}
			new_vex = new_vex /  vf_num * 1.0; // 对于新顶点就是 8 个cell 的 重心坐标 的和 / cell的个数
			changed_vex[vex_num] = new_vex;
		}
	}
	for (uint64_t vex_num = 0; vex_num < mesh.vertices_.size(); vex_num++)
	{
		mesh.vertices_[vex_num] = changed_vex[vex_num];
	}
	return mesh;
}



/**************************************************
@brief   : MLCA细分PipeLine
@author  : lee
@input   ：none
@output  ：none
@time    : none
**************************************************/
//HexV3fMesh& MLCASubdivision::Subdivising() {
//	//*mesh = Spliting();
//	//Averaging();
//	return NULL;
//}

}