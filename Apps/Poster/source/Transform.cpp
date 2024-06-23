#include"Transform.h"


//输出三角化模型


void Transform::Quad2Triangle() {

	triangleMesh tm;
	std::ofstream out(OutputPath);
	int bdyv = 0;
	int bdyf = 0;
	out << "OFF" << std::endl;
	HexModel->examine_bdy();
	std::unordered_map<int, int> vv;
	for (int i = 0; i < HexModel->vertices_.size(); ++i) {
		if (HexModel->is_bdy(VertexHandle(i))) {
			//out << "v " << HexModel->vertices_[i].x << " " << HexModel->vertices_[i].y << " "<< HexModel->vertices_[i].z << std::endl;
			vv[i] = bdyv;
			bdyv++;
		}
	}
	for (int i = 0; i < HexModel->faces_.size(); ++i) {
		if (HexModel->is_bdy(FaceHandle(i))) {
			bdyf++;
			bdyf++;
		}
	}
	out << bdyv << " " << bdyf << " 0" << std::endl;
	bdyv = 0;
	for (int i = 0; i < HexModel->vertices_.size(); ++i) {
		if (HexModel->is_bdy(VertexHandle(i))) {
			tm.vertices.push_back(HexModel->vertices_[i]);
		}
	}
	std::map<std::vector<int>,int> eID;
	for (int i = 0; i < HexModel->faces_.size(); ++i) {
		if (HexModel->is_bdy(FaceHandle(i))) {
			std::vector<int> f1{ vv[HexModel->faces_[i].vs_f[0]] ,vv[HexModel->faces_[i].vs_f[1]] ,vv[HexModel->faces_[i].vs_f[2]] };
			std::vector<int> f2{ vv[HexModel->faces_[i].vs_f[0]] ,vv[HexModel->faces_[i].vs_f[2]] ,vv[HexModel->faces_[i].vs_f[3]] };
			tm.faces.push_back(f1);
			for (int kk = 0; kk < 3; ++kk) {
				std::vector<int> e{ std::min(f1[kk], f1[(kk + 1) % 3]), std::max(f1[kk], f1[(kk + 1) % 3]) };
				if (eID[e] == 0){
					eID[e] = eID.size();
				}
				tm.e[eID[e]] = e;
				tm.fe[tm.faces.size() - 1].push_back(eID[e]);
				tm.ef[eID[e]].push_back(tm.faces.size() - 1);
			}
			tm.faces.push_back(f2);
			for (int kk = 0; kk < 3; ++kk) {
				std::vector<int> e{ std::min(f2[kk], f2[(kk + 1) % 3]), std::max(f2[kk], f2[(kk + 1) % 3]) };
				if (eID[e] == 0) {
					eID[e] = eID.size();
				}
				tm.e[eID[e]] = e;
				tm.fe[tm.faces.size() - 1].push_back(eID[e]);
				tm.ef[eID[e]].push_back(tm.faces.size() - 1);
			}

		}
	}
	tm.ReorderMesh();
	//std::cout << "tm.faces.size() = " << tm.faces.size() << std::endl;
	for (int i = 0; i < HexModel->vertices_.size(); ++i) {
		if (HexModel->is_bdy(VertexHandle(i))) {
			out << HexModel->vertices_[i].x << " " << HexModel->vertices_[i].y << " "<< HexModel->vertices_[i].z << std::endl;
		}
	}
	for (int i = 0; i < tm.faces.size(); i += 2) {
		/*out << "f" << " " << tm.faces[i + 1][2] + 1 << " " << tm.faces[i + 1][1] + 1 << " " << tm.faces[i + 1][0] + 1 << std::endl;
		out << "f" << " " << tm.faces[i][2]+1 << " " << tm.faces[i][1]+1 << " " << tm.faces[i][0]+1 ;*/
		out << 3 << " " << tm.faces[i + 1][2]  << " " << tm.faces[i + 1][1]  << " " << tm.faces[i + 1][0]  << std::endl;
		out << 3 << " " << tm.faces[i][2] << " " << tm.faces[i][1]  << " " << tm.faces[i][0] ;
		out << std::endl;
	}
	out.close();
}


void Transform::Hex2IGEs() {
	;
}