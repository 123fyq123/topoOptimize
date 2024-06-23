#include <queue>
#include "edge_subdivision.h"
using namespace std;

void EdgeSubdivision::Spliting(MeshKernel::VolumeMesh* hexMesh, int eid, int split_num) {
	vector< vector<VH> > newCellsVhs;

	unordered_map<EH, vector<VH>> deleteEdges_v;	// ±ßµã
	//vector<bool> edge_visit(hexMesh->EdgeSize(), false);	// queue visit
	vector<bool> cell_visit(hexMesh->CellSize(), false);	// solve cell finish is true
	
	auto lineSpliting = [&](EH eh, VH kown_first_vh) {
		vector<VH> one_edge;
		VH first_vh, second_vh;
		if (kown_first_vh == hexMesh->edges(eh).vh1()) {
			cout << "yes" << endl;
			first_vh = hexMesh->edges(eh).vh1();
			second_vh = hexMesh->edges(eh).vh2();
		}
		else {
			cout << "no" << endl;
			first_vh = hexMesh->edges(eh).vh2();
			second_vh = hexMesh->edges(eh).vh1();
		}
		Vex& first = hexMesh->vertices(first_vh);
		Vex& second = hexMesh->vertices(second_vh);
		one_edge.push_back(first_vh);
		for (int i = 1; i < split_num; i++) {
			Vex nowV = first + (second - first) * i / split_num;
			one_edge.push_back(hexMesh->AddVertex(nowV));
		}
		one_edge.push_back(second_vh);
		// cout << one_edge.size() << endl;
		assert(one_edge.size() == split_num+1);
		deleteEdges_v[eh] = one_edge;
	};

	queue<EH> bfs_queue;
	bfs_queue.push((EH)eid);
	lineSpliting((EH)eid, hexMesh->edges((EH)eid).vh1());

	while (!bfs_queue.empty()) {
		EH eh_first = bfs_queue.front();
		bfs_queue.pop();
		for (auto ch : hexMesh->NeighborCh(eh_first)) {
			if (!cell_visit[ch.idx()]) {
				std::cout << ch.idx() << std::endl;
				// find two face
				FH first_fh, second_fh;
				for (auto fh : hexMesh->cells(ch).getFaceHandle()) {
					if (hexMesh->NeighborFh(eh_first).count(fh)) continue;
					if (hexMesh->NeighborFh(deleteEdges_v[eh_first][0]).count(fh)) first_fh = fh;
					else if (hexMesh->NeighborFh(deleteEdges_v[eh_first][split_num]).count(fh)) second_fh = fh;
				}
				MeshKernel::iGameFace& first_face = hexMesh->faces(first_fh);
				MeshKernel::iGameFace& second_face = hexMesh->faces(second_fh);
				// perprocess
				unordered_set<EH> cell_edge_set;
				for (auto eh : hexMesh->cells(ch).getEdgeHandle())
					cell_edge_set.insert(eh);
				unordered_set<VH> face_vertex_set;
				for (auto vh : second_face.getVertexHandle())
					face_vertex_set.insert(vh);
				// find four line
				vector<EH> four_edge;
				for (auto vh : first_face.getVertexHandle()) {
					for (auto eh : hexMesh->NeighborEh(vh)) {
						if (cell_edge_set.count(eh) == 0) continue;
						VH to_vh = hexMesh->edges(eh).vh1() == vh ? hexMesh->edges(eh).vh2() : hexMesh->edges(eh).vh1();
						if (face_vertex_set.count(to_vh) == 0) continue;
						// find one edge
						four_edge.push_back(eh);
						if (deleteEdges_v.count(eh) == 0) {
							bfs_queue.push(eh);
							lineSpliting(eh, vh);
						}
						break;
					}
				}
				//std::cout << "edge size: " << four_edge.size() << std::endl;
				assert(four_edge.size() == 4);
				// solve cell
				std::vector<VH> one_cell;
				for (int i = 0; i <= split_num; i++) {
					for (auto eh : four_edge)
						one_cell.push_back(deleteEdges_v[eh][i]);
				}
				newCellsVhs.push_back(one_cell);
				cell_visit[ch.idx()] = true;
			}
		}
	}

	for (auto oneCell : newCellsVhs) {
		assert(oneCell.size() % 4 == 0);
		for (int i = 1; i < oneCell.size() / 4; i++) {
			std::vector<VH> new_cell;
			for (int j = 0; j < 4; j++)
				new_cell.push_back(oneCell[(i - 1) * 4 + j]);
			for (int j = 0; j < 4; j++)
				new_cell.push_back(oneCell[i * 4 + j]);
			hexMesh->AddCell(new_cell);
		}
	}
	for (auto eh : deleteEdges_v)
		hexMesh->DeleteEdge(eh.first);

	hexMesh->updateAllHandles();
	std::cout << "vsize: " << hexMesh->vsize() << std::endl;
	std::cout << "fsize: " << hexMesh->fsize() << std::endl;
}