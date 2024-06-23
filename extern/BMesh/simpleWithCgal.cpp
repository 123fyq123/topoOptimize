#include "simpleWithCgal.h"
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Mean_curvature_flow_skeletonization.h>
#include <CGAL/extract_mean_curvature_flow_skeleton.h>
#include <CGAL/boost/graph/split_graph_into_polylines.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include "Common/GeometryKernel.h"
#include <fstream>

typedef CGAL::Simple_cartesian<double>                        Kernel;
typedef Kernel::Point_3                                       Point;
typedef CGAL::Polyhedron_3<Kernel>                            Polyhedron;

typedef boost::graph_traits<Polyhedron>::vertex_descriptor    vertex_descriptor;

typedef CGAL::Mean_curvature_flow_skeletonization<Polyhedron> Skeletonization;
typedef Skeletonization::Skeleton                             SkeletonL;

typedef SkeletonL::vertex_descriptor                           Skeleton_vertex;
typedef SkeletonL::edge_descriptor                             Skeleton_edge;


struct Display_polylines {
	const Skeletonization::Skeleton& skeleton;
	std::ofstream& out;
	int polyline_size;
	std::stringstream sstr;
	Display_polylines(const Skeletonization::Skeleton& skeleton, std::ofstream& out)
		: skeleton(skeleton), out(out)
	{}
	void start_new_polyline() {
		polyline_size = 0;
		sstr.str("");
		sstr.clear();
	}
	void add_node(Skeleton_vertex v) {
		++polyline_size;
		sstr << " " << skeleton[v].point;
	}
	void end_polyline()
	{
		out << polyline_size << sstr.str() << "\n";
	}
};

int CgalMcfskel::Run() {
	std::ifstream input(NowProcessPath);
	Polyhedron tmesh;
	input >> tmesh;
	if (!CGAL::is_triangle_mesh(tmesh)) {
		std::cerr << "[ERROR] Input geometry is not triangulated." << std::endl;
		return EXIT_FAILURE;
	}
	SkeletonL skeleton;
	Skeletonization mcs(tmesh);

	// 1. Contract the mesh by mean curvature flow.
	mcs.contract_geometry();
	// 2. Collapse short edges and split bad triangles.
	mcs.collapse_edges();
	mcs.split_faces();
	// 3. Fix degenerate vertices.
	mcs.detect_degeneracies();
	// Perform the above three steps in one iteration.
	mcs.contract();
	// Iteratively apply step 1 to 3 until convergence.
	mcs.contract_until_convergence();
	// Convert the contracted mesh into a curve skeleton and
	// get the correspondent surface points
	mcs.convert_to_skeleton(skeleton);

	//std::cout << "[DEBUG] Number of vertices of the skeleton: " << boost::num_vertices(skeleton) << "\n";
	//std::cout << "[DEBUG] Number of edges of the skeleton: " << boost::num_edges(skeleton) << "\n";
	/*Test*/
	//std::ofstream output("skel-poly.polylines.txt");
	//Display_polylines display(skeleton, output);
	//CGAL::split_graph_into_polylines(skeleton, display);
	//output.close();
	//// Output skeleton points and the corresponding surface points
	//output.open("correspondance-poly.polylines.txt");
	//for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton)))
	//	for (vertex_descriptor vd : skeleton[v].vertices)
	//		output << "2 " << skeleton[v].point << " "
	//		<< get(CGAL::vertex_point, tmesh, vd) << "\n";
	/*Test*/
	std::vector<Point> SkelPoint;
	std::vector<E2I> SkelEdge;
	std::vector<Point>::iterator it1, it2;//寻找是否存过
	int nVertices = 0, nEdges = 0;
	int left, right;
	for (Skeleton_edge e : CGAL::make_range(edges(skeleton))) {
		const Point& s = skeleton[source(e, skeleton)].point;
		const Point& t = skeleton[target(e, skeleton)].point;
		it1 = find(SkelPoint.begin(), SkelPoint.end(), s);
		if (it1 == SkelPoint.end()) {
			SkelPoint.push_back(s);
			left = nVertices;
			nVertices++;
		}
		else {
			left = it1 - SkelPoint.begin();
		}
		it2 = find(SkelPoint.begin(), SkelPoint.end(), t);
		if (it2 == SkelPoint.end()) {
			SkelPoint.push_back(t);
			right = nVertices;
			nVertices++;
		}
		else {
			right = it2 - SkelPoint.begin();
		}
		E2I temp(left + 1, right + 1);
		SkelEdge.push_back(temp);
		nEdges++;
	}
	//骨架点邻接
	std::vector<std::vector<int>> point_nei(nVertices);
	for (auto iter : SkelEdge) {
		point_nei[iter.begin - 1].push_back(iter.end - 1);
		point_nei[iter.end - 1].push_back(iter.begin - 1);
	}
	//以上存好了点、边，现在单独处理，把对应半径提取出来
	std::vector<double> SkelRadius(nVertices, 0.0);
	for (Skeleton_vertex v : CGAL::make_range(vertices(skeleton))) {
		int count = 0;
		double TotalRadius = 0.0;
		double minRadius = 1e10;
		for (vertex_descriptor vd : skeleton[v].vertices) { //vertices
			auto x = sqrt((get(CGAL::vertex_point, tmesh, vd) - skeleton[v].point).squared_length());
			if (x < minRadius) {
				minRadius = x;
			}
			count++;
		}
		double resultRadius = 0.0;
		if (0 == count) {
			resultRadius = -1.0;
		}
		else {
			resultRadius = minRadius * 2 / sqrt(3);
		}
		auto y = skeleton[v].point;
		std::vector<Point>::iterator it2;//找骨架点在vector里的下标
		it2 = find(SkelPoint.begin(), SkelPoint.end(), y);
		int idx = it2 - SkelPoint.begin();
		SkelRadius[idx] = resultRadius;
	}
	// some radius of point have not update, then let his ring radius of point to update it;
	for (int i = 0; i < nVertices; i++) {
		if (SkelRadius[i] != -1.0)	continue;
		if (SkelRadius[i] == -1.0) {
			for (auto iter : point_nei[i]) {
				if (SkelRadius[iter] == -1.0)	continue;
				else {
					SkelRadius[i] = SkelRadius[iter];
				}
			}
		}
	}
	OurSkel newskel;
	newskel.SkelAll.resize(nVertices);
	double avgRadius = 0;
	for (int i = 0; i < nVertices; i++) {
		newskel.SkelAll[i].index = i;
		newskel.SkelAll[i].pos = V3f(SkelPoint[i][0], SkelPoint[i][1], SkelPoint[i][2]);
		newskel.SkelAll[i].radius = SkelRadius[i] * 0.1;
		avgRadius += newskel.SkelAll[i].radius;
		newskel.SkelAll[i].is_valid = true;
		newskel.SkelAll[i].Indegree = point_nei[i].size();
		for (auto iter : point_nei[i])
			if (iter != 0)
				newskel.SkelAll[i].neighbor.push_back(iter);
	}
	avgRadius /= newskel.SkelAll.size();
	SkelResult = newskel;
	// 避免生成的骨架半径相交，对点集做预处理


	//开始输出skel文件
	std::ofstream output2(outProcessPath);
	output2 << "ID Cx Cy Cz RADIUS #NEIGHBORS NEIGHBORS_LIST\n";
	output2 << nVertices << "\n";
	for (int i = 0; i < nVertices; ++i) {
		output2 << i << " " << newskel.SkelAll[i].pos.x << " " << newskel.SkelAll[i].pos.y << " " << newskel.SkelAll[i].pos.z << " "
			<< avgRadius << " " << point_nei[i].size();
		for (auto iter : point_nei[i]) {
			output2 << " " << iter;
		}
		output2 << "\n";
	}
	output2.close();

	return EXIT_SUCCESS;
}