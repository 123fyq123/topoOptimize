#include "AABBTree.h"
#include<queue>
#include<algorithm>
#include"head.h"
AABBTree::AABBTree(const std::vector<std::vector<CPoint>>& V, const std::vector<std::vector<CPoint>>& P_N, std::vector<int>& blockIndices)
	:m_V(V),m_P(P_N)
{
	m_EachNode.resize(blockIndices.size());
	CreateEachNodeBox(V, blockIndices);
	m_Root = BuildNode(blockIndices);   // 建树
}

AABBTree::AABBTree(const std::vector<std::vector<CPoint>>& V, std::vector<int>& blockIndices)
{
	m_EachNode.resize(blockIndices.size());
	CreateEachNodeBox(V, blockIndices);
	m_Root = BuildNode(blockIndices);   // 建树
}


std::vector<pair<int, int>> AABBTree::blockFindCollisionCandidates()
{
	std::vector<pair<int, int>> candidates;

	std::queue<std::shared_ptr<AABBTreeNode>> searchQueue;

	for (int i = 0; i < m_EachNode.size(); ++i)
	{
		searchQueue.push(m_Root);

		while (!searchQueue.empty()) {
			auto cur = searchQueue.front();
			searchQueue.pop();
			if (cur->is_Leaf) {
				if (cur->ind > i) {
					candidates.push_back({ i, cur->ind });
				}
				continue;
			}
			if (IsCollision(m_EachNode[i], cur->left->box))
				searchQueue.push(cur->left);
			if (IsCollision(m_EachNode[i], cur->right->box))
				searchQueue.push(cur->right);
		}
	}

	return candidates;
}

std::vector<pair<int, int>> AABBTree::refineFindCollisionCandidates(const std::vector<std::vector<CPoint>>& V)
{
	std::vector<pair<int, int>> candidates;
	
	std::queue<std::shared_ptr<AABBTreeNode>> searchQueue;
	
	for (int i = 0; i < m_EachNode.size(); ++i)
	{
		searchQueue.push(m_Root);
	
		while (!searchQueue.empty()) {
			auto cur = searchQueue.front();
			searchQueue.pop();
			if (cur->is_Leaf) {
				if (cur->ind > i && !IsConnect(i,cur->ind,V)) {
					candidates.push_back({ i, cur->ind });
				}
				continue;
			}
			if (IsCollision(m_EachNode[i], cur->left->box))
				searchQueue.push(cur->left);
			if (IsCollision(m_EachNode[i], cur->right->box))
				searchQueue.push(cur->right);
		}
	}
	
	return candidates;
}

bool AABBTree::refineCollideDetect(const AABBTree& other,  ContactPoint& collid1, ContactPoint& collid2)
{
	std::queue<std::shared_ptr<AABBTreeNode>> searchQueue;
	bool isCollide = false;
	for (int i = 0; i < m_EachNode.size(); ++i)
	{
		searchQueue.push(other.m_Root);

		while (!searchQueue.empty()) {
			auto cur = searchQueue.front();
			searchQueue.pop();
			if (cur->is_Leaf) {
				int pi = 0, pj = 0;
				double distTemp = m_V[i][pi].getDistance(other.m_V[cur->ind][pj]);
				for (int piIndex = 0; piIndex < m_V[i].size()-1; ++piIndex) {   // 计算两个包围盒的距离
					for (int pjIndex = 0; pjIndex < other.m_V[cur->ind].size()-1; ++pjIndex) {
						double tempTemp = m_V[i][piIndex].getDistance(other.m_V[cur->ind][pjIndex]);
						if ( tempTemp< distTemp) {
							pi = piIndex;
							pj = pjIndex;
							distTemp = tempTemp;
						}
					}
				}
				if (distTemp > collid1.dist) {   // 找到最大距离的两个点 记录起点 终点 距离 起点的参数坐标
					collid1.src = m_V[i][pi]; collid1.des = other.m_V[cur->ind][pj];
					collid1.dist = distTemp; collid1.parm = m_P[i][pi];
				}
				if (distTemp > collid2.dist) {  //  找到最大距离的两个点 记录起点 终点 距离 起点的参数坐标
					collid2.des = m_V[i][pi]; collid2.src = other.m_V[cur->ind][pj];
					collid2.dist = distTemp; collid2.parm = other.m_P[cur->ind][pj];
				}
				isCollide = true;
				continue;
			}
			if (IsCollision(m_EachNode[i], cur->left->box))
				searchQueue.push(cur->left);
			if (IsCollision(m_EachNode[i], cur->right->box))
				searchQueue.push(cur->right);
		}
	}
	return isCollide;
}

void AABBTree::drawAABBTree(const string filename)
{
	int count = 0;

	std::ofstream myfile;
	myfile.open("D:/Resource/FileDataOutput/AABB/" + filename + ".pvd");
	myfile << "<?xml version = \"1.0\"?>\n"
		<< "<VTKFile type = \"Collection\" version = \"0.1\"><Collection>\n";


	for (int i = 0; i < m_EachNode.size(); ++i)
	{
		std::shared_ptr<AABBTreeNode> it = std::make_shared<AABBTreeNode>(m_EachNode[i], nullptr, nullptr);
		printFile("D:/Resource/FileDataOutput/AABB/" + filename + "_" + to_string(count) + ".vts", it);
		myfile << "<DataSet part = \"" + to_string(count) + "\" file = \"" + filename + "_" + to_string(count) + ".vts" + "\"/>\n";
		count++;
	}

	//queue<std::shared_ptr<AABBTreeNode>> node;
	//node.push(m_Root);
	//while (!node.empty())
	//{
	//	auto it = node.front();
	//	node.pop();
	//	if (it->is_Leaf) {
	//		printFile(filename + "_" + to_string(count) + ".vts", it);
	//		myfile << "<DataSet part = \""+to_string(count)+"\" file = \"" + filename + "_" + to_string(count) + ".vts" + "\"/>\n";
	//		count++;
	//	}
	//	if (it->left) node.push(it->left);
	//	if (it->right) node.push(it->right);
	//}
	myfile << "</Collection>\n"
		<< "</VTKFile>\n";
	myfile.close();
}

void AABBTree::printFile(std::string finename, const std::shared_ptr<AABBTreeNode>& node)
{
	double lx = node->box.lowerBound[0], ly = node->box.lowerBound[1], lz = node->box.lowerBound[2];
	double ux = node->box.upperBound[0], uy = node->box.upperBound[1], uz = node->box.upperBound[2];
	std::ofstream myfile;
	myfile.open(finename);
	myfile << "<?xml version=\"1.0\"?>\n"
		<< "<VTKFile type = \"StructuredGrid\" version = \"0.1\">\n"
		<< "<StructuredGrid WholeExtent = \"0 1 0 1 0 1\">\n"
		<< "<Piece Extent = \"0 1 0 1 0 1\">\n"
		<< "<Points>\n"
		<< "<DataArray type = \"Float32\" NumberOfComponents = \"3\">\n"
		<< to_string(lx) + " " + to_string(ly) + " " + to_string(lz) + " " + to_string(ux) + " " + to_string(ly) + " " + to_string(lz) + "\n"
		<< to_string(lx) + " " + to_string(uy) + " " + to_string(lz) + " " + to_string(ux) + " " + to_string(uy) + " " + to_string(lz) + "\n"
		<< to_string(lx) + " " + to_string(ly) + " " + to_string(uz) + " " + to_string(ux) + " " + to_string(ly) + " " + to_string(uz) + "\n"
		<< to_string(lx) + " " + to_string(uy) + " " + to_string(uz) + " " + to_string(ux) + " " + to_string(uy) + " " + to_string(uz) + "\n"
		<< "</DataArray>\n"
		<< "</Points>\n"
		<< "</Piece>\n"
		<< "</StructuredGrid>\n"
		<< "</VTKFile>\n";
	myfile.close();
}

std::shared_ptr<AABBTreeNode> AABBTree::BuildNode(std::vector<int>& blockIndices)
{
	AABB box = FindBoundary(blockIndices);
	std::shared_ptr<AABBTreeNode> root = std::make_shared<AABBTreeNode>(box, nullptr, nullptr);
	if (blockIndices.size() == 1) {
		root->is_Leaf = true;
		root->ind = blockIndices.back();
		return root;
	}

	double x_Dist = box.upperBound[0] - box.lowerBound[0];
	double y_Dist = box.upperBound[1] - box.lowerBound[1];
	double z_Dist = box.upperBound[2] - box.lowerBound[2];
	int maxIndex = x_Dist > y_Dist ? (x_Dist > z_Dist ? 0 : 2) : (y_Dist > z_Dist ? 1 : 2);
	auto& node = m_EachNode;
	std::sort(blockIndices.begin(), blockIndices.end(), [&maxIndex, &node](const int a, const int b) {
		return node[a].lowerBound[maxIndex] < node[b].lowerBound[maxIndex];
		});
	int half = blockIndices.size() / 2;
	std::vector<int> firstInd, secondInd;
	for (int i = 0; i < half; ++i) firstInd.emplace_back(blockIndices[i]);
	for (int i = half; i < blockIndices.size(); ++i) secondInd.emplace_back(blockIndices[i]);
	root->left = BuildNode(firstInd);
	root->right = BuildNode(secondInd);
	return root;
}

void AABBTree::CreateEachNodeBox(const std::vector<std::vector<CPoint>>& V, const std::vector<int>& blockIndices)
{
	for (size_t i = 0; i < V.size(); ++i)
		m_EachNode[blockIndices[i]] = BoxInArrayPoints(V[i]);
}

AABB AABBTree::FindBoundary(const std::vector<int>& blockIndices)
{
	if (blockIndices.size() < 1)
	{
		std::cout << " Error bounding\n" << std::endl;
	}
	if (blockIndices.size() == 1) return m_EachNode[blockIndices[0]];
	double minX = m_EachNode[blockIndices[0]].lowerBound[0], minY = m_EachNode[blockIndices[0]].lowerBound[1],
		minZ = m_EachNode[blockIndices[0]].lowerBound[2];
	double maxX = m_EachNode[blockIndices[0]].upperBound[0], maxY = m_EachNode[blockIndices[0]].upperBound[1],
		maxZ = m_EachNode[blockIndices[0]].upperBound[2];

	for (size_t i = 0; i < blockIndices.size(); ++i)
	{
		int ind = blockIndices[i];
		minX = min(m_EachNode[ind].lowerBound[0], minX); maxX = max(m_EachNode[ind].upperBound[0], maxX);
		minY = min(m_EachNode[ind].lowerBound[1], minY); maxY = max(m_EachNode[ind].upperBound[1], maxY);
		minZ = min(m_EachNode[ind].lowerBound[2], minZ); maxZ = max(m_EachNode[ind].upperBound[2], maxZ);
	}
	AABB box{
		minX, minY,minZ,
		maxX,maxY,maxZ
	};
	return box;
}

AABB AABBTree::BoxInArrayPoints(const std::vector<CPoint>& Points) const
{
	if (Points.size() < 1) {
		std::cout << "empty points\n" << std::endl;
		exit(1);
	}
	double minX = Points[0].getX(), minY = Points[0].getY(), minZ = Points[0].getZ();
	double maxX = Points[0].getX(), maxY = Points[0].getY(), maxZ = Points[0].getZ();
	for (size_t i = 0; i < Points.size(); ++i)
	{
		minX = min(Points[i].getX(), minX); maxX = max(Points[i].getX(), maxX);
		minY = min(Points[i].getY(), minY); maxY = max(Points[i].getY(), maxY);
		minZ = min(Points[i].getZ(), minZ); maxZ = max(Points[i].getZ(), maxZ);
	}
	AABB box{
		minX,minY,minZ,
		maxX,maxY,maxZ
	};
	return box;
}

bool AABBTree::IsConnect(int i, int j, const std::vector<std::vector<CPoint>>& V)
{
	for (int m = 0; m < V[i].size(); ++m) {
		LocInd temp1; temp1.x = V[i][m].getX();
		temp1.y = V[i][m].getY(); temp1.z = V[i][m].getZ();
		for (int n = 0; n < V[j].size(); ++n) {
			LocInd temp2; temp2.x = V[j][n].getX();
			temp2.y = V[j][n].getY(); temp2.z = V[j][n].getZ();
			if (temp1 == temp2) return true;
		}
	}
	return false;
}

bool AABBTree::IsCollision(const AABB& block1, const AABB& block2)
{
	double xOver = calOverOne(block1.upperBound[0], block1.lowerBound[0], block2.upperBound[0], block2.lowerBound[0]);
	double yOver = calOverOne(block1.upperBound[1], block1.lowerBound[1], block2.upperBound[1], block2.lowerBound[1]);
	double zOver = calOverOne(block1.upperBound[2], block1.lowerBound[2], block2.upperBound[2], block2.lowerBound[2]);

	return xOver > 0 && yOver > 0 && zOver > 0;
}

double AABBTree::calOverOne(float amax, float amin, float bmax, float bmin)
{
	float minMax = 0;
	float maxMin = 0;
	if (amax < bmax)//a物体在b物体左侧
	{
		minMax = amax;
		maxMin = bmin;
	}
	else //a物体在b物体右侧
	{
		minMax = bmax;
		maxMin = amin;
	}

	if (minMax > maxMin) {
		return minMax - maxMin;
	}
	else {
		return 0;
	}
}


