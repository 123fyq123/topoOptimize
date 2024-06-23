#pragma once
#include<vector>
#include"Geom/CPoint.h"
#include<memory>
#include<unordered_set>
#include<set>
struct AABB
{
	double lowerBound[3];  // 左下
	double upperBound[3];  // 右上
};

struct AABBTreeNode
{
	AABB box;
	std::shared_ptr<AABBTreeNode> left;
	std::shared_ptr<AABBTreeNode> right;
	bool is_Leaf;        // 是否叶子节点
	int ind;             // 包围盒索引
	AABBTreeNode(AABB _box, std::shared_ptr<AABBTreeNode> _left, std::shared_ptr<AABBTreeNode> _right)
		:box(_box), left(_left), right(_right), is_Leaf(false), ind(-1)
	{}
};

struct ContactPoint {
	CPoint src;
	CPoint des;
	double dist = 0.0;
	CPoint parm;
};

class AABBTree {
	std::vector<AABB> m_EachNode;
	std::vector<std::vector<CPoint>> m_V;
	std::vector<std::vector<CPoint>> m_P;
	std::shared_ptr<AABBTreeNode> m_Root;
public:
	AABBTree(const std::vector<std::vector<CPoint>>& V, const std::vector<std::vector<CPoint>>& P_N, std::vector<int>& blockIndices);
	AABBTree(const std::vector<std::vector<CPoint>>& V,  std::vector<int>& blockIndices);

	AABBTree():m_Root(nullptr){}
	//param: 点集
	//		 
	//碰撞候选
	std::vector<pair<int, int>> blockFindCollisionCandidates();

	std::vector<pair<int, int>> refineFindCollisionCandidates(const std::vector<std::vector<CPoint>>& V);

	bool refineCollideDetect(const AABBTree& other,ContactPoint& collid1, ContactPoint& collid2);


	void drawAABBTree(const string filename);

	void printFile(std::string finename, const std::shared_ptr<AABBTreeNode>& node);

	int getLength() { return (int)m_EachNode.size(); }
	//private:
		// 递归构造AABB树
	std::shared_ptr<AABBTreeNode> BuildNode(std::vector<int>& blockIndices);

	//构造所有块的包围盒
	void CreateEachNodeBox(const std::vector<std::vector<CPoint>>& V, const std::vector<int>& blockIndices);

	// 控制块的包围盒
	AABB FindBoundary(const std::vector<int>& blockIndices);

	//一组控制点下的包围盒
	AABB BoxInArrayPoints(const std::vector<CPoint>& Points) const;

	bool IsConnect(int i, int j, const std::vector<std::vector<CPoint>>& V);

	//碰撞
	bool IsCollision(const AABB& block1, const AABB& block2);
	double calOverOne(float amax, float amin, float bmax, float bmin);
};
