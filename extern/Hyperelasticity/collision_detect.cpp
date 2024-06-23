#include "collision_detect.h"
#include<set>
#include<unordered_map>


CollisionDetect::~CollisionDetect()
{
}

bool CollisionDetect::ObjectDetect(std::vector<CBSplineSolid>& geo, std::vector<std::vector<CollideInfo>>& collidInfo, std::vector<pair<int, int>> &patchIds)
{
	bool res = false;
	std::vector<std::vector<CPoint>> V;
	 // patch  boundary
	for (int i = 0; i < geo.size(); ++i) {
		boundCofes(geo[i], V, i, patchIds);
	}
	collidInfo.resize(patchIds.size());
	auto collisonBlock = findCollisions(V);
	std::set<int> isDraw;
	for (int i = 0; i < collisonBlock.size(); ++i)
	{
		int p = collisonBlock[i].first, q = collisonBlock[i].second;
		if (patchIds[p].first != patchIds[q].first && !m_ConnectPatch[patchIds[p].first].count(patchIds[q].first))
		{
			if (m_Tree.find(p) == m_Tree.end())
			{
				std::vector<std::vector<CPoint>> V_N, P_N;  // 物理点和参数点
				uniformGetPoints(geo[patchIds[p].first], V_N, P_N, 1000, patchIds[p].second);
				m_Tree[p] = CreateTree(V_N,P_N);
			}
			if (m_Tree.find(q) == m_Tree.end())
			{
				std::vector<std::vector<CPoint>> V_N, P_N; // 物理点和参数点          
				uniformGetPoints(geo[patchIds[q].first], V_N, P_N,1000, patchIds[q].second);
				m_Tree[q] = CreateTree(V_N, P_N);
			}
			ContactPoint collid1, collid2;
			
			bool isCollide = m_Tree[p].refineCollideDetect(m_Tree[q], collid1, collid2);
			//if (isCollide) {
			//	if(collidInfo[p].size() == 0) m_Tree[p].drawAABBTree("block_" + to_string(p) + "_box");
			//	if(collidInfo[q].size() == 0) m_Tree[q].drawAABBTree("block_" + to_string(q) + "_box");
			//}
			 
			if (isCollide) {
				collidInfo[p].push_back({ collid1.parm, collid1.des - collid1.src });
				collidInfo[q].push_back({ collid2.parm, collid2.des - collid2.src });
			}
			//if (collid2.size() > 0 && !isDraw.count(q)) {
			//	//m_Tree[q].drawAABBTree("block_" + to_string(q) + "_box");
			//	std::cout << q << std::endl;
			//}
			res |= isCollide;
		}
	}
	return res;
}

std::vector<std::pair<int, int>> CollisionDetect::findCollisions(const std::vector<std::vector<CPoint>>& V)
{
    std::vector<int> indices(V.size());
    for (int i = 0; i < V.size(); ++i) indices[i] = i;
    AABBTree root(V, indices);
    //std::string filename = "beam_cube";
    //root.drawAABBTree(filename);
    return root.blockFindCollisionCandidates();
}
bool CollisionDetect::refineFindCollisions(const std::vector<std::vector<CPoint>>& V)
{
    std::vector<int> indices(V.size());
    for (int i = 0; i < V.size(); ++i) indices[i] = i;
    AABBTree root(V, indices);
    std::string filename = "beam_refine_cube";
    //root.drawAABBTree(filename);
    return root.refineFindCollisionCandidates(V).size()  > 0;
}
void CollisionDetect::boundCofes(const CBSplineSolid& geo, std::vector<std::vector<CPoint>>& V, int patchID, std::vector<pair<int,int>>& patchIds)
{
	int usize = geo.getSizeUControlPoints();
	int vsize = geo.getSizeVControlPoints();
	int wsize = geo.getSizeWControlPoints();
	auto coefs = geo.getControlPoints();
	std::vector<std::vector<CPoint>> VTemp(6);
	for (int i = 0; i < usize; ++i)
	{
		for (int j = 0; j < vsize; ++j)
		{
			for (int w = 0; w < wsize; ++w)
			{
				if (i == 0 && m_Inter[patchID][0].first == -1) VTemp[0].emplace_back(coefs[i][j][w]);
				if (i == usize - 1 && m_Inter[patchID][1].first == -1) VTemp[1].emplace_back(coefs[i][j][w]);
				if (j == 0 && m_Inter[patchID][2].first == -1) VTemp[2].emplace_back(coefs[i][j][w]);
				if (j == vsize - 1 && m_Inter[patchID][3].first == -1) VTemp[3].emplace_back(coefs[i][j][w]);
				if (w == 0 && m_Inter[patchID][4].first == -1) VTemp[4].emplace_back(coefs[i][j][w]);
				if (w == wsize - 1 && m_Inter[patchID][5].first == -1) VTemp[5].emplace_back(coefs[i][j][w]);
			}
		}
	}
	for (int i = 0; i < 6; ++i) {
		if (VTemp[i].size() > 0) {
			V.emplace_back(VTemp[i]);
			patchIds.emplace_back(patchID,i+1);
		}
	}
}

void CollisionDetect::uniformGetPoints(CBSplineSolid& geo, std::vector<std::vector<CPoint>>& VSet, std::vector<std::vector<CPoint>>& P_N, int numpoints, int dir)
{
	std::vector<std::vector<CPoint>> V, P, GaV;
	getParam(geo, V,P,GaV, numpoints, dir);

	for (int i = 0; i < V.size() - 1; ++i)
	{
		for (int j = 0; j < V[0].size() - 1; ++j)
		{
			std::vector<CPoint> temp, temp_P;
			temp.emplace_back(V[i][j]); temp_P.emplace_back(P[i][j]);
			temp.emplace_back(V[i + 1][j]); temp_P.emplace_back(P[i+1][j]);
			temp.emplace_back(V[i][j + 1]); temp_P.emplace_back(P[i][j+1]);
			temp.emplace_back(V[i + 1][j + 1]); temp_P.emplace_back(P[i+1][j+1]);
			temp.emplace_back(GaV[i][j]);
			VSet.emplace_back(temp); P_N.emplace_back(temp_P);
		}
	}
}

void CollisionDetect::getParam(CBSplineSolid& geo, std::vector<std::vector<CPoint>>& V, std::vector<std::vector<CPoint>>& P, std::vector<std::vector<CPoint>>& GaV, int numpoints, int dir)
{
	double low = 0.0, up = 1.0;
	const double h = pow(1.0 / numpoints, 1.0 / 2);
	if (dir == 1)
	{
		for (double w = low; w < up + h; w += h) {
			std::vector<CPoint> temp, temp_P, temp_GaV;
			for (double v = low; v < up + h; v += h) {
				temp.emplace_back(geo.getPointFromParam(low, min(v, up), min(w, up)));
				temp_P.emplace_back(CPoint{ low, min(v, up), min(w, up)});
				temp_GaV.emplace_back(geo.getPointFromParam(low + Ga, min(v, up), min(w, up)));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
	else if (dir == 2)
	{
		for (double w = low; w < up + h; w += h) {
			std::vector<CPoint> temp, temp_P, temp_GaV;
			for (double v = low; v < up + h; v += h) {
				temp.emplace_back(geo.getPointFromParam(up, min(v, up), min(w, up)));
				temp_P.emplace_back(CPoint{ up, min(v, up), min(w, up) });
				temp_GaV.emplace_back(geo.getPointFromParam(up - Ga, min(v, up), min(w, up)));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
	else if (dir == 3)
	{
		for (double w = low; w < up + h; w += h) {
			std::vector<CPoint> temp, temp_P, temp_GaV;
			for (double u = low; u < up + h; u += h) {
				temp.emplace_back(geo.getPointFromParam(min(u, up), low, min(w, up)));
				temp_P.emplace_back(CPoint{ min(u, up), low, min(w, up) });
				temp_GaV.emplace_back(geo.getPointFromParam(min(u, up), low+Ga, min(w, up)));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
	else if (dir == 4)
	{
		for (double w = low; w < up + h; w += h) {
			std::vector<CPoint> temp, temp_P,temp_GaV;
			for (double u = low; u < up + h; u += h) {
				temp.emplace_back(geo.getPointFromParam(min(u, up), up, min(w, up)));
				temp_P.emplace_back(CPoint{ min(u, up), up, min(w, up) });
				temp_GaV.emplace_back(geo.getPointFromParam(min(u, up), up - Ga, min(w, up)));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
	else if (dir == 5)
	{
		for (double v = low; v < up + h; v += h) {
			std::vector<CPoint> temp, temp_P, temp_GaV;
			for (double u = low; u < up + h; u += h) {
				temp.emplace_back(geo.getPointFromParam(min(u, up), min(v, up), low));
				temp_P.emplace_back(CPoint{ min(u, up), min(v, up), low });
				temp_GaV.emplace_back(geo.getPointFromParam(min(u, up), min(v, up), low+Ga));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
	else
	{
		for (double v = low; v < up + h; v += h) {
			std::vector<CPoint> temp, temp_P, temp_GaV;
			for (double u = low; u < up + h; u += h) {
				temp.emplace_back(geo.getPointFromParam(min(u, up), min(v, up), up));
				temp_P.emplace_back(CPoint{ min(u, up), min(v, up), up });
				temp_GaV.emplace_back(geo.getPointFromParam(min(u, up), min(v, up), up-Ga));
			}
			V.emplace_back(temp);
			P.emplace_back(temp_P);
			GaV.emplace_back(temp_GaV);
		}
	}
}

AABBTree CollisionDetect::CreateTree(const std::vector<std::vector<CPoint>>& V, const std::vector<std::vector<CPoint>>& P_N)
{
	std::vector<int> indices(V.size());
	for (int i = 0; i < V.size(); ++i) indices[i] = i;
	AABBTree root(V,P_N, indices);
	return root;
}
