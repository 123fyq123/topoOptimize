#pragma once
#include<vector>
#include<Eigen/Core>
#include <Eigen/Dense>
#include <Geom/CPoint.h>
#include"CBSplineSolid.h"
#include<map>
#include<unordered_set>
#include"AABBTree.h"
#include"ToolFunc.h"

const double Ga = 0.1;

class CollisionDetect {
    std::vector<std::unordered_set<int>> m_ConnectPatch;
    std::vector<std::vector<pair<int, int>>> m_Inter;

    std::map<int, AABBTree> m_Tree;
public:
    CollisionDetect(std::vector<std::unordered_set<int>>& connect, std::vector<std::vector<pair<int, int>>>& inter)
        :m_ConnectPatch(connect),m_Inter(inter){}

	~CollisionDetect();

    bool ObjectDetect(std::vector<CBSplineSolid>& geo, std::vector<std::vector<CollideInfo>>& collidInfo, std::vector<pair<int, int>>& patchIds);

private:
    std::vector<std::pair<int, int>> findCollisions(const std::vector<std::vector<CPoint>>& V);

    bool refineFindCollisions(const std::vector<std::vector<CPoint>>& V);

    void boundCofes(const CBSplineSolid& geo, std::vector<std::vector<CPoint>>& V, int patchID, std::vector<pair<int,int>>& patchIds);

    void uniformGetPoints(CBSplineSolid& geo, std::vector<std::vector<CPoint>>& VSet, std::vector<std::vector<CPoint>>& P_N,int numpoints, int dir);

    void getParam(CBSplineSolid& geo, std::vector<std::vector<CPoint>>& V, std::vector<std::vector<CPoint>>& P,std::vector<std::vector<CPoint>>&GaV, int numpoints, int dir);

    AABBTree CreateTree(const std::vector<std::vector<CPoint>>& V, const std::vector<std::vector<CPoint>>& P_N);
};