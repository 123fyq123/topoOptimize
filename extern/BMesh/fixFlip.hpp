#pragma once
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <map>
#include <vector>
#include <queue>
#include <set>
#include <algorithm>
#include "BMesh/point.hpp"
#include "BMesh/local_vector.hpp"
#include "BMesh/geometry.hpp"
using namespace std;
struct triangle {
    int pv1;
    int pv2;
    int pv3;
    triangle(int pv1_, int pv2_, int pv3_) :pv1(pv1_), pv2(pv2_), pv3(pv3_) {};
    bool operator == (const triangle& rhs) {
        return pv1 == rhs.pv1 && pv2 == rhs.pv2 && pv3 == rhs.pv3;
    }
    bool operator < (const triangle& rhs) {
        return (pv1 + pv2 + pv3) < (rhs.pv1 + rhs.pv2 + rhs.pv3);
    }
    bool operator < (const triangle& rhs) const {
        return (pv1 + pv2 + pv3) < (rhs.pv1 + rhs.pv2 + rhs.pv3);
    }
};
class FixFlip {
public:
	std::vector<BMesh::Point3D> p;
    int pointNum;
    int faceNum;
    std::queue<int> triQ;
    std::vector<bool> check;
    std::vector<triangle> allTriangle;
    std::map<int, std::vector<int> > mm;
public:
    bool isFlip(triangle b, int a1, int a2);
    void fixFlip(triangle& a);
    bool refine(std::set<int>& s1, std::set<int>& s2, int a1, int a2, int& falseTri);
    void BFS(int i);
    void BFSTraverse();
    void run(std::string filenameIN, std::string filenameOut);
    bool IsFurtherAway(BMesh::Point3D& center, triangle& t);
};