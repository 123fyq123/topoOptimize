#include "BMesh/fixFlip.hpp"


bool FixFlip::isFlip(triangle b, int a1, int a2) {
    if (b.pv1 == a1 && b.pv2 == a2) {
        return true;
    }
    if (b.pv2 == a1 && b.pv3 == a2) {
        return true;
    }
    if (b.pv3 == a1 && b.pv1 == a2) {
        return true;
    }
    return false;
}

void FixFlip::fixFlip(triangle& a) {
    int a1 = a.pv1;
    int a2 = a.pv2;
    int a3 = a.pv3;
    a.pv1 = a3;
    a.pv3 = a1;
}


/**
 * @description:
 * @param {set<int> &s1, std::set<int>} &s2
 * @param {int} a1
 * @param {int} a2
 * @return 是否修复了, 或者除了一种情况两个三角形都是true的情况
 */
bool FixFlip::refine(std::set<int>& s1, std::set<int>& s2, int a1, int a2, int& falseTri) {
    std::set<int> c;
    std::set_intersection(s1.begin(), s1.end(), s2.begin(), s2.end(), inserter(c, c.begin()));
    if (c.size() != 2) { // 对于封闭的曲面而言是对的
        std::cout << "[ERROR] c size not equal 2 " << c.size() << "\n";
        for (auto it : c) {
            std::cout << allTriangle[it].pv1 << " " << allTriangle[it].pv2 << " " << allTriangle[it].pv3 << "\n====\n";
            std::cout << p[allTriangle[it].pv1].x << " " << p[allTriangle[it].pv1].y << " " << p[allTriangle[it].pv1].z << "\n";
            std::cout << p[allTriangle[it].pv2].x << " " << p[allTriangle[it].pv2].y << " " << p[allTriangle[it].pv2].z << "\n";
            std::cout << p[allTriangle[it].pv3].x << " " << p[allTriangle[it].pv3].y << " " << p[allTriangle[it].pv3].z << "\n";
        }
        return false;
    }
    std::vector<int> cv;
    for (auto it : c) {
        cv.push_back(it);
    }
    int trueTri = -1;
    falseTri = -1;
    if (check[cv[0]] == true && check[cv[1]] == false) {
        trueTri = cv[0];
        falseTri = cv[1];
    }
    else if (check[cv[0]] == false && check[cv[1]] == true) {
        trueTri = cv[1];
        falseTri = cv[0];
    }
    else {
        //std::cout << "[DEBUG] maybe this two triangle is true??\n";
        return false;
    }
    if (falseTri != -1) {
        if (isFlip(allTriangle[falseTri], a1, a2)) {
            fixFlip(allTriangle[falseTri]);
        }
    }
    return true;
}


void FixFlip::BFS(int i) {
    check[i] = true;
    triQ.push(i);
    while (!triQ.empty()) {
        int front = triQ.front();
        triQ.pop();
        triangle t = allTriangle[front];
        std::set<int> s1;
        std::set<int> s2;
        std::set<int> s3;
        std::set<int> c;
        for (auto it : mm[t.pv1]) {
            s1.insert(it);
        }
        for (auto it : mm[t.pv2]) {
            s2.insert(it);
        }
        for (auto it : mm[t.pv3]) {
            s3.insert(it);
        }
        int falseTri;
        if (refine(s1, s2, t.pv1, t.pv2, falseTri)) {
            triQ.push(falseTri);
            check[falseTri] = true;
        }
        if (refine(s2, s3, t.pv2, t.pv3, falseTri)) {
            triQ.push(falseTri);
            check[falseTri] = true;
        }
        if (refine(s3, s1, t.pv3, t.pv1, falseTri)) {
            triQ.push(falseTri);
            check[falseTri] = true;
        }
    }
}

void FixFlip::BFSTraverse() {
    // 访问数组初始化
    for (int i = 0; i < faceNum; i++) {
        if (!check[i]) {
            BFS(i);
        }
    }
}

void FixFlip::run(std::string filenameIN, std::string filenameOut) {
    std::ifstream iff = std::ifstream(filenameIN.c_str(), std::ios::in);
    std::ofstream off = std::ofstream(filenameOut.c_str(), std::ios::out);
    off << "OFF" << std::endl;
    if (!off.good()) {
        std::cout << "[Error] Could not open file " << filenameOut << " for writing!" << std::endl;
        off.close();
        return;
    }
    if (!iff.good()) {
        std::cout << "[ERROR] could not open file " << filenameIN << " for reading!" << std::endl;
        iff.close();
        return;
    }
    std::string line;
    getline(iff, line);
    std::string los("OFF");
    if (line != los) {
        std::cout << "[ERROR] format not supported " << std::endl;
        std::cout << "los " << los << std::endl;
        std::cout << "line " << line << std::endl;
        return;
    }
    getline(iff, line);
    stringstream sstr;
    sstr.str(line);

    sstr >> pointNum;
    sstr >> faceNum;
    sstr.clear();

    off << line << "\n";
    BMesh::Point3D center(0, 0, 0);
    for (int i = 0; i < pointNum; i++) {
        getline(iff, line);
        sstr.str(line);
        double b1, b2, b3;
        sstr >> b1;
        sstr >> b2;
        sstr >> b3;
        p.push_back(BMesh::Point3D(b1, b2, b3));
        center += BMesh::Point3D(b1, b2, b3);
        off << line << "\n";
    }
    center /= pointNum;


    std::map<triangle, int> mi;
    int tmp;
    int pv1;
    int pv2;
    int pv3;
    sstr.clear();
    // 
    for (int i = 0; i < faceNum; i++) {
        getline(iff, line);
        sstr.str(line);
        sstr >> tmp;
        if (tmp != 3) {
            std::cout << "[ERROR] format not supported_" << std::endl;
            return;
        }
        sstr >> pv1;
        sstr >> pv2;
        sstr >> pv3;
        triangle t(pv1, pv2, pv3);
        if (i == 0) {
            if (!IsFurtherAway(center, t)) {
                std::swap(t.pv1, t.pv2);
            }
            std::swap(t.pv1, t.pv2);
        }
        sstr.clear();
        mm[pv1].push_back(i);
        mm[pv2].push_back(i);
        mm[pv3].push_back(i);
        allTriangle.push_back(t);
        mi[t] = i;
    }
    check.resize(allTriangle.size(), false);
    BFSTraverse();
    for (auto it : allTriangle) {
        off << "3 " << it.pv1 << " " << it.pv2 << " " << it.pv3 << "\n";
    }
    off.close();
    iff.close();
}


/**
 * @description: Determine whether the vertex of the first face extends further away from the center point
 * if furtherAway , return true, else return false
 * @param {type}
 * @return {type}
 */
bool FixFlip::IsFurtherAway(BMesh::Point3D& center, triangle& t) {
    BMesh::Lvector l1(p[t.pv1], p[t.pv2]);
    BMesh::Lvector l2(p[t.pv1], p[t.pv3]);
    BMesh::Lvector l3 = l1 % l2;
    BMesh::Point3D fpoint = l3 * 0.1 + p[t.pv1];
    double distance1 = BMesh::Geometry::distanceTwoPoint(fpoint, center);
    double distance2 = BMesh::Geometry::distanceTwoPoint(p[t.pv1], center);
    if (distance1 > distance2) {
        return true;
    }
    return false;
}