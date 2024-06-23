/*
 * @Author: your name
 * @Date: 2020-10-05 21:02:30
 * @LastEditTime: 2022-05-10 23:18:08
 * @LastEditors: lishaohsuai 1049188593@qq.com
 * @Description: In User Settings Edit
 * @FilePath: /my_skel/src/line.cc
 */
#include "rays.hpp"
using namespace BMesh;
/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays::Rays(Point3D a, Point3D b) {
    Lvector tt(a, b);
    l = tt;
    l.normalize();
    v.push_back(a);
    v.push_back(b);
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays::~Rays() {}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays::Rays() { l = Lvector(0, 0, 0); }

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays::Rays(Rays& rhs) {
    l = rhs.l;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays::Rays(const Rays& rhs) {
    l = rhs.l;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Rays& Rays::operator=(const Rays& rhs) {
    this->l = rhs.l;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
    return *this;
}
Rays& Rays::operator=(Rays& rhs) {
    this->l = rhs.l;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
    return *this;
}

/**
 * @description: �ж������Ƿ�����
 * @param {type}
 * @return {type}
 */
bool Rays::operator==(Rays& rhs) {
    for (int i = 0; i < this->v.size(); i++) {
        if (rhs.v[i] != this->v[i]) {
            return false;
        }
    }
    return true;
}

/**
 * @description: Rays not equal
 * @param {type}
 * @return {type}
 */
bool Rays::operator!=(Rays& rhs) {
    if ((*this) == rhs) {
        return false;
    }
    return true;
}

/**
 * @description: Is the point online
 * If there is a line, then the vector formed by a point and any point on the line is either equal or opposite
 * @param {type}
 * @return {type}
 */
bool Rays::isPointInRays(const Point3D& p) {
    Lvector tt(v[0], p);
    tt.normalize();
    if (fabs(angleTwoVector(tt, l)) < 5) {
        return true;
    }
    return false;
    //if (tt != l) {
    //    return false;
    //}
    //return true;
}

double Rays::angleTwoVector(Lvector& a, Lvector& b) {
    double cosTheta = (a * b) / (a.norm() * b.norm());
    if (cosTheta > 1) {
        cosTheta = 1;
    }
    if (cosTheta < -1) {
        cosTheta = -1;
    }
    double theta = std::acos(cosTheta);
    theta = theta / 3.1415926 * 180.0;
    return theta;
}
