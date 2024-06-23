/*
 * @Author: your name
 * @Date: 2020-10-05 21:33:43
 * @LastEditTime: 2022-05-10 23:16:55
 * @LastEditors: lishaohsuai 1049188593@qq.com
 * @Description: In User Settings Edit
 * @FilePath: /my_skel/src/face.cc
 */
#include "triangle.hpp"
#include <fstream>
#include <sstream>
using namespace BMesh;
/**
 * @description:
 * @param {type}
 * @return {type}
 */
Triangle::Triangle(Point3D a, Point3D b, Point3D c) {
    Lvector u(a, b);
    Lvector vv(a, c);
    n = (u % vv).normalize();
    v.push_back(a);
    v.push_back(b);
    v.push_back(c);
    id = -1;
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Triangle::~Triangle() {
    //std::cout << "[EEE] \n";
}

/**
 * @description: Set face name
 * @param {type}
 * @return {type}
 */
void Triangle::setName(std::string name) { name = name; }

/**
 * @description: The constructor
 * @param {type}
 * @return {type}
 */
Triangle::Triangle() {
    v.resize(3);
    v[0] = v[1] = v[2] = Point3D(0, 0, 0);
    n = Lvector(0, 0, 0);
    id = -1;
}

/**
 * @description: Copy constructor
 * @param {type}
 * @return {type}
 */
Triangle::Triangle(Triangle& rhs) {
    if (rhs.v.size() == v.size() && rhs.v.size() == 3) {
        v[0] = rhs.v[0];
        v[1] = rhs.v[1];
        v[2] = rhs.v[2];
    }
    else {
        v.clear();
        for (int i = 0; i < rhs.v.size(); i++) {
            v.push_back(rhs.v[i]);
        }
    }

    n = rhs.n;
    id = rhs.id;
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Triangle& Triangle::operator=(const Triangle& rhs) {
    this->n = rhs.n;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
    this->id = rhs.id;
    return *this;
}
Triangle& Triangle::operator=(Triangle& rhs) {
    this->n = rhs.n;
    this->id = rhs.id;
    v.clear();
    for (int i = 0; i < rhs.v.size(); i++) {
        v.push_back(rhs.v[i]);
    }
    return *this;
}

/**
 * @description:
 * @param {type}
 * @return {type}
 */
Triangle::Triangle(const Triangle& rhs) {
    if (rhs.v.size() == v.size() && rhs.v.size() == 3) {
        v[0] = rhs.v[0];
        v[1] = rhs.v[1];
        v[2] = rhs.v[2];
    }
    else {
        v.clear();
        for (int i = 0; i < rhs.v.size(); i++) {
            v.push_back(rhs.v[i]);
        }
    }
    n = rhs.n;
    this->id = rhs.id;
}

/**
 * @description: Returns the point description of the surface
 * @param {type}
 * @return {type}
 */
std::string Triangle::getPointsString() {
    std::string str;
    for (int i = 0; i < v.size(); i++) {
        str += v[i].getXYZ();
    }
    return str;
}

/**
 * @description: You get the mean points of the faces
 * @param {type}
 * @return {type}
 */
Point3D Triangle::getMidFacePoint() {
    Point3D midPoint(0, 0, 0);
    int number = 0;
    for (int i = 0; i < v.size(); i++) {
        midPoint += v[i];
        number++;
    }
    midPoint /= (double)number;
    return midPoint;
}

/**
 * @description: usemidPoint then one triangle to be three three triangle
 * @param {type}
 * @return {type}
 */
std::vector<Triangle> Triangle::oneToThree() {
    Point3D mid = getMidFacePoint();
    std::vector<Triangle> rlt;
    rlt.push_back(Triangle(v[0], v[1], mid));
    rlt.push_back(Triangle(v[1], v[2], mid));
    rlt.push_back(Triangle(v[2], v[0], mid));
    return rlt;
}

/**
 * @description: Determine if the faces are equal
 * @param {type}
 * @return {type}
 */
bool Triangle::operator==(Triangle& rhs) {
    for (int i = 0; i < this->v.size(); i++) {
        if (std::find(rhs.v.begin(), rhs.v.end(), this->v[i]) == rhs.v.end()) {
            return false;
        }
    }
    return true;
}

/**
 * @description: You get the area of a surface
 *               Take the area of two triangles for a quadrilateral
 *               For the time being only tetrahedral
 * @param {type}
 * @return {type}
 */
double Triangle::getArea() {
    if (v.size() != 3) {
        std::cout << "[ERROR] triangle's v number is not 3\n";
        return 0;
    }
    double sumArea = 0;
    sumArea += calcAreaOfThreePoints(v[0], v[1], v[2]);
    return sumArea;
}

/**
 * @description: �²���������Ϊ�໥�������µ�����
 * @param {type}
 * @return {type}
 */
double Triangle::calcAreaOfThreePoints(Point3D& a, Point3D& b, Point3D& c) {
    //Ӧ�ú��׹�ʽ   S=1/4sqrt[(a+b+c)(a+b-c)(a+c-b)(b+c-a)]
    double lenA =
        std::sqrt(std::pow(b.x - c.x, 2) + std::pow(b.y - c.y, 2) + std::pow(b.z - c.z, 2)); // b - c ����������
    double lenB =
        std::sqrt(std::pow(a.x - c.x, 2) + std::pow(a.y - c.y, 2) + std::pow(a.z - c.z, 2)); // a - c ����������
    double lenC =
        std::sqrt(std::pow(b.x - a.x, 2) + std::pow(b.y - a.y, 2) + std::pow(b.z - a.z, 2)); // a - b ����������
    double Area = 1.0 / 4.0 *
        std::sqrt((lenA + lenB + lenC) * (lenA + lenB - lenC) * (lenA + lenC - lenB) * (lenB + lenC - lenA));
    return Area;
}

/**
 * @description: I'm going to give you a bunch of points, and then I'm going to get points that aren't on the surface
 * @param {type}
 * @return {type}
 */
std::vector<Point3D> Triangle::getOtherPointFromAFace(std::vector<Point3D>& rhs) {
    std::unordered_map<std::string, bool> um;
    for (int i = 0; i < rhs.size(); i++) {
        um[rhs[i].getXYZ()] = true;
    }
    std::vector<Point3D> rlt;
    for (int i = 0; i < v.size(); i++) {
        if (um[v[i].getXYZ()] == false) {
            rlt.push_back(v[i]);
        }
    }
    return rlt;
}

/**
 * @description: Calculate the Angle between the two vectors
 * @param {type}
 * @return {type}
 */
double Triangle::angleTwoVector(Lvector& a, Lvector& b) {
    double cosTheta = (a * b) / (a.norm() * b.norm());
    if (cosTheta > 1) {
        cosTheta = 1;
    }
    if (cosTheta < -1) {
        cosTheta = -1;
    }
    double theta = std::acos(cosTheta);
    theta = theta / PI * 180.0;
    return theta;
}

/**
 * @description: Determine if the point is within the range of the surface
 * @param {type}
 * @return {type}
 */
bool Triangle::isInTheFace(const Point3D& a) {
    // The sum of the angles around point A is 360 degrees
    Lvector t1(a, v[0]);
    Lvector t2(a, v[1]);
    Lvector t3(a, v[2]);
    double angle1 = std::fabs(angleTwoVector(t1, t2)) + std::fabs(angleTwoVector(t2, t3)) +
        std::fabs(angleTwoVector(t3, t1));
    if (std::fabs(angle1 - 360) > 10) {
        return false;
    }
    return true;
}


/**
 * @description: Get bounding box
 * @param {type}
 * @return {type}
 */
void Triangle::genFrame() {
    smallPoint = Point3D(MAXDATA, MAXDATA, MAXDATA);
    bigPoint = Point3D(MINDATA, MINDATA, MINDATA);
    for (int i = 0; i < v.size(); i++) {
        if (v[i].x < smallPoint.x) {
            smallPoint.x = v[i].x;
        }
        if (v[i].y < smallPoint.y) {
            smallPoint.y = v[i].y;
        }
        if (v[i].z < smallPoint.z) {
            smallPoint.z = v[i].z;
        }
        // ---
        if (v[i].x > bigPoint.x) {
            bigPoint.x = v[i].x;
        }
        if (v[i].y > bigPoint.y) {
            bigPoint.y = v[i].y;
        }
        if (v[i].z > bigPoint.z) {
            bigPoint.z = v[i].z;
        }
    }
}

