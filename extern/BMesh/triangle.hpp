/*
 * @Author: your name
 * @Date: 2020-10-05 21:30:52
 * @LastEditTime: 2020-11-26 21:05:38
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /my_skel/header/face.hpp
 */
#ifndef __TRIANGLE_HPP__
#define __TRIANGLE_HPP__
#include "line.hpp"
#include "local_vector.hpp"
#include "point.hpp"
#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

namespace BMesh {
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
    class Triangle {
        const double MAXDATA = 10000000;
        const double MINDATA = -MAXDATA;
#ifndef PI
        const double PI = 3.1415926;
#endif
        const double ESP = 0.000001;
    public:
        Triangle(Point3D a, Point3D b, Point3D c);
        explicit Triangle();
        Triangle(Triangle& rhs);
        Triangle(const Triangle& rhs);
        ~Triangle();

    public:
        Triangle& operator=(const Triangle& rhs);
        Triangle& operator=(Triangle& rhs);
        bool operator==(Triangle& rhs);

    public:
        
        void setName(std::string name);
        std::string getPointsString();
        Point3D getMidFacePoint();
        std::vector<Triangle> oneToThree();
        double getArea();
        std::vector<Point3D> getOtherPointFromAFace(std::vector<Point3D>& rhs);
        bool isInTheFace(const Point3D& a);
        void genFrame();

    private:
        double calcAreaOfThreePoints(Point3D& a, Point3D& b, Point3D& c);
        double angleTwoVector(Lvector& a, Lvector& b);

    public:
        Lvector n;              // Normal vector of surface
        std::vector<Point3D> v; // The point on the surface
        std::string name;
        int id;
        Point3D smallPoint;
        Point3D bigPoint;
    };
} // namespace BMesh

#endif