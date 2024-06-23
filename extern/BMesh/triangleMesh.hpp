/*
 * @Author: your name
 * @Date: 2020-10-05 21:30:52
 * @LastEditTime: 2020-11-26 21:05:38
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /my_skel/header/face.hpp
 */
#ifndef __TRIANGLEMESH_HPP__
#define __TRIANGLEMESH_HPP__
#include "line.hpp"
#include "local_vector.hpp"
#include "point.hpp"
#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>
#include "triangle.hpp"
#include <map>

namespace BMesh {
    class TriangleMesh {
        const double MAXDATA = 10000000;
        const double MINDATA = -MAXDATA;
        const double PI = 3.1415926;
        const double ESP = 0.000001;
    public:
        bool openOff(std::string filename);
    private:
        std::vector<Point3D> v; // save all point
        std::vector<BMesh::triangle> allTriangle;
        int pointNum;
        int faceNum;
        std::map<int, std::vector<int> > mm; // point in triangle
    };
} // namespace BMesh

#endif