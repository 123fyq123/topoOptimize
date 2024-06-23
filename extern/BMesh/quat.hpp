/*
 * @Author: your name
 * @Date: 2020-10-05 08:25:37
 * @LastEditTime: 2020-11-26 21:15:03
 * @LastEditors: Please set LastEditors
 * @Description: In User Settings Edit
 * @FilePath: /my_skel/header/point.hpp
 */
#pragma once

#include <cmath>
#include <iostream>
#include <string>
#include <vector>
#include <glm/gtc/quaternion.hpp>
#include <glm/gtx/quaternion.hpp>
#include "BMesh/local_vector.hpp"
namespace BMesh {
    class Quat {
    public:
        explicit Quat();
        ~Quat();

        Lvector operator*(const Lvector& v);
        void QuatWithtwoVec(Lvector& a, Lvector& b);
        void QuatWithtwoVec(const Lvector& a, const Lvector& b);
        void setAxisAngle( Lvector& axis, double angle);
        void setAxisAngle(const Lvector& axis, double angle);
        Lvector rotate(const Lvector& v) const;
        


    public:
        double q[4];
        int id;

    private:
        const double ESP = 1.0e-7;
    };
} // namespace BMesh
