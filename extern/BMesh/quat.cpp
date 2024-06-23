#include "BMesh/quat.hpp"
using namespace BMesh;

Quat::Quat(){
}

Quat::~Quat() {}

void Quat::QuatWithtwoVec(Lvector& from, Lvector& to) {
    const double epsilon = 1E-10;

    const double fromSqNorm = from.squaredNorm();
    const double toSqNorm = to.squaredNorm();
    // Identity Quaternion when one vector is null
    if ((fromSqNorm < epsilon) || (toSqNorm < epsilon)) {
        q[0] = q[1] = q[2] = 0.0;
        q[3] = 1.0;
    }
    else {
        Lvector axis = from % to;
        const double axisSqNorm = axis.squaredNorm();

        // Aligned vectors, pick any axis, not aligned with from or to
        if (axisSqNorm < epsilon)
            axis = from.orthogonalVec();

        double angle = asin(sqrt(axisSqNorm / (fromSqNorm * toSqNorm)));

        if (from * to < 0.0)
            angle = M_PI - angle;

        setAxisAngle(axis, angle);
    }
}

void Quat::QuatWithtwoVec(const Lvector& a, const Lvector& b) {
    Lvector aa = a;
    Lvector bb = b;
    QuatWithtwoVec(aa, bb);
}


void Quat::setAxisAngle(Lvector& axis, double angle) {
    double norm = axis.norm();
    if (norm < 1E-8) {
        // Null rotation
        q[0] = 0.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 1.0;
    }
    else {
        const double sin_half_angle = sin(angle / 2.0);
        q[0] = sin_half_angle * axis.x / norm;
        q[1] = sin_half_angle * axis.y / norm;
        q[2] = sin_half_angle * axis.z / norm;
        q[3] = cos(angle / 2.0);
    }
}

void Quat::setAxisAngle(const Lvector& axis, double angle) {
    double norm = axis.norm();
    if (norm < 1E-8) {
        // Null rotation
        q[0] = 0.0;
        q[1] = 0.0;
        q[2] = 0.0;
        q[3] = 1.0;
    }
    else {
        const double sin_half_angle = sin(angle / 2.0);
        q[0] = sin_half_angle * axis.x / norm;
        q[1] = sin_half_angle * axis.y / norm;
        q[2] = sin_half_angle * axis.z / norm;
        q[3] = cos(angle / 2.0);
    }
}


Lvector Quat::rotate(const Lvector& v) const {
    const double q00 = 2.0 * q[0] * q[0];
    const double q11 = 2.0 * q[1] * q[1];
    const double q22 = 2.0 * q[2] * q[2];

    const double q01 = 2.0 * q[0] * q[1];
    const double q02 = 2.0 * q[0] * q[2];
    const double q03 = 2.0 * q[0] * q[3];

    const double q12 = 2.0 * q[1] * q[2];
    const double q13 = 2.0 * q[1] * q[3];

    const double q23 = 2.0 * q[2] * q[3];

    return Lvector((1.0 - q11 - q22) * v.x + (q01 - q23) * v.y + (q02 + q13) * v.z,
        (q01 + q23) * v.x + (1.0 - q22 - q00) * v.y + (q12 - q03) * v.z,
        (q02 - q13) * v.x + (q12 + q03) * v.y +
        (1.0 - q11 - q00) * v.z);
}


Lvector  Quat::operator*(const Lvector& v) {
    return rotate(v);
}