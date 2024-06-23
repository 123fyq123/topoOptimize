#include "triangleMesh.hpp"
#include <fstream>
#include <sstream>

using namespace BMesh;

/**
 * @description: open off triangle mesh
 * @param {type}
 * @return {type}
 */
bool TriangleMesh::openOff(std::string filename) {
    std::ifstream iff = std::ifstream(filename.c_str(), std::ios::in);
    if (iff.good()) {
        std::cout << "[ERROR] Could not open file " << filename << " for reading!" << std::endl;
        iff.close();
        return false;
    }
    std::string line;
    getline(iff, line);
    std::string los("OFF");
    if (line != los) {
        std::cout << "[ERROR] format not supported " << std::endl;
        std::cout << "los " << los << std::endl;
        std::cout << "line " << line << std::endl;
        return false;
    }
    getline(iff, line);
    std::stringstream sstr;
    sstr.str(line);
    sstr >> pointNum;
    sstr >> faceNum;
    sstr.clear();

    for (int i = 0; i < pointNum; i++) {
        getline(iff, line);
        sstr.str(line);
        double b1, b2, b3;
        sstr >> b1;
        sstr >> b2;
        sstr >> b3;
        v.push_back(Point3D(b1, b2, b3));
    }

    int tmp;
    int pv1;
    int pv2;
    int pv3;

    sstr.clear();
    for (int i = 0; i < faceNum; i++) {
        getline(iff, line);
        sstr.str(line);
        sstr >> tmp;
        if (tmp != 3) {
            std::cout << "[ERROR] format not supported_ " << std::endl;
            return false;
        }
        sstr >> pv1;
        sstr >> pv2;
        sstr >> pv3;
        triangle t(pv1, pv3, pv3);
        sstr.clear();
        allTriangle.push_back(t);
    }
    iff.close();
    return true;
    //sstr >>

}