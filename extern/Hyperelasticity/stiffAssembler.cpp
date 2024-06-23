#include "stiffAssembler.h"
#include "collision_detect.h"
#include "OutputParaview.h"
#include <omp.h>

stiffAssembler::stiffAssembler(const std::vector<CBSplineSolid> &m_patch,
                               const std::vector<boundary_condition> &m_dirichletSides,
                               const std::vector<boundary_condition> &m_neumannSides)
        : assembler(), geometry(m_patch), dirichletSides(m_dirichletSides), neumanSides(m_neumannSides) {
    m_system = SparseSystem();
    lambda = E * pr / ((1. + pr) * (1. - 2. * pr));    //拉梅常数
    mu = E / (2. * (1. + pr));
    refresh_sides();
    if (isPreCompute) preComputeBasis();
}

stiffAssembler::stiffAssembler(const std::vector<CBSplineSolid> &m_patch,
                               const std::vector<boundary_point> &m_dirichletPoints,
                               const std::vector<boundary_point> &m_neumannPoints)
        : assembler(), geometry(m_patch), dirichletPoints(m_dirichletPoints), neumannPoints(m_neumannPoints) {

    lambda = E * pr / ((1. + pr) * (1. - 2. * pr));    //拉梅常数
    mu = E / (2. * (1. + pr));
    refresh_points();
    if (isPreCompute) preComputeBasis();
}

void stiffAssembler::refresh_sides() {
    m_system = SparseSystem(geometry, dirichletSides, dim);
    extraRhs.setZero(m_system.rhs().rows(), 1);
    int num = 0;
    for (int i = 0; i < geometry.size(); ++i) {
        num += geometry[i].getSizeControlPoints();
        controlPointsSumNum.push_back(num);
    }

}

void stiffAssembler::refresh_points() {
    auto res = SparseSystem(geometry, dirichletPoints, dim);
    m_system = res;
    extraRhs.setZero(m_system.rhs().rows(), 1);
    int num = 0;
    for (int i = 0; i < geometry.size(); ++i) {
        num += geometry[i].getSizeControlPoints();
        controlPointsSumNum.push_back(num);
    }
}

int stiffAssembler::numDofs() {
    return m_system.getControlNums() * dim;
}

void stiffAssembler::stiffMatrix() {
    //遍历所有面片
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim), C;
    ToolFunc::setC(C, identity, lambda, mu);
#pragma omp parallel for num_threads(24)
    for (int i = 0; i < geometry.size(); ++i) {
        // 某一片模型
        int N_D;
        std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
        std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
        CBSplineSolid *geo = &geometry[i];
        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        //std::vector<Eigen::MatrixXd> &basisValue = basisValueDer;
        std::vector<Eigen::MatrixXd> basisValue(2);
        std::vector<int> localIndex(gaussNum);
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
        {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
            {
                //高斯点选择
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                    int cnt =
                            w_interval * (vknots.size() - 1) * (uknots.size() - 1) + v_interval * (uknots.size() - 1) +
                            u_interval;
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //局部索引号
                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    //基函数
                    if (!isPreCompute)
                        geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                    else
                        basisValue = basisValueDer[cnt];
                    Eigen::MatrixXd localMats;
                    localMats.setZero(dim * N_D, dim * N_D);
                    Eigen::MatrixXd localMats_truth;
                    localMats_truth.setZero(dim * N_D, dim * N_D);

                    bool bVerify = true;

                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                Eigen::MatrixXd jac, P;
                                //geo->getJacobi(0.000, 0.000, 0.000, localIndex, jac);
                                //cout << jac;
                                geo->getJacobi(basisValue[1].col(q), localIndex, jac);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());
                                Eigen::MatrixXd physGrad;
                                Eigen::MatrixXd B_i, B_j, tempK, K;
                                ToolFunc::transformGradients(jac, q, basisValue[1], physGrad, dim);
                                // 局部矩阵计算 装配
                                double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
                                for (int j = 0; j < N_D; j++) {
                                    m_Rx = physGrad.col(j)(0);
                                    m_Ry = physGrad.col(j)(1);
                                    m_Rz = physGrad.col(j)(2);

                                    ToolFunc::setB(B_i, identity, physGrad.col(j));
                                    tempK = B_i.transpose() * C;
                                    for (int k = 0; k < N_D; k++) {
                                        ToolFunc::setB(B_j, identity, physGrad.col(k));
                                        m_Qx = physGrad.col(k)(0);
                                        m_Qy = physGrad.col(k)(1);
                                        m_Qz = physGrad.col(k)(2);

                                        if (bVerify) {
                                            K = tempK * B_j;
                                            for (int di = 0; di < dim; di++) {
                                                for (int dj = 0; dj < dim;
                                                     dj++) {
                                                    localMats_truth(di * N_D + j, dj * N_D + k) +=
                                                            weight * K(di, dj);
                                                }
                                            }
                                        }

                                        //localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
                                        localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx +
                                                                     C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
                                        //localMats(j, k + N_D) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
                                        localMats(j, k + N_D) +=
                                                weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
                                        //localMats(j, k + 2 * N_D) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        localMats(j, k + 2 * N_D) +=
                                                weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx);
                                        //localMats(j + N_D, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
                                        localMats(j + N_D, k) +=
                                                weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
                                        //localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
                                        localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy +
                                                                                 C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
                                        //localMats(j + N_D, k + 2 * N_D) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        localMats(j + N_D, k + 2 * N_D) +=
                                                weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        //localMats(j + 2 * N_D, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
                                        localMats(j + 2 * N_D, k) +=
                                                weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
                                        //localMats(j + 2 * N_D, k + N_D) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
                                        localMats(j + 2 * N_D, k + N_D) +=
                                                weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
                                        //localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));
                                        localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz +
                                                                                         C(3, 3) *
                                                                                         (m_Ry * m_Qy + m_Rx * m_Qx));
                                    }
                                }

                            }
                        }
                    }

                    if (bVerify) {
                        for (int i = 0; i < 192; ++i) {
                            for (int j = 0; j < 192; ++j) {
                                if (abs(localMats(i, j) - localMats_truth(i, j)) > 0.01) {
                                    printf("error in (%d, %d), with mat-%.5f || truth-%.5f", i, j, localMats(i, j),
                                           localMats_truth(i, j));
                                }

                                assert(abs(localMats(i, j) - localMats_truth(i, j)) < 1e-6);
                            }
                        }
                    }

                    std::cout << "Matrix Test Pass!!!!" << std::endl;
                    //cout << localMats << endl;
#pragma omp critical
                    m_system.pushToMatrix(i, localIndex, localMats);

                }
            }
        }
    }
}

void stiffAssembler::test_stiffMatrix(std::vector<Eigen::MatrixXd> &localMatsVec) {
    //遍历所有面片
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim), C;
    ToolFunc::setC(C, identity, lambda, mu);

    localMatsVec.resize(geometry.size());

#pragma omp parallel for num_threads(24)
    for (int i = 0; i < geometry.size(); ++i) {
        // 某一片模型
        int N_D;
        std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
        std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
        CBSplineSolid *geo = &geometry[i];
        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        //std::vector<Eigen::MatrixXd> &basisValue = basisValueDer;
        std::vector<Eigen::MatrixXd> basisValue(2);
        std::vector<int> localIndex(gaussNum);
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
        {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
            {
                //高斯点选择
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                    int cnt =
                            w_interval * (vknots.size() - 1) * (uknots.size() - 1) + v_interval * (uknots.size() - 1) +
                            u_interval;
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //局部索引号
                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    //基函数
                    if (!isPreCompute)
                        geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                    else
                        basisValue = basisValueDer[cnt];
                    Eigen::MatrixXd localMats;
                    localMats.setZero(dim * N_D, dim * N_D);
                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                Eigen::MatrixXd jac, P;
                                //geo->getJacobi(0.000, 0.000, 0.000, localIndex, jac);
                                //cout << jac;
                                geo->getJacobi(basisValue[1].col(q), localIndex, jac);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());
                                Eigen::MatrixXd physGrad;
                                // B_i, B_j, tempK, K;
                                ToolFunc::transformGradients(jac, q, basisValue[1], physGrad, dim);
                                // 局部矩阵计算 装配
                                double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;
                                for (int j = 0; j < N_D; j++) {
                                    m_Rx = physGrad.col(j)(0);
                                    m_Ry = physGrad.col(j)(1);
                                    m_Rz = physGrad.col(j)(2);

                                    //ToolFunc::setB(B_i, identity, physGrad.col(j));
                                    //tempK = B_i.transpose() * C;
                                    for (int k = 0; k < N_D; k++) {
                                        //ToolFunc::setB(B_j, identity, physGrad.col(k));
                                        m_Qx = physGrad.col(k)(0);
                                        m_Qy = physGrad.col(k)(1);
                                        m_Qz = physGrad.col(k)(2);

                                        //K = tempK * B_j;
                                        //for (int di = 0; di < dim; di++) {
                                        //	for (int dj = 0; dj < dim; dj++) {
                                        //		localMats(di * N_D + j, dj * N_D + k) += weight * K(di, dj);
                                        //	}
                                        //}
                                        //localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
                                        localMats(j, k) += weight * (C(0, 0) * m_Rx * m_Qx +
                                                                     C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
                                        //localMats(j, k + N_D) += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
                                        localMats(j, k + N_D) +=
                                                weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
                                        //localMats(j, k + 2 * N_D) += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        localMats(j, k + 2 * N_D) +=
                                                weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx);
                                        //localMats(j + N_D, k) += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
                                        localMats(j + N_D, k) +=
                                                weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
                                        //localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
                                        localMats(j + N_D, k + N_D) += weight * (C(0, 0) * m_Ry * m_Qy +
                                                                                 C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
                                        //localMats(j + N_D, k + 2 * N_D) += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        localMats(j + N_D, k + 2 * N_D) +=
                                                weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                        //localMats(j + 2 * N_D, k) += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
                                        localMats(j + 2 * N_D, k) +=
                                                weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
                                        //localMats(j + 2 * N_D, k + N_D) += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
                                        localMats(j + 2 * N_D, k + N_D) +=
                                                weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
                                        //localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));
                                        localMats(j + 2 * N_D, k + 2 * N_D) += weight * (C(0, 0) * m_Rz * m_Qz +
                                                                                         C(3, 3) *
                                                                                         (m_Ry * m_Qy + m_Rx * m_Qx));

                                    }
                                }

                            }
                        }
                    }
                    //cout << localMats << endl;
#pragma omp critical
                    localMatsVec[i] = localMats;
                    m_system.pushToMatrix(i, localIndex, localMats);
                }
            }
        }
    }
}


void stiffAssembler::stiffNonlinearMatrix(const Eigen::MatrixXd &displacement, bool isAssembleMatrix) {
    //遍历所有面片
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim);
    //C;
    double CD[6][6];
    ToolFunc::setC<6>(CD, identity, lambda, mu);
#pragma omp parallel for num_threads(24)
    for (int i = 0; i < geometry.size(); ++i) {
        int N_D;
        std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
        std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
        // 某一片模型
        CBSplineSolid *geo = &geometry[i];
        int num = geometry[i].getSizeControlPoints();
        Eigen::MatrixXd blockU = displacement.block(controlPointsSumNum[i] - num, 0, num, 3);

        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        //std::vector<Eigen::MatrixXd> &basisValue = basisValueDer;
        std::vector<Eigen::MatrixXd> basisValue(2);
        std::vector<int> localIndex(gaussNum);
        Eigen::MatrixXd forceValues;
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
        {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
            {
                //高斯点选择
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                    int cnt =
                            w_interval * (vknots.size() - 1) * (uknots.size() - 1) + v_interval * (uknots.size() - 1) +
                            u_interval;
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //局部索引号
                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    //基函数
                    if (!isPreCompute)
                        geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                    else
                        basisValue = basisValueDer[cnt];
                    //TODO
                    bodyforce.eval(p_u, p_v, p_w, forceValues);
                    Eigen::MatrixXd localMats, localRhs, jacGeo, jacDis, F, RCG, E, S;
                    // materialTangent;
                    Eigen::VectorXd Svec, geometricTangentTemp;
                    // localResidual;
                    Eigen::MatrixXd physGrad, phyDispJac;
                    // materialTangentTemp, B_i, B_j;
                    double Bi[3][6], Bj[6][3], materialTangentTemp1[3][6], materialTangent1[3][3];
                    double Bi_L[3][3], Bi_R[3][3], Bj_L[3][3], Bj_R[3][3], materialTangentTemp1L[3][3], materialTangentTemp1R[3][3];
                    if (isAssembleMatrix)
                        localMats.setZero(dim * N_D, dim * N_D);
                    localRhs.setZero(dim * N_D, 1);
                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                //geo->getJacobiU(p_u[ii], p_v[jj], p_w[kk], localIndex,blockU, jacGeo, jacDis);  // 几何雅可比X/ξ  位移雅可比 U/ξ
                                geo->getJacobiU(basisValue[1].col(q), localIndex, blockU, jacGeo, jacDis);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jacGeo.determinant()); //权值
                                ToolFunc::transformGradients(jacGeo, q, basisValue[1], physGrad, dim);
                                phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());
                                F = identity + phyDispJac;   // F = I + du/dx, deformation gradient
                                double J = F.determinant();

                                RCG = F.transpose() * F;
                                E = 0.5 * (RCG - identity);

                                if (m_materLaw.first == material_law::saint_venant_kirchhoff) {
                                    S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
                                    ToolFunc::voigtStressToK(S, Svec);                              // 应力张量
                                } else {    // 其它材料模型

                                }
                                // 局部矩阵计算 装配
                                for (size_t j = 0; j < N_D; j++) {
                                    //B = BL + BN 即小位移矩阵+大位移矩阵
                                    //ToolFunc::setB(B_i, F, physGrad.col(j));
                                    ToolFunc::setBT(Bi, F, physGrad.col(j));

                                    if (isAssembleMatrix) {
                                        //materialTangentTemp = B_i.transpose() * C;
                                        ToolFunc::matMutly<3, 6, 6>(Bi, CD, materialTangentTemp1);

                                        geometricTangentTemp = S * physGrad.col(j);
                                        //loop over active basis functions (v_j)
                                        for (size_t k = 0; k < N_D; ++k) {
                                            //ToolFunc::setB(B_j, F, physGrad.col(k));
                                            ToolFunc::setB(Bj, F, physGrad.col(k));

                                            //materialTangent = materialTangentTemp * B_j;

                                            ToolFunc::matMutly<3, 6, 3>(materialTangentTemp1, Bj, materialTangent1);


                                            double geometricTangent =
                                                    geometricTangentTemp.transpose() * physGrad.col(k);

                                            for (short d = 0; d < dim; ++d)
                                                materialTangent1[d][d] += geometricTangent;

                                            for (short di = 0; di < dim; ++di)
                                                for (short dj = 0; dj < dim; ++dj)
                                                    localMats(di * N_D + j, dj * N_D + k) +=
                                                            weight * materialTangent1[di][dj];
                                            //localMats(di * N_D + j, dj * N_D + k) += weight * materialTangent(di, dj);
                                        }
                                    }

                                    double localResiduals[3] = {0};
                                    for (int d = 0; d < 6; ++d) {
                                        localResiduals[0] += Bi[0][d] * Svec(d);
                                        localResiduals[1] += Bi[1][d] * Svec(d);
                                        localResiduals[2] += Bi[2][d] * Svec(d);
                                    }
                                    //localResidual = B_i.transpose() * Svec;
                                    for (short d = 0; d < dim; d++)
                                        localRhs(d * N_D + j) -= weight * localResiduals[d];
                                }
                                for (short d = 0; d < dim; ++d)
                                    localRhs.middleRows(d * N_D, N_D).noalias() +=
                                            weight * forceScaling * forceValues(d, q) * basisValue[0].col(q);
                            }
                        }
                    }
                    // 单元刚度矩阵组装总刚度矩阵
#pragma omp critical
                    m_system.pushToSystem(i, localIndex, localRhs, localMats, isAssembleMatrix);

                }
            }
        }
    }
}

void stiffAssembler::test_stiffNonlinearMatrix(const Eigen::MatrixXd &displacement,
                                               bool isAssembleMatrix,
                                               std::vector<Eigen::MatrixXd> &localMatsVec) {
    localMatsVec.resize(geometry.size());

    //遍历所有面片
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim);
    //C;
    double CD[6][6];
    ToolFunc::setC<6>(CD, identity, lambda, mu);
#pragma omp parallel for num_threads(24)
    for (int i = 0; i < geometry.size(); ++i) {
        int N_D;
        std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
        std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
        // 某一片模型
        CBSplineSolid *geo = &geometry[i];
        int num = geometry[i].getSizeControlPoints();
        Eigen::MatrixXd blockU = displacement.block(controlPointsSumNum[i] - num, 0, num, 3);

        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        //std::vector<Eigen::MatrixXd> &basisValue = basisValueDer;
        std::vector<Eigen::MatrixXd> basisValue(2);
        std::vector<int> localIndex(gaussNum);
        Eigen::MatrixXd forceValues;
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
        {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
            {
                //高斯点选择
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                    int cnt =
                            w_interval * (vknots.size() - 1) * (uknots.size() - 1) + v_interval * (uknots.size() - 1) +
                            u_interval;
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //局部索引号
                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    //基函数
                    if (!isPreCompute)
                        geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                    else
                        basisValue = basisValueDer[cnt];
                    //TODO
                    bodyforce.eval(p_u, p_v, p_w, forceValues);
                    Eigen::MatrixXd localMats, localRhs, jacGeo, jacDis, F, RCG, E, S;
                    // materialTangent;
                    Eigen::VectorXd Svec, geometricTangentTemp;
                    // localResidual;
                    Eigen::MatrixXd physGrad, phyDispJac;
                    // materialTangentTemp, B_i, B_j;
                    double Bi[3][6], Bj[6][3], materialTangentTemp1[3][6], materialTangent1[3][3];
                    double Bi_L[3][3], Bi_R[3][3], Bj_L[3][3], Bj_R[3][3], materialTangentTemp1L[3][3], materialTangentTemp1R[3][3];
                    if (isAssembleMatrix)
                        localMats.setZero(dim * N_D, dim * N_D);
                    localRhs.setZero(dim * N_D, 1);
                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                //geo->getJacobiU(p_u[ii], p_v[jj], p_w[kk], localIndex,blockU, jacGeo, jacDis);  // 几何雅可比X/ξ  位移雅可比 U/ξ
                                geo->getJacobiU(basisValue[1].col(q), localIndex, blockU, jacGeo, jacDis);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jacGeo.determinant()); //权值
                                ToolFunc::transformGradients(jacGeo, q, basisValue[1], physGrad, dim);
                                phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());
                                F = identity + phyDispJac;   // F = I + du/dx
                                double J = F.determinant();

                                RCG = F.transpose() * F;
                                E = 0.5 * (RCG - identity);

                                if (m_materLaw.first == material_law::saint_venant_kirchhoff) {
                                    S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
                                    ToolFunc::voigtStressToK(S, Svec);                              // 应力张量
                                } else {    // 其它材料模型

                                }
                                // 局部矩阵计算 装配
                                for (size_t j = 0; j < N_D; j++) {
                                    //B = BL + BN 即小位移矩阵+大位移矩阵
                                    //ToolFunc::setB(B_i, F, physGrad.col(j));
                                    ToolFunc::setBT(Bi, F, physGrad.col(j));

                                    if (isAssembleMatrix) {
                                        //materialTangentTemp = B_i.transpose() * C;
                                        ToolFunc::matMutly<3, 6, 6>(Bi, CD, materialTangentTemp1);

                                        geometricTangentTemp = S * physGrad.col(j);
                                        //loop over active basis functions (v_j)
                                        for (size_t k = 0; k < N_D; ++k) {
                                            //ToolFunc::setB(B_j, F, physGrad.col(k));
                                            ToolFunc::setB(Bj, F, physGrad.col(k));

                                            //materialTangent = materialTangentTemp * B_j;

                                            ToolFunc::matMutly<3, 6, 3>(materialTangentTemp1, Bj, materialTangent1);


                                            double geometricTangent =
                                                    geometricTangentTemp.transpose() * physGrad.col(k);

                                            for (short d = 0; d < dim; ++d)
                                                materialTangent1[d][d] += geometricTangent;

                                            for (short di = 0; di < dim; ++di)
                                                for (short dj = 0; dj < dim; ++dj)
                                                    localMats(di * N_D + j, dj * N_D + k) +=
                                                            weight * materialTangent1[di][dj];
                                            //localMats(di * N_D + j, dj * N_D + k) += weight * materialTangent(di, dj);
                                        }
                                    }

                                    double localResiduals[3] = {0};
                                    for (int d = 0; d < 6; ++d) {
                                        localResiduals[0] += Bi[0][d] * Svec(d);
                                        localResiduals[1] += Bi[1][d] * Svec(d);
                                        localResiduals[2] += Bi[2][d] * Svec(d);
                                    }
                                    //localResidual = B_i.transpose() * Svec;
                                    for (short d = 0; d < dim; d++)
                                        localRhs(d * N_D + j) -= weight * localResiduals[d];
                                }
                                for (short d = 0; d < dim; ++d)
                                    localRhs.middleRows(d * N_D, N_D).noalias() +=
                                            weight * forceScaling * forceValues(d, q) * basisValue[0].col(q);
                            }
                        }
                    }
                    // 单元刚度矩阵组装总刚度矩阵
#pragma omp critical
                    localMatsVec[i] = localMats;
                    m_system.pushToSystem(i, localIndex, localRhs, localMats, isAssembleMatrix);
                }
            }
        }
    }
}

void stiffAssembler::test_stiffNonlinearMatrix(const Eigen::MatrixXd &displacement,
                                               bool isAssembleMatrix,
                                               std::vector<Eigen::MatrixXd> &localMatsVec,
                                               std::vector<Eigen::MatrixXd> &localRhsVec) {
    localMatsVec.resize(geometry.size());
    localRhsVec.resize(geometry.size());

    //遍历所有面片
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim);
    //C;
    double CD[6][6];
    ToolFunc::setC<6>(CD, identity, lambda, mu);
#pragma omp parallel for num_threads(24)
    for (int i = 0; i < geometry.size(); ++i) {
        int N_D;
        std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
        std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
        // 某一片模型
        CBSplineSolid *geo = &geometry[i];
        int num = geometry[i].getSizeControlPoints();
        Eigen::MatrixXd blockU = displacement.block(controlPointsSumNum[i] - num, 0, num, 3);

        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        //std::vector<Eigen::MatrixXd> &basisValue = basisValueDer;
        std::vector<Eigen::MatrixXd> basisValue(2);
        std::vector<int> localIndex(gaussNum);
        Eigen::MatrixXd forceValues;
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
        {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
            {
                //高斯点选择
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                    int cnt =
                            w_interval * (vknots.size() - 1) * (uknots.size() - 1) + v_interval * (uknots.size() - 1) +
                            u_interval;
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //局部索引号
                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    //基函数
                    if (!isPreCompute)
                        geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                    else
                        basisValue = basisValueDer[cnt];
                    //TODO
                    bodyforce.eval(p_u, p_v, p_w, forceValues);
                    Eigen::MatrixXd localMats, localRhs, jacGeo, jacDis, F, RCG, E, S;
                    // materialTangent;
                    Eigen::VectorXd Svec, geometricTangentTemp;
                    // localResidual;
                    Eigen::MatrixXd physGrad, phyDispJac;
                    // materialTangentTemp, B_i, B_j;
                    double Bi[3][6], Bj[6][3], materialTangentTemp1[3][6], materialTangent1[3][3];
                    double Bi_L[3][3], Bi_R[3][3], Bj_L[3][3], Bj_R[3][3], materialTangentTemp1L[3][3], materialTangentTemp1R[3][3];
                    if (isAssembleMatrix)
                        localMats.setZero(dim * N_D, dim * N_D);
                    localRhs.setZero(dim * N_D, 1);
                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                //geo->getJacobiU(p_u[ii], p_v[jj], p_w[kk], localIndex,blockU, jacGeo, jacDis);  // 几何雅可比X/ξ  位移雅可比 U/ξ
                                geo->getJacobiU(basisValue[1].col(q), localIndex, blockU, jacGeo, jacDis);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jacGeo.determinant()); //权值
                                ToolFunc::transformGradients(jacGeo, q, basisValue[1], physGrad, dim);
                                phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());
                                F = identity + phyDispJac;   // F = I + du/dx
                                double J = F.determinant();

                                RCG = F.transpose() * F;
                                E = 0.5 * (RCG - identity);

                                if (m_materLaw.first == material_law::saint_venant_kirchhoff) {
                                    S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
                                    ToolFunc::voigtStressToK(S, Svec);                              // 应力张量
                                } else {    // 其它材料模型

                                }
                                // 局部矩阵计算 装配
                                for (size_t j = 0; j < N_D; j++) {
                                    //B = BL + BN 即小位移矩阵+大位移矩阵
                                    //ToolFunc::setB(B_i, F, physGrad.col(j));
                                    ToolFunc::setBT(Bi, F, physGrad.col(j));

                                    if (isAssembleMatrix) {
                                        //materialTangentTemp = B_i.transpose() * C;
                                        ToolFunc::matMutly<3, 6, 6>(Bi, CD, materialTangentTemp1);

                                        geometricTangentTemp = S * physGrad.col(j);
                                        //loop over active basis functions (v_j)
                                        for (size_t k = 0; k < N_D; ++k) {
                                            //ToolFunc::setB(B_j, F, physGrad.col(k));
                                            ToolFunc::setB(Bj, F, physGrad.col(k));

                                            //materialTangent = materialTangentTemp * B_j;

                                            ToolFunc::matMutly<3, 6, 3>(materialTangentTemp1, Bj, materialTangent1);


                                            double geometricTangent =
                                                    geometricTangentTemp.transpose() * physGrad.col(k);

                                            for (short d = 0; d < dim; ++d)
                                                materialTangent1[d][d] += geometricTangent;

                                            for (short di = 0; di < dim; ++di)
                                                for (short dj = 0; dj < dim; ++dj)
                                                    localMats(di * N_D + j, dj * N_D + k) +=
                                                            weight * materialTangent1[di][dj];
                                            //localMats(di * N_D + j, dj * N_D + k) += weight * materialTangent(di, dj);
                                        }
                                    }

                                    double localResiduals[3] = {0};
                                    for (int d = 0; d < 6; ++d) {
                                        localResiduals[0] += Bi[0][d] * Svec(d);
                                        localResiduals[1] += Bi[1][d] * Svec(d);
                                        localResiduals[2] += Bi[2][d] * Svec(d);
                                    }
                                    //localResidual = B_i.transpose() * Svec;
                                    for (short d = 0; d < dim; d++)
                                        localRhs(d * N_D + j) -= weight * localResiduals[d];
                                }
                                for (short d = 0; d < dim; ++d)
                                    localRhs.middleRows(d * N_D, N_D).noalias() +=
                                            weight * forceScaling * forceValues(d, q) * basisValue[0].col(q);
                            }
                        }
                    }
                    // 单元刚度矩阵组装总刚度矩阵
#pragma omp critical
                    localMatsVec[i] = localMats;
                    localRhsVec[i] = localRhs;
                    m_system.pushToSystem(i, localIndex, localRhs, localMats, isAssembleMatrix);
                }
            }
        }
    }
}

void stiffAssembler::rhsVector() {
    std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
    std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
    int N_D;
    auto it = neumanSides.begin();
    for (; it != neumanSides.end(); ++it) {
        int patchIndex = it->patchSide();
        boxSide side = it->side();
        CBSplineSolid *geo = &geometry[patchIndex];
        std::vector<double> uknots = geo->getUMKnots(),
                vknots = geo->getVMKnots(), wknots = geo->getWMKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;

        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_w.resize(gaussNum_w);
        //确定单元范围
        int u_interval_min = 1, u_interval_max = uknots.size() - 2;
        int v_interval_min = 1, v_interval_max = vknots.size() - 2;
        int w_interval_min = 1, w_interval_max = wknots.size() - 2;
        int u_add = 0, v_add = 0, w_add = 0;
        if (side.index() == boundary::west) {
            u_interval_min = 0;
            u_interval_max = 1;
            u_add = 1;
        } else if (side.index() == boundary::east) {
            u_interval_min = u_interval_max;
            u_interval_max += 1;
            u_add = -1;
        } else if (side.index() == boundary::south) {
            v_interval_min = 0;
            v_interval_max = 1;
            v_add = 1;
        } else if (side.index() == boundary::north) {
            v_interval_min = v_interval_max;
            v_interval_max += 1;
            v_add = -1;
        } else if (side.index() == boundary::front) {
            w_interval_min = 0;
            w_interval_max = 1;
            w_add = 1;
        } else {
            w_interval_min = w_interval_max;
            w_interval_max += 1;
            w_add = -1;
        }

        std::vector<int> localIndex(gaussNum);

        for (int w_interval = w_interval_min; w_interval < w_interval_max; ++w_interval) {
            ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
            for (int v_interval = v_interval_min; v_interval < v_interval_max; ++v_interval) {
                ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                for (int u_interval = u_interval_min; u_interval < u_interval_max; ++u_interval) {
                    ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                    //高斯点选择
                    geo->active(u_interval - 1 + u_add, v_interval - 1 + v_add, w_interval - 1 + w_add, localIndex);
                    N_D = localIndex.size();
                    Eigen::MatrixXd neumanValues, basisValue, localRhs;
                    it->eval(p_u, p_v, p_w, neumanValues);
                    geo->eval_into(p_u, p_v, p_w, basisValue);
                    localRhs.setZero(dim * N_D, 1);
                    for (int kk = 0; kk < gaussNum_w; ++kk) {
                        for (int jj = 0; jj < gaussNum_v; ++jj) {
                            for (int ii = 0; ii < gaussNum_u; ++ii) {
                                int q = kk * gaussNum_v * gaussNum_u + jj * gaussNum_u + ii;
                                Eigen::MatrixXd jac;
                                geo->getJacobi(p_u[ii], p_v[jj], p_w[kk], localIndex, jac);
                                Eigen::VectorXd unormal;
                                ToolFunc::outerNormal(jac, side, unormal);
                                double weight = w_u[ii] * w_v[jj] * w_w[kk] * unormal.norm();
                                for (short d = 0; d < dim; ++d)
                                    localRhs.middleRows(d * N_D, N_D).noalias() +=
                                            weight * neumanValues(d, q) * basisValue.col(q);
                            }
                        }
                    }
                    m_system.pushToRhs(patchIndex, localIndex, localRhs);
                }
            }
        }
    }
}

void stiffAssembler::rhsVector_NeumannPoint()
{
    for(auto bd_point : neumannPoints)
    {
        assert(bd_point.m_condition_type == 1);
        uint32_t global_id = bd_point.m_cpoint_global_id;

        m_system.rhs()(global_id + m_system.numControlPoints() * 0) += bd_point.m_val.getX();
        m_system.rhs()(global_id + m_system.numControlPoints() * 1) += bd_point.m_val.getY();
        m_system.rhs()(global_id + m_system.numControlPoints() * 2) += bd_point.m_val.getZ();
    }
}

void stiffAssembler::assemble() {
    m_system.matrix().setZero();
    m_system.rhs().setZero();
    stiffMatrix();
    rhsVector();
    //if (!isCollision)
    m_system.dataPush();
    //else
    //	m_system.dataPush(extraRhs);
}

void stiffAssembler::assemble_points() {
    m_system.matrix().setZero();
    m_system.rhs().setZero();
    stiffMatrix();
    rhsVector_NeumannPoint();
    //if (!isCollision)
    m_system.dataPush_Dirichlet();
    //else
    //	m_system.dataPush(extraRhs);
}

void stiffAssembler::assemble(const Eigen::VectorXd &solution, bool isAssemblerMatrix) {
    Eigen::MatrixXd displacement;
    constructCoefs(solution, displacement);  // 构造初始解

    if (isAssemblerMatrix)
        m_system.matrix().setZero();
    m_system.rhs().setZero();
    stiffNonlinearMatrix(displacement, isAssemblerMatrix);
    rhsVector();
    if (!isCollision)
        m_system.dataPush_Dirichlet();
    else
        m_system.dataPush(extraRhs);
}

void stiffAssembler::assemble_points(const Eigen::VectorXd &solution, bool isAssemblerMatrix) {
    Eigen::MatrixXd displacement;
    constructCoefs(solution, displacement);  // 构造初始解

    if (isAssemblerMatrix)
        m_system.matrix().setZero();
    m_system.rhs().setZero();
    stiffNonlinearMatrix(displacement, isAssemblerMatrix);
    rhsVector_NeumannPoint();
    if (!isCollision)
        m_system.dataPush_Dirichlet();
    else
        m_system.dataPush(extraRhs);
}

void stiffAssembler::assembleNoRhs(const Eigen::VectorXd &solution, std::vector<pair<int, double>> f,
                                   bool isAssemblerMatrix) {
    Eigen::MatrixXd displacement;
    constructCoefs(solution, displacement);  // 构造初始解

    if (isAssemblerMatrix)
        m_system.matrix().setZero();
    m_system.rhs().setZero();
    stiffNonlinearMatrix(displacement, isAssemblerMatrix);
    m_system.pushToRhs(f);
    m_system.dataPush();
}

void stiffAssembler::setMaterialLaw(const pair<int, vector<double>> material) {
    m_materLaw = material;
}

void stiffAssembler::setYangsPossion(double Yangs, double Possion)
{
    E = Yangs;
    pr = Possion;

    lambda = E * pr / ((1. + pr) * (1. - 2. * pr));    //拉梅常数
    mu = E / (2. * (1. + pr));
}

void stiffAssembler::updateCoefs(Eigen::MatrixXd &coefs) {
    int begin = 0;
    for (size_t i = 0; i < geometry.size(); ++i) {
        int num = geometry[i].getSizeControlPoints();
        geometry[i].updateCoefs(coefs.block(begin, 0, num, 3));
        begin += num;
    }
}

void stiffAssembler::constructCoefs(const Eigen::VectorXd &solution, Eigen::MatrixXd &displacement) {
    if (solution.rows() != dim * m_system.getControlNums()) {
        cout << "solution length is not Equal to the size of contropoints" << endl;
        return;
    }
    m_system.constructCoefs(solution, displacement);
}

const Eigen::SparseMatrix<double> &stiffAssembler::getMatrix() const {
    return m_system.matrix();
}

const Eigen::MatrixXd &stiffAssembler::getRhs() const {
    return m_system.rhs();
}

void stiffAssembler::preComputeBasis() {
    int N_D;
    //遍历所有面片  // 某一片模型
    std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
    std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
    CBSplineSolid *geo = &geometry[0];
    std::vector<double> uknots = geo->getUKnots(),
            vknots = geo->getVKnots(), wknots = geo->getWKnots();
    int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
            gaussNum_w = geo->getWDegree() + 1;
    int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
    p_u.resize(gaussNum_u);
    p_v.resize(gaussNum_v);
    p_v.resize(gaussNum_w);
    std::vector<Eigen::MatrixXd> basisValue(2);
    std::vector<int> localIndex(gaussNum);
    //ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[0], wknots[1], p_w, w_w);
    //ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[0], vknots[1], p_v, w_v);
    //ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[0], uknots[1], p_u, w_u);
    //基函数
    for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
    {
        ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
        for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
        {
            //高斯点选择
            ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
            for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                geo->eval_AllDers(p_u, p_v, p_w, basisValue);
                basisValueDer.emplace_back(basisValue);
            }
        }
    }
}

bool stiffAssembler::isCollide(const Eigen::MatrixXd &dispVector, Eigen::SparseMatrix<double> &newMatrix,
                               Eigen::MatrixXd &extRhs) {
    Eigen::MatrixXd displacement;
    constructCoefs(dispVector, displacement);
    updateCoefs(displacement);           // 更新位置

    auto inter = getInterface();
    auto connect = getConnectPatch();
    CollisionDetect detect(connect, inter);
    std::vector<std::vector<CollideInfo>> collidInfo;   // 记录每一边界面的碰撞情况
    std::vector<pair<int, int>> patchIds;              // 记录边界面信息
    bool flag = detect.ObjectDetect(geometry, collidInfo, patchIds);   // 碰撞检测
    displacement *= -1;
    updateCoefs(displacement);                    // 复原位置

    if (!flag) {
        extraRhs.setZero(m_system.rhs().rows(), 1);
        isCollision = false;
        return false;
    }
    std::cout << "Yes" << std::endl;

    //isCollision = true;
    Eigen::VectorXd newRhs;
    newRhs.setZero(m_system.rhs().rows());

    for (int i = 0; i < collidInfo.size(); ++i) {
        if (!collidInfo[i].size()) continue;
        int patchId = patchIds[i].first;
        int faceId = patchIds[i].second;
        int contropoints = m_system.numControlPoints();

        for (int j = 0; j < collidInfo[i].size(); ++j) {

            CPoint para = collidInfo[i][j].src;
            std::vector<int> localIndex;
            Eigen::VectorXd basisValue;
            geometry[patchId].eval(para.getX(), para.getY(), para.getZ(), basisValue);
            localIndex.resize(basisValue.rows());
            geometry[patchId].activeFromParam(para.getX(), para.getY(), para.getZ(), localIndex);


            double sum = 0;
            for (int ids = 0; ids < localIndex.size(); ++ids) {
                if (geometry[patchId].isBelongToTheFace(faceId, localIndex[ids]) &&
                    !m_system.IsDirichlet(patchId, localIndex[ids]))
                    sum += basisValue(ids);
            }

            for (int ids = 0; ids < localIndex.size(); ++ids) {
                if (geometry[patchId].isBelongToTheFace(faceId, localIndex[ids]) &&
                    !m_system.IsDirichlet(patchId, localIndex[ids])) {
                    CPoint disp = collidInfo[i][j].move * (basisValue(ids) / sum);
                    int globalLocal = m_system.getIndexById(patchId, localIndex[ids]);

                    newRhs(globalLocal) = fabs(newRhs(globalLocal)) < fabs(disp.getX()) ?
                                          disp.getX() : newRhs(globalLocal);
                    //else newRhs(globalLocal) += disp.getX();

                    newRhs(globalLocal + contropoints) =
                            fabs(newRhs(globalLocal + contropoints)) < fabs(disp.getY()) ? disp.getY() : newRhs(
                                    globalLocal + contropoints);
                    //else newRhs(globalLocal + contropoints) += disp.getY();

                    newRhs(globalLocal + 2 * contropoints) =
                            fabs(newRhs(globalLocal + 2 * contropoints)) < fabs(disp.getZ()) ? disp.getZ() : newRhs(
                                    globalLocal + 2 * contropoints);
                    //else newRhs(globalLocal + 2 * contropoints) += disp.getZ();

                }
            }
        }
    }

    for (int i = 0; i < newRhs.rows(); ++i) {
        if (fabs(newRhs(i)) > 1e-8) {
            newRhs(i) = dispVector(i, 0) + newRhs(i);
        }
    }
    extraRhs = newRhs;
    //std::cout << extraRhs << "\n";
    m_system.constructNewMatrix(newRhs, newMatrix, extRhs);

    return true;
}


void stiffAssembler::WriteParaview(std::vector<CBSplineSolid> &patch, const std::string &fn, Eigen::MatrixXd &solution,
                                   unsigned npts, bool isStress) {
    const unsigned numP = patch.size();
    Eigen::MatrixXd displacement;
    constructCoefs(solution, displacement);
    std::vector<CBSplineSolid> solVector;
    int begin = 0;
    for (size_t i = 0; i < numP; ++i) {
        solVector.push_back(geometry[i]);
        int num = geometry[i].getSizeControlPoints();
        solVector[i].setCoefs(displacement.block(begin, 0, num, 3));
        begin += num;
    }


    ParaviewCollection collection(fn);
    std::string fileName = fn.substr(fn.find_last_of("/\\") + 1);

    for (unsigned i = 0; i < numP; ++i) {
        //cout << "output" << " " << i << " " << "all " << numP << endl;

        WriteSingleParaview(i, fn + to_string(i), npts, patch[i], solVector[i], displacement, isStress);

        collection.addPart(fileName + to_string(i), ".vts");

    }
    collection.save();


}

void stiffAssembler::WriteSingleParaview(const unsigned patchNum, const std::string &fn, unsigned npts,
                                         CBSplineSolid &displacement, CBSplineSolid &dispVector,
                                         Eigen::MatrixXd &solution, bool isStress) {

    const CBSplineSolid &geo = displacement;
    const short n = dim;
    const short d = dim;

    Eigen::MatrixXd ab = geo.support();
    Eigen::VectorXd a = ab.col(0);
    Eigen::VectorXd b = ab.col(1);

    auto np = ToolFunc::uniformSampleCount(a, b, npts);
    int prod = 1;
    for (int i = 0; i < np.size(); ++i) prod *= np[i];
    Eigen::MatrixXd pts(a.size(), prod);

    Eigen::VectorXd h = b - a;

    for (int i = 0; i < np.size(); ++i) h(i) /= (np[i] - 1);
    unsigned cnt = 0;
    for (int k = 0; k < np[2]; ++k) {
        for (int j = 0; j < np[1]; ++j) {
            for (int i = 0; i < np[0]; ++i) {
                pts(0, cnt) = min(a(0) + h(0) * i, b(0));
                pts(1, cnt) = min(a(1) + h(1) * j, b(1));
                pts(2, cnt++) = min(a(2) + h(2) * k, b(2));
            }
        }
    }
    Eigen::MatrixXd eval_disp = dispVector.eval(pts);
    Eigen::MatrixXd eval_geo = displacement.eval(pts);

    std::map<std::string, Eigen::MatrixXd> data;

    data["Deformation"] = eval_disp;
    //auto stress = eval_stress(pts, solution, patchNum);
    if (isStress)
        data["Stress"] = eval_stress(pts, solution, patchNum);
    //data["Stress"] = stress;
    //cout << eval_disp.rows() << " " << eval_disp.cols() << endl;
    //cout << data["Stress"].rows() << " " << data["Stress"].cols() << endl;
    OutputParaview::ParaviewMulti(eval_geo, data, np, fn);
}

Eigen::MatrixXd
stiffAssembler::eval_stress(const Eigen::MatrixXd &u, const Eigen::MatrixXd &displacement, int patchId) {
    auto &geo = geometry[patchId];
    int gaussNum = (geo.getUDegree() + 1) * (geo.getVDegree() + 1) * (geo.getWDegree() + 1);
    Eigen::MatrixXd jacGeo, jacDis, dispGrad;
    int num = geo.getSizeControlPoints();
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(dim, dim);
    Eigen::MatrixXd blockU = displacement.block(controlPointsSumNum[patchId] - num, 0, num, 3);
    Eigen::MatrixXd result;
    result.setZero(3, u.cols());
    for (int q = 0; q < u.cols(); ++q) {
        std::vector<Eigen::VectorXd> basisValue;
        std::vector<int> localIndex(gaussNum);
        geo.active(0, 0, 0, localIndex);
        int N_D = localIndex.size();
        //geo.eval_2Ders(u(0,q), u(1,q), u(2,q), basisValue);
        geo.getJacobiU(u(0, q), u(1, q), u(2, q), localIndex, blockU, jacGeo, jacDis);
        dispGrad = jacDis.transpose() * (jacGeo.transpose().inverse());
        //cout << jacGeo << " \n" << jacDis;

        double s11 = lambda * (dispGrad(0, 0) + dispGrad(1, 1) + dispGrad(2, 2)) + 2 * mu * dispGrad(0, 0);
        double s22 = lambda * (dispGrad(0, 0) + dispGrad(1, 1) + dispGrad(2, 2)) + 2 * mu * dispGrad(1, 1);
        double s33 = lambda * (dispGrad(0, 0) + dispGrad(1, 1) + dispGrad(2, 2)) + 2 * mu * dispGrad(2, 2);
        double s12 = mu * (dispGrad(0, 1) + dispGrad(1, 0));
        double s13 = mu * (dispGrad(0, 2) + dispGrad(2, 0));
        double s23 = mu * (dispGrad(1, 2) + dispGrad(2, 1));
        result(0, q) = sqrt(0.5 * ((s11 - s22) * (s11 - s22) + (s11 - s33) * (s11 - s33) + (s22 - s33) * (s22 - s33) +
                                   6 * (s12 * s12 + s13 * s13 + s23 * s23)));

        //F = identity + phyDispJac;   // F = I + du/dx

        //RCG = F.transpose() * F;
        //E = 0.5 * (RCG - identity);
        //S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
        //if (patchId == 11) {
        //		cout << "\n\n";
        //		cout << dispGrad << "\n";
        //		cout << s11 << " " << s22 << " " << s33 << endl;
        //
        //}
        //result(0, q) = sqrt(0.5 * ((S(0, 0) - S(1, 1)) * (S(0, 0) - S(1, 1)) + (S(0, 0) - S(2, 2)) * (S(0, 0) - S(2, 2)) +
        //	(S(1, 1) - S(2, 2)) * (S(1, 1) - S(2, 2)) + 6 * (S(0, 1) * S(0, 1) + S(0, 2) * S(0, 2) + S(1, 2) * S(1, 2))));
    }
    return result;
}

