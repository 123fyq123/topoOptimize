#include "massAssembler.h"

namespace Poster_Needed {
    massAssembler::massAssembler(std::vector<Poster_Needed::CBSplineSolid> _geometry, double _density,
                                 sparseSystem _system) {
        geometry = _geometry;
        density = _density;
        m_system = _system;
    }

    void massAssembler::assemble() {
        m_system.matrix().setZero();
        //遍历所有面片
        int dim = m_system.getDim();
#pragma omp parallel for num_threads(8)
        for (int i = 0; i < geometry.size(); ++i) {
            // 某一片模型
            Poster_Needed::CBSplineSolid *geo = &geometry[i];
            std::vector<double> uknots = geo->getUKnots(),
                    vknots = geo->getVKnots(), wknots = geo->getWKnots();
            int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                    gaussNum_w = geo->getWDegree() + 1;
            int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
            std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
            std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值
            int N_D;
            p_u.resize(gaussNum_u);
            p_v.resize(gaussNum_v);
            p_v.resize(gaussNum_w);
            Eigen::MatrixXd basisValue;
            std::vector<int> localIndex(gaussNum);
            for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // 单元u方向循环
            {
                ToolFunc::gaussIntegalPoint(gaussNum_w, wknots[w_interval], wknots[w_interval + 1], p_w, w_w);
                for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // 单元v方向循环
                {
                    //高斯点选择
                    ToolFunc::gaussIntegalPoint(gaussNum_v, vknots[v_interval], vknots[v_interval + 1], p_v, w_v);
                    for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval) {
                        ToolFunc::gaussIntegalPoint(gaussNum_u, uknots[u_interval], uknots[u_interval + 1], p_u, w_u);
                        //局部索引号
                        geo->active(u_interval, v_interval, w_interval, localIndex);
                        N_D = localIndex.size();
                        //基函数
                        geo->eval_into(p_u, p_v, p_w, basisValue);
                        Eigen::MatrixXd localMats, block;
                        localMats.setZero(dim * N_D, dim * N_D);
                        for (int kk = 0; kk < gaussNum_w; ++kk) {
                            for (int jj = 0; jj < gaussNum_v; ++jj) {
                                for (int ii = 0; ii < gaussNum_u; ++ii) {
                                    int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                    Eigen::MatrixXd jac, P;
                                    geo->getJacobi(p_u[ii], p_v[jj], p_w[kk], localIndex, jac);
                                    double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant()) * density;
                                    block = weight * basisValue.col(q) * basisValue.col(q).transpose();
                                    for (short d = 0; d < dim; ++d)
                                        localMats.block(d * N_D, d * N_D, N_D, N_D) += block.block(0, 0, N_D, N_D);

                                }
                            }
                        }
#pragma omp critical
                        m_system.pushToMatrix(i, localIndex, localMats);
                    }
                }
            }
        }
        m_system.dataPush();
        m_matrix = m_system.matrix();
        m_matrix.makeCompressed();
    }

    const Eigen::SparseMatrix<double> massAssembler::matrix() const {
        return m_matrix;
    }
}