#include "GP_IGA_CPU_Kernel.h"
#include <tbb/parallel_for.h>
#include <tbb/parallel_reduce.h>

template<typename T>
void linear_stiff_assembly(const vector<CBSplineSolid> &patch, const gpmesh::AssemblyContextCPU<T> &assembly_context,
                           const double lambda, const double mu) {
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(3, 3);
    Eigen::MatrixXd C;
    ToolFunc::setC(C, identity, lambda, mu);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, patch.size()),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); ++i) {

                              T *localMats = assembly_context.get_assembly_info(i).local_mats;

                              std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
                              std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值

                              const CBSplineSolid &geo = patch[i];

                              std::vector<double> uknots = geo.getUKnots(),
                                      vknots = geo.getVKnots(), wknots = geo.getWKnots();
                              int gaussNum_u = geo.getUDegree() + 1, gaussNum_v = geo.getVDegree() + 1,
                                      gaussNum_w = geo.getWDegree() + 1;
                              int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;

                              p_u.resize(gaussNum_u);
                              p_v.resize(gaussNum_v);
                              p_v.resize(gaussNum_w);

                              std::vector<Eigen::MatrixXd> basisValue(2);
                              std::vector<int> localIndex(gaussNum);
                              size_t N_D = localIndex.size();

                              ToolFunc::gaussIntegalPoint(gaussNum_w, 0, 1, p_w, w_w);
                              ToolFunc::gaussIntegalPoint(gaussNum_v, 0, 1, p_v, w_v);
                              ToolFunc::gaussIntegalPoint(gaussNum_u, 0, 1, p_u, w_u);

                              geo.active(0, 0, 0, localIndex);
                              geo.eval_AllDers(p_u, p_v, p_w, basisValue);

                              for (int kk = 0; kk < gaussNum_w; ++kk) {
                                  for (int jj = 0; jj < gaussNum_v; ++jj) {
                                      for (int ii = 0; ii < gaussNum_u; ++ii) {
                                          int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;
                                          Eigen::MatrixXd jac, P;
                                          geo.getJacobi(basisValue[1].col(q), localIndex, jac);
                                          double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jac.determinant());

                                          Eigen::MatrixXd physGrad;
                                          Eigen::MatrixXd B_i, B_j, tempK, K;
                                          ToolFunc::transformGradients(jac, q, basisValue[1], physGrad, 3);

                                          double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;

                                          for (size_t j = 0; j < N_D; ++j) {
                                              m_Rx = physGrad.col(j)(0);
                                              m_Ry = physGrad.col(j)(1);
                                              m_Rz = physGrad.col(j)(2);

                                              ToolFunc::setB(B_i, identity, physGrad.col(j));
                                              tempK = B_i.transpose() * C;

                                              for (size_t k = 0; k < N_D; ++k) {
                                                  ToolFunc::setB(B_j, identity, physGrad.col(k));
                                                  m_Qx = physGrad.col(k)(0);
                                                  m_Qy = physGrad.col(k)(1);
                                                  m_Qz = physGrad.col(k)(2);

                                                  localMats[(j + 0 * N_D) + 192 * (k + 0 * N_D)] += weight *
                                                                                                    (C(0, 0) * m_Rx *
                                                                                                     m_Qx + C(3, 3) *
                                                                                                            (m_Ry *
                                                                                                             m_Qy +
                                                                                                             m_Rz *
                                                                                                             m_Qz));
                                                  localMats[(j + 0 * N_D) + 192 * (k + 1 * N_D)] +=
                                                          weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
                                                  localMats[(j + 0 * N_D) + 192 * (k + 2 * N_D)] +=
                                                          weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx);
                                                  localMats[(j + 1 * N_D) + 192 * (k + 0 * N_D)] +=
                                                          weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
                                                  localMats[(j + 1 * N_D) + 192 * (k + 1 * N_D)] += weight *
                                                                                                    (C(0, 0) * m_Ry *
                                                                                                     m_Qy + C(3, 3) *
                                                                                                            (m_Rx *
                                                                                                             m_Qx +
                                                                                                             m_Rz *
                                                                                                             m_Qz));
                                                  localMats[(j + 1 * N_D) + 192 * (k + 2 * N_D)] +=
                                                          weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
                                                  localMats[(j + 2 * N_D) + 192 * (k + 0 * N_D)] +=
                                                          weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
                                                  localMats[(j + 2 * N_D) + 192 * (k + 1 * N_D)] +=
                                                          weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
                                                  localMats[(j + 2 * N_D) + 192 * (k + 2 * N_D)] += weight *
                                                                                                    (C(0, 0) * m_Rz *
                                                                                                     m_Qz + C(3, 3) *
                                                                                                            (m_Ry *
                                                                                                             m_Qy +
                                                                                                             m_Rx *
                                                                                                             m_Qx));
                                              }
                                          }
                                      }
                                  }
                              }
                          }
                      });
}

template void linear_stiff_assembly<float>(const vector<CBSplineSolid> &patch,
                                           const gpmesh::AssemblyContextCPU<float> &assembly_context,
                                           const double lambda, const double mu);

template void linear_stiff_assembly<double>(const vector<CBSplineSolid> &patch,
                                            const gpmesh::AssemblyContextCPU<double> &assembly_context,
                                            const double lambda, const double mu);


template<typename T>
void hyperelasticity_stiff_assembly(const vector<CBSplineSolid> &patch,
                                    const gpmesh::AssemblyContextCPU<T> &assembly_context,
                                    const gpmesh::PhysicsContextCPU<T> &physics_context,
                                    const double lambda, const double mu) {
    Eigen::MatrixXd identity = Eigen::MatrixXd::Identity(3, 3);
    double CD[6][6];
    ToolFunc::setC<6>(CD, identity, lambda, mu);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, patch.size()),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t i = r.begin(); i < r.end(); ++i) {

                              T *localMats = assembly_context.get_assembly_info(i).local_mats;

                              std::vector<double> p_u, p_v, p_w; // 存储单方向高斯点
                              std::vector<double> w_u, w_v, w_w; // 存储单方向高斯点权值

                              const CBSplineSolid &geo = patch[i];

                              Eigen::Matrix<T, 64, 3> &blockU = physics_context.get_physics_info(
                                      i).m_cpoints_displacement;

                              std::vector<double> uknots = geo.getUKnots(),
                                      vknots = geo.getVKnots(), wknots = geo.getWKnots();
                              int gaussNum_u = geo.getUDegree() + 1, gaussNum_v = geo.getVDegree() + 1,
                                      gaussNum_w = geo.getWDegree() + 1;
                              int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;

                              p_u.resize(gaussNum_u);
                              p_v.resize(gaussNum_v);
                              p_v.resize(gaussNum_w);

                              std::vector<Eigen::MatrixXd> basisValue(2);
                              std::vector<int> localIndex(gaussNum);
                              size_t N_D = localIndex.size();

                              ToolFunc::gaussIntegalPoint(gaussNum_w, 0, 1, p_w, w_w);
                              ToolFunc::gaussIntegalPoint(gaussNum_v, 0, 1, p_v, w_v);
                              ToolFunc::gaussIntegalPoint(gaussNum_u, 0, 1, p_u, w_u);

                              geo.active(0, 0, 0, localIndex);
                              geo.eval_AllDers(p_u, p_v, p_w, basisValue);

                              Eigen::MatrixXd jacGeo, jacDis, F, RCG, E, S;
                              Eigen::VectorXd Svec, geometricTangentTemp;
                              Eigen::MatrixXd physGrad, phyDispJac;
                              double Bi[3][6], Bj[6][3], materialTangentTemp1[3][6], materialTangent1[3][3];
                              double Bi_L[3][3], Bi_R[3][3], Bj_L[3][3], Bj_R[3][3], materialTangentTemp1L[3][3], materialTangentTemp1R[3][3];

                              for (int kk = 0; kk < gaussNum_w; ++kk) {
                                  for (int jj = 0; jj < gaussNum_v; ++jj) {
                                      for (int ii = 0; ii < gaussNum_u; ++ii) {
                                          int q = kk * gaussNum_u * gaussNum_v + jj * gaussNum_u + ii;


                                          geo.getJacobiU(basisValue[1].col(q), localIndex,
                                                         blockU.template cast<double>(), jacGeo, jacDis);
                                          double weight = w_u[ii] * w_v[jj] * w_w[kk] * fabs(jacGeo.determinant()); //权值
                                          ToolFunc::transformGradients(jacGeo, q, basisValue[1], physGrad, 3);
                                          phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());
                                          F = identity + phyDispJac;   // F = I + du/dx
                                          double J = F.determinant();

                                          RCG = F.transpose() * F;
                                          E = 0.5 * (RCG - identity);

                                          S = lambda * E.trace() * identity + 2 * mu * E;       // 应力
                                          ToolFunc::voigtStressToK(S, Svec);                              // 应力张量

                                          for (size_t j = 0; j < N_D; ++j) {

                                              ToolFunc::setBT(Bi, F, physGrad.col(j));

                                              ToolFunc::matMutly<3, 6, 6>(Bi, CD, materialTangentTemp1);

                                              geometricTangentTemp = S * physGrad.col(j);

                                              for (size_t k = 0; k < N_D; ++k) {

                                                  ToolFunc::setB(Bj, F, physGrad.col(k));

                                                  ToolFunc::matMutly<3, 6, 3>(materialTangentTemp1, Bj,
                                                                              materialTangent1);

                                                  double geometricTangent =
                                                          geometricTangentTemp.transpose() * physGrad.col(k);

                                                  for (short d = 0; d < 3; ++d)
                                                      materialTangent1[d][d] += geometricTangent;

                                                  for (size_t di = 0; di < 3; ++di) {
                                                      for (size_t dj = 0; dj < 3; ++dj) {
                                                          localMats[(j + di * N_D) + 192 * (k + dj * N_D)] +=
                                                                  weight * materialTangent1[di][dj];
                                                      }
                                                  }
                                              }
                                          }
                                      }
                                  }
                              }
                          }
                      });
}


template void hyperelasticity_stiff_assembly<double>(const vector<CBSplineSolid> &patch,
                                                     const gpmesh::AssemblyContextCPU<double> &assembly_context,
                                                     const gpmesh::PhysicsContextCPU<double> &physics_context,
                                                     const double lambda, const double mu);

template void hyperelasticity_stiff_assembly<float>(const vector<CBSplineSolid> &patch,
                                                    const gpmesh::AssemblyContextCPU<float> &assembly_context,
                                                    const gpmesh::PhysicsContextCPU<float> &physics_context,
                                                    const double lambda, const double mu);

template<typename T>
void add_dirichlet_point_info_to_assembly(const gpmesh::AssemblyContextCPU<T> &assembly_context,
                                          const gpmesh::DirichletPointContextCPU<T> &dirichlet_point_context) {
    tbb::parallel_for(tbb::blocked_range<size_t>(0, dirichlet_point_context.get_num_dirichlet_point()),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t id = r.begin(); id < r.end(); ++id) {
                              auto &dirichlet_point_info = dirichlet_point_context.get_dirichlet_point_info(id);

                              uint32_t geometry_id = dirichlet_point_info.geometry_id;
                              uint16_t local_id = dirichlet_point_info.local_id;


                              auto &geometry_info = assembly_context.get_assembly_info(geometry_id);
                              assert(geometry_id == geometry_info.geometry_id);

                              auto local_mat = geometry_info.local_mats;

                              for (uint16_t j = 0; j < 192; ++j) {
                                  for (uint16_t d = 0; d < 3; ++d) {
                                      local_mat[(local_id + d * 64) + 192 * (j)] = 0;
                                  }
                              }

                              for (uint16_t d = 0; d < 3; ++d) {
                                  local_mat[(local_id + d * 64) + 192 * (local_id + d * 64)] += 1;
                              }
                          }
                      }
    );
}

template void add_dirichlet_point_info_to_assembly<float>(const gpmesh::AssemblyContextCPU<float> &assembly_context,
                                                          const gpmesh::DirichletPointContextCPU<float> &dirichlet_point_context);

template void add_dirichlet_point_info_to_assembly<double>(const gpmesh::AssemblyContextCPU<double> &assembly_context,
                                                           const gpmesh::DirichletPointContextCPU<double> &dirichlet_point_context);

template<typename T>
void matvec(const gpmesh::AssemblyContextCPU<T> &assembly_context,
            uint32_t num_geometry,
            std::shared_ptr<gpmesh::SplineAttributeCPU<T>> in,
            std::shared_ptr<gpmesh::SplineAttributeCPU<T>> out) {

    tbb::task_arena ta(16);

    ta.execute([&] {

        tbb::parallel_for(tbb::blocked_range<size_t>(0, num_geometry),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t geometry_id = r.begin(); geometry_id < r.end(); ++geometry_id) {
                                  //out = local_mats * in
                                  auto local_mat = assembly_context.get_assembly_info(geometry_id).local_mats;
                                  for (uint16_t i = 0; i < 192; ++i) {
                                      for (uint16_t j = 0; j < 192; ++j) {
                                          out->RawAccessLocal(geometry_id, i) +=
                                                  local_mat[i + 192 * j] * in->RawAccessLocal(geometry_id, j);
                                      }
                                      //cout << out->RawAccessLocal(geometry_id, i) << endl;
                                  }
                              }
                          }
        );

//        tbb::parallel_for(tbb::blocked_range<size_t>(0, num_geometry, num_geometry / (2*numproc)),
//                          [&](tbb::blocked_range<size_t> r) {
//                              for (size_t geometry_id = r.begin(); geometry_id < r.end(); ++geometry_id) {
//                                  //out = local_mats * in
//                                  auto local_mat = assembly_context.get_assembly_info(geometry_id).local_mats;
//                                  for (uint16_t i = 0; i < 192; ++i) {
//                                      for (uint16_t j = 0; j < 192; ++j) {
//                                          out->RawAccessLocal(geometry_id, i) +=
//                                                  local_mat[i + 192 * j] * in->RawAccessLocal(geometry_id, j);
//                                      }
//                                      //cout << out->RawAccessLocal(geometry_id, i) << endl;
//                                  }
//                              }
//                          }, tbb::simple_partitioner{}
//        );
    });

}


template void matvec<float>(const gpmesh::AssemblyContextCPU<float> &assembly_context,
                            uint32_t num_geometry,
                            std::shared_ptr<gpmesh::SplineAttributeCPU<float>> in,
                            std::shared_ptr<gpmesh::SplineAttributeCPU<float>> out);

template void matvec<double>(const gpmesh::AssemblyContextCPU<double> &assembly_context,
                             uint32_t num_geometry,
                             std::shared_ptr<gpmesh::SplineAttributeCPU<double>> in,
                             std::shared_ptr<gpmesh::SplineAttributeCPU<double>> out);


template<typename T>
void init_P_R(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> S,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> P,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> R) {
    uint32_t globalSize = B->m_num_cpoints * B->m_num_attributes;
    assert(globalSize == S->m_num_cpoints * S->m_num_attributes);
    assert(globalSize == P->m_num_cpoints * P->m_num_attributes);
    assert(globalSize == R->m_num_cpoints * R->m_num_attributes);
    tbb::parallel_for(tbb::blocked_range<size_t>(0, globalSize),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t id = r.begin(); id < r.end(); ++id) {
                              P->RawAccessGlobal(id) = B->RawAccessGlobal(id) - S->RawAccessGlobal(id);
                              R->RawAccessGlobal(id) = P->RawAccessGlobal(id);
                          }
                      }
    );
}

template void init_P_R<double>(std::shared_ptr<gpmesh::SplineAttributeCPU<double>> B,
                               std::shared_ptr<gpmesh::SplineAttributeCPU<double>> S,
                               std::shared_ptr<gpmesh::SplineAttributeCPU<double>> P,
                               std::shared_ptr<gpmesh::SplineAttributeCPU<double>> R);

template void init_P_R<float>(std::shared_ptr<gpmesh::SplineAttributeCPU<float>> B,
                              std::shared_ptr<gpmesh::SplineAttributeCPU<float>> S,
                              std::shared_ptr<gpmesh::SplineAttributeCPU<float>> P,
                              std::shared_ptr<gpmesh::SplineAttributeCPU<float>> R);

template<typename T>
void Spline_Attribute_CWise_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> res,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> other) {
    uint32_t globalSize = res->m_num_cpoints * res->m_num_attributes;
    assert(globalSize == other->m_num_cpoints * other->m_num_attributes);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, globalSize),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t id = r.begin(); id < r.end(); ++id) {
                              res->RawAccessGlobal(id) *= other->RawAccessGlobal(id);
                          }
                      }
    );
}

template void Spline_Attribute_CWise_dot<double>(std::shared_ptr<gpmesh::SplineAttributeCPU<double>> res,
                                                 std::shared_ptr<gpmesh::SplineAttributeCPU<double>> other);

template void Spline_Attribute_CWise_dot<float>(std::shared_ptr<gpmesh::SplineAttributeCPU<float>> res,
                                                std::shared_ptr<gpmesh::SplineAttributeCPU<float>> other);

template<typename T>
T Spline_Attribute_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> A,
                       std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B) {
    uint32_t globalSize = A->m_num_cpoints * A->m_num_attributes;
    assert(globalSize == B->m_num_cpoints * B->m_num_attributes);
    T res = tbb::parallel_reduce(tbb::blocked_range<size_t>(0, globalSize), static_cast<T>(0),
                                 [&](tbb::blocked_range<size_t> r, T local_res) -> T {
                                     for (size_t id = r.begin(); id < r.end(); ++id) {
                                         local_res += A->RawAccessGlobal(id) * B->RawAccessGlobal(id);
                                     }

                                     return local_res;
                                 },
                                 [](T x, T y) -> T {
                                     return x + y;
                                 }
    );

    return res;
}

template double Spline_Attribute_dot<double>(std::shared_ptr<gpmesh::SplineAttributeCPU<double>> A,
                                             std::shared_ptr<gpmesh::SplineAttributeCPU<double>> B);

template float Spline_Attribute_dot<float>(std::shared_ptr<gpmesh::SplineAttributeCPU<float>> A,
                                           std::shared_ptr<gpmesh::SplineAttributeCPU<float>> B);

//Y = alpha * X + beta * Y
template<typename T>
void Spline_Attribute_axpy(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> Y,
                           std::shared_ptr<gpmesh::SplineAttributeCPU<T>> X,
                           T alpha,
                           T beta) {
    uint32_t globalSize = X->m_num_cpoints * X->m_num_attributes;
    assert(globalSize == Y->m_num_cpoints * Y->m_num_attributes);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, globalSize),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t id = r.begin(); id < r.end(); ++id) {
                              Y->RawAccessGlobal(id) = alpha * X->RawAccessGlobal(id) + beta * Y->RawAccessGlobal(id);
                          }
                      }
    );
}

template void Spline_Attribute_axpy<double>(std::shared_ptr<gpmesh::SplineAttributeCPU<double>> Y,
                                            std::shared_ptr<gpmesh::SplineAttributeCPU<double>> X,
                                            double alpha,
                                            double beta);

template void Spline_Attribute_axpy<float>(std::shared_ptr<gpmesh::SplineAttributeCPU<float>> Y,
                                           std::shared_ptr<gpmesh::SplineAttributeCPU<float>> X,
                                           float alpha,
                                           float beta);


template<typename T>
void Spline_Attribute_CWise_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> res,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> A,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B) {
    uint32_t globalSize = res->m_num_cpoints * res->m_num_attributes;
    assert(globalSize == A->m_num_cpoints * A->m_num_attributes);
    assert(globalSize == B->m_num_cpoints * B->m_num_attributes);

    tbb::parallel_for(tbb::blocked_range<size_t>(0, globalSize),
                      [&](tbb::blocked_range<size_t> r) {
                          for (size_t id = r.begin(); id < r.end(); ++id) {
                              res->RawAccessGlobal(id) = A->RawAccessGlobal(id) * B->RawAccessGlobal(id);
                          }
                      }
    );
}

template void Spline_Attribute_CWise_dot<double>(std::shared_ptr<gpmesh::SplineAttributeCPU<double>> res,
                                                 std::shared_ptr<gpmesh::SplineAttributeCPU<double>> A,
                                                 std::shared_ptr<gpmesh::SplineAttributeCPU<double>> B);


template void Spline_Attribute_CWise_dot<float>(std::shared_ptr<gpmesh::SplineAttributeCPU<float>> res,
                                                std::shared_ptr<gpmesh::SplineAttributeCPU<float>> A,
                                                std::shared_ptr<gpmesh::SplineAttributeCPU<float>> B);


template<typename T>
T Spline_Attribute_squared_norm(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> in) {
    return Spline_Attribute_dot<T>(in, in);
}

template double Spline_Attribute_squared_norm(shared_ptr<gpmesh::SplineAttributeCPU<double>> in);

template float Spline_Attribute_squared_norm(shared_ptr<gpmesh::SplineAttributeCPU<float>> in);