#pragma once

#include <Hyperelasticity/ToolFunc.h>
#include <Hyperelasticity/CBSplineSolid.h>
#include "GP_Context.h"
#include "GP_Spline_Attribute_CPU.h"

template<typename T>
void linear_stiff_assembly(const vector<CBSplineSolid> &patch,
                           const gpmesh::AssemblyContextCPU<T> &assembly_context,
                           const double lambda, const double mu);

template<typename T>
void hyperelasticity_stiff_assembly(const vector<CBSplineSolid> &patch,
                                    const gpmesh::AssemblyContextCPU<T> &assembly_context,
                                    const gpmesh::PhysicsContextCPU<T> &physics_context,
                                    const double lambda, const double mu);

template<typename T>
void add_dirichlet_point_info_to_assembly(const gpmesh::AssemblyContextCPU<T> &assembly_context,
                                          const gpmesh::DirichletPointContextCPU<T> &dirichlet_point_context);

template<typename T>
void matvec(const gpmesh::AssemblyContextCPU<T> &assembly_context,
            uint32_t num_geometry,
            std::shared_ptr<gpmesh::SplineAttributeCPU<T>> in,
            std::shared_ptr<gpmesh::SplineAttributeCPU<T>> out);

template<typename T>
void init_P_R(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> S,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> P,
              std::shared_ptr<gpmesh::SplineAttributeCPU<T>> R);

template<typename T>
void Spline_Attribute_CWise_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> res,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> other);

template<typename T>
T Spline_Attribute_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> A,
                       std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B);

template<typename T>
void Spline_Attribute_axpy(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> Y,
                           std::shared_ptr<gpmesh::SplineAttributeCPU<T>> X,
                           T alpha,
                           T beta);

template<typename T>
void Spline_Attribute_CWise_dot(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> res,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> A,
                                std::shared_ptr<gpmesh::SplineAttributeCPU<T>> B);

template<typename T>
T Spline_Attribute_squared_norm(std::shared_ptr<gpmesh::SplineAttributeCPU<T>> in);


//Y = alpha * X + beta * Y;
template<typename T>
inline void axpy(Eigen::Matrix<T, -1, 1> &Y,
                 Eigen::Matrix<T, -1, 1> &X,
                 T alpha,
                 T beta) {
    Y = alpha * X + beta * Y;
}

template<typename T, typename U>
inline void axpy(T &Y, T &X, U alpha, U beta) {
    Y = alpha * X + beta * Y;
}


template<typename T>
inline void matvec_serial(const gpmesh::AssemblyContextCPU<T> &assembly_context,
                          uint32_t num_geometry,
                          std::shared_ptr<gpmesh::SplineAttributeCPU<T>> in,
                          std::shared_ptr<gpmesh::SplineAttributeCPU<T>> out) {
    for (size_t geometry_id = 0; geometry_id < num_geometry; ++geometry_id) {
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