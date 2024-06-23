#pragma once

#include <stdint.h>
#include "GP_Patch_Info.h"
#include "GP_Geometry_info.h"
#include "Util/GP_Macros.h"
#include <Eigen/Core>
#include <Hyperelasticity/SparseSystem.h>

namespace gpmesh {

/**
 * @brief context for the mesh parameters and pointers. Everything is allocated
 * on and managed by gpmesh. This class is meant to be a vehicle to copy various
 * parameters to the device kernels.
 */
class Context
{
   public:
    friend class GPMesh;
    friend class GPMeshDynamic;
    friend class GPMeshStatic;

    /**
     * @brief Default constructor
     */
    Context()
        : m_d_num_edges(nullptr),
          m_d_num_faces(nullptr),
          m_d_num_vertices(nullptr),
          m_d_num_patches(nullptr),
          m_patches_info(nullptr)
    {
    }

    Context(const Context&) = default;

    /**
     * @brief Total number of edges in mesh
     */
    __device__ __forceinline__ uint32_t get_num_edges() const
    {
        return *m_d_num_edges;
    }

    /**
     * @brief Total number of faces in mesh
     */
    __device__ __forceinline__ uint32_t get_num_faces() const
    {
        return *m_d_num_faces;
    }

    /**
     * @brief Total number of vertices in mesh
     */
    __device__ __forceinline__ uint32_t get_num_vertices() const
    {
        return *m_d_num_vertices;
    }

    /**
     * @brief Total number of patches in mesh
     */
    __device__ __forceinline__ uint32_t get_num_patches() const
    {
        return *m_d_num_patches;
    }

    /**
     * @brief A pointer to device PatchInfo used to store various information
     * about the patches
     */
    __device__ __forceinline__ PatchInfo* get_patches_info() const
    {
        return m_patches_info;
    }

    /**
     * @brief Unpack an edge to its edge ID and direction
     * @param edge_dir The input packed edge as stored in PatchInfo and
     * internally in gpmesh
     * @param edge The unpacked edge ID
     * @param dir The unpacked edge direction
     */
    static __device__ __host__ __forceinline__ void
    unpack_edge_dir(const uint16_t edge_dir, uint16_t& edge, flag_t& dir)
    {
        dir  = (edge_dir & 1) != 0;
        edge = edge_dir >> 1;
    }

   private:
    /**
     * @brief initialize various members
     * @param num_edges total number of edges in the mesh
     * @param num_faces total number of faces in the mesh
     * @param num_vertices total number of vertices in the mesh
     * @param num_patches number of patches
     * @param patches pointer to PatchInfo that contains different info about
     * the patches
     */
    void init(const uint32_t num_edges,
              const uint32_t num_faces,
              const uint32_t num_vertices,
              const uint32_t num_patches,
              PatchInfo*     patches)
    {
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_vertices, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_edges, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_faces, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_patches, sizeof(uint32_t)));

        CUDA_ERROR(cudaMemcpy(m_d_num_vertices,
                              &num_vertices,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(
            m_d_num_edges, &num_edges, sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(
            m_d_num_faces, &num_faces, sizeof(uint32_t), cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_num_patches,
                              &num_patches,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));
        m_patches_info = patches;
    }

    void release()
    {
//        CUDA_ERROR(cudaFree(m_d_num_edges));
//        CUDA_ERROR(cudaFree(m_d_num_faces));
//        CUDA_ERROR(cudaFree(m_d_num_vertices));
//        CUDA_ERROR(cudaFree(m_d_num_patches));

        GPU_FREE(m_d_num_edges);
        GPU_FREE(m_d_num_faces);
        GPU_FREE(m_d_num_vertices);
        GPU_FREE(m_d_num_patches);
    }


    uint32_t * m_d_num_edges, *m_d_num_faces, *m_d_num_vertices, *m_d_num_patches;
    PatchInfo* m_patches_info;
};

class SplineContext
{
   public:
    friend class GPSplineSolid;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineV0;
    friend class GPSplineSolidV0;

    /**
     * @brief Default constructor
     */
    SplineContext()
        : m_d_num_cpoints(nullptr),
          m_d_num_geometry(nullptr),
          m_d_geometry_info(nullptr)
    {
    }

    SplineContext(const SplineContext&) = default;

    /**
     * @brief Total number of edges in mesh
     */
    __device__ __forceinline__ uint32_t get_num_cpoins() const
    {
        return *m_d_num_cpoints;
    }

    __device__ __forceinline__ uint32_t get_num_geometry() const
    {
        return *m_d_num_geometry;
    }

    /**
     * @brief A pointer to device PatchInfo used to store various information
     * about the patches
     */
    __device__ __forceinline__ GeometryInfo* get_geometry_info() const
    {
        return m_d_geometry_info;
    }


   private:
    /**
     * @brief initialize various members
     * @param num_edges total number of edges in the mesh
     * @param num_faces total number of faces in the mesh
     * @param num_vertices total number of vertices in the mesh
     * @param num_patches number of patches
     * @param patches pointer to PatchInfo that contains different info about
     * the patches
     */
    void init(const uint32_t num_cpoints,
              const uint32_t num_geometry,
              GeometryInfo*     patches)
    {
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_cpoints, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_num_geometry, sizeof(uint32_t)));

        CUDA_ERROR(cudaMemcpy(m_d_num_cpoints,
                              &num_cpoints,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));

        CUDA_ERROR(cudaMemcpy(m_d_num_geometry,
                              &num_geometry,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));
        m_d_geometry_info = patches;
    }

    void release()
    {
        //        CUDA_ERROR(cudaFree(m_d_num_edges));
        //        CUDA_ERROR(cudaFree(m_d_num_faces));
        //        CUDA_ERROR(cudaFree(m_d_num_vertices));
        //        CUDA_ERROR(cudaFree(m_d_num_patches));

        GPU_FREE(m_d_num_cpoints);
        GPU_FREE(m_d_num_geometry);

    }

    uint32_t * m_d_num_cpoints, *m_d_num_geometry;
    //m_d_geometry_info[geometry_id]
    GeometryInfo* m_d_geometry_info;

};

class GaussianIntegralContext
{
    friend class GPSplineSolid;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineV0;
    friend class GPSplineSolidV0;

   private:
    uint32_t *m_d_gaussian_num_u, *m_d_gaussian_num_v, *m_d_gaussian_num_w;
    double *m_d_points_u, *m_d_points_v, *m_d_points_w;
    double *m_d_weights_u, *m_d_weights_v, *m_d_weights_w;
    double *m_d_basis_value_0, *m_d_basis_value_1;
    Eigen::Matrix<double, 64, 64>* m_d_basis_value_0_eigen_matrix;
    Eigen::Matrix<double, 192, 64>* m_d_basis_value_1_eigen_matrix;

   public:

    __device__ __host__ __forceinline__ uint32_t get_gaussian_num_u() const
    {
        return *m_d_gaussian_num_u;
    }

    __device__ __host__ __forceinline__ uint32_t get_gaussian_num_v() const
    {
        return *m_d_gaussian_num_v;
    }

    __device__ __host__ __forceinline__ uint32_t get_gaussian_num_w() const
    {
        return *m_d_gaussian_num_w;
    }

    __device__ __host__ __forceinline__ double* get_points_u() const
    {
        return m_d_points_u;
    }

    __device__ __host__ __forceinline__ double* get_points_v() const
    {
        return m_d_points_v;
    }

    __device__ __host__ __forceinline__ double* get_points_w()
    {
        return m_d_points_w;
    }

    __device__ __host__ __forceinline__ double* get_weights_u() const
    {
        return m_d_weights_u;
    }

    __device__ __host__ __forceinline__ double* get_weights_v() const
    {
        return m_d_weights_v;
    }

    __device__ __host__ __forceinline__ double* get_weights_w() const
    {
        return m_d_weights_w;
    }

    __device__ __host__ __forceinline__ double* get_basis_value_0() const
    {
        return m_d_basis_value_0;
    }

    __device__ __host__ __forceinline__ double* get_basis_value_1() const
    {
        return m_d_basis_value_1;
    }

    __device__ __host__ __forceinline__ Eigen::Matrix<double, 64, 64>* get_basis_value_0_matrix() const
    {
        return m_d_basis_value_0_eigen_matrix;
    }

    __device__ __host__ __forceinline__ Eigen::Matrix<double, 192, 64>* get_basis_value_1_matrix() const
    {
        return m_d_basis_value_1_eigen_matrix;
    }

    void init(uint32_t gaussian_num_u, uint32_t gaussian_num_v, uint32_t gaussian_num_w,
              double* points_u, double* points_v, double* points_w,
              double* weights_u,double* weights_v, double* weights_w)
    {
        CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_u, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_v, sizeof(uint32_t)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_w, sizeof(uint32_t)));

        CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_u,
                              &gaussian_num_u,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_v,
                              &gaussian_num_v,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_w,
                              &gaussian_num_w,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));

        CUDA_ERROR(cudaMalloc((void**)&m_d_points_u, gaussian_num_u * sizeof(double)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_points_v, gaussian_num_v * sizeof(double)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_points_w, gaussian_num_w * sizeof(double)));

        CUDA_ERROR(cudaMalloc((void**)&m_d_weights_u, gaussian_num_u * sizeof(double)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_weights_v, gaussian_num_v * sizeof(double)));
        CUDA_ERROR(cudaMalloc((void**)&m_d_weights_w, gaussian_num_w * sizeof(double)));

        CUDA_ERROR(cudaMemcpy(m_d_points_u,
                              points_u,
                              gaussian_num_u * sizeof(double),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_points_v,
                              points_v,
                              gaussian_num_v * sizeof(double),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_points_w,
                              points_w,
                              gaussian_num_w * sizeof(double),
                              cudaMemcpyHostToDevice));

        CUDA_ERROR(cudaMemcpy(m_d_weights_u,
                              weights_u,
                              gaussian_num_u * sizeof(double),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_weights_v,
                              weights_v,
                              gaussian_num_v * sizeof(double),
                              cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_weights_w,
                              weights_w,
                              gaussian_num_w * sizeof(double),
                              cudaMemcpyHostToDevice));

    }

    void init_derivatives(double* basis_value_0, size_t length_0, double* basis_value_1, size_t length_1)
    {
        CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_0, length_0 * sizeof(double)));
        CUDA_ERROR(cudaMemcpy(m_d_basis_value_0, basis_value_0,
                              length_0 * sizeof(double),
                              cudaMemcpyHostToDevice));

        CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_1, length_1 * sizeof(double)));
        CUDA_ERROR(cudaMemcpy(m_d_basis_value_1, basis_value_1,
                              length_1 * sizeof(double),
                              cudaMemcpyHostToDevice));
    }

    void init_derivatives_matrix(Eigen::Matrix<double, 64, 64> basis_value_0_matrix,
                                 Eigen::Matrix<double, 192, 64> basis_value_1_matrix)
    {
        CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_0_eigen_matrix, sizeof(Eigen::Matrix<double, 64, 64>)));
        CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_1_eigen_matrix, sizeof(Eigen::Matrix<double, 192, 64>)));
        CUDA_ERROR(cudaMemcpy(m_d_basis_value_0_eigen_matrix, &basis_value_0_matrix,
                              sizeof(Eigen::Matrix<double, 64, 64>), cudaMemcpyHostToDevice));
        CUDA_ERROR(cudaMemcpy(m_d_basis_value_1_eigen_matrix, &basis_value_1_matrix,
                              sizeof(Eigen::Matrix<double, 192, 64>), cudaMemcpyHostToDevice));
    }

    void release()
    {
        GPU_FREE(m_d_gaussian_num_u);
        GPU_FREE(m_d_gaussian_num_v);
        GPU_FREE(m_d_gaussian_num_w);
        GPU_FREE(m_d_points_u);
        GPU_FREE(m_d_points_v);
        GPU_FREE(m_d_points_w)
        GPU_FREE(m_d_weights_u);
        GPU_FREE(m_d_weights_v);
        GPU_FREE(m_d_weights_w);
        GPU_FREE(m_d_basis_value_0);
        GPU_FREE(m_d_basis_value_1);
        GPU_FREE(m_d_basis_value_0_eigen_matrix);
        GPU_FREE(m_d_basis_value_1_eigen_matrix);
    }


};

class PhysicContext
{
    friend class GPSplineSolid;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineV0;
    friend class GPSplineSolidV0;

   public:
    __device__ __host__ Eigen::Matrix<double, 64, 3>* get_cpoints_displacement_array() const
    {
        return m_d_cpoints_displacement;
    }

    __device__ __host__ Eigen::Matrix<double, 64, 3>* get_cpoints_displacement(uint32_t geometry_id) const
    {
        return &(m_d_cpoints_displacement[geometry_id]);
    }

    __device__ __host__ Eigen::Matrix<double, 6, 6>* get_C() const
    {
        return m_d_C;
    }

   private:
    Eigen::Matrix<double, 64, 3>* m_d_cpoints_displacement;

    Eigen::Matrix<double, 6, 6>* m_d_C;

    uint32_t* m_d_num_geometry;

    void init(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement, SparseSystem& system, uint32_t num_geometry)
    {
        assert(solution.rows() == 3 * system.getControlNums());

        system.constructCoefs(solution, displacement);

        CUDA_ERROR(cudaMalloc((void**)&m_d_num_geometry, sizeof(uint32_t)));
        CUDA_ERROR(cudaMemcpy(m_d_num_geometry,
                              &num_geometry,
                              sizeof(uint32_t),
                              cudaMemcpyHostToDevice));

        std::vector<Eigen::Matrix<double, 64, 3>>
            m_h_cpoints_displacement_per_geometry(num_geometry, Eigen::Matrix<double, 64, 3>::Zero());

        for (int i = 0; i < num_geometry; ++i) {
            // default situation : each geometry have 64 cpoints
            m_h_cpoints_displacement_per_geometry[i] = displacement.block(
                system.get_cpoints_inclusive_sum()[i] - 64, 0, 64, 3);
        }

        CUDA_ERROR(
            cudaMalloc((void**)&m_d_cpoints_displacement,
                       num_geometry * sizeof(Eigen::Matrix<double, 64, 3>)));
        CUDA_ERROR(
            cudaMemcpy(m_d_cpoints_displacement,
                       m_h_cpoints_displacement_per_geometry.data(),
                       num_geometry * sizeof(Eigen::Matrix<double, 64, 3>),
                       cudaMemcpyHostToDevice));
    }

    void init_C(const Eigen::Matrix<double, 6, 6> C)
    {
        CUDA_ERROR(cudaMalloc((void**)&m_d_C,
                              sizeof(Eigen::Matrix<double, 6, 6>)));
        CUDA_ERROR(cudaMemcpy(m_d_C, &C, sizeof(Eigen::Matrix<double, 6, 6>),
                              cudaMemcpyHostToDevice));
    }
    //Release Strategy
    void release()
    {
        GPU_FREE(m_d_cpoints_displacement);
        GPU_FREE(m_d_C);
        GPU_FREE(m_d_num_geometry);
    }

public:
    //update Strategy
    void update_host_info(Eigen::Matrix<double, 64, 3>* m_h_cpoints_displacement_array, uint32_t num_geometry)
    {
        Eigen::Matrix<double, 64, 3>* m_host_map_displacement_info = new Eigen::Matrix<double, 64, 3>[num_geometry];

        CUDA_ERROR(cudaMemcpy(m_host_map_displacement_info,
                              m_d_cpoints_displacement,
                              num_geometry * sizeof(Eigen::Matrix<double, 64, 3>),
                              cudaMemcpyDeviceToHost));

        for (int i = 0;  i < num_geometry; ++i)
        {
            m_h_cpoints_displacement_array[i] = m_host_map_displacement_info[i];
        }

        delete[] m_host_map_displacement_info;
    }



};

class AssemblyContext
{
    friend class GPSplineSolid;
    friend class GPSplineSolidTest;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineSolidV0;

   public:

    __host__ __device__ AssemblyInfo*  get_assembly_info() const
    {
        return m_d_assembly_info;
    }

    __host__ __device__ AssemblyInfo& get_assembly_info(uint32_t geometry_id) const
    {
        return m_d_assembly_info[geometry_id];
    }

    void get_host_local_matrix(double** &local_matrix, uint32_t geometry_size) const{
        auto m_host_map_assembly_info = (AssemblyInfo*) malloc(geometry_size * sizeof(AssemblyInfo));
        CUDA_ERROR(cudaMemcpy(m_host_map_assembly_info,
                              m_d_assembly_info,
                              geometry_size * sizeof(AssemblyInfo),
                              cudaMemcpyDeviceToHost));
        local_matrix = (double**) malloc(geometry_size * sizeof(double*));
        for(int i = 0; i < geometry_size; ++i){
            local_matrix[i] = (double*) malloc(192 * 192 * sizeof(double));
            CUDA_ERROR(cudaMemcpy(local_matrix[i], m_host_map_assembly_info[i].local_mats,
                                  192 * 192 * sizeof(double), cudaMemcpyDeviceToHost));
        }
        free(m_host_map_assembly_info);
    }

   private:

    AssemblyInfo* m_d_assembly_info;

    void init(AssemblyInfo* assembly_info)
    {
        m_d_assembly_info = assembly_info;
    }

    void release()
    {
        GPU_FREE(m_d_assembly_info);
    }

    void update_host_info(AssemblyInfo* m_h_assembly_info, uint32_t geometry_size)
    {
        AssemblyInfo* m_host_map_assembly_info = (AssemblyInfo*) malloc(geometry_size * sizeof(AssemblyInfo));

        CUDA_ERROR(cudaMemcpy(m_host_map_assembly_info,
                              m_d_assembly_info,
                              geometry_size * sizeof(AssemblyInfo),
                              cudaMemcpyDeviceToHost));
        for (int i = 0;  i < geometry_size; ++i)
        {
            m_h_assembly_info[i].geometry_id = m_host_map_assembly_info[i].geometry_id;
            CUDA_ERROR(cudaMemcpy(m_h_assembly_info[i].local_mats, m_host_map_assembly_info[i].local_mats,
                                  192 * 192 * sizeof(double), cudaMemcpyDeviceToHost));
            CUDA_ERROR(cudaMemcpy(m_h_assembly_info[i].local_rhs, m_host_map_assembly_info[i].local_rhs,
                                  192 * sizeof(double), cudaMemcpyDeviceToHost));
        }

        //CUDA_ERROR(cudaStreamSynchronize(0));
        free(m_host_map_assembly_info);
    }
};

class DirichletContext
{
    friend class GPSplineSolid;
    friend class GPSplineSolidTest;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineSolidV0;

   private:

    DirichletInfo* m_d_dirichlet_info = nullptr;
    uint32_t m_num_dirichlet_boundary = 0;

   public:

    __host__ __device__ DirichletInfo* get_dirichlet_info()
    {
        return m_d_dirichlet_info;
    }

    __host__ __device__ DirichletInfo get_dirichlet_info(uint32_t info_id) const
    {
        assert(info_id < m_num_dirichlet_boundary);

        return m_d_dirichlet_info[info_id];
    }

    __host__ __device__ uint32_t get_num_dirichlet_boundary() const
    {
        return m_num_dirichlet_boundary;
    }

   public:

    void init(DirichletInfo* dirichlet_info, uint32_t num_dirichlet_boundary_num)
    {
        m_num_dirichlet_boundary = num_dirichlet_boundary_num;

        CUDA_ERROR(cudaMalloc((void**)&m_d_dirichlet_info,
                              num_dirichlet_boundary_num * sizeof(DirichletInfo)));
        CUDA_ERROR(cudaMemcpy(m_d_dirichlet_info, dirichlet_info,
                              num_dirichlet_boundary_num * sizeof(DirichletInfo),
                              cudaMemcpyHostToDevice));
    }

    void release()
    {
        GPU_FREE(m_d_dirichlet_info);
    }

};

class DirichletPointContext

{
    friend class GPSplineSolid;
    friend class GPSplineSolidTest;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineSolidV0;

   private:

    DirichletPointInfo* m_d_dirichlet_point_info = nullptr;
    uint32_t m_num_dirichlet_point = 0;

   public:

    __host__ __device__ DirichletPointInfo* get_dirichlet_point_info()
    {
        return m_d_dirichlet_point_info;
    }

    __host__ __device__ DirichletPointInfo get_dirichlet_point_info(uint32_t info_id) const
    {
        assert(info_id < m_num_dirichlet_point);

        return m_d_dirichlet_point_info[info_id];
    }

    __host__ __device__ uint32_t get_num_dirichlet_point() const
    {
        return m_num_dirichlet_point;
    }

   public:

    void init(DirichletPointInfo* h_dirichlet_point_info, uint32_t num_dirichlet_point)
    {
        m_num_dirichlet_point = num_dirichlet_point;

        CUDA_ERROR(cudaMalloc((void**)&m_d_dirichlet_point_info,
                              m_num_dirichlet_point * sizeof(DirichletPointInfo)));
        CUDA_ERROR(cudaMemcpy(m_d_dirichlet_point_info, h_dirichlet_point_info,
                              m_num_dirichlet_point * sizeof(DirichletPointInfo),
                              cudaMemcpyHostToDevice));
    }

    void release()
    {
        GPU_FREE(m_d_dirichlet_point_info);
    }

};

class MassContext
{
    friend class GPSplineSolid;
    friend class GPSplineSolidTest;
    friend class GPSplineDynamic;
    friend class GPSplineStatic;
    friend class GPSplineTest;
    friend class GPSplineSolidV0;

   public:

    __host__ __device__ MassInfo*  get_mass_info() const
    {
        return m_d_mass_info;
    }

    __host__ __device__ MassInfo& get_mass_info(uint32_t geometry_id) const
    {
        return m_d_mass_info[geometry_id];
    }

   private:

    MassInfo* m_d_mass_info;

    void init(MassInfo* d_mass_info)
    {
        m_d_mass_info = d_mass_info;
    }

    void release()
    {
        GPU_FREE(m_d_mass_info);
    }

    void update_host_info(MassInfo* h_mass_info, uint32_t geometry_size)
    {
        MassInfo* m_host_map_mass_info = (MassInfo*) malloc(geometry_size * sizeof(MassInfo));

        CUDA_ERROR(cudaMemcpy(m_host_map_mass_info,
                              m_d_mass_info,
                              geometry_size * sizeof(MassInfo),
                              cudaMemcpyDeviceToHost));

        for (int i = 0;  i < geometry_size; ++i)
        {
            h_mass_info[i].geometry_id = m_host_map_mass_info[i].geometry_id;
            CUDA_ERROR(cudaMemcpy(h_mass_info[i].mass_local_mats, m_host_map_mass_info[i].mass_local_mats,
                                  192 * 192 * sizeof(double), cudaMemcpyDeviceToHost));
        }

        free(m_host_map_mass_info);
    }
};

    class SplineAttributeGlobalMappingContext
    {
        friend class GPSplineSolid;
        friend class GPSplineSolidTest;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

    public:

        __host__ __device__ SplineAttributeGlobalMappingInfo*  get_attribute_global_mapping_info() const
        {
            return m_d_attr_mapping_info;
        }

        __host__ __device__ SplineAttributeGlobalMappingInfo& get_attribute_global_mapping_info(uint32_t geometry_id) const
        {
            return m_d_attr_mapping_info[geometry_id];
        }

    private:

        SplineAttributeGlobalMappingInfo* m_d_attr_mapping_info;

        void init(SplineAttributeGlobalMappingInfo* d_mapping_info)
        {
            m_d_attr_mapping_info = d_mapping_info;
        }

        void release()
        {
            //may cause memory leak
            GPU_FREE(m_d_attr_mapping_info);
        }
    };

    template<typename T>
    class SplineContext_t
    {
    public:
        friend class GPSplineSolid;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

        /**
         * @brief Default constructor
         */
        SplineContext_t()
                : m_d_num_cpoints(nullptr),
                  m_d_num_geometry(nullptr),
                  m_d_geometry_info(nullptr)
        {
        }

        SplineContext_t(const SplineContext_t&) = default;

        /**
         * @brief Total number of edges in mesh
         */
        __device__ __forceinline__ uint32_t get_num_cpoins() const
        {
            return *m_d_num_cpoints;
        }

        __device__ __forceinline__ uint32_t get_num_geometry() const
        {
            return *m_d_num_geometry;
        }

        /**
         * @brief A pointer to device PatchInfo used to store various information
         * about the patches
         */
        __device__ __forceinline__ GeometryInfo_t<T>* get_geometry_info() const
        {
            return m_d_geometry_info;
        }


    private:
        /**
         * @brief initialize various members
         * @param num_edges total number of edges in the mesh
         * @param num_faces total number of faces in the mesh
         * @param num_vertices total number of vertices in the mesh
         * @param num_patches number of patches
         * @param patches pointer to PatchInfo that contains different info about
         * the patches
         */

        void init(const uint32_t num_cpoints,
                  const uint32_t num_geometry,
                  GeometryInfo_t<T>*     patches)
        {
            CUDA_ERROR(cudaMalloc((void**)&m_d_num_cpoints, sizeof(uint32_t)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_num_geometry, sizeof(uint32_t)));

            CUDA_ERROR(cudaMemcpy(m_d_num_cpoints,
                                  &num_cpoints,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_num_geometry,
                                  &num_geometry,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));
            m_d_geometry_info = patches;
        }

        void release()
        {
            GPU_FREE(m_d_num_cpoints);
            GPU_FREE(m_d_num_geometry);
        }

        uint32_t * m_d_num_cpoints, *m_d_num_geometry;
        GeometryInfo_t<T> m_d_geometry_info;
    };

    template<typename T>
    class GaussianIntegralContext_t
    {
        friend class GPSplineSolid;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

    private:
        uint32_t *m_d_gaussian_num_u, *m_d_gaussian_num_v, *m_d_gaussian_num_w;
        T *m_d_points_u, *m_d_points_v, *m_d_points_w;
        T *m_d_weights_u, *m_d_weights_v, *m_d_weights_w;
        T *m_d_basis_value_0, *m_d_basis_value_1;
        Eigen::Matrix<T, 64, 64>* m_d_basis_value_0_eigen_matrix;
        Eigen::Matrix<T, 192, 64>* m_d_basis_value_1_eigen_matrix;

    public:

        __device__ __host__ __forceinline__ uint32_t get_gaussian_num_u() const
        {
            return *m_d_gaussian_num_u;
        }

        __device__ __host__ __forceinline__ uint32_t get_gaussian_num_v() const
        {
            return *m_d_gaussian_num_v;
        }

        __device__ __host__ __forceinline__ uint32_t get_gaussian_num_w() const
        {
            return *m_d_gaussian_num_w;
        }

        __device__ __host__ __forceinline__ T* get_points_u() const
        {
            return m_d_points_u;
        }

        __device__ __host__ __forceinline__ T* get_points_v() const
        {
            return m_d_points_v;
        }

        __device__ __host__ __forceinline__ T* get_points_w()
        {
            return m_d_points_w;
        }

        __device__ __host__ __forceinline__ T* get_weights_u() const
        {
            return m_d_weights_u;
        }

        __device__ __host__ __forceinline__ T* get_weights_v() const
        {
            return m_d_weights_v;
        }

        __device__ __host__ __forceinline__ T* get_weights_w() const
        {
            return m_d_weights_w;
        }

        __device__ __host__ __forceinline__ T* get_basis_value_0() const
        {
            return m_d_basis_value_0;
        }

        __device__ __host__ __forceinline__ T* get_basis_value_1() const
        {
            return m_d_basis_value_1;
        }

        __device__ __host__ __forceinline__ Eigen::Matrix<T, 64, 64>* get_basis_value_0_matrix() const
        {
            return m_d_basis_value_0_eigen_matrix;
        }

        __device__ __host__ __forceinline__ Eigen::Matrix<T, 192, 64>* get_basis_value_1_matrix() const
        {
            return m_d_basis_value_1_eigen_matrix;
        }

        void init(uint32_t gaussian_num_u, uint32_t gaussian_num_v, uint32_t gaussian_num_w,
                  T* points_u, T* points_v, T* points_w,
                  T* weights_u,T* weights_v, T* weights_w)
        {
            CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_u, sizeof(uint32_t)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_v, sizeof(uint32_t)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_gaussian_num_w, sizeof(uint32_t)));

            CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_u,
                                  &gaussian_num_u,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_v,
                                  &gaussian_num_v,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_gaussian_num_w,
                                  &gaussian_num_w,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMalloc((void**)&m_d_points_u, gaussian_num_u * sizeof(T)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_points_v, gaussian_num_v * sizeof(T)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_points_w, gaussian_num_w * sizeof(T)));

            CUDA_ERROR(cudaMalloc((void**)&m_d_weights_u, gaussian_num_u * sizeof(T)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_weights_v, gaussian_num_v * sizeof(T)));
            CUDA_ERROR(cudaMalloc((void**)&m_d_weights_w, gaussian_num_w * sizeof(T)));

            CUDA_ERROR(cudaMemcpy(m_d_points_u,
                                  points_u,
                                  gaussian_num_u * sizeof(T),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_points_v,
                                  points_v,
                                  gaussian_num_v * sizeof(T),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_points_w,
                                  points_w,
                                  gaussian_num_w * sizeof(T),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_weights_u,
                                  weights_u,
                                  gaussian_num_u * sizeof(T),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_weights_v,
                                  weights_v,
                                  gaussian_num_v * sizeof(T),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_weights_w,
                                  weights_w,
                                  gaussian_num_w * sizeof(T),
                                  cudaMemcpyHostToDevice));

        }

        void init_derivatives(T* basis_value_0, size_t length_0, T* basis_value_1, size_t length_1)
        {
            CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_0, length_0 * sizeof(double)));
            CUDA_ERROR(cudaMemcpy(m_d_basis_value_0, basis_value_0,
                                  length_0 * sizeof(T),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_1, length_1 * sizeof(T)));
            CUDA_ERROR(cudaMemcpy(m_d_basis_value_1, basis_value_1,
                                  length_1 * sizeof(T),
                                  cudaMemcpyHostToDevice));
        }

        void init_derivatives_matrix(Eigen::Matrix<T, 64, 64> basis_value_0_matrix,
                                     Eigen::Matrix<T, 192, 64> basis_value_1_matrix)
        {
            CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_0_eigen_matrix, sizeof(Eigen::Matrix<T, 64, 64>)));
            CUDA_ERROR(cudaMalloc((void**) &m_d_basis_value_1_eigen_matrix, sizeof(Eigen::Matrix<T, 192, 64>)));
            CUDA_ERROR(cudaMemcpy(m_d_basis_value_0_eigen_matrix, &basis_value_0_matrix,
                                  sizeof(Eigen::Matrix<T, 64, 64>), cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMemcpy(m_d_basis_value_1_eigen_matrix, &basis_value_1_matrix,
                                  sizeof(Eigen::Matrix<T, 192, 64>), cudaMemcpyHostToDevice));
        }

        void release()
        {
            GPU_FREE(m_d_gaussian_num_u);
            GPU_FREE(m_d_gaussian_num_v);
            GPU_FREE(m_d_gaussian_num_w);
            GPU_FREE(m_d_points_u);
            GPU_FREE(m_d_points_v);
            GPU_FREE(m_d_points_w)
            GPU_FREE(m_d_weights_u);
            GPU_FREE(m_d_weights_v);
            GPU_FREE(m_d_weights_w);
            GPU_FREE(m_d_basis_value_0);
            GPU_FREE(m_d_basis_value_1);
            GPU_FREE(m_d_basis_value_0_eigen_matrix);
            GPU_FREE(m_d_basis_value_1_eigen_matrix);
        }


    };

    template<typename T>
    class PhysicContext_t
    {
        friend class GPSplineSolid;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

    public:
        __device__ __host__ Eigen::Matrix<T, 64, 3>* get_cpoints_displacement() const
        {
            return m_d_cpoints_displacement;
        }

        __device__ __host__ Eigen::Matrix<T, 64, 3>* get_cpoints_displacement(uint32_t geometry_id) const
        {
            return &(m_d_cpoints_displacement[geometry_id]);
        }

        __device__ __host__ Eigen::Matrix<T, 6, 6>* get_C() const
        {
            return m_d_C;
        }

    private:
        Eigen::Matrix<T, 64, 3>* m_d_cpoints_displacement;

        Eigen::Matrix<T, 6, 6>* m_d_C;

        uint32_t* m_d_num_geometry;

        void init(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement, SparseSystem& system, uint32_t num_geometry)
        {
            assert(solution.rows() == 3 * system.getControlNums());

            system.constructCoefs(solution, displacement);

            CUDA_ERROR(cudaMalloc((void**)&m_d_num_geometry, sizeof(uint32_t)));
            CUDA_ERROR(cudaMemcpy(m_d_num_geometry,
                                  &num_geometry,
                                  sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));

            std::vector<Eigen::Matrix<T, 64, 3>>
                    m_h_cpoints_displacement_per_geometry(
                    num_geometry, Eigen::Matrix<T, 64, 3>::Zero());

            for (int i = 0; i < num_geometry; ++i) {
                // default situation : each geometry have 64 cpoints
                m_h_cpoints_displacement_per_geometry[i] = displacement.block(
                        system.get_cpoints_inclusive_sum()[i] - 64, 0, 64, 3);
            }

            CUDA_ERROR(
                    cudaMalloc((void**)&m_d_cpoints_displacement,
                               num_geometry * sizeof(Eigen::Matrix<T, 64, 3>)));
            CUDA_ERROR(
                    cudaMemcpy(m_d_cpoints_displacement,
                               m_h_cpoints_displacement_per_geometry.data(),
                               num_geometry * sizeof(Eigen::Matrix<T, 64, 3>),
                               cudaMemcpyHostToDevice));
        }

        void init_C(const Eigen::Matrix<T, 6, 6> C)
        {
            CUDA_ERROR(cudaMalloc((void**)&m_d_C,
                                  sizeof(Eigen::Matrix<T, 6, 6>)));
            CUDA_ERROR(cudaMemcpy(m_d_C, &C, sizeof(Eigen::Matrix<T, 6, 6>),
                                  cudaMemcpyHostToDevice));
        }
        //update Strategy

        //Release Strategy
        void release()
        {
            GPU_FREE(m_d_cpoints_displacement);
            GPU_FREE(m_d_C);
            GPU_FREE(m_d_num_geometry);
        }
    };

    template<typename T>
    class AssemblyContext_t
    {
        friend class GPSplineSolid;
        friend class GPSplineSolidTest;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

        using m_AssemblyInfo = AssemblyInfo_t<T>;

    public:

        __host__ __device__ m_AssemblyInfo*  get_assembly_info() const
        {
            return m_d_assembly_info;
        }

        __host__ __device__ m_AssemblyInfo& get_assembly_info(uint32_t geometry_id) const
        {
            return m_d_assembly_info[geometry_id];
        }

    private:

        m_AssemblyInfo* m_d_assembly_info;

        void init(m_AssemblyInfo* assembly_info)
        {
            m_d_assembly_info = assembly_info;
        }

        void release()
        {
            GPU_FREE(m_d_assembly_info);
        }

        void update_host_info(m_AssemblyInfo* m_h_assembly_info, uint32_t geometry_size)
        {
            m_AssemblyInfo* m_host_map_assembly_info = (m_AssemblyInfo*) malloc(geometry_size * sizeof(m_AssemblyInfo));

            CUDA_ERROR(cudaMemcpy(m_host_map_assembly_info,
                                  m_d_assembly_info,
                                  geometry_size * sizeof(m_AssemblyInfo),
                                  cudaMemcpyDeviceToHost));
            for (int i = 0;  i < geometry_size; ++i)
            {
                m_h_assembly_info[i].geometry_id = m_host_map_assembly_info[i].geometry_id;
                CUDA_ERROR(cudaMemcpy(m_h_assembly_info[i].local_mats, m_host_map_assembly_info[i].local_mats,
                                      192 * 192 * sizeof(T), cudaMemcpyDeviceToHost));
                CUDA_ERROR(cudaMemcpy(m_h_assembly_info[i].local_rhs, m_host_map_assembly_info[i].local_rhs,
                                      192 * sizeof(T), cudaMemcpyDeviceToHost));
            }

            free(m_host_map_assembly_info);
        }
    };

    template<typename T>
    class DirichletContext_t
    {
        friend class GPSplineSolid;
        friend class GPSplineSolidTest;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

        using m_DirchletInfo = DirichletInfo_t<T>;

    private:

        m_DirchletInfo* m_d_dirichlet_info;
        uint32_t m_num_dirichlet_boundary;

    public:

        __host__ __device__ m_DirchletInfo* get_dirichlet_info()
        {
            return m_d_dirichlet_info;
        }

        __host__ __device__ m_DirchletInfo get_dirichlet_info(uint32_t info_id) const
        {
            assert(info_id < m_num_dirichlet_boundary);

            return m_d_dirichlet_info[info_id];
        }

        __host__ __device__ uint32_t get_num_dirichlet_boundary() const
        {
            return m_num_dirichlet_boundary;
        }

    public:

        void init(m_DirchletInfo* dirichlet_info, uint32_t num_dirichlet_boundary_num)
        {
            m_num_dirichlet_boundary = num_dirichlet_boundary_num;

            CUDA_ERROR(cudaMalloc((void**)&m_d_dirichlet_info,
                                  num_dirichlet_boundary_num * sizeof(m_DirchletInfo)));
            CUDA_ERROR(cudaMemcpy(m_d_dirichlet_info, dirichlet_info,
                                  num_dirichlet_boundary_num * sizeof(m_DirchletInfo),
                                  cudaMemcpyHostToDevice));
        }

        void release()
        {
            GPU_FREE(m_d_dirichlet_info);
        }

    };

    template<typename T>
    class MassContext_t
    {
        friend class GPSplineSolid;
        friend class GPSplineSolidTest;
        friend class GPSplineDynamic;
        friend class GPSplineStatic;
        friend class GPSplineTest;
        friend class GPSplineSolidV0;

        using m_MassInfo = MassInfo_t<T>;

    public:

        __host__ __device__ m_MassInfo*  get_mass_info() const
        {
            return m_d_mass_info;
        }

        __host__ __device__ m_MassInfo& get_mass_info(uint32_t geometry_id) const
        {
            return m_d_mass_info[geometry_id];
        }

    private:

        m_MassInfo* m_d_mass_info;

        void init(m_MassInfo* d_mass_info)
        {
            m_d_mass_info = d_mass_info;
        }

        void release()
        {
            GPU_FREE(m_d_mass_info);
        }

        void update_host_info(m_MassInfo* h_mass_info, uint32_t geometry_size)
        {
            m_MassInfo* m_host_map_mass_info = (m_MassInfo*) malloc(geometry_size * sizeof(m_MassInfo));

            CUDA_ERROR(cudaMemcpy(m_host_map_mass_info,
                                  m_d_mass_info,
                                  geometry_size * sizeof(m_MassInfo),
                                  cudaMemcpyDeviceToHost));

            for (int i = 0;  i < geometry_size; ++i)
            {
                h_mass_info[i].geometry_id = m_host_map_mass_info[i].geometry_id;
                CUDA_ERROR(cudaMemcpy(h_mass_info[i].mass_local_mats, m_host_map_mass_info[i].mass_local_mats,
                                      192 * 192 * sizeof(T), cudaMemcpyDeviceToHost));
            }

            free(m_host_map_mass_info);
        }
    };

    template<typename T>
    class MassContextCPU
    {
    private:
        using m_MassInfo = MassInfoCPU<T>;
        m_MassInfo* m_h_mass_info = nullptr;

    public:

        void init(m_MassInfo* h_mass_info)
        {
            m_h_mass_info = h_mass_info;
        }

        void release()
        {
            if(m_h_mass_info != nullptr) delete[] m_h_mass_info;
        }

        __host__ m_MassInfo*  get_mass_info() const
        {
            return m_h_mass_info;
        }

        __host__ m_MassInfo& get_mass_info(uint32_t geometry_id) const
        {
            return m_h_mass_info[geometry_id];
        }


    };

    template<typename T>
    class AssemblyContextCPU
    {
        using m_AssemblyInfo = AssemblyInfoCPU<T>;

    private:

        m_AssemblyInfo* m_h_assembly_info = nullptr;

    public:

        void init(m_AssemblyInfo* h_assembly_info)
        {
            m_h_assembly_info = h_assembly_info;
        }

        void release()
        {
            if(m_h_assembly_info!= nullptr) delete[] m_h_assembly_info;
        }

        m_AssemblyInfo*  get_assembly_info() const
        {
            return m_h_assembly_info;
        }

        m_AssemblyInfo& get_assembly_info(uint32_t geometry_id) const
        {
            return m_h_assembly_info[geometry_id];
        }
    };

    template<typename T>
    class PhysicsContextCPU
    {
    private:
        using m_PhysicsInfo = PhysicsInfoCPU<T>;

        m_PhysicsInfo* m_h_physics_info = nullptr;

    public:
        void init(m_PhysicsInfo* h_physics_info)
        {
            m_h_physics_info = h_physics_info;
        }

        void release()
        {
            if(m_h_physics_info!= nullptr) delete[] m_h_physics_info;
        }

        m_PhysicsInfo*  get_physics_info() const
        {
            return m_h_physics_info;
        }

        m_PhysicsInfo& get_physics_info(uint32_t geometry_id) const
        {
            return m_h_physics_info[geometry_id];
        }
    };

    template<typename T>
    class DirichletPointContextCPU

    {

    private:
        using m_DirichletInfo = DirichletPointInfoCPU<T>;
        m_DirichletInfo* m_h_dirichlet_point_info = nullptr;
        uint32_t m_num_dirichlet_point = 0;

    public:

        m_DirichletInfo* get_dirichlet_point_info() const
        {
            return m_h_dirichlet_point_info;
        }

        m_DirichletInfo& get_dirichlet_point_info(uint32_t info_id) const
        {
            assert(info_id < m_num_dirichlet_point);

            return m_h_dirichlet_point_info[info_id];
        }

        uint32_t get_num_dirichlet_point() const
        {
            return m_num_dirichlet_point;
        }

    public:

        void init(m_DirichletInfo* h_dirichlet_point_info, uint32_t num_dirichlet_point)
        {
            m_num_dirichlet_point = num_dirichlet_point;

            m_h_dirichlet_point_info = h_dirichlet_point_info;
        }

        void release()
        {
            if(m_h_dirichlet_point_info != nullptr)
            {
                delete[] m_h_dirichlet_point_info;
            }
        }

    };


}  // namespace gpmesh