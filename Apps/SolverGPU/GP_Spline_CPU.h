#pragma once

#include <memory>
#include <unordered_map>
#include <vector>
#include "GP_Context.h"
#include "GP_Geometry_info.h"
#include "GP_Types.h"
#include "Util/GP_Log.h"
#include "Util/GP_Macros.h"
#include "Util/GP_Util.h"
#include "GP_Spline_Attribute_CPU.h"
#include <Hyperelasticity/CBSplineSolid.h>
#include <Hyperelasticity/SparseSystem.h>

#include "GP_BoundaryCPoint.h"


namespace gpmesh {

/**
 * @brief The main class for creating GPSplineSolid data structure. It takes an input
 * mesh on the host, computes the patches, and creates the data structure on the
 * GPU. It is not mean to be used directly by the user. Users should use
 * GPSplineStatic instead
 */

    template<typename T>
    class GPSplineSolidCPU {
    private:

        bool m_quite = true;
        bool bInitDirichlet = false;
        bool bInitPatch = false;

        uint32_t m_num_geometry = INVALID32;
        uint32_t m_num_cpoints = INVALID32;
        uint16_t m_cpoint_num_per_geometry = 64;

        //装配矩阵
        AssemblyContextCPU<T> m_assembly_context;

        //质量矩阵
        MassContextCPU<T> m_mass_context;

        //控制点位移
        PhysicsContextCPU<T> m_physics_context;

        //Dirichlet边界条件
        DirichletPointContextCPU<T> m_dirichlet_point_context;

        SplineAttributeContainerCPU<T> m_attr_container;

    public:
        //建立映射关系
        SparseSystem m_sparse_system;

    public:
        virtual ~GPSplineSolidCPU();

        GPSplineSolidCPU() = default;

        GPSplineSolidCPU(const GPSplineSolidCPU &) = delete;

        void init(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                  std::vector<boundary_point> &neumann_points, const bool quite);

        void init_patch(const vector<CBSplineSolid> &patch);

        void init_dirichlet_points(std::vector<boundary_point> &dirichlet_points);

        void build_assembly_context(const vector<CBSplineSolid> &patch);

        void build_mass_context(const vector<CBSplineSolid> &patch);

        void build_physics_context(const vector<CBSplineSolid> &patch);

        void build_dirichlet_point_context(const std::vector<boundary_point> &dirichlet_points);

        void build_rhs(shared_ptr<SplineAttributeCPU<T>> rhs, std::vector<boundary_point> &dirichlet_points,
                       std::vector<boundary_point> &neumann_points);

        void debug_init(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                        std::vector<boundary_point> &neumann_points);

        void debug_assembly(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                            std::vector<boundary_point> &neumann_points);

        void debug_attribute_l2g(const vector<CBSplineSolid> &patch, vector<boundary_point> &dirichlet_points,
                                 vector<boundary_point> &neumann_points);

        void debug_attribute_g2l(const vector<CBSplineSolid> &patch, vector<boundary_point> &dirichlet_points,
                                 vector<boundary_point> &neumann_points);

        void debug_solving(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                           std::vector<boundary_point> &neumann_points,
                           Eigen::Matrix<T, -1, -1> &U);

        void linear_elasticity_pcg_solving(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                                  std::vector<boundary_point> &neumann_points,
                                                  Eigen::Matrix<T, -1, -1> &U,
                                                  T Yangs, T possion);

        void linear_elasticity_pardiso(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                       std::vector<boundary_point> &neumann_points,
                                       Eigen::Matrix<T, -1, -1> &U,
                                       T Yangs, T possion);

        void linear_elasticity_Eigen(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                       std::vector<boundary_point> &neumann_points,
                                       Eigen::Matrix<T, -1, -1> &U,
                                       T Yangs, T possion, bool bUseCG);

        inline void helper_set_neumann_point_in_sparse_system_rhs(std::vector<boundary_point>& neumann_points)
        {
            for (auto bd_point : neumann_points)
            {
                assert(bd_point.m_condition_type == 1);

                uint32_t global_id = bd_point.m_cpoint_global_id;
                m_sparse_system.rhs()(global_id + m_num_cpoints * 0) = bd_point.m_val.getX();
                m_sparse_system.rhs()(global_id + m_num_cpoints * 1) = bd_point.m_val.getY();
                m_sparse_system.rhs()(global_id + m_num_cpoints * 2) = bd_point.m_val.getZ();
            }
        }

        inline void helper_set_dirichlet_point_in_sparse_system_rhs(std::vector<boundary_point>& dirichlet_points)
        {
            const double       MAXNUM_DIRICHLET = 1e10;
            for (auto bd_point : dirichlet_points)
            {
                assert(bd_point.m_condition_type == 0);

                uint32_t global_id = bd_point.m_cpoint_global_id;
                m_sparse_system.rhs()(global_id + m_num_cpoints * 0) = bd_point.m_val.getX() * MAXNUM_DIRICHLET;
                m_sparse_system.rhs()(global_id + m_num_cpoints * 1) = bd_point.m_val.getY() * MAXNUM_DIRICHLET;
                m_sparse_system.rhs()(global_id + m_num_cpoints * 2) = bd_point.m_val.getZ() * MAXNUM_DIRICHLET;
            }
        }

        void helper_PCG_Solver_JacobiPrecondition(Eigen::Matrix<T, -1, -1> &result, float &time, T threshold = 1e-13,
                                                  bool bQuiet = true);

        void helper_pardiso_solver(Eigen::Matrix<T, -1, -1> &result, float &time, T threshold = 1e-13,
                                   bool bQuiet = true);

        void helper_eigen_solver(Eigen::Matrix<T, -1, -1> &result, float &time, bool bUseCG = true, T threshold = 1e-13,
                                 bool bQuiet = true);

        std::shared_ptr<SplineAttributeCPU<T>>
        helper_add_spline_attribute(const std::string &name, const uint32_t num_attributes, layoutT layout = SoA);

        void
        helper_Comp_PCG_Solver(Eigen::Matrix<T, -1, -1> &result, float &time, T threshold = 1e-13, bool bQuiet = true);

        void helper_build_jacobi_preconditioner_from_assembly(std::shared_ptr<SplineAttributeCPU<T>> Preconditioner);

        void helper_check_assembly_context();

        inline uint32_t get_num_cpoints()
        {
            return m_num_cpoints;
        }

        void fill_assembly_context(T fill_data) const;

        void fill_physics_context(T fill_data) const;

        //accessor of spline context
        const AssemblyContextCPU<T> &get_assembly_context() const {
            return m_assembly_context;
        }

        const MassContextCPU<T> &get_mass_context() const {
            return m_mass_context;
        }

        const PhysicsContextCPU<T> &get_physics_context() const {
            return m_physics_context;
        }

        const DirichletPointContextCPU<T> &get_dirichlet_piont_context() const {
            return m_dirichlet_point_context;
        };
    };


    //explicit specification of template, for def, imp separation
    template
    class GPSplineSolidCPU<double>;

    template
    class GPSplineSolidCPU<float>;

}  // namespace gpmesh


