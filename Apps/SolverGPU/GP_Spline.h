#pragma once
#include <memory>
#include <unordered_map>
#include <vector>
#include "GP_Context.h"
#include "GP_Geometry_info.h"
#include "Patcher/GP_Patcher.h"
#include "GP_Types.h"
#include "Util/GP_Log.h"
#include "Util/GP_Macros.h"
#include "Util/GP_Util.h"
#include <Hyperelasticity/CBSplineSolid.h>
#include <Hyperelasticity/SparseSystem.h>

#include "GP_BoundaryCPoint.h"

class GPSplineTest;
class GPSplineSolidV0;

namespace gpmesh {

/**
 * @brief The main class for creating GPSplineSolid data structure. It takes an input
 * mesh on the host, computes the patches, and creates the data structure on the
 * GPU. It is not mean to be used directly by the user. Users should use
 * GPSplineStatic instead
 */
class GPSplineSolid
{
   public:
    /**
     * @brief Total number of vertices in the mesh
     */
    uint32_t get_num_cpoints() const
    {
        return m_num_cpoints;
    }

    /**
     * @brief Return a context that store various information about the mesh on
     * the GPU
     */
    const SplineContext& get_spline_context() const
    {
        return m_gpspline_context;
    }

    const GaussianIntegralContext& get_gaussian_integral_context() const
    {
        return m_gs_integral_context;
    }

    const PhysicContext& get_physic_context() const
    {
        return m_physic_context;
    }

    const AssemblyContext& get_assembly_context() const
    {
        return m_assembly_context;
    }

    const DirichletContext& get_dirichlet_context() const
    {
        return m_dirichlet_context;
    }

    const DirichletPointContext& get_dirichelt_point_context() const
    {
        return m_dirichlet_point_context;
    }

    const MassContext& get_mass_context() const
    {
        return m_mass_context;
    }

    /**
     * @brief returns the patch size used during partitioning the input mesh
     */
    uint32_t get_geometry_size() const
    {
        return m_geometry_size;
    }

    /**
     * @brief Total number of patches of the input mesh
     */
    uint32_t get_num_geometry() const
    {
        return m_num_geometry;
    }

    /**
     * @brief Maximum number of control points in a patch
     */
    uint32_t get_per_patch_max_cpoints() const
    {
        return m_max_cp_per_geometry;
    }


    GPSplineSolid();

   public:
    virtual ~GPSplineSolid();

    GPSplineSolid(const GPSplineSolid&) = delete;

    void init(const vector<CBSplineSolid>& patch, std::vector<boundary_condition>& m_dirichletSides,
         std::vector<boundary_condition>& m_neumannSides, const bool quite);

    void init(const vector<CBSplineSolid>& patch, std::vector<boundary_point>& dirichlet_points,
              std::vector<boundary_point>& neumann_points, const bool quite);

    //two steps initialization
    void init_patch(const vector<CBSplineSolid>& patch);
    void init_dirichlet_points(const std::vector<boundary_point>& dirichlet_points);

    void setYoungPossion(double _Young, double _Poisson);
    /**
     * @brief Calculate various statistics for the input mesh
     *
     * Calculate max valence, max edge incident faces, max face adjacent faces,
     * if the input is closed, if the input is edge manifold, and max number of
     * vertices/edges/faces per patch
     *
     * @param fv input face incident vertices
     * @param ef input edge incident faces
     */
    //void calc_statistics();

    //void calc_max_not_owned_elements();

    void build(const std::vector<CBSplineSolid>& patch);

    void build_single_geometry(const std::vector<CBSplineSolid>& patch,
                               const uint32_t geometry_id);
//    void build_all_patch_parallel(const std::vector<std::vector<uint32_t>>& fv);
//    void build_single_patch(const std::vector<std::vector<uint32_t>>& fv,
//                            const uint32_t                            patch_id);
//
//    void build_single_patch_ltog(const std::vector<std::vector<uint32_t>>& fv,
//                                 const uint32_t patch_id);

    void build_device(const std::vector<CBSplineSolid>& patch);
    void build_gaussian_integral_context(const std::vector<CBSplineSolid>& patch);
    void build_physic_context(const std::vector<CBSplineSolid>& patch);
    void build_assembly_context(const std::vector<CBSplineSolid>& patch);
    void build_dirichlet_context(const std::vector<boundary_condition>& m_dirichletSides);
    void build_dirichlet_point_context(const vector<boundary_point>& dirichlet_points);
    void build_attr_mapping_context();

    inline double get_lambda() const { return lambda; }
    inline double get_mu() const { return mu; }

    friend class ::GPSplineTest;
    friend class ::GPSplineSolidV0;

    SplineContext m_gpspline_context;
    GaussianIntegralContext m_gs_integral_context;
    PhysicContext m_physic_context;
    AssemblyContext m_assembly_context;
    DirichletContext m_dirichlet_context;
    DirichletPointContext m_dirichlet_point_context;
    MassContext m_mass_context;
    SplineAttributeGlobalMappingContext m_attr_mapping_context;

    uint32_t m_num_cpoints;

    uint32_t m_max_cp_per_geometry;

    uint32_t m_max_not_owned_cp;

    uint32_t       m_num_geometry;
    const uint32_t m_geometry_size;

    bool           m_quite;

    int m_dims = 3;
    // the number of owned mesh elements per geometry
    std::vector<uint16_t> m_h_num_owned_cp;

    // mappings
    // local to global map for (cp) control points
    std::vector<std::vector<uint32_t>> m_h_geometry_ltog_cp;
    // global to local map for control points
    std::vector<std::vector<uint32_t>> m_h_geometry_gtol_cp;

    GeometryInfo *m_d_geometry_info, *m_h_geometry_info;
    AssemblyInfo *m_d_assembly_info, *m_h_assembly_info;
    MassInfo *m_d_mass_info, *m_h_mass_info;

    Eigen::Matrix<double, 64, 3>* m_h_cpoints_displacement_array;

    SparseSystem m_sparse_system;

    //cpoints_num * dims
    Eigen::VectorXd m_h_cpoints_displacement;
    //coefs_num * dims
    Eigen::MatrixXd m_h_coefs_displacement;
    //double* m_d_cpoints_displacement_per_geometry; //[geometry_id][local_cp_id * dims + x/y/z]

    //flag for initialization
    bool bInitPatch = false;
    bool bInitDirichlet = false;

    double Young = 1000;
    double Poisson = 0.3;


    double lambda = 0;
    double mu = 0;
};
}  // namespace gpmesh
