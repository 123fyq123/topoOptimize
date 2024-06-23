#pragma once

#include <cuda_runtime_api.h>
#include <stdint.h>
#include <string>
#include <utility>
#include "GP_Local.h"
#include "Util/GP_Macros.h"
#include <Eigen/Core>

namespace gpmesh {

/**
 * @brief GeometryInfo stores the information needed for query operations in a
 * patch
 */
    struct ALIGN(16) GeometryInfo {

        // Non-owned mesh elements patch ID
        uint32_t *not_owned_Geometry_cp;

        // Non-owned mesh elements local ID
        LocalControlPointT *not_owned_id_cp;

        // Number of mesh elements in the patch
        uint16_t num_CPoints;

        // Number of mesh elements owned by this patch
        uint16_t num_owned_CPoints;

        //Spline Geometry Information
        double *u_knots, *v_knots, *w_knots;
        int u_size, v_size, w_size;
        int u_degree, v_degree, w_degree;

        double *coefs;

        double *cpoints_displacement;

        // The index of this patch (relative to all other patches)
        uint32_t geometry_id;
    };

    struct ALIGN(16) AssemblyInfo {
        double *local_mats;
        double *local_rhs;
        uint32_t geometry_id;
    };

    struct ALIGN(16) MassInfo {
        double *mass_local_mats;
        uint32_t geometry_id;
    };

    struct ALIGN(16) SplineAttributeGlobalMappingInfo
    {
        uint32_t* ltog_mapping;
        uint32_t geometry_id;
    };

    enum class BoundarySide {
        NONE = 0,
        WEST = 1,
        EAST = 2,
        SOUTH = 3,
        NORTH = 4,
        FRONT = 5,
        BACK = 6
    };

    struct ALIGN(16) DirichletInfo {
        BoundarySide side;
        uint16_t dim;
        double value;

        //belonging geometry id
        uint32_t geometry_id;

        uint32_t dirichlet_boundary_id;
    };

    //try to use only one patch-local info if this dirichlet point is on boundary
    struct ALIGN(16) DirichletPointInfo
    {
        double valueX, valueY, valueZ;

        uint32_t global_id;

        //if this point is on boundary and has many patch-local info, choose first one
        uint32_t geometry_id;
        uint16_t local_id;
    };

    template<typename T>
    struct ALIGN(16) DirichletPointInfoCPU
    {
        T valueX = 0, valueY = 0, valueZ = 0;

        uint32_t global_id = INVALID32;

        uint32_t geometry_id = INVALID32;
        uint16_t local_id = INVALID16;
    };

    template<typename T>
    struct ALIGN(16) GeometryInfo_t {

        // Non-owned mesh elements patch ID
        uint32_t *not_owned_Geometry_cp;

        // Non-owned mesh elements local ID
        LocalControlPointT *not_owned_id_cp;

        // Number of mesh elements in the patch
        uint16_t num_CPoints;

        // Number of mesh elements owned by this patch
        uint16_t num_owned_CPoints;

        //Spline Geometry Information
        T *u_knots, *v_knots, *w_knots;
        int u_size, v_size, w_size;
        int u_degree, v_degree, w_degree;

        T *coefs;

        T *cpoints_displacement;

        // The index of this patch (relative to all other patches)
        uint32_t geometry_id;
    };

    template<typename T>
    struct ALIGN(16) AssemblyInfo_t {
        T *local_mats;
        T *local_rhs;
        uint32_t geometry_id;
    };

    template<typename T>
    struct ALIGN(16) AssemblyInfoCPU {
        T *local_mats = nullptr;
        T *local_rhs = nullptr;
        uint32_t geometry_id = INVALID64;

        void release()
        {
            if(local_mats != nullptr) delete[] local_mats;
            if(local_rhs != nullptr) delete[] local_rhs;
        }
    };

    template<typename T>
    struct ALIGN(16) PhysicsInfoCPU
    {
        Eigen::Matrix<T, 64, 3> m_cpoints_displacement;
        uint32_t geometry_id;
    };

    template<typename T>
    struct ALIGN(16) MassInfo_t {
        T *mass_local_mats;
        uint32_t geometry_id;
    };

    template<typename T>
    struct ALIGN(16) MassInfoCPU {
        T *mass_local_mats = nullptr;
        uint32_t geometry_id = INVALID64;

        void release()
        {
            if(mass_local_mats != nullptr) delete[] mass_local_mats;
        }
    };

    template<typename T>
    struct ALIGN(16) DirichletInfo_t {
        BoundarySide side;
        uint16_t dim;
        T value;

        //belonging geometry id
        uint32_t geometry_id;

        uint32_t dirichlet_boundary_id;
    };
}  // namespace gpmesh