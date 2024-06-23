#pragma once
#include <stdint.h>
#include <string>
#include "Util/GP_Macros.h"

namespace gpmesh {

/**
 * @brief Local vertex type (wrapped around uint16_t)
 */
struct LocalVertexT
{
    /**
     * @brief Default constructor
     */
    __device__ __host__ LocalVertexT() : id(INVALID16)
    {
    }

    /**
     * @brief Constructor using local index
     * @param id vertex local index in the owner patch
     * @return
     */
    __device__ __host__ LocalVertexT(uint16_t id) : id(id)
    {
    }
    uint16_t id;
};

/**
 * @brief Local edge type (wrapped around uint16_t)
 */
struct LocalEdgeT
{
    /**
     * @brief Default constructor
     */
    __device__ __host__ LocalEdgeT() : id(INVALID16)
    {
    }

    /**
     * @brief Constructor using local index
     * @param id edge local index in the owner patch
     * @return
     */
    __device__ __host__ LocalEdgeT(uint16_t id) : id(id)
    {
    }
    uint16_t id;
};

/**
 * @brief Local face type (wrapped around uint16_t)
 */
struct LocalFaceT
{
    /**
     * @brief Default constructor
     */
    __device__ __host__ LocalFaceT() : id(INVALID16)
    {
    }

    /**
     * @brief Constructor using local index
     * @param id face local index in the owner patch
     * @return
     */
    __device__ __host__ LocalFaceT(uint16_t id) : id(id)
    {
    }
    uint16_t id;
};

struct LocalControlPointT
{
    __device__ __host__ LocalControlPointT() : id(INVALID16) {};

    __device__ __host__ LocalControlPointT(uint16_t _id) : id(_id) {};

    uint16_t id;
};

}  // namespace gpmesh