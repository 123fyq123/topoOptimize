#pragma once

#include "GP_Spline_Attribute.h"
#include "Kernels/GP_CUDA_Attribute.cuh"

namespace gpmesh
{
    template <typename T>
    class ReduceHandle
    {

    public:
        ReduceHandle()                    = delete;
        ReduceHandle(const ReduceHandle&) = default;

        /**
         * @brief Constructor which allocates internal memory used in all reduce
         * operations
         * @param attr one of Attribute used for subsequent reduction
         * operations
         */
        ReduceHandle(const SplineAttributeGlobal<T>& attr) : m_num_blocks(32),
              m_num_threads(attr.m_block_size)
        {
            uint32_t elements_num = attr.m_num_cpoints * 3;
            uint32_t tmp_block_num = DIVIDE_UP(elements_num, m_num_threads);

            m_num_blocks = std::min(tmp_block_num, m_num_blocks);

            CUDA_ERROR(cudaMalloc(&m_d_reduce_1st_stage, m_num_blocks * sizeof(T)));

            CUDA_ERROR(cudaMalloc(&m_d_reduce_2nd_stage, sizeof(T)));

            m_d_reduce_temp_storage = NULL;
            cub::DeviceReduce::Sum(m_d_reduce_temp_storage,
                                   m_reduce_temp_storage_bytes,
                                   m_d_reduce_1st_stage,
                                   m_d_reduce_2nd_stage,
                                   m_num_blocks);

            CUDA_ERROR(
                cudaMalloc(&m_d_reduce_temp_storage, m_reduce_temp_storage_bytes));
        }

        ~ReduceHandle()
        {
            GPU_FREE(m_d_reduce_1st_stage);
            GPU_FREE(m_d_reduce_2nd_stage);
            GPU_FREE(m_d_reduce_temp_storage);
            m_reduce_temp_storage_bytes = 0;
        }


        /**
         * @brief compute L2 norm between two input attributes and return the output
         * on the host
         * @param attr input attribute
         * @param attribute_id specific attribute ID to compute its norm2. Default
         * is INVALID32 which compute norm2 for all attributes
         * @param stream stream to run the computation on
         * @return the output of L2 norm on the host
         */
        T squaredNorm(const SplineAttributeGlobal<T>& attr,
                cudaStream_t        stream       = NULL)
        {

            if((attr.get_global_allocated() & DEVICE) != DEVICE)
            {
                GPMESH_ERROR(
                    "ReduceHandle::norm2() input attribute to should be "
                    "allocated on the device"
                    );

                return -1;
            }

            detail::squaredNorm_kernel_spline_attribute_global<T, attr.m_block_size><<<m_num_blocks, attr.m_block_size, 0, stream>>>
                (attr, attr.m_num_cpoints * 3, m_d_reduce_1st_stage);

            //print_arr_host<T><<<1,1>>>(m_num_blocks, m_d_reduce_1st_stage);

            return reduce_2nd_stage(stream, cub::Sum(), 0.0);
        }

        T dot(const SplineAttributeGlobal<T>& attr1,
              const SplineAttributeGlobal<T>& attr2,
              cudaStream_t stream = NULL)
        {
            if((attr1.get_global_allocated() & DEVICE) != DEVICE)
            {
                GPMESH_ERROR(
                    "ReduceHandle::dot() input attr1 to should be "
                    "allocated on the device"
                );

                assert(0);
            }

            if((attr2.get_global_allocated() & DEVICE) != DEVICE)
            {
                GPMESH_ERROR(
                    "ReduceHandle::dot() input attr2 to should be "
                    "allocated on the device"
                );
                assert(0);
            }

            assert(attr1.m_block_size == attr2.m_block_size);
            assert(attr1.m_num_cpoints == attr2.m_num_cpoints);

//            detail::dot_kernel_spline_attribute_global<T, attr1.m_block_size><<<m_num_blocks, attr1.m_block_size, 0, stream>>>
//                (attr1, attr2, attr1.m_num_cpoints * 3, m_d_reduce_1st_stage);


            detail::dot_kernel_spline_attribute_global<T, attr1.m_block_size><<<m_num_blocks, attr1.m_block_size, 0, stream>>>
                            (attr1, attr2, attr1.m_num_cpoints * 3, m_d_reduce_1st_stage);

//            detail::dot_kernel_spline_attribute_global<T, attr1.m_block_size><<<m_num_blocks, attr1.m_block_size, 0, stream>>>
//                (attr1.m_d_global_attr, attr2.m_d_global_attr, attr1.m_num_cpoints * 3, m_d_reduce_1st_stage);

            return reduce_2nd_stage(stream, cub::Sum(), 0.0);
        }

        uint32_t get_num_block() const
        {
            return m_num_blocks;
        }

        uint32_t  get_num_thread() const
        {
            return m_num_threads;
        }

    private:
        template <typename ReductionOp>
        T reduce_2nd_stage(cudaStream_t stream, ReductionOp reduction_op, T init)
        {
            T h_output = 0;

            cub::DeviceReduce::Reduce(m_d_reduce_temp_storage,
                                      m_reduce_temp_storage_bytes,
                                      m_d_reduce_1st_stage,
                                      m_d_reduce_2nd_stage,
                                      m_num_blocks,
                                      reduction_op,
                                      init,
                                      stream);

            CUDA_ERROR(cudaMemcpy(&h_output,
                                       m_d_reduce_2nd_stage,
                                       sizeof(T),
                                       cudaMemcpyDeviceToHost));

            return h_output;
        }

        size_t   m_reduce_temp_storage_bytes;
        T*       m_d_reduce_1st_stage;
        T*       m_d_reduce_2nd_stage;
        void*    m_d_reduce_temp_storage;
        uint32_t m_num_blocks = 32;
        const uint32_t m_num_threads = 256;
    };
}