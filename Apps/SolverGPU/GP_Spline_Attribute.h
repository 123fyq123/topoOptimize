#pragma once

#include <assert.h>
#include <utility>

#include "GP_Attribute.h"
#include "GP_Handle.h"
#include "Kernels/GP_CUDA_Attribute.cuh"
#include "Kernels/GP_Collective.cuh"
#include "Kernels/GP_CUDA_Util.cuh"
#include "GP_Patch_Info.h"
#include "GP_Mesh.h"
#include "GP_Types.h"
#include "Util/GP_Cuda_Query.h"
#include "Util/GP_Log.h"
#include "Util/GP_Util.h"
#include "Util/GP_Vector.h"

class GPMeshTest;


namespace gpmesh {

template<typename T>
class SplineAttributeBase : public AttributeBase
{
   private:
    //default is shallow copy
    char*     m_name;
    uint32_t  m_num_attributes;
    locationT m_allocated;
    T**       m_h_attr;
    T**       m_h_ptr_on_device;
    T**       m_d_attr;

    uint32_t  m_num_geometry;

    //here 64 per geometry is default situation
    constexpr static uint32_t m_block_size = 256;

    uint16_t  m_cpoint_num_per_geometry;
    //in iGame spline, it is not used
    //uint16_t* m_d_element_per_patch;
    //uint16_t* m_h_element_per_patch;

    layoutT   m_layout;



   public:

    SplineAttributeBase(const char*            name,
                        uint32_t               num_geometry,
                        uint16_t               cpoint_num_per_geometry,
                        uint32_t               num_attributes,
                        locationT              location,
                        layoutT                layout)
    : m_name(nullptr), m_allocated(LOCATION_NONE), m_h_attr(nullptr), m_h_ptr_on_device(nullptr), m_d_attr(nullptr)
    {
        if (name != nullptr) {
            //strlen`count exclude \0
            this->m_name = (char*)malloc(sizeof(char) * (strlen(name) + 1));
            strcpy(this->m_name, name);
        }

        this->init(num_attributes, num_geometry, cpoint_num_per_geometry, location, layout);
    }

    const char* get_name() const override
    {
        return m_name;
    }

    void init(const uint32_t               num_attributes,
              const uint32_t               num_geometry,
              const uint16_t               cpoint_num_per_geometry = 64,
              locationT                    location = LOCATION_ALL,
              const layoutT                layout   = AoS)
    {
        release();
        m_num_geometry   = num_geometry;
        m_num_attributes = num_attributes;
        m_cpoint_num_per_geometry = cpoint_num_per_geometry;
        m_layout         = layout;

        if (m_num_geometry == 0) {
            return;
        }

        allocate(location);
    }

    void release(locationT location = LOCATION_ALL) override
    {
        if (((location & HOST) == HOST) && ((m_allocated & HOST) == HOST)) {
            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                free(m_h_attr[p]);
            }
            free(m_h_attr);
            m_h_attr = nullptr;
            //free(m_h_element_per_patch);
            //m_h_element_per_patch = nullptr;
            m_allocated           = m_allocated & (~HOST);
        }

        if (((location & DEVICE) == DEVICE) &&
            ((m_allocated & DEVICE) == DEVICE)) {
            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                GPU_FREE(m_h_ptr_on_device[p]);
            }
            GPU_FREE(m_d_attr);
            //GPU_FREE(m_d_element_per_patch);
            m_allocated = m_allocated & (~DEVICE);
        }
    }

    void allocate(locationT location)
    {

        if (m_num_geometry != 0) {

            if ((location & HOST) == HOST) {
                release(HOST);

                //Store T* array per patch
                m_h_attr = static_cast<T**>(malloc(sizeof(T*) * m_num_geometry));

                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    m_h_attr[p] = static_cast<T*>(malloc(
                        sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes));
                }

                m_allocated = m_allocated | HOST;
            }

            if ((location & DEVICE) == DEVICE) {
                release(DEVICE);

                CUDA_ERROR(cudaMalloc((void**)&(m_d_attr),
                                      sizeof(T*) * m_num_geometry));
                m_h_ptr_on_device =
                    static_cast<T**>(malloc(sizeof(T*) * m_num_geometry));


                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMalloc((void**)&(m_h_ptr_on_device[p]),
                                          sizeof(T) * m_cpoint_num_per_geometry *
                                              m_num_attributes));
                }
                CUDA_ERROR(cudaMemcpy(m_d_attr,
                                      m_h_ptr_on_device,
                                      sizeof(T*) * m_num_geometry,
                                      cudaMemcpyHostToDevice));

                m_allocated = m_allocated | DEVICE;
            }
        }
    }

    void reset(const T value, locationT location, cudaStream_t stream = NULL)
    {
        if ((location & DEVICE) == DEVICE) {

            assert((m_allocated & DEVICE) == DEVICE);

            const int threads = 64;
            detail::template memset_spline_attribute<T>
                <<<m_num_geometry, threads, 0, stream>>>(*this,
                                                        value,
                                                        m_cpoint_num_per_geometry,
                                                        m_num_geometry,
                                                        m_num_attributes);
        }


        if ((location & HOST) == HOST) {
            assert((m_allocated & HOST) == HOST);
#pragma omp parallel for
            for (int p = 0; p < static_cast<int>(m_num_geometry); ++p) {
                for (int e = 0; e < m_cpoint_num_per_geometry; ++e) {
                    m_h_attr[p][e] = value;
                }
            }
        }
    }

    __host__ __device__ __forceinline__ bool is_empty() const
    {
        return m_num_geometry == 0;
    }

    __host__ __device__ __forceinline__ T& RawOperator(const uint32_t geometry_id ,
                                                       const uint32_t id) const
    {
#ifdef __CUDA_ARCH__
        return m_d_attr[geometry_id][id];
#else
        return m_h_attr[geometry_id][id];
#endif
    }

    __host__ __device__ __forceinline__ T& operator()(const uint32_t geometry_id,
                                                      const uint16_t local_id,
                                                      const uint32_t attr) const
    {
        assert(geometry_id < m_num_geometry);
        assert(attr < m_num_attributes);

        const uint32_t pitch_x = (m_layout == AoS) ? m_num_attributes : 1;
#ifdef __CUDA_ARCH__
        const uint32_t pitch_y =
            (m_layout == AoS) ? 1 : m_cpoint_num_per_geometry;
        return m_d_attr[geometry_id][local_id * pitch_x + attr * pitch_y];
#else
        const uint32_t pitch_y =
            (m_layout == AoS) ? 1 : m_cpoint_num_per_geometry;
        return m_h_attr[geometry_id][local_id * pitch_x + attr * pitch_y];
#endif
    }



    void copy_from(SplineAttributeBase<T>& source,
                   locationT     source_flag,
                   locationT     dst_flag,
                   cudaStream_t  stream = NULL)
    {


        if (source.m_layout != m_layout) {
            GPMESH_ERROR(
                "Attribute::copy_from() does not support copy from "
                "source of different layout!");
        }

        if ((source_flag & LOCATION_ALL) == LOCATION_ALL &&
            (dst_flag & LOCATION_ALL) != LOCATION_ALL) {
            GPMESH_ERROR("Attribute::copy_from() Invalid configuration!");
        }

        if (m_num_attributes != source.get_num_attributes()) {
            GPMESH_ERROR(
                "Attribute::copy_from() number of attributes is "
                "different!");
        }

        if (this->is_empty() || this->m_num_patches == 0) {
            return;
        }

        // 1) copy from HOST to HOST
        if ((source_flag & HOST) == HOST && (dst_flag & HOST) == HOST) {
            if ((source_flag & source.m_allocated) != source_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because it was not allocated on host");
            }
            if ((dst_flag & m_allocated) != dst_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because location (this) was not allocated on host");
            }

            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                //���
                std::memcpy(
                    m_h_attr[p],
                    source.m_h_attr[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes);
            }
        }


        // 2) copy from DEVICE to DEVICE
        if ((source_flag & DEVICE) == DEVICE && (dst_flag & DEVICE) == DEVICE) {
            if ((source_flag & source.m_allocated) != source_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because it was not allocated on device");
            }
            if ((dst_flag & m_allocated) != dst_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because location (this) was not allocated on device");
            }

            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                CUDA_ERROR(cudaMemcpyAsync(
                    m_h_ptr_on_device[p],
                    source.m_h_ptr_on_device[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                    cudaMemcpyDeviceToDevice,
                    stream));
            }
        }


        // 3) copy from DEVICE to HOST
        if ((source_flag & DEVICE) == DEVICE && (dst_flag & HOST) == HOST) {
            if ((source_flag & source.m_allocated) != source_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because it was not allocated on host");
            }
            if ((dst_flag & m_allocated) != dst_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because location (this) was not allocated on device");
            }


            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                CUDA_ERROR(cudaMemcpyAsync(
                    m_h_attr[p],
                    source.m_h_ptr_on_device[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                    cudaMemcpyDeviceToHost,
                    stream));
            }
        }


        // 4) copy from HOST to DEVICE
        if ((source_flag & HOST) == HOST && (dst_flag & DEVICE) == DEVICE) {
            if ((source_flag & source.m_allocated) != source_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because it was not allocated on device");
            }
            if ((dst_flag & m_allocated) != dst_flag) {
                GPMESH_ERROR(
                    "Attribute::copy() copying source is not valid"
                    " because location (this) was not allocated on host");
            }


            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                CUDA_ERROR(cudaMemcpyAsync(
                    m_h_ptr_on_device[p],
                    source.m_h_attr[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                    cudaMemcpyHostToDevice,
                    stream));
            }
        }
    }

    void move(locationT source, locationT target, cudaStream_t stream = NULL)
    {
        if (source == target) {
            GPMESH_WARN(
                "Attribute::move() source ({}) and target ({}) "
                "are the same.",
                location_to_string(source),
                location_to_string(target));
            return;
        }

        if ((source == HOST || source == DEVICE) &&
            ((source & m_allocated) != source)) {
            GPMESH_ERROR(
                "Attribute::move() moving source is not valid"
                " because it was not allocated on source i.e., {}",
                location_to_string(source));
        }

        if (((target & HOST) == HOST || (target & DEVICE) == DEVICE) &&
            ((target & m_allocated) != target)) {
            GPMESH_WARN(
                "Attribute::move() allocating target before moving to {}",
                location_to_string(target));
            allocate(target);
        }

        if (this->m_num_geometry == 0) {
            return;
        }

        if (source == HOST && target == DEVICE) {
            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                CUDA_ERROR(cudaMemcpyAsync(
                    m_h_ptr_on_device[p],
                    m_h_attr[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                    cudaMemcpyHostToDevice,
                    stream));
            }
        } else if (source == DEVICE && target == HOST) {
            for (uint32_t p = 0; p < m_num_geometry; ++p) {
                CUDA_ERROR(cudaMemcpyAsync(
                    m_h_attr[p],
                    m_h_ptr_on_device[p],
                    sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                    cudaMemcpyDeviceToHost,
                    stream));
            }
        }
    }


};

    template<typename T>
    class SplineAttributeGlobal : public AttributeBase
    {
    private:
        //default is shallow copy
        char*     m_name;
        uint32_t  m_num_attributes = -1;
        locationT m_allocated;

        //local_attr (per geoemtry)
        T**       m_h_attr;
        T**       m_h_ptr_on_device;
        T**       m_d_attr;

        uint32_t  m_num_geometry = 100;

        //here 64 per geometry is default situation
        //constexpr static uint32_t m_block_size = 64;

        uint16_t  m_cpoint_num_per_geometry;


        //in iGame spline, it is not used
        //uint16_t* m_d_element_per_patch;
        //uint16_t* m_h_element_per_patch;

        layoutT   m_layout = 100;

    public:
        //global_attr (require mapping information)
        locationT m_global_allocated;
        T* m_h_global_attr;
        T* m_d_global_attr;
        uint32_t m_num_cpoints;

        constexpr static uint32_t m_block_size = 256;

    public:

        SplineAttributeGlobal(const char* name)
                : m_name(nullptr), m_allocated(LOCATION_NONE), m_global_allocated(LOCATION_NONE),
                m_h_attr(nullptr), m_h_ptr_on_device(nullptr), m_d_attr(nullptr),
                m_h_global_attr(nullptr), m_d_global_attr(nullptr)
        {
            if (name != nullptr) {
                //strlen`count exclude \0
                this->m_name = (char*)malloc(sizeof(char) * (strlen(name) + 1));
                strcpy(this->m_name, name);
            }
        }

        const char* get_name() const override
        {
            return m_name;
        }

        __host__ __device__ __forceinline__ locationT get_local_allocated() const
        {
            return this->m_allocated;
        }

        __host__ __device__ __forceinline__ locationT get_global_allocated() const
        {
            return this->m_global_allocated;
        }

        void allocate(locationT location)
        {
            allocate_global_part(location);
            allocate_local_part(location);
        }

        void release(locationT location) override
        {
            release_global_part(location);
            release_local_part(location);
        }

        //not finish
        void move_local_to_global(locationT location)
        {
            if(location == gpmesh::LOCATION_NONE || location == gpmesh::LOCATION_ALL)
            {
                //
                assert(0);
                return;
            }
        }

        void load_global_part(Eigen::MatrixX<T>& global_attr)
        {
            assert(global_attr.rows() == m_num_cpoints * m_num_attributes);

            if ((m_global_allocated & gpmesh::DEVICE) == gpmesh::DEVICE)
            {
                CUDA_ERROR(cudaMemcpy(m_d_global_attr,
                                      global_attr.data(),
                                      sizeof(T) * m_num_cpoints * m_num_attributes,
                                      cudaMemcpyHostToDevice));
            }
            else
            {
                assert(0);
            }
        }

        void move_local_to_global_gpu(SplineAttributeGlobalMappingContext attr_mapping_context)
        {

            if((m_allocated & DEVICE) == DEVICE &&
               (m_global_allocated & DEVICE) == DEVICE)
            {
                detail::spline_attr_local_to_global<T>
                        <<<m_num_geometry, m_cpoint_num_per_geometry>>>
                (this->m_d_global_attr, this->m_d_attr, attr_mapping_context, m_num_cpoints);

//                detail::spline_attr_local_to_global<T>
//                <<<m_num_geometry, m_cpoint_num_per_geometry>>>
//                (*this, attr_mapping_context, m_num_cpoints);
            }
            else
            {
                assert(0);
            }

            return;
        }

        void move_global_to_local_gpu(SplineAttributeGlobalMappingContext attr_mapping_context)
        {
            if((m_allocated & DEVICE) == DEVICE &&
                (m_global_allocated & DEVICE) == DEVICE)
            {
                detail::spline_attr_global_to_local<T>
                    <<<m_num_geometry, m_cpoint_num_per_geometry>>>
                    (this->m_d_global_attr, this->m_d_attr, attr_mapping_context, m_num_cpoints);

                //                detail::spline_attr_local_to_global<T>
                //                <<<m_num_geometry, m_cpoint_num_per_geometry>>>
                //                (*this, attr_mapping_context, m_num_cpoints);
            }
            else
            {
                assert(0);
            }

            return;
        }

        void init_global_part(const uint32_t num_attributes,
                              const uint32_t num_cpoints,
                              locationT location = LOCATION_ALL)
        {
            release_global_part();
            m_num_attributes = num_attributes;

            m_num_cpoints = num_cpoints;

            allocate_global_part(location);
        }

        void release_global_part(locationT location = LOCATION_ALL)
        {
            if (((location & HOST) == HOST) && ((m_global_allocated & HOST) == HOST)) {
                free(m_h_global_attr);
                m_h_global_attr = nullptr;
                //free(m_h_element_per_patch);
                //m_h_element_per_patch = nullptr;
                m_global_allocated  = m_global_allocated & (~HOST);
            }

            if (((location & DEVICE) == DEVICE) &&
                ((m_global_allocated & DEVICE) == DEVICE)) {
                GPU_FREE(m_d_global_attr);
                m_d_global_attr = nullptr;
                m_global_allocated = m_global_allocated & (~DEVICE);
            }
        }

        void allocate_global_part(locationT location)
        {
            if ((location & HOST) == HOST) {
                release_global_part(HOST);

                //Store T* array per patch
                m_h_global_attr = static_cast<T*>(malloc(sizeof(T) * m_num_cpoints * m_num_attributes));

                m_global_allocated = m_global_allocated | HOST;
            }

            if ((location & DEVICE) == DEVICE) {
                release_global_part(DEVICE);

                CUDA_ERROR(cudaMalloc((void**)&(m_d_global_attr),
                                      sizeof(T) * m_num_cpoints * m_num_attributes));

                m_global_allocated = m_global_allocated | DEVICE;
            }
        }

        void reset_global_part(const T value, locationT location, cudaStream_t stream = NULL)
        {
            if ((location & DEVICE) == DEVICE) {

                assert((m_global_allocated & DEVICE) == DEVICE);

                uint32_t num_elements = m_num_cpoints * m_num_attributes;
                const size_t threads_num = 256;
                size_t blocks_num = DIVIDE_UP(num_elements, threads_num);
                detail::memset_spline_attribute_global_part<double><<<blocks_num, threads_num>>>(*this, value, num_elements);
            }

            if ((location & HOST) == HOST) {
                if(m_h_global_attr == nullptr)
                {
                    throw std::runtime_error("m_h_global_attr is nullptr");
                }
                assert((m_global_allocated & HOST) == HOST);
#pragma omp parallel for
                for (int p = 0; p < static_cast<int>(m_num_cpoints * m_num_attributes); ++p) {
                    m_h_global_attr[p] = value;
                }
            }
        }

        void move_global_part(locationT source, locationT target, cudaStream_t stream = NULL)
        {
            if (source == target) {
                GPMESH_WARN(
                        "Attribute::move_global_part() source ({}) and target ({}) "
                        "are the same.",
                        location_to_string(source),
                        location_to_string(target));
                return;
            }

            if ((source == HOST || source == DEVICE) &&
                ((source & m_global_allocated) != source)) {
                GPMESH_ERROR(
                        "Attribute::move_global_part() moving source is not valid"
                        " because it was not allocated on source i.e., {}",
                        location_to_string(source));
            }

            if (((target & HOST) == HOST || (target & DEVICE) == DEVICE) &&
                ((target & m_global_allocated) != target)) {
                GPMESH_WARN(
                        "Attribute::move_global_part() allocating target before moving to {}",
                        location_to_string(target));
                allocate_global_part(target);
            }

            if (this->m_num_geometry == 0) {
                return;
            }

            if (source == HOST && target == DEVICE) {
                CUDA_ERROR(cudaMemcpy(
                            m_d_global_attr,
                            m_h_global_attr,
                            sizeof(T) * m_num_cpoints * m_num_attributes,
                            cudaMemcpyHostToDevice));

            } else if (source == DEVICE && target == HOST) {
                CUDA_ERROR(cudaMemcpy(
                            m_h_global_attr,
                            m_d_global_attr,
                            sizeof(T) * m_num_cpoints * m_num_attributes,
                            cudaMemcpyDeviceToHost));
            }
        }

        void init_local_part(const uint32_t    num_attributes,
                  const uint32_t               num_geometry,
                  const uint16_t               cpoint_num_per_geometry = 64,
                  locationT                    location = LOCATION_ALL,
                  const layoutT                layout   = AoS)
        {
            release_local_part();
            m_num_geometry   = num_geometry;
            m_num_attributes = num_attributes;
            m_cpoint_num_per_geometry = cpoint_num_per_geometry;
            m_layout         = layout;

            if (m_num_geometry == 0) {
                return;
            }

            allocate_local_part(location);
        }

        void release_local_part(locationT location = LOCATION_ALL)
        {
            if (((location & HOST) == HOST) && ((m_allocated & HOST) == HOST)) {
                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    free(m_h_attr[p]);
                }
                free(m_h_attr);
                m_h_attr = nullptr;
                //free(m_h_element_per_patch);
                //m_h_element_per_patch = nullptr;
                m_allocated           = m_allocated & (~HOST);
            }

            if (((location & DEVICE) == DEVICE) &&
                ((m_allocated & DEVICE) == DEVICE)) {
                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    GPU_FREE(m_h_ptr_on_device[p]);
                }
                GPU_FREE(m_d_attr);
                //GPU_FREE(m_d_element_per_patch);
                m_allocated = m_allocated & (~DEVICE);
            }
        }

        void allocate_local_part(locationT location)
        {

            if (m_num_geometry != 0) {

                if ((location & HOST) == HOST) {
                    release_local_part(HOST);

                    //Store T* array per patch
                    m_h_attr = static_cast<T**>(malloc(sizeof(T*) * m_num_geometry));

                    for (uint32_t p = 0; p < m_num_geometry; ++p) {
                        m_h_attr[p] = static_cast<T*>(malloc(
                                sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes));
                    }

                    m_allocated = m_allocated | HOST;
                }

                if ((location & DEVICE) == DEVICE) {
                    release_local_part(DEVICE);

                    CUDA_ERROR(cudaMalloc((void**)&(m_d_attr),
                                          sizeof(T*) * m_num_geometry));
                    m_h_ptr_on_device =
                            static_cast<T**>(malloc(sizeof(T*) * m_num_geometry));


                    for (uint32_t p = 0; p < m_num_geometry; ++p) {
                        CUDA_ERROR(cudaMalloc((void**)&(m_h_ptr_on_device[p]),
                                              sizeof(T) * m_cpoint_num_per_geometry *
                                              m_num_attributes));
                    }
                    CUDA_ERROR(cudaMemcpy(m_d_attr,
                                          m_h_ptr_on_device,
                                          sizeof(T*) * m_num_geometry,
                                          cudaMemcpyHostToDevice));

                    m_allocated = m_allocated | DEVICE;
                }
            }
        }

        void reset_local_part(const T value, locationT location, cudaStream_t stream = NULL)
        {
            if ((location & DEVICE) == DEVICE) {

                assert((m_allocated & DEVICE) == DEVICE);

                const int threads = 64;
                detail::template memset_spline_attribute<T>
                <<<m_num_geometry, threads, 0, stream>>>(*this,
                        value,
                        m_cpoint_num_per_geometry,
                        m_num_geometry,
                        m_num_attributes);
            }


            if ((location & HOST) == HOST) {
                assert((m_allocated & HOST) == HOST);
#pragma omp parallel for
                for (int p = 0; p < static_cast<int>(m_num_geometry); ++p) {
                    for (int e = 0; e < static_cast<int>(m_cpoint_num_per_geometry); ++e) {
                        m_h_attr[p][e] = value;
                    }
                }
            }
        }

        __host__ __device__ __forceinline__ bool is_empty() const
        {
            return m_num_geometry == 0;
        }

        //local part access
        __host__ __device__ __forceinline__ T& RawOperator(const uint32_t geometry_id ,
                                                           const uint32_t id) const
        {
#ifdef __CUDA_ARCH__
            return m_d_attr[geometry_id][id];
#else
            return m_h_attr[geometry_id][id];
#endif
        }

        __host__ __device__ __forceinline__ T& operator()(const uint32_t geometry_id,
                                                          const uint16_t local_id,
                                                          const uint32_t attr) const
        {
            assert(geometry_id < m_num_geometry);
            assert(attr < m_num_attributes);

            const uint32_t pitch_x = (m_layout == AoS) ? m_num_attributes : 1;
#ifdef __CUDA_ARCH__
            const uint32_t pitch_y =
                    (m_layout == AoS) ? 1 : m_cpoint_num_per_geometry;
            return m_d_attr[geometry_id][local_id * pitch_x + attr * pitch_y];
#else
            const uint32_t pitch_y =
            (m_layout == AoS) ? 1 : m_cpoint_num_per_geometry;
        return m_h_attr[geometry_id][local_id * pitch_x + attr * pitch_y];
#endif
        }

        __host__ __device__ __forceinline__ T* get_global_attr()
        {
#ifdef __CUDA_ARCH__
            return m_d_global_attr;
#else
            return m_h_global_attr;
#endif
        }

        void copy_from_local_part(SplineAttributeBase<T>& source,
                       locationT     source_flag,
                       locationT     dst_flag,
                       cudaStream_t  stream = NULL)
        {


            if (source.m_layout != m_layout) {
                GPMESH_ERROR(
                        "Attribute::copy_from() does not support copy from "
                        "source of different layout!");
            }

            if ((source_flag & LOCATION_ALL) == LOCATION_ALL &&
                (dst_flag & LOCATION_ALL) != LOCATION_ALL) {
                GPMESH_ERROR("Attribute::copy_from() Invalid configuration!");
            }

            if (m_num_attributes != source.get_num_attributes()) {
                GPMESH_ERROR(
                        "Attribute::copy_from() number of attributes is "
                        "different!");
            }

            if (this->is_empty() || this->m_num_patches == 0) {
                return;
            }

            // 1) copy from HOST to HOST
            if ((source_flag & HOST) == HOST && (dst_flag & HOST) == HOST) {
                if ((source_flag & source.m_allocated) != source_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because it was not allocated on host");
                }
                if ((dst_flag & m_allocated) != dst_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because location (this) was not allocated on host");
                }

                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    //���
                    std::memcpy(
                            m_h_attr[p],
                            source.m_h_attr[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes);
                }
            }


            // 2) copy from DEVICE to DEVICE
            if ((source_flag & DEVICE) == DEVICE && (dst_flag & DEVICE) == DEVICE) {
                if ((source_flag & source.m_allocated) != source_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because it was not allocated on device");
                }
                if ((dst_flag & m_allocated) != dst_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because location (this) was not allocated on device");
                }

                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMemcpyAsync(
                            m_h_ptr_on_device[p],
                            source.m_h_ptr_on_device[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                            cudaMemcpyDeviceToDevice,
                            stream));
                }

                CUDA_ERROR(cudaStreamSynchronize(0));
            }


            // 3) copy from DEVICE to HOST
            if ((source_flag & DEVICE) == DEVICE && (dst_flag & HOST) == HOST) {
                if ((source_flag & source.m_allocated) != source_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because it was not allocated on host");
                }
                if ((dst_flag & m_allocated) != dst_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because location (this) was not allocated on device");
                }


                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMemcpyAsync(
                            m_h_attr[p],
                            source.m_h_ptr_on_device[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                            cudaMemcpyDeviceToHost,
                            stream));
                }
                CUDA_ERROR(cudaStreamSynchronize(0));
            }


            // 4) copy from HOST to DEVICE
            if ((source_flag & HOST) == HOST && (dst_flag & DEVICE) == DEVICE) {
                if ((source_flag & source.m_allocated) != source_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because it was not allocated on device");
                }
                if ((dst_flag & m_allocated) != dst_flag) {
                    GPMESH_ERROR(
                            "Attribute::copy() copying source is not valid"
                            " because location (this) was not allocated on host");
                }


                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMemcpyAsync(
                            m_h_ptr_on_device[p],
                            source.m_h_attr[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                            cudaMemcpyHostToDevice,
                            stream));
                }
                CUDA_ERROR(cudaStreamSynchronize(0));
            }
        }

        void move_local_part(locationT source, locationT target, cudaStream_t stream = NULL)
        {
            if (source == target) {
                GPMESH_WARN(
                        "Attribute::move_local_part() source ({}) and target ({}) "
                        "are the same.",
                        location_to_string(source),
                        location_to_string(target));
                return;
            }

            if ((source == HOST || source == DEVICE) &&
                ((source & m_allocated) != source)) {
                GPMESH_ERROR(
                        "Attribute::move_local_part() moving source is not valid"
                        " because it was not allocated on source i.e., {}",
                        location_to_string(source));
            }

            if (((target & HOST) == HOST || (target & DEVICE) == DEVICE) &&
                ((target & m_allocated) != target)) {
                GPMESH_WARN(
                        "Attribute::move_local_part() allocating target before moving to {}",
                        location_to_string(target));
                allocate_local_part(target);
            }

            if (this->m_num_geometry == 0) {
                return;
            }

            if (source == HOST && target == DEVICE) {
                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMemcpyAsync(
                            m_h_ptr_on_device[p],
                            m_h_attr[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                            cudaMemcpyHostToDevice,
                            stream));
                }
            } else if (source == DEVICE && target == HOST) {
                for (uint32_t p = 0; p < m_num_geometry; ++p) {
                    CUDA_ERROR(cudaMemcpyAsync(
                            m_h_attr[p],
                            m_h_ptr_on_device[p],
                            sizeof(T) * m_cpoint_num_per_geometry * m_num_attributes,
                            cudaMemcpyDeviceToHost,
                            stream));
                }
            }
            CUDA_ERROR(cudaStreamSynchronize(stream));
        }
    };
/**
 * @brief Attribute container used to manage a collection of attributes by
 * RXMeshStatic
 */
class SplineAttributeContainer
{
   public:
    /**
     * @brief Default constructor
     */
    SplineAttributeContainer() = default;

    /**
     * @brief Destructor which releases all attribute managed by this container
     */
    virtual ~SplineAttributeContainer()
    {
        while (!m_attr_container.empty()) {
            m_attr_container.back()->release();
            m_attr_container.pop_back();
        }

        for (auto iter = m_attr_hash_container.begin(); iter != m_attr_hash_container.end(); ++iter)
        {
            (*iter).second->release();
        }
    }

    /**
     * @brief Number of attribute managed by this container
     */
    size_t size()
    {
        //return m_attr_container.size();
            return m_attr_hash_container.size();
    }

    /**
     * @brief get a list of name of the attributes managed by this container
     * @return
     */
    std::vector<std::string> get_attribute_names() const
    {
        std::vector<std::string> names;
//        for (size_t i = 0; i < m_attr_container.size(); ++i) {
//            names.push_back(m_attr_container[i]->get_name());
//        }

        for (auto iter = m_attr_hash_container.begin(); iter != m_attr_hash_container.end(); ++iter)
        {
            names.emplace_back((*iter).first);
        }
        return names;
    }

    /**
     * @brief add a new attribute to be managed by this container
     * @tparam AttrT attribute type
     * @param name unique name given to the attribute
     * @param element_per_patch number of mesh element owned by each patch
     * @param num_attributes number of attributes per mesh element
     * @param location where the attributes will be allocated
     * @param layout memory layout in case of num_attributes > 1
     * @return a shared pointer to the attribute
     */

    //add SplineAttributeBase
    template <typename AttrT>
    std::shared_ptr<AttrT> add(const char*            name,
                               uint32_t               num_geometry,
                               uint16_t               cpoint_num_per_geometry,
                               uint32_t               num_attributes,
                               locationT              location,
                               layoutT                layout)
    {
        if (does_exist(name)) {
            GPMESH_WARN(
                "AttributeContainer::add() adding an attribute with "
                "name {} already exists!",
                std::string(name));
        }

        auto new_attr = std::make_shared<AttrT>(
            name, num_geometry, cpoint_num_per_geometry, num_attributes, location, layout);

//        m_attr_container.push_back(
//            std::dynamic_pointer_cast<AttributeBase>(new_attr));

        std::string key(name);

        //check in compiling state
        m_attr_hash_container[key] = std::dynamic_pointer_cast<AttributeBase>(new_attr);

        return new_attr;
    }


    template<typename T>
    std::shared_ptr<SplineAttributeGlobal<T>> add_global_attr(const char*            name,
                                           uint32_t               num_geometry,
                                           uint16_t               cpoint_num_per_geometry,
                                           uint32_t               num_attributes,
                                           uint32_t               num_cpoints,
                                           locationT              location,
                                           layoutT                layout)
    {
        if (does_exist(name)) {
//            GPMESH_WARN(
//                    "AttributeContainer::add() adding an attribute with "
//                    "name {} already exists!",
//                    std::string(name));
//
//            GPMESH_WARN(
//                    "this is a danger operation, should be cautious!"
//                    );

            return std::dynamic_pointer_cast<SplineAttributeGlobal<T>>(m_attr_hash_container[string(name)]);
        }

        auto new_attr = std::make_shared<SplineAttributeGlobal<T>>(name);
        new_attr->init_local_part(num_attributes, num_geometry, cpoint_num_per_geometry, location, layout);
        new_attr->init_global_part(num_attributes, num_cpoints, location);


        std::string key = name;

        m_attr_hash_container[key] = std::dynamic_pointer_cast<AttributeBase>(new_attr);

        return new_attr;
    }

    /**
     * @brief Check if an attribute exists
     * @param name of the attribute
     */
    bool does_exist(const char* name)
    {

        return m_attr_hash_container.find(string(name)) != m_attr_hash_container.end();

        for (size_t i = 0; i < m_attr_container.size(); ++i) {
            if (!strcmp(m_attr_container[i]->get_name(), name)) {
                return true;
            }
        }
        return false;
    }

    /**
     * @brief remove an attribute and release its memory
     * @param name of the attribute
     */
    void remove(const char* name)
    {
        if (does_exist(name))
        {
            string key(name);
            m_attr_hash_container[key]->release();
            m_attr_hash_container.erase(key);
        }

        return;

        for (auto it = m_attr_container.begin(); it != m_attr_container.end();
             ++it) {

            if (!strcmp((*it)->get_name(), name)) {
                (*it)->release(LOCATION_ALL);
                m_attr_container.erase(it);
                break;
            }
        }
    }

   private:
    std::vector<std::shared_ptr<AttributeBase>> m_attr_container;
    //ʹ��unordered_map���Ӻ���
    std::unordered_map<std::string, std::shared_ptr<AttributeBase>> m_attr_hash_container;

};



}  // namespace gpmesh