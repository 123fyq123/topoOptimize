#pragma once

#include <assert.h>
#include <utility>
#include "Util/GP_Macros.h"
#include "GP_Types.h"

#include <tbb/parallel_for.h>

class GPMeshTest;


namespace gpmesh {

    template<typename T>
    class SplineAttributeCPU {
    private:
        //default is shallow copy
        char *m_name;


        //local_attr (per geoemtry)
        T **m_h_local_attr;
        T *m_h_global_attr;

        layoutT m_layout = 100;

    public:
        //global_attr (require mapping information)

        uint32_t m_num_geometry = INVALID32;
        uint32_t m_num_attributes = INVALID32;
        uint16_t m_cpoint_num_per_geometry = INVALID16;
        uint32_t m_num_cpoints = INVALID32;

        constexpr static uint32_t m_block_size = 256;

    public:

        SplineAttributeCPU(const char *name)
                : m_name(nullptr),
                  m_h_local_attr(nullptr),
                  m_h_global_attr(nullptr) {
            if (name != nullptr) {
                //strlen`count exclude \0
                this->m_name = (char *) malloc(sizeof(char) * (strlen(name) + 1));
                strcpy(this->m_name, name);
            }
        }

        ~SplineAttributeCPU() {
            release_local_part();
            release_global_part();
        }

        const char *get_name() const {
            return m_name;
        }

        void allocate() {
            allocate_global_part();
            allocate_local_part();
        }

        void release() {
            release_global_part();
            release_local_part();
        }

        void load_global_part(Eigen::MatrixX<T> &global_attr) {
            assert(global_attr.rows() == m_num_cpoints * m_num_attributes);
            assert(0);
        }

        void move_local_to_global(const vector<vector<int>> &g2l) {

            assert(m_num_geometry != 0 && m_num_geometry != INVALID32);
            assert(m_cpoint_num_per_geometry != 0 && m_cpoint_num_per_geometry != INVALID16);
            assert(m_num_attributes != 0 && m_num_attributes != INVALID32);
            assert(m_num_cpoints != 0 && m_num_cpoints != INVALID32);
            assert(g2l.size() == m_num_cpoints);

            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_cpoints),
                              [&](tbb::blocked_range<size_t> r) {
                                  for (size_t global_id = r.begin(); global_id < r.end(); ++global_id) {
                                      uint32_t num_sub = g2l[global_id].size();
                                      for (uint32_t sub = 0; sub < num_sub; ++sub) {
                                          uint32_t patch_id = g2l[global_id][sub] / m_cpoint_num_per_geometry;
                                          uint16_t local_id =
                                                  g2l[global_id][sub] - m_cpoint_num_per_geometry * patch_id;
                                          assert(local_id >= 0 && local_id < m_cpoint_num_per_geometry);
                                          assert(patch_id < m_num_geometry);

                                          for (uint32_t d = 0; d < m_num_attributes; ++d) {
                                              m_h_global_attr[global_id + d * m_num_cpoints] +=
                                                      m_h_local_attr[patch_id][local_id +
                                                                               d * m_cpoint_num_per_geometry];
                                          }
                                      }
                                  }
                              });
        }

        void move_local_to_global_serial(const vector<vector<int>> &g2l)
        {
            for (size_t global_id = 0; global_id < m_num_cpoints; ++global_id) {
                uint32_t num_sub = g2l[global_id].size();
                for (uint32_t sub = 0; sub < num_sub; ++sub) {
                    uint32_t patch_id = g2l[global_id][sub] / m_cpoint_num_per_geometry;
                    uint16_t local_id =
                            g2l[global_id][sub] - m_cpoint_num_per_geometry * patch_id;
                    assert(local_id >= 0 && local_id < m_cpoint_num_per_geometry);
                    assert(patch_id < m_num_geometry);

                    for (uint32_t d = 0; d < m_num_attributes; ++d) {
                        m_h_global_attr[global_id + d * m_num_cpoints] +=
                                m_h_local_attr[patch_id][local_id +
                                                         d * m_cpoint_num_per_geometry];
                    }
                }
            }
        }

        void move_global_to_local(const vector<vector<int>> &l2g) {
            assert(m_num_geometry != 0 && m_num_geometry != INVALID32);
            assert(m_cpoint_num_per_geometry != 0 && m_cpoint_num_per_geometry != INVALID16);
            assert(m_num_attributes != 0 && m_num_attributes != INVALID32);
            assert(m_num_cpoints != 0 && m_num_cpoints != INVALID32);
            assert(l2g.size() == m_num_geometry);
            assert(l2g[0].size() == 64);

            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_geometry),
                              [&](tbb::blocked_range<size_t> r) {
                                  for (size_t geometry_id = r.begin(); geometry_id < r.end(); ++geometry_id) {
                                      for (uint16_t local_id = 0; local_id < m_cpoint_num_per_geometry; ++local_id) {
                                          uint32_t global_id = l2g[geometry_id][local_id];

                                          for (uint32_t d = 0; d < m_num_attributes; ++d) {
                                              m_h_local_attr[geometry_id][local_id + d * m_cpoint_num_per_geometry] =
                                                      m_h_global_attr[global_id + d * m_num_cpoints];
                                          }
                                      }
                                  }
                              }
            );
        }

        void move_global_to_local_serial(const vector<vector<int>> &l2g) {
            for (size_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
                for (uint16_t local_id = 0; local_id < m_cpoint_num_per_geometry; ++local_id) {
                    uint32_t global_id = l2g[geometry_id][local_id];

                    for (uint32_t d = 0; d < m_num_attributes; ++d) {
                        m_h_local_attr[geometry_id][local_id + d * m_cpoint_num_per_geometry] =
                                m_h_global_attr[global_id + d * m_num_cpoints];
                    }
                }
            }
        }

        void init_global_part(const uint32_t num_attributes,
                              const uint32_t num_cpoints) {

            release_global_part();
            m_num_attributes = num_attributes;
            m_num_cpoints = num_cpoints;
            allocate_global_part();
        }

        void release_global_part() {
            if (m_h_global_attr != nullptr) {
                delete[] m_h_global_attr;
                m_h_global_attr = nullptr;
            }
        }

        void allocate_global_part() {
            m_h_global_attr = new T[m_num_attributes * m_num_cpoints]();
        }

        void reset_global_part(const T value) {

            if (m_h_global_attr == nullptr) {
                throw std::runtime_error("m_h_global_attr is nullptr");
            }
            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_cpoints * m_num_attributes),
                              [&](tbb::blocked_range<size_t> r) {
                                  for (size_t i = r.begin(); i < r.end(); ++i) {
                                      m_h_global_attr[i] = value;
                                  }
                              }
            );
        }

        void init_local_part(const uint32_t num_attributes,
                             const uint32_t num_geometry,
                             const uint16_t cpoint_num_per_geometry = 64,
                             const layoutT layout = SoA) {
            release_local_part();
            m_num_geometry = num_geometry;
            m_num_attributes = num_attributes;
            m_cpoint_num_per_geometry = cpoint_num_per_geometry;
            m_layout = layout;

            if (m_num_geometry == 0) {
                return;
            }

            allocate_local_part();
        }

        void release_local_part() {
            if (m_h_local_attr == nullptr) return;

            for (uint32_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
                if (m_h_local_attr[geometry_id] != nullptr) {
                    delete[] m_h_local_attr[geometry_id];
                }

            }

            delete[] m_h_local_attr;
            m_h_local_attr = nullptr;
        }

        void allocate_local_part() {
            assert(m_num_geometry != 0 && m_num_geometry != INVALID32);
            assert(m_cpoint_num_per_geometry != 0 && m_cpoint_num_per_geometry != INVALID16);
            assert(m_num_attributes != 0 && m_num_attributes != INVALID32);

            m_h_local_attr = new T *[m_num_geometry]();

            for (uint32_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
                m_h_local_attr[geometry_id] = new T[m_cpoint_num_per_geometry * m_num_attributes]();
            }
        }

        void reset_local_part(const T value) {

            if (m_h_local_attr == nullptr) {
                throw std::runtime_error("m_h_local_attr is nullptr");
            }
            tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_geometry),
                              [&](tbb::blocked_range<size_t> r) {
                                  for (size_t i = r.begin(); i < r.end(); ++i) {
                                      for (uint32_t e = 0; e < m_cpoint_num_per_geometry * m_num_attributes; ++e) {
                                          m_h_local_attr[i][e] = value;
                                      }
                                  }
                              }
            );
        }

        __forceinline__ bool is_empty() const {
            return m_num_geometry == 0;
        }

        __forceinline__ T &RawAccessLocal(const uint32_t geometry_id, const uint16_t i) const {
            return m_h_local_attr[geometry_id][i];
        }

        __forceinline__ T &RawAccessGlobal(uint32_t i) {
            return m_h_global_attr[i];
        }

        __forceinline__ T &operator()(const uint32_t geometry_id,
                                      const uint16_t local_id,
                                      const uint32_t attr) const {
            assert(geometry_id < m_num_geometry);
            assert(attr < m_num_attributes);

            const uint32_t pitch_x = (m_layout == AoS) ? m_num_attributes : 1;

            const uint32_t pitch_y = (m_layout == AoS) ? 1 : m_cpoint_num_per_geometry;
            return m_h_local_attr[geometry_id][local_id * pitch_x + attr * pitch_y];
        }

        __forceinline__ T &operator()(const uint32_t cpoint_id, const uint32_t attr_id) const {
            return m_h_global_attr[cpoint_id + attr_id * m_num_cpoints];
        }


    };

    template<typename T>
    class SplineAttributeContainerCPU {
    public:
        /**
         * @brief Default constructor
         */
        SplineAttributeContainerCPU() = default;

        /**
         * @brief Destructor which releases all attribute managed by this container
         */
        virtual ~SplineAttributeContainerCPU() {
            while (!m_attr_container.empty()) {
                m_attr_container.back()->release();
                m_attr_container.pop_back();
            }

            for (auto iter = m_attr_hash_container.begin(); iter != m_attr_hash_container.end(); ++iter) {
                (*iter).second->release();
            }
        }

        /**
         * @brief Number of attribute managed by this container
         */
        size_t size() {
            //return m_attr_container.size();
            return m_attr_hash_container.size();
        }

        /**
         * @brief get a list of name of the attributes managed by this container
         * @return
         */
        std::vector<std::string> get_attribute_names() const {
            std::vector<std::string> names;
//        for (size_t i = 0; i < m_attr_container.size(); ++i) {
//            names.push_back(m_attr_container[i]->get_name());
//        }

            for (auto iter = m_attr_hash_container.begin(); iter != m_attr_hash_container.end(); ++iter) {
                names.emplace_back((*iter).first);
            }
            return names;
        }

        std::shared_ptr<SplineAttributeCPU<T>> add_spline_attr(
                const char *name,
                uint32_t
                num_geometry,
                uint16_t
                cpoint_num_per_geometry,
                uint32_t
                num_attributes,
                uint32_t num_cpoints,
                layoutT layout
        ) {
            if (does_exist(name)) {
//                GPMESH_WARN(
//                        "AttributeContainer::add() adding an attribute with "
//                        "name {} already exists!",
//                        std::string(name));
//
//                GPMESH_WARN(
//                        "this is a danger operation, should be cautious!"
//                );

                return m_attr_hash_container[name];
            }

            auto new_attr = std::make_shared<SplineAttributeCPU<T>>(name);
            new_attr->init_local_part(num_attributes, num_geometry, cpoint_num_per_geometry, layout);
            new_attr->init_global_part(num_attributes, num_cpoints);

            std::string key = name;

            m_attr_hash_container[key] = new_attr;

            return new_attr;
        }

        /**
         * @brief Check if an attribute exists
         * @param name of the attribute
         */
        bool does_exist(const char *name) {

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
        void remove(const char *name) {
            if (does_exist(name)) {
                string key = name;
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

        std::shared_ptr<SplineAttributeCPU<T>> get_attr(string name) {
            if (m_attr_hash_container.find(name) == m_attr_hash_container.end())
                throw std::runtime_error("no attr find in attribute container");

            return m_attr_hash_container[name];
        }

    private:
        std::vector<std::shared_ptr<SplineAttributeCPU<T>>> m_attr_container;
        std::unordered_map<std::string, std::shared_ptr<SplineAttributeCPU<T>>> m_attr_hash_container;

    };


}  // namespace gpmesh