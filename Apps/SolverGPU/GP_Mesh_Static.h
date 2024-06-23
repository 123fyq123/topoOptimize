#pragma once
#include <assert.h>
#include <fstream>
#include <memory>
#include <thread>

#include <cuda_profiler_api.h>

#include "GP_Attribute.h"
#include "GP_Handle.h"
#include "Kernels/GP_For_Each.cuh"
#include "GP_Launch_Box.h"
#include "GP_Mesh.h"
#include "GP_Types.h"
#include "Util/GP_Import_Obj.h"
#include "Util/GP_Log.h"
#include "Util/GP_Timer.h"

#include <filesystem>

#include <thread>

#pragma once
#include <assert.h>
#include <fstream>
#include <memory>


//template include for debug
#include "GP_vertex_normal_kernel.cuh"
#include "GP_vertex_normal_ref.h"
#include "GP_shallow_wave_kernel.cuh"
#include "GP_mean_curvature_flow_kernel.cuh"
#include "GP_cloth_debug_kernel.cuh"
#include "GP_debug_kernel.cuh"

#include "Util/GP_Report.h"
#include "GP_Debug.h"

//MeshKernel
#include "gpmesh_MeshKernel_helper.h"

namespace gpmesh {

/**
 * @brief This class is responsible for query operations of static meshes. It
 * extends GPMesh with methods needed to launch kernel and do computation on the
 * mesh as well as managing mesh attributes
 */
class GPMeshStatic : public GPMesh
{
   public:
    GPMeshStatic(const GPMeshStatic&) = delete;

    // debug part
    explicit GPMeshStatic(const bool quite = false) : GPMesh()
    {
    }

    void testOutputSinglePatchWithRibbon(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices, std::string& file_path, uint16_t patch_id)
    {
        std::string output_name = file_path + "_" + std::to_string(patch_id) + ".obj";

        // patch start and end
        const uint32_t p_start =
            (patch_id == 0) ? 0 : m_patcher->get_patches_offset()[patch_id - 1];
        const uint32_t p_end = m_patcher->get_patches_offset()[patch_id];

        // ribbon start and end
        const uint32_t r_start =
            (patch_id == 0) ? 0 :
                              m_patcher->get_external_ribbon_offset()[patch_id - 1];
        const uint32_t r_end = m_patcher->get_external_ribbon_offset()[patch_id];


        const uint32_t total_patch_num_faces =
            (p_end - p_start) + (r_end - r_start);

        int n_local_vh = 0;
        std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
        std::vector<std::vector<float>> local_vertices;

        for (uint32_t f = p_start; f < p_end; ++f) {
            uint32_t face_id = m_patcher->get_patches_val()[f];
            auto vhs = fv[face_id];

            //Todo: read vh once per time, and add to hashTable and local_vertices
            for (int j = 0; j < vhs.size(); ++j) {
                auto global_vh = vhs[j];
                if (m_hash_g2l.find(global_vh) == m_hash_g2l.end()) {
                    m_hash_g2l[global_vh] = n_local_vh++;
                    local_vertices.push_back(vertices[global_vh]);
                }
            }
        }

        for (uint32_t f = r_start; f < r_end; ++f) {
            uint32_t face_id = m_patcher->get_external_ribbon_val()[f];

            auto vhs = fv[face_id];

            //Todo: read vh once per time, and add to hashTable and local_vertices
            for (int j = 0; j < vhs.size(); ++j) {
                auto global_vh = vhs[j];
                if (m_hash_g2l.find(global_vh) == m_hash_g2l.end()) {
                    m_hash_g2l[global_vh] = n_local_vh++;
                    local_vertices.push_back(vertices[global_vh]);
                }
            }
        }


        //Todo: Output Patch File
        std::ofstream output(output_name);

        for (int j = 0 ; j < local_vertices.size(); ++j)
        {
            output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
        }


        for (uint32_t f = p_start; f < p_end; ++f) {
            uint32_t face_id = m_patcher->get_patches_val()[f];
            auto& vhs = fv[face_id];

            output << "f" ;
            for (auto vh : vhs)
            {
                output << " " << m_hash_g2l[vh] + 1;
            }

            output << std::endl;
        }

        for (uint32_t f = r_start; f < r_end; ++f) {
            uint32_t face_id = m_patcher->get_external_ribbon_val()[f];

            auto& vhs = fv[face_id];

            output << "f" ;
            for (auto vh : vhs)
            {
                output << " " << m_hash_g2l[vh] + 1;
            }

            output << std::endl;
        }

        output.close();
    }

    void testOutputPatcherWithRibbonParallel(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices, std::string& file_path)
    {
        std::vector<std::thread> threads_pool;

        for (int i = 0; i < m_patcher->get_num_patches(); ++i)
        {
            threads_pool.push_back(std::thread(&GPMeshStatic::testOutputSinglePatchWithRibbon, this, std::ref(fv), std::ref(vertices), std::ref(file_path), i));
        }

        for (int i = 0; i < m_patcher->get_num_patches(); ++i)
        {
            threads_pool[i].join();
        }

        return;
    }

    void testOutputPatcherWithRibbon(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices, std::string& file_path)
    {
        testOutputSinglePatchWithRibbon(fv, vertices, file_path, 0);
        return;
    }

    void testOutputSinglePatchTopology(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices,
                                       std::string& file_path, uint16_t patch_id)
    {
        //m_h_patches_ltog_v: get global_vertex_index from local_vertex_index, then read vertex position from vertices

        //m_h_patches_fe, get four local_edge_index from local_face_index
        //m_h_patches_ev, get two local_vertices_index from local_edge_index
        //so, we can get local_vertices_index from local_face_index
        //then output the local_vertices_index in Obj file

        //m_h_patches_ltog_f: can use to verify the accuracy of patch topology

        int n_local_vh = 0;
        std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
        std::vector<std::vector<float>> local_vertices;

        //according to the local_vertex_id ordering
        for (auto& global_vertex_id : m_h_patches_ltog_v[patch_id])
        {
            local_vertices.push_back(vertices[global_vertex_id]);
        }


        int patch_num_faces = m_h_patches_ltog_f[patch_id].size();

        std::vector<std::vector<uint32_t>> local_faces;
        local_faces.resize(patch_num_faces);

        for(int local_face_id = 0; local_face_id < patch_num_faces; ++local_face_id)
        {
            std::vector<uint32_t> local_face;
            for (int j = 0; j < fv[0].size(); ++j) {
                uint32_t local_edge_id =
                    m_h_patches_fe[patch_id][local_face_id * 4 + j];
                int dir       = local_edge_id & 1;
                local_edge_id = local_edge_id >> 1;

                uint32_t offset;
                if (dir == 1)
                {
                    offset = 1;
                }
                else
                {
                    offset = 0;
                }
                auto ek = std::make_pair(m_h_patches_ev[patch_id][local_edge_id * 2 + 0], m_h_patches_ev[patch_id][local_edge_id * 2 + 1]);
                local_face.push_back(m_h_patches_ev[patch_id][local_edge_id*2 + offset]);
            }
            local_faces[local_face_id] = local_face;
        }



        std::string output_name = file_path + "_" + std::to_string(patch_id) + ".obj";


        //Todo: Output Patch File
        std::ofstream output(output_name);

        for (int j = 0 ; j < local_vertices.size(); ++j)
        {
            output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
        }


        for (uint32_t j = 0; j < local_faces.size(); ++j) {
            auto& vhs = local_faces[j];

            output << "f" ;

            for (auto vh : vhs)
            {
                output << " " << vh + 1;
            }

            output << std::endl;
        }

        output.close();
    }

    void testOutputPathNotOwnedRing(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices,
                                    std::string& file_path, uint16_t patch_id)
    {
        int n_local_vh = 0;
        std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
        std::vector<std::vector<float>> local_vertices;

        uint16_t patch_f_num = m_h_patches_ltog_f[patch_id].size();
        uint16_t patch_owned_f_num = m_h_num_owned_f[patch_id];
        uint16_t patch_not_owned_f_num = patch_f_num - patch_owned_f_num;

        for (int i = 0; i < patch_not_owned_f_num; ++i)
        {
            uint16_t local_id     = i + m_h_num_owned_f[patch_id];
            uint32_t global_id    = m_h_patches_ltog_f[patch_id][local_id];

            uint32_t owning_patch = m_patcher->get_face_patch_id(global_id);
            //confirm the ribbon face
            bool isFind = false;
            for (auto global_f_id_tmp : m_h_patches_ltog_f[owning_patch])
            {
                if (global_f_id_tmp == global_id)
                {
                    isFind = !isFind;
                }
            }

            if (!isFind)
            {
                GPMESH_ERROR("ribbon face from patch-{} is not found in patch-{}", patch_id, owning_patch);
            }
        }

        //output not owned patch id
        int count = 0;

        std::vector<std::vector<uint32_t>> local_faces;
        local_faces.resize(patch_not_owned_f_num);

        for (int i = 0; i < patch_not_owned_f_num; ++i)
        {
            uint16_t patch_not_owned_fid =  m_h_patches_info[patch_id].not_owned_id_f[i].id;
            uint32_t patch_not_owned_f_pid = m_h_patches_info[patch_id].not_owned_patch_f[i];

            uint32_t global_fid = m_h_patches_ltog_f[patch_not_owned_f_pid][patch_not_owned_fid];

            auto& local_face = local_faces[i];

            for (auto global_vid : fv[global_fid])
            {
                if (m_hash_g2l.find(global_vid) == m_hash_g2l.end())
                {
                    m_hash_g2l[global_vid] = count++;
                }

                local_face.push_back(m_hash_g2l[global_vid]);
            }
        }

        local_vertices = std::vector<std::vector<float>> (m_hash_g2l.size(), std::vector<float>(3, -1));

        for(auto& gl_pair : m_hash_g2l)
        {
            for (int i = 0; i < 3; ++i)
            {
                local_vertices[gl_pair.second][i] = vertices[gl_pair.first][i];
            }
        }

        std::string output_name = file_path + "_" + std::to_string(patch_id) + ".obj";

        //Todo: Output Patch File
        std::ofstream output(output_name);

        for (int j = 0 ; j < local_vertices.size(); ++j)
        {
            output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
        }


        for (uint32_t j = 0; j < local_faces.size(); ++j) {
            auto& vhs = local_faces[j];

            output << "f" ;

            for (auto vh : vhs)
            {
                output << " " << vh + 1;
            }

            output << std::endl;
        }

        output.close();
    }

    void testOutputPatchTopology(std::vector<std::vector<uint32_t>>& fv, std::vector<std::vector<float>>& vertices, std::string& file_path)
    {
        for (int i = 0; i < m_num_patches; ++i)
        {
            testOutputSinglePatchTopology(fv, vertices, file_path, i);
        }
        return;
    }

    std::vector<std::vector<float>> m_vertices;
    std::vector<std::vector<uint32_t>> m_fv;

    void debug()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = "/run/media/whalesea/Projects/ManjaroProjects/RXMesh/input/QuadMesh/plane-quad3.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        // Output Patchers in Obj Format for Debug
        bool testOutput = true;
        if(testOutput)
        {
            std::vector<uint32_t>& face_patch = m_patcher->get_face_patch();
            size_t num_patch = m_patcher->get_num_patches();
            //output different patch into obj format
            std::vector<std::vector<uint32_t>> faces_per_patch;
            faces_per_patch.resize(num_patch);
            for (int i = 0; i < face_patch.size(); ++i)
            {
                faces_per_patch[face_patch[i]].push_back(i);
            }

            for (int i = 0; i < num_patch; ++i)
            {
                std::string base_file_name=  file_path.substr(file_path.find_last_of("/")+1, file_path.find_last_of('.') - file_path.find_last_of('/') - 1);
                std::string output_file_name = "Patches/patch_" + std::to_string(i) + "_from_" + base_file_name + ".obj";

                //we have fv and vertices
                int n_local_vh = 0;
                std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
                std::vector<std::vector<float>> local_vertices;

                auto local_faces = faces_per_patch[i];

                for (auto local_face : local_faces)
                {
                    auto vhs = fv[local_face];

                    //Todo: read vh once per time, and add to hashTable and local_vertices
                    for (int j = 0; j < vhs.size(); ++j)
                    {
                        auto global_vh = vhs[j];
                        if (m_hash_g2l.find(global_vh) == m_hash_g2l.end())
                        {
                            m_hash_g2l[global_vh] = n_local_vh++;
                            local_vertices.push_back(vertices[global_vh]);
                        }
                    }
                }

                //Todo: Output Patch File
                std::ofstream output(output_file_name);

                for (int j = 0 ; j < local_vertices.size(); ++j)
                {
                    output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
                }

                for (auto local_face : local_faces)
                {
                    output << "f" ;
                    for (auto vh : fv[local_face])
                    {
                        output << " " << m_hash_g2l[vh] + 1;
                    }

                    output << std::endl;
                }

                output.close();
            }

            std::string output_file_path = "Patches/patch_with_ribbon_from_kiten_quad";
            testOutputPatcherWithRibbonParallel(fv, vertices, output_file_path);
        }

        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        //bool test_build_patch_output = true;
        if (testOutput)
        {
            std::string output_file_path = "PatchesTopologyDebug/patch_topology_debug_from_kiten_quad";
            testOutputPatchTopology(fv, vertices, output_file_path);
        }

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);
        
        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        auto coords = add_vertex_attribute<float>(vertices, "coordinates");

        auto v_normals = add_vertex_attribute<float>("v_normals", 3, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box;

        // compute<float, blockThreads><<<>>>
        // calculate shared memory required for each block
        prepare_launch_box({gpmesh::Op::FV},
                           launch_box,
                           (void*) compute_vertex_normal<float, blockThreads>
                           );

        TestData testData;
        testData.test_name = "VertexNormal";
        testData.num_threads = launch_box.num_threads;
        testData.num_blocks  = launch_box.blocks;
        testData.dyn_smem    = launch_box.smem_bytes_dyn;
        testData.static_smem = launch_box.smem_bytes_static;
        testData.num_reg     = launch_box.num_registers_per_thread;

        float vn_time = 0;

        int num_runs = 1;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            v_normals->reset(0, gpmesh::DEVICE);

            GPUTimer timer;

            timer.start();

            compute_vertex_normal<float, blockThreads><<<launch_box.blocks,
                                                         launch_box.num_threads,
                                                         launch_box.smem_bytes_dyn>>>
                (get_context(), *coords, *v_normals);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            testData.time_ms.push_back(timer.elapsed_millis());
            vn_time += timer.elapsed_millis();
        }

        v_normals->move(gpmesh::DEVICE, gpmesh::HOST);

        GPMESH_TRACE("vertex_normal_gpmesh vertex normal kernel took {} (ms)",
                     vn_time / num_runs);

        //calculate vertex normal ref in CPU

        std::vector<float> vertices_normal_ref(3*m_num_vertices, 0);

        vertex_normal_ref(fv, vertices, vertices_normal_ref);


        //可视化法向:

        std::vector<std::vector<float>> vertices_color = std::vector<std::vector<float>>(m_num_vertices, std::vector<float>(3, -1));

        std::vector<float> vertices_normal_gpu(3*m_num_vertices, 0);

        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            for (uint32_t i = 0; i < 3; ++i) {
                //assert(std::abs(std::abs(vertices_normal_ref[v_id * 3 + i]) - std::abs((*v_normals)(vh, i))) < 0.0001);
                vertices_color[v_id][i] = ((*v_normals)(vh, i) + 4.0f) / 6.0f;
                vertices_normal_gpu[3 * v_id + i] = (*v_normals)(vh, i);
            }
            //GPMESH_TRACE("v_id {} - normal ({}, {}, {})", v_id, (*v_normals)(vh, 0), (*v_normals)(vh, 1), (*v_normals)(vh, 2));
        });

        for (int i = 0; i < m_num_vertices; ++i)
        {
            bool flagSame = true;

            for (int j = 0; j < 3; ++j)
            {
                if (std::abs(vertices_normal_ref[3 * i + j] - vertices_normal_gpu[3 * i + j]) >= 0.1)
                {
                    flagSame = false;
                }
            }

            if (flagSame)
            {
                //GPMESH_TRACE("vid {} has same normal in gpu and cpu", i);
            }
        }

        bool use_cpu_normal = true;
        bool use_gpu_normal = true;

        if(use_cpu_normal)
        {
            outputVerticesInfo2File<float, 3>("NormalDebug/vertices_normal_cpu.txt", vertices_normal_ref);

            if(testOutput)
            {

                std::vector<std::vector<float>> vertices_color_cpu = std::vector<std::vector<float>>(m_num_vertices, std::vector<float>(3, -1));

                for (int i = 0; i < m_num_vertices; ++i)
                {
                    for (int j = 0; j < 3; ++j)
                    {
                        vertices_color_cpu[i][j] = (vertices_normal_ref[3 * i + j] + 4.0f) / 6.0f;
                    }
                }

                std::string output_filename = "NormalDebug/NormalColor_CPU.obj";

                std::ofstream output(output_filename);

                for (int i = 0; i < m_num_vertices; ++i)
                {
                    output << "v";

                    for (int j = 0; j < 3; ++j)
                    {
                        output << " " << vertices[i][j];
                    }

                    for (int j = 0; j < 3; ++j)
                    {
                        output << " " << vertices_color_cpu[i][j];
                    }

                    output << std::endl;
                }

                for (int i = 0; i < m_num_faces; ++i)
                {
                    output << "f";

                    for (int j = 0; j < 4; ++j)
                    {
                        output << " " << fv[i][j] + 1;
                    }

                    output << std::endl;
                }

                output.close();

                GPMESH_TRACE("Output obj with cpu normal as color");
            }
        }

        //output Obj mesh with vertices color
        if (testOutput && use_gpu_normal)
        {
            outputVerticesInfo2File<float, 3>("NormalDebug/vertices_normal_gpu.txt", vertices_normal_gpu);

            std::string output_filename = "NormalDebug/NormalColor_GPU.obj";

            std::ofstream output(output_filename);

            for (int i = 0; i < m_num_vertices; ++i)
            {
                output << "v";

                for (int j = 0; j < 3; ++j)
                {
                    output << " " << vertices[i][j];
                }

                for (int j = 0; j < 3; ++j)
                {
                    output << " " << vertices_color[i][j];
                }

                output << std::endl;
            }

            for (int i = 0; i < m_num_faces; ++i)
            {
                output << "f";

                for (int j = 0; j < 4; ++j)
                {
                    output << " " << fv[i][j] + 1;
                }

                output << std::endl;
            }

            output.close();

            GPMESH_TRACE("Output obj with normal as color");
        }
    }

    void debug_neighborSize()
    {
        auto coords = add_vertex_attribute<float>(m_vertices, "coordinates for neighbor size");
        auto v_neighborSize = add_vertex_attribute<float>("v_neighborSize", 1, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box;

        prepare_launch_box({Op::FV},
                           launch_box,
                           (void *)compute_vertex_neighbor_size<float, blockThreads>);


        TestData testData;
        testData.test_name = "Neighbor Size";
        testData.num_threads = launch_box.num_threads;
        testData.num_blocks  = launch_box.blocks;
        testData.dyn_smem    = launch_box.smem_bytes_dyn;
        testData.static_smem = launch_box.smem_bytes_static;
        testData.num_reg     = launch_box.num_registers_per_thread;

        v_neighborSize->reset(0, gpmesh::DEVICE);


        float ns_time = 0.0f;

        GPUTimer timer;

        timer.start();

        compute_vertex_neighbor_size<float, blockThreads><<<launch_box.blocks,
                                                            launch_box.num_threads,
                                                            launch_box.smem_bytes_dyn>>>
            (get_context(), *coords, *v_neighborSize);

        timer.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        testData.time_ms.push_back(timer.elapsed_millis());

        ns_time += timer.elapsed_millis();

        GPMESH_TRACE("compute_vertex_neighbor_size kenrel took {}ms",
                     ns_time);

        v_neighborSize->move(gpmesh::DEVICE, gpmesh::HOST);

        //move data from attribute to vector
        std::vector<float> vertices_neighbor_size_gpu(m_num_vertices,  0);
        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            vertices_neighbor_size_gpu[v_id] = (*v_neighborSize)(vh, 0);
        });

        // calculate vertex length ref in cpu;
        std::vector<float> vertices_neighbor_size_cpu(m_num_vertices,  0);

        vertex_neighbor_sizes_ref(m_fv, m_vertices, vertices_neighbor_size_cpu);

        for (int i = 0; i < vertices_neighbor_size_gpu.size(); ++i)
        {
            float diff = std::abs(vertices_neighbor_size_cpu[i] - vertices_neighbor_size_gpu[i]);
            if ( diff >= 0.00001)
            {
                GPMESH_ERROR("Validation Error for vid-{} with diff {}", i, diff);
            }
        }

        ns_time = 0.0f;

        auto v_fvids = add_vertex_attribute<uint32_t>("v_fvids", 1, gpmesh::LOCATION_ALL);

        LaunchBox<blockThreads> launch_box_fv;

        prepare_launch_box({Op::FV},
                           launch_box_fv,
                           (void *)compute_face_vertex_order<uint32_t , blockThreads>);


        TestData testData_fv;
        testData_fv.test_name = "Topology Info";
        testData_fv.num_threads = launch_box.num_threads;
        testData_fv.num_blocks  = launch_box.blocks;
        testData_fv.dyn_smem    = launch_box.smem_bytes_dyn;
        testData_fv.static_smem = launch_box.smem_bytes_static;
        testData_fv.num_reg     = launch_box.num_registers_per_thread;

        v_fvids->reset(0, gpmesh::DEVICE);


        float ns_time_fv = 0.0f;

        GPUTimer timer_fv;

        timer_fv.start();

        compute_face_vertex_order<uint32_t, blockThreads><<<launch_box_fv.blocks,
                                                            launch_box_fv.num_threads,
                                                            launch_box_fv.smem_bytes_dyn
                                                            >>>(get_context(), *v_fvids);

        timer_fv.stop();

        CUDA_ERROR(cudaDeviceSynchronize());;
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        testData_fv.time_ms.push_back(timer_fv.elapsed_millis());

        ns_time += timer.elapsed_millis();

        GPMESH_TRACE("compute_vertex_fvids kenrel took {}ms",
                     ns_time);

        v_fvids->move(gpmesh::DEVICE, gpmesh::HOST);

        //calculate info in cpu
        std::vector<uint32_t> vertices_fvids(m_num_vertices,  0);

        for (int i = 0; i < m_num_vertices; ++i)
        {
            auto patch_id = m_patcher->get_vertex_patch_id(i);

            uint32_t local_id = INVALID32;

            for (int j = 0; j < m_h_patches_ltog_v[patch_id].size(); ++j)
            {
                if (m_h_patches_ltog_v[patch_id][j] == i)
                {
                    local_id = j;
                    break;
                }
            }

            if (local_id == INVALID32)
            {
                GPMESH_ERROR("global vid-{} is not found in patch-{}", i, patch_id);
                break;
            }
            else
            {
                vertices_fvids[i] += patch_id * 100000 + local_id;
            }
        }

        bool successFlag4fvid = true;

        std::vector<uint32_t> vertices_neighbor_face_count(m_num_vertices, 0);
        for (auto f : m_fv)
        {
            for (auto f_vid : f)
            {
                vertices_neighbor_face_count[f_vid]++;
            }
        }

        for_each_vertex(gpmesh::HOST,[&](const VertexHandle& vh) {
            uint32_t global_v_id = map_to_global(vh);

            if (vertices_fvids[global_v_id] != (*v_fvids)(vh, 0) / vertices_neighbor_face_count[global_v_id])
            {
                GPMESH_ERROR("Validation Error for global vid-{}, with cpu value {} against gpu value {}", global_v_id, vertices_fvids[global_v_id], (*v_fvids)(vh, 0) / 4);
                successFlag4fvid = false;
            }
            else
            {
                //GPMESH_TRACE("Validation Success for global vid-{}, with cpu value {} against gpu value {}", global_v_id, vertices_fvids[global_v_id], (*v_fvids)(vh, 0) / 4);
            }
        });

        if (successFlag4fvid)
        {
            GPMESH_INFO("Validation Success for fvid check!");
        }
        else
        {
            GPMESH_ERROR("Validation Error for fvid check!");
        }

        ns_time = 0.0f;

        LaunchBox<blockThreads> launch_box_fvids;

        prepare_launch_box({Op::FV},
                           launch_box_fvids,
                           (void *)compute_fvids<uint32_t , blockThreads>);

        auto f_vids = add_face_attribute<uint32_t>("f_vids", 4, gpmesh::LOCATION_ALL);

        TestData testData_fvids;
        testData_fvids.test_name = "Topology Info";
        testData_fvids.num_threads = launch_box.num_threads;
        testData_fvids.num_blocks  = launch_box.blocks;
        testData_fvids.dyn_smem    = launch_box.smem_bytes_dyn;
        testData_fvids.static_smem = launch_box.smem_bytes_static;
        testData_fvids.num_reg     = launch_box.num_registers_per_thread;

        f_vids->reset(0, gpmesh::DEVICE);


        float ns_time_fvids = 0.0f;

        GPUTimer timer_fvids;

        timer_fvids.start();

        compute_fvids<uint32_t, blockThreads><<<launch_box_fvids.blocks,
                                                launch_box_fvids.num_threads,
                                                launch_box_fvids.smem_bytes_dyn>>>(
            get_context(),
            *f_vids
            );

        timer_fvids.stop();

        CUDA_ERROR(cudaDeviceSynchronize());;
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        testData_fvids.time_ms.push_back(timer_fvids.elapsed_millis());

        ns_time += timer.elapsed_millis();

        GPMESH_TRACE("compute_fvids_kenrel took {}ms",
                     ns_time);


        f_vids->move(gpmesh::DEVICE, gpmesh::HOST);
        bool successFvidCheckFlag = true;
        for_each_face(gpmesh::HOST, [&](const FaceHandle& fh) -> void
                      {
                          uint32_t global_f_id = map_to_global(fh);
                          bool local_success_flag = true;
                          for (int i = 0; i < 4; ++i)
                          {
                            if(vertices_fvids[m_fv[global_f_id][i]] != (*f_vids)(fh, i)){
                                local_success_flag = false;
                            }
                          }

                          if(local_success_flag == false)
                          {
                              GPMESH_ERROR("Validation Error for fid-{}", global_f_id);
                              for (int i = 0; i < 4; ++i)
                              {
                                  GPMESH_ERROR("#### {} - {}", vertices_fvids[m_fv[global_f_id][i]], (*f_vids)(fh, i));
                              }
                              successFvidCheckFlag = false;
                          }
                      });

        if (successFvidCheckFlag)
        {
            GPMESH_INFO("fvid check success!");
        }
        else
        {
            GPMESH_ERROR("fvid check failed !!!");
        }
    }

    void debug_check_patch_info()
    {
        auto coords = add_vertex_attribute<float>(m_vertices, "coordinates for neighbor size");
        auto v_topology_info = add_vertex_attribute<float>("v_topology_info", 1, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box;

        prepare_launch_box({Op::FV},
                           launch_box,
                           (void *)compute_vertex_topology<float, blockThreads>);


        TestData testData;
        testData.test_name = "Topology Info";
        testData.num_threads = launch_box.num_threads;
        testData.num_blocks  = launch_box.blocks;
        testData.dyn_smem    = launch_box.smem_bytes_dyn;
        testData.static_smem = launch_box.smem_bytes_static;
        testData.num_reg     = launch_box.num_registers_per_thread;

        v_topology_info->reset(0, gpmesh::DEVICE);


        float ns_time = 0.0f;

        GPUTimer timer;

        timer.start();

        compute_vertex_topology<float, blockThreads><<<launch_box.blocks,
                                                            launch_box.num_threads,
                                                            launch_box.smem_bytes_dyn>>>
            (get_context(), *coords, *v_topology_info);

        timer.stop();

        CUDA_ERROR(cudaDeviceSynchronize());
        CUDA_ERROR(cudaGetLastError());
        CUDA_ERROR(cudaProfilerStop());

        testData.time_ms.push_back(timer.elapsed_millis());

        ns_time += timer.elapsed_millis();

        GPMESH_TRACE("compute_vertex_neighbor_size kenrel took {}ms",
                     ns_time);

        v_topology_info->move(gpmesh::DEVICE, gpmesh::HOST);

        //move data from attribute to vector
        std::vector<float> vertices_topology_info_gpu(m_num_vertices,  0);
        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            vertices_topology_info_gpu[v_id] = (*v_topology_info)(vh, 0);
        });

        // calculate vertex length ref in cpu;
        std::vector<float> vertices_topology_info_cpu(m_num_vertices,  0);

        vertex_neighbor_sizes_ref(m_fv, m_vertices, vertices_topology_info_cpu);

        for (int i = 0; i < vertices_topology_info_gpu.size(); ++i)
        {
            float diff = std::abs(vertices_topology_info_cpu[i] - vertices_topology_info_gpu[i]);
            if ( diff >= 0.0001)
            {
                GPMESH_TRACE("Validation Error for vid-{} with diff {}", i, diff);
            }
        }
    }

    void debug_check_patch_info_v2()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = "/run/media/whalesea/Projects/ManjaroProjects/RXMesh/input/QuadMesh/dragon-quad.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        // Output Patchers in Obj Format for Debug
        bool testOutput = true;
        if(testOutput)
        {
            std::vector<uint32_t>& face_patch = m_patcher->get_face_patch();
            size_t num_patch = m_patcher->get_num_patches();
            //output different patch into obj format
            std::vector<std::vector<uint32_t>> faces_per_patch;
            faces_per_patch.resize(num_patch);
            for (int i = 0; i < face_patch.size(); ++i)
            {
                faces_per_patch[face_patch[i]].push_back(i);
            }

            for (int i = 0; i < num_patch; ++i)
            {
                std::string base_file_name=  file_path.substr(file_path.find_last_of("/")+1, file_path.find_last_of('.') - file_path.find_last_of('/') - 1);
                std::string output_file_name = "Patches/patch_" + std::to_string(i) + "_from_" + base_file_name + ".obj";

                //we have fv and vertices
                int n_local_vh = 0;
                std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
                std::vector<std::vector<float>> local_vertices;

                auto local_faces = faces_per_patch[i];

                for (auto local_face : local_faces)
                {
                    auto vhs = fv[local_face];

                    //Todo: read vh once per time, and add to hashTable and local_vertices
                    for (int j = 0; j < vhs.size(); ++j)
                    {
                        auto global_vh = vhs[j];
                        if (m_hash_g2l.find(global_vh) == m_hash_g2l.end())
                        {
                            m_hash_g2l[global_vh] = n_local_vh++;
                            local_vertices.push_back(vertices[global_vh]);
                        }
                    }
                }

                //Todo: Output Patch File
                std::ofstream output(output_file_name);

                for (int j = 0 ; j < local_vertices.size(); ++j)
                {
                    output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
                }

                for (auto local_face : local_faces)
                {
                    output << "f" ;
                    for (auto vh : fv[local_face])
                    {
                        output << " " << m_hash_g2l[vh] + 1;
                    }

                    output << std::endl;
                }

                output.close();
            }

            std::string output_file_path = "Patches/patch_with_ribbon_from_kiten_quad";
            testOutputPatcherWithRibbonParallel(fv, vertices, output_file_path);
        }

        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        //bool test_build_patch_output = true;
        if (testOutput)
        {
            std::string output_file_path = "PatchesTopologyDebug/patch_topology_debug_from_kiten_quad";
            testOutputPatchTopology(fv, vertices, output_file_path);
        }

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        bool isTestOutput = true;

        if(isTestOutput)
        {
            std::string output_not_owned = "NotOwnedPatchPart/NotOwned";
            for (int i = 0; i < m_num_patches; ++i)
            {
                testOutputPathNotOwnedRing(fv, vertices, output_not_owned, 0);
            }

            GPMESH_TRACE("Output Ribbon Region For Every Patch");
        }
    }


    void debug_normal()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = "/run/media/whalesea/Projects/ManjaroProjects/RXMesh/input/QuadMesh/plane-quad3.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        // Output Patchers in Obj Format for Debug
        bool testOutput = true;
        if(testOutput)
        {
            std::vector<uint32_t>& face_patch = m_patcher->get_face_patch();
            size_t num_patch = m_patcher->get_num_patches();
            //output different patch into obj format
            std::vector<std::vector<uint32_t>> faces_per_patch;
            faces_per_patch.resize(num_patch);
            for (int i = 0; i < face_patch.size(); ++i)
            {
                faces_per_patch[face_patch[i]].push_back(i);
            }

            for (int i = 0; i < num_patch; ++i)
            {
                std::string base_file_name=  file_path.substr(file_path.find_last_of("/")+1, file_path.find_last_of('.') - file_path.find_last_of('/') - 1);
                std::string output_file_name = "Patches/patch_" + std::to_string(i) + "_from_" + base_file_name + ".obj";

                //we have fv and vertices
                int n_local_vh = 0;
                std::unordered_map<uint32_t, uint32_t> m_hash_g2l;
                std::vector<std::vector<float>> local_vertices;

                auto local_faces = faces_per_patch[i];

                for (auto local_face : local_faces)
                {
                    auto vhs = fv[local_face];

                    //Todo: read vh once per time, and add to hashTable and local_vertices
                    for (int j = 0; j < vhs.size(); ++j)
                    {
                        auto global_vh = vhs[j];
                        if (m_hash_g2l.find(global_vh) == m_hash_g2l.end())
                        {
                            m_hash_g2l[global_vh] = n_local_vh++;
                            local_vertices.push_back(vertices[global_vh]);
                        }
                    }
                }

                //Todo: Output Patch File
                std::ofstream output(output_file_name);

                for (int j = 0 ; j < local_vertices.size(); ++j)
                {
                    output << "v " << local_vertices[j][0] << " " << local_vertices[j][1] << " " << local_vertices[j][2] << std::endl;
                }

                for (auto local_face : local_faces)
                {
                    output << "f" ;
                    for (auto vh : fv[local_face])
                    {
                        output << " " << m_hash_g2l[vh] + 1;
                    }

                    output << std::endl;
                }

                output.close();
            }

            std::string output_file_path = "Patches/patch_with_ribbon_from_kiten_quad";
            testOutputPatcherWithRibbonParallel(fv, vertices, output_file_path);
        }

        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        //bool test_build_patch_output = true;
        if (testOutput)
        {
            std::string output_file_path = "PatchesTopologyDebug/patch_topology_debug_from_kiten_quad";
            testOutputPatchTopology(fv, vertices, output_file_path);
        }

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        auto coords = add_vertex_attribute<float>(vertices, "coordinates");

        auto v_normals = add_vertex_attribute<float>("v_normals", 3, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box;

        // compute<float, blockThreads><<<>>>
        // calculate shared memory required for each block
        prepare_launch_box({gpmesh::Op::FV},
                           launch_box,
                           (void*) compute_vertex_normal_v2<float, blockThreads>
        );

        TestData testData;
        testData.test_name = "VertexNormal_v2";
        testData.num_threads = launch_box.num_threads;
        testData.num_blocks  = launch_box.blocks;
        testData.dyn_smem    = launch_box.smem_bytes_dyn;
        testData.static_smem = launch_box.smem_bytes_static;
        testData.num_reg     = launch_box.num_registers_per_thread;

        float vn_time = 0;

        int num_runs = 1;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            v_normals->reset(0, gpmesh::DEVICE);

            GPUTimer timer;

            timer.start();

            compute_vertex_normal_v2<float, blockThreads><<<launch_box.blocks,
                                                         launch_box.num_threads,
                                                         launch_box.smem_bytes_dyn>>>
                (get_context(), *coords, *v_normals);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            testData.time_ms.push_back(timer.elapsed_millis());
            vn_time += timer.elapsed_millis();
        }

        v_normals->move(gpmesh::DEVICE, gpmesh::HOST);

        GPMESH_TRACE("vertex_normal_gpmesh_v2e vertex normal kernel took {} (ms)",
                     vn_time / num_runs);

        //calculate vertex normal ref in CPU

        std::vector<float> vertices_normal_ref(3*m_num_vertices, 0);

        vertex_normal_ref_v2(fv, vertices, vertices_normal_ref);


        //可视化法向:

        std::vector<float> vertices_normal_gpu(3*m_num_vertices, 0);

        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            for (uint32_t i = 0; i < 3; ++i) {
                //assert(std::abs(std::abs(vertices_normal_ref[v_id * 3 + i]) - std::abs((*v_normals)(vh, i))) < 0.0001);
                //vertices_color[v_id][i] = ((*v_normals)(vh, i) + 4.0f) / 6.0f;
                vertices_normal_gpu[3 * v_id + i] = (*v_normals)(vh, i);
            }

            float diff = 0;

            for (int i = 0; i < 3; ++i)
            {
                float local_diff  = vertices_normal_ref[3 * v_id + i] - vertices_normal_gpu[3 * v_id +i];
                diff += local_diff * local_diff;
            }

            if(diff > 0.00000001f)
            {
                GPMESH_ERROR("NormaL for v_id-{} error : diff-{}", v_id, diff);
            }

            //GPMESH_TRACE("v_id {} - normal ({}, {}, {})", v_id, (*v_normals)(vh, 0), (*v_normals)(vh, 1), (*v_normals)(vh, 2));
        });

        for (int i = 0; i < m_num_vertices; ++i)
        {
            bool flagSame = true;

            for (int j = 0; j < 3; ++j)
            {
                if (std::abs(vertices_normal_ref[3 * i + j] - vertices_normal_gpu[3 * i + j]) >= 0.1)
                {
                    flagSame = false;
                }
            }

            if (flagSame)
            {
                //GPMESH_TRACE("vid {} has same normal in gpu and cpu", i);
            }
        }

        bool use_cpu_normal = true;
        bool use_gpu_normal = true;
    }

    void debug_MCF()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = STRINGIFY(INPUT_DIR) "QuadMesh/dragon-quad.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);

        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        if (!is_closed())
        {
            GPMESH_ERROR("mean curvature flow only take watertight mesh without boundaries");
        }

        auto input_coord = add_vertex_attribute<float>(vertices, "coordinates");

        //CG Solver
        auto S = add_vertex_attribute<float>("S", 3, gpmesh::DEVICE, gpmesh::SoA);
        S->reset(0.0, gpmesh::DEVICE);
        auto D = add_vertex_attribute<float>("D", 3, gpmesh::DEVICE, gpmesh::SoA);
        D->reset(0.0, gpmesh::DEVICE);
        auto R = add_vertex_attribute<float>("R", 3, gpmesh::DEVICE, gpmesh::SoA);
        R->reset(0.0, gpmesh::DEVICE);
        auto B = add_vertex_attribute<float>("B", 3, gpmesh::DEVICE, gpmesh::SoA);
        B->reset(0.0, gpmesh::DEVICE);
        auto X = add_vertex_attribute<float>("X", 3, gpmesh::DEVICE, gpmesh::SoA);
        X->copy_from(*input_coord, gpmesh::DEVICE, gpmesh::DEVICE);

        auto surround_size = add_vertex_attribute<float>("surround_size", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_size->reset(0, gpmesh::DEVICE);

        auto surround_info = add_vertex_attribute<uint32_t>("surround_info", 16, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_info->reset(0, gpmesh::DEVICE);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box_init_B;
        LaunchBox<blockThreads> launch_box_mat_vec_Multiply;
//        prepare_launch_box({gpmesh::Op::VV},
//                           launch_box_init_B,
//                           (void*)init_B<float, blockThreads>);
//        prepare_launch_box({gpmesh::Op::VV},
//                            launch_box_mat_vec_Multiply,
//                            (void*)mat_vec<float, blockThreads>);

        LaunchBox<blockThreads> launch_box_debug_kernel;
        LaunchBox<blockThreads> launch_box_debug_surround_kernel;
        LaunchBox<blockThreads> launch_box_debug_surround_kernel_v2;

        prepare_launch_box({gpmesh::Op::V_SURROUND_V},
                           launch_box_debug_kernel,
                           (void*)debug_kernel<float, blockThreads>,
                               true);

        prepare_launch_box({gpmesh::Op::V_SURROUND_V},
                           launch_box_debug_surround_kernel,
                           (void*)debug_surround_kernel<float, blockThreads>,
                               true);

        prepare_launch_box({gpmesh::Op::V_SURROUND_V},
                           launch_box_debug_surround_kernel_v2,
                           (void*)debug_surround_kerne_v2<float, blockThreads>,
                               true);


        TestData testData;
        testData.test_name = "Debug V_S_V";
        testData.num_threads = launch_box_debug_kernel.num_threads;
        testData.num_blocks  = launch_box_debug_kernel.blocks;
        testData.dyn_smem    = launch_box_debug_kernel.smem_bytes_dyn;
        testData.static_smem = launch_box_debug_kernel.smem_bytes_static;
        testData.num_reg     = launch_box_debug_kernel.num_registers_per_thread;

        TestData data_debug_surround_size;
        data_debug_surround_size.test_name = "Debug V_S_V";
        data_debug_surround_size.num_threads = launch_box_debug_surround_kernel.num_threads;
        data_debug_surround_size.num_blocks  = launch_box_debug_surround_kernel.blocks;
        data_debug_surround_size.dyn_smem    = launch_box_debug_surround_kernel.smem_bytes_dyn;
        data_debug_surround_size.static_smem = launch_box_debug_surround_kernel.smem_bytes_static;
        data_debug_surround_size.num_reg     = launch_box_debug_surround_kernel.num_registers_per_thread;

        TestData data_debug_surround_info;

        data_debug_surround_info.test_name = "Debug surround info";
        data_debug_surround_info.num_threads = launch_box_debug_surround_kernel_v2.num_threads;
        data_debug_surround_info.num_blocks  = launch_box_debug_surround_kernel_v2.blocks;
        data_debug_surround_info.dyn_smem    = launch_box_debug_surround_kernel_v2.smem_bytes_dyn;
        data_debug_surround_info.static_smem = launch_box_debug_surround_kernel_v2.smem_bytes_static;
        data_debug_surround_info.num_reg     = launch_box_debug_surround_kernel_v2.num_registers_per_thread;

        float debug_time = 0.0f;

        int num_runs = 1;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            GPUTimer timer;

            timer.start();

//            debug_kernel<float, blockThreads><<<launch_box_debug_kernel.blocks,
//                                                launch_box_debug_kernel.num_threads,
//                                                launch_box_debug_kernel.smem_bytes_dyn>>>
//                (get_context());

//            debug_surround_kernel<float, blockThreads><<<launch_box_debug_surround_kernel.blocks,
//                                                         launch_box_debug_surround_kernel.num_threads,
//                                                         launch_box_debug_surround_kernel.smem_bytes_dyn>>>
//                (get_context(), *surround_size);

            debug_surround_kerne_v2<float, blockThreads><<<launch_box_debug_surround_kernel_v2.blocks,
                                                           launch_box_debug_surround_kernel_v2.num_threads,
                                                           launch_box_debug_surround_kernel_v2.smem_bytes_dyn>>>
                (get_context(), *surround_size, *surround_info);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            testData.time_ms.push_back(timer.elapsed_millis());

            debug_time += timer.elapsed_millis();

        }

        GPMESH_TRACE("debug kernel {} (ms)",
                     debug_time / num_runs);

        surround_size->move(gpmesh::DEVICE, gpmesh::HOST);
        surround_info->move(gpmesh::DEVICE, gpmesh::HOST);

        //test for real data;
        MeshKernel::SurfaceMesh mesh_kernel;
        convert_to_mesh_kernel(vertices, fv, mesh_kernel);

        std::vector<uint32_t> neighbor_size(vertices.size(), 0);

        uint32_t max_neighbor_size = 0;

        for (auto& vp : mesh_kernel.allvertices())
        {
            auto vh = vp.first;
            auto v_id =vh.idx();

            neighbor_size[v_id] = mesh_kernel.NeighborVh(vh).size();
            max_neighbor_size = std::max(max_neighbor_size, neighbor_size[v_id]);
        }

        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh) {
           uint32_t v_id = map_to_global(vh);
           assert(neighbor_size[v_id] * 2 == (*surround_size)(vh));
        });

        std::vector<std::vector<uint32_t>> neighbor_info_from_gpu(vertices.size(), std::vector<uint32_t>());

        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh)
                        {
                            uint32_t v_id = map_to_global(vh);

                            auto [patch_id, local_id] = vh.unpack();

                            uint32_t size = (*surround_size)(vh);

                            for (int i = 0; i < size; ++i)
                            {
                                uint32_t info = (*surround_info)(vh, i);

                                uint32_t surr_patch_id = info / 100000;
                                uint32_t surr_local_id = info - 100000 * surr_patch_id;

                                uint32_t global_surround_info = m_h_patches_ltog_v[surr_patch_id][surr_local_id];

                                neighbor_info_from_gpu[v_id].push_back(global_surround_info);
                            }
                        });


        //test the neighbor_info_from_gpu with MeshKernel
        std::vector<uint32_t> error_v_ids;
        for (uint32_t v_id = 0; v_id < vertices.size(); ++v_id)
        {
            VH center(v_id);

            auto neighborVhs = mesh_kernel.NeighborVh(center);
            auto& neighborEhs = mesh_kernel.NeighborEh(center);

            int n = neighbor_info_from_gpu[v_id].size();

            for (uint32_t j = 0; j < n; ++j)
            {
                 VH surr_vh(neighbor_info_from_gpu[v_id][j]);
                 VH surr_next_vh(neighbor_info_from_gpu[v_id][(j+1)%n]);

                 MeshKernel::iGameEdge surr_e(surr_vh, surr_next_vh);
                 EH surr_eh = mesh_kernel.edgehandle(surr_e);

                 bool find_surr_eh = neighborEhs.find(surr_eh) != neighborEhs.end();
                 assert(find_surr_eh == false);

                 bool find_surr_vh = neighborVhs.find(surr_vh) != neighborVhs.end();
                 bool find_surr_next_vh = neighborVhs.find(surr_next_vh) != neighborVhs.end();
                 uint32_t test_find = static_cast<uint32_t>(find_surr_vh) + static_cast<uint32_t>(find_surr_next_vh);
                 if(test_find != 1)
                 {
                     error_v_ids.push_back(v_id);
                 }
            }
        }

        assert(error_v_ids.size() == 0);

        GPMESH_INFO("v_surround_v test Success!");
    }

    void debug_shallowWave()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = "/run/media/whalesea/Projects/ManjaroProjects/RXMesh/input/QuadMesh/plane-quad3.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        auto coords = add_vertex_attribute<float>(vertices, "coordinates");

        auto old_h = add_vertex_attribute<float>("height0", 1, gpmesh::LOCATION_ALL);
        auto h = add_vertex_attribute<float>("height1", 1, gpmesh::LOCATION_ALL);
        auto new_h = add_vertex_attribute<float>("height2", 1, gpmesh::LOCATION_ALL);


        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box;

        // compute<float, blockThreads><<<>>>
        // calculate shared memory required for each block
        prepare_launch_box({gpmesh::Op::VV},
                           launch_box,
                           (void*) shallow_wave<float, blockThreads>
        );

        TestData testData;
        testData.test_name = "Shallow Wave V1";
        testData.num_threads = launch_box.num_threads;
        testData.num_blocks  = launch_box.blocks;
        testData.dyn_smem    = launch_box.smem_bytes_dyn;
        testData.static_smem = launch_box.smem_bytes_static;
        testData.num_reg     = launch_box.num_registers_per_thread;

        float shallow_wave_time = 0;



        old_h->reset(0, gpmesh::DEVICE);
        h->reset(0, gpmesh::LOCATION_ALL);
        new_h->reset(0, gpmesh::DEVICE);

        float v = 0.2f * rand_FloatRange(0.5f, 1.0f) * 1.0f;

        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            if (v_id == 201 * 100 + 100)
            {
               h->operator()(vh) += v;
            }
            if (v_id == 201 * 99 + 100 || v_id == 201 * 101 + 100 || v_id == 201 * 100 + 100 + 1 || v_id == 201 * 100 + 100 - 1)
            {
                h->operator()(vh) -= v / 4.0f;
            }
        });

        for_each_vertex(HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);

            vertices[v_id][2] = (*h)(vh);
        });

        std::vector<std::thread> thread_pool;

        //multi-thread output
        thread_pool.push_back(std::thread(write_mesh, vertices, fv, "ShallowWaveObj/shallow_wave_frame_start.obj"));

        h->move(gpmesh::HOST, gpmesh::DEVICE);

        int num_runs = 1000;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            GPUTimer timer;

            timer.start();

            shallow_wave<float, blockThreads><<<launch_box.blocks,
                                                            launch_box.num_threads,
                                                            launch_box.smem_bytes_dyn>>>
                (get_context(), *old_h, *h, *new_h);

            timer.stop();

            std::swap(new_h, h);
            std::swap(old_h, new_h);

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            testData.time_ms.push_back(timer.elapsed_millis());
            shallow_wave_time += timer.elapsed_millis();

            //output mesh per frame
            h->move(gpmesh::DEVICE, gpmesh::HOST);

            for_each_vertex(HOST, [&](const VertexHandle& vh) {
                uint32_t v_id = map_to_global(vh);

                vertices[v_id][2] = (*h)(vh);
            });

            //multi-thread output

            if (iter % 30 == 0)
            {
                thread_pool.push_back(std::thread(write_mesh, vertices, fv, "ShallowWaveObjV2/shallow_wave_frame_" + std::to_string(iter) + ".obj"));
            }
        }

        for (auto& t : thread_pool)
        {
            t.join();
        }

        old_h->move(gpmesh::DEVICE, gpmesh::HOST);
        h->move(gpmesh::DEVICE, gpmesh::HOST);

        GPMESH_TRACE("shallow wave kernel took {} (ms)",
                     shallow_wave_time / num_runs);

    }

    void debug_cloth_simulation()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = STRINGIFY(INPUT_DIR) "QuadMesh/plane-quad.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        //test for real data;
        MeshKernel::SurfaceMesh mesh_kernel;

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        if (is_closed())
        {
            GPMESH_ERROR("cloth simulation demo only take mesh with boundaries");
        }

        auto input_coord = add_vertex_attribute<float>(vertices, "coordinates");

        auto surround_size = add_vertex_attribute<float>("surround_size", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_size->reset(0, gpmesh::DEVICE);

        auto surround_info = add_vertex_attribute<uint32_t>("surround_info", 16, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_info->reset(0, gpmesh::DEVICE);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box_debug_cloth;

        prepare_launch_box({Op::V_V_SURROUND_UNCLOSED},
                           launch_box_debug_cloth,
                           (void*)cloth_debug_kernel<float, blockThreads>,
                               true
                           );


        TestData data_cloth_debug_kernel;
        data_cloth_debug_kernel.test_name = "Debug Cloth Kernel";
        data_cloth_debug_kernel.num_threads = launch_box_debug_cloth.num_threads;
        data_cloth_debug_kernel.num_blocks  = launch_box_debug_cloth.blocks;
        data_cloth_debug_kernel.dyn_smem    = launch_box_debug_cloth.smem_bytes_dyn;
        data_cloth_debug_kernel.static_smem = launch_box_debug_cloth.smem_bytes_static;
        data_cloth_debug_kernel.num_reg     = launch_box_debug_cloth.num_registers_per_thread;

        float debug_time = 0.0f;

        int num_runs = 1;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            GPUTimer timer;

            timer.start();


            cloth_debug_kernel<float, blockThreads><<<launch_box_debug_cloth.blocks,
                                                      launch_box_debug_cloth.num_threads,
                                                      launch_box_debug_cloth.smem_bytes_dyn>>>
                //(get_context());
                (get_context(), *surround_size, *surround_info);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            //CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            data_cloth_debug_kernel.time_ms.push_back(timer.elapsed_millis());

            debug_time += timer.elapsed_millis();

        }

        GPMESH_TRACE("debug kernel {} (ms)",
                     debug_time / num_runs);


        surround_size->move(gpmesh::DEVICE, gpmesh::HOST);
        surround_info->move(gpmesh::DEVICE, gpmesh::HOST);

        convert_to_mesh_kernel(vertices, fv, mesh_kernel);

        std::vector<uint32_t> neighbor_size(vertices.size(), 0);

        uint32_t max_neighbor_size = 0;

        for (auto& vp : mesh_kernel.allvertices())
        {
            auto vh = vp.first;
            auto v_id =vh.idx();

            neighbor_size[v_id] = mesh_kernel.NeighborVh(vh).size();
            max_neighbor_size = std::max(max_neighbor_size, neighbor_size[v_id]);
        }

        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh) {
            uint32_t v_id = map_to_global(vh);
            VH mk_vh(v_id);
            uint32_t surr_size = (*surround_size)(vh);
            if (mesh_kernel.isOnBoundary(mk_vh))
            {
                assert(neighbor_size[v_id] * 2 - 1 == surr_size);
            }
            else
            {
                assert(neighbor_size[v_id] * 2 == surr_size);
            }
        });

        std::vector<std::vector<uint32_t>> neighbor_info_from_gpu(vertices.size(), std::vector<uint32_t>());

        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh)
                        {
                            uint32_t v_id = map_to_global(vh);

                            auto [patch_id, local_id] = vh.unpack();

                            uint32_t size = (*surround_size)(vh);

                            for (int i = 0; i < size; ++i)
                            {
                                uint32_t info = (*surround_info)(vh, i);

                                uint32_t surr_patch_id = info / 100000;
                                uint32_t surr_local_id = info - 100000 * surr_patch_id;

                                uint32_t global_surround_info = m_h_patches_ltog_v[surr_patch_id][surr_local_id];

                                neighbor_info_from_gpu[v_id].push_back(global_surround_info);
                            }
                        });

        //test the neighbor_info_from_gpu with MeshKernel
        std::vector<uint32_t> error_v_ids;
        for (uint32_t v_id = 0; v_id < vertices.size(); ++v_id)
        {
            VH center(v_id);

            auto neighborVhs = mesh_kernel.NeighborVh(center);
            auto& neighborEhs = mesh_kernel.NeighborEh(center);

            int n = neighbor_info_from_gpu[v_id].size();

            uint32_t iter_end = n;

            if (mesh_kernel.isOnBoundary(center))
            {
                iter_end = n - 1;
            }

            for (uint32_t j = 0; j < iter_end; ++j)
            {
                VH surr_vh(neighbor_info_from_gpu[v_id][j]);
                VH surr_next_vh(neighbor_info_from_gpu[v_id][(j+1)%n]);

                MeshKernel::iGameEdge surr_e(surr_vh, surr_next_vh);
                EH surr_eh = mesh_kernel.edgehandle(surr_e);
                assert(surr_eh.idx() != -1);

                bool find_surr_eh = neighborEhs.find(surr_eh) != neighborEhs.end();
                assert(find_surr_eh == false);

                bool find_surr_vh = neighborVhs.find(surr_vh) != neighborVhs.end();
                bool find_surr_next_vh = neighborVhs.find(surr_next_vh) != neighborVhs.end();
                uint32_t test_find = static_cast<uint32_t>(find_surr_vh) + static_cast<uint32_t>(find_surr_next_vh);
                if(test_find != 1)
                {
                    error_v_ids.push_back(v_id);
                }
            }
        }

        assert(error_v_ids.size() == 0);

        GPMESH_INFO("v_v_surround_unclosed test Success!");
    }

    void debug_cloth_simulation_explicit_integration()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = STRINGIFY(INPUT_DIR) "QuadMesh/Plane-5-21.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        //test for real data;
        MeshKernel::SurfaceMesh mesh_kernel;

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        if (is_closed())
        {
            GPMESH_ERROR("cloth simulation demo only take mesh with boundaries");
        }

        auto input_coord = add_vertex_attribute<float>(vertices, "coordinates");

        auto surround_size = add_vertex_attribute<float>("surround_size", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_size->reset(0, gpmesh::DEVICE);

        auto surround_info = add_vertex_attribute<uint32_t>("surround_info", 16, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_info->reset(0, gpmesh::DEVICE);

        auto x = add_vertex_attribute<float>(vertices, "v position");
        auto ox = add_vertex_attribute<float>(vertices, "original v positions");
        auto v = add_vertex_attribute<float>("velocity", 3, gpmesh::LOCATION_ALL, gpmesh::SoA);
        auto f = add_vertex_attribute<float>("force", 3, gpmesh::LOCATION_ALL, gpmesh::SoA);
        auto flag = add_vertex_attribute<bool>("fix flag", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);

        v->reset(0.0f, gpmesh::LOCATION_ALL);
        f->reset(0.0f, gpmesh::LOCATION_ALL);
        flag->reset(false, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box_debug_cloth;
        LaunchBox<blockThreads> launch_box_cloth_explicit_integration;

        prepare_launch_box({Op::V_V_SURROUND_UNCLOSED},
                           launch_box_debug_cloth,
                           (void*)cloth_debug_kernel<float, blockThreads>,
                           true
        );
        prepare_launch_box( {Op::V_V_SURROUND_UNCLOSED},
                           launch_box_cloth_explicit_integration,
                           (void*)cloth_explicit_integration_kernel<float, blockThreads>,
                           true
                           );

        TestData data_cloth_debug_kernel;
        data_cloth_debug_kernel.test_name = "Debug Cloth Kernel";
        data_cloth_debug_kernel.num_threads = launch_box_debug_cloth.num_threads;
        data_cloth_debug_kernel.num_blocks  = launch_box_debug_cloth.blocks;
        data_cloth_debug_kernel.dyn_smem    = launch_box_debug_cloth.smem_bytes_dyn;
        data_cloth_debug_kernel.static_smem = launch_box_debug_cloth.smem_bytes_static;
        data_cloth_debug_kernel.num_reg     = launch_box_debug_cloth.num_registers_per_thread;

        TestData data_cloth_explicit_integration_kernel;
        data_cloth_explicit_integration_kernel.test_name = "Debug Cloth Kernel";
        data_cloth_explicit_integration_kernel.num_threads =
            launch_box_cloth_explicit_integration.num_threads;
        data_cloth_explicit_integration_kernel.num_blocks  =
            launch_box_cloth_explicit_integration.blocks;
        data_cloth_explicit_integration_kernel.dyn_smem    =
            launch_box_cloth_explicit_integration.smem_bytes_dyn;
        data_cloth_explicit_integration_kernel.static_smem =
            launch_box_cloth_explicit_integration.smem_bytes_static;
        data_cloth_explicit_integration_kernel.num_reg     =
            launch_box_cloth_explicit_integration.num_registers_per_thread;

        float debug_time = 0.0f;

        float mass = 1.0f;
        float stiffness = 5e5;
        float damping = 5;
        float dt = 2e-4;
        float eps = 1e-6;


        uint32_t flag_count = 0;
        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh) -> void
                        {
                           uint32_t global_v_id = map_to_global(vh);
                           if (global_v_id == 0 || global_v_id == 20 ||
                                global_v_id == 420 || global_v_id == 440) {
                               (*flag)(vh) = true;
                               flag_count++;
                           }
                        });

        assert(flag_count == 4);
        flag->move(gpmesh::HOST, gpmesh::DEVICE);

        auto advance_lambda = [&]() -> void {
            for_each_vertex(gpmesh::DEVICE,
                            [Pos=*x,
                             Vel=*v,
                             F=*f,
                             flag=*flag,
                             dt = dt] __device__ (const VertexHandle& vh) -> void
                            {

                                if(flag(vh)) return;

                                for (int i = 0; i < 3; ++i)
                                {
                                    Vel(vh, i) += dt * F(vh, i);
                                }

                                for (int i = 0; i < 3; ++i)
                                {
                                    Pos(vh, i) += dt * Vel(vh, i);
                                }
                            });
        };

        //test advance_lambda
        {
            GPUTimer timer;

            timer.start();
            advance_lambda();
            timer.stop();
            GPMESH_TRACE("advance test takes {} (ms)",
                         timer.elapsed_millis());
        }


        int num_runs = 600;

        std::vector<std::thread> thread_pool;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            GPUTimer timer;

            timer.start();

            advance_lambda();

            cloth_explicit_integration_kernel<float, blockThreads><<<
                launch_box_cloth_explicit_integration.blocks,
                launch_box_cloth_explicit_integration.num_threads,
                launch_box_cloth_explicit_integration.smem_bytes_dyn
                >>>
                (get_context(), *x, *ox, *v, *f, dt, damping, stiffness);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            data_cloth_debug_kernel.time_ms.push_back(timer.elapsed_millis());

            debug_time += timer.elapsed_millis();

            if ((iter+1) % 100 == 0)
            {
                v->move(gpmesh::DEVICE, gpmesh::HOST);
                for_each_vertex(HOST, [&](const VertexHandle& vh) {
                    uint32_t v_id = map_to_global(vh);
                    vertices[v_id][0] = (*v)(vh, 0);
                    vertices[v_id][1] = (*v)(vh, 1);
                    vertices[v_id][2] = (*v)(vh, 2);
                });
                std::string output_path = STRINGIFY(OUTPUT_DIR) "ClothSimulation/plane3-frame";
                output_path += std::to_string(iter+1) + ".obj";
                thread_pool.push_back(std::thread(write_mesh, vertices, fv, output_path));
            }
        }

        GPMESH_TRACE("Cloth Explicit Integration kernel takes {} (ms)",
                     debug_time / num_runs);

        for (auto& t : thread_pool)
        {
            t.join();
        }

    }

    void debug_cloth_simulation_explicit_integration_cmp_with_unity()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>> vertices;
        std::string file_path = STRINGIFY(INPUT_DIR) "QuadMesh/Plane-5-21.obj";
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        //test for real data;
        MeshKernel::SurfaceMesh mesh_kernel;

        m_fv = fv;
        m_vertices = vertices;

        std::vector<uint32_t>              ff_values;
        std::vector<uint32_t>              ff_offset;
        std::vector<std::vector<uint32_t>> ef;

        build_supporting_structures(fv, ef, ff_offset, ff_values);

        //run Lloy Algorithm throught Patcher
        m_patcher = std::make_unique<patcher::Patcher>(m_patch_size,
                                                       ff_offset,
                                                       ff_values,
                                                       fv,
                                                       m_edges_map,
                                                       m_num_vertices,
                                                       m_num_edges,
                                                       m_quite);


        m_num_patches = m_patcher->get_num_patches();

        m_h_patches_ltog_f.resize(m_num_patches);
        m_h_patches_ltog_e.resize(m_num_patches);
        m_h_patches_ltog_v.resize(m_num_patches);
        m_h_num_owned_f.resize(m_num_patches);
        m_h_num_owned_v.resize(m_num_patches);
        m_h_num_owned_e.resize(m_num_patches);
        m_h_patches_fe.resize(m_num_patches);
        m_h_patches_ev.resize(m_num_patches);

        // Then build Single Patch Topology Structure

        //test
        build_all_patch_parallel(fv);

        calc_statistics(fv, ef);

        //store the patch information in GPU
        build_device();

        calc_max_not_owned_elements();

        // a vehicle to copy various parameters to the device kernels.
        // number of edge, face, vertex, patch
        // pathch information
        m_gpmesh_context.init(m_num_edges,
                              m_num_faces,
                              m_num_vertices,
                              m_num_patches,
                              m_d_patches_info);

        if (!m_quite) {
            GPMESH_TRACE("#Vertices = {}, #Faces= {}, #Edges= {}",
                         m_num_vertices,
                         m_num_faces,
                         m_num_edges);
            GPMESH_TRACE("Input is {} edge manifold",
                         ((m_is_input_edge_manifold) ? "" : " Not"));
            GPMESH_TRACE("Input is {} closed", ((m_is_input_closed) ? "" : " Not"));
            GPMESH_TRACE("max valence = {}", m_max_valence);
            GPMESH_TRACE("max edge incident faces = {}", m_max_edge_incident_faces);
            GPMESH_TRACE("max face adjacent faces = {}", m_max_face_adjacent_faces);
            GPMESH_TRACE("per-patch maximum face count = {}",
                         m_max_faces_per_patch);
            GPMESH_TRACE("per-patch maximum edge count = {}",
                         m_max_edges_per_patch);
            GPMESH_TRACE("per-patch maximum vertex count = {}",
                         m_max_vertices_per_patch);
            GPMESH_TRACE("per-patch maximum not-owned face count = {}",
                         m_max_not_owned_faces);
            GPMESH_TRACE("per-patch maximum not-owned edge count = {}",
                         m_max_not_owned_edges);
            GPMESH_TRACE("per-patch maximum not-owned vertex count = {}",
                         m_max_not_owned_vertices);
        }

        m_attr_container = std::make_shared<AttributeContainer>();

        if (is_closed())
        {
            GPMESH_ERROR("cloth simulation demo only take mesh with boundaries");
        }

        auto input_coord = add_vertex_attribute<float>(vertices, "coordinates");

        auto surround_size = add_vertex_attribute<float>("surround_size", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_size->reset(0, gpmesh::DEVICE);

        auto surround_info = add_vertex_attribute<uint32_t>("surround_info", 16, gpmesh::LOCATION_ALL, gpmesh::SoA);
        surround_info->reset(0, gpmesh::DEVICE);

        auto x = add_vertex_attribute<float>(vertices, "v position");
        auto ox = add_vertex_attribute<float>(vertices, "original v positions");
        auto v = add_vertex_attribute<float>("velocity", 3, gpmesh::LOCATION_ALL, gpmesh::SoA);
        auto f = add_vertex_attribute<float>("force", 3, gpmesh::LOCATION_ALL, gpmesh::SoA);
        auto flag = add_vertex_attribute<bool>("fix flag", 1, gpmesh::LOCATION_ALL, gpmesh::SoA);

        v->reset(0.0f, gpmesh::LOCATION_ALL);
        f->reset(0.0f, gpmesh::LOCATION_ALL);
        flag->reset(false, gpmesh::LOCATION_ALL);

        constexpr uint32_t blockThreads = 256;

        LaunchBox<blockThreads> launch_box_debug_cloth;
        LaunchBox<blockThreads> launch_box_cloth_explicit_integration;

        prepare_launch_box({Op::V_V_SURROUND_UNCLOSED},
                           launch_box_debug_cloth,
                           (void*)cloth_debug_kernel<float, blockThreads>,
                           true
        );
        prepare_launch_box( {Op::V_V_SURROUND_UNCLOSED},
                           launch_box_cloth_explicit_integration,
                           (void*)cloth_explicit_integration_kernel<float, blockThreads>,
                           true
        );

        TestData data_cloth_debug_kernel;
        data_cloth_debug_kernel.test_name = "Debug Cloth Kernel";
        data_cloth_debug_kernel.num_threads = launch_box_debug_cloth.num_threads;
        data_cloth_debug_kernel.num_blocks  = launch_box_debug_cloth.blocks;
        data_cloth_debug_kernel.dyn_smem    = launch_box_debug_cloth.smem_bytes_dyn;
        data_cloth_debug_kernel.static_smem = launch_box_debug_cloth.smem_bytes_static;
        data_cloth_debug_kernel.num_reg     = launch_box_debug_cloth.num_registers_per_thread;

        TestData data_cloth_explicit_integration_kernel;
        data_cloth_explicit_integration_kernel.test_name = "Debug Cloth Kernel";
        data_cloth_explicit_integration_kernel.num_threads =
            launch_box_cloth_explicit_integration.num_threads;
        data_cloth_explicit_integration_kernel.num_blocks  =
            launch_box_cloth_explicit_integration.blocks;
        data_cloth_explicit_integration_kernel.dyn_smem    =
            launch_box_cloth_explicit_integration.smem_bytes_dyn;
        data_cloth_explicit_integration_kernel.static_smem =
            launch_box_cloth_explicit_integration.smem_bytes_static;
        data_cloth_explicit_integration_kernel.num_reg     =
            launch_box_cloth_explicit_integration.num_registers_per_thread;

        float debug_time = 0.0f;

        float mass = 1.0f;
        float spring_k = 8000.0f;
        float damping = 0.99;
        float rho = 0.995f;
        float dt = 0.0333f;
        float eps = 0.2f;


        uint32_t flag_count = 0;
        for_each_vertex(gpmesh::HOST, [&](const VertexHandle& vh) -> void
                        {
                            uint32_t global_v_id = map_to_global(vh);
                            if (global_v_id == 0 || global_v_id == 20 ||
                                global_v_id == 420 || global_v_id == 440) {
                                (*flag)(vh) = true;
                                flag_count++;
                            }
                        });

        assert(flag_count == 4);
        flag->move(gpmesh::HOST, gpmesh::DEVICE);

        auto advance_lambda = [&]() -> void {
            for_each_vertex(gpmesh::DEVICE,
                            [Pos=*x,
                             Vel=*v,
                             F=*f,
                             flag=*flag,
                             dt = dt] __device__ (const VertexHandle& vh) -> void
                            {

                                if(flag(vh)) return;

                                for (int i = 0; i < 3; ++i)
                                {
                                    Vel(vh, i) += dt * F(vh, i);
                                }

                                for (int i = 0; i < 3; ++i)
                                {
                                    Pos(vh, i) += dt * Vel(vh, i);
                                }
                            });
        };

        //test advance_lambda
        {
            GPUTimer timer;

            timer.start();
            advance_lambda();
            timer.stop();
            GPMESH_TRACE("advance test takes {} (ms)",
                         timer.elapsed_millis());
        }


        int num_runs = 2;

        std::vector<std::thread> thread_pool;

        for (uint32_t iter = 0; iter < num_runs; ++iter)
        {
            GPUTimer timer;

            timer.start();

            advance_lambda();

            cloth_explicit_integration_kernel<float, blockThreads><<<
                launch_box_cloth_explicit_integration.blocks,
                launch_box_cloth_explicit_integration.num_threads,
                launch_box_cloth_explicit_integration.smem_bytes_dyn
                >>>
                (get_context(), *x, *ox, *v, *f, dt, damping, spring_k);

            timer.stop();

            CUDA_ERROR(cudaDeviceSynchronize());
            CUDA_ERROR(cudaGetLastError());
            CUDA_ERROR(cudaProfilerStop());
            data_cloth_debug_kernel.time_ms.push_back(timer.elapsed_millis());

            debug_time += timer.elapsed_millis();

            //if ((iter+1) % 100 == 0)
            {
                v->move(gpmesh::DEVICE, gpmesh::HOST);
                for_each_vertex(HOST, [&](const VertexHandle& vh) {
                    uint32_t v_id = map_to_global(vh);
                    vertices[v_id][0] = (*v)(vh, 0);
                    vertices[v_id][1] = (*v)(vh, 1);
                    vertices[v_id][2] = (*v)(vh, 2);
                });
                std::string output_path = STRINGIFY(OUTPUT_DIR) "ClothSimulation/Plane-5-21-frame-";
                output_path += std::to_string(iter) + ".obj";
                thread_pool.push_back(std::thread(write_mesh, vertices, fv, output_path));
            }
        }

        GPMESH_TRACE("Cloth Explicit Integration kernel takes {} (ms)",
                     debug_time / num_runs);

        for (auto& t : thread_pool)
        {
            t.join();
        }

    }
    /**
     * @brief Constructor using path to obj file
     * @param file_path path to an obj file
     * @param quite run in quite mode
     */
    GPMeshStatic(const std::string file_path, const bool quite = false)
        : GPMesh()
    {
        std::vector<std::vector<uint32_t>> fv;
        std::vector<std::vector<float>>    vertices;
        if (!import_obj(file_path, vertices, fv)) {
            GPMESH_ERROR(
                "GPMeshStatic::GPMeshStatic could not read the input file {}",
                file_path);
            exit(EXIT_FAILURE);
        }

        this->init(fv, quite);

        m_attr_container = std::make_shared<AttributeContainer>();

        m_input_vertex_coordinates =
            this->add_vertex_attribute<float>(vertices, "RX:vertices");

#if USE_POLYSCOPE
        m_polyscope_mesh = polyscope::registerSurfaceMesh(
            polyscope::guessNiceNameFromPath(file_path),
            *m_input_vertex_coordinates,
            fv);
#endif
    };

    /**
     * @brief Constructor using triangles and vertices
     * @param fv Face incident vertices as read from an obj file
     * @param quite run in quite mode
     */
    GPMeshStatic(std::vector<std::vector<uint32_t>>& fv,
                 const bool                          quite = false)
        : GPMesh(), m_input_vertex_coordinates(nullptr)
    {
        this->init(fv, quite);
        m_attr_container = std::make_shared<AttributeContainer>();
    };

    virtual ~GPMeshStatic()
    {
    }

#if USE_POLYSCOPE
    /**
     * @brief return a pointer to polyscope surface which has been registered
     * with this instance
     */
    polyscope::SurfaceMesh* get_polyscope_mesh()
    {
        return m_polyscope_mesh;
    }

    /**
     * @brief add the face's patch scalar quantity to the polyscope instance
     * associated GPMeshStatic
     * @return pointer to polyscope's face scalar quantity
     */
    polyscope::SurfaceFaceScalarQuantity* polyscope_render_face_patch()
    {
        return m_polyscope_mesh->addFaceScalarQuantity(
            "rx:FPatch", this->m_patcher->get_face_patch());
    }

    /**
     * @brief add the vertex's patch scalar quantity to the polyscope instance
     * associated GPMeshStatic
     * @return pointer to polyscope's vertex scalar quantity
     */
    polyscope::SurfaceVertexScalarQuantity* polyscope_render_vertex_patch()
    {
        return m_polyscope_mesh->addVertexScalarQuantity(
            "rx:VPatch", this->m_patcher->get_vertex_patch());
    }

    /**
     * @brief add the edge's patch scalar quantity to the polyscope instance
     * associated GPMeshStatic
     * @return pointer to polyscope's edge scalar quantity
     */
    polyscope::SurfaceEdgeScalarQuantity* polyscope_render_edge_patch()
    {
        return m_polyscope_mesh->addEdgeScalarQuantity(
            "rx:EPatch", this->m_patcher->get_edge_patch());
    }
#endif

    /**
     * @brief Apply a lambda function on all vertices in the mesh
     * @tparam LambdaT type of the lambda function (inferred)
     * @param location the execution location
     * @param apply lambda function to be applied on all vertices. The lambda
     * function signature takes a VertexHandle
     * @param stream the stream used to run the kernel in case of DEVICE
     * execution location
     */
    template <typename LambdaT>
    void for_each_vertex(locationT    location,
                         LambdaT      apply,
                         cudaStream_t stream = NULL)
    {
        if ((location & HOST) == HOST) {
            const int num_patches = this->get_num_patches();
#pragma omp parallel for
            for (int p = 0; p < num_patches; ++p) {
                for (uint16_t v = 0;
                     v < this->m_h_patches_info[p].num_owned_vertices;
                     ++v) {
                    const VertexHandle v_handle(static_cast<uint32_t>(p), v);
                    apply(v_handle);
                }
            }
        }

        if ((location & DEVICE) == DEVICE) {

            if constexpr (IS_HD_LAMBDA(LambdaT) || IS_D_LAMBDA(LambdaT)) {

                const int num_patches = this->get_num_patches();
                const int threads     = 256;
                detail::for_each_vertex<<<num_patches, threads, 0, stream>>>(
                    num_patches, this->m_d_patches_info, apply);
            } else {
                GPMESH_ERROR(
                    "GPMeshStatic::for_each_vertex() Input lambda function "
                    "should be annotated with  __device__ for execution on "
                    "device");
            }
        }
    }

    /**
     * @brief Apply a lambda function on all edges in the mesh
     * @tparam LambdaT type of the lambda function (inferred)
     * @param location the execution location
     * @param apply lambda function to be applied on all edges. The lambda
     * function signature takes a EdgeHandle
     * @param stream the stream used to run the kernel in case of DEVICE
     * execution location
     */
    template <typename LambdaT>
    void for_each_edge(locationT    location,
                       LambdaT      apply,
                       cudaStream_t stream = NULL)
    {
        if ((location & HOST) == HOST) {
            const int num_patches = this->get_num_patches();
#pragma omp parallel for
            for (int p = 0; p < num_patches; ++p) {
                for (uint16_t e = 0;
                     e < this->m_h_patches_info[p].num_owned_edges;
                     ++e) {
                    const EdgeHandle e_handle(static_cast<uint32_t>(p), e);
                    apply(e_handle);
                }
            }
        }

        if ((location & DEVICE) == DEVICE) {
            if constexpr (IS_HD_LAMBDA(LambdaT) || IS_D_LAMBDA(LambdaT)) {

                const int num_patches = this->get_num_patches();
                const int threads     = 256;
                detail::for_each_edge<<<num_patches, threads, 0, stream>>>(
                    num_patches, this->m_d_patches_info, apply);
            } else {
                GPMESH_ERROR(
                    "GPMeshStatic::for_each_edge() Input lambda function "
                    "should be annotated with  __device__ for execution on "
                    "device");
            }
        }
    }

    /**
     * @brief Apply a lambda function on all faces in the mesh
     * @tparam LambdaT type of the lambda function (inferred)
     * @param location the execution location
     * @param apply lambda function to be applied on all faces. The lambda
     * function signature takes a FaceHandle
     * @param stream the stream used to run the kernel in case of DEVICE
     * execution location
     */
    template <typename LambdaT>
    void for_each_face(locationT    location,
                       LambdaT      apply,
                       cudaStream_t stream = NULL)
    {
        if ((location & HOST) == HOST) {
            const int num_patches = this->get_num_patches();
#pragma omp parallel for
            for (int p = 0; p < num_patches; ++p) {
                for (int f = 0; f < this->m_h_patches_info[p].num_owned_faces;
                     ++f) {
                    const FaceHandle f_handle(static_cast<uint32_t>(p), f);
                    apply(f_handle);
                }
            }
        }

        if ((location & DEVICE) == DEVICE) {
            if constexpr (IS_HD_LAMBDA(LambdaT) || IS_D_LAMBDA(LambdaT)) {

                const int num_patches = this->get_num_patches();
                const int threads     = 256;
                detail::for_each_face<<<num_patches, threads, 0, stream>>>(
                    num_patches, this->m_d_patches_info, apply);
            } else {
                GPMESH_ERROR(
                    "GPMeshStatic::for_each_face() Input lambda function "
                    "should be annotated with  __device__ for execution on "
                    "device");
            }
        }
    }

    /**
     * @brief populate the launch_box with grid size and dynamic shared memory
     * needed for kernel launch
     * @param op List of query operations done inside this the kernel
     * @param launch_box input launch box to be populated
     * @param kernel The kernel to be launched
     * @param oriented if the query is oriented. Valid only for Op::VV queries
     */
    template <uint32_t blockThreads>
    void prepare_launch_box(const std::vector<Op>    op,
                            LaunchBox<blockThreads>& launch_box,
                            const void*              kernel,
                            const bool               oriented = false) const
    {

        launch_box.blocks         = this->m_num_patches;
        launch_box.smem_bytes_dyn = 0;

        for (auto o : op) {
            launch_box.smem_bytes_dyn = std::max(
                launch_box.smem_bytes_dyn,
                this->template calc_shared_memory<blockThreads>(o, oriented));
        }

        if (!this->m_quite) {
            GPMESH_TRACE(
                "GPMeshStatic::calc_shared_memory() launching {} blocks with "
                "{} threads on the device",
                launch_box.blocks,
                blockThreads);
        }

        check_shared_memory(launch_box.smem_bytes_dyn,
                            launch_box.smem_bytes_static,
                            launch_box.num_registers_per_thread,
                            kernel);
    }


    /**
     * @brief Adding a new face attribute
     * @tparam T type of the attribute
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param num_attributes number of the attributes
     * @param location where to allocate the attributes
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<FaceAttribute<T>> add_face_attribute(
        const std::string& name,
        uint32_t           num_attributes,
        locationT          location = LOCATION_ALL,
        layoutT            layout   = SoA)
    {
        return m_attr_container->template add<FaceAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_f,
            num_attributes,
            location,
            layout,
            this);
    }

    /**
     * @brief Adding a new face attribute by reading values from a host buffer
     * f_attributes where the order of faces is the same as the order of
     * faces given to the constructor.The attributes are populated on device
     * and host
     * @tparam T type of the attribute
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<FaceAttribute<T>> add_face_attribute(
        const std::vector<std::vector<T>>& f_attributes,
        const std::string&                 name,
        layoutT                            layout = SoA)
    {
        if (f_attributes.empty()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_face_attribute() input attribute is empty");
        }

        if (f_attributes.size() != get_num_faces()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_face_attribute() input attribute size ({}) "
                "is not the same as number of faces in the input mesh ({})",
                f_attributes.size(),
                get_num_faces());
        }

        uint32_t num_attributes = f_attributes[0].size();

        auto ret = m_attr_container->template add<FaceAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_f,
            num_attributes,
            LOCATION_ALL,
            layout,
            this);

        // populate the attribute before returning it
        const int num_patches = this->get_num_patches();
#pragma omp parallel for
        for (int p = 0; p < num_patches; ++p) {
            for (uint16_t f = 0; f < this->m_h_num_owned_f[p]; ++f) {

                const FaceHandle f_handle(static_cast<uint32_t>(p), f);

                uint32_t global_f = m_h_patches_ltog_f[p][f];

                for (uint32_t a = 0; a < num_attributes; ++a) {
                    (*ret)(f_handle, a) = f_attributes[global_f][a];
                }
            }
        }

        // move to device
        ret->move(gpmesh::HOST, gpmesh::DEVICE);
        return ret;
    }

    /**
     * @brief Adding a new face attribute by reading values from a host buffer
     * f_attributes where the order of faces is the same as the order of
     * faces given to the constructor.The attributes are populated on device
     * and host
     * @tparam T type of the attribute
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<FaceAttribute<T>> add_face_attribute(
        const std::vector<T>& f_attributes,
        const std::string&    name,
        layoutT               layout = SoA)
    {
        if (f_attributes.empty()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_face_attribute() input attribute is empty");
        }

        if (f_attributes.size() != get_num_faces()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_face_attribute() input attribute size ({}) "
                "is not the same as number of faces in the input mesh ({})",
                f_attributes.size(),
                get_num_faces());
        }

        uint32_t num_attributes = 1;

        auto ret = m_attr_container->template add<FaceAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_f,
            num_attributes,
            LOCATION_ALL,
            layout,
            this);

        // populate the attribute before returning it
        const int num_patches = this->get_num_patches();
#pragma omp parallel for
        for (int p = 0; p < num_patches; ++p) {
            for (uint16_t f = 0; f < this->m_h_num_owned_f[p]; ++f) {

                const FaceHandle f_handle(static_cast<uint32_t>(p), f);

                uint32_t global_f = m_h_patches_ltog_f[p][f];

                (*ret)(f_handle, 0) = f_attributes[global_f];
            }
        }

        // move to device
        ret->move(gpmesh::HOST, gpmesh::DEVICE);
        return ret;
    }

    /**
     * @brief Adding a new edge attribute
     * @tparam T type of the attribute
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param num_attributes number of the attributes
     * @param location where to allocate the attributes
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<EdgeAttribute<T>> add_edge_attribute(
        const std::string& name,
        uint32_t           num_attributes,
        locationT          location = LOCATION_ALL,
        layoutT            layout   = SoA)
    {
        return m_attr_container->template add<EdgeAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_e,
            num_attributes,
            location,
            layout,
            this);
    }

    /**
     * @brief Adding a new vertex attribute
     * @tparam T type of the attribute
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param num_attributes number of the attributes
     * @param location where to allocate the attributes
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<VertexAttribute<T>> add_vertex_attribute(
        const std::string& name,
        uint32_t           num_attributes,
        locationT          location = LOCATION_ALL,
        layoutT            layout   = SoA)
    {
        return m_attr_container->template add<VertexAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_v,
            num_attributes,
            location,
            layout,
            this);
    }

    /**
     * @brief Adding a new vertex attribute by reading values from a host buffer
     * v_attributes where the order of vertices is the same as the order of
     * vertices given to the constructor. The attributes are populated on device
     * and host
     * @tparam T type of the attribute
     * @param v_attributes attributes to read
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<VertexAttribute<T>> add_vertex_attribute(
        const std::vector<std::vector<T>>& v_attributes,
        const std::string&                 name,
        layoutT                            layout = SoA)
    {
        if (v_attributes.empty()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_vertex_attribute() input attribute is "
                "empty");
        }

        if (v_attributes.size() != get_num_vertices()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_vertex_attribute() input attribute size "
                "({}) is not the same as number of vertices in the input mesh "
                "({})",
                v_attributes.size(),
                get_num_vertices());
        }

        uint32_t num_attributes = v_attributes[0].size();

        //allocate host and device memory
        auto ret = m_attr_container->template add<VertexAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_v,
            num_attributes,
            LOCATION_ALL,
            layout,
            this);

        // populate the attribute before returning it
        const int num_patches = this->get_num_patches();
#pragma omp parallel for
        for (int p = 0; p < num_patches; ++p) {
            for (uint16_t v = 0; v < this->m_h_num_owned_v[p]; ++v) {

                const VertexHandle v_handle(static_cast<uint32_t>(p), v);

                uint32_t global_v = m_h_patches_ltog_v[p][v];

                for (uint32_t a = 0; a < num_attributes; ++a) {
                    (*ret)(v_handle, a) = v_attributes[global_v][a];
                }
            }
        }

        // move to device
        ret->move(gpmesh::HOST, gpmesh::DEVICE);
        return ret;
    }

    /**
     * @brief Adding a new vertex attribute by reading values from a host buffer
     * v_attributes where the order of vertices is the same as the order of
     * vertices given to the constructor. The attributes are populated on device
     * and host
     * @tparam T type of the attribute
     * @param v_attributes attributes to read
     * @param name of the attribute. Should not collide with other attributes
     * names
     * @param layout as SoA or AoS
     * operations
     * @return shared pointer to the created attribute
     */
    template <class T>
    std::shared_ptr<VertexAttribute<T>> add_vertex_attribute(
        const std::vector<T>& v_attributes,
        const std::string&    name,
        layoutT               layout = SoA)
    {
        if (v_attributes.empty()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_vertex_attribute() input attribute is "
                "empty");
        }

        if (v_attributes.size() != get_num_vertices()) {
            GPMESH_ERROR(
                "GPMeshStatic::add_vertex_attribute() input attribute size "
                "({}) is not the same as number of vertices in the input mesh "
                "({})",
                v_attributes.size(),
                get_num_vertices());
        }

        uint32_t num_attributes = 1;

        auto ret = m_attr_container->template add<VertexAttribute<T>>(
            name.c_str(),
            this->m_h_num_owned_v,
            num_attributes,
            LOCATION_ALL,
            layout,
            this);

        // populate the attribute before returning it
        const int num_patches = this->get_num_patches();
#pragma omp parallel for
        for (int p = 0; p < num_patches; ++p) {
            for (uint16_t v = 0; v < this->m_h_num_owned_v[p]; ++v) {

                const VertexHandle v_handle(static_cast<uint32_t>(p), v);

                uint32_t global_v = m_h_patches_ltog_v[p][v];

                (*ret)(v_handle, 0) = v_attributes[global_v];
            }
        }

        // move to device
        ret->move(gpmesh::HOST, gpmesh::DEVICE);
        return ret;
    }

    /**
     * @brief Checks if an attribute exists given its name
     * @param name the attribute name
     * @return True if the attribute exists. False otherwise.
     */
    bool does_attribute_exist(const std::string& name)
    {
        return m_attr_container->does_exist(name.c_str());
    }

    /**
     * @brief Remove an attribute. Could be vertex, edge, or face attribute
     * @param name the attribute name
     */
    void remove_attribute(const std::string& name)
    {
        if (!this->does_attribute_exist(name)) {
            GPMESH_WARN(
                "GPMeshStatic::remove_attribute() trying to remove an "
                "attribute that does not exit with name {}",
                name);
            return;
        }

        m_attr_container->remove(name.c_str());
    }


    /**
     * @brief return a shared pointer the input vertex position
     */
    std::shared_ptr<VertexAttribute<float>> get_input_vertex_coordinates()
    {
        if (!m_input_vertex_coordinates) {
            GPMESH_ERROR(
                "GPMeshStatic::get_input_vertex_coordinates input vertex was "
                "not initialized. Call GPMeshStatic with constructor to the "
                "obj file path");
            exit(EXIT_FAILURE);
        }
        return m_input_vertex_coordinates;
    }

    /**
     * @brief Map a vertex handle into a global index as seen in the input
     * to GPMeshStatic
     * @param vh input vertex handle
     * @return the global index of vh
     */
    uint32_t map_to_global(const VertexHandle vh) const
    {
        auto pl = vh.unpack();
        return m_h_patches_ltog_v[pl.first][pl.second];
    }

    /**
     * @brief Map an edge handle into a global index
     * @param eh input edge handle
     * @return the global index of eh
     */
    uint32_t map_to_global(const EdgeHandle eh) const
    {
        auto pl = eh.unpack();
        return m_h_patches_ltog_e[pl.first][pl.second];
    }

    /**
     * @brief Map a face handle into a global index as seen in the input
     * to GPMeshStatic
     * @param vh input face handle
     * @return the global index of fh
     */
    uint32_t map_to_global(const FaceHandle fh) const
    {
        auto pl = fh.unpack();
        return m_h_patches_ltog_f[pl.first][pl.second];
    }

    /**
     * @brief Export the mesh to obj file
     * @tparam T type of vertices coordinates
     * @param filename the output file
     * @param coords vertices coordinates
     */
    template <typename T>
    void export_obj(const std::string&        filename,
                    const VertexAttribute<T>& coords)
    {
        std::string  fn = filename;
        std::fstream file(fn, std::ios::out);
        file.precision(30);

        uint32_t num_v = 0;
        for (uint32_t p = 0; p < this->m_num_patches; ++p) {

            const uint32_t p_num_vertices =
                this->m_h_patches_info[p].num_vertices;

            for (uint16_t v = 0; v < p_num_vertices; ++v) {
                uint16_t v_id = v;
                uint32_t p_id = p;
                if (v >= this->m_h_patches_info[p].num_owned_vertices) {
                    uint16_t l =
                        v - this->m_h_patches_info[p].num_owned_vertices;
                    v_id = this->m_h_patches_info[p].not_owned_id_v[l].id;
                    p_id = this->m_h_patches_info[p].not_owned_patch_v[l];
                }
                VertexHandle vh(p_id, {v_id});
                file << "v " << coords(vh, 0) << " " << coords(vh, 1) << " "
                     << coords(vh, 2) << std::endl;
            }

            const uint32_t p_num_faces =
                this->m_h_patches_info[p].num_owned_faces;

            for (uint32_t f = 0; f < p_num_faces; ++f) {

                file << "f ";
                for (uint32_t e = 0; e < 3; ++e) {
                    uint16_t edge = this->m_h_patches_info[p].fe[3 * f + e].id;
                    flag_t   dir(0);
                    Context::unpack_edge_dir(edge, edge, dir);
                    uint16_t e_id = (2 * edge) + dir;
                    uint16_t v    = this->m_h_patches_info[p].ev[e_id].id;
                    file << v + num_v + 1 << " ";
                }
                file << std::endl;
            }

            num_v += p_num_vertices;
        }
    }

   protected:
    template <uint32_t blockThreads>
    size_t calc_shared_memory(const Op op, const bool oriented = false) const
    {
        // Operations that uses matrix transpose needs a template parameter
        // that is by default TRANSPOSE_ITEM_PER_THREAD. Here we check if
        // this default parameter is valid otherwise, it needs to be increased.
        if (op == Op::VV || op == Op::VE || op == Op::V_SURROUND_V || op == Op::V_V_SURROUND_UNCLOSED) {
            if (2 * this->m_max_edges_per_patch >
                blockThreads * TRANSPOSE_ITEM_PER_THREAD) {
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() "
                    "TRANSPOSE_ITEM_PER_THREAD = {} needs "
                    "to be increased for op = {}",
                    TRANSPOSE_ITEM_PER_THREAD,
                    op_to_string(op));
            }
        } else if (op == Op::VE || op == Op::EF || op == Op::FF) {
            if (4 * this->m_max_faces_per_patch >
                blockThreads * TRANSPOSE_ITEM_PER_THREAD) {
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() "
                    "TRANSPOSE_ITEM_PER_THREAD = {} needs "
                    "to be increased for op = {}",
                    TRANSPOSE_ITEM_PER_THREAD,
                    op_to_string(op));
            }
        }


        if (oriented && (op != Op::VV && op!= Op::V_SURROUND_V && op != Op::V_V_SURROUND_UNCLOSED)) {
            GPMESH_ERROR(
                "GPMeshStatic::calc_shared_memory() Oriented is only "
                "allowed on VV and V_SURROUND_V. The input op is {}",
                op_to_string(op));
        }

        if (oriented && (!(op == Op::VV && is_closed())
                         && !(op == Op::VV && is_closed())
                         && !(op == Op::V_V_SURROUND_UNCLOSED && !is_closed()))
            )
        {
            GPMESH_ERROR(
                "GPMeshStatic::calc_shared_memory() Can't generate oriented "
                "output for input with boundaries "
                "or Can`t generate oriented v_v_surround_unclosed for input without boundaryes"
                );
        }

        size_t dynamic_smem = 0;

        if (op == Op::FE) {
            // only FE will be loaded
            dynamic_smem = 4 * this->m_max_faces_per_patch * sizeof(uint16_t);
            // to load not-owned edges local and patch id
            dynamic_smem += this->m_max_not_owned_edges * sizeof(uint32_t);
            dynamic_smem += this->m_max_not_owned_edges * sizeof(uint16_t);
        } else if (op == Op::EV) {
            // only EV will be loaded
            dynamic_smem = 2 * this->m_max_edges_per_patch * sizeof(uint16_t);
            // to load not-owned vertices local and patch id
            dynamic_smem += this->m_max_not_owned_vertices * sizeof(uint32_t);
            dynamic_smem += this->m_max_not_owned_vertices * sizeof(uint16_t);
        } else if (op == Op::FV) {
            // We load both FE and EV. We don't change EV.
            // FE are updated to contain FV instead of FE by reading from
            // EV
            dynamic_smem = 2 * this->m_max_edges_per_patch * sizeof(uint16_t);
            dynamic_smem += 4 * this->m_max_faces_per_patch * sizeof(uint16_t);
            // no need for extra memory to load not-owned vertices local and
            // patch id. We load them and overwrite EV.
            const uint32_t not_owned_v_bytes =
                this->m_max_not_owned_vertices *
                (sizeof(uint16_t) + sizeof(uint32_t));
            const uint32_t edges_bytes =
                2 * this->m_max_edges_per_patch * sizeof(uint16_t);
            if (not_owned_v_bytes > edges_bytes) {
                // dynamic_smem += not_owned_v_bytes - edges_bytes;
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() FV query might fail!");
            }
        } else if (op == Op::VE) {
            // load EV and then transpose it in place
            // The transpose needs two buffer; one for prefix sum and another
            // for the actual output
            // The prefix sum will be stored in place (where EV are loaded)
            // The output will be stored in another buffer with size equal to
            // the EV (i.e., 2*#edges) since this output buffer will stored the
            // nnz and the nnz of a matrix the same before/after transpose
            dynamic_smem =
                (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);

            // to load the not-owned edges local and patch id
            dynamic_smem += this->m_max_not_owned_edges * sizeof(uint32_t);
            dynamic_smem += this->m_max_not_owned_edges * sizeof(uint16_t);
        } else if (op == Op::EF) {
            // same as Op::VE but with faces
            dynamic_smem =
                (2 * 4 * this->m_max_faces_per_patch) * sizeof(uint16_t) +
                sizeof(uint16_t) + sizeof(uint16_t);

            // to load the not-owned faces local and patch id
            dynamic_smem += this->m_max_not_owned_faces * sizeof(uint32_t);
            dynamic_smem += this->m_max_not_owned_faces * sizeof(uint16_t);
        } else if (op == Op::VF) {
            // load EV and FE simultaneously. changes FE to FV using EV. Then
            // transpose FV in place and use EV to store the values/output while
            // using FV to store the prefix sum. Thus, the space used to store
            // EV should be max(3*#faces, 2*#edges)
            dynamic_smem = 3 * this->m_max_faces_per_patch * sizeof(uint16_t);
            dynamic_smem += std::max(4 * this->m_max_faces_per_patch,
                                     2 * this->m_max_edges_per_patch) *
                                sizeof(uint16_t) +
                            sizeof(uint16_t);

            // to load the not-owned faces local and patch id
            dynamic_smem += this->m_max_not_owned_faces * sizeof(uint32_t);
            dynamic_smem += this->m_max_not_owned_faces * sizeof(uint16_t);
        } else if (op == Op::VV) {
            // similar to VE but we also need to store the EV even after
            // we do the transpose
            dynamic_smem =
                (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            dynamic_smem +=
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            // no need for extra memory to load not-owned local and patch id.
            // We load them and overwrite the extra EV
            if (this->m_max_not_owned_vertices *
                    (sizeof(uint16_t) + sizeof(uint32_t)) >
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t)) {
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() VV query might fail!");
            }

        } else if (op == Op::V_SURROUND_V)
        {
            dynamic_smem =
                (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            dynamic_smem +=
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            // no need for extra memory to load not-owned local and patch id.
            // We load them and overwrite the extra EV
            if (this->m_max_not_owned_vertices *
                    (sizeof(uint16_t) + sizeof(uint32_t)) >
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t)) {
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() V_V_Surround query might fail!");
            }
        }
        else if (op == Op::V_V_SURROUND_UNCLOSED)
        {
            dynamic_smem =
                (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            dynamic_smem +=
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
            // no need for extra memory to load not-owned local and patch id.
            // We load them and overwrite the extra EV
            if (this->m_max_not_owned_vertices *
                    (sizeof(uint16_t) + sizeof(uint32_t)) >
                (2 * this->m_max_edges_per_patch) * sizeof(uint16_t)) {
                GPMESH_ERROR(
                    "GPMeshStatic::calc_shared_memory() V_V_Surround_unclosed query might fail!");
            }
        }
        else if (op == Op::FF) {
            // FF needs to store FE and EF along side with the output itself
            // FE needs 3*max_num_faces
            // EF is FE transpose
            // FF is max_num_faces + (on average) 3*max_num_faces
            // Since we have so many boundary faces (due to ribbons), they will
            // make up this averaging

            dynamic_smem = (3 * this->m_max_faces_per_patch +        // FE
                            2 * (3 * this->m_max_faces_per_patch) +  // EF
                            4 * this->m_max_faces_per_patch) *       // FF
                           sizeof(uint16_t);
            // no need for extra memory to load not-owned faces local and
            // patch id. We load them and overwrite FE.
        }

        if (op == Op::VV && oriented) {
            // For oriented VV, we additionally need to store FE and EF along
            // with the (transposed) VE
            // FE needs 4*max_num_faces
            // Since oriented is only done on manifold, EF needs only
            // 2*max_num_edges since every edge is neighbor to maximum of two
            // faces (which we write on the same place as the extra EV)
            dynamic_smem +=
                (4 * this->m_max_faces_per_patch) * sizeof(uint16_t);
        }

        if(op == Op::V_SURROUND_V && oriented)
        {
            dynamic_smem +=
                (4 * this->m_max_faces_per_patch) * sizeof(uint16_t);
            dynamic_smem += (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
        }

        if (op == Op::V_V_SURROUND_UNCLOSED && oriented)
        {
            dynamic_smem +=
                (4 * this->m_max_faces_per_patch) * sizeof(uint16_t);
            dynamic_smem += (2 * 2 * this->m_max_edges_per_patch) * sizeof(uint16_t);
        }

        //return 40000; //3070 laptop max-49152
        return dynamic_smem;
    }

    void check_shared_memory(const uint32_t smem_bytes_dyn,
                             size_t&        smem_bytes_static,
                             uint32_t&      num_reg_per_thread,
                             const void*    kernel) const
    {
        // check if total shared memory (static + dynamic) consumed by
        // k_base_query are less than the max shared per block

        //here especially set shared memory per block to 96kb
        int maxbytes =  98304;
        cudaFuncSetAttribute(kernel, cudaFuncAttributeMaxDynamicSharedMemorySize, maxbytes);

        cudaFuncAttributes func_attr = cudaFuncAttributes();
        CUDA_ERROR(cudaFuncGetAttributes(&func_attr, kernel));


        smem_bytes_static  = func_attr.sharedSizeBytes;
        num_reg_per_thread = static_cast<uint32_t>(func_attr.numRegs);
        int device_id;
        CUDA_ERROR(cudaGetDevice(&device_id));
        cudaDeviceProp devProp;
        CUDA_ERROR(cudaGetDeviceProperties(&devProp, device_id));

        if (!this->m_quite) {
            GPMESH_TRACE(
                "GPMeshStatic::check_shared_memory() user function requires "
                "shared memory = {} (dynamic) + {} (static) = {} (bytes) and "
                "{} registers per thread",
                smem_bytes_dyn,
                smem_bytes_static,
                smem_bytes_dyn + smem_bytes_static,
                num_reg_per_thread);

            GPMESH_TRACE(
                "GPMeshStatic::check_shared_memory() available total shared "
                "memory per block = {} (bytes) = {} (Kb)",
                devProp.sharedMemPerBlock,
                float(devProp.sharedMemPerBlock) / 1024.0f);
        }

        if (smem_bytes_static + smem_bytes_dyn > devProp.sharedMemPerBlock) {
            GPMESH_ERROR(
                " GPMeshStatic::check_shared_memory() shared memory needed for"
                " input function ({} bytes) exceeds the max shared memory "
                "per block on the current device ({} bytes)",
                smem_bytes_static + smem_bytes_dyn,
                devProp.sharedMemPerBlock);
            exit(EXIT_FAILURE);
        }
    }

    std::shared_ptr<AttributeContainer>     m_attr_container;
    std::shared_ptr<VertexAttribute<float>> m_input_vertex_coordinates;
};
}  // namespace gpmesh