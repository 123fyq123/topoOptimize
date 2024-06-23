#include "GP_Spline_CPU.h"
#include "GP_IGA_CPU_Kernel.h"
#include "Util/GP_Timer.h"
#include "GP_Spline_Attribute_CPU.h"
#include "Hyperelasticity/stiffAssembler.h"
#include <tbb/parallel_for.h>
#define EIGEN_USE_MKL_ALL
#include <Eigen/PardisoSupport>
#include <iostream>

namespace gpmesh {

    template<typename T>
    void GPSplineSolidCPU<T>::init(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                   std::vector<boundary_point> &neumann_points, const bool quite) {
        m_sparse_system = SparseSystem(patch, dirichlet_points, 3);
        m_num_geometry = patch.size();
        m_num_cpoints = m_sparse_system.numControlPoints();
        assert(m_num_cpoints == m_sparse_system.get_g2l_map().size());

        build_assembly_context(patch);
        GPMESH_INFO("cpu assembly context init success !");

        build_mass_context(patch);
        GPMESH_INFO("cpu mass context init success !");

        build_physics_context(patch);
        GPMESH_INFO("cpu physics context init success !");

        build_dirichlet_point_context(dirichlet_points);
        GPMESH_INFO("cpu dirichlet point context init success !");
    }
    template<typename T>
    void GPSplineSolidCPU<T>::init_patch(const vector<CBSplineSolid> &patch)
    {
        m_sparse_system = SparseSystem();
        m_sparse_system.init_geometry(patch, 3);
        m_num_cpoints = m_sparse_system.numControlPoints();

        m_num_geometry = patch.size();

        assert(m_num_cpoints == m_sparse_system.get_g2l_map().size());

        build_assembly_context(patch);
        GPMESH_INFO("cpu assembly context init success !");

        build_mass_context(patch);
        GPMESH_INFO("cpu mass context init success !");

        build_physics_context(patch);
        GPMESH_INFO("cpu physics context init success !");

        bInitPatch = true;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::init_dirichlet_points(std::vector<boundary_point> &dirichlet_points)
    {
        m_sparse_system.init_dirichlet(dirichlet_points);
        build_dirichlet_point_context(dirichlet_points);
        GPMESH_INFO("cpu dirichlet point context init success !");

        bInitDirichlet = true;
    }

    template<typename T>
    GPSplineSolidCPU<T>::~GPSplineSolidCPU() {
        m_mass_context.release();
        m_assembly_context.release();
        m_physics_context.release();
        m_dirichlet_point_context.release();
    }

    template<typename T>
    void GPSplineSolidCPU<T>::build_assembly_context(const vector<CBSplineSolid> &patch) {
        AssemblyInfoCPU<T> *m_h_assembly_info;

        assert(m_num_geometry > 0);
        m_h_assembly_info = new AssemblyInfoCPU<T>[m_num_geometry]();

        for (size_t g = 0; g < m_num_geometry; ++g) {
            m_h_assembly_info[g].local_mats = new T[64 * 3 * 64 * 3]();
            m_h_assembly_info[g].local_rhs = new T[64 * 3]();
            m_h_assembly_info[g].geometry_id = g;
        }

        m_assembly_context.init(m_h_assembly_info);
    };

    template<typename T>
    void GPSplineSolidCPU<T>::build_mass_context(const vector<CBSplineSolid> &patch) {
        MassInfoCPU<T> *m_h_mass_info;

        assert(m_num_geometry > 0);
        m_h_mass_info = new MassInfoCPU<T>[m_num_geometry]();

        for (size_t g = 0; g < m_num_geometry; ++g) {
            m_h_mass_info[g].mass_local_mats = new T[64 * 3 * 64 * 3]();
            m_h_mass_info[g].geometry_id = g;
        }

        m_mass_context.init(m_h_mass_info);
    }

    template<typename T>
    void GPSplineSolidCPU<T>::build_physics_context(const vector<CBSplineSolid> &patch) {
        PhysicsInfoCPU<T> *m_h_physics_info;
        assert(m_num_geometry > 0);
        m_h_physics_info = new PhysicsInfoCPU<T>[m_num_geometry]();

        for (size_t g = 0; g < m_num_geometry; ++g) {
            m_h_physics_info[g].m_cpoints_displacement = Eigen::Matrix<T, 64, 3>::Zero();
            m_h_physics_info[g].geometry_id = g;
        }

        m_physics_context.init(m_h_physics_info);
    }

    template<typename T>
    void GPSplineSolidCPU<T>::build_dirichlet_point_context(const vector<boundary_point> &dirichlet_points) {
        //the process of building dirichlet point context depends on SparseSystem of Patches
        assert(m_sparse_system.getControlNums() > 0);

        uint32_t num_dirichlet_point = dirichlet_points.size();

        vector<DirichletPointInfoCPU<T>> dirichletInfoVec;
        auto &g2l_map = m_sparse_system.get_g2l_map();

        unordered_set<uint32_t> gid_set;

        for (uint32_t id = 0; id < num_dirichlet_point; ++id) {
            auto &bd_point = dirichlet_points[id];

            assert(bd_point.m_condition_type == 0);

            uint32_t bp_global_id = bd_point.m_cpoint_global_id;

            uint32_t num_local = g2l_map[bp_global_id].size();

            T valueX = bd_point.m_val.getX();
            T valueY = bd_point.m_val.getY();
            T valueZ = bd_point.m_val.getZ();

            for (uint32_t local = 0; local < num_local; ++local) {
                uint32_t patch_id = g2l_map[bp_global_id][local] / 64;
                gid_set.insert(patch_id);
                uint16_t local_id = g2l_map[bp_global_id][local] - 64 * patch_id;
                assert(patch_id < m_num_geometry);
                assert(local_id >= 0 && local_id < 64);
                assert(m_sparse_system.get_l2g_map()[patch_id][local_id] == bp_global_id);

                dirichletInfoVec.push_back({});
                dirichletInfoVec.back().valueX = valueX;
                dirichletInfoVec.back().valueY = valueY;
                dirichletInfoVec.back().valueZ = valueZ;
                dirichletInfoVec.back().global_id = bp_global_id;
                dirichletInfoVec.back().geometry_id = patch_id;
                dirichletInfoVec.back().local_id = local_id;
            }
        }

        DirichletPointInfoCPU<T> *m_h_dirichlet_point_info = new DirichletPointInfoCPU<T>[dirichletInfoVec.size()]();

        for (uint32_t id = 0; id < dirichletInfoVec.size(); ++id) {
            m_h_dirichlet_point_info[id] = dirichletInfoVec[id];
        }

        m_dirichlet_point_context.init(m_h_dirichlet_point_info, dirichletInfoVec.size());
    }

    template<typename T>
    void
    GPSplineSolidCPU<T>::build_rhs(shared_ptr<SplineAttributeCPU<T>> rhs, std::vector<boundary_point> &dirichlet_points,
                                   std::vector<boundary_point> &neumann_points) {

        tbb::parallel_for(tbb::blocked_range<size_t>(0, neumann_points.size()),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t i = r.begin(); i < r.end(); ++i) {
                                  const auto &n_point = neumann_points[i];

                                  assert(n_point.m_condition_type == 1);

                                  uint32_t global_id = n_point.m_cpoint_global_id;
                                  rhs->operator()(global_id, 0) = n_point.m_val.getX();
                                  rhs->operator()(global_id, 1) = n_point.m_val.getY();
                                  rhs->operator()(global_id, 2) = n_point.m_val.getZ();
                              }
                          }
        );

        //dirichlet boundary
        tbb::parallel_for(tbb::blocked_range<size_t>(0, dirichlet_points.size()),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t i = r.begin(); i < r.end(); ++i) {
                                  const auto &d_point = dirichlet_points[i];

                                  assert(d_point.m_condition_type == 0);

                                  uint32_t global_id = d_point.m_cpoint_global_id;

                                  uint32_t bd_size = m_sparse_system.get_g2l_map()[global_id].size();

                                  rhs->operator()(global_id, 0) = d_point.m_val.getX() * bd_size;
                                  rhs->operator()(global_id, 1) = d_point.m_val.getY() * bd_size;
                                  rhs->operator()(global_id, 2) = d_point.m_val.getZ() * bd_size;
                              }
                          }
        );


        rhs->move_global_to_local(m_sparse_system.get_l2g_map());
    }

    template<typename T>
    void GPSplineSolidCPU<T>::debug_init(const vector<CBSplineSolid> &patch, vector<boundary_point> &dirichlet_points,
                                         vector<boundary_point> &neumann_points) {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;
    }

    template<typename T>
    void
    GPSplineSolidCPU<T>::debug_assembly(const vector<CBSplineSolid> &patch, vector<boundary_point> &dirichlet_points,
                                        vector<boundary_point> &neumann_points) {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = 1000.0000000000000;
        double pr = 0.3;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);
        //linear_stiff_assembly<T>(patch, get_assembly_context(), lambda, mu);
        hyperelasticity_stiff_assembly<T>(patch, get_assembly_context(), get_physics_context(), lambda, mu);

        add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());

        stiffAssembler stiff(patch, dirichlet_points, neumann_points);
        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
                displacement, true, test_localMatsVec, test_localRhsVec);

        for (uint32_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
            auto local_mat_eigen = test_localMatsVec[geometry_id];
            auto local_mat_cpu = get_assembly_context().get_assembly_info(geometry_id).local_mats;

            for (uint16_t i = 0; i < 192; ++i) {
                for (uint16_t j = 0; j < 192; ++j) {
                    if (abs(local_mat_eigen.coeff(i, j) - local_mat_cpu[i + j * 192]) > 1e-10) {
                        GPMESH_INFO("Error for {},{} = {}", i, j, local_mat_cpu[i + j * 192]);
                    }
                }
            }
        }
        //add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());
    }

    template<typename T>
    void
    GPSplineSolidCPU<T>::debug_attribute_l2g(const vector<CBSplineSolid> &patch,
                                             vector<boundary_point> &dirichlet_points,
                                             vector<boundary_point> &neumann_points) {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = 1000.0000000000000;
        double pr = 0.3;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);

        auto rhs = m_attr_container.add_spline_attr("rhs", m_num_geometry, 64, 3, m_num_cpoints, gpmesh::SoA);
        rhs->reset_global_part(0);
        rhs->move_global_to_local(m_sparse_system.get_l2g_map());
        rhs->move_local_to_global(m_sparse_system.get_g2l_map());
        rhs->reset_local_part(0);

        rhs->reset_local_part(1.0);
        rhs->move_local_to_global_serial(m_sparse_system.get_g2l_map());
        for (uint32_t global_id = 0; global_id < m_num_cpoints; ++global_id) {
            if (m_sparse_system.get_g2l_map()[global_id].size() > 1) {
                for (int d = 0; d < 3; ++d) {
                    if (abs(rhs->RawAccessGlobal(global_id + d * m_num_cpoints) -
                            m_sparse_system.get_g2l_map()[global_id].size()) > 1e-13) {
                        GPMESH_INFO("Error Dirichlet idx={},d={} : {}", global_id, d,
                                    rhs->RawAccessGlobal(global_id + d * m_num_cpoints));
                    }
                }
            } else {
                for (int d = 0; d < 3; ++d) {
                    if (abs(rhs->RawAccessGlobal(global_id + d * m_num_cpoints) - 1) > 1e-13) {
                        GPMESH_INFO("Error idx={},d={} : {}", global_id, d,
                                    rhs->RawAccessGlobal(global_id + d * m_num_cpoints));
                    }
                }
            }
        }
    }

    template<typename T>
    void
    GPSplineSolidCPU<T>::debug_attribute_g2l(const vector<CBSplineSolid> &patch,
                                             vector<boundary_point> &dirichlet_points,
                                             vector<boundary_point> &neumann_points) {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = 1000.0000000000000;
        double pr = 0.3;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        auto rhs = m_attr_container.add_spline_attr("rhs", m_num_geometry, 64, 3, m_num_cpoints, gpmesh::SoA);

        Eigen::VectorXd test = Eigen::VectorXd::Random(m_num_cpoints * 3);

        rhs->reset_global_part(0);

        for (int i = 0; i < m_num_cpoints * 3; ++i) {
            rhs->RawAccessGlobal(i) = test(i);
        }

        GPMESH_INFO("squared norm = {}", Spline_Attribute_squared_norm(rhs));

        rhs->move_global_to_local(m_sparse_system.get_l2g_map());

        for (uint32_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
            for (uint16_t p = 0; p < 192; ++p) {
                uint16_t local_id = p % 64;
                uint16_t d = p / 64;

                uint32_t global_id = m_sparse_system.get_l2g_map()[geometry_id][local_id];

                if (abs(rhs->RawAccessLocal(geometry_id, p) - rhs->RawAccessGlobal(global_id + d * m_num_cpoints)) >
                    1e-13) {
                    GPMESH_INFO("Error for g={} l={} d={} : ", geometry_id, local_id, d,
                                rhs->RawAccessLocal(geometry_id, p));
                }
            }
        }

        T val = 0;
        //T val = 0.0;
        rhs->reset_global_part(val);

        for (uint32_t global_id = 0; global_id < m_num_cpoints * 3; ++global_id) {
            if (abs(rhs->RawAccessGlobal(global_id) - val) > 1e-13) {
                GPMESH_WARN("Error id={} : {}", global_id, rhs->RawAccessGlobal(global_id));
            }
        }

        rhs->move_local_to_global(m_sparse_system.get_g2l_map());

        for (uint32_t global_id = 0; global_id < m_num_cpoints; ++global_id) {
            for (int d = 0; d < 3; ++d) {
                T sum = val;
                for (auto sub: m_sparse_system.get_g2l_map()[global_id]) {
                    uint32_t patch_id = sub / 64;
                    uint32_t local_id = sub - 64 * patch_id;

                    sum += rhs->RawAccessLocal(patch_id, local_id + d * 64);
                }

                if (abs(rhs->RawAccessGlobal(global_id + d * m_num_cpoints) - sum) > 1e-13) {
                    GPMESH_INFO("ERROR g={}, d={}, sum={}, raw={}", global_id, d, sum,
                                rhs->RawAccessGlobal(global_id + d * m_num_cpoints));
                }
            }
        }

        for (uint32_t global_id = 0; global_id < m_num_cpoints; ++global_id) {
            int sub_count = m_sparse_system.get_g2l_map()[global_id].size();

            for (int d = 0; d < 3; ++d) {
                uint32_t id = global_id + d * m_num_cpoints;

                if (abs(rhs->RawAccessGlobal(id) - test[id] * sub_count - val) > 1e-13) {
                    GPMESH_INFO("gid={}, d={} \n raw={}, test={} \n mul={}, sub={}",
                                global_id, d, rhs->RawAccessGlobal(id), test[id], rhs->RawAccessGlobal(id) / test[id],
                                sub_count
                    );
                }
            }
        }

    }

    template<typename T>
    void
    GPSplineSolidCPU<T>::debug_solving(const vector<CBSplineSolid> &patch,
                                       std::vector<boundary_point> &dirichlet_points,
                                       std::vector<boundary_point> &neumann_points,
                                       Eigen::Matrix<T, -1, -1> &U) {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = 1000.0000000000000;
        double pr = 0.3;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);
        //linear_stiff_assembly<T>(patch, get_assembly_context(), lambda, mu);
        hyperelasticity_stiff_assembly<T>(patch, get_assembly_context(), get_physics_context(), lambda, mu);

        auto rhs = m_attr_container.add_spline_attr("rhs", m_num_geometry, 64, 3, m_num_cpoints, gpmesh::SoA);
        rhs->reset_global_part(0);
        rhs->move_global_to_local(m_sparse_system.get_l2g_map());
        rhs->move_local_to_global(m_sparse_system.get_g2l_map());
        rhs->reset_local_part(0);

        build_rhs(rhs, dirichlet_points, neumann_points);

        stiffAssembler stiff(patch, dirichlet_points, neumann_points);
        auto dispSol = Eigen::MatrixXd::Zero(stiff.numDofs(), 1);
        Eigen::MatrixXd displacement;
        stiff.constructCoefs(dispSol, displacement);
        stiff.system().matrix().setZero();
        stiff.system().rhs().setZero();
        std::vector<Eigen::MatrixXd> test_localMatsVec;
        std::vector<Eigen::MatrixXd> test_localRhsVec;
        stiff.test_stiffNonlinearMatrix(
                displacement, true, test_localMatsVec, test_localRhsVec);

        //local_mat test
        for (uint32_t geometry_id = 0; geometry_id < m_num_geometry; ++geometry_id) {
            auto local_mat_eigen = test_localMatsVec[geometry_id];
            auto local_mat_cpu = get_assembly_context().get_assembly_info(geometry_id).local_mats;

            for (uint16_t i = 0; i < 192; ++i) {
                for (uint16_t j = 0; j < 192; ++j) {
                    if (abs(local_mat_eigen.coeff(i, j) - local_mat_cpu[i + j * 192]) > 1e-10) {
                        GPMESH_WARN("Error for {},{} = {}", i, j, local_mat_cpu[i + j * 192]);
                    }
                }
            }
        }


        bool addDirichlet = true;
        if (addDirichlet)
            stiff.system().dataPush_Dirichlet_V2();
        else
            stiff.system().dataPush();

        //get matrix from stiff assembler
        m_sparse_system.matrix() = stiff.system().matrix();

        if (addDirichlet)
            add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());

        //build Sparse Matrix From assembly context
        stiff.system().matrix().setZero();
        for (uint32_t g = 0; g < m_num_geometry; ++g) {
            stiff.system().pushToMatrix(g, get_assembly_context().get_assembly_info(g).local_mats);
        }
        stiff.system().dataPush();

        //using sparse matrix from assembly context
        //m_sparse_system.matrix() = stiff.system().matrix();

        GPMESH_INFO("Matrix Diff = {}", (m_sparse_system.matrix() - stiff.system().matrix()).norm());

        bool local_debug_matrix = false;

        if (local_debug_matrix) {
            GPMESH_INFO("num = {}", get_dirichlet_piont_context().get_num_dirichlet_point());
            for (uint32_t d = 0; d < get_dirichlet_piont_context().get_num_dirichlet_point(); ++d) {
                auto info = get_dirichlet_piont_context().get_dirichlet_point_info(d);

                auto gid = info.geometry_id;
                auto lid = info.local_id;

                //GPMESH_INFO("gid={}, lid={}, d = {}", gid, lid, d);

                auto assembly_info = get_assembly_context().get_assembly_info(gid);

                for (uint16_t iter = 0; iter < 192; ++iter) {
                    for (uint32_t d1 = 0; d1 < 3; ++d1) {
                        auto clid = lid + d1 * 64;

                        if (abs(assembly_info.local_mats[clid + 192 * iter]) > 1e-10)
                            //if(clid == iter)
                        {
                            //GPMESH_INFO("Geometry = {}, Matrix[{}][{}] = {}", gid, clid, iter, assembly_info.local_mats[clid + 192 * iter]);
                        }
                    }
                }
            }

            for (uint32_t cpoint_id = 0; cpoint_id < m_num_cpoints; ++cpoint_id) {
                for (uint16_t d = 0; d < 3; ++d) {
                    uint32_t global_id = cpoint_id + d * m_num_cpoints;

                    if (m_sparse_system.IsDirichlet(global_id)) {
                        for (uint32_t j = 0; j < 3 * m_num_cpoints; ++j) {
                            if (abs(stiff.system().matrix().coeff(global_id, j)) > 1e-13) {
                                GPMESH_INFO("Matrix[{}][{}] = {}", global_id, j,
                                            stiff.system().matrix().coeff(global_id, j));
                            }
                        }
                    }
                }
            }
        }

        Eigen::Matrix<T, -1, -1> comp_result;
        float comp_pcg_time = 0.0;

        helper_Comp_PCG_Solver(comp_result, comp_pcg_time);

        Eigen::Matrix<T, -1, -1> cpu_result;
        float cpu_pcg_time = 0.0;
        helper_PCG_Solver_JacobiPrecondition(cpu_result, cpu_pcg_time);

        U = cpu_result;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::linear_elasticity_pcg_solving(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                       std::vector<boundary_point> &neumann_points,
                                       Eigen::Matrix<T, -1, -1> &U,
                                       T Yangs, T possion)
    {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = Yangs;
        double pr = possion;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);
        //linear_stiff_assembly<T>(patch, get_assembly_context(), lambda, mu);
        hyperelasticity_stiff_assembly<T>(patch, get_assembly_context(), get_physics_context(), lambda, mu);

        auto rhs = m_attr_container.add_spline_attr("rhs", m_num_geometry, 64, 3, m_num_cpoints, gpmesh::SoA);
        rhs->reset_global_part(0);
        rhs->move_global_to_local(m_sparse_system.get_l2g_map());
        rhs->move_local_to_global(m_sparse_system.get_g2l_map());
        rhs->reset_local_part(0);

        build_rhs(rhs, dirichlet_points, neumann_points);

        add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());

        Eigen::Matrix<T, -1, -1> cpu_result;
        float cpu_pcg_time = 0.0;
        helper_PCG_Solver_JacobiPrecondition(cpu_result, cpu_pcg_time);

        U = cpu_result;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::linear_elasticity_pardiso(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                                            std::vector<boundary_point> &neumann_points,
                                                            Eigen::Matrix<T, -1, -1> &U,
                                                            T Yangs, T possion)
    {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = Yangs;
        double pr = possion;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);
        //linear_stiff_assembly<T>(patch, get_assembly_context(), lambda, mu);
        hyperelasticity_stiff_assembly<T>(patch, get_assembly_context(), get_physics_context(), lambda, mu);

        add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());

        for (uint32_t g = 0; g < m_num_geometry; ++g)
        {
            m_sparse_system.pushToMatrix(g, get_assembly_context().get_assembly_info(g).local_mats);
        }

        m_sparse_system.dataPush();
        m_sparse_system.rhs().setZero();

        helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);
        helper_set_neumann_point_in_sparse_system_rhs(neumann_points);


        Eigen::Matrix<T, -1, -1> cpu_result;
        float cpu_pcg_time = 0.0;
        helper_pardiso_solver(cpu_result, cpu_pcg_time, 1e-13, false);

        U = cpu_result;
    }
    template<typename T>
    void GPSplineSolidCPU<T>::linear_elasticity_Eigen(const vector<CBSplineSolid> &patch, std::vector<boundary_point> &dirichlet_points,
                                 std::vector<boundary_point> &neumann_points,
                                 Eigen::Matrix<T, -1, -1> &U,
                                 T Yangs, T possion, bool bUseCG)

    {
        init(patch, dirichlet_points, neumann_points, false);

        double timeSpan = 1;
        double timeStep = 0.1;
        double density = 3320;

        double e = Yangs;
        double pr = possion;

        double lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
        double mu = e / (2. * (1. + pr));

        double Beta = 0.25;
        double Gamma = 0.5;
        double tStep = 0.1;

        fill_assembly_context(0);
        fill_physics_context(0);
        //linear_stiff_assembly<T>(patch, get_assembly_context(), lambda, mu);
        hyperelasticity_stiff_assembly<T>(patch, get_assembly_context(), get_physics_context(), lambda, mu);

        add_dirichlet_point_info_to_assembly<T>(get_assembly_context(), get_dirichlet_piont_context());

        for (uint32_t g = 0; g < m_num_geometry; ++g)
        {
            m_sparse_system.pushToMatrix(g, get_assembly_context().get_assembly_info(g).local_mats);
        }

        m_sparse_system.dataPush();
        m_sparse_system.rhs().setZero();

        helper_set_dirichlet_point_in_sparse_system_rhs(dirichlet_points);
        helper_set_neumann_point_in_sparse_system_rhs(neumann_points);


        Eigen::Matrix<T, -1, -1> cpu_result;
        float cpu_pcg_time = 0.0;
        helper_eigen_solver(cpu_result, cpu_pcg_time, bUseCG, 1e-13, false);

        U = cpu_result;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::helper_PCG_Solver_JacobiPrecondition(Eigen::Matrix<T, -1, -1> &result, float &time,
                                                                   T threshold, bool bQuiet) {

        CPUTimer timer;

        auto B = m_attr_container.get_attr("rhs");

        if (!bQuiet) GPMESH_INFO("rhs norm = {}", Spline_Attribute_squared_norm(B));

        auto X = helper_add_spline_attribute("X", 3, gpmesh::SoA);
        X->reset_global_part(0);
        X->reset_local_part(0);

        auto S = helper_add_spline_attribute("S", 3, gpmesh::SoA);
        S->reset_global_part(0);
        S->reset_local_part(0);

        auto R = helper_add_spline_attribute("R", 3, gpmesh::SoA);
        R->reset_global_part(0);
        R->reset_local_part(0);

        auto P = helper_add_spline_attribute("P", 3, gpmesh::SoA);
        P->reset_global_part(0);
        P->reset_local_part(0);

        auto Preconditioner = helper_add_spline_attribute("Preconditioner", 3, gpmesh::SoA);
        Preconditioner->reset_global_part(0);
        Preconditioner->reset_local_part(0);

        helper_build_jacobi_preconditioner_from_assembly(Preconditioner);

        if (!bQuiet) GPMESH_INFO("Preconditioner norm = {}", Spline_Attribute_squared_norm(Preconditioner));

        matvec(get_assembly_context(), m_num_geometry, X, S);

        S->move_local_to_global(m_sparse_system.get_g2l_map());

        init_P_R(B, S, P, R);


        Spline_Attribute_CWise_dot(P, Preconditioner);
        P->move_global_to_local(m_sparse_system.get_l2g_map());

        T Delta_new = Spline_Attribute_dot(R, P);
        T Delta_0 = Delta_new;

        if (!bQuiet) GPMESH_INFO("Delta_0 = {}", Delta_0);

        uint32_t iter_num = 0;
        uint32_t iter_limit = 10000;

        T Alpha(0), Beta(0), Delta_old(0);

        timer.start();

        while (iter_num < iter_limit) {
            S->reset_local_part(0);
            S->reset_global_part(0);

            if (!bQuiet) GPMESH_INFO("S squared norm before = {}", Spline_Attribute_squared_norm(S));
            matvec(get_assembly_context(), m_num_geometry, P, S);
            S->move_local_to_global(m_sparse_system.get_g2l_map());

            if (!bQuiet) GPMESH_INFO("P squared norm = {}", Spline_Attribute_squared_norm(P));
            if (!bQuiet) GPMESH_INFO("S squared norm = {}", Spline_Attribute_squared_norm(S));

            Alpha = Spline_Attribute_dot(S, P);
            if (!bQuiet) GPMESH_INFO("dot(S,P) = {}", Alpha);

            Alpha = Delta_new / Alpha;
            if (!bQuiet) GPMESH_INFO("Alpha={}", Alpha);

            Spline_Attribute_axpy<T>(X, P, Alpha, 1.0);
            Spline_Attribute_axpy<T>(R, S, -Alpha, 1.0);

            Delta_old = Delta_new;

            Spline_Attribute_CWise_dot(S, Preconditioner, R);

            Delta_new = Spline_Attribute_dot(S, R);
            if (!bQuiet) GPMESH_INFO("Delta New = {}", Delta_new);

            if (Delta_new < threshold * Delta_0) {
                break;
            }

            Beta = Delta_new / Delta_old;
            if (!bQuiet) GPMESH_INFO("Beta = {}", Beta);

            Spline_Attribute_axpy<T>(P, S, 1.0, Beta);
            if (!bQuiet) GPMESH_INFO("1 P squared norm = {}", Spline_Attribute_squared_norm(P));
            P->move_global_to_local(m_sparse_system.get_l2g_map());

            ++iter_num;
        }
        timer.stop();

        time = timer.elapsed_millis();
        GPMESH_INFO("Parallel PCG Solver Takes {}(ms)", time);

        GPMESH_INFO("PCG CPU Takes {} rounds iterations", iter_num);

        X->move_local_to_global(m_sparse_system.get_g2l_map());

        Eigen::Matrix<T, -1, -1> updateVector_cps;
        updateVector_cps.setZero(m_num_cpoints * 3, 1);

        for (int d = 0; d < m_num_cpoints * 3; ++d) {
            updateVector_cps(d, 0) = X->RawAccessGlobal(d);
        }

        result = updateVector_cps;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::helper_Comp_PCG_Solver(Eigen::Matrix<T, -1, -1> &result, float &time, T threshold,
                                                     bool bQuiet) {
        Eigen::SparseMatrix<double> &A = m_sparse_system.matrix();

        auto rhs = m_attr_container.get_attr("rhs");

        uint32_t globalSize = rhs->m_num_attributes * rhs->m_num_cpoints;
        Eigen::VectorXd rhsVec = Eigen::VectorXd::Zero(globalSize);

        Eigen::VectorXd X = Eigen::VectorXd::Zero(globalSize);

        for (uint32_t i = 0; i < rhsVec.rows(); ++i) {
            rhsVec(i) = rhs->RawAccessGlobal(i);
        }

        if (!bQuiet) GPMESH_INFO("rhs norm Comp = {}", rhsVec.squaredNorm());

        Eigen::VectorXd Preconditioner = Eigen::VectorXd::Zero(rhs->m_num_attributes * rhs->m_num_cpoints);

        {
            for (uint32_t i = 0; i < globalSize; ++i) {
                T val = A.coeff(i, i);
                if (abs(val) < 1e-10) {
                    val = 1;
                } else {
                    val = static_cast<T>(1.0) / val;
                }

                Preconditioner(i) = val;
            }
        }

        if (!bQuiet) GPMESH_INFO("Preconditioner squaredNorm Comp = {}", Preconditioner.squaredNorm());

        Eigen::VectorXd S = A * X;

        //init P R
        Eigen::VectorXd R = rhsVec - S;
        Eigen::VectorXd P = R;

        //P = (P * Preconditioner.transpose()).diagonal();
        P = P.cwiseProduct(Preconditioner);

        double Delta_new = R.dot(P);
        double Delta_0 = Delta_new;
        if (!bQuiet) GPMESH_INFO("Delta_0 Comp = {}", Delta_0);

        uint32_t iter_num = 0;
        uint32_t iter_limit = 10000;

        double Alpha(0), Beta(0), Delta_old(0);
        CPUTimer timer;
        timer.start();
        while (iter_num < iter_limit) {
            S = A * P;

            if (!bQuiet) GPMESH_INFO("P squaredNorm Comp = {}", P.squaredNorm());
            if (!bQuiet) GPMESH_INFO("S squaredNorm Comp = {}", S.squaredNorm());

            Alpha = S.dot(P);

            if (!bQuiet) GPMESH_INFO("dot(S,P) Comp = {}", Alpha);

            Alpha = Delta_new / Alpha;

            if (!bQuiet) GPMESH_INFO("Alpha Comp ={} ", Alpha);

            axpy<decltype(X), double>(X, P, Alpha, 1.0);
            axpy<decltype(R), double>(R, S, -Alpha, 1.0);

            Delta_old = Delta_new;

            S = Preconditioner.cwiseProduct(R);

            Delta_new = S.dot(R);

            if (!bQuiet) GPMESH_INFO("Delta New Comp = {}", Delta_new);

            if (Delta_new < threshold * Delta_0) {
                break;
            }

            Beta = Delta_new / Delta_old;
            if (!bQuiet) GPMESH_INFO("Beta Comp = {}", Beta);

            axpy<decltype(P), double>(P, S, 1.0, Beta);

            ++iter_num;
        }
        timer.stop();
        GPMESH_INFO("Serial PCG Solver Takes {}(ms)", timer.elapsed_millis());
        GPMESH_INFO("Serial PCG CPU Takes {} rounds iterations", iter_num);
    }

    template<typename T>
    void GPSplineSolidCPU<T>::helper_pardiso_solver(Eigen::Matrix<T, -1, -1> &result, float &time,
                                                                   T threshold, bool bQuiet)
    {
        CPUTimer timer;
        timer.start();
        Eigen::PardisoLDLT<Eigen::SparseMatrix<double>> pardisoSolve;

        // cg.setMaxIterations(50);
        pardisoSolve.compute(m_sparse_system.matrix());
        result = pardisoSolve.solve(m_sparse_system.rhs()).cast<T>();
        if (pardisoSolve.info() == Eigen::NoConvergence ||
            pardisoSolve.info() == Eigen::NumericalIssue) {
            if (!bQuiet)
                GPMESH_TRACE("Pardiso LDLT Solve Fail!");
        } else {
            if (!bQuiet)
                GPMESH_TRACE("Pardiso LDLT Solve Success!");
        }

        timer.stop();
        time = timer.elapsed_millis();
        GPMESH_INFO("Pardiso LDLT Takes {}(ms)", time);
    }
    template<typename T>
    void GPSplineSolidCPU<T>::helper_eigen_solver(Eigen::Matrix<T, -1, -1> &result, float &time, bool bUseCG, T threshold,
                             bool bQuiet)

    {
        Eigen::initParallel();
        Eigen::setNbThreads(24);

        CPUTimer timer;
        timer.start();

        Eigen::ConjugateGradient<Eigen::SparseMatrix<double>> eigen_CG_Solver;
        Eigen::SimplicialLDLT<Eigen::SparseMatrix<double>> eigen_LDLT_Solver;

        // cg.setMaxIterations(50);
        if(bUseCG)
        {
            eigen_CG_Solver.compute(m_sparse_system.matrix());
            result = eigen_CG_Solver.solve(m_sparse_system.rhs()).cast<T>();

            if (eigen_CG_Solver.info() == Eigen::NoConvergence ||
                    eigen_CG_Solver.info() == Eigen::NumericalIssue) {
                if (!bQuiet) GPMESH_ERROR("Eigen CG Solve Fail!");
            } else {
                if (!bQuiet) GPMESH_TRACE("Eigen CG Solve Success!");
            }
        }
        else
        {
            eigen_LDLT_Solver.compute(m_sparse_system.matrix());
            result = eigen_LDLT_Solver.solve(m_sparse_system.rhs()).cast<T>();

            if (eigen_LDLT_Solver.info() == Eigen::NoConvergence ||
                    eigen_LDLT_Solver.info() == Eigen::NumericalIssue) {
                if (!bQuiet) GPMESH_ERROR("Eigen LDLT Solve Fail!");
            } else {
                if (!bQuiet) GPMESH_TRACE("Eigen LDLT Solve Success!");
            }
        }

        timer.stop();
        time = timer.elapsed_millis();
        GPMESH_INFO("Eigen Takes {}(ms)", time);
    }

    template<typename T>
    std::shared_ptr<SplineAttributeCPU<T>>
    GPSplineSolidCPU<T>::helper_add_spline_attribute(const std::string &name, const uint32_t num_attributes,
                                                     layoutT layout) {
        auto attr_ptr = m_attr_container.add_spline_attr(name.c_str(),
                                                         m_num_geometry,
                                                         64,
                                                         num_attributes,
                                                         m_num_cpoints,
                                                         layout);

        return attr_ptr;
    }

    template<typename T>
    void GPSplineSolidCPU<T>::helper_build_jacobi_preconditioner_from_assembly(
            std::shared_ptr<SplineAttributeCPU<T>> Preconditioner) {

        uint32_t arraySize = Preconditioner->m_num_attributes *
                             Preconditioner->m_cpoint_num_per_geometry;
        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_geometry),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t geometry_id = r.begin(); geometry_id < r.end(); ++geometry_id) {
                                  for (uint32_t i = 0; i < arraySize; ++i) {
                                      Preconditioner->RawAccessLocal(geometry_id, i) =
                                              get_assembly_context().get_assembly_info(geometry_id).local_mats[i +
                                                                                                               arraySize *
                                                                                                               i];
                                  }
                              }
                          }
        );

        Preconditioner->move_local_to_global(m_sparse_system.get_g2l_map());

        tbb::parallel_for(
                tbb::blocked_range<size_t>(0, Preconditioner->m_num_cpoints * Preconditioner->m_num_attributes),
                [&](tbb::blocked_range<size_t> r) {
                    for (size_t raw_id = r.begin(); raw_id < r.end(); ++raw_id) {
                        T val = Preconditioner->RawAccessGlobal(raw_id);
                        if (abs(val) < 1e-10) {
                            val = 1;
                        } else {
                            val = static_cast<T>(1.0) / val;
                        }

                        Preconditioner->RawAccessGlobal(raw_id) = val;
                    }
                }
        );

        Preconditioner->move_global_to_local(m_sparse_system.get_l2g_map());

    }


    template<typename T>
    void GPSplineSolidCPU<T>::fill_assembly_context(T fill_data) const {
        assert(m_num_geometry > 0);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_geometry),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t i = r.begin(); i < r.end(); ++i) {
                                  T *mat_data = m_assembly_context.get_assembly_info(i).local_mats;
                                  std::fill(mat_data, mat_data + 192 * 192, fill_data);
                              }
                          }
        );
    }

    template<typename T>
    void GPSplineSolidCPU<T>::fill_physics_context(T fill_data) const {
        assert(m_num_geometry > 0);

        tbb::parallel_for(tbb::blocked_range<size_t>(0, m_num_geometry),
                          [&](tbb::blocked_range<size_t> r) {
                              for (size_t i = r.begin(); i < r.end(); ++i) {
                                  m_physics_context.get_physics_info(i).m_cpoints_displacement =
                                          Eigen::Matrix<T, 64, 3>::Ones() * fill_data;
                              }
                          }
        );
    }


}  // namespace gpmesh
