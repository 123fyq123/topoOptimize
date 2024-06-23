#include <assert.h>
//#include <omp.h>
#include <iostream>
#include <memory>
#include <numeric>
#include <queue>

#include "Util/GP_Util.h"
#include "GP_Context.h"
#include "GP_Spline.h"

#include <Hyperelasticity/ToolFunc.h>

#include "Eigen_helper.h"

namespace gpmesh {

    GPSplineSolid::GPSplineSolid()
            : m_num_cpoints(0),
              m_max_cp_per_geometry(0),
              m_max_not_owned_cp(0),
              m_num_geometry(0),
              m_geometry_size(512),
              m_quite(false),
              m_h_geometry_info(nullptr),
              m_d_geometry_info(nullptr) {
    }

    void GPSplineSolid::init(const vector<CBSplineSolid> &patch, std::vector<boundary_condition> &m_dirichletSides,
                             std::vector<boundary_condition> &m_neumannSides, const bool quite) {
        m_quite = quite;

        if (patch.empty()) {
            GPMESH_ERROR(
                    "GPSpline::inti input patch is empty"
            );
        }

        m_sparse_system = SparseSystem(patch, m_dirichletSides, 3);
        m_num_cpoints = m_sparse_system.getControlNums();

        build(patch);
        build_device(patch);

        m_gpspline_context.init(m_num_cpoints, m_num_geometry, m_d_geometry_info);

        build_dirichlet_context(m_dirichletSides);

        GPMESH_INFO("GP Spline Init Success!");

        bInitPatch = true;
        bInitDirichlet = true;
    }

    void GPSplineSolid::init(const vector<CBSplineSolid> &patch,
                             vector<boundary_point> &dirichlet_points,
                             vector<boundary_point> &neumann_points,
                             const bool quite) {
        m_quite = quite;

        if (patch.empty()) {
            GPMESH_ERROR(
                    "GPSpline::inti input patch is empty"
            );
        }

        //m_sparse_system = SparseSystem(patch, m_dirichletSides, 3);
        m_sparse_system = SparseSystem(patch, dirichlet_points, 3);
        m_num_cpoints = m_sparse_system.getControlNums();
        std::cout << "ControlNum = " << m_num_cpoints << "\n";
        build(patch);
        build_device(patch);

        m_gpspline_context.init(m_num_cpoints, m_num_geometry, m_d_geometry_info);

        //build_dirichlet_context(m_dirichletSides);
        build_dirichlet_point_context(dirichlet_points);

        GPMESH_INFO("GP Spline Init Success!");

        bInitPatch = true;
        bInitDirichlet = true;
    }


    void GPSplineSolid::init_patch(const vector<CBSplineSolid> &patch) {
        m_sparse_system = SparseSystem();
        m_sparse_system.init_geometry(patch, 3);
        m_num_cpoints = m_sparse_system.getControlNums();
        build(patch);
        build_device(patch);

        m_gpspline_context.init(m_num_cpoints, m_num_geometry, m_d_geometry_info);
        GPMESH_INFO("GP Spline Geometry Init Success!");

        lambda = Young * Poisson / ((1. + Poisson) * (1. - 2. * Poisson));    //拉梅常数
        mu = Young / (2. * (1. + Poisson));

        bInitPatch = true;
    }

    void GPSplineSolid::init_dirichlet_points(const std::vector<boundary_point>& dirichlet_points)
    {
        m_sparse_system.init_dirichlet(dirichlet_points);

        build_dirichlet_point_context(dirichlet_points);
        GPMESH_INFO("GP Spline Dirichlet point Init Success!");

        bInitDirichlet = true;
    }

    void GPSplineSolid::setYoungPossion(double _Young, double _Poisson) {
        Young = _Young;
        Poisson = _Poisson;
        lambda = Young * Poisson / ((1. + Poisson) * (1. - 2. * Poisson));    //拉梅常数
        mu = Young / (2. * (1. + Poisson));
    }

    void GPSplineSolid::build(const vector<CBSplineSolid> &patch) {
        m_num_geometry = patch.size();
        m_h_geometry_ltog_cp.resize(m_num_geometry);
        m_h_num_owned_cp.resize(m_num_geometry);

        uint32_t cp_num = m_sparse_system.getControlNums();
        m_h_geometry_gtol_cp.resize(cp_num);

        for (int i = 0; i < cp_num; ++i) {
            for (int j = 0; j < m_sparse_system.get_g2l_map()[i].size(); ++j) {
                m_h_geometry_gtol_cp[i].push_back(m_sparse_system.get_g2l_map()[i][j]);
            }
        }

        for (int i = 0; i < m_num_geometry; ++i) {
            build_single_geometry(patch, i);
        }

        m_h_cpoints_displacement.resize(m_num_cpoints * m_dims, 1);
        m_h_cpoints_displacement.setZero();
    }

    void GPSplineSolid::build_single_geometry(const vector<CBSplineSolid> &patch,
                                              const uint32_t geometry_id) {
        auto &m_sparse_system_map = m_sparse_system.get_l2g_map();

        uint32_t geo_size = m_sparse_system_map[geometry_id].size();

        m_h_geometry_ltog_cp[geometry_id].resize(geo_size);

        for (int i = 0; i < geo_size; ++i) {
            m_h_geometry_ltog_cp[geometry_id][i] = static_cast<uint32_t>(m_sparse_system_map[geometry_id][i]);
        }
    }

    void GPSplineSolid::build_device(const std::vector<CBSplineSolid> &patch) {
        //gpu store every geometry cpoints information
        CUDA_ERROR(cudaMalloc((void **) &m_d_geometry_info,
                              m_num_geometry * sizeof(GeometryInfo)));

        m_h_geometry_info = (GeometryInfo *) malloc(m_num_geometry * sizeof(GeometryInfo));


#pragma omp parallel for
        for (int g = 0; g < m_num_geometry; ++g) {

            //GPU side
            GeometryInfo d_geometry_info;

            d_geometry_info.geometry_id = g;
            d_geometry_info.num_CPoints = patch[g].getSizeControlPoints();

            d_geometry_info.u_size = patch[g].getUSize();
            d_geometry_info.v_size = patch[g].getVSize();
            d_geometry_info.w_size = patch[g].getWSize();

            d_geometry_info.u_degree = patch[g].getUDegree();
            d_geometry_info.v_degree = patch[g].getVDegree();
            d_geometry_info.w_degree = patch[g].getWDegree();

            vector<double> h_uknots = patch[g].getUAllKnots();
            vector<double> h_vknots = patch[g].getVAllKnots();
            vector<double> h_wknots = patch[g].getWAllKnots();
            assert(h_uknots.size() == d_geometry_info.u_size);
            CUDA_ERROR(cudaMalloc((void **) &d_geometry_info.u_knots,
                                  h_uknots.size() * sizeof(double)));
            CUDA_ERROR(cudaMemcpy(d_geometry_info.u_knots,
                                  h_uknots.data(), h_uknots.size() * sizeof(double),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMalloc((void **) &d_geometry_info.v_knots,
                                  h_vknots.size() * sizeof(double)));
            CUDA_ERROR(cudaMemcpy(d_geometry_info.v_knots,
                                  h_vknots.data(), h_vknots.size() * sizeof(double),
                                  cudaMemcpyHostToDevice));
            CUDA_ERROR(cudaMalloc((void **) &d_geometry_info.w_knots,
                                  h_wknots.size() * sizeof(double)));
            CUDA_ERROR(cudaMemcpy(d_geometry_info.w_knots,
                                  h_wknots.data(), h_wknots.size() * sizeof(double),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMalloc((void **) &d_geometry_info.coefs,
                                  patch[g].getSizeControlPoints() * m_dims * sizeof(double)));

            vector<double> h_coefs(patch[g].getSizeControlPoints() * m_dims, 0);

            for (int i = 0; i < patch[g].getSizeControlPoints(); ++i) {
                h_coefs[3 * i + 0] = patch[g].getCoefs()[i].getX();
                h_coefs[3 * i + 1] = patch[g].getCoefs()[i].getY();
                h_coefs[3 * i + 2] = patch[g].getCoefs()[i].getZ();
            }

            CUDA_ERROR(cudaMemcpy(d_geometry_info.coefs,
                                  h_coefs.data(),
                                  patch[g].getSizeControlPoints() * m_dims * sizeof(double),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_geometry_info + g,
                                  &d_geometry_info,
                                  sizeof(GeometryInfo),
                                  cudaMemcpyHostToDevice));


            vector<double> h_cps_disp_per_geometry(patch[g].getSizeControlPoints() * m_dims, 0);


            //cpu side
            m_h_geometry_info[g].geometry_id = g;
            m_h_geometry_info[g].num_CPoints = patch[g].getSizeControlPoints();

            m_h_geometry_info[g].u_size = patch[g].getUSize();
            m_h_geometry_info[g].v_size = patch[g].getVSize();
            m_h_geometry_info[g].w_size = patch[g].getWSize();

            m_h_geometry_info[g].u_degree = patch[g].getUDegree();
            m_h_geometry_info[g].v_degree = patch[g].getVDegree();
            m_h_geometry_info[g].w_degree = patch[g].getWDegree();

            //not fill data yet
            m_h_geometry_info[g].u_knots = (double *) malloc(patch[g].getUSize() * sizeof(double));
            m_h_geometry_info[g].v_knots = (double *) malloc(patch[g].getVSize() * sizeof(double));
            m_h_geometry_info[g].w_knots = (double *) malloc(patch[g].getWSize() * sizeof(double));

            m_h_geometry_info[g].coefs = (double *) malloc(h_coefs.size() * sizeof(double));

            //Gaussian Integral Context;

        }

        build_gaussian_integral_context(patch);
        build_physic_context(patch);
        build_assembly_context(patch);
        build_attr_mapping_context();
    }

    GPSplineSolid::~GPSplineSolid() {
        for (int g = 0; g < m_num_geometry; ++g) {
            free(m_h_geometry_info[g].u_knots);
            free(m_h_geometry_info[g].v_knots);
            free(m_h_geometry_info[g].w_knots);
            free(m_h_geometry_info[g].coefs);
        }

        CUDA_ERROR(cudaMemcpy(m_h_geometry_info,
                              m_d_geometry_info,
                              m_num_geometry * sizeof(GeometryInfo),
                              cudaMemcpyDeviceToHost));

        for (int g = 0; g < m_num_geometry; ++g) {
            GPU_FREE(m_h_geometry_info[g].u_knots);
            GPU_FREE(m_h_geometry_info[g].v_knots);
            GPU_FREE(m_h_geometry_info[g].w_knots);
            GPU_FREE(m_h_geometry_info[g].coefs);
        }

        GPU_FREE(m_d_geometry_info);
        free(m_h_geometry_info);

        delete[] m_h_cpoints_displacement_array;

        m_gpspline_context.release();
        m_gs_integral_context.release();
        m_physic_context.release();
        m_dirichlet_context.release();
    }

    void GPSplineSolid::build_gaussian_integral_context(const std::vector<CBSplineSolid> &patch) {
        std::vector<double> p_u, p_v, p_w; // �洢�������˹��
        std::vector<double> w_u, w_v, w_w; // �洢�������˹��Ȩֵ

        int gaussNum_u = m_h_geometry_info[0].u_degree + 1, gaussNum_v = m_h_geometry_info[0].u_degree + 1,
                gaussNum_w = m_h_geometry_info[0].w_degree + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;

        p_u.resize(gaussNum_u);
        p_v.resize(gaussNum_v);
        p_v.resize(gaussNum_w);
        w_u.resize(gaussNum_u);
        w_v.resize(gaussNum_v);
        w_v.resize(gaussNum_w);

        ToolFunc::gaussIntegalPoint(gaussNum_u, 0, 1, p_u, w_u);
        ToolFunc::gaussIntegalPoint(gaussNum_v, 0, 1, p_v, w_v);
        ToolFunc::gaussIntegalPoint(gaussNum_w, 0, 1, p_w, w_w);

        m_gs_integral_context.init(gaussNum_u, gaussNum_v, gaussNum_w,
                                   p_u.data(), p_v.data(), p_w.data(),
                                   w_u.data(), w_v.data(), w_w.data());

        std::vector<Eigen::MatrixXd> basisValue;

        CBSplineSolid *geo = const_cast<CBSplineSolid *>(&patch[0]);
        geo->eval_AllDers(p_u, p_v, p_w, basisValue);

        std::vector<double> basis_value_0, basis_value_1;

        //basis_value_0.resize(basisValue[0].rows() * basisValue[0].cols());
        //basis_value_1.resize(basisValue[1].rows() * basisValue[1].cols());

        convert_eigen_matrix_to_vector(basisValue[0], basis_value_0);
        convert_eigen_matrix_to_vector(basisValue[1], basis_value_1);
        //memory leak may cause strange bug
        m_gs_integral_context.init_derivatives(basis_value_0.data(), basis_value_0.size(), basis_value_1.data(),
                                               basis_value_1.size());
        m_gs_integral_context.init_derivatives_matrix(basisValue[0], basisValue[1]);
    }

    void GPSplineSolid::build_physic_context(const std::vector<CBSplineSolid> &patch) {
        double E = 1e3;
        double pr = 0.3;

        double lambda = E * pr / ((1. + pr) * (1. - 2. * pr));    //��÷����
        double mu = E / (2. * (1. + pr));
        Eigen::MatrixXd identity = Eigen::Matrix<double, 3, 3>::Identity();
        Eigen::MatrixXd C;

        ToolFunc::setC(C, identity, lambda, mu);
        Eigen::Matrix<double, 6, 6> C_copy = C;

        m_physic_context.init_C(C_copy);
        m_physic_context.init(m_h_cpoints_displacement, m_h_coefs_displacement, m_sparse_system, m_num_geometry);

        m_h_cpoints_displacement_array = new Eigen::Matrix<double, 64, 3>[m_num_geometry];
        for (int g = 0; g < m_num_geometry; ++g) {
            m_h_cpoints_displacement_array[g] = Eigen::Matrix<double, 64, 3>::Zero();
        }
    }

    void GPSplineSolid::build_assembly_context(const std::vector<CBSplineSolid> &patch) {
        AssemblyInfo *m_d_assembly_info;
        MassInfo *m_d_mass_info;

        CUDA_ERROR(cudaMalloc((void **) &m_d_assembly_info,
                              m_num_geometry * sizeof(AssemblyInfo)));
        m_h_assembly_info = (AssemblyInfo *) malloc(m_num_geometry * sizeof(AssemblyInfo));

        CUDA_ERROR(cudaMalloc((void **) &m_d_mass_info,
                              m_num_geometry * sizeof(MassInfo)));

        m_h_mass_info = (MassInfo *) malloc(m_num_geometry * sizeof(MassInfo));

        for (int g = 0; g < m_num_geometry; ++g) {
            AssemblyInfo d_assembly_info;

            d_assembly_info.geometry_id = g;

            CUDA_ERROR(cudaMalloc((void **) &d_assembly_info.local_mats,
                                  192 * 192 * sizeof(double)));

            std::vector<double> zero_array(192 * 192, 0.0);

            CUDA_ERROR(cudaMemcpy(d_assembly_info.local_mats,
                                  zero_array.data(), 192 * 192 * sizeof(double), cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMalloc((void **) &d_assembly_info.local_rhs,
                                  64 * 3 * sizeof(double)));

            std::vector<double> zero_rhs(192, 0);

            CUDA_ERROR(cudaMemcpy(d_assembly_info.local_rhs,
                                  zero_rhs.data(), 64 * 3 * sizeof(double), cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_assembly_info + g, &d_assembly_info,
                                  sizeof(AssemblyInfo),
                                  cudaMemcpyHostToDevice));

            m_h_assembly_info[g].geometry_id = g;
            m_h_assembly_info[g].local_mats = (double *) malloc(192 * 192 * sizeof(double));
            m_h_assembly_info[g].local_rhs = (double *) malloc(192 * sizeof(double));

            MassInfo d_mass_info;
            d_mass_info.geometry_id = g;

            CUDA_ERROR(cudaMalloc((void **) &d_mass_info.mass_local_mats,
                                  192 * 192 * sizeof(double)));

            CUDA_ERROR(cudaMemcpy(d_mass_info.mass_local_mats,
                                  zero_array.data(), 192 * 192 * sizeof(double), cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_mass_info + g, &d_mass_info,
                                  sizeof(MassInfo),
                                  cudaMemcpyHostToDevice));

            m_h_mass_info[g].geometry_id = g;
            m_h_mass_info[g].mass_local_mats = (double *) malloc(192 * 192 * sizeof(double));
        }

        m_assembly_context.init(m_d_assembly_info);
        m_mass_context.init(m_d_mass_info);
    }

    void GPSplineSolid::build_dirichlet_context(
            const vector<boundary_condition> &m_dirichletSides) {
        uint32_t num_dirichlet_boundary = m_dirichletSides.size();

        DirichletInfo *m_h_dirichlet_info = (DirichletInfo *) malloc(num_dirichlet_boundary * sizeof(DirichletInfo));

        for (uint32_t dirichlet_id = 0; dirichlet_id < num_dirichlet_boundary; ++dirichlet_id) {
            m_h_dirichlet_info[dirichlet_id].geometry_id = m_dirichletSides[dirichlet_id].side().patch();

            BoundarySide side;

            switch (m_dirichletSides[dirichlet_id].side().index()) {
                case boundary::south:
                    side = BoundarySide::SOUTH;
                    break;
                case boundary::north:
                    side = BoundarySide::NORTH;
                    break;
                case boundary::west:
                    side = BoundarySide::WEST;
                    break;
                case boundary::east:
                    side = BoundarySide::EAST;
                    break;
                case boundary::front:
                    side = BoundarySide::FRONT;
                    break;
                case boundary::back:
                    side = BoundarySide::BACK;
                    break;
            }

            m_h_dirichlet_info[dirichlet_id].side = side;

            m_h_dirichlet_info[dirichlet_id].dim = m_dirichletSides[dirichlet_id].getDir();

            m_h_dirichlet_info[dirichlet_id].value = m_dirichletSides[dirichlet_id].getValue();

            m_h_dirichlet_info[dirichlet_id].dirichlet_boundary_id = dirichlet_id;
        }

        m_dirichlet_context.init(m_h_dirichlet_info, num_dirichlet_boundary);
    }

    void GPSplineSolid::build_dirichlet_point_context(const vector<boundary_point> &dirichlet_points) {
        uint32_t num_dirichlet_point = dirichlet_points.size();

        DirichletPointInfo *m_h_dirichlet_point_info = (DirichletPointInfo *) malloc(
                num_dirichlet_point * sizeof(DirichletPointInfo));

        for (uint32_t dirichlet_point_id = 0; dirichlet_point_id < num_dirichlet_point; ++dirichlet_point_id) {
            auto &bd_point = dirichlet_points[dirichlet_point_id];

            assert(bd_point.m_condition_type == 0);

            uint32_t bp_global_id = bd_point.m_cpoint_global_id;

            m_h_dirichlet_point_info[dirichlet_point_id].valueX = bd_point.m_val.getX();
            m_h_dirichlet_point_info[dirichlet_point_id].valueY = bd_point.m_val.getY();
            m_h_dirichlet_point_info[dirichlet_point_id].valueZ = bd_point.m_val.getZ();

            m_h_dirichlet_point_info[dirichlet_point_id].global_id = bp_global_id;

            //获得local信息
            uint32_t patch_id = m_h_geometry_gtol_cp[bp_global_id][0] / 64;
            uint16_t local_id = m_h_geometry_gtol_cp[bp_global_id][0] - 64 * patch_id;
            assert(m_h_geometry_ltog_cp[patch_id][local_id] == bp_global_id);

            m_h_dirichlet_point_info[dirichlet_point_id].geometry_id = patch_id;
            m_h_dirichlet_point_info[dirichlet_point_id].local_id = local_id;
        }

        //m_dirichlet_context.init(m_h_dirichlet_info, num_dirichlet_boundary);
        m_dirichlet_point_context.init(m_h_dirichlet_point_info, num_dirichlet_point);
    }


    void GPSplineSolid::build_attr_mapping_context() {
        SplineAttributeGlobalMappingInfo *m_d_attr_mapping_info;
        CUDA_ERROR(cudaMalloc((void **) &m_d_attr_mapping_info,
                              m_num_geometry * sizeof(SplineAttributeGlobalMappingInfo)));

        for (uint32_t g = 0; g < m_num_geometry; ++g) {
            SplineAttributeGlobalMappingInfo d_mapping_info;
            d_mapping_info.geometry_id = g;

            CUDA_ERROR(cudaMalloc((void **) &d_mapping_info.ltog_mapping,
                                  64 * sizeof(uint32_t)));

            assert(m_h_geometry_ltog_cp[g].size() == 64);

            CUDA_ERROR(cudaMemcpy(d_mapping_info.ltog_mapping,
                                  m_h_geometry_ltog_cp[g].data(),
                                  64 * sizeof(uint32_t),
                                  cudaMemcpyHostToDevice));

            CUDA_ERROR(cudaMemcpy(m_d_attr_mapping_info + g, &d_mapping_info,
                                  sizeof(SplineAttributeGlobalMappingInfo),
                                  cudaMemcpyHostToDevice));
        }

        m_attr_mapping_context.init(m_d_attr_mapping_info);
    }

}