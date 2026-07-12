#pragma once

#include <assert.h>
#include "GP_Attribute.h"
#include "GP_Context.h"
#include "Kernels/GP_Query_Dispatcher.cuh"
#include "Util/GP_Vector.h"
#include <Eigen/Core>
#include "Kernels/GP_Debug.cuh"
#include "Kernels/GP_Debug_Eigen.cuh"


template<typename T>
__host__ __device__ void get_jacobi_u(Eigen::Vector<T, 192>& basis_value,
	const Eigen::Matrix<T, 64, 3>& displacement,
	T* coefs,
	Eigen::Matrix<T, 3, 3>& jacGeo,
	Eigen::Matrix<T, 3, 3>& jacDis) {
	Eigen::Map<Eigen::Matrix<T, 3, 64>> Nder(basis_value.data(), 3, 64);

	Eigen::Matrix<T, 64, 3> U, X;

	for (uint32_t i = 0; i < 64; ++i) {
		X(i, 0) = coefs[3 * i + 0];
		X(i, 1) = coefs[3 * i + 1];
		X(i, 2) = coefs[3 * i + 2];
	}

	U = displacement;

	jacGeo = Nder * X; //jacX
	jacDis = Nder * U; //jacU
}


template<typename T>
__host__ __device__ void transform_gradients(const Eigen::Matrix<T, 3, 3>& jac,
	const Eigen::Vector<T, 3 * 64>& grads,
	Eigen::Matrix<double, 3, 64>& trfGradsK,
	int dim = 3) {
	const Eigen::Map<const Eigen::Matrix<double, 3, 64>> grads_k(grads.data(), 3, 64);
	//declain no alias name confuse
	trfGradsK.noalias() = jac.inverse() * grads_k;
}

template<typename T>
__host__ __device__ void setBT(T B[][6], const Eigen::Matrix<T, 3, 3>& F, const Eigen::Vector<T, 3> bGrad) {
	short a_dim = F.cols();
	for (short j = 0; j < a_dim; ++j) {
		for (short i = 0; i < a_dim; ++i)
			B[j][i] = F(j, i) * bGrad(i);

		if (a_dim == 2)
			B[j][2] = F(j, 0) * bGrad(1) + F(j, 1) * bGrad(0);

		if (a_dim == 3)
			for (short i = 0; i < a_dim; ++i) {
				short k = (i + 1) % a_dim;
				B[j][i + a_dim] = F(j, i) * bGrad(k) + F(j, k) * bGrad(i);
			}
	}
}

template<typename T>
__host__ __device__ void setBT(T B[3][6], const T F[3][3], const T bGradX, const T bGradY, const T bGradZ) {
	short a_dim = 3;
	for (short j = 0; j < a_dim; ++j) {

		B[j][0] = F[j][0] * bGradX;
		B[j][1] = F[j][1] * bGradY;
		B[j][2] = F[j][2] * bGradZ;
		B[j][3] = F[j][0] * bGradY + F[j][1] * bGradX;
		B[j][4] = F[j][1] * bGradZ + F[j][2] * bGradY;
		B[j][5] = F[j][2] * bGradX + F[j][0] * bGradZ;
	}
}

template<typename T>
__host__ __device__ void setB(T B[][3], const Eigen::Matrix<T, 3, 3>& F, const Eigen::Vector<T, 3>& bGrad) {
	short a_dim = F.cols();

	for (short j = 0; j < a_dim; ++j) {
		for (short i = 0; i < a_dim; ++i)
			B[i][j] = F(j, i) * bGrad(i);

		if (a_dim == 2)
			B[2][j] = F(j, 0) * bGrad(1) + F(j, 1) * bGrad(0);

		if (a_dim == 3)
			for (short i = 0; i < a_dim; ++i) {
				short k = (i + 1) % a_dim;
				B[i + a_dim][j] = F(j, i) * bGrad(k) + F(j, k) * bGrad(i);
			}
	}
}

template<typename T>
__host__ __device__ void setB(T B[][3], const T F[3][3], const T bGradX, const T bGradY, const T bGradZ) {
	short a_dim = 3;

	for (short j = 0; j < a_dim; ++j) {

		B[0][j] = F[j][0] * bGradX;
		B[1][j] = F[j][1] * bGradY;
		B[2][j] = F[j][2] * bGradZ;
		B[3][j] = F[j][0] * bGradY + F[j][1] * bGradX;
		B[4][j] = F[j][1] * bGradZ + F[j][2] * bGradY;
		B[5][j] = F[j][2] * bGradX + F[j][0] * bGradZ;
	}
}

template<int m, int n, int k>
__host__ __device__ static void matMutly(double A[m][n], double B[n][k], double res[m][k]) {

	for (int i = 0; i < m; ++i) {
		for (int j = 0; j < k; ++j)
			res[i][j] = 0.0;
	}

	for (int aa = 0; aa < n; aa++) {
		for (int mm = 0; mm < m; mm++) {
			for (int nn = 0; nn < k; nn++) {
				res[mm][nn] += A[mm][aa] * B[aa][nn];
			}
		}
	}
}

template<int m>
__host__ __device__ static void setC(double C[m][m], const Eigen::Matrix<double, 3, 3>& R, double lambda, double mu) {
	double a = lambda, b = mu;
	short a_dim = R.cols();
	short dimTensor = m;
	// voigt indices
	Eigen::Matrix<int, 6, 2> v(dimTensor, 2);
	if (a_dim == 2)
		v << 0, 0,
		1, 1,
		0, 1;
	if (a_dim == 3)
		v << 0, 0,
		1, 1,
		2, 2,
		0, 1,
		1, 2,
		0, 2;
	// componentwise definition
	for (short i = 0; i < dimTensor; ++i)
		for (short j = 0; j < dimTensor; ++j)
			C[i][j] = a * R(v(i, 0), v(i, 1)) * R(v(j, 0), v(j, 1)) +
			b * (R(v(i, 0), v(j, 0)) * R(v(i, 1), v(j, 1)) +
				R(v(i, 0), v(j, 1)) * R(v(i, 1), v(j, 0)));
}

template<typename T>
__host__ __device__ void voigtStressToK(const Eigen::Matrix<T, 3, 3>& S,
	Eigen::Vector<double, 6>& Svec) {
	//matrix to tensor
	Svec(0) = S(0, 0);
	Svec(1) = S(1, 1);
	Svec(2) = S(2, 2);
	Svec(3) = S(0, 1);
	Svec(4) = S(1, 2);
	Svec(5) = S(0, 2);
}

inline __host__ __device__ void get_gs_3d_id(uint16_t gs_3d_id[3], uint16_t gs_point_id) {
	gs_3d_id[0] = gs_point_id / 16;
	gs_3d_id[1] = (gs_point_id - gs_3d_id[0] * 16) / 4;
	gs_3d_id[2] = gs_point_id - gs_3d_id[0] * 16 - gs_3d_id[1] * 4;

	assert(gs_point_id == 16 * gs_3d_id[0] + 4 * gs_3d_id[1] + gs_3d_id[2]);
}

///<<<dirichlet_num, surf_point_num=16>>>
__global__ static void IGA_add_dirichlet_to_local_mats(
	const gpmesh::SplineContext context,
	const gpmesh::AssemblyContext assembly_context,
	const gpmesh::DirichletContext dirichlet_context
) {
	using namespace gpmesh;

	uint32_t dirichlet_id = blockIdx.x;
	uint16_t surf_id = threadIdx.x;

	uint16_t i = surf_id / 4;
	uint16_t j = surf_id - 4 * i;

	if (dirichlet_id >= dirichlet_context.get_num_dirichlet_boundary()) {
		return;
	}

	DirichletInfo dirichlet_info = dirichlet_context.get_dirichlet_info(dirichlet_id);

	assert(dirichlet_info.dirichlet_boundary_id == dirichlet_id);

	uint32_t geometry_id = dirichlet_info.geometry_id;
	uint16_t dim = dirichlet_info.dim;
	double value = dirichlet_info.value;

	AssemblyInfo& assembly_info = assembly_context.get_assembly_info(geometry_id);

	uint16_t local_index = INVALID16;

	assert(dirichlet_info.side != BoundarySide::NONE);

	uint16_t local_index_v2 = INVALID16;
	uint16_t i_multiply = INVALID16;
	uint16_t j_multiply = INVALID16;
	uint16_t index_add = INVALID16;

	switch (dirichlet_info.side) {
	case BoundarySide::EAST:
		i_multiply = 16;
		j_multiply = 4;
		index_add = 3;
		break;
	case BoundarySide::WEST:
		i_multiply = 16;
		j_multiply = 4;
		index_add = 0;
		break;
	case BoundarySide::NORTH:
		i_multiply = 16;
		j_multiply = 1;
		index_add = 12;
		break;
	case BoundarySide::SOUTH:
		i_multiply = 16;
		j_multiply = 1;
		index_add = 0;
		break;
	case BoundarySide::BACK:
		i_multiply = 4;
		j_multiply = 1;
		index_add = 48;
		break;
	case BoundarySide::FRONT:
		i_multiply = 4;
		j_multiply = 1;
		index_add = 0;
		break;
	}

	local_index_v2 = i_multiply * i + j_multiply * j + index_add;

	if (dirichlet_info.side == BoundarySide::WEST) {
		local_index = 16 * i + 4 * j;
	}
	else if (dirichlet_info.side == BoundarySide::EAST) {
		local_index = 16 * i + 4 * j + 3;
	}
	else if (dirichlet_info.side == BoundarySide::SOUTH) {
		local_index = 16 * i + j;
	}
	else if (dirichlet_info.side == BoundarySide::NORTH) {
		local_index = 16 * i + j + 12;
	}
	else if (dirichlet_info.side == BoundarySide::FRONT) {
		local_index = 4 * i + j;
	}
	else if (dirichlet_info.side == BoundarySide::BACK) {
		local_index = 4 * i + j + 48;
	}

	assert(local_index != INVALID16);
	assert(local_index == local_index_v2);

	constexpr double DirchletMaxNum = 1e10;

	//assembly_info.local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)] += DirchletMaxNum;
	atomicAdd(&assembly_info.local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)], DirchletMaxNum);
}

//<<<DIVIDE_UP(num_dirichlet_point, 256),  256>>>
__global__ static void IGA_add_dirichlet_point_info_to_local_mats(
	const gpmesh::SplineContext context,
	const gpmesh::AssemblyContext assembly_context,
	const gpmesh::DirichletPointContext dirichlet_point_context,
	uint32_t num_dirichlet_point
) {
	using namespace gpmesh;

	uint32_t dirichlet_id = blockIdx.x * blockDim.x + threadIdx.x;

	if (dirichlet_id >= num_dirichlet_point) return;

	DirichletPointInfo dirichlet_point_info = dirichlet_point_context.get_dirichlet_point_info(dirichlet_id);


	uint32_t geometry_id = dirichlet_point_info.geometry_id;
	uint16_t local_id = dirichlet_point_info.local_id;

	AssemblyInfo& assembly_info = assembly_context.get_assembly_info(geometry_id);

	constexpr double DirchletMaxNum = 1e10;

	//assembly_info.local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)] += DirchletMaxNum;
	//atomicAdd(&assembly_info.local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)], DirchletMaxNum);
	atomicAdd(&assembly_info.local_mats[(local_id + 0 * 64) + 192 * (local_id + 0 * 64)], DirchletMaxNum);
	atomicAdd(&assembly_info.local_mats[(local_id + 1 * 64) + 192 * (local_id + 1 * 64)], DirchletMaxNum);
	atomicAdd(&assembly_info.local_mats[(local_id + 2 * 64) + 192 * (local_id + 2 * 64)], DirchletMaxNum);
}

template<typename T, size_t blockThreads>
__global__ static void
IGA_debug_kernel(const gpmesh::SplineContext context, const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	//load all information in shrd_mem
	uint32_t geometry_id = blockIdx.x;
	uint32_t gs_point_id = threadIdx.x;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	//    if (blockIdx.x == 0 && threadIdx.x < integral_context.get_gaussian_num_u())
	//    {
	//        printf("point-%f with index-%d has weight-%f\n", integral_context.get_points_u()[threadIdx.x], threadIdx.x, integral_context.get_weights_u()[threadIdx.x]);
	//    }

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();


	uint32_t q = gs_point_id;
	uint32_t kk = q / 16;
	uint32_t jj = (q - kk * 16) / 4;
	uint32_t ii = q - kk * 16 - jj * 4;
	assert(q == kk * 16 + jj * 4 + ii);

	//calculate Displacement Jacobi, and Geometry Jacobi

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);

	//    if (blockIdx.x == 0 && threadIdx.x == 0)
	//    {
	//        for (int i = 0; i < 64; ++i)
	//        {
	//            for (int j = 0; j < 3; ++j)
	//            {
	//                printf("row[%d],col[%d]=%f\n", i, j, disp_block_geo->operator()(i, j));
	//            }
	//        }
	//    }
	//
		//Eigen
		//Eigen::Matrix<double, 3 * 64, 3 * 64> localMats;

	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();
	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(q);
	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;
	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[ii] * integral_context.get_weights_v()[jj] *
		integral_context.get_weights_w()[kk] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;
	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	double m_Rx, m_Ry, m_Rz, m_Qx, m_Qy, m_Qz;

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	if (geometry_id == 0 && (ii + jj + kk) == 0) {
		//printf("start: localMas(%d,%d)=%f\n", 0, 0, localMats[0 + 192 * 0]);
	}

	for (int j = 0; j < 64; ++j) {
		m_Rx = physGrad.col(j)(0);
		m_Ry = physGrad.col(j)(1);
		m_Rz = physGrad.col(j)(2);

		for (int k = 0; k < 64; ++k) {
			m_Qx = physGrad.col(k)(0);
			m_Qy = physGrad.col(k)(1);
			m_Qz = physGrad.col(k)(2);

			uint32_t N_D = 64;

			//            if (blockIdx.x == 0 && j == 0 && k == 0 && (ii + jj + kk) == 0)
			//            {
			//                printf("Rx=%.15f, Ry=%.15f, Rz=%.15f \n Qx=%.15f, Qy=%.15f, Qz=%.15f \n weight=%.20f \n",
			//                       m_Rx, m_Ry, m_Rz,
			//                       m_Qx, m_Qy, m_Qz,
			//                       weight);
			//            }

			//            //test C
			//            if (blockIdx.x == 0 && j == 0 && k == 0 && (ii + jj + kk) == 0)
			//            {
			//                printf("Matrix C\n");
			//
			//                for (int i = 0; i < 6; ++i)
			//                {
			//                    for (int j = 0; j < 6; ++j)
			//                    {
			//                        printf("C(%d,%d)=%.5f ", i, j, C(i, j));
			//                    }
			//
			//                    printf("\n");
			//                }
			//            }

			//            double test = weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
			//
			//            if (blockIdx.x == 0 && (ii + jj + kk) == 0 && j == 0 && k == 0)
			//            {
			//                printf("localMas(%d,%d)=%f\n", j, k, localMats[j + 192 * k]);
			//                printf("test=%f\n", test);
			//            }

			atomicAdd(&localMats[(j + 0 * N_D) + 192 * (k + 0 * N_D)],
				weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz)));
			atomicAdd(&localMats[(j + 0 * N_D) + 192 * (k + 1 * N_D)],
				weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx));
			atomicAdd(&localMats[(j + 0 * N_D) + 192 * (k + 2 * N_D)],
				weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qx));

			atomicAdd(&localMats[(j + 1 * N_D) + 192 * (k + 0 * N_D)],
				weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy));
			atomicAdd(&localMats[(j + 1 * N_D) + 192 * (k + 1 * N_D)],
				weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz)));
			atomicAdd(&localMats[(j + 1 * N_D) + 192 * (k + 2 * N_D)],
				weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy));

			atomicAdd(&localMats[(j + 2 * N_D) + 192 * (k + 0 * N_D)],
				weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz));
			atomicAdd(&localMats[(j + 2 * N_D) + 192 * (k + 1 * N_D)],
				weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz));
			atomicAdd(&localMats[(j + 2 * N_D) + 192 * (k + 2 * N_D)],
				weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx)));

			//            localMats[192 * j + k] += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
			//            localMats[192 * j + k + N_D] += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
			//            localMats[192 * j + k + 2 * N_D] += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
			//            localMats[192 * (j + N_D) + k] += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
			//            localMats[192 * (j + N_D) + k + N_D] += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
			//            localMats[192 * (j + N_D) + k + 2 * N_D] += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
			//            localMats[192 * (j + 2 * N_D) + k] += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
			//            localMats[192 * (j + 2 * N_D) + k + N_D] += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
			//            localMats[192 * (j + 2 * N_D) + k + 2 * N_D] += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));

			//            localMats[j][k] += weight * (C(0, 0) * m_Rx * m_Qx + C(3, 3) * (m_Ry * m_Qy + m_Rz * m_Qz));
			//            localMats[j][k + N_D] += weight * (C(0, 1) * m_Rx * m_Qy + C(3, 3) * m_Ry * m_Qx);
			//            localMats[j][k + 2 * N_D] += weight * (C(0, 1) * m_Rx * m_Qz + C(3, 3) * m_Rz * m_Qy);
			//            localMats[j + N_D][ k] += weight * (C(0, 1) * m_Ry * m_Qx + C(3, 3) * m_Rx * m_Qy);
			//            localMats[j + N_D][k + N_D] += weight * (C(0, 0) * m_Ry * m_Qy + C(3, 3) * (m_Rx * m_Qx + m_Rz * m_Qz));
			//            localMats[j + N_D][k + 2 * N_D] += weight * (C(0, 1) * m_Ry * m_Qz + C(3, 3) * m_Rz * m_Qy);
			//            localMats[j + 2 * N_D][ k] += weight * (C(0, 1) * m_Rz * m_Qx + C(3, 3) * m_Rx * m_Qz);
			//            localMats[j + 2 * N_D][k + N_D] += weight * (C(0, 1) * m_Rz * m_Qy + C(3, 3) * m_Ry * m_Qz);
			//            localMats[j + 2 * N_D][k + 2 * N_D] += weight * (C(0, 0) * m_Rz * m_Qz + C(3, 3) * (m_Ry * m_Qy + m_Rx * m_Qx));
		}
	}

}

template<typename T, uint16_t blockThreads>
__global__ static void IGA_debug_matvec_A_Multiply(const gpmesh::SplineContext context,
	const gpmesh::AssemblyContext assembly_context,
	const gpmesh::SplineAttributeBase<T> in,
	const gpmesh::SplineAttributeBase<T> out
) {
	using namespace gpmesh;

	//load all information in shrd_mem
	uint32_t geometry_id = blockIdx.x;
	uint16_t local_cpoint_id = threadIdx.x;

	if (geometry_id >= context.get_num_geometry() || local_cpoint_id >= blockThreads) {
		return;
	}

	AssemblyInfo assembly_info = assembly_context.get_assembly_info()[geometry_id];
	double* localMats = assembly_info.local_mats;

	assert(assembly_info.geometry_id == geometry_id);

	for (uint32_t d = 0; d < 3; ++d) {
		double sum = 0.0;

		for (uint16_t ci = 0; ci < 64; ++ci) {
			for (int cd = 0; cd < 3; ++cd) {
				sum += localMats[(local_cpoint_id + d * blockThreads) + 192 * (ci + cd * blockThreads)] *
					in(geometry_id, ci, cd);
			}
		}

		out(geometry_id, local_cpoint_id, d) = sum;
	}
}

//<<<(geometry_num), (gaussian_num)>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assemble_debug_kernel(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	//load all information in shrd_mem
	//(num_geometry, 1, 1)
	uint32_t geometry_id = blockIdx.x;
	//(64, 64, 64)
	uint32_t gs_point_id = threadIdx.x;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint32_t q = gs_point_id;
	uint32_t kk = q / 16;
	uint32_t jj = (q - kk * 16) / 4;
	uint32_t ii = q - kk * 16 - jj * 4;
	assert(q == kk * 16 + jj * 4 + ii);

	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));    //��÷����
	mu = e / (2. * (1. + pr));

	//calculate Displacement Jacobi, and Geometry Jacobi

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);

	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();
	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(q);
	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[ii] * integral_context.get_weights_v()[jj] *
		integral_context.get_weights_w()[kk] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	double J = F.determinant();

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	Eigen::Vector<double, 6> Svec;

	voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z
	uint32_t j = threadIdx.y;
	uint32_t k = threadIdx.z;

	//�ֲ�����װ��
	double B_i[3][6], B_j[6][3];

	setBT<double>(B_i, F, physGrad.col(j));

	double materialTangentTemp1[3][6];

	double CD[6][6];

	double materialTangent1[3][3];

	setC<6>(CD, identity, lambda, mu);

	matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);

	setB<double>(B_j, F, physGrad.col(k));

	matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

	Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

	double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

	for (uint16_t d = 0; d < 3; ++d) {
		materialTangent1[d][d] += geometricTangent;
	}

	constexpr int N_D = 64;
	for (uint16_t di = 0; di < 3; ++di) {
		for (uint16_t dj = 0; dj < 3; ++dj) {
			atomicAdd(&localMats[(j + di * N_D) + 192 * (k + dj * N_D)], weight * materialTangentTemp1[di][dj]);
		}
	}
}

//dim3 numBlocks(num_geometry, 64, 1) numBlocks = blocksInGrid, Candidate choice = (num_geometry, 8, 8) but depends on memory hit situation
//v1: dim3 threadsPerBlock(64, 1, 1)  (max threads 1024, 每个thread 负责64 iters)
//v2: dim3 threadsPerBlock(32, 32, 1) (max threads 1024, 每个thread 负责2x2 iters)

//v1版本
//<<<dim3(11, 64, 1), dim3(64, 1, 1)>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assembly_debug_kernel_v1(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint32_t q = gs_point_id;
	uint32_t kk = q / 16;
	uint32_t jj = (q - kk * 16) / 4;
	uint32_t ii = q - kk * 16 - jj * 4;
	assert(q == kk * 16 + jj * 4 + ii);


	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));    //��÷����
	mu = e / (2. * (1. + pr));

	//calculate Displacement Jacobi, and Geometry Jacobi

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	//可以修改成并行化读取
	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();
	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(q);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	//可以修改成并行化计算
	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[ii] * integral_context.get_weights_v()[jj] *
		integral_context.get_weights_w()[kk] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;

	//可以改成并行化计算
	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	double J = F.determinant();

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	Eigen::Vector<double, 6> Svec;

	voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	uint16_t j = threadIdx.x;

	setBT<double>(B_i, F, physGrad.col(j));

	double materialTangentTemp1[3][6];

	double CD[6][6];

	double materialTangent1[3][3];

	setC<6>(CD, identity, lambda, mu);

	matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);

	constexpr int N_D = 64;

	for (uint16_t k = 0; k < N_D; ++k) {
		setB<double>(B_j, F, physGrad.col(k));

		matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

		Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

		double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

		for (uint16_t d = 0; d < 3; ++d) {
			materialTangent1[d][d] += geometricTangent;
		}

		for (uint16_t di = 0; di < 3; ++di) {
			for (uint16_t dj = 0; dj < 3; ++dj) {
				atomicAdd(&localMats[(j + di * N_D) + 192 * (k + dj * N_D)], weight * materialTangentTemp1[di][dj]);
			}
		}
	}
}

//v2版本
//<<<dim3(m_num_geometry, m_num_gaussian_point=64, 1), dim3(16, 16, 1)>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assembly_debug_kernel_v2(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//for debug
	//print_matrix3d(jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	double J = F.determinant();

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	Eigen::Vector<double, 6> Svec;

	voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}


//v2 debug版本
//<<<dim3(m_num_geometry, m_num_gaussian_point=64, 1), dim3(16, 16, 1)>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assembly_debug_kernel_v2_debug(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;


	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//print_double(5, weight);
	//Debug jacGeo
	//print_matrix3d(5, jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	//Debug Show physGrad
	//print_RowxCol_matrix3d<uint16_t, 3, 64>(5, physGrad);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	//Debug phyDispJac
	//print_matrix3d(5, phyDispJac);

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	//double J = F.determinant();

	//Debug
	//print_matrix3d(5, F);

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	//Debug
	//print_matrix3d(5, S);

	//Eigen::Vector<double, 6> Svec;

	//voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}

template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assembly_debug_kernel_v3(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context,
	double lambda, double mu) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
//    double lambda = 576.92307692307691;
//    double mu = 384.61538461538458;
//    double e = 1000.0000000000000;
//    double pr = 0.3;
//
//
//    lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
//    mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//print_double(5, weight);
	//Debug jacGeo
	//print_matrix3d(5, jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	//Debug Show physGrad
	//print_RowxCol_matrix3d<uint16_t, 3, 64>(5, physGrad);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	//Debug phyDispJac
	//print_matrix3d(5, phyDispJac);

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	//double J = F.determinant();

	//Debug
	//print_matrix3d(5, F);

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	//Debug
	//print_matrix3d(5, S);

	//Eigen::Vector<double, 6> Svec;

	//voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}

template<typename T, size_t blockThreads>
__global__ static void Assembly_Local_Matrix_With_Yangs(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context,
	double* d_Yangs,
	double* d_Poisson) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	e = *d_Yangs;
	pr = *d_Poisson;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//print_double(5, weight);
	//Debug jacGeo
	//print_matrix3d(5, jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	//Debug Show physGrad
	//print_RowxCol_matrix3d<uint16_t, 3, 64>(5, physGrad);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	//Debug phyDispJac
	//print_matrix3d(5, phyDispJac);

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	//double J = F.determinant();

	//Debug
	//print_matrix3d(5, F);

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	//Debug
	//print_matrix3d(5, S);

	//Eigen::Vector<double, 6> Svec;

	//voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}


template<typename T, size_t blockThreads>
__global__ static void TopOpt_Assembly_Local_Matrix(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context,
	double* top_density) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;
	//printf("density = %f", top_density[geometry_id]);

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//print_double(5, weight);
	//Debug jacGeo
	//print_matrix3d(5, jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	//Debug Show physGrad
	//print_RowxCol_matrix3d<uint16_t, 3, 64>(5, physGrad);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	//Debug phyDispJac
	//print_matrix3d(5, phyDispJac);

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	//double J = F.determinant();

	//Debug
	//print_matrix3d(5, F);

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	//Debug
	//print_matrix3d(5, S);

	//Eigen::Vector<double, 6> Svec;

	//voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

			double top_extra_mul = top_density[geometry_id] * top_density[geometry_id] * top_density[geometry_id];
			//double top_extra_mul = top_density[geometry_id];
			//std::cout <<"density = " << top_density[geometry_id] << "\n";
#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)],
						weight * materialTangent1[di][dj] * top_extra_mul);
				}
			}
		}
	}
}




template<typename T, size_t blockThreads>
__global__ static void TopOpt_Assembly_Local_Matrix_MMC(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context,
	double* top_density) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;
	//printf("density = %f", top_density[geometry_id]);

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	//print_double(5, weight);
	//Debug jacGeo
	//print_matrix3d(5, jacGeo);

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	//Debug Show physGrad
	//print_RowxCol_matrix3d<uint16_t, 3, 64>(5, physGrad);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	//Debug phyDispJac
	//print_matrix3d(5, phyDispJac);

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	//double J = F.determinant();

	//Debug
	//print_matrix3d(5, F);

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	//Debug
	//print_matrix3d(5, S);

	//Eigen::Vector<double, 6> Svec;

	//voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

			//double top_extra_mul = top_density[geometry_id] * top_density[geometry_id] * top_density[geometry_id];
			double top_extra_mul = top_density[geometry_id];
			//std::cout <<"density = " << top_density[geometry_id] << "\n";
#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)],
						weight * materialTangent1[di][dj] * top_extra_mul);
				}
			}
		}
	}
}

template<typename T, size_t blockThreads>
__global__ static void TopOpt_Assembly_Local_Matrix_CellScale(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context,
	double* cell_scale) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry() || blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;
	uint16_t gs_point_id = blockIdx.y;
	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];
	extern __shared__ double shrd_coefs[];

	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];
	get_gs_3d_id(gs_kkjjii, gs_point_id);

	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0;
	double pr = 0.3;
	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();
	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();
	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;
	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;
	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());
	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;
	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;
	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();
	double B_i[3][6], B_j[6][3];
	constexpr uint16_t N_D = 4;
	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;
	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;
		setBT<double>(B_i, F, physGrad.col(j));

		double materialTangentTemp1[3][6];
		double CD[6][6];
		double materialTangent1[3][3];
		setC<6>(CD, identity, lambda, mu);
		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);

		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;
			setB<double>(B_j, F, physGrad.col(k));
			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);
			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

			double scale = cell_scale[geometry_id];
#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)],
						weight * materialTangent1[di][dj] * scale);
				}
			}
		}
	}
}


//v2版本
//<<<dim3(m_num_geometry, m_num_gaussian_point=64, 1), dim3(16, 16, 1)>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_hyperelasticity_assembly_debug_kernel_v3(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	double J = F.determinant();

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	Eigen::Vector<double, 6> Svec;

	voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}

//<<<(num_geometry,num_gauss_point), (16, 16)>>>
template<typename T>
__device__ void load_basis_value1_to_shared_memory(const Eigen::Matrix<T, 192, 64>* basis_value1, T gradients[3][64]) {
	uint16_t gauss_id = blockIdx.y;
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	if (idx < 192) {
		uint16_t r = idx % 3;
		uint16_t c = idx / 3;

		gradients[r][c] = basis_value1->operator()(idx, gauss_id);
	}

}

//<<<(num_geometry,num_gauss_point), (16, 16)>>>
template<typename T>
__device__ void load_physical_cpoints_position_to_shared_memory(const gpmesh::SplineContext* ptr_spline_context,
	T coefficients[3][64]) {
	uint32_t geometry_id = blockIdx.x;

	uint16_t gauss_id = blockIdx.y;
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	if (idx < 192) {
		uint16_t r = idx % 3;
		uint16_t c = idx / 3;

		coefficients[r][c] = ptr_spline_context->get_geometry_info()[geometry_id].coefs[idx];
	}
}

//<<<(num_geometry,num_gauss_point), (16, 16)>>>
template<typename T>
__device__ void
load_cpoints_displacement_to_shared_memory(const gpmesh::PhysicContext* ptr_physic_context, T coefficients[3][64]) {
	uint32_t geometry_id = blockIdx.x;
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	Eigen::Matrix<T, 64, 3>* local_cpoint_disp_mat = ptr_physic_context->get_cpoints_displacement(geometry_id);

	if (idx < 192) {
		uint16_t r = idx % 3;
		uint16_t c = idx / 3;

		coefficients[r][c] = local_cpoint_disp_mat->operator()(c, r);
	}
}

template<typename T>
__device__ void set_jacobian_matrix(T jacobian[3][3], T value) {
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	if (idx < 9) {
		uint16_t r = idx % 3;
		uint16_t c = idx / 3;

		jacobian[r][c] = value;
	}
}

//<<<(num_geometry, 64, 1), (16, 16, 1)>>>
template<typename T>
__device__ void calculate_jacobian_matrix(T gradients[3][64], T coefficients[3][64], T jacobian[3][3]) {
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	uint16_t working_id = idx % 64;

	uint16_t group_id[2];
	group_id[0] = idx / 64;
	group_id[0] *= 2;
	group_id[1] = group_id[0] + 1;


	//Handling [0-7]
#pragma unroll
	for (int i = 0; i < 2; ++i) {
		uint16_t r = group_id[i] / 3;
		uint16_t c = group_id[i] % 3;

		atomicAdd(&jacobian[r][c], gradients[r][working_id] * coefficients[c][working_id]);
	}

	//Handling [8]
	if (idx < 64) {
		atomicAdd(&jacobian[2][2], gradients[2][idx] * coefficients[2][idx]);
	}
}

template<typename T>
__device__ void caluclate_jac_inverse(T jac[3][3], T jac_inv[3][3]) {
	double det = jac[0][0] * (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) -
		jac[0][1] * (jac[1][0] * jac[2][2] - jac[1][2] * jac[2][0]) +
		jac[0][2] * (jac[1][0] * jac[2][1] - jac[1][1] * jac[2][0]);

	double invdet = 1 / det;

	jac_inv[0][0] = (jac[1][1] * jac[2][2] - jac[2][1] * jac[1][2]) * invdet;
	jac_inv[0][1] = (jac[0][2] * jac[2][1] - jac[0][1] * jac[2][2]) * invdet;
	jac_inv[0][2] = (jac[0][1] * jac[1][2] - jac[0][2] * jac[1][1]) * invdet;
	jac_inv[1][0] = (jac[1][2] * jac[2][0] - jac[1][0] * jac[2][2]) * invdet;
	jac_inv[1][1] = (jac[0][0] * jac[2][2] - jac[0][2] * jac[2][0]) * invdet;
	jac_inv[1][2] = (jac[1][0] * jac[0][2] - jac[0][0] * jac[1][2]) * invdet;
	jac_inv[2][0] = (jac[1][0] * jac[2][1] - jac[2][0] * jac[1][1]) * invdet;
	jac_inv[2][1] = (jac[2][0] * jac[0][1] - jac[0][0] * jac[2][1]) * invdet;
	jac_inv[2][2] = (jac[0][0] * jac[1][1] - jac[1][0] * jac[0][1]) * invdet;
}

template<typename T>
//phyDispjac = transpose(jacDis) * transpose(jacGeo_inv)
//正确性需要验证
__host__ __device__ void calculate_phyDispJac(T jacDis[3][3], T jacGeo_inv[3][3], T phyDispJac[3][3]) {
	phyDispJac[0][0] =
		jacGeo_inv[0][0] * jacDis[0][0] + jacGeo_inv[0][1] * jacDis[1][0] + jacGeo_inv[0][2] * jacDis[2][0];
	phyDispJac[0][1] =
		jacGeo_inv[1][0] * jacDis[0][0] + jacGeo_inv[1][1] * jacDis[1][0] + jacGeo_inv[1][2] * jacDis[2][0];
	phyDispJac[0][2] =
		jacGeo_inv[2][0] * jacDis[0][0] + jacGeo_inv[2][1] * jacDis[1][0] + jacGeo_inv[2][2] * jacDis[2][0];

	phyDispJac[1][0] =
		jacGeo_inv[0][0] * jacDis[0][1] + jacGeo_inv[0][1] * jacDis[1][1] + jacGeo_inv[0][2] * jacDis[2][1];
	phyDispJac[1][1] =
		jacGeo_inv[1][0] * jacDis[0][1] + jacGeo_inv[1][1] * jacDis[1][1] + jacGeo_inv[1][2] * jacDis[2][1];
	phyDispJac[1][2] =
		jacGeo_inv[2][0] * jacDis[0][1] + jacGeo_inv[2][1] * jacDis[1][1] + jacGeo_inv[2][2] * jacDis[2][1];

	phyDispJac[2][0] =
		jacGeo_inv[0][0] * jacDis[0][2] + jacGeo_inv[0][1] * jacDis[1][2] + jacGeo_inv[0][2] * jacDis[2][2];
	phyDispJac[2][1] =
		jacGeo_inv[1][0] * jacDis[0][2] + jacGeo_inv[1][1] * jacDis[1][2] + jacGeo_inv[1][2] * jacDis[2][2];
	phyDispJac[2][2] =
		jacGeo_inv[2][0] * jacDis[0][2] + jacGeo_inv[2][1] * jacDis[1][2] + jacGeo_inv[2][2] * jacDis[2][2];
}

template<typename T>
__device__ void calculate_physic_gradients(T jacGeo[3][3], T gradients[3][64], T physic_gradients[3][64]) {
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	uint16_t col_id = idx % 64;
	//[0, 3]
	uint16_t row_id = idx / 64;

	//Todo: 冗余计算，是否可以尝试并行化
	//calculate inverse of jacGeo
	T jacGeo_inc[3][3];

	// determin of matrix

	double det = jacGeo[0][0] * (jacGeo[1][1] * jacGeo[2][2] - jacGeo[2][1] * jacGeo[1][2]) -
		jacGeo[0][1] * (jacGeo[1][0] * jacGeo[2][2] - jacGeo[1][2] * jacGeo[2][0]) +
		jacGeo[0][2] * (jacGeo[1][0] * jacGeo[2][1] - jacGeo[1][1] * jacGeo[2][0]);

	double invdet = 1 / det;

	jacGeo_inc[0][0] = (jacGeo[1][1] * jacGeo[2][2] - jacGeo[2][1] * jacGeo[1][2]) * invdet;
	jacGeo_inc[0][1] = (jacGeo[0][2] * jacGeo[2][1] - jacGeo[0][1] * jacGeo[2][2]) * invdet;
	jacGeo_inc[0][2] = (jacGeo[0][1] * jacGeo[1][2] - jacGeo[0][2] * jacGeo[1][1]) * invdet;
	jacGeo_inc[1][0] = (jacGeo[1][2] * jacGeo[2][0] - jacGeo[1][0] * jacGeo[2][2]) * invdet;
	jacGeo_inc[1][1] = (jacGeo[0][0] * jacGeo[2][2] - jacGeo[0][2] * jacGeo[2][0]) * invdet;
	jacGeo_inc[1][2] = (jacGeo[1][0] * jacGeo[0][2] - jacGeo[0][0] * jacGeo[1][2]) * invdet;
	jacGeo_inc[2][0] = (jacGeo[1][0] * jacGeo[2][1] - jacGeo[2][0] * jacGeo[1][1]) * invdet;
	jacGeo_inc[2][1] = (jacGeo[2][0] * jacGeo[0][1] - jacGeo[0][0] * jacGeo[2][1]) * invdet;
	jacGeo_inc[2][2] = (jacGeo[0][0] * jacGeo[1][1] - jacGeo[1][0] * jacGeo[0][1]) * invdet;

	if (row_id < 3) {
		physic_gradients[row_id][col_id] = jacGeo_inc[row_id][0] * gradients[0][col_id] +
			jacGeo_inc[row_id][1] * gradients[1][col_id] +
			jacGeo_inc[row_id][2] * gradients[2][col_id];
	}
}

template<typename T>
__device__ void
calculate_physic_gradients_optimized(T jacGeo_inv[3][3], T gradients[3][64], T physic_gradients[3][64]) {
	uint16_t idx = threadIdx.x + blockDim.x * threadIdx.y;

	uint16_t col_id = idx % 64;
	//[0, 3]
	uint16_t row_id = idx / 64;

	if (row_id < 3) {
		physic_gradients[row_id][col_id] = jacGeo_inv[row_id][0] * gradients[0][col_id] +
			jacGeo_inv[row_id][1] * gradients[1][col_id] +
			jacGeo_inv[row_id][2] * gradients[2][col_id];
	}
}

template<typename T>
__device__ static void calculate_F_and_S(T phyDispJac[3][3], T lambda, T mu, T F[3][3], T S[3][3]) {

#pragma unroll
	for (int i = 0; i < 3; ++i) {
#pragma unroll
		for (int j = 0; j < 3; ++j) {
			F[i][j] = phyDispJac[i][j];
		}
	}

	F[0][0] += 1.0;
	F[1][1] += 1.0;
	F[2][2] += 1.0;

	//Calculate RCG * 0.5
	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			S[i][j] = (F[0][i] * F[0][j] + F[1][i] * F[1][j] + F[2][i] * F[2][j]) * 0.5;
		}
	}

	T trace = 0;
	//Calculate E and it`s trace
	S[0][0] -= 0.5;
	trace += S[0][0];
	S[1][1] -= 0.5;
	trace += S[1][1];
	S[2][2] -= 0.5;
	trace += S[2][2];

	for (int i = 0; i < 3; ++i) {
		for (int j = 0; j < 3; ++j) {
			S[i][j] *= 2 * mu;
		}
	}

	S[0][0] += lambda * trace;
	S[1][1] += lambda * trace;
	S[2][2] += lambda * trace;
}

template<typename T>
__device__ static T calculate_determinant(T matrix[3][3]) {
	T x = ((matrix[1][1] * matrix[2][2]) - (matrix[2][1] * matrix[1][2]));
	T y = ((matrix[1][0] * matrix[2][2]) - (matrix[2][0] * matrix[1][2]));
	T z = ((matrix[1][0] * matrix[2][1]) - (matrix[2][0] * matrix[1][1]));

	return ((matrix[0][0] * x) - (matrix[0][1] * y) + (matrix[0][2] * z));
}

template<typename T>
//<<<(num_geometry, num_gaussian=64, 1), (16, 16, 1)>>>
__global__ static void IGA_he_assembly_kernel_optimized_v0(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	__shared__ T gradients[3][64];
	__shared__ T coefficients[3][64];

	//load basisValue1 to gradients
	load_basis_value1_to_shared_memory(integral_context.get_basis_value_1_matrix(), gradients);

	//load cpoint physical position to coefficients
	load_physical_cpoints_position_to_shared_memory(&context, coefficients);

	__syncthreads();

	__shared__ T jacGeo[3][3];
	__shared__ T jacDisp[3][3];

	//Empty Matrix
	set_jacobian_matrix<double>(jacGeo, 0);

	__syncthreads();

	calculate_jacobian_matrix<double>(gradients, coefficients, jacGeo);

	__syncthreads();

	//Debug jacGeo
	//print_3x3_double_array(5, jacGeo);

	//load cpoint displacement to coefficients;
	//如果是第一次运行，控制点并没有发生位移
	load_cpoints_displacement_to_shared_memory<double>(&physic_context, coefficients);

	__syncthreads();

	//Empty Matrix
	set_jacobian_matrix<double>(jacDisp, 0);

	__syncthreads();

	calculate_jacobian_matrix<double>(gradients, coefficients, jacDisp);

	T jacGeo_inv[3][3];

	caluclate_jac_inverse(jacGeo, jacGeo_inv);

	__syncthreads();

	__shared__ T physGrad[3][64];

	calculate_physic_gradients_optimized(jacGeo_inv, gradients, physGrad);

	T phyDispJac[3][3];

	calculate_phyDispJac(jacDisp, jacGeo_inv, phyDispJac);

	__syncthreads();

	//Debug PhyDispJac
	//print_3x3_double_array(5, phyDispJac);

	//Debug Show physGrad
	//print_RowxCol_double_array<uint16_t, 3, 64>(5, physGrad);

	uint16_t gs_kkjjii[3];

	//Calculate Gaussian Point Index
	get_gs_3d_id(gs_kkjjii, gs_point_id);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(calculate_determinant<T>(jacGeo));

	//print_double(5, weight);

	//physics properties sent to PhysicContext
	const T e = 1000.0000000000000;
	const T pr = 0.3;

	const T lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	const T mu = e / (2. * (1. + pr));

	T* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	Eigen::Matrix<T, 6, 6>& C = *physic_context.get_C();

	T F[3][3];
	T S[3][3];

	calculate_F_and_S<T>(phyDispJac, lambda, mu, F, S);

	//print_3x3_double_array(5, F);
	//print_3x3_double_array(5, S);

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad[0][j], physGrad[1][j], physGrad[2][j]);

		//print_RowxCol_double_array<uint16_t, 3, 6>(5, B_i);

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		Eigen::Matrix<T, 3, 3> identity = Eigen::Matrix<T, 3, 3>::Identity();

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);

		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad[0][k], physGrad[1][k], physGrad[2][k]);

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			T geometricTangentTemp[3];

			geometricTangentTemp[0] = S[0][0] * physGrad[0][j] + S[0][1] * physGrad[1][j] + S[0][2] * physGrad[2][j];
			geometricTangentTemp[1] = S[1][0] * physGrad[0][j] + S[1][1] * physGrad[1][j] + S[1][2] * physGrad[2][j];
			geometricTangentTemp[2] = S[2][0] * physGrad[0][j] + S[2][1] * physGrad[1][j] + S[2][2] * physGrad[2][j];

			double geometricTangent = geometricTangentTemp[0] * physGrad[0][k]
				+ geometricTangentTemp[1] * physGrad[1][k]
				+ geometricTangentTemp[2] * physGrad[2][k];

			if (local_j == 0 && local_k == 0) {
				//print_RowxCol_double_array<uint16_t, 6, 3>(5, B_j);
				//print_double(5, geometricTangent);
				//print_3x3_double_array(5, materialTangent1);
			}

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}
	}
}

//<<<dim3(m_num_geometry, 1, 1), dim3(16, 16, 1)>>>
__global__ static void IGA_clear_assembly_context(const gpmesh::AssemblyContext assembly_context) {

	uint32_t geometry_id = blockIdx.x;

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

#pragma unroll
			for (uint16_t di = 0; di < 3; ++di) {
#pragma unroll
				for (uint16_t dj = 0; dj < 3; ++dj) {
					localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)] = 0;
				}
			}
		}
	}
}


template<typename T, size_t blockThreads>
__global__ static void
IGA_hyperelasticity_assembly_debug_kernel_v2_add_residual_to_rhs(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::AssemblyContext assembly_context) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 192> basis_value = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant());

	Eigen::Matrix<double, 3, 64> physGrad;

	transform_gradients<double>(jacGeo, basis_value, physGrad, 3);

	Eigen::Matrix<double, 3, 3> phyDispJac = jacDis.transpose() * (jacGeo.transpose().inverse());

	Eigen::Matrix<double, 3, 3> F = identity + phyDispJac;
	double J = F.determinant();

	Eigen::Matrix<double, 3, 3> RCG = F.transpose() * F;

	Eigen::Matrix<double, 3, 3> E = 0.5 * (RCG - identity);

	//using stvk material
	Eigen::Matrix<double, 3, 3> S = lambda * E.trace() * identity + 2 * mu * E;

	Eigen::Vector<double, 6> Svec;

	voigtStressToK(S, Svec);

	double* localMats = assembly_context.get_assembly_info()[geometry_id].local_mats;
	double* localRhs = assembly_context.get_assembly_info()[geometry_id].local_rhs;

	Eigen::Matrix<double, 6, 6>& C = *physic_context.get_C();

	//replace j,k iteration with threadIdx.y, threadIdx.z)

	double B_i[3][6], B_j[6][3];

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		setBT<double>(B_i, F, physGrad.col(j));

		double materialTangentTemp1[3][6];

		double CD[6][6];

		double materialTangent1[3][3];

		setC<6>(CD, identity, lambda, mu);

		matMutly<3, 6, 6>(B_i, CD, materialTangentTemp1);


		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			setB<double>(B_j, F, physGrad.col(k));

			matMutly<3, 6, 3>(materialTangentTemp1, B_j, materialTangent1);

			Eigen::Matrix<double, 3, 1> geometricTangentTemp = S * physGrad.col(j);

			double geometricTangent = geometricTangentTemp.transpose() * physGrad.col(k);

			for (uint16_t d = 0; d < 3; ++d) {
				materialTangent1[d][d] += geometricTangent;
			}

			for (uint16_t di = 0; di < 3; ++di) {
				for (uint16_t dj = 0; dj < 3; ++dj) {
					atomicAdd(&localMats[(j + di * GS_D) + 192 * (k + dj * GS_D)], weight * materialTangent1[di][dj]);
				}
			}
		}

		double localResiduals[3] = { 0 };
		for (int d = 0; d < 6; ++d) {
			localResiduals[0] += B_i[0][d] * Svec(d);
			localResiduals[1] += B_i[1][d] * Svec(d);
			localResiduals[2] += B_i[2][d] * Svec(d);
		}

		for (uint16_t d = 0; d < 3; ++d) {
			atomicAdd(&localRhs[d * 64 + j], -weight * localResiduals[d]);
		}
		//此外，需要对每个localRhs施加外力项 需要执行一个额外的kernel函数
	}
}

//rhs装配


//mass装配
//v2版本
//<<<dim3(11, 64, 1), dim3(16, 16, 1)>>>
//非常差的并行化版本
template<typename T, size_t blockThreads>
__global__ static void IGA_mass_assembly_debug_kernel_v0(const gpmesh::SplineContext context,
	const gpmesh::GaussianIntegralContext integral_context,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::MassContext mass_context,
	double density
) {
	using namespace gpmesh;

	if (blockIdx.x >= context.get_num_geometry()) {
		return;
	}

	if (blockIdx.y >= 64) {
		return;
	}

	uint32_t geometry_id = blockIdx.x;

	uint16_t gs_point_id = blockIdx.y;

	GeometryInfo& m_geometry_info = context.get_geometry_info()[geometry_id];

	extern __shared__ double shrd_coefs[];

	//load coefs to shared memory
	detail::load_async<double, uint16_t>(m_geometry_info.coefs,
		m_geometry_info.num_CPoints * 3,
		shrd_coefs,
		true);

	uint16_t gs_kkjjii[3];

	get_gs_3d_id(gs_kkjjii, gs_point_id);

	//physics properties sent to PhysicContext
	double lambda = 576.92307692307691;
	double mu = 384.61538461538458;
	double e = 1000.0000000000000;
	double pr = 0.3;

	lambda = e * pr / ((1. + pr) * (1. - 2. * pr));
	mu = e / (2. * (1. + pr));

	Eigen::Matrix<double, 3, 3> identity = Eigen::Matrix<double, 3, 3>::Identity();

	Eigen::Matrix<double, 64, 3>* disp_block_geo = physic_context.get_cpoints_displacement(geometry_id);
	Eigen::Matrix<double, 64, 64>* basisValue0_matrix = integral_context.get_basis_value_0_matrix();
	Eigen::Matrix<double, 192, 64>* basisValue1_matrix = integral_context.get_basis_value_1_matrix();

	Eigen::Vector<double, 64> basis_value0 = basisValue0_matrix->col(gs_point_id);
	Eigen::Vector<double, 192> basis_value1 = basisValue1_matrix->col(gs_point_id);

	Eigen::Matrix<double, 3, 3> jacGeo, jacDis;

	get_jacobi_u<double>(basis_value1, *disp_block_geo, shrd_coefs, jacGeo, jacDis);

	double weight = integral_context.get_weights_u()[gs_kkjjii[2]] * integral_context.get_weights_v()[gs_kkjjii[1]] *
		integral_context.get_weights_w()[gs_kkjjii[0]] * fabs(jacGeo.determinant()) * density;

	double* massLocalMats = mass_context.get_mass_info()[geometry_id].mass_local_mats;

	Eigen::Matrix<double, 64, 64> block = weight * basis_value0 * basis_value0.transpose();

	constexpr uint16_t N_D = 4;

	uint16_t base_j = N_D * threadIdx.x;
	uint16_t base_k = N_D * threadIdx.y;

	constexpr uint16_t GS_D = 64;

	for (uint16_t local_j = 0; local_j < N_D; ++local_j) {
		uint16_t j = base_j + local_j;

		for (uint16_t local_k = 0; local_k < N_D; ++local_k) {
			uint16_t k = base_k + local_k;

			for (uint16_t d = 0; d < 3; ++d) {
				atomicAdd(&massLocalMats[(j + d * GS_D) + 192 * (k + d * GS_D)], block.coeff(j, k));
			}
		}
	}
}

__global__ static void IGA_add_dirichlet_to_mass_local_mats(
	const gpmesh::SplineContext context,
	const gpmesh::MassContext mass_context,
	const gpmesh::DirichletContext dirichlet_context
) {
	using namespace gpmesh;

	uint32_t dirichlet_id = blockIdx.x;
	uint16_t surf_id = threadIdx.x;

	uint16_t i = surf_id / 4;
	uint16_t j = surf_id - 4 * i;

	if (dirichlet_id >= dirichlet_context.get_num_dirichlet_boundary()) {
		return;
	}

	DirichletInfo dirichlet_info = dirichlet_context.get_dirichlet_info(dirichlet_id);

	assert(dirichlet_info.dirichlet_boundary_id == dirichlet_id);

	uint32_t geometry_id = dirichlet_info.geometry_id;
	uint16_t dim = dirichlet_info.dim;
	double value = dirichlet_info.value;

	MassInfo& mass_info = mass_context.get_mass_info(geometry_id);

	uint16_t local_index = INVALID16;

	assert(dirichlet_info.side != BoundarySide::NONE);

	uint16_t local_index_v2 = INVALID16;
	uint16_t i_multiply = INVALID16;
	uint16_t j_multiply = INVALID16;
	uint16_t index_add = INVALID16;

	switch (dirichlet_info.side) {
	case BoundarySide::EAST:
		i_multiply = 16;
		j_multiply = 4;
		index_add = 3;
		break;
	case BoundarySide::WEST:
		i_multiply = 16;
		j_multiply = 4;
		index_add = 0;
		break;
	case BoundarySide::NORTH:
		i_multiply = 16;
		j_multiply = 1;
		index_add = 12;
		break;
	case BoundarySide::SOUTH:
		i_multiply = 16;
		j_multiply = 1;
		index_add = 0;
		break;
	case BoundarySide::BACK:
		i_multiply = 4;
		j_multiply = 1;
		index_add = 48;
		break;
	case BoundarySide::FRONT:
		i_multiply = 4;
		j_multiply = 1;
		index_add = 0;
		break;
	}

	local_index_v2 = i_multiply * i + j_multiply * j + index_add;

	if (dirichlet_info.side == BoundarySide::WEST) {
		local_index = 16 * i + 4 * j;
	}
	else if (dirichlet_info.side == BoundarySide::EAST) {
		local_index = 16 * i + 4 * j + 3;
	}
	else if (dirichlet_info.side == BoundarySide::SOUTH) {
		local_index = 16 * i + j;
	}
	else if (dirichlet_info.side == BoundarySide::NORTH) {
		local_index = 16 * i + j + 12;
	}
	else if (dirichlet_info.side == BoundarySide::FRONT) {
		local_index = 4 * i + j;
	}
	else if (dirichlet_info.side == BoundarySide::BACK) {
		local_index = 4 * i + j + 48;
	}

	assert(local_index != INVALID16);
	assert(local_index == local_index_v2);

	constexpr double DirchletMaxNum = 1e10;

	//mass_info.mass_local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)] += DirchletMaxNum;
	atomicAdd(&(mass_info.mass_local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)]), DirchletMaxNum);
}

__global__ static void IGA_add_dirichlet_point_to_mass_local_mats(
	const gpmesh::SplineContext context,
	const gpmesh::MassContext mass_context,
	const gpmesh::DirichletPointContext dirichlet_point_context,
	uint32_t num_dirichlet_point
) {
	using namespace gpmesh;

	uint32_t dirichlet_id = blockIdx.x * blockDim.x + threadIdx.x;

	if (dirichlet_id >= num_dirichlet_point) return;


	DirichletPointInfo dirichlet_point_info = dirichlet_point_context.get_dirichlet_point_info(dirichlet_id);


	uint32_t geometry_id = dirichlet_point_info.geometry_id;
	uint16_t local_id = dirichlet_point_info.local_id;


	MassInfo& mass_info = mass_context.get_mass_info(geometry_id);

	constexpr double DirchletMaxNum = 1e10;

	//mass_info.mass_local_mats[(local_index + dim * 64) + 192 * (local_index + dim * 64)] += DirchletMaxNum;
	atomicAdd(&(mass_info.mass_local_mats[(local_id + 0 * 64) + 192 * (local_id + 0 * 64)]), DirchletMaxNum);
	atomicAdd(&(mass_info.mass_local_mats[(local_id + 1 * 64) + 192 * (local_id + 1 * 64)]), DirchletMaxNum);
	atomicAdd(&(mass_info.mass_local_mats[(local_id + 2 * 64) + 192 * (local_id + 2 * 64)]), DirchletMaxNum);
}

//<<<geometry_size, 256>>>
template<typename T, size_t blockThreads>
__global__ static void IGA_add_mass_matrix_to_stiff_matrix(
	const gpmesh::AssemblyContext assembly_context,
	const gpmesh::MassContext mass_context,
	double coeff) {
	uint32_t geometry_id = blockIdx.x;

	double* massLocalMats = mass_context.get_mass_info(geometry_id).mass_local_mats;
	double* assemblyLocalMats = assembly_context.get_assembly_info(geometry_id).local_mats;

	for (uint32_t i = threadIdx.x; i < 192 * 192; i += blockThreads) {
		assemblyLocalMats[i] += coeff * massLocalMats[i];
	}
}


//<<<geometry_size, 192>>>
template<typename T>
__global__ static void IGA_update_physics_context_by_spline_attr_global(
	const gpmesh::SplineAttributeGlobal<T> update,
	const gpmesh::PhysicContext physic_context,
	const gpmesh::SplineAttributeGlobalMappingContext mapping_context
) {
	uint32_t geometry_id = blockIdx.x;
	uint16_t col_id = threadIdx.x;

	uint16_t dim = col_id / 64;
	uint16_t local_id = col_id - dim * 64;

	Eigen::Matrix<T, 64, 3>* local_cpoint_disp_mat = physic_context.get_cpoints_displacement(geometry_id);

	local_cpoint_disp_mat->operator()(local_id, dim) += update(geometry_id, local_id, dim);
}

//<<<num_geometry, 64*3>>>
template<typename T>
__global__ static void
IGA_kernel_build_preconditioner_from_system_matrix(const gpmesh::SplineAttributeGlobal<T> preconditioner,
	const gpmesh::AssemblyContext assembly_context
) {
	uint32_t geometry_id = blockIdx.x;
	uint32_t thread_id = threadIdx.x;
	assert(thread_id < 192);

	//add to local part
	preconditioner.RawOperator(geometry_id, thread_id) = assembly_context.get_assembly_info(geometry_id).local_mats[
		thread_id * 192 + thread_id];
}

//<<<DIVE_UP(threads, 256), 256)>>>
template<typename T>
__global__ static void IGA_kernel_convert_preconditioner_reciprocal(
	const gpmesh::SplineAttributeGlobal<T> preconditioner,
	uint32_t num_threads
) {
	uint32_t id = blockDim.x * blockIdx.x + threadIdx.x;

	if (id >= num_threads) return;

	assert(id < num_threads);

	T val = preconditioner.m_d_global_attr[id];

	if (abs(val) < 1e-10) {
		val = 1;
	}
	else {
		val = 1.0 / val;
	}

	preconditioner.m_d_global_attr[id] = val;
}

//<<<DIVE_UP(threads, 256), 256)>>>
//res CWise* = other
template<typename T>
__global__ static void IGA_kernel_CWise_operation_dot(const gpmesh::SplineAttributeGlobal<T> res,
	const gpmesh::SplineAttributeGlobal<T> other,
	uint32_t num_threads) {
	uint32_t id = blockDim.x * blockIdx.x + threadIdx.x;

	if (id >= num_threads) return;

	assert(id < num_threads);

	res.m_d_global_attr[id] = res.m_d_global_attr[id]
		* other.m_d_global_attr[id];
}

//<<<DIVE_UP(threads, 256), 256)>>>
//res = A CWise* B
template<typename T>
__global__ static void IGA_kernel_CWise_operation_dot(const gpmesh::SplineAttributeGlobal<T> res,
	const gpmesh::SplineAttributeGlobal<T> A,
	const gpmesh::SplineAttributeGlobal<T> B,
	uint32_t num_threads) {
	uint32_t id = blockDim.x * blockIdx.x + threadIdx.x;

	if (id >= num_threads) return;

	assert(id < num_threads);

	res.m_d_global_attr[id] = A.m_d_global_attr[id]
		* B.m_d_global_attr[id];
}

__global__ static void IGA_clear_physics_displacement(const gpmesh::PhysicContext physic_context)
{
	uint32_t geometry_id = blockIdx.x;
	uint32_t thread_id = threadIdx.x;
	uint16_t row = thread_id / 3;
	uint16_t col = thread_id % 3;
	physic_context.get_cpoints_displacement(geometry_id)->operator()(row, col) = 0;
}
