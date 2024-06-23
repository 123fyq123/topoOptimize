#include "BSplineSolid.h"
#include<iomanip>
//求点x在T中的位置
int BSplineSolid::X_LOCATION(int solid_id, int uvw, double x)
{
	int n;
	n = Local_Node[solid_id][uvw].size() - Local_Order[solid_id][uvw] - 2;//控制顶点个数-1
	if (fabs(x - Local_Node[solid_id][uvw][n + 1]) < 1.e-6) return n;
	int low, high;
	low = Local_Order[solid_id][uvw]; high = n + 1;
	int mid = (low + high) / 2;
	while (x < Local_Node[solid_id][uvw][mid] || x >= Local_Node[solid_id][uvw][mid + 1])
	{
		if (x < Local_Node[solid_id][uvw][mid]) high = mid;
		else low = mid;
		mid = (low + high) / 2;
	}
	return mid;
}

//m表示节点向量的最大下标，即节点向量共有m+1个节点，i表示基Ni（i从0开始算），
//n表示求导到n次，结果在ders中（ders共n+1个元素）。
//void DersOneBasisFun(int degree,int m,vector<double> U, int i,double u,int n,double ders[])
void BSplineSolid::DersOneBasisFun(int solid_id, int uvw, int i, double u, int n, double ders[])
{
	double** N, * ND;
	int j, k;
	int degree = Local_Order[solid_id][uvw];                        //3
	N = new double* [degree + 1];
	ND = new double[degree + 1];
	for (j = 0; j <= degree; j++)
		N[j] = new double[degree + 1];

	double saved, temp, Uleft, Uright;


	// 参数u 不在该单元的节点
	if (u < Local_Node[solid_id][uvw][i] || u > Local_Node[solid_id][uvw][i + degree + 1])
	{
		for (k = 0; k <= n; k++)
			ders[k] = 0.0;
		return;
	}

	double m_u;
	m_u = u;
	int m = Local_Node[solid_id][uvw].size() - 1;
	if (u == Local_Node[solid_id][uvw][m])
	{
		i = m - degree - 1 - i;
		u = Local_Node[solid_id][uvw][degree];
	}
	for (j = 0; j <= degree; j++)
		if (u >= Local_Node[solid_id][uvw][i + j] && u < Local_Node[solid_id][uvw][i + j + 1])  N[j][0] = 1.0;
		else N[j][0] = 0.0;
	for (k = 1; k <= degree; k++)
	{
		if (N[0][k - 1] == 0.0) saved = 0.0;
		else saved = ((u - Local_Node[solid_id][uvw][i]) * N[0][k - 1]) / (Local_Node[solid_id][uvw][i + k] - Local_Node[solid_id][uvw][i]);
		for (j = 0; j < degree - k + 1; j++)
		{
			Uleft = Local_Node[solid_id][uvw][i + j + 1];
			Uright = Local_Node[solid_id][uvw][i + j + k + 1];
			if (N[j + 1][k - 1] == 0.0)
			{
				N[j][k] = saved; saved = 0.0;
			}
			else
			{
				temp = N[j + 1][k - 1] / (Uright - Uleft);
				N[j][k] = saved + (Uright - u) * temp;
				saved = (u - Uleft) * temp;
			}
		}
	}
	ders[0] = N[0][degree];

	for (k = 1; k <= n; k++)
	{
		for (j = 0; j <= k; j++)
			ND[j] = N[j][degree - k];
		for (int jj = 1; jj <= k; jj++)
		{
			if (ND[0] == 0.0) saved = 0.0;
			else saved = ND[0] / (Local_Node[solid_id][uvw][i + degree - k + jj] - Local_Node[solid_id][uvw][i]);
			for (j = 0; j < k - jj + 1; j++)
			{
				Uleft = Local_Node[solid_id][uvw][i + j + 1];
				Uright = Local_Node[solid_id][uvw][i + j + degree - k + jj + 1];
				if (ND[j + 1] == 0.0)
				{
					ND[j] = (degree - k + jj) * saved; saved = 0.0;
				}
				else
				{
					temp = ND[j + 1] / (Uright - Uleft);
					ND[j] = (degree - k + jj) * (saved - temp);
					saved = temp;
				}
			}
		}
		ders[k] = ND[0];
	}

	//仅考虑一阶导数的情形。奇数阶导数变号，偶数阶不变号
	if (m_u == Local_Node[solid_id][uvw][m])
	{
		ders[1] = -ders[1];
	}
	for (j = 0; j <= degree; j++)
		delete[] N[j];
	delete[] N;
	delete[] ND;
}

//计算几何变换.X=G(u,v),也就是曲面在u,v的值。
V3f BSplineSolid::Geometry(int solid_id, double u, double v, double w)
{
	int m_uspan, m_vspan, m_wspan;
	double m_x, m_y, m_z, u_ders[2], v_ders[2], w_ders[2];
	int u_order = Local_Order[solid_id][0] + 1;
	int v_order = Local_Order[solid_id][1] + 1;
	int w_order = Local_Order[solid_id][2] + 1;
	V3f vertex(0, 0, 0);
	m_x = m_y = m_z = 0.0;
	m_uspan = X_LOCATION(solid_id, 0, u);
	m_vspan = X_LOCATION(solid_id, 1, v);
	m_wspan = X_LOCATION(solid_id, 2, w);
	for (int l = m_uspan - u_order + 1; l <= m_uspan; l++)
	{
		DersOneBasisFun(solid_id, 0, l, u, 0, u_ders);
		for (int k = m_vspan - v_order + 1; k <= m_vspan; k++)
		{
			DersOneBasisFun(solid_id, 1, k, v, 0, v_ders);
			for (int p = m_wspan - w_order + 1; p <= m_wspan; p++)
			{
				DersOneBasisFun(solid_id, 2, p, w, 0, w_ders);
				/*int m_id = l*(Local_Node[solid_id][0].size() - v_order) + k +
					p*(Local_Node[solid_id][0].size() - v_order)*(Local_Node[solid_id][1].size() - u_order);*/
				vertex = vertex + Local_BSplineSolid[solid_id].vertices_[Control_Point_id[solid_id][p][k][l]] * u_ders[0] * v_ders[0] * w_ders[0];
				//m_y += CPts[m_id].y*u_ders[0] * v_ders[0];
			}
		}
	}
	//X = m_x; Y = m_y;
	return vertex;
}

V3f BSplineSolid::Geometry_new(int solid_id, double u, double v, double w)
{
	int m_uspan, m_vspan, m_wspan;
	double m_x, m_y, m_z, u_ders[2], v_ders[2], w_ders[2];
	int u_order = Local_Order[solid_id][0] + 1;
	int v_order = Local_Order[solid_id][1] + 1;
	int w_order = Local_Order[solid_id][2] + 1;
	V3f vertex(0, 0, 0);
	m_x = m_y = m_z = 0.0;
	m_uspan = X_LOCATION(solid_id, 0, u);
	m_vspan = X_LOCATION(solid_id, 1, v);
	m_wspan = X_LOCATION(solid_id, 2, w);
	for (int l = m_uspan - u_order + 1; l <= m_uspan; l++)
	{
		DersOneBasisFun(solid_id, 0, l, u, 0, u_ders);
		for (int k = m_vspan - v_order + 1; k <= m_vspan; k++)
		{
			DersOneBasisFun(solid_id, 1, k, v, 0, v_ders);
			for (int p = m_wspan - w_order + 1; p <= m_wspan; p++)
			{
				DersOneBasisFun(solid_id, 2, p, w, 0, w_ders);
				vertex = vertex + Local_V3f[solid_id][Control_Point[solid_id][p][k][l]] * u_ders[0] * v_ders[0] * w_ders[0];
			}
		}
	}
	return vertex;
}

V3f BSplineSolid::Geometry_new(int solid_id, int uid, int vid, int wid, const std::vector<double>& PreCalDerBasisFunction)
{
	int m_uspan, m_vspan, m_wspan;
	double m_x, m_y, m_z, u_ders[2], v_ders[2], w_ders[2];
	int u_order = Local_Order[solid_id][0] + 1;
	int v_order = Local_Order[solid_id][1] + 1;
	int w_order = Local_Order[solid_id][2] + 1;
	V3f vertex(0, 0, 0);
	m_x = m_y = m_z = 0.0;
	m_uspan = 3; // X_LOCATION(solid_id, 0, u);
	m_vspan = 3; // X_LOCATION(solid_id, 1, v);
	m_wspan = 3; // X_LOCATION(solid_id, 2, w);
	for (int l = m_uspan - u_order + 1; l <= m_uspan; l++)
	{
		//DersOneBasisFun(solid_id, 0, l, u, 0, u_ders);
		u_ders[0] = PreCalDerBasisFunction[0 * 32 + (l * 1 + uid * 4) * 2];
		for (int k = m_vspan - v_order + 1; k <= m_vspan; k++)
		{
			//DersOneBasisFun(solid_id, 1, k, v, 0, v_ders);
			v_ders[0] = PreCalDerBasisFunction[1 * 32 + (k * 1 + vid * 4) * 2];
			for (int p = m_wspan - w_order + 1; p <= m_wspan; p++)
			{
				//DersOneBasisFun(solid_id, 2, p, w, 0, w_ders);
				w_ders[0] = PreCalDerBasisFunction[2 * 32 + (p * 1 + wid * 4) * 2];
				vertex = vertex + Local_V3f[solid_id][Control_Point[solid_id][p][k][l]] * u_ders[0] * v_ders[0] * w_ders[0];
			}
		}
	}
	return vertex;
}

//对边界函数近似求基的系数，m_bdy[0]=1表示u方向，m_bdy[1]=1表示v方向，m_bdy[2]=1表示w方向, m_para表示另一个方向的参数值
void BSplineSolid::SolveCoef(int solid_id, std::vector<double> m_scoef, std::vector<int> m_para,
	std::vector<int> m_bdy, std::vector<double>& m_coef)
{
	std::vector< std::vector<double> > m_matrixvec;
	std::vector<double> m_rhsvec;
	double m_der[2], m_value;
	int m_nptnum, m_cptnum;
	V3f m_pt;
	//判断是什么方向的边 以后不用重复判断
	int uvw;
	if (m_bdy[0] == 1) uvw = 0;
	else if (m_bdy[1] == 1) uvw = 1;
	else if (m_bdy[2] == 1) uvw = 2;
	else return;
	int m_uptnum = Local_Node[solid_id][0].size() - Local_Order[solid_id][0] - 1;
	int m_vptnum = Local_Node[solid_id][1].size() - Local_Order[solid_id][1] - 1;
	int m_wptnum = Local_Node[solid_id][2].size() - Local_Order[solid_id][2] - 1;
	int sameuvw = m_wptnum;
	m_coef.clear();
	m_cptnum = Local_Node[solid_id][uvw].size() - Local_Order[solid_id][uvw] - 1;
	m_nptnum = 6 * m_cptnum;
	m_matrixvec.resize(m_nptnum - 1);
	for (int j = 1; j < m_nptnum; j++)
	{
		for (int k = 0; k < m_cptnum - 2; k++)
			m_matrixvec[j - 1].push_back(0);//m_nptnum-1的size的m_matrixvec，每个元素有m_cptnum-2个
	}
	for (int j = 1; j < m_nptnum; j++)
	{
		double m_tmp = 1.0 * j / m_nptnum;
		for (int k = 1; k < m_cptnum - 1; k++)
		{
			DersOneBasisFun(solid_id, uvw, k, m_tmp, 1, m_der);
			m_matrixvec[j - 1][k - 1] = m_der[0];
		}
		if (uvw == 0)
			m_pt = Geometry(solid_id, m_tmp, m_para[1], m_para[2]);
		else if (uvw == 1)
			m_pt = Geometry(solid_id, m_para[0], m_tmp, m_para[2]);
		else if (uvw == 2)
			m_pt = Geometry(solid_id, m_para[0], m_para[1], m_tmp);
		else return;
		m_value = ExactSol(m_pt.x, m_pt.y, m_pt.z);
		DersOneBasisFun(solid_id, uvw, 0, m_tmp, 1, m_der);
		if (uvw == 0)
		{
			if (m_para[1] == 0 && m_para[2] == 0)//0,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][0]];
			else if (m_para[1] == 0 && m_para[2] == 1)//0,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][0][0]];
			else if (m_para[1] == 1 && m_para[2] == 0)//0,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][0]];
			else if (m_para[1] == 1 && m_para[2] == 1)//0,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][0]];
		}
		else if (uvw == 1)
		{
			if (m_para[0] == 0 && m_para[2] == 0)//0,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][0]];
			else if (m_para[0] == 0 && m_para[2] == 1)//0,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][0][0]];
			else if (m_para[0] == 1 && m_para[2] == 0)//1,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][sameuvw - 1]];
			else if (m_para[0] == 1 && m_para[2] == 1)//1,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][0][sameuvw - 1]];
		}
		else if (uvw == 2)
		{
			if (m_para[0] == 0 && m_para[1] == 0)//0,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][0]];
			else if (m_para[0] == 0 && m_para[1] == 1)//0,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][0]];
			else if (m_para[0] == 1 && m_para[1] == 0)//1,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][sameuvw - 1]];
			else if (m_para[0] == 1 && m_para[1] == 1)//1,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][sameuvw - 1]];
		}
		DersOneBasisFun(solid_id, uvw, m_cptnum - 1, m_tmp, 1, m_der);
		if (uvw == 0)
		{
			if (m_para[1] == 0 && m_para[2] == 0)//1,0,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][sameuvw - 1]];
			else if (m_para[1] == 0 && m_para[2] == 1)//1,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][0][sameuvw - 1]];
			else if (m_para[1] == 1 && m_para[2] == 0)//1,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][sameuvw - 1]];
			else if (m_para[1] == 1 && m_para[2] == 1)//1,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][sameuvw - 1]];
		}
		else if (uvw == 1)
		{
			if (m_para[0] == 0 && m_para[2] == 0)//0,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][0]];
			else if (m_para[0] == 0 && m_para[2] == 1)//0,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][0]];
			else if (m_para[0] == 1 && m_para[2] == 0)//1,1,0
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][sameuvw - 1][sameuvw - 1]];
			else if (m_para[0] == 1 && m_para[2] == 1)//1,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][sameuvw - 1]];
		}
		else if (uvw == 2)
		{
			if (m_para[0] == 0 && m_para[1] == 0)//0,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][0][0][sameuvw - 1]];
			else if (m_para[0] == 0 && m_para[1] == 1)//0,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][0]];
			else if (m_para[0] == 1 && m_para[1] == 0)//1,0,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][0][sameuvw - 1]];
			else if (m_para[0] == 1 && m_para[1] == 1)//1,1,1
				m_value = m_value - m_der[0] * m_scoef[Control_Point_id[solid_id][sameuvw - 1][sameuvw - 1][sameuvw - 13]];
		}
		m_rhsvec.push_back(m_value);
	}

	double* m_matrix, * m_rvec;
	int m_size = m_cptnum - 2;
	m_matrix = new double[m_size * m_size];
	m_rvec = new double[m_size];
	for (int k = 0; k < m_size; k++)
	{
		double m_tmp;
		for (int j = 0; j < m_size; j++)
		{
			m_tmp = 0;
			for (int nn = 0; nn < m_nptnum - 1; nn++)
				m_tmp = m_tmp + m_matrixvec[nn][k] * m_matrixvec[nn][j];
			m_matrix[k * m_size + j] = m_tmp;
		}
		m_tmp = 0;
		for (int j = 0; j < m_nptnum - 1; j++)
			m_tmp = m_tmp + m_matrixvec[j][k] * m_rhsvec[j];
		m_rvec[k] = m_tmp;
	}


	int M = m_size;
	int N = 1;
	int lda, ldb;
	int* ipiv = new int[m_size];
	int INFO;
	lda = M;
	ldb = M;
	dgesv_((int*)&M, (int*)&N, m_matrix, (int*)&lda, (int*)ipiv, m_rvec, (int*)&ldb, (int*)&INFO);
	for (int k = 1; k < m_cptnum - 1; k++)
		m_coef.push_back(m_rvec[k - 1]);
	delete[]ipiv;
	delete[]m_matrix;
	delete[]m_rvec;
}

//寻找结点向量中所含的区域个数
void BSplineSolid::IntegralArea(int solid_id, std::vector<double>& u_interval, std::vector<double>& v_interval,
	std::vector<double>& w_interval)
{
	int i;
	double uvw;
	u_interval.clear(); v_interval.clear(); w_interval.clear();

	i = 0;
	do
	{
		uvw = Local_Node[solid_id][0][i];
		if (uvw == Local_Node[solid_id][0][i + 1])
		{
			if (i == (int)Local_Node[solid_id][0].size() - 2) u_interval.push_back(Local_Node[solid_id][0][i + 1]);
			i++;
		}
		else
		{
			u_interval.push_back(uvw);
			if (i == (int)Local_Node[solid_id][0].size() - 2) u_interval.push_back(Local_Node[solid_id][0][i + 1]);
			i++;
		}
	} while (i < (int)Local_Node[solid_id][0].size() - 1);

	i = 0;
	do
	{
		uvw = Local_Node[solid_id][1][i];
		if (uvw == Local_Node[solid_id][1][i + 1])
		{
			if (i == (int)Local_Node[solid_id][1].size() - 2) v_interval.push_back(Local_Node[solid_id][1][i + 1]);
			i++;
		}
		else
		{
			v_interval.push_back(uvw);
			if (i == (int)Local_Node[solid_id][1].size() - 2) v_interval.push_back(Local_Node[solid_id][1][i + 1]);
			i++;
		}
	} while (i < (int)Local_Node[solid_id][1].size() - 1);

	i = 0;
	do
	{
		uvw = Local_Node[solid_id][2][i];
		if (uvw == Local_Node[solid_id][2][i + 1])
		{
			if (i == (int)Local_Node[solid_id][2].size() - 2) w_interval.push_back(Local_Node[solid_id][2][i + 1]);
			i++;
		}
		else
		{
			w_interval.push_back(uvw);
			if (i == (int)Local_Node[solid_id][2].size() - 2) w_interval.push_back(Local_Node[solid_id][2][i + 1]);
			i++;
		}
	} while (i < (int)Local_Node[solid_id][2].size() - 1);
}

//计算一系列点上的雅可比矩阵及行列式
void BSplineSolid::JacobiMatrix(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
	std::vector<double>& jacobidet, std::vector<Matrix>& jacobi)
{
	int m_umaxindex = Local_Node[solid_id][0].size() - 1;
	int m_vmaxindex = Local_Node[solid_id][1].size() - 1;
	int m_wmaxindex = Local_Node[solid_id][2].size() - 1;
	int m_uptnum = m_umaxindex - Local_Order[solid_id][0];
	int m_vptnum = m_vmaxindex - Local_Order[solid_id][1];
	int m_wptnum = m_wmaxindex - Local_Order[solid_id][2];
	double m_ubasevalue[2], m_vbasevalue[2], m_wbasevalue[2];

	jacobidet.clear();
	jacobi.clear();

	Matrix m_jacobi;
	double m_jacobidet;

	m_jacobi.resize(3);
	m_jacobi[0].resize(3); m_jacobi[1].resize(3); m_jacobi[2].resize(3);

	for (int m = 0; m < w.size(); m++)
	{
		for (int l = 0; l < v.size(); l++)
		{
			for (int k = 0; k < u.size(); k++)
			{
				double m_value[3] = { 0,0,0 };
				double m_uder[3] = { 0,0,0 };
				double m_vder[3] = { 0,0,0 };
				double m_wder[3] = { 0,0,0 };
				int m_uIndex = X_LOCATION(solid_id, 0, u[k]);
				int m_vIndex = X_LOCATION(solid_id, 1, v[l]);
				int m_wIndex = X_LOCATION(solid_id, 2, w[m]);
				for (int i = m_uIndex - Local_Order[solid_id][0]; i <= m_uIndex; i++)
				{
					DersOneBasisFun(solid_id, 0, i, u[k], 1, m_ubasevalue);
					for (int j = m_vIndex - Local_Order[solid_id][1]; j <= m_vIndex; j++)
					{
						DersOneBasisFun(solid_id, 1, j, v[l], 1, m_vbasevalue);
						for (int n = m_wIndex - Local_Order[solid_id][2]; n <= m_wIndex; n++)
						{
							DersOneBasisFun(solid_id, 2, n, w[m], 1, m_wbasevalue);

							V3f vertex = Local_BSplineSolid[solid_id].vertices_[Control_Point_id[solid_id][n][j][i]];
							m_value[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];

							m_uder[0] += vertex.x * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[1] += vertex.y * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[2] += vertex.z * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];

							m_vder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];

							m_wder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
						}
					}
				}
				/*double check1 = m_uder[0] * (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]);
				double check2 = -m_vder[0] * (m_uder[1] * m_wder[2] - m_wder[1] * m_uder[2]);
				double check3 = m_wder[0] * (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]);*/
				m_jacobidet = m_uder[0] * (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) - m_vder[0] * (m_uder[1] * m_wder[2] - m_wder[1] * m_uder[2])
					+ m_wder[0] * (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]);
				m_jacobi[0][0] = (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) / m_jacobidet;
				m_jacobi[0][1] = (m_wder[1] * m_uder[2] - m_uder[1] * m_wder[2]) / m_jacobidet;
				m_jacobi[0][2] = (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]) / m_jacobidet;
				m_jacobi[1][0] = (m_wder[0] * m_vder[2] - m_vder[0] * m_wder[2]) / m_jacobidet;
				m_jacobi[1][1] = (m_uder[0] * m_wder[2] - m_wder[0] * m_uder[2]) / m_jacobidet;
				m_jacobi[1][2] = (m_vder[0] * m_uder[2] - m_uder[0] * m_vder[2]) / m_jacobidet;
				m_jacobi[2][0] = (m_vder[0] * m_wder[1] - m_wder[0] * m_vder[1]) / m_jacobidet;
				m_jacobi[2][1] = (m_wder[0] * m_uder[1] - m_uder[0] * m_wder[1]) / m_jacobidet;
				m_jacobi[2][2] = (m_uder[0] * m_vder[1] - m_vder[0] * m_uder[1]) / m_jacobidet;
				jacobidet.push_back(m_jacobidet);
				jacobi.push_back(m_jacobi);
			}
		}
	}
}

void BSplineSolid::JacobiMatrix_new(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
	std::vector<double>& jacobidet, std::vector<Matrix>& jacobi)
{
	int m_umaxindex = Local_Node[solid_id][0].size() - 1;
	int m_vmaxindex = Local_Node[solid_id][1].size() - 1;
	int m_wmaxindex = Local_Node[solid_id][2].size() - 1;
	int m_uptnum = m_umaxindex - Local_Order[solid_id][0];
	int m_vptnum = m_vmaxindex - Local_Order[solid_id][1];
	int m_wptnum = m_wmaxindex - Local_Order[solid_id][2];
	double m_ubasevalue[2], m_vbasevalue[2], m_wbasevalue[2];

	jacobidet.clear();
	jacobi.clear();

	Matrix m_jacobi;
	double m_jacobidet;

	m_jacobi.resize(3);
	m_jacobi[0].resize(3); m_jacobi[1].resize(3); m_jacobi[2].resize(3);

	for (int m = 0; m < w.size(); m++)
	{
		for (int l = 0; l < v.size(); l++)
		{
			for (int k = 0; k < u.size(); k++)
			{
				double m_value[3] = { 0,0,0 };
				double m_uder[3] = { 0,0,0 };
				double m_vder[3] = { 0,0,0 };
				double m_wder[3] = { 0,0,0 };
				int m_uIndex = X_LOCATION(solid_id, 0, u[k]);
				int m_vIndex = X_LOCATION(solid_id, 1, v[l]);
				int m_wIndex = X_LOCATION(solid_id, 2, w[m]);
				for (int i = m_uIndex - Local_Order[solid_id][0]; i <= m_uIndex; i++)
				{
					DersOneBasisFun(solid_id, 0, i, u[k], 1, m_ubasevalue);
					for (int j = m_vIndex - Local_Order[solid_id][1]; j <= m_vIndex; j++)
					{
						DersOneBasisFun(solid_id, 1, j, v[l], 1, m_vbasevalue);
						for (int n = m_wIndex - Local_Order[solid_id][2]; n <= m_wIndex; n++)
						{
							DersOneBasisFun(solid_id, 2, n, w[m], 1, m_wbasevalue);

							V3f vertex = Local_V3f[solid_id][Control_Point[solid_id][n][j][i]];
							m_value[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];

							m_uder[0] += vertex.x * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[1] += vertex.y * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[2] += vertex.z * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];

							m_vder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];

							m_wder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
						}
					}
				}
				m_jacobidet = m_uder[0] * (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) - m_vder[0] * (m_uder[1] * m_wder[2] - m_wder[1] * m_uder[2])
					+ m_wder[0] * (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]);
				if (m_jacobidet == 0) std::cout << "m_jacobidet == 0" << std::endl;
				m_jacobi[0][0] = (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) / m_jacobidet;
				m_jacobi[0][1] = (m_wder[1] * m_uder[2] - m_uder[1] * m_wder[2]) / m_jacobidet;
				m_jacobi[0][2] = (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]) / m_jacobidet;
				m_jacobi[1][0] = (m_wder[0] * m_vder[2] - m_vder[0] * m_wder[2]) / m_jacobidet;
				m_jacobi[1][1] = (m_uder[0] * m_wder[2] - m_wder[0] * m_uder[2]) / m_jacobidet;
				m_jacobi[1][2] = (m_vder[0] * m_uder[2] - m_uder[0] * m_vder[2]) / m_jacobidet;
				m_jacobi[2][0] = (m_vder[0] * m_wder[1] - m_wder[0] * m_vder[1]) / m_jacobidet;
				m_jacobi[2][1] = (m_wder[0] * m_uder[1] - m_uder[0] * m_wder[1]) / m_jacobidet;
				m_jacobi[2][2] = (m_uder[0] * m_vder[1] - m_vder[0] * m_uder[1]) / m_jacobidet;
				jacobidet.push_back(m_jacobidet);
				jacobi.push_back(m_jacobi);
			}
		}
	}
}

void BSplineSolid::JacobiMatrix_new(int solid_id, std::vector<double> u, std::vector<double> v, std::vector<double> w,
	std::vector<double>& jacobidet, std::vector<Matrix>& jacobi, std::vector<double>& PreCalDerBasisFunction, bool isPreCal)
{
	int m_umaxindex = Local_Node[solid_id][0].size() - 1;
	int m_vmaxindex = Local_Node[solid_id][1].size() - 1;
	int m_wmaxindex = Local_Node[solid_id][2].size() - 1;
	int m_uptnum = m_umaxindex - Local_Order[solid_id][0];
	int m_vptnum = m_vmaxindex - Local_Order[solid_id][1];
	int m_wptnum = m_wmaxindex - Local_Order[solid_id][2];

	double m_ubasevalue[2], m_vbasevalue[2], m_wbasevalue[2];

	jacobidet.clear();
	jacobi.clear();

	Matrix m_jacobi;
	double m_jacobidet;

	m_jacobi.resize(3);
	m_jacobi[0].resize(3); m_jacobi[1].resize(3); m_jacobi[2].resize(3);
	PreCalDerBasisFunction.resize(24 * 4);
	for (int m = 0; m < w.size(); m++)
	{
		for (int l = 0; l < v.size(); l++)
		{
			for (int k = 0; k < u.size(); k++)
			{

				// 每一个采样点

				double m_value[3] = { 0,0,0 };
				double m_uder[3] = { 0,0,0 };
				double m_vder[3] = { 0,0,0 };
				double m_wder[3] = { 0,0,0 };
				int m_uIndex = X_LOCATION(solid_id, 0, u[k]);     //3
				int m_vIndex = X_LOCATION(solid_id, 1, v[l]);
				int m_wIndex = X_LOCATION(solid_id, 2, w[m]);

				for (int i = m_uIndex - Local_Order[solid_id][0]; i <= m_uIndex; i++)
				{
					int id = (i * 1 + k * 4) * 2 + 0 * 32;
					if (!isPreCal) {
						DersOneBasisFun(solid_id, 0, i, u[k], 1, m_ubasevalue);
						PreCalDerBasisFunction[id] = m_ubasevalue[0];
						PreCalDerBasisFunction[id + 1] = m_ubasevalue[1];
					}
					else {
						//DersOneBasisFun(solid_id, 0, i, u[k], 1, m_ubasevalue);
						m_ubasevalue[0] = PreCalDerBasisFunction[id];
						m_ubasevalue[1] = PreCalDerBasisFunction[id + 1];
					}

					for (int j = m_vIndex - Local_Order[solid_id][1]; j <= m_vIndex; j++)
					{
						int id = (j * 1 + l * 4) * 2 + 1 * 32;
						if (!isPreCal) {

							DersOneBasisFun(solid_id, 1, j, v[l], 1, m_vbasevalue);
							PreCalDerBasisFunction[id] = m_vbasevalue[0];
							PreCalDerBasisFunction[id + 1] = m_vbasevalue[1];
						}
						else {
							//DersOneBasisFun(solid_id, 1, j, v[l], 1, m_vbasevalue);
							m_vbasevalue[0] = PreCalDerBasisFunction[id];
							m_vbasevalue[1] = PreCalDerBasisFunction[id + 1];
						}

						for (int n = m_wIndex - Local_Order[solid_id][2]; n <= m_wIndex; n++)
						{
							int id = (n * 1 + m * 4) * 2 + 2 * 32;
							if (!isPreCal) {
								DersOneBasisFun(solid_id, 2, n, w[m], 1, m_wbasevalue);
								PreCalDerBasisFunction[id] = m_wbasevalue[0];
								PreCalDerBasisFunction[id + 1] = m_wbasevalue[1];
							}
							else {
								//DersOneBasisFun(solid_id, 2, n, w[m], 1, m_wbasevalue);
								m_wbasevalue[0] = PreCalDerBasisFunction[id];
								m_wbasevalue[1] = PreCalDerBasisFunction[id + 1];
							}


							V3f vertex = Local_V3f[solid_id][Control_Point[solid_id][n][j][i]];

							m_value[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];
							m_value[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[0];

							m_uder[0] += vertex.x * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[1] += vertex.y * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];
							m_uder[2] += vertex.z * m_ubasevalue[1] * m_vbasevalue[0] * m_wbasevalue[0];

							m_vder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];
							m_vder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[1] * m_wbasevalue[0];

							m_wder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
							m_wder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0] * m_wbasevalue[1];
						}
					}
				}
				m_jacobidet = m_uder[0] * (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) - m_vder[0] * (m_uder[1] * m_wder[2] - m_wder[1] * m_uder[2])
					+ m_wder[0] * (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]);
				//if (m_jacobidet <= 0) std::cout <<"m_jacobidet " << m_jacobidet << std::endl;
				m_jacobi[0][0] = (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) / m_jacobidet;
				m_jacobi[0][1] = (m_wder[1] * m_uder[2] - m_uder[1] * m_wder[2]) / m_jacobidet;
				m_jacobi[0][2] = (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]) / m_jacobidet;

				m_jacobi[1][0] = (m_wder[0] * m_vder[2] - m_vder[0] * m_wder[2]) / m_jacobidet;
				m_jacobi[1][1] = (m_uder[0] * m_wder[2] - m_wder[0] * m_uder[2]) / m_jacobidet;
				m_jacobi[1][2] = (m_vder[0] * m_uder[2] - m_uder[0] * m_vder[2]) / m_jacobidet;

				m_jacobi[2][0] = (m_vder[0] * m_wder[1] - m_wder[0] * m_vder[1]) / m_jacobidet;
				m_jacobi[2][1] = (m_wder[0] * m_uder[1] - m_uder[0] * m_wder[1]) / m_jacobidet;
				m_jacobi[2][2] = (m_uder[0] * m_vder[1] - m_vder[0] * m_uder[1]) / m_jacobidet;
				/*m_jacobi[0][0] = (m_vder[1] * m_wder[2] - m_wder[1] * m_vder[2]) / m_jacobidet;
				m_jacobi[1][0] = (m_wder[1] * m_uder[2] - m_uder[1] * m_wder[2]) / m_jacobidet;
				m_jacobi[2][0] = (m_uder[1] * m_vder[2] - m_vder[1] * m_uder[2]) / m_jacobidet;

				m_jacobi[0][1] = (m_wder[0] * m_vder[2] - m_vder[0] * m_wder[2]) / m_jacobidet;
				m_jacobi[1][1] = (m_uder[0] * m_wder[2] - m_wder[0] * m_uder[2]) / m_jacobidet;
				m_jacobi[2][1] = (m_vder[0] * m_uder[2] - m_uder[0] * m_vder[2]) / m_jacobidet;

				m_jacobi[0][2] = (m_vder[0] * m_wder[1] - m_wder[0] * m_vder[1]) / m_jacobidet;
				m_jacobi[1][2] = (m_wder[0] * m_uder[1] - m_uder[0] * m_wder[1]) / m_jacobidet;
				m_jacobi[2][2] = (m_uder[0] * m_vder[1] - m_vder[0] * m_uder[1]) / m_jacobidet;*/
				jacobidet.push_back(m_jacobidet);
				jacobi.push_back(m_jacobi);
			}
		}
	}
}

double BSplineSolid::ExactSol(double x, double y, double z)
{
	return x * (1 - x) * y * (1 - y) * (1 + z) * (1 - z);
	//return x * x + y * y + z * z;
}

double BSplineSolid::RHS(double x, double y, double z)
{
	//return 2 * x*y*(x - 1)*(y - 1) + 2 * x*(x - 1)*(z - 1)*(z + 1) + 2 * y*(y - 1)*(z - 1)*(z + 1);
	return 1;
	//return -8;
}

double BSplineSolid::RHS1(double u, double v, double w)
{
	double x;
	x = sqrt((u - 0.3) * (u - 0.3) + (v - 0.3) * (v - 0.3) + (w - 0.3) * (w - 0.3));
	return (1 - x) * 10;
}

void BSplineSolid::JacobiMatrixFace(int solid_id, int uvw, int uvw_id, std::vector<double> u, std::vector<double> v,
	std::vector<double>& jacobidet, std::vector<Matrix>& jacobi)
{
	int m_umaxindex = Local_Node[solid_id][0].size() - 1;
	int m_vmaxindex = Local_Node[solid_id][1].size() - 1;
	int m_uptnum = m_umaxindex - Local_Order[solid_id][0];
	int m_vptnum = m_vmaxindex - Local_Order[solid_id][1];
	double m_ubasevalue[2], m_vbasevalue[2];

	jacobidet.clear();
	jacobi.clear();

	Matrix m_jacobi;
	double m_jacobidet;

	m_jacobi.resize(2);
	m_jacobi[0].resize(2); m_jacobi[1].resize(2);

	int uvw_size;
	if (uvw_id == 0)
		uvw_size = 0;
	else if (uvw_id == 1)
		uvw_size = m_uptnum - 1;

	for (int l = 0; l < v.size(); l++)
	{
		for (int k = 0; k < u.size(); k++)
		{
			double m_value[3] = { 0,0,0 };
			double m_uder[3] = { 0,0,0 };
			double m_vder[3] = { 0,0,0 };
			int m_uIndex = X_LOCATION(solid_id, 0, u[k]);
			int m_vIndex = X_LOCATION(solid_id, 1, v[l]);
			for (int i = m_uIndex - Local_Order[solid_id][0]; i <= m_uIndex; i++)
			{
				DersOneBasisFun(solid_id, 0, i, u[k], 1, m_ubasevalue);
				for (int j = m_vIndex - Local_Order[solid_id][1]; j <= m_vIndex; j++)
				{
					DersOneBasisFun(solid_id, 1, j, v[l], 1, m_vbasevalue);

					V3f vertex;
					if (uvw == U)
						vertex = Local_BSplineSolid[solid_id].vertices_[Control_Point_id[solid_id][j][i][uvw_size]];
					else if (uvw == V)
						vertex = Local_BSplineSolid[solid_id].vertices_[Control_Point_id[solid_id][j][uvw_size][i]];
					else if (uvw == W)
						vertex = Local_BSplineSolid[solid_id].vertices_[Control_Point_id[solid_id][uvw_size][j][i]];

					m_value[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[0];
					m_value[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[0];
					m_value[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[0];

					m_uder[0] += vertex.x * m_ubasevalue[1] * m_vbasevalue[0];
					m_uder[1] += vertex.y * m_ubasevalue[1] * m_vbasevalue[0];
					m_uder[2] += vertex.z * m_ubasevalue[1] * m_vbasevalue[0];

					m_vder[0] += vertex.x * m_ubasevalue[0] * m_vbasevalue[1];
					m_vder[1] += vertex.y * m_ubasevalue[0] * m_vbasevalue[1];
					m_vder[2] += vertex.z * m_ubasevalue[0] * m_vbasevalue[1];
				}
			}
			m_jacobidet = m_uder[0] * m_vder[1] - m_vder[0] * m_uder[1];
			m_jacobi[0][0] = m_vder[1] / m_jacobidet; m_jacobi[0][1] = -m_uder[1] / m_jacobidet;
			m_jacobi[1][0] = -m_vder[0] / m_jacobidet; m_jacobi[1][1] = m_uder[0] / m_jacobidet;
			jacobidet.push_back(m_jacobidet);
			jacobi.push_back(m_jacobi);
		}
	}
}
