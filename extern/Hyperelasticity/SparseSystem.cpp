#include "SparseSystem.h"

const double maxNum = 1e10;

void SparseSystem::init_geometry(const std::vector<CBSplineSolid>& geometry, int dims)
{
    dim = dims;
    m_globalIndex.resize(geometry.size());
    computerIndex(geometry);
    m_gl_coefs_index.resize(controlpoints);
    int count = 0;
    for (int p = 0; p < m_globalIndex.size(); ++p) {
        const vector<CPoint>& coefs = geometry[p].coefs();
        coefsNum += coefs.size();
        for (size_t m = 0; m < coefs.size(); ++m)
        {
            m_gl_coefs_index[m_globalIndex[p][m]].push_back(count);   // ȫ�� �� ���Ƶ���� ������
            assert(count == m + p * coefs.size());
            count++;
        }
    }
    m_globalDirichlet.resize(dim * controlpoints, false);
    m_globalDirichlet_count.resize(dim * controlpoints, false);
    m_dirichletValue.resize(dim * controlpoints, 0);

    m_bigmatrix.resize(dim * controlpoints);
    m_colids.resize(dim * controlpoints);
    m_rhs.resize(dim * controlpoints,1);
    m_matrix.resize(dim * controlpoints, dim * controlpoints);

    createStore(geometry);

    int num = 0;
    for (int i = 0; i < geometry.size();++i) {
        num += geometry[i].getSizeControlPoints();
        m_cpoints_inclusive_sum.push_back(num);
    }
}

void SparseSystem::init_dirichlet(const std::vector<boundary_point>& dirichletPoints)
{
    auto it = dirichletPoints.begin();
    for (; it != dirichletPoints.end(); ++it) {
        assert(it->m_condition_type == 0);
        uint32_t global_idx = it->m_cpoint_global_id;

        m_globalDirichlet[global_idx + controlpoints * 0] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 0] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 0] = it->m_val.getX();

        m_globalDirichlet[global_idx + controlpoints * 1] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 1] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 1] = it->m_val.getY();

        m_globalDirichlet[global_idx + controlpoints * 2] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 2] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 2] = it->m_val.getZ();
    }
}

SparseSystem::SparseSystem(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_condition>& dirichletSides, int dims)
{
    dim = dims;
    m_globalIndex.resize(geometry.size());
    computerIndex(geometry);
    m_gl_coefs_index.resize(controlpoints);
    int count = 0;
    for (int p = 0; p < m_globalIndex.size(); ++p) {
        const vector<CPoint>& coefs = geometry[p].coefs();
        coefsNum += coefs.size();
        for (size_t m = 0; m < coefs.size(); ++m)
        {
            m_gl_coefs_index[m_globalIndex[p][m]].push_back(count);   // ȫ�� �� ���Ƶ���� ������
            assert(count == m + p * coefs.size());
            count++;
        }
    }
    m_globalDirichlet.resize(dim * controlpoints, false);
    m_globalDirichlet_count.resize(dim * controlpoints, false);
    m_dirichletValue.resize(dim * controlpoints, 0);
    setDirichlet(geometry, dirichletSides);
    m_bigmatrix.resize(dim * controlpoints);
    m_colids.resize(dim * controlpoints);
    m_rhs.resize(dim * controlpoints,1);
    m_matrix.resize(dim * controlpoints, dim * controlpoints);
    createStore(geometry);
    computerInterface(geometry);

    int num = 0;
    for (int i = 0; i < geometry.size();++i) {
        num += geometry[i].getSizeControlPoints();
        m_cpoints_inclusive_sum.push_back(num);
    }
}

SparseSystem::SparseSystem(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_point>& dirichletPoints, int dims)
{
    dim = dims;
    m_globalIndex.resize(geometry.size());
    computerIndex(geometry);
    m_gl_coefs_index.resize(controlpoints);
    int count = 0;
    for (int p = 0; p < m_globalIndex.size(); ++p) {
        const vector<CPoint>& coefs = geometry[p].coefs();
        coefsNum += coefs.size();
        for (size_t m = 0; m < coefs.size(); ++m)
        {
            m_gl_coefs_index[m_globalIndex[p][m]].push_back(count);   //  p = count / 64, m = count - 64 * p;
            assert(count == m + p * coefs.size());
            count++;
        }
    }
    m_globalDirichlet.resize(dim * controlpoints, false);
    m_globalDirichlet_count.resize(dim * controlpoints, false);
    m_dirichletValue.resize(dim * controlpoints, 0);
    setDirichlet(geometry, dirichletPoints);
    m_bigmatrix.resize(dim * controlpoints);
    m_colids.resize(dim * controlpoints);
    m_rhs.resize(dim * controlpoints,1);
    m_matrix.resize(dim * controlpoints, dim * controlpoints);
    createStore(geometry);
    computerInterface(geometry);

    int num = 0;
    for (int i = 0; i < geometry.size();++i) {
        num += geometry[i].getSizeControlPoints();
        m_cpoints_inclusive_sum.push_back(num);
    }
}


void SparseSystem::computerInterface(const std::vector<CBSplineSolid>& geometry)
{
    std::vector<CFace> allface(geometry.size()*6);
    for (size_t k = 0; k < geometry.size(); ++k)
    {
        int usize = geometry[k].getSizeUControlPoints();
        int vsize = geometry[k].getSizeVControlPoints();
        int wsize = geometry[k].getSizeWControlPoints();
        auto coefs = geometry[k].getControlPoints();
        allface[k * 6].p = k; allface[k * 6].dir = 1;
        allface[k * 6+1].p = k; allface[k * 6+1].dir = 2;
        allface[k * 6+2].p = k; allface[k * 6+2].dir = 3;
        allface[k * 6+3].p = k; allface[k * 6+3].dir = 4;
        allface[k * 6+4].p = k; allface[k * 6+4].dir = 5;
        allface[k * 6+5].p = k; allface[k * 6+5].dir = 6;
        for (int i = 0; i < usize; ++i)
        {
            for (int j = 0; j < vsize; ++j)
            {
                for (int w = 0; w < wsize; ++w)
                {
                    LocInd temp;
                    temp.p = i; temp.x = coefs[i][j][w].getX();
                    temp.y = coefs[i][j][w].getY(); temp.z = coefs[i][j][w].getZ();

                    if (i == 0) allface[k * 6].locInfo.emplace_back(temp);
                    if (i == usize - 1) allface[k * 6+1].locInfo.emplace_back(temp);
                    if (j == 0) allface[k * 6+2].locInfo.emplace_back(temp);
                    if (j == vsize - 1) allface[k * 6+3].locInfo.emplace_back(temp);
                    if (w == 0) allface[k * 6+4].locInfo.emplace_back(temp);
                    if (w == wsize - 1) allface[k * 6+5].locInfo.emplace_back(temp);

                }
            }
        }
        std::sort(allface[k * 6].locInfo.begin(), allface[k * 6].locInfo.end());
        std::sort(allface[k * 6+1].locInfo.begin(), allface[k * 6+1].locInfo.end());
        std::sort(allface[k * 6+2].locInfo.begin(), allface[k * 6+2].locInfo.end());
        std::sort(allface[k * 6+3].locInfo.begin(), allface[k * 6+3].locInfo.end());
        std::sort(allface[k * 6+4].locInfo.begin(), allface[k * 6+4].locInfo.end());
        std::sort(allface[k * 6+5].locInfo.begin(), allface[k * 6+5].locInfo.end());
    }
    m_inter.resize(geometry.size(), std::vector<pair<int,int>>(6,{-1,-1}));
    std::sort(allface.begin(), allface.end());
    for (int i = 0; i < allface.size()-1; ++i)
    {
        if (allface[i] == allface[i + 1]) {
            m_inter[allface[i].p][allface[i].dir - 1] = { allface[i + 1].p,allface[i + 1].dir };
            m_inter[allface[i+1].p][allface[i+1].dir - 1] = { allface[i].p,allface[i].dir };
        }
    }
}

void SparseSystem::computerIndex(const std::vector<CBSplineSolid>& geometry)
{
    std::vector<LocInd> indexSelect;
    for (size_t i = 0; i < geometry.size(); ++i)
    {
        auto coefs = geometry[i].coefs();
        m_globalIndex[i].resize((int)coefs.size());
        int index = 0;
        for (CPoint& it : coefs) {
            LocInd temp;
            temp.p = i; temp.x = it.getX();
            temp.y = it.getY(); temp.z = it.getZ();
            temp.index = index++;
            indexSelect.emplace_back(temp);
        }
    }
    m_ConnectPatch.resize(geometry.size());
    std::sort(indexSelect.begin(), indexSelect.end());
    int p = indexSelect[0].p, loc = indexSelect[0].index;
    m_globalIndex[p][loc] = controlpoints++;
    for (size_t i = 1; i < indexSelect.size(); ++i) {
        p = indexSelect[i].p; loc = indexSelect[i].index;
        if (indexSelect[i] == indexSelect[i - 1]) {
            if (indexSelect[i - 1].p != p) {
                m_ConnectPatch[indexSelect[i - 1].p].insert(p);
                m_ConnectPatch[p].insert(indexSelect[i - 1].p);
            }
            m_globalIndex[p][loc] = controlpoints - 1;
        }
        else
            m_globalIndex[p][loc] = controlpoints++;
    }

    //for debug
    if(geometry.size() == 1)
    {
        for (int i = 0; i < 64; ++i)
        {
            m_globalIndex[0][i] = i;
        }
    }
}

void SparseSystem::setDirichlet(const std::vector<CBSplineSolid>& geometry, const std::vector<boundary_condition>& dirichletSides)
{

    auto it = dirichletSides.begin();
    for (; it != dirichletSides.end(); ++it) {
        int patchIndex = it->patchSide();
        boxSide side = it->side();
        const CBSplineSolid* geo = &geometry[patchIndex];
        int m = geo->getSizeUControlPoints(), n = geo->getSizeVControlPoints(),
                o = geo->getSizeWControlPoints();
        if (side.index() == boundary::east || side.index() == boundary::west) {
            int u = side.index() == boundary::east ? 1 : 0;
            for (int k = 0; k < o; ++k) {
                for (int j = 0; j < n; ++j) {
                    int localIndex = k * m * n + j * m + u * (m - 1);
                    m_globalDirichlet[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = true;
                    m_globalDirichlet_count[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()]++;
                    m_dirichletValue[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
                }
            }
        }
        else if (side.index() == boundary::south || side.index() == boundary::north) {
            int v = side.index() == boundary::north ? 1 : 0;
            for (int k = 0; k < o; ++k) {
                for (int i = 0; i < m; ++i) {
                    int localIndex = k * m * n + i + v * ((n - 1) * m);
                    m_globalDirichlet[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = true;
                    m_globalDirichlet_count[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()]++;
                    m_dirichletValue[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
                }
            }
        }
        else {
            int w = side.index() == boundary::back ? 1 : 0;
            for (int j = 0; j < n; ++j) {
                for (int i = 0; i < m; ++i) {
                    int localIndex = j * m + i + w * ((o - 1) * m * n);
                    m_globalDirichlet[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = true;
                    m_globalDirichlet_count[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()]++;
                    m_dirichletValue[m_globalIndex[patchIndex][localIndex] + controlpoints * it->getDir()] = it->getValue();
                }
            }
        }
    }
    //m_globalDirichlet[m_globalIndex[0][20] + controlpoints * 2] = true;
    //m_globalDirichlet[m_globalIndex[0][20] + controlpoints] = true;
    //m_globalDirichlet[m_globalIndex[0][20]] = true;
    //m_dirichletValue[m_globalIndex[0][20] + controlpoints * 2] = 0.5;
    //m_dirichletValue[m_globalIndex[0][20] + controlpoints] = 0;
    //m_dirichletValue[m_globalIndex[0][20]] = 0;
}

void SparseSystem::setDirichlet(const vector<CBSplineSolid>&  geometry,
                                const vector<boundary_point>& dirichletPoints)
{

    for (auto it = dirichletPoints.begin(); it != dirichletPoints.end(); ++it) {
        assert(it->m_condition_type == 0);
        uint32_t global_idx = it->m_cpoint_global_id;

        m_globalDirichlet[global_idx + controlpoints * 0] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 0] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 0] = it->m_val.getX();

        m_globalDirichlet[global_idx + controlpoints * 1] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 1] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 1] = it->m_val.getY();

        m_globalDirichlet[global_idx + controlpoints * 2] = true;
        m_globalDirichlet_count[global_idx + controlpoints * 2] = get_g2l_map()[global_idx].size();
        m_dirichletValue[global_idx + controlpoints * 2] = it->m_val.getZ();
    }
}


void SparseSystem::createStore(const std::vector<CBSplineSolid>& geometry)
{
    for (int i = 0; i < geometry.size(); ++i) {

        const CBSplineSolid* geo = &geometry[i];
        std::vector<double> uknots = geo->getUKnots(),
                vknots = geo->getVKnots(), wknots = geo->getWKnots();
        int gaussNum_u = geo->getUDegree() + 1, gaussNum_v = geo->getVDegree() + 1,
                gaussNum_w = geo->getWDegree() + 1;
        int gaussNum = gaussNum_u * gaussNum_v * gaussNum_w;
        std::vector<int> localIndex(gaussNum);
        int N_D;
        for (int w_interval = 0; w_interval < wknots.size() - 1; ++w_interval)  // ��Ԫu����ѭ��
        {
            for (int v_interval = 0; v_interval < vknots.size() - 1; ++v_interval) // ��Ԫv����ѭ��
            {
                for (int u_interval = 0; u_interval < uknots.size() - 1; ++u_interval)
                {

                    geo->active(u_interval, v_interval, w_interval, localIndex);
                    N_D = localIndex.size();
                    for (int m = 0; m < N_D; ++m) {
                        for (int n = 0; n < N_D; ++n) {
                            for (int di = 0; di < dim; di++) {
                                for (int dj = 0; dj < dim; dj++) {

                                    int r = m_globalIndex[i][localIndex[m]] + di * controlpoints,
                                            c = m_globalIndex[i][localIndex[n]] + dj * controlpoints;
                                    //if (r != c && (m_globalDirichlet[r] || m_globalDirichlet[c])) continue;

                                    if (!m_colids[r].count(c)) {
                                        m_colids[r][c] = m_bigmatrix[r].size();
                                        m_bigmatrix[r].push_back(0.0);
                                    }
                                }
                            }
                        }
                    }
                }
            }
        }
    }
}

void SparseSystem::constructCoefs(const Eigen::VectorXd& solution, Eigen::MatrixXd& displacement)
{
    displacement.setZero(coefsNum, dim);
    for (int i = 0; i < controlpoints; ++i) {
        for (auto index : m_gl_coefs_index[i]) {
            for (int di = 0; di < dim; ++di) {
                displacement(index, di) = solution(i + di * controlpoints);
            }
        }
    }
}

void SparseSystem::pushToMatrix(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localMats)
{

    for (int i = 0; i < 192; ++i)
    {
        for (int j = 0; j < 192; ++j)
        {
            if (abs(localMats(i, j) - localMats(j, i)) > 0.0001 * abs(localMats(i, j)) && abs(localMats(i, j) - localMats(j, i)) > 0.00001)
            {
                std:: cout << "noequal for i-" << i << ", j-" << j << " with A(i,j) =" << localMats(i, j) << "  , A(j, i) =" << localMats(j, i) << std::endl;
                std::cout << "difference is " << abs(localMats(i, j) - localMats(j, i)) << endl;
            }
        }
    }

    int N_D = localIndex.size();
    for (int m = 0; m < N_D; ++m) {
        for (int n = 0; n < N_D; ++n) {
            for (int di = 0; di < dim; di++) {
                for (int dj = 0; dj < dim; dj++) {
                    //if (fabs(localMats(N_D * di + m, N_D * dj + n)) < tol) continue;
                    int r = m_globalIndex[patch][localIndex[m]] + di * controlpoints,
                            c = m_globalIndex[patch][localIndex[n]] + dj * controlpoints;
                    //if (m_globalDirichlet[r] || m_globalDirichlet[c] ) continue;

                    if (!m_colids[r].count(c)) {
                        cout << "Index error" << endl;
                        exit(1);
                    }
                    m_bigmatrix[r][m_colids[r][c]] += localMats(N_D * di + m, N_D * dj + n);
                }
            }
        }
    }
}

void SparseSystem::pushToMatrix(uint32_t patch_id, void* data)
{
    double* array = reinterpret_cast<double*>(data);

    for(uint16_t m = 0; m < 64; ++m)
    {
        for(uint16_t n = 0; n < 64; ++n)
        {
            for(uint16_t di = 0; di < 3; ++di)
            {
                for(uint16_t dj = 0; dj < 3; ++dj)
                {
                    auto r = m_globalIndex[patch_id][m] + di * controlpoints;
                    auto c = m_globalIndex[patch_id][n] + dj * controlpoints;

                    if (!m_colids[r].count(c)) {
                        cout << "Index error" << endl;
                        exit(1);
                    }

                    m_bigmatrix[r][m_colids[r][c]] += array[(m + 64 * di) + 192 * (n + 64 * dj)];
                }
            }
        }
    }
}

void SparseSystem::pushToRhs(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localRhs)
{
    int N_D = localIndex.size();
    for (int m = 0; m < N_D; ++m) {
        for (int di = 0; di < dim; di++) {
            int r = m_globalIndex[patch][localIndex[m]] + di * controlpoints;
            //if (m_globalDirichlet[r] == 1) continue;
            m_rhs(m_globalIndex[patch][localIndex[m]] + di * controlpoints, 0) += localRhs(N_D * di + m, 0);
        }
    }
}

void SparseSystem::pushToRhs(const std::vector<pair<int, double>>& f)
{
    for (int i = 0; i < f.size(); ++i)
    {
        //if (m_globalDirichlet[f[i].first] == 1) continue;
        m_rhs(f[i].first, 0) += f[i].second;
    }
}



void SparseSystem::pushToSystem(int patch, const std::vector<int>& localIndex, const Eigen::MatrixXd& localRhs, const Eigen::MatrixXd& localMats, bool isAssembleMatrix)
{
    if(isAssembleMatrix)
        pushToMatrix(patch, localIndex, localMats);
    pushToRhs(patch, localIndex, localRhs);
}

void SparseSystem::dataPush(bool isOne)
{
    for (int i = 0; i < m_globalDirichlet.size(); ++i) {
        if (m_globalDirichlet[i]) {
            //m_rhs(i, 0) = m_dirichletValue[i];

            //if(isOne)
            //	m_bigmatrix[i][m_colids[i][i]] = 1.0;
            //else
            //	m_bigmatrix[i][m_colids[i][i]] = 0.0;

            //test cancel Dirichlet Boundary Condition

            bool bCGTest = false;

            if (bCGTest)
            {
                m_rhs(i, 0) += m_dirichletValue[i] * maxNum * m_globalDirichlet_count[i];
                m_bigmatrix[i][m_colids[i][i]] += maxNum * m_globalDirichlet_count[i];
            }
        }
    }



    //clear m_bidmatrix data when dataPush finish
    std::vector< Eigen::Triplet<double> > tripletList;
    for (int begin = 0; begin < m_colids.size(); ++begin) {
        for (auto& col : m_colids[begin]) {
            tripletList.emplace_back(Eigen::Triplet<double>(begin, col.first, m_bigmatrix[begin][col.second]));
            m_bigmatrix[begin][col.second] = 0.0;
        }
    }

    m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    m_matrix.makeCompressed();
}

void SparseSystem::dataPush_Dirichlet(bool isOne)
{
    for (int i = 0; i < m_globalDirichlet.size(); ++i) {
        if (m_globalDirichlet[i]) {
            //m_rhs(i, 0) = m_dirichletValue[i];

            //if(isOne)
            //	m_bigmatrix[i][m_colids[i][i]] = 1.0;
            //else
            //	m_bigmatrix[i][m_colids[i][i]] = 0.0;

            //test cancel Dirichlet Boundary Condition

            bool bCGTest = true;

            if (bCGTest)
            {
                m_rhs(i, 0) += m_dirichletValue[i] * maxNum * m_globalDirichlet_count[i];
                m_bigmatrix[i][m_colids[i][i]] += maxNum * m_globalDirichlet_count[i];
            }
        }
    }

    std::vector< Eigen::Triplet<double> > tripletList;
    for (int begin = 0; begin < m_colids.size(); ++begin) {
        for (auto& col : m_colids[begin]) {
            tripletList.emplace_back(Eigen::Triplet<double>(begin, col.first, m_bigmatrix[begin][col.second]));
            m_bigmatrix[begin][col.second] = 0.0;
        }
    }

    m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    m_matrix.makeCompressed();
}


void SparseSystem::dataPush_Dirichlet_V2() {
    for (int i = 0; i < m_globalDirichlet.size(); ++i) {
        if (m_globalDirichlet[i]) {

            m_rhs(i, 0) += m_dirichletValue[i] * m_globalDirichlet_count[i];

            for(auto& j : m_colids[i])
            {
                //j.second = 列号
                m_bigmatrix[i][j.second] = 0;
            }
            m_bigmatrix[i][m_colids[i][i]] += m_globalDirichlet_count[i];
        }
    }

    std::vector< Eigen::Triplet<double> > tripletList;
    for (int begin = 0; begin < m_colids.size(); ++begin) {
        for (auto& col : m_colids[begin]) {
            tripletList.emplace_back(Eigen::Triplet<double>(begin, col.first, m_bigmatrix[begin][col.second]));
            m_bigmatrix[begin][col.second] = 0.0;
        }
    }

    m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    m_matrix.makeCompressed();
}




void SparseSystem::dataPush(Eigen::MatrixXd& extraRhs)
{
    for (int i = 0; i < m_globalDirichlet.size(); ++i) {
        if (m_globalDirichlet[i]) {
            m_rhs(i, 0) += m_dirichletValue[i] * maxNum;
            m_bigmatrix[i][m_colids[i][i]] += maxNum;
        }
    }
    std::cout << "--------------" << std::endl;
    for (int i = 0; i < extraRhs.rows(); ++i)
    {
        if (fabs(extraRhs(i, 0)) > 1e-8 && !m_globalDirichlet[i]) {
            //if (i >= controlpoints && i < 2*controlpoints) {
            //	for (int j = 0; j < m_gl_coefs_index[i - controlpoints].size(); ++j)
            //		std::cout << m_gl_coefs_index[i - controlpoints][j] / 64 << " " << m_gl_coefs_index[i - controlpoints][j] % 64 << " ";
            //	std::cout << extraRhs(i, 0) << std::endl;
            //}
            m_rhs(i, 0) += extraRhs(i, 0) * maxNum;
            m_bigmatrix[i][m_colids[i][i]] += maxNum;
        }
    }



    std::vector< Eigen::Triplet<double> > tripletList;
    for (int begin = 0; begin < m_colids.size(); ++begin) {
        for (auto& col : m_colids[begin]) {
            tripletList.emplace_back(Eigen::Triplet<double>(begin, col.first, m_bigmatrix[begin][col.second]));
            m_bigmatrix[begin][col.second] = 0.0;
        }
    }

    m_matrix.setFromTriplets(tripletList.begin(), tripletList.end());
    m_matrix.makeCompressed();
}

int SparseSystem::getIndexById(int patch, int index)
{
    return m_globalIndex[patch][index];
}

bool SparseSystem::IsDirichlet(int patch, int index)
{
    return m_globalDirichlet[getIndexById(patch, index)];
}

bool  SparseSystem::IsDirichlet(uint32_t global_id)
{
    return m_globalDirichlet[global_id];
}

void SparseSystem::constructNewMatrix(Eigen::VectorXd& extraRhs, Eigen::SparseMatrix<double>& newMatrix, Eigen::MatrixXd& newRhs)
{
    newMatrix.resize(m_rhs.rows(), m_rhs.rows());
    newMatrix.setZero();
    newRhs.setZero(m_rhs.rows(), 1);

    for (int i = 0; i < extraRhs.rows(); ++i)
    {
        if (fabs(extraRhs(i)) > 1e-8 && !m_globalDirichlet[i]) {
            newRhs(i) += extraRhs(i) * 1e12;
            newMatrix.coeffRef(i,i) += 1e12;
        }
    }
}
