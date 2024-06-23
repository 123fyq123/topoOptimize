//
// Created by XieJin on 10/25/2023.
//

#include "IGA_CC_Helper.h"

vector<CBSplineSurface> IGA_CC_Helper::muti_bezier_batch_ = {};
vector<CBSplineSurface> IGA_CC_Helper::muti_bezier_batch_dir_ = {};

void IGA_CC_Helper::write2Bezier(const std::string& _filename, BSplineSolid& m_solids) {
    std::ofstream off(_filename.c_str(), std::ios::out | std::ios::trunc);

    int n = m_solids.Local_V3f.size();
    off << "<xml>" << std::endl;
    for (int i = 0; i < n; ++i) {
        off << "<Geometry type=\"TensorBSpline3\" id=\" " << i << "\">" << std::endl;

        off << "   <Basis type=\"TensorBSplineBasis3\" parDim=\"3\">" << std::endl;
        off << "    <Basis type = \"BSplineBasis\" index = \"0\">" << std::endl;
        off << "     <KnotVector degree = \"3\">0.00000   0.00000   0.00000   0.00000    1.00000   1.00000   1.00000    1.00000 </KnotVector>" << std::endl;
        off << "   </Basis>" << std::endl;
        off << "    <Basis type = \"BSplineBasis\" index = \"1\">" << std::endl;
        off << "     <KnotVector degree = \"3\">0.00000   0.00000   0.00000   0.00000    1.00000   1.00000   1.00000    1.00000 </KnotVector>" << std::endl;
        off << "   </Basis>" << std::endl;
        off << "    <Basis type = \"BSplineBasis\" index = \"2\">" << std::endl;
        off << "     <KnotVector degree = \"3\">0.00000   0.00000   0.00000   0.00000    1.00000   1.00000   1.00000    1.00000 </KnotVector>" << std::endl;
        off << "   </Basis>" << std::endl;
        off << "  </Basis>" << std::endl;
        off << " <coefs geoDim = \"3\">";
        for (const auto& f : m_solids.Control_Point[i]) {
            for (const auto& e : f) {
                for (const auto& v : e) {
                    off << std::endl << std::fixed << m_solids.Local_V3f[i][v].x << " "
                        << m_solids.Local_V3f[i][v].y << " "
                        << m_solids.Local_V3f[i][v].z;
                }
            }
        }
        off << "  </coefs>" << std::endl;
        off << " </Geometry>" << std::endl;
    }
    off << "</xml>" << std::endl;
}

void IGA_CC_Helper::write2bdy(const std::string& _filename, IGA::CCBezierBdy& ccbdy, BSplineSolid& m_solids) {

    std::ofstream off(_filename.c_str(), std::ios::out | std::ios::trunc);
    int n = ccbdy.size();
    off << n << "\n";

    for (int i = 0; i < n; ++i) {
        int c_id = ccbdy[i].cell_id;
        int l_id = ccbdy[i].local_base_id;

        off << m_solids.Local_V3f[c_id][l_id].x << " " << m_solids.Local_V3f[c_id][l_id].y << " " << m_solids.Local_V3f[c_id][l_id].z << " "
            << ccbdy[i].bdy_type << " " << ccbdy[i].bdy_value[0] << " " << ccbdy[i].bdy_value[1] << " " << ccbdy[i].bdy_value[2] << "\n";
    }
}

void IGA_CC_Helper::readUU(const std::string& _filename, std::vector<std::vector<CPoint>>& uu) {
    int n = uu.size();

    std::ifstream in(_filename);

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < 64; ++j) {
            double x, y, z;
            in >> x >> y >> z;
            uu[i][j] = CPoint(x, y, z);
        }
    }
}

void IGA_CC_Helper::write2prot(const std::string& _filename, const std::vector<int>& protect_ids){
    std::ofstream off(_filename.c_str(), std::ios::out | std::ios::trunc);
    size_t n = protect_ids.size();
    off << n << "\n";

    for (auto i = 0; i < n; ++i) {
        off << protect_ids[i] <<" ";
    }
    off <<"\n";
    off.close();
}
void IGA_CC_Helper::write2neighbor(const std::string& _filename,
                                   const std::vector<std::vector<CellHandle>>& aroundcell,
                                   const std::vector<std::vector<double>>& aroundcell_weight){
    std::ofstream off(_filename.c_str(), std::ios::out | std::ios::trunc);
    size_t n = aroundcell.size();
    off << n << "\n";

    for (auto i = 0; i < n; ++i) {
        off << aroundcell[i].size() <<"\n";
        for(int j = 0; i < aroundcell[i].size(); ++j){
            off << aroundcell[i][j] <<" ";
        }
        off << "\n";
        for(int j = 0; i < aroundcell_weight[i].size(); ++j){
            off << aroundcell_weight[i][j] <<" ";
        }
        off << "\n";
    }
    off <<"\n";
    off.close();
}

std::string IGA_CC_Helper::get_model_name(std::string& file_path){
    size_t pos = file_path.find_last_of('/');
    std::string file_name = file_path.substr(pos + 1);
    pos = file_name.find_last_of('.');
    std::string file_prefix = file_name.substr(0, pos);
    std::string file_suffix = file_name.substr(pos);
    return file_prefix;
}


bool IGA_CC_Helper::read2BezierOnly(const std::string& _filename) {
    // 边界面对应数字计算
    int cube[4][4][4] = { 0 };
    int to_face[6][4][4] = { 0 };
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                cube[i][j][k] = i * 16 + j * 4 + k;

    for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
            to_face[0][j][k] = cube[0][j][k];
            to_face[1][j][k] = cube[3][j][k];
            to_face[2][j][k] = cube[j][0][k];
            to_face[3][j][k] = cube[j][3][k];
            to_face[4][j][k] = cube[j][k][0];
            to_face[5][j][k] = cube[j][k][3];
        }
    }

    vector<vector<vector<vector<int>>>> control_point;
    vector<V3f> global_V3f;
    vector<int> global_V3f_count;
    function<int(double, double, double)> addVectex = [&](double x, double y, double z) {
        V3f point(x, y, z);
        for (int i = 0; i < global_V3f.size(); ++i)
        {
            if (distance(point, global_V3f[i]) <= 1e-5) {
                ++global_V3f_count[i];
                return i;
            }
        }
        global_V3f.emplace_back(point);
        global_V3f_count.emplace_back(1);
        return (int)global_V3f.size() - 1;
    };

    std::ifstream iff(_filename.c_str(), std::ios::in);
    if (!iff.good())
    {
        std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        iff.close();
        return false;
    }
    std::string line;
    bool is_number = false;
    int soild_num = 0;
    int control_num = 0;
    vector<vector<vector<int>>> oneSoild(4, vector<vector<int>>(4, vector<int>(4)));
    while (std::getline(iff, line))
    {
        // 逐行读取文件内容
        vector<string> ss;
        if (std::isdigit(line[0]) || line[0] == '-') {
            //std::cout << line[0] << std::endl;
            is_number = true;
            int left = 0;
            for (int i = 0; i < line.length(); ++i) {
                if (line[i] == ' ') {
                    ss.emplace_back(line.substr(left, i - left));
                    left = i + 1;
                }
            }
            ss.emplace_back(line.substr(left, line.length() - left));
            cout << "point " << control_num << endl;
            double value[3] = { 0 };
            for (int i = 0; i < ss.size(); i++) {
                if (ss[i].length() == 0) continue;
                if (std::isdigit(ss[i][0]) || ss[i][0] == '-') {
                    value[i] = std::stod(ss[i]);
                    cout << value[i] << " ";
                }
            }
            int vid = addVectex(value[0], value[1], value[2]);
            int level_num = control_num / 16;
            int row_num = (control_num % 16) / 4;
            int column_num = control_num % 4;
            oneSoild[level_num][row_num][column_num] = vid;
            std::cout << std::endl;
            control_num++;
        }
        else {
            if (is_number) {
                ++soild_num;
                control_point.emplace_back(oneSoild);
                control_num = 0;
                cout << soild_num << endl;
                is_number = false;
            }
        }
    }

    vector<vector<vector<Vec>>> faces_control_points;
    vector<vector<vector<Vec>>> faces_control_points_dir;
    // 边界面
    for (int soild_i = 0; soild_i < soild_num; ++soild_i)
    {
        for (int face_i = 0; face_i < 6; ++face_i)
        {
            vector<vector<Vec>> oneface_control_points(4, vector<Vec>(4));
            vector<vector<Vec>> oneface_control_points_dir(4, vector<Vec>(4));
            bool is_boundary = false;
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    int control_num = to_face[face_i][i][j];
                    int level_num = control_num / 16;
                    int row_num = (control_num % 16) / 4;
                    int column_num = control_num % 4;
                    int vid = control_point[soild_i][level_num][row_num][column_num];
                    oneface_control_points[i][j] = Vec(global_V3f[vid].x, global_V3f[vid].y, global_V3f[vid].z);
                    if (global_V3f_count[vid] == 1)
                        is_boundary = true;
                }
            }
            if (is_boundary) {
                faces_control_points.emplace_back(oneface_control_points);
            }

        }
    }
    cad_restruction(faces_control_points);
}

bool IGA_CC_Helper::read2Bezier(const std::string& _filename) {
    // 边界面对应数字计算
    int cube[4][4][4] = { 0 };
    int to_face[6][4][4] = { 0 };
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            for (int k = 0; k < 4; ++k)
                cube[i][j][k] = i * 16 + j * 4 + k;

    for (int j = 0; j < 4; ++j) {
        for (int k = 0; k < 4; ++k) {
            to_face[0][j][k] = cube[0][j][k];
            to_face[1][j][k] = cube[3][j][k];
            to_face[2][j][k] = cube[j][0][k];
            to_face[3][j][k] = cube[j][3][k];
            to_face[4][j][k] = cube[j][k][0];
            to_face[5][j][k] = cube[j][k][3];
        }
    }

    vector<vector<vector<vector<int>>>> control_point;
    vector<V3f> global_V3f;
    vector<int> global_V3f_count;
    function<int(double, double, double)> addVectex = [&](double x, double y, double z) {
        V3f point(x, y, z);
        for (int i = 0; i < global_V3f.size(); ++i)
        {
            if (distance(point, global_V3f[i]) <= 1e-5) {
                ++global_V3f_count[i];
                return i;
            }
        }
        global_V3f.emplace_back(point);
        global_V3f_count.emplace_back(1);
        return (int)global_V3f.size() - 1;
    };

    std::ifstream iff(_filename.c_str(), std::ios::in);
    if (!iff.good())
    {
        std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        iff.close();
        return false;
    }
    std::string line;
    bool is_number = false;
    int soild_num = 0;
    int control_num = 0;
    vector<vector<vector<int>>> oneSoild(4, vector<vector<int>>(4, vector<int>(4)));
    while (std::getline(iff, line))
    {
        // 逐行读取文件内容
        vector<string> ss;
        if (std::isdigit(line[0]) || line[0] == '-') {
            //std::cout << line[0] << std::endl;
            is_number = true;
            int left = 0;
            for (int i = 0; i < line.length(); ++i){
                if (line[i] == ' ') {
                    ss.emplace_back(line.substr(left, i - left));
                    left = i + 1;
                }
            }
            ss.emplace_back(line.substr(left, line.length() - left));
            cout << "point " << control_num << endl;
            double value[3] = {0};
            for (int i = 0; i < ss.size(); i++) {
                if (ss[i].length() == 0) continue;
                if (std::isdigit(ss[i][0]) || ss[i][0] == '-') {
                    value[i] = std::stod(ss[i]);
                    cout << value[i] << " ";
                }
            }
            int vid = addVectex(value[0], value[1], value[2]);
            int level_num = control_num / 16;
            int row_num = (control_num % 16) / 4;
            int column_num = control_num % 4;
            oneSoild[level_num][row_num][column_num] = vid;
            std::cout << std::endl;
            control_num++;
        }
        else {
            if (is_number) {
                ++soild_num;
                control_point.emplace_back(oneSoild);
                control_num = 0;
                cout << soild_num << endl;
                is_number = false;
            }
        }
    }

    vector<Vec> global_V3f_dir(global_V3f.size());
    std::string filename = "C:\\Users\\Administrator\\Desktop\\long beam_skel_fithex_contorl.u";
    std::ifstream in(filename);

    double min_x = DBL_MAX;
    double max_x = DBL_MIN;
    double min_y = DBL_MAX;
    double max_y = DBL_MIN;
    double min_z = DBL_MAX;
    double max_z = DBL_MIN;
    for (int i = 0; i < soild_num; ++i) {
        for (int control_num = 0; control_num < 64; ++control_num) {
            int level_num = control_num / 16;
            int row_num = (control_num % 16) / 4;
            int column_num = control_num % 4;
            int vid = control_point[i][level_num][row_num][column_num];
            double x, y, z;
            in >> x >> y >> z;
            min_x = min(x, min_x);
            min_y = min(y, min_y);
            min_z = min(z, min_z);
            max_x = max(x, max_x);
            max_y = max(y, max_y);
            max_z = max(z, max_z);
            global_V3f_dir[vid] = Vec(x, y, z);
        }
    }

    double range_x = max_x - min_x;
    cout << "x: " << min_x << " - " << max_x << endl;
    double range_y = max_y - min_y;
    cout << "y: " << min_y << " - " << max_y << endl;
    double range_z = max_z - min_z;
    cout << "z: " << min_z << " - " << max_z << endl;
    for (int i = 0; i < global_V3f_dir.size(); ++i)
    {
        global_V3f_dir[i].setX( (global_V3f_dir[i].x() - min_x) / range_x);
        global_V3f_dir[i].setY((global_V3f_dir[i].y() - min_y) / range_y);
        global_V3f_dir[i].setZ((global_V3f_dir[i].z() - min_z) / range_z);
    }

    vector<vector<vector<Vec>>> faces_control_points;
    vector<vector<vector<Vec>>> faces_control_points_dir;
    // 边界面
    for (int soild_i = 0; soild_i < soild_num; ++soild_i)
    {
        for (int face_i = 0; face_i < 6; ++face_i)
        {
            vector<vector<Vec>> oneface_control_points(4, vector<Vec>(4));
            vector<vector<Vec>> oneface_control_points_dir(4, vector<Vec>(4));
            bool is_boundary = false;
            for (int i = 0; i < 4; ++i)
            {
                for (int j = 0; j < 4; ++j)
                {
                    int control_num = to_face[face_i][i][j];
                    int level_num = control_num / 16;
                    int row_num = (control_num % 16) / 4;
                    int column_num = control_num % 4;
                    int vid = control_point[soild_i][level_num][row_num][column_num];
                    oneface_control_points[i][j] = Vec(global_V3f[vid].x, global_V3f[vid].y, global_V3f[vid].z);

                    oneface_control_points_dir[i][j] = global_V3f_dir[vid];
                    if (global_V3f_count[vid] == 1)
                        is_boundary = true;
                }
            }
            if (is_boundary) {
                faces_control_points.emplace_back(oneface_control_points);
                faces_control_points_dir.emplace_back(oneface_control_points_dir);
            }

        }
    }
    cad_restruction(faces_control_points, faces_control_points_dir);
}

bool IGA_CC_Helper::cad_restruction(vector<vector<vector<Vec>>>& control_points, vector<vector<vector<Vec>>>& control_points_dir) {
    int udegree = 3;
    int vdegree = 3;
    vector<double> uknots{ 0,1 };
    vector<double> vknots{ 0,1 };
    vector<int> umultis{ 4,4 };
    vector<int> vmultis{ 4,4 };

    muti_bezier_batch_.clear();
    muti_bezier_batch_dir_.clear();
    for (int i = 0; i < control_points.size(); i++) {
        const auto& single_patch = control_points[i];
        vector<vector<CPoint>> control_surface;
        for (int j = 0; j < 4; ++j) {
            vector<CPoint> control_line;
            for (int k = 0; k < 4; ++k) {
                control_line.push_back(CPoint(single_patch[j][k].x(), single_patch[j][k].y(), single_patch[j][k].z()));
            }
            control_surface.push_back(control_line);
        }
        CBSplineSurface bsplineSurface(udegree, vdegree, control_surface, uknots, vknots, umultis, vmultis);
        muti_bezier_batch_.emplace_back(bsplineSurface);
    }
    for (int i = 0; i < control_points_dir.size(); i++) {
        const auto& single_patch = control_points_dir[i];
        vector<vector<CPoint>> control_surface;
        for (int j = 0; j < 4; ++j) {
            vector<CPoint> control_line;
            for (int k = 0; k < 4; ++k) {
                control_line.push_back(CPoint(single_patch[j][k].x(), single_patch[j][k].y(), single_patch[j][k].z()));
            }
            control_surface.push_back(control_line);
        }
        CBSplineSurface bsplineSurface(udegree, vdegree, control_surface, uknots, vknots, umultis, vmultis);
        muti_bezier_batch_dir_.emplace_back(bsplineSurface);
    }
    return true;
}

bool IGA_CC_Helper::cad_restruction(vector<vector<vector<Vec>>>& control_points)
{
    int udegree = 3;
    int vdegree = 3;
    vector<double> uknots{ 0,1 };
    vector<double> vknots{ 0,1 };
    vector<int> umultis{ 4,4 };
    vector<int> vmultis{ 4,4 };

    muti_bezier_batch_.clear();
    muti_bezier_batch_dir_.clear();
    for (int i = 0; i < control_points.size(); i++) {
        const auto& single_patch = control_points[i];
        vector<vector<CPoint>> control_surface;
        for (int j = 0; j < 4; ++j) {
            vector<CPoint> control_line;
            for (int k = 0; k < 4; ++k) {
                control_line.push_back(CPoint(single_patch[j][k].x(), single_patch[j][k].y(), single_patch[j][k].z()));
            }
            control_surface.push_back(control_line);
        }
        CBSplineSurface bsplineSurface(udegree, vdegree, control_surface, uknots, vknots, umultis, vmultis);
        muti_bezier_batch_.emplace_back(bsplineSurface);
    }
    return true;
}