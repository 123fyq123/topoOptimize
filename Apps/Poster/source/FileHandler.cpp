#include "FileHandler.h"

using namespace std;
FileHandler::FileHandler(string str_file_path, string str_file_name)
{
	if (str_file_path.empty()) {
		std::cout << "file path is empty!";
		isReading = false;
		checkIsFilled = false;
		return;
	}
	else if (str_file_name.empty()) {
		std::cout << "file name is empty!";
		isReading = false;
		checkIsFilled = false;
		return;
	}

	if (str_file_name.find("mesh") == string::npos&&
		str_file_name.find("MESH") == string::npos) {
		std::cout << "file type is wrong!";
		isReading = false;
		checkIsFilled = false;
		return;
	}

		
	std::stringstream stream;
	std::string result;


	stream << str_file_path;
	if (str_file_path.rfind("\\") != str_file_path.length() - 1 &&
		str_file_path.rfind("/") != str_file_path.length() - 1) {

		stream << "\\";

	}
	stream << str_file_name;
	stream >> result;

	filepath_checked = result;
	isReading = true;
	checkIsFilled = false;
}

FileHandler::FileHandler(std::string path)
{
	if (path.empty()) 
	{
		std::cout << "file path is empty!";
		isReading = false;
		checkIsFilled = false;
		return;
	}

	filepath_checked = path;
	isReading = true;
	checkIsFilled = false;
}

void FileHandler::extractQuotedText(string& _string) const {

	// Trim Both leading and trailing quote marks
	size_t start = _string.find_first_of("\""); ++start;
	size_t end = _string.find_last_not_of("\"");

	if ((string::npos == start) || (string::npos == end)) {
		_string = "";
	}
	else {
		_string = _string.substr(start, end - start + 1);
	}
}

void FileHandler::trimString(string& _string) const {

	// Trim Both leading and trailing spaces
	size_t start = _string.find_first_not_of(" \t\r\n");
	size_t end = _string.find_last_not_of(" \t\r\n");

	if ((string::npos == start) || (string::npos == end)) {
		_string = "";
	}
	else {
		_string = _string.substr(start, end - start + 1);
	}
}

bool FileHandler::getCleanLine(istream& _ifs, string& _string, bool _skipEmptyLines) const {

	// While we are not at the end of the file
	while (true) {

		// Get the current line:
		getline(_ifs, _string);

		// Remove whitespace at beginning and end
		trimString(_string);

		// Check if string is not empty ( otherwise we continue
		if (_string.size() != 0) {

			// Check if string is a comment ( starting with # )
			if (_string[0] != '#') {
				return true;
			}

		}
		else {
			if (!_skipEmptyLines)
				return true;
		}

		if (_ifs.eof() || _string.compare("End")) {
			cerr << "End of file reached while searching for input!" << endl;
			return false;
		}
	}

	return false;
}
//read mesh file to HexV3fMesh
bool FileHandler::read2HexMeshObj(HexV3fMesh *hexmesh, bool isWrite2File, bool _computBottom, bool _topologyCheck)
{
	
	if (!isReading) {
		return false;
	}
	this->pm_convert_from = hexmesh;

	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;
	string s_tmp;
	uint64_t c = 0u;

		
	int v_arr[8];

	pm_convert_from = hexmesh;

	//*************************
	//header
	bool header_found = true;

	getCleanLine(mesh_iff, line);
	sstr.str(line);

	// Check header
	sstr >> s_tmp;
	transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
	if (s_tmp != "MESHVERSIONFORMATTED") {
		//iff.close();
		header_found = false;
		mesh_iff.close();
		std::cerr << "The specified file might not be in .mesh format!" << std::endl;
		return false;
	}

	sstr.clear();
	getCleanLine(mesh_iff, line);
	sstr.str(line);

	// Get Dimension
	string s_dim;
	sstr >> s_tmp;
		
	sstr >> s_dim;
	transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
	if (s_tmp != "DIMENSION" || s_dim != "3") {
		mesh_iff.close();
		std::cerr << "Other dimension for files are not supported at the moment!" << std::endl;
		return false;
	}

	if (!header_found) {
		sstr.clear();
		sstr.str(line);
	}
	else {
		sstr.clear();
		getCleanLine(mesh_iff, line);
		sstr.str(line);
	}

	sstr >> s_tmp;
	std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
	if (s_tmp != "VERTICES") {
		mesh_iff.close();
		std::cerr << "No vertex section defined!" << std::endl;
		return false;
	}
	else {

        if(!sstr.eof())
        {
            sstr >> c;
            if(c == 0)
            {
                getCleanLine(mesh_iff, line);
                sstr.clear();
                sstr.str(line);
                sstr >> c;
            }
        }
        else
        {
            getCleanLine(mesh_iff, line);
            sstr.clear();
            sstr.str(line);
            sstr >> c;
        }

		for (uint64_t i = 0u; i < c; ++i) {
			/*if(i%10000 == 0)
			cout << i << ":" << c << endl;*/
			int check;
			getCleanLine(mesh_iff, line);
			sstr.clear();
			sstr.str(line);
			V3f v = V3f(0.00, 0.00, 0.00);
			sstr >> v.x;
			sstr >> v.y;
			sstr >> v.z;
			//v.z -= 1.5;
			sstr >> check;
			if (check == 0)
				v.is_valid = false;
			pm_convert_from->add_vertex(v);
		}
	}

    c = 0;

	while (true)
	{
        sstr.clear();
		if(!getCleanLine(mesh_iff, line))
        {
            throw std::runtime_error("No Hexahedra Information");
        }
        sstr.str(line);
        sstr >> s_tmp;
        std::transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
        if (s_tmp == "HEXAHEDRA") {

            if(!sstr.eof())
            {
                sstr >> c;
                if(c == 0)
                {
                    getCleanLine(mesh_iff, line);
                    sstr.clear();
                    sstr.str(line);
                    sstr >> c;
                }
            }
            else
            {
                getCleanLine(mesh_iff, line);
                sstr.clear();
                sstr.str(line);
                sstr >> c;
            }

            break;
        }
	}

	int temp_inde = 0;
	for (unsigned int f = 0; f < c; ++f) {
		/*if(f%10000 == 0)
			cout << f << ":" << c << endl;*/
		vector<VertexHandle> temp_vetex;
			
		getCleanLine(mesh_iff, line);
		sstr.clear();
		sstr.str(line);

		for (int i = 0; i < 8; ++i) {
			sstr >> temp_inde;
			temp_vetex.push_back(VertexHandle(temp_inde - 1));
		}

			
		CellHandle cell_now = pm_convert_from->add_cell(temp_vetex);

			
	}
	mesh_iff.close();
	hexmesh->examine_bdy();
	hexmesh->cal_cell_cen();
	hexmesh->cal_cen();
	return true;
}

bool FileHandler::writeFile(const std::string& _filename,HexV3fMesh* _mesh) const {

	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}

	// Write header
	off << "MeshVersionFormatted 1" << std::endl;
	off << "Dimension 3" << std::endl;
	uint64_t n_vertices(_mesh->vertices_.size());
	off << "Vertices" << std::endl;
	off << n_vertices << std::endl;

	// write vertices
	for (uint64_t v_it = 0; v_it<n_vertices ; ++v_it) {

		off << _mesh->vertices_[v_it].x << " " 
			<< _mesh->vertices_[v_it].y << " "
			<< _mesh->vertices_[v_it].z << " "
			<<"-1"<<std::endl;
	}
	uint64_t n_faces(_mesh->faces_.size());
	off << "Quadrilaterals" << std::endl;
	off << n_faces << std::endl;

	for (uint64_t c_it = 0; c_it<n_faces; ++c_it) {


		std::vector<VertexHandle> ves = _mesh->faces_[c_it].vs_f;

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx()) + 1) << " ";

		}
		off << "1" << std::endl;
	}
	
		uint64_t n_cells(_mesh->cells_.size());
		off << "Hexahedra" << std::endl;
		//off << n_cells << std::endl;
		int cnt = 0;
		for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
			if (!_mesh->cells_[c_it].is_valid) continue;
			cnt++;
		}
		off << cnt << std::endl;

		for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
			if (!_mesh->cells_[c_it].is_valid) continue;


			std::vector<VertexHandle> ves = _mesh->cells_[c_it].vertex_;

			for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
				!= ves.end(); ++it_ves) {
				off << ((it_ves->idx()) + 1) << " ";

			}
			off << "1" << std::endl;
		}
	
	
	off << "End" << std::endl;
	off.close();

	return true;
}

bool FileHandler::write2OffFile(const std::string& _filename, HexV3fMesh* _mesh) const {

	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}

	// Write header
	off << "OFF" << std::endl;
	off << _mesh->vertices_.size()<<" "<<_mesh->faces_.size()<<" 0"<< std::endl;
	uint64_t n_vertices(_mesh->vertices_.size());
	
	// write vertices
	for (uint64_t v_it = 0; v_it<n_vertices; ++v_it) {

		off << std::fixed << std::setprecision(6)<<
			_mesh->vertices_[v_it].x << " "
			<< _mesh->vertices_[v_it].y << " "
			<< _mesh->vertices_[v_it].z << " "
			<< std::endl;
	}
	//清除点数模式
	off.unsetf(ios::fixed);
	uint64_t n_faces(_mesh->faces_.size());
	
	for (uint64_t c_it = 0; c_it<n_faces; ++c_it) {
		
		off << "4 ";
		std::vector<VertexHandle> ves = _mesh->faces_[c_it].vs_f;
		
		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << (it_ves->idx()) << " ";

		}
		off <<std::endl;
	}

	off << "End" << std::endl;
	off.close();

	return true;
}

bool FileHandler::write2DatFile(const std::string & _filename, std::vector<HexV3fMesh>& _bc_mesh_vec) const
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}

	// Write header
	off << "TITLE = \"title\" VARIABLES = \"X\", \"Y\", \"Z\", \"elementcolor\"" << "\n";
	uint16_t bc_index = 0;
	//循环写n个bc网格
	uint32_t color[4] = { 0xD2B48C,
	0xB22222,
	0x87CEEB,
	0x228B22
	};
	for (bc_index = 0; bc_index < _bc_mesh_vec.size(); ++bc_index) {
		off << "ZONE T =" << "zone_" << (int)bc_index
			<< ", N = " << (int)_bc_mesh_vec[bc_index].vertices_.size()
			<< ", E = " << (int)_bc_mesh_vec[bc_index].cells_.size() << ", ET = BRICK, F = FEPOINT\n";

		uint64_t n_vertices(_bc_mesh_vec[bc_index].vertices_.size());
	
		// write vertices
		for (uint64_t v_it = 0; v_it<n_vertices; ++v_it) {
			/*if (bc_index == 686 || bc_index == 435 || bc_index == 750|| bc_index == 609 || bc_index == 627)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[0]) << "\n";
			}
			else if (bc_index == 674 || bc_index == 434 || bc_index == 740 || bc_index == 610 || bc_index == 620 || bc_index == 652)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[1]) << "\n";
			}
			else if (bc_index == 672 || bc_index == 756 || bc_index == 611 || bc_index == 485 || bc_index == 655)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[2]) << "\n";
			}
			else if (bc_index == 638|| bc_index == 664 || bc_index == 670 || bc_index == 432 || bc_index == 768)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[3]) << "\n";
			}
			else
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[bc_index % 4]) << "\n";
			}*/
			if (bc_index == 1)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[0]) << "\n";
			}
			else if (bc_index == 2)
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[3]) << "\n";
			}
			else
			{
				off << std::fixed << std::setprecision(6)
					<< _bc_mesh_vec[bc_index].vertices_[v_it].x << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
					<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
					<< (double)(color[bc_index % 4]) << "\n";
			}
			/*off << std::fixed << std::setprecision(6)
				<<_bc_mesh_vec[bc_index].vertices_[v_it].x << " "
				<< _bc_mesh_vec[bc_index].vertices_[v_it].y << " "
				<< _bc_mesh_vec[bc_index].vertices_[v_it].z << " "
				<< (double)(color[bc_index%4]) << "\n";*/
			
		}

		//清除点数模式
		off.unsetf(ios::fixed);
		//write hex

		uint64_t n_cells(_bc_mesh_vec[bc_index].cells_.size());

		for (uint64_t v_it = 0; v_it<n_cells; ++v_it) {
			off << (int)(_bc_mesh_vec[bc_index].cells_[v_it].vertex_[0]+1) ;
			for (uint8_t v_index = 1; v_index < 8; ++v_index) {
				off << " "<<(int)(_bc_mesh_vec[bc_index].cells_[v_it].vertex_[v_index]+1);
			}
			off << "\n";
		}

	}

	off.close();

	return true;

}

bool FileHandler::readVtk2HexMeshObj(HexV3fMesh * hexmesh)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;
	string s_tmp;
	uint64_t c = 0u;


	int v_arr[8];

	pm_convert_from = hexmesh;

	//*************************
	//header
	bool header_found = true;

	getCleanLine(mesh_iff, line);
	getCleanLine(mesh_iff, line);
	getCleanLine(mesh_iff, line);
	sstr.str(line);

	// Check header
	sstr >> s_tmp;
	transform(s_tmp.begin(), s_tmp.end(), s_tmp.begin(), ::toupper);
	if (s_tmp != "ASCII") {
		//iff.close();
		header_found = false;
		mesh_iff.close();
		std::cerr << "The specified file might not be in .vtk format!" << std::endl;
		return false;
	}
	sstr.clear();

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	sstr.clear();

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	std::string type;
	uint64_t count_point;
	sstr >> type;
	sstr >> count_point;
	sstr >> type;

	for (uint64_t i = 0u; i < count_point; ++i) {

		getCleanLine(mesh_iff, line);
		sstr.str("");
		sstr.clear();
		sstr.str(line);
		V3f v = V3f(0.00, 0.00, 0.00);
		sstr >> v.x;
		sstr >> v.y;
		sstr >> v.z;
		pm_convert_from->add_vertex(v);
	}


	sstr.clear();
	getCleanLine(mesh_iff, line);
	sstr.str(line);
	//std::string type;
	uint64_t cell_point;
	sstr >> type;
	sstr >> cell_point;
	sstr >> type;

	int temp_inde;


	for (unsigned int f = 0; f < cell_point; ++f) {
		//cout << f << ":" << 62858 << endl;
		vector<VertexHandle> temp_vetex;
		getCleanLine(mesh_iff, line);
		sstr.clear();
		sstr.str(line);
		sstr >> temp_inde;
		for (int i = 0; i < 8; ++i) {
			sstr >> temp_inde;
			temp_vetex.push_back(VertexHandle(temp_inde));
		}


		CellHandle cell_now = pm_convert_from->add_cell(temp_vetex);


	}
	mesh_iff.close();

	return true;
}

bool FileHandler::write2AXLFile(const std::string& _filename, std::vector<std::vector<std::vector<std::vector<V3f>>>> _SolidData,
	std::vector<std::vector<double>> _bezier_color)
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}
	off << "<axl>" << "\n";
	int i = 0;
	int v_sum = 4 * 4 * 4;
	for (std::vector<std::vector<std::vector<std::vector<V3f>>>>::iterator _SolidData_it = _SolidData.begin();_SolidData_it != _SolidData.end();
		_SolidData_it++)
	{
		off << "<volume type=" << '"' << "bspline" << '"' << " name=" << '"' << "bv_0" << '"' << " color=" << '"';
		off << _bezier_color[i][0] << " " << _bezier_color[i][1] << " " << _bezier_color[i][2] << '"' << '>' << "\n";
		i++;
		off << "  <number>4 4 4</number>" << "\n";
		off << "  <order>4 4 4</order>" << "\n";
		off << "  <knots>0 0 0 0 1 1 1 1 </knots>" << "\n";
		off << "  <knots>0 0 0 0 1 1 1 1 </knots>" << "\n";
		off << "  <knots>0 0 0 0 1 1 1 1 </knots>" << "\n";
		off << "  <points>";
		int j = 0;
		std::vector<std::vector<std::vector<V3f>>> solid = *_SolidData_it;
		for (std::vector<std::vector<std::vector<V3f>>>::iterator solid_it = solid.begin();solid_it != solid.end();solid_it++)
		{
			std::vector<std::vector<V3f>> face = *solid_it;
			for (std::vector<std::vector<V3f>>::iterator face_it = face.begin();face_it != face.end();face_it++)
			{
				std::vector<V3f> line = *face_it;
				for (std::vector<V3f>::iterator line_it = line.begin();line_it != line.end();line_it++)
				{
					V3f vertex = *line_it;
					j++;
					if (j == v_sum)
					{
						off << std::fixed << std::setprecision(6)
							<< vertex.x << " "
							<< vertex.y << " "
							<< vertex.z;
					}
					else
					{
						off << std::fixed << std::setprecision(6)
							<< vertex.x << " "
							<< vertex.y << " "
							<< vertex.z << " "
							<< "\n";
					}	
				}
			}
		}
		off << "</points>" << "\n";
		off << "</volume>" << "\n";
	}
	off << "<script></script>" << "\n";
	off << "<camera>" << "\n";
	off << " <Parameters fieldOfView=" << '"' << "0.785398" << '"' << " Type=" << '"' << "PERSPECTIVE" << '"' << "\n";
	off << "zNearCoefficient=" << '"' << "0.005" << '"' << " orthoCoef=" << '"' << "0.414214" << '"' << " zClippingCoefficient=" << '"' << "1.73205" << '"' << "\n";
	off << "/>" << "\n";
	off << " <Stereo distToScreen=" << '"' << "0.5" << '"' << " focusDistance=" << '"' << "241.421" << '"' << " physScreenWidth=" << '"' << "0.4" << '"' << "\n";
	off << "IODist=" << '"' << "0.062" << '"' << " />" << "\n";
	off << " <ManipulatedCameraFrame>" << "\n";
	off << "  <position x=" << '"' << "-8.36938" << '"' << " y=" << '"' << "0.880508" << '"' << " z=" << '"' << "7.02155" << '"' << " />" << "\n";
	off << "  <orientation q0=" << '"' << "-0.498947" << '"' << " q1=" << '"' << "0.060419" << '"' << " q2=" << '"' << "0.842227" << '"' << " q3=" << '"' << "-0.195076" << '"' << " />" << "\n";
	off << "  <ManipulatedParameters transSens=" << '"' << "1" << '"' << " rotSens=" << '"' << "1" << '"' << " wheelSens=" << '"' << "1" << '"' << "\n";
	off << "spinSens=" << '"' << "0.3" << '"' << " />" << "\n";
	off << "  <ManipulatedCameraParameters flySpeed=" << '"' << "1" << '"' << " >" << "\n";
	off << "   <flyUpVector x=" << '"' << "0" << '"' << " y=" << '"' << "1" << '"' << " z=" << '"' << "0" << '"' << " />" << "\n";
	off << "  </ ManipulatedCameraParameters>" << "\n";
	off << " </ManipulatedCameraFrame>" << "\n";
	off << "</camera>"<<"\n";
	off << "</axl>" << "\n";
}

bool FileHandler::writeFile1(const std::string& _filename, HexV3fMesh* _mesh) const
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << "\r\n";
		off.close();
		return false;
	}

	// Write header
	off << "MeshVersionFormatted 1" << "\r\n";
	off << "Dimension 3" << "\r\n";
	uint64_t n_vertices(_mesh->vertices_.size());
	off << "Vertices" << "\r\n";
	off << n_vertices << "\r\n";

	// write vertices
	for (uint64_t v_it = 0; v_it<n_vertices; ++v_it) {

		off << _mesh->vertices_[v_it].x << " "
			<< _mesh->vertices_[v_it].y << " "
			<< _mesh->vertices_[v_it].z << " ";
		if (_mesh->vertices_[v_it].is_valid == true)
			off << 1 << " ";
		else
			off << 0 << " ";
		off << "-1" << "\r\n";
	}
	uint64_t n_faces(_mesh->faces_.size());
	off << "Quadrilaterals" << "\r\n";
	off << n_faces << "\r\n";

	for (uint64_t c_it = 0; c_it<n_faces; ++c_it) {


		std::vector<VertexHandle> ves = _mesh->faces_[c_it].vs_f;

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx()) + 1) << " ";

		}
		off << "1" << "\r\n";
	}

	uint64_t n_cells(_mesh->cells_.size());
	off << "Hexahedra" << "\r\n";
	uint32_t cout_ = 0;
	for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
		if (_mesh->cells_[c_it].is_valid) ++cout_;
	}
	off << cout_ << "\r\n";

	for (uint64_t c_it = 0; c_it<n_cells; ++c_it) {
		if (!_mesh->cells_[c_it].is_valid) continue;

		std::vector<VertexHandle> ves = _mesh->cells_[c_it].vertex_;

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx()) + 1) << " ";

		}
		off << "1" << "\r\n";
	}


	off << "End" << "\r\n";
	off.close();

	return true;
}

bool FileHandler::readpickedpoint(std::vector<std::vector<VertexHandle>> &pickpoint_vtoline)
{
	std::cout<<"Start read picked point"<<std::endl;
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int pickpoint_vtolin_size;
	sstr >> pickpoint_vtolin_size;
	for (int i = 0; i < pickpoint_vtolin_size; i++)
	{
		std::cout << "Schedule of read picked point " << i + 1 << ":" << pickpoint_vtolin_size << std::endl;
		getCleanLine(mesh_iff, line);
		sstr.clear();
		sstr.str(line);
		int pickpoint_size;
		sstr >> pickpoint_size;
		std::vector<VertexHandle> temp_v;

		getCleanLine(mesh_iff, line);
		sstr.clear();
		sstr.str(line);
		for (int j = 0; j < pickpoint_size; j++)
		{
			int v_id;
			sstr >> v_id;
			VertexHandle v(v_id);
			temp_v.push_back(v);
		}
		pickpoint_vtoline.push_back(temp_v);
	}

	mesh_iff.close();

	std::cout<<"Read picked point completely"<<std::endl;
	return true;
}

bool FileHandler::readpickedpoint(std::vector<std::vector<V3f>> &pickpoint_vtoline)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int pickpoint_vtolin_size;
	sstr >> pickpoint_vtolin_size;

	for (int i = 0; i < pickpoint_vtolin_size; i++)
	{
		getCleanLine(mesh_iff, line);
		sstr.clear();
		sstr.str(line);
		int pickpoint_size;
		sstr >> pickpoint_size;
		std::vector<V3f> temp_v;

		for (int j = 0; j < pickpoint_size; j++)
		{
			getCleanLine(mesh_iff, line);
			sstr.str("");
			sstr.clear();
			sstr.str(line);
			V3f v = V3f(0.00, 0.00, 0.00);
			sstr >> v.x;
			sstr >> v.y;
			sstr >> v.z;
			temp_v.push_back(v);
		}
		pickpoint_vtoline.push_back(temp_v);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::readpickededge(std::vector<EdgeHandle>& pickededge)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int pickededge_size;
	sstr >> pickededge_size;

	getCleanLine(mesh_iff, line);
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickededge_size; j++)
	{
		int e_id;
		sstr >> e_id;
		EdgeHandle e(e_id);
		pickededge.push_back(e);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::readpickedface(std::vector<FaceHandle>& pickedface)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int pickedface_size;
	sstr >> pickedface_size;

	getCleanLine(mesh_iff, line);
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickedface_size; j++)
	{
		int f_id;
		sstr >> f_id;
		FaceHandle f(f_id);
		pickedface.push_back(f);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::readpickedall(std::vector<VertexHandle> &pickedvertex, std::vector<EdgeHandle> &pickededge,
	std::vector<FaceHandle> &pickedface, std::vector<CellHandle> &pickedcell)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	int pickedvertex_size;
	sstr >> pickedvertex_size;
	
	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickedvertex_size; j++)
	{
		int v_id;
		sstr >> v_id;
		VertexHandle v(v_id);
		pickedvertex.push_back(v);
	}

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	int pickededge_size;
	sstr >> pickededge_size;

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickededge_size; j++)
	{
		int e_id;
		sstr >> e_id;
		EdgeHandle e(e_id);
		pickededge.push_back(e);
	}

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	int pickedface_size;
	sstr >> pickedface_size;
	
	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickedface_size; j++)
	{
		int f_id;
		sstr >> f_id;
		FaceHandle f(f_id);
		pickedface.push_back(f);
	}

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	int pickedcell_size;
	sstr >> pickedcell_size;

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < pickedcell_size; j++)
	{
		int c_id;
		sstr >> c_id;
		CellHandle c(c_id);
		pickedcell.push_back(c);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::writebdy(const std::string& _filename, const std::vector<std::vector<std::pair<int, int>>> &m_bdyquick)const
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}


	// write picked point
	off << m_bdyquick.size() << "\r\n";

	for (int i = 0; i < m_bdyquick.size(); i++)
	{
		std::vector<std::pair<int, int>> temp = m_bdyquick[i];
		off << temp.size() << "\r\n";
		for (int j = 0; j < temp.size(); ++j)
		{
			off << temp[j].first << " " << temp[j].second << "\r\n";
		}
		
	}
	off << "\r\n";


	off.close();

	return true;


}

bool FileHandler::readbdy(std::vector<std::vector<std::pair<int, int>>> &m_bdyquick)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str("");
	sstr.clear();
	sstr.str(line);
	int bdyquick_size;
	sstr >> bdyquick_size;

	
	for (int j = 0; j < bdyquick_size; j++)
	{
		getCleanLine(mesh_iff, line);
		sstr.str("");
		sstr.clear();
		sstr.str(line);
		int tempsize;
		sstr >> tempsize;
		std::vector<std::pair<int, int>> temp;
		for (int i = 0; i < tempsize; ++i)
		{
			getCleanLine(mesh_iff, line);
			sstr.str("");
			sstr.clear();
			sstr.str(line);
			int a, b;
			sstr >> a >> b;
			temp.push_back(pair<int, int>(a, b));
		}
		m_bdyquick.push_back(temp);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::write2VTKFile(const std::string & _filename, HexV3fMesh* _mesh) const
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}

	// Write header
	off << "# vtk DataFile Version 2.0" << std::endl;
	off << "mesh vtk data - converted from .off" << std::endl;
	off << "ASCII" << std::endl;
	off << "DATASET UNSTRUCTURED_GRID" << std::endl;
	uint64_t n_vertices(_mesh->vertices_.size());
	off << "POINTS ";
	off << n_vertices <<" double" <<std::endl;

	// write vertices
	for (uint64_t v_it = 0; v_it < n_vertices; ++v_it) {

		off << _mesh->vertices_[v_it].x << " "
			<< _mesh->vertices_[v_it].y << " "
			<< _mesh->vertices_[v_it].z << std::endl;
	}
	
	
	uint64_t n_cells(_mesh->cells_.size());
	
	uint32_t cout_ = 0;
	for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
		if (_mesh->cells_[c_it].is_valid) ++cout_;
	}
	off << "CELLS ";
	off << cout_ << " " << cout_ * 9 << std::endl;

	for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
		if (!_mesh->cells_[c_it].is_valid) continue;

		std::vector<VertexHandle> ves = _mesh->cells_[c_it].vertex_;

		off << 8 << " ";

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx())) << " ";

		}
		off << std::endl;
	}

	off << "CELL_TYPES  ";
	off << cout_ << std::endl;

	for(int i = 0;i< cout_;++i)
		off << 12 << std::endl;

	off.close();

	return true;
}

bool FileHandler::readcsv(std::vector<VertexHandle>& csv)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int csv_size;
	sstr >> csv_size;

	getCleanLine(mesh_iff, line);
	sstr.clear();
	sstr.str(line);
	for (int j = 0; j < csv_size; j++)
	{
		int csv_id;
		sstr >> csv_id;
		VertexHandle c(csv_id);
		csv.push_back(c);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::readpickededge(std::vector<V3f>& edgemid)
{
	mesh_iff = ifstream(filepath_checked.c_str(), ios::in);

	if (!mesh_iff.good()) {
		std::cerr << "Error: Could not open file " << filepath_checked << " for reading!" << std::endl;
		mesh_iff.close();
		return false;
	}

	stringstream sstr;
	string line;

	getCleanLine(mesh_iff, line);
	sstr.str(line);
	int pickededge_size;
	sstr >> pickededge_size;

	for (int j = 0; j < pickededge_size; j++)
	{
		getCleanLine(mesh_iff, line);
		sstr.str("");
		sstr.clear();
		sstr.str(line);
		V3f v = V3f(0.00, 0.00, 0.00);
		sstr >> v.x;
		sstr >> v.y;
		sstr >> v.z;
		edgemid.push_back(v);
	}

	mesh_iff.close();

	return true;
}

bool FileHandler::writeXJ(const std::string& _filename, HexV3fMesh * _mesh, std::vector<double> &m_color)
{
	std::ofstream off(_filename.c_str(), std::ios::out);

	if (!off.good()) {
		std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
		off.close();
		return false;
	}

	// Write header
	off << "MeshVersionFormatted 1" << std::endl;
	off << "Dimension 3" << std::endl;
	uint64_t n_vertices(_mesh->vertices_.size());
	off << "Vertices" << std::endl;
	off << n_vertices << std::endl;

	// write vertices
	for (uint64_t v_it = 0; v_it < n_vertices; ++v_it) {

		off << _mesh->vertices_[v_it].x << " "
			<< _mesh->vertices_[v_it].y << " "
			<< _mesh->vertices_[v_it].z << " "
			<< m_color[v_it] << " ";
		if (_mesh->vertices_[v_it].is_valid == true)
			off << 1 << " ";
		else
			off << 0 << " ";
		off << "-1" << endl;
	}
	uint64_t n_faces(_mesh->faces_.size());
	off << "Quadrilaterals" << std::endl;
	off << n_faces << std::endl;

	for (uint64_t c_it = 0; c_it < n_faces; ++c_it) {


		std::vector<VertexHandle> ves = _mesh->faces_[c_it].vs_f;

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx()) + 1) << " ";

		}
		off << "1" << std::endl;
	}

	uint64_t n_cells(_mesh->cells_.size());
	off << "Hexahedra" << std::endl;
	uint32_t cout_ = 0;
	for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
		if (_mesh->cells_[c_it].is_valid) ++cout_;
	}
	off << cout_ << std::endl;

	for (uint64_t c_it = 0; c_it < n_cells; ++c_it) {
		if (!_mesh->cells_[c_it].is_valid) continue;

		std::vector<VertexHandle> ves = _mesh->cells_[c_it].vertex_;

		for (std::vector<VertexHandle>::iterator it_ves = ves.begin(); it_ves
			!= ves.end(); ++it_ves) {
			off << ((it_ves->idx()) + 1) << " ";

		}
		off << "1" << std::endl;
	}


	off << "End" << std::endl;
	off.close();

	return true;
}