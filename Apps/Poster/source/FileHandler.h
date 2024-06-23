#ifndef FILE_HANDLE_H
#define FILE_HANDLE_H
#include "Common/ComHeader.h"
#include "Common/GeometryKernel.h"
#include <iomanip> 
class FileHandler {
	
private:
	HexV3fMesh *pm_convert_from;
	std::ifstream mesh_iff;
	std::ofstream mesh_off;
	bool checkIsFilled;
	bool isReading;
	std::string filepath_checked;
private:
	bool getCleanLine(std::istream& _ifs, std::string& _string, bool _skipEmptyLines = true) const;
	void trimString(std::string& _string) const;
	void extractQuotedText(std::string& _string) const;
public:
	FileHandler(std::string str_file_path, std::string str_file_name);
	FileHandler(std::string path);
	FileHandler() {};
	bool read2HexMeshObj(HexV3fMesh *hexmesh, bool isWrite2File, bool _computBottom, bool _topologyCheck);
	bool writeFile(const std::string& _filename, HexV3fMesh* _mesh) const;
	bool writeFile1(const std::string& _filename, HexV3fMesh* _mesh) const;
	bool write2OffFile(const std::string& _filename, HexV3fMesh* _mesh)const;
	bool write2DatFile(const std::string& _filename, std::vector<HexV3fMesh>& _bc_mesh_vec)const;
	bool readVtk2HexMeshObj(HexV3fMesh * hexmesh);
	bool write2AXLFile(const std::string& _filename,std::vector<std::vector<std::vector<std::vector<V3f>>>> _SolidData,
		std::vector<std::vector<double>> _bezier_color);
	bool readpickedpoint(std::vector<std::vector<VertexHandle>> &pickpoint_vtoline);
	bool readpickedpoint(std::vector<std::vector<V3f>> &pickpoint_vtoline);
	bool readpickededge(std::vector<EdgeHandle> &pickededge);
	bool readpickedface(std::vector<FaceHandle> &pickedface);
	bool readpickedall(std::vector<VertexHandle> &pickedvertex, std::vector<EdgeHandle> &pickededge,
		std::vector<FaceHandle> &pickedface, std::vector<CellHandle> &pickedcell);
	bool writebdy(const std::string& _filename, const std::vector<std::vector<std::pair<int, int>>> &m_bdyquick)const;
	bool readbdy(std::vector<std::vector<std::pair<int, int>>> &m_bdyquick);
	bool write2VTKFile(const std::string& _filename, HexV3fMesh* _bc_mesh_vec)const;
	bool readcsv(std::vector<VertexHandle> &pickedface);
	bool readpickededge(std::vector<V3f> &edgemid);
	bool writeXJ(const std::string& _filename, HexV3fMesh * _mesh, std::vector<double> &m_color);
};

#endif // !FILE_HANDLE_H
