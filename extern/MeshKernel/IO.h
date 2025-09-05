#pragma once
#include "Mesh.h"
#include "Geom/CBSplineSurface.h"
#include <fstream>
#include <string>
#include <sstream>
#include <array>
#include <vector>
//Step读取
#include <..\opencascade\opencascade770-install\inc\TopoDS_Shape.hxx>
#include <..\opencascade\opencascade770-install\inc\TopoDS.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepPrimAPI_MakeBox.hxx>
#include <..\opencascade\opencascade770-install\inc\TopExp_Explorer.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepAdaptor_Surface.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepBndLib.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepMesh_IncrementalMesh.hxx>
#include <..\opencascade\opencascade770-install\inc\STEPControl_Reader.hxx>
#include <..\opencascade\opencascade770-install\inc\Poly_Triangulation.hxx>
//Step导出
#include <..\opencascade\opencascade770-install\inc\BRepBuilderAPI_MakePolygon.hxx>
#include <..\opencascade\opencascade770-install\inc\BRep_Builder.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepBuilderAPI_MakeVertex.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepBuilderAPI_MakeWire.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepBuilderAPI_MakeEdge.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepBuilderAPI_MakeFace.hxx>
#include <..\opencascade\opencascade770-install\inc\STEPControl_Writer.hxx>
#include <..\opencascade\opencascade770-install\inc\TopoDS_Vertex.hxx>
#include <..\opencascade\opencascade770-install\inc\Geom_BSplineSurface.hxx>
#include <..\opencascade\opencascade770-install\inc\TColgp_Array2OfPnt.hxx>
#include <..\opencascade\opencascade770-install\inc\GeomAPI_PointsToBSplineSurface.hxx>
//Iges读取
#include <..\opencascade\opencascade770-install\inc\IGESControl_Reader.hxx>
#include <..\opencascade\opencascade770-install\inc\Standard.hxx>
#include <..\opencascade\opencascade770-install\inc\Standard_Type.hxx>
#include <..\opencascade\opencascade770-install\inc\StepRepr_RepresentationItem.hxx>
#include <..\opencascade\opencascade770-install\inc\Interface_Static.hxx>
#include <..\opencascade\opencascade770-install\inc\IGESData_IGESModel.hxx>
#include <..\opencascade\opencascade770-install\inc\TCollection_AsciiString.hxx>
//Iges导出
#include <..\opencascade\opencascade770-install\inc\IGESControl_Writer.hxx>
//Stl导出
#include <..\opencascade\opencascade770-install\inc\TopoDS_Iterator.hxx>
#include <..\opencascade\opencascade770-install\inc\TopoDS_Face.hxx>
#include <..\opencascade\opencascade770-install\inc\StlAPI_Writer.hxx>
#include <..\opencascade\opencascade770-install\inc\BRepTools.hxx>
#include <..\opencascade\opencascade770-install\inc\TopoDS_Compound.hxx>
#include <..\opencascade\opencascade770-install\inc\gp_Pnt.hxx>

namespace MeshKernel {
	class IO {
	public:

		VolumeMesh iGameReadVolumeFile(const std::string&);// 读六面体网格
		VolumeMesh iGameReadVolumeFile_TopOptMeshGeneration(const std::string&);
		void iGameWriteVolumeFile(const VolumeMesh& _mesh, const std::string&);// 写六面体网格

		TetMesh iGameReadTetMeshFile(const std::string&);// 读四面体网格

        SurfaceMesh ReadObjFile(const std::string& _InputFile, int& sides_num);
        SurfaceMesh ReadOffFile(const std::string& _InputFile, int& sides_num);
        SurfaceMesh ReadStlFile(const std::string& _InputFile, int& sides_num);
        SurfaceMesh ReadIgesFile(const std::string& _InputFile, int& sides_num);
        SurfaceMesh ReadStepFile(const std::string& _InputFile, int& sides_num);
        //TetMesh ReadTetMeshFile(const std::string& _InputFile);
        bool WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
        bool WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
        bool WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
        bool WriteStepFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
        bool WriteObjFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& outputFile);
        bool WriteStlFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile);
        bool WriteStepFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile);
        bool WriteIgesFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile);
        //bool WriteMeshFile(const VolumeMesh& _mesh, const std::string& _OutputFile);
        bool WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile);
        std::string WriteOffString(const SurfaceMesh& _mesh);
        TetMesh ReadMeshFileFromStr(const std::string& data);
        std::vector<std::string> SplitFileName(const std::string& fileName);

	private:
		void ReOrderiGameVertexHandle(const SurfaceMesh& _mesh);
		std::vector<iGameVertexHandle> reorderedvh_;                        // 重排顶点
		std::unordered_map<iGameVertexHandle, std::size_t> newvh_;          // 新的顶点handle
	};
}
