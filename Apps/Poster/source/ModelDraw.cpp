#include "ModelDraw.h"
#include "Data.h"
#include "singalFather.h"
#include "TimeIntegrator.h"
#include "ColorBar.h"
#include "CoverageAxis.h"

using namespace std;
using namespace Poster_Needed;
ModelDraw::ModelDraw(QWidget* parent) :ViewOpenGL(parent), isEnableLighting(false), isTwoSideLighting(false), editmodel(DEFULTEDIT) {
	/* 初始化 colorbar */
	genColors();
	initColorbar();
	/*Temp-lee*/
	oldHexMesh = new HexV3fMesh();
	depth = 0;
	bshowAxis = false;
	bshowGround = false;
	bshowSKEL = false;
	bdragMode = false;
	copyTime = false;
	bVertexSnap = false;
	bshowTriMesh = false;
	solutionSkel = new SolutionSkel();
	projectSkel = new Project();
	mouse_move_flag = false;
	local_coordinater_axis = new CoordinateAxis();
}
ModelDraw::~ModelDraw(void) {

}

bool ModelDraw::loadTriMesh(std::string filename) {
	ModelClear();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	//std::cout << "file format: " << fileformat << std::endl;
	tri_mesh = new iGame_OpenMesh_TriMesh();
	tri_mesh->request_face_normals();
	tri_mesh->request_vertex_normals();
	if (!OpenMesh::IO::read_mesh(*tri_mesh, filename)) {
		std::cout << "Can not read obj" << std::endl;
		return false;
	}
	if (!tri_mesh->is_trimesh()) {
		std::cout << "Support trimesh only!" << std::endl;
		return false;
	}


	SetDrawModule(TRIMESH);
	if (resetViewFlag) updateMeshModel();
	updateDrawData_TriMesh();
	update();
	return true;
}

bool ModelDraw::saveTriMesh(std::string filename) {
	if (!tri_mesh) return false;

	OpenMesh::IO::write_mesh(*tri_mesh, filename);

	return true;
}

void ModelDraw::smoothing_openmesh(size_t iter) {
	iGame_OpenMesh_Smoother app(*tri_mesh);
	app.smoothing(iter);
	updateDrawData();
	update();
}

void ModelDraw::remeshing_uniform_openmesh(double ratio) {
	iGame_OpenMesh_Remeshing app(*tri_mesh);
	app.remeshing_uniform(ratio, 10);
	updateDrawData();
	update();
}

void ModelDraw::simplification_qem_openmesh(size_t target_vcnt) {
	iGame_OpenMesh_Simplification app(*tri_mesh);
	app.simplification_qem(target_vcnt);
	updateDrawData();
	update();
}

void ModelDraw::denoising_bnf_openmesh(size_t iter = 1) {
	for (size_t i = 0; i < iter; ++i) {
		iGame_OpenMesh_Denoising app(*tri_mesh);
		app.denoising_bnf();
	}
	updateDrawData();
	update();
}

bool ModelDraw::loadSurfaceMesh(std::string filename) {
	ModelClear();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	//std::cout << "file format: " << fileformat << std::endl;

	MeshKernel::IO io;
	int tmp;
	if (fileformat == "obj") {
		this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename, tmp));
	}
	else if (fileformat == "off") {
		this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadOffFile(filename, tmp));
	}
    else if (fileformat == "stl") {
        this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadStlFile(filename, tmp));
    }
    else if (fileformat == "IGS") {
        this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadIgesFile(filename, tmp));
    }
    else if (fileformat == "step") {
        this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadStepFile(filename, tmp));
    }
    else{
        std::cout << "support obj and off files only!\n";
        return false;
    }

	///* iGameKernel 三角网格读入以及 KD 树使用示例 */
	//MeshKernel::SurfaceMesh* mesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename, tmp));
	///*if (fileformat == "obj") {
	//	this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadObjFile(filename, tmp));
	//} 	else if (fileformat == "off") {
	//	this->surface_mesh = new MeshKernel::SurfaceMesh(io.ReadOffFile(filename, tmp));
	//}*/
	//TriangleKdTree* kd_tree = new TriangleKdTree(*surface_mesh, 1);
	//for (auto& vp : mesh->allvertices()) {
	//	auto pos = vp.second;
	//	std::cout << "Input pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << "; ";
	//	auto face_nearest = kd_tree->nearest(pos);// 输入一个顶点(位置), 函数返回该曲面模型上离这个点最近的点坐标(可能在面上或边上)
	//	pos = face_nearest.nearest;
	//	std::cout << "Nearest pos: " << pos.x() << ", " << pos.y() << ", " << pos.z() << std::endl;
	//}
	bvh_tree->destoryBVH_Tree();
	bvh_tree->buildBVH_Tree(*surface_mesh);
	SetDrawModule(SURFACEMESH);
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	update();
	return true;
}

bool ModelDraw::saveSurfaceMesh(std::string filename) {

	//if (!this->surface_mesh) return false;
    if (this->surface_mesh) surface_mesh->updateAllHandles();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	MeshKernel::IO io;
	if (fileformat == "obj") {
        switch (drawmodule) {
            case DrawModule::EMPTY:
                break;

            case DrawModule::SURFACEMESH:
                io.WriteObjFile(*(this->surface_mesh), filename);
                break;

            case DrawModule::TETMESH:
                break;

            case DrawModule::HEXMESH:
                break;

            case DrawModule::SKEL:
                break;

            case DrawModule::BSPLINESURFACE:
                io.WriteObjFileWithSpline(this->surfacesDraw, filename);
                break;

            default:
                break;
        }
	}
	else if (fileformat == "off") {
		io.WriteOffFile(*(this->surface_mesh), filename);
	}
    else if (fileformat == "stl") {
        switch (drawmodule) {
            case DrawModule::EMPTY:
                break;

            case DrawModule::SURFACEMESH:
                io.WriteStlFile(*(this->surface_mesh), filename);
                break;

            case DrawModule::TETMESH:
                break;

            case DrawModule::HEXMESH:
                break;

            case DrawModule::SKEL:
                break;

            case DrawModule::BSPLINESURFACE:
                io.WriteStlFileWithSpline(this->surfacesDraw, filename);
                break;

            default:
                break;
        }
    }
    else if (fileformat == "IGS") {
        //io.WriteStepFile(*(this->surface_mesh), filename);
        switch (drawmodule) {

            case DrawModule::EMPTY:
                break;

            case DrawModule::SURFACEMESH:
                break;

            case DrawModule::TETMESH:
                break;

            case DrawModule::HEXMESH:
                break;

            case DrawModule::SKEL:
                break;

            case DrawModule::BSPLINESURFACE:
                io.WriteIgesFileWithSpline(this->surfacesDraw, filename);
                break;

            default:
                break;
        }
    }
    else if (fileformat == "step") {
        //io.WriteStepFile(*(this->surface_mesh), filename);
        switch (drawmodule) {

            case DrawModule::EMPTY:
                break;

            case DrawModule::SURFACEMESH:
                io.WriteStepFile(*(this->surface_mesh), filename);
                break;

            case DrawModule::TETMESH:
                break;

            case DrawModule::HEXMESH:
                break;

            case DrawModule::SKEL:
                break;

            case DrawModule::BSPLINESURFACE:
                io.WriteStepFileWithSpline(this->surfacesDraw, filename);
                break;

            default:
                break;
        }
    }
	else if (fileformat == "vtk") {
        io.WriteVtkFile(*(this->surface_mesh), filename);
	}
	return true;
}

bool ModelDraw::mesh_to_iges(std::string meshFileName, std::string igesFileName) {
	// TODO inplement mesh to iges
	if (meshFileName == "" || igesFileName == "") return false;

	HexV3fMesh* mesh = new HexV3fMesh();
	FileHandler file(meshFileName);
	file.read2HexMeshObj(mesh, true, true, true);// read mesh file
	mesh->examine_bdy();
	mesh->cal_cell_cen();
	mesh->cal_cen();
	mesh->examine_cs();
	mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();

	int udegree = 3;
	int vdegree = 3;
	vector<double> uknots{ 0,1 };
	vector<double> vknots{ 0,1 };
	vector<int> umultis{ 4,4 };
	vector<int> vmultis{ 4,4 };

	HexV3fMesh* SkelMesh = mesh;
	std::vector<EdgeHandle> TempEdge;
	ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

	int idx = 0;
	AppCC.approximationing();

	vector<CBSplineSurface> MutiBezierPatch;
	for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
		if (SkelMesh->is_bdy((FaceHandle)i)) {
			vector<vector<CPoint>> controlPointsSurface;
			for (int j = 0; j < 4; ++j) {
				vector<CPoint> controlPointLine;
				for (int k = 0; k < 4; ++k) {
					controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
				}
				controlPointsSurface.push_back(controlPointLine);
			}
			CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
			MutiBezierPatch.push_back(splineSurface);
			idx++;
		}
	}
	CAidFunction::getIGESSurface(igesFileName, MutiBezierPatch);// write iges file

	return true;
}

bool ModelDraw::saveIGESFromMesh(std::string igesFileName) {
	//if (!this->hexMesh || igesFileName == "") return false;
	//
	//int udegree = 3;
	//int vdegree = 3;
	//vector<double> uknots{ 0,1 };
	//vector<double> vknots{ 0,1 };
	//vector<int> umultis{ 4,4 };
	//vector<int> vmultis{ 4,4 };

	//HexV3fMesh* SkelMesh = this->hexMesh;
	//std::vector<EdgeHandle> TempEdge;
	//ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

	//int idx = 0;
	//AppCC.approximationingSharp();

	//vector<CBSplineSurface> MutiBezierPatch;
	//for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
	//	if (SkelMesh->is_bdy((FaceHandle)i)) {
	//		vector<vector<CPoint>> controlPointsSurface;
	//		for (int j = 0; j < 4; ++j) {
	//			vector<CPoint> controlPointLine;
	//			for (int k = 0; k < 4; ++k) {
	//				controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
	//			}
	//			controlPointsSurface.push_back(controlPointLine);
	//		}
	//		CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
	//		MutiBezierPatch.push_back(splineSurface);
	//		idx++;
	//	}
	//}
	//CAidFunction::getIGESSurface(igesFileName, MutiBezierPatch);// write iges file

	return true;
}

bool ModelDraw::hexMeshConvertToIGES() {
	if (!this->hexMesh) return false;

	std::string filename = "./Resources/tmp_convert.mesh";

	/*MeshKernel::IO io;
	io.iGameWriteVolumeFile(*hexMesh, filename);*/

	HexV3fMesh* _hexMesh = new HexV3fMesh();

	/* iGame_hexMesh to XIE_hexMesh */
	{
		uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();
		for (uint64_t i = 0; i < vcnt; ++i) {
			MeshKernel::iGameVertexHandle vh(i);
			auto& v = hexMesh->vertices(vh);
			V3f _v(v.x(), v.y(), v.z());
			_hexMesh->add_vertex(_v);
		}

		for (uint64_t i = 0; i < ccnt; ++i) {
			MeshKernel::iGameCellHandle ch(i);
			if (hexMesh->isValid(ch)) {
				auto& c = hexMesh->cells(ch);
				auto vhs = c.getVertexHandle();
				vector<VertexHandle> _vhs;
				for (auto& vh : vhs) {
					_vhs.push_back(VertexHandle(vh));
				}
				_hexMesh->add_cell(_vhs);
			}
		}

		/*std::cout << "iGameHexMesh convert to OldHexMesh success. OldHexMesh: #V = " << _hexMesh->vertices_.size() <<
			", #C = " << _hexMesh->cells_.size() << std::endl;*/

	}

	std::cout << "trans success" << std::endl;
	/*FileHandler file(filename);
	assert(file.read2HexMeshObj(_hexMesh, true, true, true));*/

	if (_hexMesh->vertices_.size() == 0) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
		return false;
	}

	_hexMesh->examine_bdy();
	_hexMesh->cal_cell_cen();
	_hexMesh->cal_cen();
	_hexMesh->examine_cs();
	//_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	//_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
	std::cout << "Converting... hexMesh: vertices = " << _hexMesh->vertices_.size() << ", cells = " << _hexMesh->cells_.size() << std::endl;
	int udegree = 3;
	int vdegree = 3;
	vector<double> uknots{ 0,1 };
	vector<double> vknots{ 0,1 };
	vector<int> umultis{ 4,4 };
	vector<int> vmultis{ 4,4 };

	HexV3fMesh* SkelMesh = _hexMesh;
	std::vector<EdgeHandle> TempEdge;
	ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

	int idx = 0;
	AppCC.approximationing();

	vector<CBSplineSurface> MutiBezierPatch;
	vector<CBSplineSurface> dx;

	for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
		if (SkelMesh->is_bdy((FaceHandle)i)) {
			vector<vector<CPoint>> controlPointsSurface;
			for (int j = 0; j < 4; ++j) {
				vector<CPoint> controlPointLine;
				for (int k = 0; k < 4; ++k) {
					controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
				}
				controlPointsSurface.push_back(controlPointLine);
			}
			CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
		
			MutiBezierPatch.push_back(splineSurface);
			
			idx++;
		}
	}
	//bsplineSurfaces = MutiBezierPatch;
	//saveIGESFile("./Resources/tmp_convert.iges");
	resetViewFlag = false;
	loadIGESSurface(MutiBezierPatch);
	resetViewFlag = true;
	cout << "HexMesh convert to surface mesh success.\n";
	return true;
}

bool ModelDraw::hexMeshConvertToIGES_2(std::string model_name) {
    if (!this->hexMesh) return false;

    std::string filename = model_name + ".mesh";

    /*MeshKernel::IO io;
    io.iGameWriteVolumeFile(*hexMesh, filename);*/

    HexV3fMesh* _hexMesh = new HexV3fMesh();

    /* iGame_hexMesh to XIE_hexMesh */
    {
        uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();
        for (uint64_t i = 0; i < vcnt; ++i) {
            MeshKernel::iGameVertexHandle vh(i);
            auto& v = hexMesh->vertices(vh);
            V3f _v(v.x(), v.y(), v.z());
            _hexMesh->add_vertex(_v);
        }

        for (uint64_t i = 0; i < ccnt; ++i) {
            MeshKernel::iGameCellHandle ch(i);
            if (hexMesh->isValid(ch)) {
                auto& c = hexMesh->cells(ch);
                auto vhs = c.getVertexHandle();
                vector<VertexHandle> _vhs;
                for (auto& vh : vhs) {
                    _vhs.push_back(VertexHandle(vh));
                }
                _hexMesh->add_cell(_vhs);
            }
        }

        /*std::cout << "iGameHexMesh convert to OldHexMesh success. OldHexMesh: #V = " << _hexMesh->vertices_.size() <<
            ", #C = " << _hexMesh->cells_.size() << std::endl;*/

    }

    std::cout << "trans success" << std::endl;
    /*FileHandler file(filename);
    assert(file.read2HexMeshObj(_hexMesh, true, true, true));*/

    if (_hexMesh->vertices_.size() == 0) {
        QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
        return false;
    }

    _hexMesh->examine_bdy();
    _hexMesh->cal_cell_cen();
    _hexMesh->cal_cen();
    _hexMesh->examine_cs();
    //_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
    //_hexMesh->smooth_mesh(); _hexMesh->smooth_mesh(); _hexMesh->smooth_mesh();
    std::cout << "Converting... hexMesh: vertices = " << _hexMesh->vertices_.size() << ", cells = " << _hexMesh->cells_.size() << std::endl;
    int udegree = 3;
    int vdegree = 3;
    vector<double> uknots{ 0,1 };
    vector<double> vknots{ 0,1 };
    vector<int> umultis{ 4,4 };
    vector<int> vmultis{ 4,4 };

    HexV3fMesh* SkelMesh = _hexMesh;
    std::vector<EdgeHandle> TempEdge;
    ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);

    int idx = 0;
    AppCC.approximationing();

    vector<CBSplineSurface> MutiBezierPatch;
    vector<CBSplineSurface> dx;

    for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
        if (SkelMesh->is_bdy((FaceHandle)i)) {
            vector<vector<CPoint>> controlPointsSurface;
            for (int j = 0; j < 4; ++j) {
                vector<CPoint> controlPointLine;
                for (int k = 0; k < 4; ++k) {
                    controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
                }
                controlPointsSurface.push_back(controlPointLine);
            }
            CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);

            MutiBezierPatch.push_back(splineSurface);

            idx++;
        }
    }
    //bsplineSurfaces = MutiBezierPatch;
    //saveIGESFile("./Resources/tmp_convert.iges");
    resetViewFlag = false;
    loadIGESSurface(MutiBezierPatch);
    resetViewFlag = true;
    cout << "HexMesh convert to surface mesh success.\n";
    return true;
}

bool ModelDraw::hexMeshConvertToIGES_skeletonExtraction() {
	//if (meshFileName == "" || igesFileName == "") return false;

	////std::string InputModel = "C:/iGameLab/iGame-Poster/iGame-Poster/Model/Beam_211106.mesh";
	////std::string OutputModel = "C:/iGameLab/iGame-Poster/iGame-Poster/Model/Beam_Test.iges";

	///*std::string InputModel = argv[1];
	//std::string OutputModel = argv[2];*/
	//std::string InputModel = "Beam.mesh";
	//std::string OutputModel = "tmp.iges";

	///*------------------------------ Step1:  InputModel ------------------------------------------*/
	//HexV3fMesh* mesh = new HexV3fMesh();
	//FileHandler file(InputModel);
	//file.read2HexMeshObj(mesh, true, true, true);       //read InputModel
	//mesh->examine_bdy();
	//mesh->cal_cell_cen();
	//mesh->cal_cen();
	//mesh->examine_cs();
	//mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	//mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	////mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	////mesh->examine_bdy();
	////mesh->cal_cen();
	////mesh->examine_bdfv();
	////mesh->cal_cell_cen();
	////mesh->cal_bdnorm();
	///*if (InputModel != "Beam.mesh" || mesh->vertices_.size() != 14553) {
	//	mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	//	mesh->smooth_mesh(); mesh->smooth_mesh(); mesh->smooth_mesh();
	//	mesh->examine_bdy();
	//	mesh->cal_cen();
	//	mesh->examine_bdfv();
	//	mesh->cal_cell_cen();
	//	mesh->cal_bdnorm();
	//}*/
	///*------------------------------ Step2: Skeleton Extraction(CGAL) ------------------------------------------*/

	////transform Hex Mesh to Triangle Mesh
	//Transform trans(mesh);
	//trans.OutputPath = "triangle.off";
	//trans.Quad2Triangle();

	////extract point line skeleton
	//CgalMcfskel skeler;
	//skeler.InputPath = "triangle.off";
	//skeler.OutputPath = "PointLine.skel";
	//skeler.Run();


	//if (InputModel == "Beam.mesh" && mesh->vertices_.size() == 14553) {
	//	std::ofstream out("PointLine.skel");
	//	out << skelData;
	//	out.close();
	//}

	////generate Hex skeleton
	//GenHexSkel("PointLine.skel");        //Output HexSkel.mesh
	//if (InputModel == "Beam.mesh" && mesh->vertices_.size() == 14553) {
	//	std::ofstream out("HexSkel.mesh");
	//	out << HexData;
	//	out.close();
	//}



	///*------------------------------ Step3: CC Subdivision Fitting -------------------------------------*/
	//HexV3fMesh* SkelMesh = new HexV3fMesh();
	//FileHandler file1("HexSkel.mesh");
	//file1.read2HexMeshObj(SkelMesh, true, true, true);
	//SkelMesh->examine_bdy();
	//SkelMesh->cal_cell_cen();
	//SkelMesh->cal_cen();

	//std::vector<V3f> VtoFit(SkelMesh->vertices_.size());
	//std::vector<V3f> VRight(SkelMesh->vertices_.size());

	//if (InputModel == "Beam.mesh" && mesh->vertices_.size() == 14553) {
	//	std::stringstream ssVtoFit(ppvToFit);
	//	std::stringstream ssVRight(ppvRight);
	//	for (int i = 0; i < SkelMesh->vertices_.size(); ++i) {
	//		ssVtoFit >> VtoFit[i].x >> VtoFit[i].y >> VtoFit[i].z;
	//		ssVRight >> VRight[i].x >> VRight[i].y >> VRight[i].z;
	//	}
	//} else {
	//	for (int i = 0; i < SkelMesh->vertices_.size(); ++i) {
	//		VtoFit[i] = SkelMesh->vertices_[i];
	//		double MinError = 99999999;
	//		int jIdx = -1;
	//		for (int j = 0; j < mesh->vertices_.size(); ++j) {
	//			double CurError = distance(VtoFit[i], mesh->vertices_[j]);
	//			if (CurError < MinError) {
	//				MinError = CurError;
	//				jIdx = j;
	//			}
	//		}
	//		VRight[i] = mesh->vertices_[jIdx];
	//	}
	//}


	//std::vector<EdgeHandle> TempEdge;
	//PartPIA pia(SkelMesh, VRight, VtoFit, TempEdge);
	//pia.PartUpDate();



	//SurfaceFitter sf(SkelMesh, mesh);
	//sf.CCSubdivisionFitting();

	//// 输出最终网格
	//FileHandler file2("FinalHex.mesh");
	//file2.writeFile("FinalHex.mesh", SkelMesh);

	///*------------------------------ step4: OutputIGES(Opencascade) ------------------------------------------*/
	////std::cout << "ok1" << std::endl;
	//int udegree = 3;
	//int vdegree = 3;
	//vector<double> uknots{ 0,1 };
	//vector<double> vknots{ 0,1 };
	//vector<int> umultis{ 4,4 };
	//vector<int> vmultis{ 4,4 };
	//vector<CBSplineSurface> MutiBezierPatch;
	////HexV3fMesh* SkelMesh = mesh;
	////std::vector<EdgeHandle> TempEdge;
	//ApproximationCC AppCC(SkelMesh, SkelMesh, TempEdge);


	//int idx = 0;
	//AppCC.approximationingSharp();
	//for (int i = 0; i < SkelMesh->faces_.size(); ++i) {
	//	if (SkelMesh->is_bdy((FaceHandle)i)) {
	//		vector<vector<CPoint>> controlPointsSurface;
	//		for (int j = 0; j < 4; ++j) {
	//			vector<CPoint> controlPointLine;
	//			for (int k = 0; k < 4; ++k) {
	//				controlPointLine.push_back(CPoint(AppCC.beiv[idx][j][k].x, AppCC.beiv[idx][j][k].y, AppCC.beiv[idx][j][k].z));
	//			}
	//			controlPointsSurface.push_back(controlPointLine);
	//		}
	//		CBSplineSurface splineSurface(udegree, vdegree, controlPointsSurface, uknots, vknots, umultis, vmultis);
	//		MutiBezierPatch.push_back(splineSurface);
	//		idx++;
	//	}
	//}
	//CAidFunction::getIGESSurface(OutputModel, MutiBezierPatch);
	///*------------------------------ Others: delete intermediate file ------------------------------------------*/
	///*remove("HexSkel.mesh");
	//remove("triangle.off");
	//remove("PointLine.skel");*/

	//loadIGESSurface(OutputModel);// added by Sheng
	return false;
}

bool ModelDraw::loadIGESCurve(std::string fileName) {
	ModelClear();

	vector<CBSplineCurve> geomcurve;
	if (CAidFunction::readIGESGeomCurve(fileName, geomcurve)) {
		if (geomcurve.empty()) {
			return false;
		}
		bsplineCurves.clear();
		vector<CBSplineCurve>().swap(bsplineCurves);
		for (int i = 0; i < geomcurve.size(); i++) {
			bsplineCurves.push_back(geomcurve[i]);
		}
		curvesDraw.clear();
		vector<CBSplineCurve>().swap(curvesDraw);
		curvesDraw = bsplineCurves;
		bsplineCurves.clear();
		bsplineCurveToDiscrete();
		updateCurveModel();
		SetDrawModule(BSPLINECURVE);
		update();
		return true;
	}
	return false;
}

bool ModelDraw::loadIGESSurface(std::string fileName) {
	ModelClear();

	vector<CBSplineSurface>geomsurface;

	if (CAidFunction::readIGESGeomSurface(fileName, geomsurface)) {
		if (geomsurface.empty()) {
			return false;
		}
		bsplineSurfaces.clear();
		vector<CBSplineSurface>().swap(bsplineSurfaces);
		for (int i = 0; i < geomsurface.size(); i++) {
			bsplineSurfaces.push_back(geomsurface[i]);
		}
		surfacesDraw.clear();
		vector<CBSplineSurface>().swap(surfacesDraw);
		surfacesDraw = bsplineSurfaces;
		bsplineSurfaceToDiscrete();
		//std::cout << "bsplineSurfaceToDiscrete ok" << std::endl;
		controlpointsToDiscrete();
		//std::cout << "controlpointsToDiscrete ok" << std::endl;
		if (resetViewFlag) updateSurfaceModel();
		SetDrawModule(BSPLINESURFACE);
		updateDrawData();
		update();
		return true;
	}
	return false;
}

bool ModelDraw::loadIGESSurface(std::vector<CBSplineSurface>& geomsurface) {
	ModelClear();
	if (geomsurface.empty()) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}

	bsplineSurfaces = geomsurface;

	surfacesDraw = bsplineSurfaces;

	bsplineSurfaceToDiscrete();


	std::cout << "bsplineSurfaceToDiscrete done" << std::endl;
	controlpointsToDiscrete();

	std::cout << "controlpointsToDiscrete done" << std::endl;


	if (resetViewFlag) updateSurfaceModel();
	SetDrawModule(BSPLINESURFACE);
	std::cout << "SetDrawModule done" << std::endl;
	updateDrawData();
	std::cout << "updateDrawData done" << std::endl;
	update();
	return true;
}

bool ModelDraw::loadIGESSurfaceColor(std::vector<CBSplineSurface>& geomsurface, std::vector<CBSplineSurface>& geomsurface_color)
{
    ModelClear();
	if (geomsurface.empty() || geomsurface_color.empty()) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Surface is empty! "));
		return false;
	}
	surfacesDraw = geomsurface_color;
	bsplineSurfaceToDiscrete();
	vector<Vector3d> colorInsert = { Vector3d(0,0,153) ,Vector3d(0,76,153) ,Vector3d(0,153,153) ,Vector3d(0,153,76) ,Vector3d(0,153,0),
		Vector3d(76,153,0),Vector3d(153,153,0), Vector3d(153,76,0), Vector3d(153,0,0) };

	ColorBar dirColorBar(colorInsert);
	for (int i = 0; i < data_surfaces.size(); ++i) {
		auto& singleSurface = data_surfaces[i];
		vector<vector<Vector3d>> singleSurfaceColor[3];
		for (int j = 0; j < singleSurface.size(); j++) {
			vector<Vector3d> tmpColorX; vector<Vector3d> tmpColorY; vector<Vector3d> tmpColorZ;
			for (int k = 0; k < singleSurface[j].size(); ++k) {
				tmpColorX.emplace_back(dirColorBar.getColor(data_surfaces[i][j][k].getX()));
				tmpColorY.emplace_back(dirColorBar.getColor(data_surfaces[i][j][k].getY()));
				tmpColorZ.emplace_back(dirColorBar.getColor(data_surfaces[i][j][k].getZ()));
			}
			singleSurfaceColor[0].emplace_back(tmpColorX);
			singleSurfaceColor[1].emplace_back(tmpColorY);
			singleSurfaceColor[2].emplace_back(tmpColorZ);
		}
		data_surfaces_color[0].emplace_back(singleSurfaceColor[0]);
		data_surfaces_color[1].emplace_back(singleSurfaceColor[1]);
		data_surfaces_color[2].emplace_back(singleSurfaceColor[2]);
	}

	bsplineSurfaces = geomsurface;
	surfacesDraw = bsplineSurfaces;
	SetDrawModule(BSPLINESURFACE);
	bsplineSurfaceToDiscrete();
	//controlpointsToDiscrete();
	bshowWireLine = true;
	have_light = false;
	bsplineSurfaceWireLineToDiscrete();
	if (resetViewFlag) updateSurfaceModel();
	updateDrawData_Surface_IGES_Color();
	update();
	return true;
}

bool ModelDraw::saveIGESFile(std::string fileName) {
	if (!bsplineCurves.empty()) {
		CAidFunction::getIGESCurve(fileName, bsplineCurves);
	}
	else if (!bsplineSurfaces.empty()) {
		CAidFunction::getIGESSurface(fileName, bsplineSurfaces);
	}
	else {
		QMessageBox::information(this, QString::fromLocal8Bit(" 警告 "), QString::fromLocal8Bit(" 没有保存任何的模型 "));
		return false;
	}
	return true;
}

void ModelDraw::project() {
	projectSkel->hexMesh = hexMesh;
	projectSkel->projectVolumeToOff_sec1();
	projectSkel->movePoint1();
	updateDrawData_HexMesh();
	update();
	iGameDFS();
}

void ModelDraw::rselectALine() {
	if (selected_c_indices.size() == 1) {

		iGameVertex tmp(0.0, 0.0, 0.0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices.back());
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			tmp += hexMesh->vertices(vh);
		}
		tmp = tmp / 8.0;
		// check with all node position
		BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(tmp[0], tmp[1], tmp[2]));
		std::vector<BMesh::Node*> ret;
		solutionSkel->getALine(ret);
		for (int i = 0; i < ret.size(); i++) {
			iGameVertex temp_pick(ret[i]->getPosition().x, ret[i]->getPosition().y, ret[i]->getPosition().z);
			int csize = hexMesh->csize();
			double mindis = 999; int mineid = 0;
			for (int i = 0; i < csize; i++) {
				iGameVertex meshv(0, 0, 0);
				auto& cell = hexMesh->cells((iGameCellHandle)i);
				auto vhs = cell.getVertexHandle();
				for (auto& vh : vhs) {
					meshv += hexMesh->vertices(vh);
				}
				meshv = meshv / 8;
				double dis = meshv.dist(temp_pick);
				if (dis < mindis) {
					mindis = dis;
					mineid = i;
				}
			}
			// find it
			std::vector<iGameCellHandle>::iterator it;
			if ((std::find(selected_c_indices.begin(), selected_c_indices.end(), iGameCellHandle(mineid))) != selected_c_indices.end()) {
				// pickcell.erase(it);
			}
			else {
				selected_c_indices.push_back(iGameCellHandle(mineid));
			}
		}
	}
	updateMeshModel();
	updateDrawData();
	update();
}

void ModelDraw::setRoot(BMesh::Node* _root) {
	root = _root;
}

void ModelDraw::setNodeVec(std::vector<BMesh::Node*>& v) {
	nodeVec = v;
}

void ModelDraw::setNodeVec(const std::vector<BMesh::Node*>& v) {
	nodeVec = v;
}

void ModelDraw::rconnectNodes() {
	std::vector<BMesh::Node*> vv;
	for (int i = 0; i < selected_c_indices.size(); i++) {
		iGameVertex meshv(0, 0, 0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[i]);
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = hexMesh->vertices(vh);
			meshv = meshv + v;
		}
		BMesh::Node* node1 = solutionSkel->selectWhichNode(V3f(meshv[0], meshv[1], meshv[2]));
		if (node1 != nullptr) {
			vv.push_back(node1);
		}
	}
	// connect vv
	for (int i = 0; i < vv.size() - 1 && vv.size(); i++) {
		for (int j = i + 1; j < vv.size(); j++) {
			solutionSkel->smartConnectNodes(vv[i], vv[j]);
		}
	}
	setNodeVec(solutionSkel->getNodeVec());
	update();
}

void ModelDraw::rdisconnectNodes() {
	std::vector<BMesh::Node*> vv;
	for (int i = 0; i < selected_c_indices.size(); i++) {
		iGameVertex meshv(0, 0, 0);
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[i]);
		auto vhs = cell.getVertexHandle();
		for (auto& vh : vhs) {
			meshv = meshv + hexMesh->vertices(vh);
		}
		meshv = meshv / 8.0;
		BMesh::Node* node1 = solutionSkel->selectWhichNode(V3f(meshv[0], meshv[1], meshv[2]));
		if (node1 != nullptr) {
			vv.push_back(node1);
		}
	}
	// connect vv
	for (int i = 0; i < vv.size() - 1 && vv.size(); i++) {
		for (int j = i + 1; j < vv.size(); j++) {
			solutionSkel->smartConnectNodes(vv[i], vv[j]);
		}
	}
	setNodeVec(solutionSkel->getNodeVec());
	update();
}

void ModelDraw::BFS(BMesh::Node* root) {
	Q.push(root);
	V3f v1;
	V3f v2;
	glLineWidth(3);
	shaderProgram.setUniformValue("colormode", 235);
	while (!Q.empty()) {
		BMesh::Node* a = Q.front();
		visited[a->id] = true;
		v1.x = a->root->position.x;
		v1.y = a->root->position.y;
		v1.z = a->root->position.z;
		Q.pop();
		for (int i = 0; i < a->children.size(); i++) {
			if (visited[a->children[i]->id] == false) {
				Q.push(a->children[i]);
				v2.x = a->children[i]->root->position.x;
				v2.y = a->children[i]->root->position.y;
				v2.z = a->children[i]->root->position.z;
				glBegin(GL_LINES);

#if 1
				GLfloat color[] = { 0,0,1,0.5 }; // GLfloat color[] = { 0.5,0.5,0.5,0.1 }
				GLfloat shine[] = { 120.0f };
				glMaterialfv(GL_FRONT_AND_BACK, GL_AMBIENT, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_DIFFUSE, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SPECULAR, color);
				glMaterialfv(GL_FRONT_AND_BACK, GL_SHININESS, shine);
#endif
				glVertex3d(v1.x, v1.y, v1.z); glVertex3d(v2.x, v2.y, v2.z);
				glEnd();
			}
		}
	}
	glLineWidth(1);
}

void ModelDraw::BFSTraverse(BMesh::Node* root) {
	visited.clear();
	std::queue<BMesh::Node*> empty;
	std::swap(empty, Q);

	visited[root->id] = true;
	BFS(root);
	for (int i = 0; i < nodeVec.size(); i++) {
		if (!visited[nodeVec[i]->id]) {
			BFS(nodeVec[i]);
		}
	}
}

void ModelDraw::drawAxis() {
	double len = radius * 10;
	glLineWidth(2);
	shaderProgram.setUniformValue("colormode", 232);
	//glColor3f(1.0, 0.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, len, 0);
	glEnd();
	shaderProgram.setUniformValue("colormode", 233);
	//glColor3f(0.0, 1.0, 0.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(len, 0, 0);
	glEnd();
	shaderProgram.setUniformValue("colormode", 234);
	//glColor3f(0.0, 0.0, 1.0);
	glBegin(GL_LINES);
	glVertex3f(0, 0, 0);
	glVertex3f(0, 0, len);
	glEnd();
}

void ModelDraw::drawSurfaceMesh_skel() {
	glLineWidth(2);
	shaderProgram.setUniformValue("colormode", 236);
	glBegin(GL_LINES);
	for (auto e_it = projectSkel->triMeshtv.edges_begin(); e_it != projectSkel->triMeshtv.edges_end(); ++e_it)
	{
		// 得到边所代表的半边
		OpenMesh::HalfedgeHandle heh1 = projectSkel->triMeshtv.halfedge_handle(*e_it, 0); // 默认一个方向的半边
		OpenMesh::Vec3d edgePoint(0, 0, 0);
		int edgePointsNumber = 0;
		OpenMesh::DefaultTraits::Point pointV = projectSkel->triMeshtv.point(projectSkel->triMeshtv.from_vertex_handle(heh1)); // 这条（半）边的起点
		OpenMesh::DefaultTraits::Point pointW = projectSkel->triMeshtv.point(projectSkel->triMeshtv.to_vertex_handle(heh1));   // 这条（半）边的终点
		glVertex3f(pointV[0], pointV[1], pointV[2]);
		glVertex3f(pointW[0], pointW[1], pointW[2]);
	}
	glEnd();
}

void ModelDraw::drawSkel() {
	BFSTraverse(root);
}

void  ModelDraw::drawHorizon() {
	float stepd = radius * 0.05;
	shaderProgram.setUniformValue("colormode", 231);
	//glColor3f(164 / 256.0, 128 / 256.0, 251 / 256.0);
	const int step = 1;
	const int maxLen = 100;
	glLineWidth(1);
	glBegin(GL_LINES);
	for (int i = -maxLen; i < maxLen - step; i += step) {
		for (int j = -maxLen; j < maxLen - step; j += step) {
			glVertex3f(i * stepd, j * stepd, 0);
			glVertex3f(i * stepd, (j + step) * stepd, 0);
			glVertex3f(i * stepd, j * stepd, 0);
			glVertex3f((i + step) * stepd, j * stepd, 0);
		}
	}
	glEnd();
}

int ModelDraw::loadSkeletonLee(std::string filename) {
	ModelClear();
	if (!solutionSkel->readSkel(filename)) {
		std::cout << "[DEBUG] read skel error" << std::endl;
		return -1;
	}
	std::cout << "[DEBUG] read skel success" << std::endl;
	hexMesh = new MeshKernel::VolumeMesh();
	// first gen node Hex
	solutionSkel->genNodeHex();
	// let Bmesh hex to mesh
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	updateMeshModel();
	//updateDrawData_HexMesh();
	updateDrawData();
	update();
	return 0;
}

int ModelDraw::loadHexMesh(std::string fileName, bool reset_view) {
	if (compare_flag) {
		MeshKernel::IO io;
		MeshKernel::VolumeMesh* tmp_hexMesh = new MeshKernel::VolumeMesh(io.iGameReadVolumeFile_TopOptMeshGeneration(fileName));
		compared_vertices_data.clear();
		for (auto& vp : tmp_hexMesh->allvertices()) {
			if (tmp_hexMesh->NeighborCh(vp.first).empty()) continue;
			compared_vertices_data.push_back(vp.second.x());
			compared_vertices_data.push_back(vp.second.y());
			compared_vertices_data.push_back(vp.second.z());
			compared_vertices_data.push_back(0.1);
			compared_vertices_data.push_back(0.1);
			compared_vertices_data.push_back(0.8);
		}
		delete tmp_hexMesh;
		tmp_hexMesh = nullptr;
		if (hexMesh) {
			compared_bind_vertices();
			updateDrawData();
		}
		return 0;
	}
	ModelClear();
	MeshKernel::IO io;
	hexMesh = new MeshKernel::VolumeMesh(io.iGameReadVolumeFile(fileName));
	int exterior_vcnt = 0;
	//std::cout << "exterior_vcnt" << std::endl;
	/*for (const auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}*/
	//std::cout << "isOnBoundary ok" << std::endl;
	//bvh_tree->destoryBVH_Tree();
	//bvh_tree->buildBVH_Tree(*hexMesh);
	//printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(HEXMESH);
	
	if (reset_view && resetViewFlag) updateMeshModel();
	compared_bind_vertices();
	updateDrawData();
	update();
	return iGameDFS();
}

bool ModelDraw::saveHexMesh(std::string fileName) {

	MeshKernel::IO io;
	io.iGameWriteVolumeFile(*hexMesh, fileName);
	return true;
}

void ModelDraw::loadTetMesh(std::string fileName) {
	ModelClear();
	MeshKernel::IO io;
	tetMesh = new MeshKernel::TetMesh(io.iGameReadTetMeshFile(fileName));
	SetDrawModule(TETMESH);
	if (resetViewFlag) updateMeshModel();
	updateDrawData_TetMesh();
	update();
}


void ModelDraw::saveTetMesh(std::string fileName) {

	std::cout << "The function is not implemented yet." << std::endl;

}



void ModelDraw::drawHexMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
	//shaderProgram.setUniformValue("colormode", 606);
	glDepthRange(0.01, 1.0);
	shaderProgram.setUniformValue("colormode", -1);
	glPointSize(9);
	core->glDrawArrays(GL_POINTS, quad_cnt * 4, used_last_vcnt);
	core->glDrawArrays(GL_QUADS, 0, quad_cnt * 4);


	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		// 绘制线框
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		if (editmodel == EditModel::PICKITEM && pickType == PickType::PickEdge) {
			glLineWidth(meshline_width * 2.5f);
		}
		else {
			glLineWidth(meshline_width);
		}

		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < quad_cnt; ++i) {
			core->glDrawArrays(GL_LINE_LOOP, i * 4, 4);
		}

		{
			// 绘制约束
			shaderProgram.setUniformValue("colormode", 231);
			glDepthRange(0.01, 0.99999);
			glPointSize(7);
			glBegin(GL_POINTS);
			for (auto& fp : spc_id) {
				VH vh(fp.first);
				if (!hexMesh->isValid(vh)) continue;
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
			glEnd();

			shaderProgram.setUniformValue("colormode", 232);
			glBegin(GL_POINTS);
			for (auto& fp : force_id) {
				VH vh(fp.first);
				if (!hexMesh->isValid(vh)) continue;
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
			glEnd();

		}


	}
	// 请在绘制其他物体后再绘制半透明物体
	glDepthRange(0.01, 2.0);
	shaderProgram.setUniformValue("colormode", 518);// 半透明
	core->glDrawArrays(GL_QUADS, quad_cnt * 4 + used_last_vcnt, translucence_vertices_data.size() / 6);

	update();
}

void ModelDraw::drawPolyMesh_OPT() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	if (axis.empty()) {
		shaderProgram.setUniformValue("colormode", -1);
		auto total_vcnt = vertices_data.size() / 6;
		auto f_vcnt = poly_mesh->n_faces();
		core->glDrawArrays(GL_QUADS, 0, f_vcnt * 4);
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width * 2.f);
		GLsizei l_vcnt = total_vcnt - f_vcnt;
		core->glDrawArrays(GL_LINES, f_vcnt * 4, l_vcnt);// 绘制折线

		if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
			glEnable(GL_POLYGON_OFFSET_FILL);
			glPolygonOffset(1.0, 1.0);
			glLineWidth(meshline_width);
			shaderProgram.setUniformValue("colormode", 530);
			for (int i = 0; i < f_vcnt; ++i) {
				core->glDrawArrays(GL_LINE_LOOP, i * 4, 4);
			}
		}

	}
	else {
		// draw 
		shaderProgram.setUniformValue("colormode", 530);
		glLineWidth(meshline_width);
		glBegin(GL_LINES);
		for (auto& pos_pair : axis) {
			auto& p1 = pos_pair.first;
			auto& p2 = pos_pair.second;
			glVertex3f(p1.x(), p1.y(), 0);
			glVertex3f(p2.x(), p2.y(), 0);
		}
		glEnd();

	}

	update();
}

void ModelDraw::drawSurfaceMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	shaderProgram.setUniformValue("colormode", 1);
	if (analysis_flag) shaderProgram.setUniformValue("colormode", -1);
	core->glDrawArrays(GL_TRIANGLES, 0, tri_cnt * 3);
	core->glDrawArrays(GL_QUADS, tri_cnt * 3, quad_cnt * 4);


	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width);
		shaderProgram.setUniformValue("colormode", 530);
//		for (int i = 0; i < tri_cnt; ++i) {
//			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
//		}
		for (int i = 0; i < quad_cnt; ++i) {
			core->glDrawArrays(GL_LINE_LOOP, i * 4 + tri_cnt * 3, 4);
		}


	}

	{
		/* draw streamlines */
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width * 2);
		shaderProgram.setUniformValue("colormode", 999);
		glBegin(GL_LINES);
		for (auto& streamline : streamlines) {
			auto& points = streamline.points;
			for (int i = 0; i < points.size() - 1; ++i) {
				glVertex3d(points[i].x(), points[i].y(), points[i].z());
				glVertex3d(points[i + 1].x(), points[i + 1].y(), points[i + 1].z());
			}
		}
		if (!streamlines.empty()) {
			for (auto& ep : surface_mesh->alledges()) {
				if (surface_mesh->isOnBoundary(ep.first)) {
					auto& v1 = surface_mesh->vertices(ep.second.vh1());
					auto& v2 = surface_mesh->vertices(ep.second.vh2());
					glVertex3d(v1.x(), v1.y(), v1.z());
					glVertex3d(v2.x(), v2.y(), v2.z());
				}
			}
		}

		glEnd();

		/* draw singularity */
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v : singularity_vertices) {
			auto& pos = v.pos;
			glVertex3d(pos.x(), pos.y(), pos.z());
		}
		glEnd();


	}

	shaderProgram.setUniformValue("colormode", 888);
	glPointSize(7);
	glBegin(GL_POINTS);
	for (auto& fix_idx : fixed_indices) {
		auto& v = surface_mesh->vertices(MeshKernel::iGameVertexHandle(fix_idx));
		glVertex3d(v.x(), v.y(), v.z());
	}
	glEnd();

	update();
}

void ModelDraw::drawTetMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	if (cutFlag) shaderProgram.setUniformValue("colormode", -1);
	else shaderProgram.setUniformValue("colormode", 1);
	core->glDrawArrays(GL_TRIANGLES, 0, vertices_data.size() / 6);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(meshline_width);
		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < vertices_data.size() / 18; i++) {
			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
		}
	}

	shaderProgram.setUniformValue("colormode", 888);
	glPointSize(7);
	glBegin(GL_POINTS);
	for (auto& fix_idx : fixed_indices) {
		auto& v = tetMesh->vertices(MeshKernel::iGameVertexHandle(fix_idx));
		glVertex3d(v.x(), v.y(), v.z());
	}
	glEnd();


	update();
}

void ModelDraw::drawTriMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可

	shaderProgram.setUniformValue("colormode", 1);
	core->glDrawArrays(GL_TRIANGLES, 0, vertices_data.size() / 6);

	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
		glEnable(GL_POLYGON_OFFSET_FILL);
		glPolygonOffset(1.0, 1.0);
		glLineWidth(1.f);
		shaderProgram.setUniformValue("colormode", 530);
		for (int i = 0; i < vertices_data.size() / 18; i++) {
			core->glDrawArrays(GL_LINE_LOOP, i * 3, 3);
		}
	}

	update();
}

void ModelDraw::updateDrawData_Topopt() {
    vertices_data.clear();

    auto& cells = hexMesh->allcells();
    std::unordered_map<int, bool> current_drawed;
    for (auto& cp : cells) {

        bool drawFlag = false;
        auto vhs = cp.second.getVertexHandle();
        if (hex_density[cp.first.idx()] >= 0.69) drawFlag = true;

        if (drawFlag) {
            auto fhs = cp.second.getFaceHandle();
            for (auto& fh : fhs) {
                auto f_vhs = hexMesh->faces(fh).getVertexHandle();
                for (auto& f_vh : f_vhs) {
                    auto& f_v = hexMesh->vertices(f_vh);
                    current_drawed[f_vh] = true;
                    vertices_data.push_back(f_v.x());
                    vertices_data.push_back(f_v.y());
                    vertices_data.push_back(f_v.z());
                    if (hexMesh->isOnBoundary(fh)) {
                        vertices_data.push_back(0.349);
                        vertices_data.push_back(0.6);
                        vertices_data.push_back(0.349);
                    }
                    else {
                        vertices_data.push_back(0.6588);
                        vertices_data.push_back(0.6588);
                        vertices_data.push_back(0.0862);
                    }
                }
            }
        }
    }
    quad_cnt = vertices_data.size() / 24;
    int last_vcnt = compared_vertices_data.size() / 6;
    used_last_vcnt = 0;
    for (int i = 0; i < last_vcnt; ++i) {
        int idx = i * 6;
        if (current_drawed.count(compared_to_current[i])) {
            used_last_vcnt++;
            for (int j = 0; j < 6; ++j) {
                vertices_data.push_back(compared_vertices_data[j + idx]);
            }
        }
    }

    vertices_data.insert(vertices_data.end(), translucence_vertices_data.begin(), translucence_vertices_data.end());

    update();
}

void ModelDraw::updateDrawData_DensityChange() {
	vertices_data.clear();
	double cut_density = topopt_density;

    if(hexMesh == nullptr) return;
	auto& cells = hexMesh->allcells();
	std::unordered_map<int, bool> current_drawed;
	for (auto& cp : cells) {
		
		bool drawFlag = false;
		auto vhs = cp.second.getVertexHandle();
		if (hex_density[cp.first.idx()] >= topopt_density) drawFlag = true;

		if (drawFlag) {
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				auto f_vhs = hexMesh->faces(fh).getVertexHandle();
				for (auto& f_vh : f_vhs) {
					auto& f_v = hexMesh->vertices(f_vh);
					current_drawed[f_vh] = true;
					vertices_data.push_back(f_v.x());
					vertices_data.push_back(f_v.y());
					vertices_data.push_back(f_v.z());
					if (hexMesh->isOnBoundary(fh)) {
						vertices_data.push_back(0.349);
						vertices_data.push_back(0.6);
						vertices_data.push_back(0.349);
					}
					else {
						vertices_data.push_back(0.6588);
						vertices_data.push_back(0.6588);
						vertices_data.push_back(0.0862);
					}
				}
			}
		}
	}
	quad_cnt = vertices_data.size() / 24;
	int last_vcnt = compared_vertices_data.size() / 6;
	used_last_vcnt = 0;
	for (int i = 0; i < last_vcnt; ++i) {
		int idx = i * 6;
		if (current_drawed.count(compared_to_current[i])) {
			used_last_vcnt++;
			for (int j = 0; j < 6; ++j) {
				vertices_data.push_back(compared_vertices_data[j + idx]);
			}
		}
	}

	vertices_data.insert(vertices_data.end(), translucence_vertices_data.begin(), translucence_vertices_data.end());

	update();
}

void ModelDraw::updateDrawData_HexMesh_CutPlane() {
	if (!this->hexMesh) return;
	vertices_data.clear();

	double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];
	printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	auto& cells = hexMesh->allcells();
	std::unordered_map<int, bool> current_drawed;
	for (auto& cp : cells) {
		bool drawFlag = true;
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& vex = hexMesh->vertices(vh);
			if (vex.x() < draw_min_x || vex.y() < draw_min_y || vex.z() < draw_min_z) {
				drawFlag = false;
				break;
			}
		}
		if (drawFlag) {
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				auto f_vhs = hexMesh->faces(fh).getVertexHandle();
				for (auto& f_vh : f_vhs) {
					auto& f_v = hexMesh->vertices(f_vh);
					current_drawed[f_vh] = true;
					vertices_data.push_back(f_v.x());
					vertices_data.push_back(f_v.y());
					vertices_data.push_back(f_v.z());
					if (hexMesh->isOnBoundary(fh)) {
						vertices_data.push_back(0.349);
						vertices_data.push_back(0.6);
						vertices_data.push_back(0.349);
					}
					else {
						vertices_data.push_back(0.6588);
						vertices_data.push_back(0.6588);
						vertices_data.push_back(0.0862);
					}
				}
			}
		}
		else {
			if (!hexMesh->isOnBoundary(cp.first)) continue;
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				if (!hexMesh->isOnBoundary(fh)) continue;
				auto f_vhs = hexMesh->faces(fh).getVertexHandle();
				for (auto& f_vh : f_vhs) {
					auto& f_v = hexMesh->vertices(f_vh);
					translucence_vertices_data.push_back(f_v.x());
					translucence_vertices_data.push_back(f_v.y());
					translucence_vertices_data.push_back(f_v.z());
					translucence_vertices_data.push_back(0.7725);
					translucence_vertices_data.push_back(0.7725);
					translucence_vertices_data.push_back(0.7725);
				}
			}

		}

	}
	quad_cnt = vertices_data.size() / 24;
	int last_vcnt = compared_vertices_data.size() / 6;
	used_last_vcnt = 0;
	for (int i = 0; i < last_vcnt; ++i) {
		int idx = i * 6;
		if (current_drawed.count(compared_to_current[i])) {
			used_last_vcnt++;
			for (int j = 0; j < 6; ++j) {
				vertices_data.push_back(compared_vertices_data[j + idx]);
			}
		}
	}

	vertices_data.insert(vertices_data.end(), translucence_vertices_data.begin(), translucence_vertices_data.end());

	update();

}

void ModelDraw::updateDrawData_TetMesh_CutPlane() {
	if (!this->tetMesh) return;
	vertices_data.clear();

	double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];
	printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	auto& cells = tetMesh->allcells();
	for (auto& cp : cells) {
		bool drawFlag = true;
		auto vhs = cp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = tetMesh->vertices(vh);
			if (v.x() < draw_min_x || v.y() < draw_min_y || v.z() < draw_min_z) {
				drawFlag = false;
				break;
			}
		}
		if (drawFlag) {
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				auto f_vhs = tetMesh->faces(fh).getVertexHandle();
				for (auto& f_vh : f_vhs) {
					auto& f_v = tetMesh->vertices(f_vh);
					vertices_data.push_back(f_v.x());
					vertices_data.push_back(f_v.y());
					vertices_data.push_back(f_v.z());
					if (tetMesh->isOnBoundary(fh)) {
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
					}
					else {
						vertices_data.push_back(0.8078);
						vertices_data.push_back(0.7098);
						vertices_data.push_back(0.5098);
					}
				}
			}
		}

	}

	update();

}



void ModelDraw::updateDrawData_SurfaceMesh_CutPlane() {
	if (!this->surface_mesh || quad_cnt != 0) return;
	vertices_data.clear();
	surface_mesh->genAllVerticesNormal();
	double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];
	printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	auto& faces = surface_mesh->allfaces();
	for (auto& fp : faces) {
		bool drawFlag = true;
		auto vhs = fp.second.getVertexHandle();
		for (auto& vh : vhs) {
			auto& v = surface_mesh->vertices(vh);
			if (v.x() < draw_min_x || v.y() < draw_min_y || v.z() < draw_min_z) {
				drawFlag = false;
				break;
			}
		}
		if (drawFlag) {
			for (auto& f_vh : vhs) {
				auto& f_v = surface_mesh->vertices(f_vh);
				vertices_data.push_back(f_v.x());
				vertices_data.push_back(f_v.y());
				vertices_data.push_back(f_v.z());
				vertices_data.push_back(f_v.getNormalX());
				vertices_data.push_back(f_v.getNormalY());
				vertices_data.push_back(f_v.getNormalZ());
			}
		}

	}
	update();

}

void ModelDraw::updateMeshModel() {// 一般只在重新加载模型的时候使用
	if (drawmodule == SURFACEMESH) {
		assert(surface_mesh != nullptr);
		this->surface_mesh->initBBox();
		this->bbMin = Vector3d(surface_mesh->BBoxMin.x(), surface_mesh->BBoxMin.y(), surface_mesh->BBoxMin.z());
		this->bbMax = Vector3d(surface_mesh->BBoxMax.x(), surface_mesh->BBoxMax.y(), surface_mesh->BBoxMax.z());
	}
	else if (drawmodule == HEXMESH || drawmodule == SKEL) {
		assert(hexMesh != nullptr);
		this->hexMesh->initBBox();
		this->bbMin = Vector3d(hexMesh->BBoxMin.x(), hexMesh->BBoxMin.y(), hexMesh->BBoxMin.z());
		this->bbMax = Vector3d(hexMesh->BBoxMax.x(), hexMesh->BBoxMax.y(), hexMesh->BBoxMax.z());
	}
	else if (drawmodule == TETMESH) {
		assert(tetMesh != nullptr);
		this->tetMesh->initBBox();
		this->bbMin = Vector3d(tetMesh->BBoxMin.x(), tetMesh->BBoxMin.y(), tetMesh->BBoxMin.z());
		this->bbMax = Vector3d(tetMesh->BBoxMax.x(), tetMesh->BBoxMax.y(), tetMesh->BBoxMax.z());
	}
	else if (drawmodule == SPARSEOPTIMIZATION) {
		if (axis.empty()) {
			double min_x = 999999999, min_y = 999999999, min_z = 999999999;
			double max_x = -999999999, max_y = -999999999, max_z = -999999999;
			assert(poly_mesh);
			poly_mesh->request_vertex_status();
			for (auto v_it = poly_mesh->vertices_begin(); v_it != poly_mesh->vertices_end(); ++v_it) {
				if (poly_mesh->status(*v_it).deleted()) continue;
				auto V = poly_mesh->point(*v_it);
				auto P = V.data();
				double x = P[0], y = P[1], z = P[2];
				min_x = std::min(x, min_x), min_y = std::min(y, min_y), min_z = std::min(z, min_z);
				max_x = std::max(x, max_z), max_y = std::max(y, max_y), max_z = std::max(z, max_z);
			}
			this->bbMin = Vector3d(min_x, min_y, min_z);
			this->bbMax = Vector3d(max_x, max_y, max_z);
		}
	}
	else if (drawmodule == TRIMESH) {
		double min_x = 999999999, min_y = 999999999, min_z = 999999999;
		double max_x = -999999999, max_y = -999999999, max_z = -999999999;
		tri_mesh->request_vertex_status();
		for (auto v_it = tri_mesh->vertices_begin(); v_it != tri_mesh->vertices_end(); ++v_it) {
			if (tri_mesh->status(*v_it).deleted()) continue;
			auto V = tri_mesh->point(*v_it);
			auto P = V.data();
			double x = P[0], y = P[1], z = P[2];
			min_x = std::min(x, min_x), min_y = std::min(y, min_y), min_z = std::min(z, min_z);
			max_x = std::max(x, max_z), max_y = std::max(y, max_y), max_z = std::max(z, max_z);
		}
		this->bbMin = Vector3d(min_x, min_y, min_z);
		this->bbMax = Vector3d(max_x, max_y, max_z);
	}
	else if (drawmodule == SKELMESH) {
		assert(skeletal_mesh != nullptr);
		this->skeletal_mesh->initBBox();
		this->bbMin = Vector3d(skeletal_mesh->BBoxMin.x(), skeletal_mesh->BBoxMin.y(), skeletal_mesh->BBoxMin.z());
		this->bbMax = Vector3d(skeletal_mesh->BBoxMax.x(), skeletal_mesh->BBoxMax.y(), skeletal_mesh->BBoxMax.z());
		std::cout << "Bounding Box Dig Length = " << (bbMax - bbMin).Length() << "\n";
	}


	Vector3d center = (bbMin + bbMax) / 2;
	double radius = (bbMax - bbMin).Length();
	SetScenePosition(center, radius);
	local_coordinater_axis->set_len(radius * 0.37);

}


void ModelDraw::drawNormalSkel_1() {
	ModelClear();
	hexMesh = new MeshKernel::VolumeMesh();
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	updateMeshModel();
	updateDrawData();
	update();
	return;
}

void ModelDraw::smoothing_hexahedron(int iter) {
	if (!this->hexMesh) return;

	vector<int> feature_chs = selected_c_indices;
	vector<int> feature_fhs = selected_f_indices;
	vector<int> feature_ehs = selected_e_indices;
	vector<int> feature_vhs = selected_v_indices;

	for (auto& vp : force_id) feature_vhs.push_back(vp.first);
	for (auto& vp : spc_id) feature_vhs.push_back(vp.first);

	HexMesh_Smoothing app(*this->hexMesh);
	app.Laplacian_Smoothing(iter, feature_chs, feature_fhs, feature_ehs, feature_vhs);

	updateDrawData();
	update();
}

void ModelDraw::subdivision(uint16_t iter) {

	if (!this->hexMesh) return;
	Subdivision_CC app(*this->hexMesh);
	app.subdivision_ordinary(iter);

	updateDrawData();
	update();

}

void ModelDraw::add_noise(double rate) {

	updateDrawData();
	update();

}


void ModelDraw::regularization_hexahedron(bool use_exterior_only, int max_iter) {
	if (!this->hexMesh) return;

	HexMesh_Regularization app(*this->hexMesh, use_exterior_only);
	app.regularization(max_iter);

	updateDrawData();
	update();
}


int ModelDraw::iGameDFS() {
	if (!hexMesh) return 0;
	int module_cnt = 0;
	unordered_map<int, bool> visited;
	for (auto& cp : hexMesh->allcells()) {
		if (visited.count(cp.first)) continue;
		module_cnt++;
		auto& ch = cp.first;
		queue<int> que;
		que.push(ch);
		visited[ch] = true;
		while (!que.empty()) {
			MeshKernel::iGameCellHandle cur_ch(que.front());
			que.pop();
			for (auto& adjch : hexMesh->NeighborCh(cur_ch)) {
				if (visited.count(adjch)) continue;
				que.push(adjch);
				visited[adjch] = true;
			}
		}
	}
	return module_cnt;
}


/**************************** B 样条曲面 ******************************************/

void ModelDraw::bsplineSurfaceToDiscrete() {
    time_t beg = clock();
    data_surfaces.clear();
    controlMesh.clear();
    double umin, umax, vmin, vmax;
    vector<CPoint> singlePatch; // 每个面片的离散点
    vector<vector<CPoint>> singlelSurface; // 每个曲面的离散点
    vector<vector<vector<CPoint>>> allSurface; // 所有曲面集合的离散点

    vector<CVector> singlePatch_normal;
    vector<vector<CVector>> singlelSurface_normal;
    vector<vector<vector<CVector>>> allSurface_normal;
    CVector normal; // 保存法向量

    int totalVerticesSize = 0;
    for (int i = 0; i < surfacesDraw.size(); i++) {
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
        int unumber, vnumber;
        unumber = upatchnum * (surfacesDraw[i].getSizeUControlPoints());
        vnumber = vpatchnum * (surfacesDraw[i].getSizeVControlPoints());
        unumber = unumber < 100 ? unumber : 100;
        vnumber = vnumber < 100 ? vnumber : 100;
        // For faster rendering, reduce the number of patches
        unumber *= 0.2;
        vnumber *= 0.2f;
        double sizeu = (umax - umin) / unumber - 1e-9; // 步距
        double sizev = (vmax - vmin) / vnumber - 1e-9;

        // Calculate the centroid of the surface for normal direction checking
        CPoint surfaceCentroid(0, 0, 0);
        int num = 0;
        for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u += sizeu) {
            for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v += sizev) {
                auto s = surfacesDraw[i].getPointAtParam(u, v);
                surfaceCentroid = s + surfacesDraw[i].getPointAtParam(u, v) + surfaceCentroid;
                num++;
            }
        }
        surfaceCentroid.setXYZ(surfaceCentroid.getX()/num, surfaceCentroid.getY()/num, surfaceCentroid.getZ()/num);

        for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u += sizeu) {
            for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v += sizev) {
                CPoint p0 = surfacesDraw[i].getPointAtParam(u, v);
                CAidFunction::geomNormal(surfacesDraw[i], u, v, normal);
                CPoint P0 = p0 - surfaceCentroid;
                if ((P0.getX() * normal.getX() + P0.getY() * normal.getY() + P0.getZ() * normal.getZ()) < 0) {
                    normal.setXYZ(normal.getX() * -1.0, normal.getY() * -1.0, normal.getZ() * -1.0);
                }
                singlePatch.push_back(p0);
                singlePatch_normal.push_back(normal);

                CPoint p1 = surfacesDraw[i].getPointAtParam(u, v + sizev);
                CAidFunction::geomNormal(surfacesDraw[i], u, v + sizev, normal);
                CPoint P1 = p1 - surfaceCentroid;
                if ((P1.getX() * normal.getX() + P1.getY() * normal.getY() + P1.getZ() * normal.getZ()) < 0) {
                    normal.setXYZ(normal.getX() * -1.0, normal.getY() * -1.0, normal.getZ() * -1.0);
                }
                singlePatch.push_back(p1);
                singlePatch_normal.push_back(normal);

                CPoint p2 = surfacesDraw[i].getPointAtParam(u + sizeu, v + sizev);
                CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v + sizev, normal);
                CPoint P2 = p2 - surfaceCentroid;
                if ((P2.getX() * normal.getX() + P2.getY() * normal.getY() + P2.getZ() * normal.getZ()) < 0) {
                    normal.setXYZ(normal.getX() * -1.0, normal.getY() * -1.0, normal.getZ() * -1.0);
                }
                singlePatch.push_back(p2);
                singlePatch_normal.push_back(normal);

                CPoint p3 = surfacesDraw[i].getPointAtParam(u + sizeu, v);
                CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v, normal);
                CPoint P3 = p3 - surfaceCentroid;
                if ((P3.getX() * normal.getX() + P3.getY() * normal.getY() + P3.getZ() * normal.getZ()) < 0) {
                    normal.setXYZ(normal.getX() * -1.0, normal.getY() * -1.0, normal.getZ() * -1.0);
                }
                singlePatch.push_back(p3);
                singlePatch_normal.push_back(normal);

                singlelSurface.push_back(singlePatch);
                singlelSurface_normal.push_back(singlePatch_normal);
                singlePatch.clear();
                singlePatch_normal.clear();
            }
        }

        allSurface.push_back(singlelSurface);
        if (!singlelSurface.empty()) {
            totalVerticesSize += singlelSurface.size() * singlelSurface[0].size();
        }
        singlelSurface.clear();
        allSurface_normal.push_back(singlelSurface_normal);
        singlelSurface_normal.clear();
    }

    this->data_surfaces = allSurface; // 曲面离散点集合
    this->surface_normals = allSurface_normal;
    this->data_bound_surface = getData_boundcurve(surfacesDraw); // 曲面边界曲线离散点集合

    // 获取控制点
    vector<vector<CPoint>> singleSurfaceControlPoints;
    vector<vector<vector<CPoint>>> allSurfaceControlPoints;
    for (int i = 0; i < surfacesDraw.size(); i++) {
        singleSurfaceControlPoints = surfacesDraw[i].getControlPoints();
        allSurfaceControlPoints.push_back(singleSurfaceControlPoints);
    }
    this->controlMesh = allSurfaceControlPoints;

    time_t end = clock();
    printf("bsplineSurfaceToDiscrete cost: %zdms, vertices size: %zd, surface size: %zd\n", int(end - beg), totalVerticesSize, allSurface.size());
    // updateDrawData();


//    time_t beg = clock();
//    data_surfaces.clear();
//    controlMesh.clear();
//
//    double umin, umax, vmin, vmax;
//    vector<CPoint> singlePatch; // 每个面片的离散点
//    vector<vector<CPoint>> singleSurface; // 每个曲面的离散点
//    vector<vector<vector<CPoint>>> allSurface; // 所有曲面集合的离散点
//
//    vector<CVector> singlePatch_normal;
//    vector<vector<CVector>> singleSurface_normal;
//    vector<vector<vector<CVector>>> allSurface_normal;
//    CVector normal; // 保存法向量
//
//    int totalVerticesSize = 0;
//    for (int i = 0; i < surfacesDraw.size(); i++) {
//        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
//        int unumber = upatchnum * surfacesDraw[i].getSizeUControlPoints();
//        int vnumber = vpatchnum * surfacesDraw[i].getSizeVControlPoints();
//        unumber = unumber < 100 ? unumber : 100;
//        vnumber = vnumber < 100 ? vnumber : 100;
//        unumber = static_cast<int>(unumber * 0.2);
//        vnumber = static_cast<int>(vnumber * 0.2f);
//        double sizeu = (umax - umin) / unumber - 1e-9; // 步距
//        double sizev = (vmax - vmin) / vnumber - 1e-9;
//
//        for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u += sizeu) {
//            for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v += sizev) {
//                singlePatch.push_back(surfacesDraw[i].getPointAtParam(u, v));
//                CAidFunction::geomNormal(surfacesDraw[i], u, v, normal);
//                singlePatch_normal.push_back(normal);
//
//                singlePatch.push_back(surfacesDraw[i].getPointAtParam(u, v + sizev));
//                CAidFunction::geomNormal(surfacesDraw[i], u, v + sizev, normal);
//                singlePatch_normal.push_back(normal);
//
//                singlePatch.push_back(surfacesDraw[i].getPointAtParam(u + sizeu, v + sizev));
//                CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v + sizev, normal);
//                singlePatch_normal.push_back(normal);
//
//                singlePatch.push_back(surfacesDraw[i].getPointAtParam(u + sizeu, v));
//                CAidFunction::geomNormal(surfacesDraw[i], u + sizeu, v, normal);
//                singlePatch_normal.push_back(normal);
//
//                singleSurface.push_back(singlePatch);
//                singleSurface_normal.push_back(singlePatch_normal);
//                singlePatch.clear();
//                singlePatch_normal.clear();
//            }
//        }
//
//        allSurface.push_back(singleSurface);
//        if (!singleSurface.empty()) totalVerticesSize += singleSurface.size() * singleSurface[0].size();
//        singleSurface.clear();
//        allSurface_normal.push_back(singleSurface_normal);
//        singleSurface_normal.clear();
//    }
//
//    this->data_surfaces = allSurface; // 曲面离散点集合
//    this->surface_normals = allSurface_normal;
//    this->data_bound_surface = getData_boundcurve(surfacesDraw); // 曲面边界曲线离散点集合
//
//// 获取控制点
//    vector<vector<CPoint>> singleSurfaceControlPoints;
//    vector<vector<vector<CPoint>>> allSurfaceControlPoints;
//    for (int i = 0; i < surfacesDraw.size(); i++) {
//        singleSurfaceControlPoints = surfacesDraw[i].getControlPoints();
//        allSurfaceControlPoints.push_back(singleSurfaceControlPoints);
//    }
//    this->controlMesh = allSurfaceControlPoints;
//
//    time_t end = clock();
//    printf("bsplineSurfaceToDiscrete cost: %zdms, vertices size: %zd, surface size: %zd\n", int(end - beg), totalVerticesSize, allSurface.size());

}

void ModelDraw::bsplineSurfaceWireLineToDiscrete() {
    data_surface_line.clear();
    double umin, umax, vmin, vmax;
    vector<vector<CPoint>> lines;

    int totalVerticesSize = 0;
    for (int i = 0; i < surfacesDraw.size(); i++) {
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
        int unumber, vnumber;
        unumber = upatchnum * (surfacesDraw[i].getSizeUControlPoints());
        vnumber = vpatchnum * (surfacesDraw[i].getSizeVControlPoints());
        unumber = unumber < 100 ? unumber : 100;
        vnumber = vnumber < 100 ? vnumber : 100;
        unumber *= 0.2, vnumber *= 0.2f;
        double sizeu = (umax - umin) / unumber - 1e-9;// 步距
        double sizev = (vmax - vmin) / vnumber - 1e-9;
        vector<CPoint> line;
        for (double u = umin + 1e-9; u <= umax - 1e-9; u = u + sizeu) {
            double v = vmin;
            line.push_back(surfacesDraw[i].getPointAtParam(u, v));
        }
        lines.push_back(line);
        line.clear();
        for (double u = umin + 1e-9; u <= umax - 1e-9; u = u + sizeu) {
            double v = vmax;
            line.push_back(surfacesDraw[i].getPointAtParam(u, v));
        }
        lines.push_back(line);
        line.clear();
        for (double v = vmin + 1e-9; v < vmax - 1e-9; v = v + sizev) {
            double u = umin;
            line.push_back(surfacesDraw[i].getPointAtParam(u, v));
        }
        lines.push_back(line);
        line.clear();
        for (double v = vmin + 1e-9; v < vmax - 1e-9; v = v + sizev) {
            double u = umax;
            line.push_back(surfacesDraw[i].getPointAtParam(u, v));
        }
        lines.push_back(line);
        line.clear();

    }
    this->data_surface_line = lines;
}

void ModelDraw::updateSurfaceModel() {
	if (surfacesDraw.empty()) {
		return;
	}
	bbMin = Vector3d(std::numeric_limits<double>::max(), std::numeric_limits<double>::max(), std::numeric_limits<double>::max());
	bbMax = Vector3d(-99999999, -99999999, -99999999);
	for (auto& vv_p : data_bound_surface) {
		for (auto& v_p : vv_p) {
			for (auto& p : v_p) {
				double x = p.getX(), y = p.getY(), z = p.getZ();
				bbMin[0] = std::min(bbMin[0], x), bbMin[1] = std::min(bbMin[1], y), bbMin[2] = std::min(bbMin[2], z);
				bbMax[0] = std::max(bbMax[0], x), bbMax[1] = std::max(bbMax[1], y), bbMax[2] = std::max(bbMax[2], z);
			}
		}
	}
	double u1, u2, v1, v2;
	surfacesDraw[0].getUVMinMax(u1, u2, v1, v2);
	double r = 1.5 * sqrt(
		pow((surfacecenter.getX() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getX()), 2) +
		pow((surfacecenter.getY() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getY()), 2) +
		pow((surfacecenter.getZ() - surfacesDraw[0].getPointAtParam(u1 + 1.e-9, v1 + 1.e-9).getZ()), 2)
	);
	SetScenePosition(Vector3d(surfacecenter.getX(), surfacecenter.getY(), surfacecenter.getZ()), r);
	update();
}

void ModelDraw::drawSurface() {
    core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
    // 设置面颜色
    shaderProgram.setUniformValue("colormode", -1);
    core->glDrawArrays(GL_QUADS, 0, vertices_data.size() / 6);
//
//    // 设置边界颜色
//    shaderProgram.setUniformValue("colormode", 530); // 使用黑色绘制边界
//    for (size_t i = 0; i < vertices_data.size(); i += 24) { // 4 vertices per quad, each vertex has 6 components
//        core->glDrawArrays(GL_LINE_LOOP, i / 6, 4); // 4 vertices per quad
//    }
//··················法向量绘制·····················
//    float normalLength = 1.1f; // 法向量的长度
//
//    for (size_t i = 0; i < data_surfaces.size(); ++i) {
//        auto& singleSurface = data_surfaces[i];
//        auto& normals = surface_normals[i];
//
////        glColor3f(1.0f, 1.0f, 1.0f); // 曲面的颜色
////        glBegin(GL_TRIANGLES);
////        for (size_t j = 0; j < singleSurface.size(); ++j) {
////            for (size_t k = 0; k < singleSurface[j].size(); ++k) {
////                auto& v = singleSurface[j][k];
////                glVertex3f(v.getX(), v.getY(), v.getZ());
////            }
////        }
////        glEnd();
//
//        // 绘制法向量
//        glColor3f(1.0f, 0.0f, 0.0f); // 法向量的颜色
//        glBegin(GL_LINES);
//        for (size_t j = 0; j < singleSurface.size(); ++j) {
//            for (size_t k = 0; k < singleSurface[j].size(); ++k) {
//                auto& v = singleSurface[j][k];
//                auto& n = normals[j][k];
//                glVertex3f(v.getX(), v.getY(), v.getZ());
//                glVertex3f(v.getX() + n.getX() * normalLength,
//                           v.getY() + n.getY() * normalLength,
//                           v.getZ() + n.getZ() * normalLength);
//            }
//        }
//        glEnd();
//    }
//··················九维渲染信息·····················
//    core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
//    if (have_light)
//        shaderProgram.setUniformValue("colormode", 1);
//    else
//        shaderProgram.setUniformValue("colormode", -1);
//    core->glDrawArrays(GL_QUADS, 0, vertices_data.size() / 9);

//    shaderProgram.setUniformValue("colormode", 530); // 使用黑色绘制边界
//    for (size_t i = 0; i < vertices_data.size(); i += 36) { // 4 vertices per quad, each vertex has 9 components
//        core->glDrawArrays(GL_LINE_LOOP, i / 9, 4); // 4 vertices per quad
//    }
}

void ModelDraw::drawSelectedObject() {

	shaderProgram.setUniformValue("colormode", 510);// Blue
	glDepthRange(0.01, 0.02);
	if (!selected_rectangle.empty() && editmodel == EditModel::PICKITEM) {
		glLineWidth(meshline_width);
		glBegin(GL_LINES);
		for (auto& vex : selected_rectangle) {
			glVertex3d(vex.x(), vex.y(), vex.z());
		}
		glEnd();
	}

    shaderProgram.setUniformValue("colormode", 237);// SpringGreen
    shaderProgram.setUniformValue("colormode", 238);// Salmon

	shaderProgram.setUniformValue("colormode", 999);// Red
	glDepthRange(0.01, 0.99999);
	switch (drawmodule) {

	case DrawModule::SURFACEMESH:
		drawSelectedObject_triMesh();
		break;
	case DrawModule::TETMESH:
		drawSelectedObject_tetMesh();
		break;
	case DrawModule::HEXMESH:
		drawSelectedObject_hexMesh();
		break;
	case DrawModule::SKELMESH:
		drawSelectedObject_sklMesh();
		break;
		/*Temp-lee*/
	case DrawModule::SKEL:
		drawSelectedObject_SKEL();
		break;
	default:
		break;
	}

	glDepthRange(0.01, 1.0);
	update();
}

void ModelDraw::drawSelectedObject_sklMesh() {
//	if (!skeletal_mesh) return;
//	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
//		glPointSize(7);
//		glBegin(GL_POINTS);
//		for (auto& v_idx : selected_v_indices) {
//			auto& v = skeletal_mesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
//			glVertex3d(v.x(), v.y(), v.z());
//		}
//		glEnd();
//	}
//	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
//		glLineWidth(1.5);
//		glBegin(GL_LINES);
//		for (auto& e_idx : selected_e_indices) {
//			auto& e = skeletal_mesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
//			auto& v1 = skeletal_mesh->vertices(e.vh1());
//			auto& v2 = skeletal_mesh->vertices(e.vh2());
//			glVertex3d(v1.x(), v1.y(), v1.z());
//			glVertex3d(v2.x(), v2.y(), v2.z());
//		}
//		glEnd();
//		//std::cout << "We are drawing edges.\n";
//	}
//	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
//		glBegin(GL_TRIANGLES);
//		for (auto& f_idx : selected_f_indices) {
//			auto vhs = skeletal_mesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
//			for (auto& vh : vhs) {
//				auto& v = skeletal_mesh->vertices(vh);
//				glVertex3d(v.x(), v.y(), v.z());
//			}
//		}
//		glEnd();
//	}

    if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
        glPointSize(17);
        glBegin(GL_POINTS);
        for (auto& v_idx : selected_v_indices) {
            auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
            glVertex3d(v.x(), v.y(), v.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
        glLineWidth(7);
        glBegin(GL_LINES);
        for (auto& e_idx : selected_e_indices) {
            auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
            auto& v1 = hexMesh->vertices(e.vh1());
            auto& v2 = hexMesh->vertices(e.vh2());
            glVertex3d(v1.x(), v1.y(), v1.z());
            glVertex3d(v2.x(), v2.y(), v2.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
        glBegin(GL_QUADS);
        for (auto& f_idx : selected_f_indices) {
            auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
            for (auto& vh : vhs) {
                auto& v = hexMesh->vertices(vh);
                glVertex3d(v.x(), v.y(), v.z());
            }
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
        glBegin(GL_QUADS);
        for (auto& c_idx : selected_c_indices) {
            auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
            for (auto& fh : fhs) {
                auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
                for (auto& vh : vhs) {
                    auto& v = hexMesh->vertices(vh);
                    glVertex3d(v.x(), v.y(), v.z());
                }
            }
        }
        glEnd();
    }
}

void ModelDraw::drawSelectedObject_SKEL() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(17);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(7);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = hexMesh->vertices(e.vh1());
			auto& v2 = hexMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& c_idx : selected_c_indices) {
			auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}

}

void ModelDraw::drawSelectedObject_hexMesh() {
    shaderProgram.setUniformValue("colormode", 237);// SpringGreen
    if ((drawSelectedAllFlag || pickType == PickPoint) && force_v_idx.size()) {
        glPointSize(17);
        glBegin(GL_POINTS);
        for (auto& v_idx : force_v_idx) {
            auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
            glVertex3d(v.x(), v.y(), v.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickEdge) && force_e_idx.size()) {
        glLineWidth(7);
        glBegin(GL_LINES);
        for (auto& e_idx : force_e_idx) {
            auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
            auto& v1 = hexMesh->vertices(e.vh1());
            auto& v2 = hexMesh->vertices(e.vh2());
            glVertex3d(v1.x(), v1.y(), v1.z());
            glVertex3d(v2.x(), v2.y(), v2.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickFace) && force_f_idx.size()) {
        glBegin(GL_QUADS);
        for (auto& f_idx : force_f_idx) {
            auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
            for (auto& vh : vhs) {
                auto& v = hexMesh->vertices(vh);
                glVertex3d(v.x(), v.y(), v.z());
            }
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickCell) && force_c_idx.size()) {
        glBegin(GL_QUADS);
        for (auto& c_idx : force_c_idx) {
            auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
            for (auto& fh : fhs) {
                auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
                for (auto& vh : vhs) {
                    auto& v = hexMesh->vertices(vh);
                    glVertex3d(v.x(), v.y(), v.z());
                }
            }
        }
        glEnd();
    }
    shaderProgram.setUniformValue("colormode", 238);// Salmon
    if ((drawSelectedAllFlag || pickType == PickPoint) && fix_v_idx.size()) {
        glPointSize(17);
        glBegin(GL_POINTS);
        for (auto& v_idx : fix_v_idx) {
            auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
            glVertex3d(v.x(), v.y(), v.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickEdge) && fix_e_idx.size()) {
        glLineWidth(7);
        glBegin(GL_LINES);
        for (auto& e_idx : fix_e_idx) {
            auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
            auto& v1 = hexMesh->vertices(e.vh1());
            auto& v2 = hexMesh->vertices(e.vh2());
            glVertex3d(v1.x(), v1.y(), v1.z());
            glVertex3d(v2.x(), v2.y(), v2.z());
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickFace) && fix_f_idx.size()) {
        glBegin(GL_QUADS);
        for (auto& f_idx : fix_f_idx) {
            auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
            for (auto& vh : vhs) {
                auto& v = hexMesh->vertices(vh);
                glVertex3d(v.x(), v.y(), v.z());
            }
        }
        glEnd();
    }
    if ((drawSelectedAllFlag || pickType == PickCell) && fix_c_idx.size()) {
        glBegin(GL_QUADS);
        for (auto& c_idx : fix_c_idx) {
            auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
            for (auto& fh : fhs) {
                auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
                for (auto& vh : vhs) {
                    auto& v = hexMesh->vertices(vh);
                    glVertex3d(v.x(), v.y(), v.z());
                }
            }
        }
        glEnd();
    }
    shaderProgram.setUniformValue("colormode", 999);// red
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(17);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = hexMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(7);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = hexMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = hexMesh->vertices(e.vh1());
			auto& v2 = hexMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& c_idx : selected_c_indices) {
			auto fhs = hexMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}

    shaderProgram.setUniformValue("colormode", 239);// Thistle1
    if ((drawSelectedAllFlag || pickType == PickFace) && RangeTodoFaces.size()) {
        glBegin(GL_QUADS);
        for (auto& f_idx : RangeTodoFaces) {
            auto vhs = hexMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
            for (auto& vh : vhs) {
                auto& v = hexMesh->vertices(vh);
                glVertex3d(v.x(), v.y(), v.z());
            }
        }
        glEnd();
    }
}

void ModelDraw::drawSelectedObject_tetMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = tetMesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(7);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = tetMesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = tetMesh->vertices(e.vh1());
			auto& v2 = tetMesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = tetMesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = tetMesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickCell) && selected_c_indices.size()) {
		glBegin(GL_QUADS);
		for (auto& c_idx : selected_c_indices) {
			auto fhs = tetMesh->cells(MeshKernel::iGameCellHandle(c_idx)).getFaceHandle();
			for (auto& fh : fhs) {
				auto vhs = tetMesh->faces(MeshKernel::iGameFaceHandle(fh)).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = tetMesh->vertices(vh);
					glVertex3d(v.x(), v.y(), v.z());
				}
			}
		}
		glEnd();
	}
}

void ModelDraw::drawSelectedObject_triMesh() {
	if ((drawSelectedAllFlag || pickType == PickPoint) && selected_v_indices.size()) {
		glPointSize(7);
		glBegin(GL_POINTS);
		for (auto& v_idx : selected_v_indices) {
			auto& v = surface_mesh->vertices(MeshKernel::iGameVertexHandle(v_idx));
			glVertex3d(v.x(), v.y(), v.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickEdge) && selected_e_indices.size()) {
		glLineWidth(1.5);
		glBegin(GL_LINES);
		for (auto& e_idx : selected_e_indices) {
			auto& e = surface_mesh->edges(MeshKernel::iGameEdgeHandle(e_idx));
			auto& v1 = surface_mesh->vertices(e.vh1());
			auto& v2 = surface_mesh->vertices(e.vh2());
			glVertex3d(v1.x(), v1.y(), v1.z());
			glVertex3d(v2.x(), v2.y(), v2.z());
		}
		glEnd();
	}
	if ((drawSelectedAllFlag || pickType == PickFace) && selected_f_indices.size()) {
		glBegin(GL_TRIANGLES);
		for (auto& f_idx : selected_f_indices) {
			auto vhs = surface_mesh->faces(MeshKernel::iGameFaceHandle(f_idx)).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = surface_mesh->vertices(vh);
				glVertex3d(v.x(), v.y(), v.z());
			}
		}
		glEnd();
	}
}

std::vector<vector<std::vector<CPoint>>> ModelDraw::getData_boundcurve(vector<CBSplineSurface> boundsurface) {
	vector<CPoint>singleBoundPoints;//一条边界上的离散点；
	vector<vector<CPoint>>singleSurfacePoints;//一张曲面上的离散点；
	vector<vector<vector<CPoint>>>allSurfacePoints;//所有边界离散点集合；
	CPoint p(0, 0, 0);
	int number = 0;
	for (int i = 0; i < boundsurface.size(); i++) {
		double umin, umax, vmin, vmax;
		boundsurface[i].getUVMinMax(umin, umax, vmin, vmax);
		int unumber, vnumber;
		unumber = upatchnum * (boundsurface[i].getSizeUControlPoints());
		vnumber = vpatchnum * (boundsurface[i].getSizeVControlPoints());
		number += 2 * (unumber + vnumber);
		double sizeu = (umax - umin) / unumber - 1.e-9;
		double sizev = (vmax - vmin) / vnumber - 1.e-9;
		for (double u = umin + 1.e-9; u <= umax - 1.e-9; u = u + sizeu) {//第一条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(u, vmin + 1.e-9));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		for (double u = umin + 1.e-9; u <= umax - 1.e-9; u = u + sizeu) {//第二条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(u, vmax - 1.e-9));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		for (double v = vmin + 1.e-9; v <= vmax - 1.e-9; v = v + sizev) {//第三条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(umin + 1.e-9, v));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		for (double v = vmin + 1.e-9; v <= vmax - 1.e-9; v = v + sizev) {//第四条边界
			singleBoundPoints.push_back(boundsurface[i].getPointAtParam(umax - 1.e-9, v));
			p.setXYZ(p.getX() + singleBoundPoints.back().getX(), p.getY() + singleBoundPoints.back().getY(), p.getZ() + singleBoundPoints.back().getZ());
		}
		singleSurfacePoints.push_back(singleBoundPoints);
		singleBoundPoints.clear();

		allSurfacePoints.push_back(singleSurfacePoints);
		singleSurfacePoints.clear();
	}
	surfacecenter.setXYZ(p.getX() / number, p.getY() / number, p.getZ() / number);
	return allSurfacePoints;
}

/**************************** B 样条曲线 ******************************************/

void ModelDraw::bsplineCurveToDiscrete() {
	data_curves.clear();
	controlLine.clear();
	double umin, umax;
	vector<CPoint>singleCurve;//每个曲线的离散点
	vector<vector<CPoint>>allCurve;//所有曲线的离散点
	CPoint p(0, 0, 0);
	int number = 0;
	for (int i = 0; i < curvesDraw.size(); i++) {
		umax = curvesDraw[i].getUmax();
		umin = curvesDraw[i].getUmin();
		int unumber;
		unumber = upatchnum * curvesDraw[i].getSizeControlPoints();
		number += unumber;
		double sizeu = (umax - umin) / unumber - 1.e-9;
		for (double u = umin + 1.e-9; u <= umax - 1.e-9; u = u + sizeu) {
			singleCurve.push_back(curvesDraw[i].getPointAtParam(u));
			p.setXYZ(p.getX() + singleCurve.back().getX(), p.getY() + singleCurve.back().getY(), p.getZ() + singleCurve.back().getZ());
		}
		allCurve.push_back(singleCurve);
		singleCurve.clear();
	}
	curvecenter.setXYZ(p.getX() / number, p.getY() / number, p.getZ() / number);
	this->data_curves = allCurve;
	//获取控制点
	vector<CPoint>singleCurveControlPoint;
	for (int i = 0; i < curvesDraw.size(); i++) {
		singleCurveControlPoint = curvesDraw[i].getControlPoints();
		controlLine.push_back(singleCurveControlPoint);
	}
	data_controlLine.clear();
	data_controlLine = controlLine;
}

void ModelDraw::updateCurveModel() {
	if (curvesDraw.empty()) {
		return;
	}
	double umax, umin;
	umax = curvesDraw[0].getUmax();
	umin = curvesDraw[0].getUmin();
	double r = 1.5 * sqrt(
		pow((curvecenter.getX() - curvesDraw[0].getPointAtParam(umin + 1.e-9).getX()), 2) +
		pow((curvecenter.getY() - curvesDraw[0].getPointAtParam(umin + 1.e-9).getY()), 2) +
		pow((curvecenter.getZ() - curvesDraw[0].getPointAtParam(umin + 1.e-9).getZ()), 2)
	);
	SetScenePosition(Vector3d(curvecenter.getX(), curvecenter.getY(), curvecenter.getZ()), r);
	update();
}

void ModelDraw::drawCurve(void) {
	for (int i = 0; i < data_curves.size(); i++) {
		glDepthRange(0.01, 1.0);
		glDepthFunc(GL_LEQUAL);
		glColor4f(0.5, 0.9, 0.2, 0.0);
		glLineWidth(2.0);
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < data_curves[i].size(); j++) {
			glVertex3d(data_curves[i][j].getX(), data_curves[i][j].getY(), data_curves[i][j].getZ());
		}
		glEnd();
	}
}

/**************************** 几何公共函数 ******************************************/

void ModelDraw::drawControlPoint(double xx, double yy, double zz, double radius, GLfloat M, GLfloat N) {
	float step_z = PI / M;
	float step_xy = 2 * PI / N;
	double x[4], y[4], z[4];
	float angle_z = 0.0;
	float angle_xy = 0.0;
	int i = 0, j = 0;
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	glBegin(GL_QUADS);
	for (i = 0; i < M; i++) {
		angle_z = i * step_z;

		for (j = 0; j < N; j++) {
			angle_xy = j * step_xy;

			x[0] = radius * sin(angle_z) * cos(angle_xy);
			y[0] = radius * sin(angle_z) * sin(angle_xy);
			z[0] = radius * cos(angle_z);

			x[1] = radius * sin(angle_z + step_z) * cos(angle_xy);
			y[1] = radius * sin(angle_z + step_z) * sin(angle_xy);
			z[1] = radius * cos(angle_z + step_z);

			x[2] = radius * sin(angle_z + step_z) * cos(angle_xy + step_xy);
			y[2] = radius * sin(angle_z + step_z) * sin(angle_xy + step_xy);
			z[2] = radius * cos(angle_z + step_z);

			x[3] = radius * sin(angle_z) * cos(angle_xy + step_xy);
			y[3] = radius * sin(angle_z) * sin(angle_xy + step_xy);
			z[3] = radius * cos(angle_z);
			for (int k = 0; k < 4; k++) {
				double distance = sqrt(pow(x[k], 2) + pow(y[k], 2) + pow(z[k], 2));
				glNormal3d(x[k] / distance, y[k] / distance, z[k] / distance);
				glVertex3d(xx + x[k], yy + y[k], zz + z[k]);
			}
		}
	}
	glEnd();

}

void ModelDraw::controlpointsToDiscrete() {

	data_controlMesh.clear();
	vector<CPoint>singlePatch;//每个面片的离散点
	vector<vector<CPoint>>singlelSurface;//每个曲面的离散点
	vector<vector<vector<CPoint>>>allSurface;//所有曲面集合的离散点
	for (int i = 0; i < controlMesh.size(); i++) {//controlMesh[i]是第i个曲面
		for (int j = 0; j < controlMesh[i].size() - 1; j++) {//controlMesh[i][j]是第i个曲面控制点集合的第j行
			for (int k = 0; k < controlMesh[i][j].size() - 1; k++) {//controlMesh[i][j][k]是第i个曲面的第j行第k列的控制点
				singlePatch.push_back(controlMesh[i][j][k]);
				singlePatch.push_back(controlMesh[i][j][k + 1]);
				singlePatch.push_back(controlMesh[i][j + 1][k + 1]);
				singlePatch.push_back(controlMesh[i][j + 1][k]);
				singlelSurface.push_back(singlePatch);
				singlePatch.clear();
			}
		}
		allSurface.push_back(singlelSurface);
		singlelSurface.clear();
	}
	this->data_controlMesh = allSurface;
}

void ModelDraw::drawControlLine(void) {
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	for (int i = 0; i < data_controlLine.size(); i++) {
		for (int j = 0; j < data_controlLine[i].size(); j++) {
			drawControlPoint(data_controlLine[i][j].getX(), data_controlLine[i][j].getY(), data_controlLine[i][j].getZ(), 0.1, 16, 16);
		}
	}

	glLineWidth(1.0);
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);

	for (int i = 0; i < data_controlLine.size(); i++) {
		glBegin(GL_LINE_STRIP);
		for (int j = 0; j < data_controlLine[i].size(); j++) {
			glVertex3d(data_controlLine[i][j].getX(), data_controlLine[i][j].getY(), data_controlLine[i][j].getZ());
		}
		glEnd();
	}
	//编辑模式下，使用填充样式
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

void ModelDraw::drawControlMesh(void) {
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	for (int i = 0; i < data_controlMesh.size(); i++) {
		for (int j = 0; j < data_controlMesh[i].size(); j++) {
			for (int k = 0; k < data_controlMesh[i][j].size(); k++) {
				drawControlPoint(data_controlMesh[i][j][k].getX(), data_controlMesh[i][j][k].getY(), data_controlMesh[i][j][k].getZ(), ControlPointSize, 16, 16);
			}
		}
	}
	glLineWidth(ControlLineWidth);
	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	glColor4f(0.2, 0.2, 0.7, 1.0);
	glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
	for (int i = 0; i < data_controlMesh.size(); i++) {
		glBegin(GL_QUADS);
		for (int j = 0; j < data_controlMesh[i].size(); j++) {
			for (int k = 0; k < data_controlMesh[i][j].size(); k++) {
				glVertex3d(data_controlMesh[i][j][k].getX(), data_controlMesh[i][j][k].getY(), data_controlMesh[i][j][k].getZ());
			}
		}
		glEnd();
	}
	//编辑模式下，使用填充样式
	glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
}

/**************************** 场景绘制 ******************************************/

void ModelDraw::DrawScene(void) {
    static int simulationFrame = 0;
	if (is_need_update_topopt_model) {
		updateDrawData_Topopt();
		is_need_update_topopt_model = false;
	}
    else if (is_need_update_spline_simulation_model)
    {
        cout << "Draw Simulation Frame id = " << simulationFrame << endl;
        ++simulationFrame;
        analysis_flag = true;
        SetDrawModule(ViewOpenGL::HEXMESH);
        updateDrawData();
        is_need_update_spline_simulation_model = false;
        analysis_flag = false;
    }
	DrawStyle();
	ColorStyle();
	DrawSceneModel();
}

void ModelDraw::DrawSceneModel(void) {

	loadShader();
	local_coordinater_axis->Draw3DCoord(*this);

	switch (drawmodule) {
	case EMPTY:

		break;

	case SKELMESH:// 骨架绘制
		drawSkelMesh();
		break;

	case SPARSEOPTIMIZATION:// 稀疏优化结果绘制
		drawPolyMesh_OPT();
		break;

	case TRIMESH:
		drawTriMesh();
		break;

	case SURFACEMESH:// 三角网格绘制
		drawSurfaceMesh();
		break;

	case TETMESH:// 四面体网格绘制
		drawTetMesh();
		break;

	case HEXMESH:// 六面体网格绘制
		drawHexMesh();
		break;

	case SKEL:// Lee
		drawHexMesh();
		if (bshowSKEL) drawSkel();
		if (bshowAxis) drawAxis();
		if (bshowGround) drawHorizon();
		if (bshowTriMesh) drawSurfaceMesh_skel();
		break;

	case BSPLINECURVE:// B样条曲线绘制
		if (PatchChange) {
			bsplineCurveToDiscrete();
			PatchChange = false;
		}
		if (isControlLine) {
			drawControlLine();
		}
		drawCurve();
		break;
	case BSPLINESURFACE:// B样条曲面绘制
		if (PatchChange) {
			bsplineSurfaceToDiscrete();
			PatchChange = false;
		}
		if (isControlMesh) drawControlMesh();

		drawSurface();

		break;


	case BSPLINESOLID:// B样条体绘制
		if (PatchChange) {
			bsplineSurfaceToDiscrete();
			PatchChange = false;
		}
		if (isControlMesh) bsplinesolid_drawControlMesh(); //控制点绘制
		if (isHexMesh)draw_BSplineSolid_HexMesh();//样条体绘制
		break;
	case BSPLINESOLIDS:

		draw_BSplineSolid_HexMesh();
		break;

	default:
		break;

	}


	// 绘制点, 用于剪断网格线框以保证框选准确性
	if (editmodel == EditModel::PICKITEM && mouse_move_flag) {
		shaderProgram.setUniformValue("colormode", 231);
		glDepthRange(0.01, 0.9999);
		glPointSize(meshline_width * 5.f);
		core->glDrawArrays(GL_POINTS, 0, vertices_data.size() / 6);
	}

	drawSelectedObject();
}

void ModelDraw::DrawStyle(void) {

	if (this->isEnableLighting) glEnable(GL_LIGHTING);
	else glDisable(GL_LIGHTING);

	switch (drawmode) {
	case POINTS:
		core->glPointSize(4);
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_POINT);
		shaderProgram.setUniformValue("rendermode", 1);
		break;
	case FLATLINES:
		core->glLineWidth(1);// 线宽 
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		shaderProgram.setUniformValue("rendermode", 2);
		break;
	case FILL:
		core->glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
		shaderProgram.setUniformValue("rendermode", 3);
		break;
	default:
		break;
	}
}

void ModelDraw::ColorStyle(void) {
	switch (colormode) {
	case LIGHT:
		shaderProgram.setUniformValue("colormode", 1);
		break;
	case NORMAL:
		shaderProgram.setUniformValue("colormode", 2);
		break;
	default:
		break;
	}
}

void ModelDraw::pickPoint(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	float zMax = -9999999;
	double pickRadius = ControlPointSize / 5;
	//double pickRadius = 0.01;
	CPoint tempVertex;
	bool picked = false;



	if (controlLineOrMesh == 2) {
		for (auto bp : bsplinesolidmesh) {
			controlMesh = bp.getControlPoints();
			for (int i = 0; i < controlMesh.size(); i++) {
				for (int j = 0; j < controlMesh[i].size(); j++) {
					for (int k = 0; k < controlMesh[i][j].size(); k++) {
						QVector3D point0(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ());
						double dist = point2LineDis(point1, point2, point0);
						if (dist < pickRadius) {
							QVector4D fDepth(point0, 0);
							fDepth = modelmatrix * fDepth;
							if (fDepth[2] > zMax) {
								zMax = fDepth[2];
								tempVertex = controlMesh[i][j][k];
								picked = true;
							}
						}
					}
				}
			}
			Is_Highlight[tempVertex] = 1;//changede
		}
	}
	else if (controlLineOrMesh == 1) {
		for (int i = 0; i < controlLine.size(); i++) {
			for (int j = 0; j < controlLine[i].size(); j++) {
				QVector3D point0(controlLine[i][j].getX(), controlLine[i][j].getY(), controlLine[i][j].getZ());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < pickRadius) {
					QVector4D fDepth(point0, 0);
					fDepth = modelmatrix * fDepth;
					if (fDepth[2] > zMax) {
						zMax = fDepth[2];
						tempVertex = controlLine[i][j];
						picked = true;
						curvepointIndex[0] = i;
						curvepointIndex[1] = j;
					}
				}
			}
		}
	}


	if (picked) {
		pickedItemPosition[0] = tempVertex.getX();
		pickedItemPosition[1] = tempVertex.getY();
		pickedItemPosition[2] = tempVertex.getZ();
	}
	else std::cout << "no picked !!" << std::endl;
}
void ModelDraw::pickPoint_chooseface(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	float zMax = -9999999;
	double pickRadius = ControlPointSize / 5;
	//double pickRadius = 0.01;
	CPoint tempVertex;
	bool picked = false;




	for (auto bp : bsplinesolidmesh) {
		vector<CPoint> V = bp.U_0_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 1;//new
					picked = true;
				}
			}
		}
		V = bp.U_1_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 2;//new
					picked = true;
				}
			}
		}
		V = bp.V_0_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 3;//new
					picked = true;
				}
			}
		}
		V = bp.V_1_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 4;//new
					picked = true;
				}
			}
		}
		V = bp.W_0_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 5;//new
					picked = true;
				}
			}
		}
		V = bp.W_1_discrete_point;
		for (int i = 0; i < V.size(); i++)
		{
			QVector3D point0(V[i].getX(), V[i].getY(), V[i].getZ());
			double dist = point2LineDis(point1, point2, point0);
			if (dist < pickRadius) {
				QVector4D fDepth(point0, 0);
				fDepth = modelmatrix * fDepth;
				if (fDepth[2] > zMax) {
					zMax = fDepth[2];
					bp_index = bp.gethandle();//new
					face_handle = 6;//new
					picked = true;
				}
			}
		}
	}

	if (picked) {
		chooseface_index = bp_index;
		std::cout << "bp_index=" << bp_index << std::endl;
		std::cout << "face_handle=" << face_handle << std::endl;
	}
	else std::cout << "no picked !!" << std::endl;
}

void  ModelDraw::ClearAllVector() {
	Is_Highlight.clear();
	AllControlPoint.clear();
	PointBelongToSolid.clear();
	bsplinesolidmesh.clear();
	singlepatch_vertics_number.clear();
	vertices_data.clear();


}

double ModelDraw::point2LineDis(QVector3D point1, QVector3D point2, QVector3D point0) {

	//直线的法向量(p,q,r)
	QVector3D normal = (point1 - point2).normalized();

	//两个向量
	QVector3D normal01 = (point1 - point0).normalized();
	QVector3D normal02 = (point2 - point0).normalized();

	//进行叉乘
	double cross_x = normal01.y() * normal02.z() - normal02.y() * normal01.z();
	double cross_y = normal02.x() * normal01.z() - normal02.z() * normal01.x();
	double cross_z = normal01.x() * normal02.y() - normal01.y() * normal02.x();

	double cross_length = sqrt(cross_x * cross_x + cross_y * cross_y + cross_z * cross_z);
	double dx = (point2 - point1).x();
	double dy = (point2 - point1).y();
	double dz = (point2 - point1).z();

	double qb_length = sqrt(dx * dx + dy * dy + dz * dz);

	double ds = cross_length / qb_length;
	return ds;

}

bool ModelDraw::isIntersectTriangle(QVector3D orig, QVector3D end, QVector3D v0, QVector3D v1, QVector3D v2, QVector3D& intersection) {
	// Refer to https://www.cnblogs.com/graphics/archive/2010/08/09/1795348.html
	float t, u, v;
	QVector3D dir = end - orig;
	QVector3D E1 = v1 - v0;
	QVector3D E2 = v2 - v0;
	// 叉乘
	QVector3D P(dir.y() * E2.z() - E2.y() * dir.z(),
		E2.x() * dir.z() - E2.z() * dir.x(),
		dir.x() * E2.y() - dir.y() * E2.x());
	// 点乘
	float det = E1.x() * P.x() + E1.y() * P.y() + E1.z() * P.z();
	QVector3D T;
	if (det > 0) {
		T = orig - v0;
	}
	else {
		T = v0 - orig;
		det = -det;
	}
	if (det < 0)
		return false;

	// Calculate u and make sure u <= 1
	u = T.x() * P.x() + T.y() * P.y() + T.z() * P.z();
	if (u < 0.0f || u > det)
		return false;

	// Q
	QVector3D Q(T.y() * E1.z() - E1.y() * T.z(),
		E1.x() * T.z() - E1.z() * T.x(),
		T.x() * E1.y() - T.y() * E1.x());

	// Calculate v and make sure u + v <= 1
	v = dir.x() * Q.x() + dir.y() * Q.y() + dir.z() * Q.z();
	if (v < 0.0f || u + v > det)
		return false;

	//// Calculate t, scale parameters, ray intersects triangle
	//t = E2.x() * Q.x() + E2.y() * Q.y() + E2.z() * Q.z();
	float fInvDet = 1.0f / det;
	//t *= fInvDet;
	u *= fInvDet;
	v *= fInvDet;
	intersection = u * E1 + v * E2 + v0;
	return true;
}

//void ModelDraw::clickCurveControlPoint(QPoint point) {
//	vector<double> p(3);
//	pickPoint_getDepth(point.x(), point.y(), p, depth);//获取鼠标处的深度值
//	if (depth >= 0.999999999) {
//		depth = 0.9;
//	}
//	CPoint clickpoint(p[0], p[1], p[2]);//选择的点坐标
//	double distance = 999;
//	for (int i = 0; i < controlLine.size(); i++) {
//		for (int j = 0; j < controlLine[i].size(); j++) {
//			//找到离点击点最近的控制点
//			double tmpdistance = sqrt(pow(controlLine[i][j].getX() - clickpoint.getX(), 2) +
//				pow(controlLine[i][j].getY() - clickpoint.getY(), 2) +
//				pow(controlLine[i][j].getZ() - clickpoint.getZ(), 2));
//			if (tmpdistance < distance) {//保存离点击点最近的控制点的索引
//				curvepointIndex[0] = i;
//				curvepointIndex[1] = j;
//				distance = tmpdistance;
//			}
//		}
//	}
//}

//void ModelDraw::clickSurfaceControlPoint(QPoint point) {
//	vector<double> p(3);
//	pickPoint_getDepth(point.x(), point.y(), p, depth);//获取鼠标处的深度值
//	if (depth >= 0.999999999) {
//		depth = 0.996766;
//	}
//	CPoint clickpoint(p[0], p[1], p[2]);//选择的点坐标
//	double distance = 999;
//	for (int i = 0; i < controlMesh.size(); i++) {
//		for (int j = 0; j < controlMesh[i].size(); j++) {
//			for (int k = 0; k < controlMesh[i][j].size(); k++) {
//				//找到离点击点最近的控制点
//				double tmpdistance = sqrt(pow(controlMesh[i][j][k].getX() - clickpoint.getX(), 2) +
//					pow(controlMesh[i][j][k].getY() - clickpoint.getY(), 2) +
//					pow(controlMesh[i][j][k].getZ() - clickpoint.getZ(), 2));
//				if (tmpdistance < distance) {//保存离点击点最近的控制点的索引
//					surfacepointIndex[0] = i;
//					surfacepointIndex[1] = j;
//					surfacepointIndex[2] = k;
//					distance = tmpdistance;
//				}
//			}
//		}
//	}
//}

//void ModelDraw::moveCurveControlPoint(QPoint point) {
//	//将新的位置点point设置为选中控制点的新坐标
//	vector<double> p(3);
//	pickPointDepth(point.x(), point.y(), p, depth);
//	controlLine[curvepointIndex[0]][curvepointIndex[1]].setX(p[0]);
//	controlLine[curvepointIndex[0]][curvepointIndex[1]].setY(p[1]);
//	controlLine[curvepointIndex[0]][curvepointIndex[1]].setZ(p[2]);
//
//	//重新绘制控制点和曲线
//	int uindex = curvepointIndex[1] + 1;//控制点索引从1开始
//	curvesDraw[curvepointIndex[0]].setSingleControlPoint(uindex, controlLine[curvepointIndex[0]][curvepointIndex[1]]);//更新单个控制点
//	/*curvesDraw.clear();
//	vector<CBSplineCurve>().swap(curvesDraw);
//	curvesDraw = bsplineCurves;*/
//	bsplineCurveToDiscrete();
//	update();
//}

//void ModelDraw::moveSurfaceControlPoint(QPoint point) {
//	//将新的位置点point设置为选中控制点的新坐标
//	vector<double> p(3);
//	pickPointDepth(point.x(), point.y(), p, depth);
//	controlMesh[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setX(p[0]);
//	controlMesh[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setY(p[1]);
//	controlMesh[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setZ(p[2]);
//	//重新绘制控制点和曲面
//	int uindex = surfacepointIndex[1] + 1;//控制点索引从1开始
//	int vindex = surfacepointIndex[2] + 1;
//	surfacesDraw[surfacepointIndex[0]].setSingleControlPoint(uindex, vindex, controlMesh[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]]);//将更改的曲面的控制点更新
//	/*surfacesDraw.clear();
//	vector<CBSplineSurface>().swap(surfacesDraw);
//	surfacesDraw=bsplineSurfaces;*/
//	bsplineSurfaceToDiscrete();
//	controlpointsToDiscrete();
//	update();
//}

void ModelDraw::removeSelectedOne() {
	if (!selected_v_indices.empty() && pickType == PickType::PickPoint) selected_v_indices.resize(selected_v_indices.size() - 1);
	if (!selected_e_indices.empty() && pickType == PickType::PickEdge) selected_e_indices.resize(selected_e_indices.size() - 1);
	if (!selected_f_indices.empty() && pickType == PickType::PickFace) selected_f_indices.resize(selected_f_indices.size() - 1);
	if (!selected_c_indices.empty() && pickType == PickType::PickCell) selected_c_indices.resize(selected_c_indices.size() - 1);
}

void ModelDraw::removeSelectedAll() {
	pickedge.clear();
	pickface.clear();
	this->selected_v_indices.clear();
	this->selected_e_indices.clear();
	this->selected_f_indices.clear();
	this->selected_c_indices.clear();
	vector<int>().swap(this->selected_v_indices);
	vector<int>().swap(this->selected_e_indices);
	vector<int>().swap(this->selected_f_indices);
	vector<int>().swap(this->selected_c_indices);
}

void ModelDraw::deleteSelectedAll() {
	switch (pickType) {
	case PickPoint: {
		for (auto& i : selected_v_indices) {
			MeshKernel::iGameVertexHandle vh(i);
			if (hexMesh) hexMesh->DeleteVertex(vh);
			else if (surface_mesh) surface_mesh->DeleteVertex(vh);
		}
		selected_v_indices.clear();
		break;
	}
	case PickEdge: {
		for (auto& i : selected_e_indices) {
			MeshKernel::iGameEdgeHandle eh(i);
			if (hexMesh) hexMesh->DeleteEdge(eh);
			else if (surface_mesh) surface_mesh->DeleteEdge(eh);
		}
		selected_e_indices.clear();
		break;
	}
	case PickFace: {
		for (auto& i : selected_f_indices) {
			MeshKernel::iGameFaceHandle fh(i);
			if (hexMesh) hexMesh->DeleteFace(fh);
			else if (surface_mesh) surface_mesh->DeleteFace(fh);
		}
		selected_f_indices.clear();
		break;
	}
	case PickCell: {
		for (auto& i : selected_c_indices) {
			MeshKernel::iGameCellHandle ch(i);
			hexMesh->DeleteCell(ch);
		}
		selected_c_indices.clear();
		break;
	}
	}
	printf("delete selected all success\n");
	updateDrawData();
	update();
}

void ModelDraw::add_new_cell() {
	if (!hexMesh) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Support Hex Mesh Only! "));
		return;
	}
	if (selected_f_indices.size() == 2) {
		FH fh1(selected_f_indices[0]), fh2(selected_f_indices[1]);
		add_cell_2fh(*hexMesh, fh1, fh2);
	}
	else if (selected_e_indices.size() == 1 && selected_f_indices.size() == 1) {
		EH eh(selected_e_indices[0]);
		FH fh(selected_f_indices[0]);
		add_cell_eh_fh(*hexMesh, eh, fh);
	}
	else if (selected_f_indices.size() == 3) {
		FH fh1(selected_f_indices[0]), fh2(selected_f_indices[1]), fh3(selected_f_indices[2]);
		add_cell_3fh(*hexMesh, fh1, fh2, fh3);
	}
	else if (selected_f_indices.size() == 1) {
		FH fh(selected_f_indices[0]);
		add_cell_fh(*hexMesh, fh);
	}

	updateDrawData();
}

void ModelDraw::initializePick() {
	selected_c_indices.clear();
	solutionSkel->selectNodes.clear();
	solutionSkel->selectNode = nullptr;
}

void ModelDraw::pickEvent(double winX, double winY) {
	switch (drawmodule) {
	case ViewOpenGL::EMPTY:
		break;

	case ViewOpenGL::SURFACEMESH:
#ifdef DEBUG
		printf("pick triMesh, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
		pickEvent_TriMesh_BVH(winX, winY);
		showPickInfo();
		break;

	case ViewOpenGL::TETMESH:
#ifdef DEBUG
		printf("pick tetMesh, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
		pickEvent_TetMesh_BVH(winX, winY);
		showPickInfo();
		break;

	case ViewOpenGL::HEXMESH:
#ifdef DEBUG
		printf("pick hexMesh, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
		pickEvent_HexMesh_BVH(winX, winY);
		//showPickInfo();
		break;

	case ViewOpenGL::SKEL:
#ifdef DEBUG
		printf("pick SKEL, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
		pickEvent_SKEL(winX, winY);
		showPickInfo();
		break;

    case ViewOpenGL::SKELMESH:
#ifdef DEBUG
        printf("pick SKELMESH, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
        pickEvent_Skeleton(winX, winY);
        break;
	case ViewOpenGL::HYBRIDMESH:
		break;
	case ViewOpenGL::BSPLINECURVE:
		break;
	case ViewOpenGL::BSPLINESURFACE:
		break;
#ifdef DEBUG
		printf("pick BSPLINESURFACE, mouse at (%.2f, %.2f)\n", winX, winY);
#endif // DEBUG
		pickEvent_Bsplinesurface(winX, winY);
		break;
	default:
		break;
	}
}

void ModelDraw::pickEvent_Skeleton(double winX, double winY) {
    if (!skeletal_mesh) return;

    QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
    QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
    QVector4D pointEnd = point + 2 * toward;
    QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
    QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
    QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
    QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());


    Vector3d orig(point1.x(), point1.y(), point1.z());
    auto _dir = point2 - point1;
    Vector3d dir(_dir.x(), _dir.y(), _dir.z());
    Ray ray(orig, -dir);
    auto isec = bvh_tree->getIntersection(ray);
    MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
    QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);


    if (tempFaceHandle != -1) {
        double min_distance = numeric_limits<double>::max();
        auto& tempFace = skeletal_mesh->faces(tempFaceHandle);
        int nearest_handle(-1);
        MeshKernel::iGameVertex pos(-1, -1, -1);
        switch (pickType) {

            case PickType::PickPoint:
                for (auto& vh : tempFace.getVertexHandle()) {
                    auto& v = skeletal_mesh->vertices(vh);
                    QVector3D point0(v.x(), v.y(), v.z());
                    double dist = point2LineDis(point1, point2, point0);
                    if (dist < min_distance) {
                        min_distance = dist;
                        nearest_handle = vh;
                    }
                }
                if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
                    selected_v_indices.push_back(nearest_handle);
                }
                break;

            case PickType::PickEdge:
                for (auto& eh : tempFace.getEdgeHandle()) {
                    auto& e = skeletal_mesh->edges(eh);
                    auto& v1 = skeletal_mesh->vertices(e.vh1());
                    auto& v2 = skeletal_mesh->vertices(e.vh2());
                    QVector3D ev1(v1.x(), v1.y(), v1.z());
                    QVector3D ev2(v2.x(), v2.y(), v2.z());
                    double dist = point2LineDis(ev1, ev2, tempIntersection);
                    if (dist < min_distance) {
                        min_distance = dist;
                        nearest_handle = eh;
                    }
                }
                if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
                    selected_e_indices.push_back(nearest_handle);
                }
                break;

            case PickType::PickFace:
                nearest_handle = tempFaceHandle;
                if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
                    selected_f_indices.push_back(nearest_handle);
                }
                break;

            default:
                break;
        }
    }
}

void ModelDraw::pickEvent_Bsplinesurface(double winX, double winY){

}

void ModelDraw::pickEvent_TriMesh(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	double zMax = -numeric_limits<double>::max();

	std::vector<QVector3D> tempPoints;
	MeshKernel::iGameFaceHandle tempFaceHandle(-1);
	MeshKernel::iGameFace tempFace;
	QVector3D tempIntersection;
	for (auto& fp : surface_mesh->allfaces()) {
		MeshKernel::iGameFaceHandle fh = fp.first;
		auto& f = fp.second;
		auto vhs = f.getVertexHandle();
		vector<QVector3D> face_vertices;
		QVector4D face_depth(0, 0, 0, 0);
		for (auto& vh : vhs) {
			auto& v = surface_mesh->vertices(vh);
			auto qv = QVector3D(v.x(), v.y(), v.z());
			face_vertices.push_back(qv);
			face_depth += QVector4D(qv, 1);
		}
		face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
		face_depth = modelmatrix * face_depth;
		bool intersection_flag = false;
		QVector3D intersection_point;
		if (face_vertices.size() == 3) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
		}
		else if (face_vertices.size() == 4) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
				| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
		}
		if (intersection_flag && zMax < face_depth[2]) {
			zMax = face_depth[2];
			tempFaceHandle = fh;
			tempFace = f;
			tempIntersection = intersection_point;
			tempPoints = face_vertices;
		}
	}

	/*Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	auto isec = bvh_tree->getIntersection(ray);

	std::cout << "\nOld-Method fh = " << tempFaceHandle << "; " <<
		" BVH fh = " << isec.index << "\n";
	printf("Old-IsecPos = (%.4f, %.4f, %.4f), BVH IsecPos = (%.4f, %.4f, %.4f)\n", tempIntersection.x(), tempIntersection.y(), tempIntersection.z(),
		isec.pos[0], isec.pos[1], isec.pos[2]);*/


	if (tempFaceHandle != -1) {
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		MeshKernel::iGameVertex pos(-1, -1, -1);
		switch (pickType) {

		case PickType::PickPoint:
			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = surface_mesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			pos = surface_mesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			//printf("the selected vertex handle is %d: (%.2f, %.2f, %.2f)\n", nearest_handle, pos.x(), pos.y(), pos.z());
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
#endif // DEBUG
			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = surface_mesh->edges(eh);
				auto& v1 = surface_mesh->vertices(e.vh1());
				auto& v2 = surface_mesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("the selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
			break;
		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("the selected face handle is %d\n", nearest_handle);
#endif // DEBUG


			break;
		default:
			break;
		}
	}
	else {
#ifdef DEBUG
		printf("no intersection!!!\n");
#endif // DEBUG
	}


}

void ModelDraw::pickEvent_TriMesh_BVH(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());


	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	auto isec = bvh_tree->getIntersection(ray);
	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);


	if (tempFaceHandle != -1) {
		double min_distance = numeric_limits<double>::max();
		auto& tempFace = surface_mesh->faces(tempFaceHandle);
		int nearest_handle(-1);
		MeshKernel::iGameVertex pos(-1, -1, -1);
		switch (pickType) {

		case PickType::PickPoint:
			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = surface_mesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			pos = surface_mesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			//printf("the selected vertex handle is %d: (%.2f, %.2f, %.2f)\n", nearest_handle, pos.x(), pos.y(), pos.z());
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
#endif // DEBUG
			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = surface_mesh->edges(eh);
				auto& v1 = surface_mesh->vertices(e.vh1());
				auto& v2 = surface_mesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("the selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
			break;
		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("the selected face handle is %d\n", nearest_handle);
#endif // DEBUG


			break;
		default:
			break;
		}
	}
	else {
#ifdef DEBUG
		printf("no intersection!!!\n");
#endif // DEBUG
	}


}

void ModelDraw::pickEvent_TetMesh(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);	// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	double zMax = -numeric_limits<double>::max();

	std::vector<QVector3D> tempPoints;
	MeshKernel::iGameFaceHandle tempFaceHandle(-1);
	MeshKernel::iGameFace tempFace;
	QVector3D tempIntersection;
	for (auto& fp : tetMesh->allfaces()) {// this->hexMesh->boundary_face
		if (!tetMesh->isOnBoundary(fp.first)) continue;
		MeshKernel::iGameFaceHandle fh = fp.first;
		auto& f = fp.second;
		auto vhs = f.getVertexHandle();
		vector<QVector3D> face_vertices;
		QVector4D face_depth(0, 0, 0, 0);
		for (auto& vh : vhs) {
			auto& v = tetMesh->vertices(vh);
			auto qv = QVector3D(v.x(), v.y(), v.z());
			face_vertices.push_back(qv);
			face_depth += QVector4D(qv, 1);
		}
		face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
		face_depth = modelmatrix * face_depth;
		bool intersection_flag = false;
		QVector3D intersection_point;
		if (face_vertices.size() == 3) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
		}
		else if (face_vertices.size() == 4) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
				| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
		}
		if (intersection_flag && zMax < face_depth[2]) {
			zMax = face_depth[2];
			tempFaceHandle = fh;
			tempFace = f;
			tempIntersection = intersection_point;
			tempPoints = face_vertices;
		}
	}
	if (tempFaceHandle != -1) {

		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = tetMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) {
		case PickType::PickPoint:

			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = tetMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
			//printf("the selected vertex handle is %d\n", nearest_handle);
			pos = tetMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = tetMesh->edges(eh);
				auto& v1 = tetMesh->vertices(e.vh1());
				auto& v2 = tetMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
			printf("the selected edge handle is %d\n", nearest_handle);
			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
			printf("the selected face handle is %d\n", nearest_handle);
			break;

		case PickType::PickCell:
			for (auto& cp : tetMesh->allcells()) {
				if (tetMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
			printf("the selected cell handle is %d\n", nearest_handle);
			break;

		default:
			break;
		}
	}
	else {
		printf("no intersection!!!\n");
	}
}

void ModelDraw::pickEvent_TetMesh_BVH(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);	// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	auto isec = bvh_tree->getIntersection(ray);
	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);

	if (tempFaceHandle != -1) {

		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = tetMesh->CellSize();
		MeshKernel::iGameVertex pos;
		auto& tempFace = tetMesh->faces(tempFaceHandle);

		switch (pickType) {
		case PickType::PickPoint:

			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = tetMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
			//printf("the selected vertex handle is %d\n", nearest_handle);
			pos = tetMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = tetMesh->edges(eh);
				auto& v1 = tetMesh->vertices(e.vh1());
				auto& v2 = tetMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
			printf("the selected edge handle is %d\n", nearest_handle);
			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
			printf("the selected face handle is %d\n", nearest_handle);
			break;

		case PickType::PickCell:
			for (auto& cp : tetMesh->allcells()) {
				if (tetMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
			printf("the selected cell handle is %d\n", nearest_handle);
			break;

		default:
			break;
		}
	}
	else {
		printf("no intersection!!!\n");
	}
}

void ModelDraw::pick_point_with_depth(float x, float y, double(&p)[3], double itemPosition[3]) {
	//if (!this->surface_mesh && !this->tetMesh) return;
	QVector4D point(x, y, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(x, y, -1, 1.0);	// 指向toward与投影矩阵的fov相关
	QVector4D pointView(0, 0, 1, 1.0);
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector4D tpointView = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointView;
	// point1为射线起点，point2为射线终点
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	QVector3D faceNormal(tpointView.x() / tpointView.w(), tpointView.y() / tpointView.w(), tpointView.z() / tpointView.w());

	// 射线向量
	double e[3]{ point2.x() - point1.x(),point2.y() - point1.y(), point2.z() - point1.z() };
	// 射线起点
	double PL[3]{ point1.x(),point1.y(), point1.z() };
	// 平面向量
	double q[3]{ itemPosition[0] - faceNormal.x(),itemPosition[1] - faceNormal.y(), itemPosition[2] - faceNormal.z() };

	double t = (q[0] * (itemPosition[0] - PL[0]) + q[1] * (itemPosition[1] - PL[1]) + q[2] * (itemPosition[2] - PL[2])) / (q[0] * e[0] + q[1] * e[1] + q[2] * e[2]);
	p[0] = PL[0] + e[0] * t;
	p[1] = PL[1] + e[1] * t;
	p[2] = PL[2] + e[2] * t;
}


bool ModelDraw::pickEvent_SKEL(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0);		// point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0);		// 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());
	double zMax = -numeric_limits<double>::max();
	double item[3]{ 0, 0, 0 };

	std::vector<QVector3D> selectPoints;
	MeshKernel::iGameCell selectCell;
	MeshKernel::iGameCellHandle selectCellHandle(-1);
	MeshKernel::iGameFace selectFace;
	MeshKernel::iGameFaceHandle selectFaceHandle(-1);
	QVector3D selectIntersection;
	for (auto& cp : hexMesh->allcells()) {		// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		MeshKernel::iGameCell c = cp.second;
		for (auto& fh : c.getFaceHandle()) {
			auto& f = hexMesh->faces(fh);
			auto vhs = f.getVertexHandle();
			vector<QVector3D> face_vertices;
			QVector4D face_depth(0, 0, 0, 0);
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				auto qv = QVector3D(v.x(), v.y(), v.z());
				face_vertices.push_back(qv);
				face_depth += QVector4D(qv, 1);
			}
			face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
			face_depth = modelmatrix * face_depth;
			bool intersection_flag = false;
			QVector3D intersection_point;
			if (face_vertices.size() == 3) {
				intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
			}
			else if (face_vertices.size() == 4) {
				intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
					| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
			}
			if (intersection_flag && zMax < face_depth[2]) {
				zMax = face_depth[2];
				selectCellHandle = cp.first;
				selectCell = cp.second;
				selectFaceHandle = fh;
				selectFace = f;
				selectPoints = face_vertices;
				selectIntersection = intersection_point;
			}
		}
	}

	// 如果有选中的体
	if (selectCellHandle != -1) {
		MeshKernel::iGameVertex pos;
		switch (pickType) {
		case PickType::PickPoint:
			/*
			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = hexMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
			//printf("the selected vertex handle is %d\n", nearest_handle);
			pos = hexMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
			break;
			*/
		case PickType::PickEdge:
			/*
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = hexMesh->edges(eh);
				auto& v1 = hexMesh->vertices(e.vh1());
				auto& v2 = hexMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
			}
			printf("the selected edge handle is %d\n", nearest_handle);
			break;
			*/
		case PickType::PickFace:
			/*
			nearest_handle = (int)tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
			printf("the selected face handle is %d\n", nearest_handle);
			break;
			*/

		case PickType::PickCell:
			px = selectIntersection[0];
			py = selectIntersection[1];
			pz = selectIntersection[2];
			if ((selected_c_indices.size() == 1 && selected_c_indices[0] != selectCellHandle)
				|| selected_c_indices.size() == 0) {
				selected_c_indices.clear();
				selected_c_indices.push_back(selectCellHandle);
				pos = hexMesh->getCellCenter(selectCellHandle);
				// check with all node position
				solutionSkel->selectNodes.clear();
				BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(pos[0], pos[1], pos[2]));
				SingalFather::getInstance()->getData(rlt->getPosition().x, rlt->getPosition().y, rlt->getPosition().z,
					rlt->root->getRotateX(), rlt->root->getRotateY(), rlt->root->getRotateZ(),
					rlt->root->getLength(), rlt->root->getWidth(), rlt->root->getHeight());

				// 坐标系设置
				local_coordinater_axis->set_orign(pos.x(), pos.y(), pos.z());
				local_coordinater_axis->set_mode_type(select_coordinater_mode);
				printf("the selected cell mid is %lf %lf %lf\n", pos.x(), pos.y(), pos.z());
			}
			printf("the selected cell handle is %d\n", selectCellHandle);

			break;

			/*
			case PickType::PickCells:
				for (auto& cp : hexMesh->allcells()) {
					if (hexMesh->isOnBoundary(cp.first)) {
						auto cell_fhs = cp.second.getFaceHandle();
						if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
							nearest_handle = cp.first;
							break;
						}
					}
				}
				px = tempIntersection[0];
				py = tempIntersection[1];
				pz = tempIntersection[2];

				if ((it = std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle)) == selected_c_indices.end()) {
					selected_c_indices.push_back(nearest_handle);
				}
				else {
					selected_c_indices.erase(it);
				}
				// update node
				solutionSkel->selectNodes.clear();
				for (int index = 0; index < selected_c_indices.size(); index++) {
					iGameVertex tmp(0.0, 0.0, 0.0);
					auto& cell = hexMesh->cells((iGameCellHandle)index);
					auto vhs = cell.getVertexHandle();
					for (auto& vh : vhs) {
						tmp += hexMesh->vertices(vh);
					}
					tmp = tmp / 8.0;
					// check with all node position
					BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(tmp[0], tmp[1], tmp[2]));
					if (rlt) {
						solutionSkel->selectNodes.push_back(rlt);
					}
				}

				printf("the selected cell handle is %d\n", nearest_handle);
				break;

			*/
		default:
			break;
		}
		return true;
	}
	else {
		printf("no intersection!!!\n");
		return false;
	}
}


void ModelDraw::rupdateGL() {
	if (selected_c_indices.size() == 1) {
		auto& cell = hexMesh->cells((iGameCellHandle)selected_c_indices[0]);
		auto vhs = cell.getVertexHandle();
		for (int i = 0; i < 8; i++) {
			iGameVertex t(solutionSkel->selectNode->root->p[i].x, solutionSkel->selectNode->root->p[i].y,
				solutionSkel->selectNode->root->p[i].z);
			auto& v = hexMesh->vertices(vhs[i]);
			v = t;
		}
	}
	SetDrawModule(SKEL);
	updateDrawData();
	//updateDrawData_HexMesh();
	//update();
	//iGameDFS();
	return;
}

void ModelDraw::drawNormalSkel() {
	ModelClear();
	hexMesh = new MeshKernel::VolumeMesh();
	solutionSkel->genViewHex(*hexMesh);
	int exterior_vcnt = 0;
	for (auto& vp : hexMesh->allvertices()) {
		if (hexMesh->isOnBoundary(vp.first)) exterior_vcnt++;
	}
	setRoot(solutionSkel->getRoot());
	setNodeVec(solutionSkel->getNodeVec());
	printf("HEXMESH: vertices = %d, cells = %d, body_vcnt = %d\n", hexMesh->VertexSize(), hexMesh->allcells().size(), exterior_vcnt);
	SetDrawModule(SKEL);
	//updateMeshModel();
	updateDrawData();
	update();
	return;
}


void ModelDraw::pickEvent_HexMesh(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	double zMax = -numeric_limits<double>::max();
	std::vector<QVector3D> tempPoints;
	MeshKernel::iGameFaceHandle tempFaceHandle(-1);
	MeshKernel::iGameFace tempFace;
	QVector3D tempIntersection;

	for (auto& fp : hexMesh->allfaces()) {// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(fp.first)) continue;
		MeshKernel::iGameFaceHandle fh = fp.first;
		auto& f = fp.second;
		auto vhs = f.getVertexHandle();
		vector<QVector3D> face_vertices;
		QVector4D face_depth(0, 0, 0, 0);
		for (auto& vh : vhs) {
			auto& v = hexMesh->vertices(vh);
			auto qv = QVector3D(v.x(), v.y(), v.z());
			face_vertices.push_back(qv);
			face_depth += QVector4D(qv, 1);
		}
		face_depth = QVector4D(face_depth.x() / face_depth.w(), face_depth.y() / face_depth.w(), face_depth.z() / face_depth.w(), 1);
		face_depth = modelmatrix * face_depth;
		bool intersection_flag = false;
		QVector3D intersection_point;
		if (face_vertices.size() == 3) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point);
		}
		else if (face_vertices.size() == 4) {
			intersection_flag = isIntersectTriangle(point1, point2, face_vertices[0], face_vertices[1], face_vertices[2], intersection_point)
				| isIntersectTriangle(point1, point2, face_vertices[2], face_vertices[3], face_vertices[0], intersection_point);
		}
		if (intersection_flag && zMax < face_depth[2]) {
			zMax = face_depth[2];
			tempFaceHandle = fh;
			tempFace = f;
			tempIntersection = intersection_point;
			tempPoints = face_vertices;
		}
	}

	/*Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	int time_bvh = clock() - time2;

	std::cout << "\nOld-Method fh = " << tempFaceHandle << ", cost: " << time_old << "ms. " <<
		" BVH fh = " << isec.index << ", cost: " << time_bvh << "ms.\n";
	printf("Old-IsecPos = (%.4f, %.4f, %.4f), BVH IsecPos = (%.4f, %.4f, %.4f)\n", tempIntersection.x(), tempIntersection.y(), tempIntersection.z(),
		isec.pos[0], isec.pos[1], isec.pos[2]);*/

	if (tempFaceHandle != -1) {

		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = hexMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) {
		case PickType::PickPoint:

			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = hexMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected vertex handle is %d\n", nearest_handle);
			pos = hexMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
#endif // DEBUG


			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = hexMesh->edges(eh);
				auto& v1 = hexMesh->vertices(e.vh1());
				auto& v2 = hexMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
				// added for Top Opt
				/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
				pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
			}
#ifdef DEBUG
			printf("The selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected face handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickCell:
			for (auto& cp : hexMesh->allcells()) {
				if (hexMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected cell handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		default:
			break;
		}
	}
	else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}

void ModelDraw::pickEvent_HexMesh_BVH(double winX, double winY) {
	QVector4D point(winX, winY, 0.9, 1.0); // point为鼠标点击的位置
	QVector4D toward(winX, winY, -1, 1.0); // 指向toward与投影矩阵的fov相关
	QVector4D pointEnd = point + 2 * toward;
	QVector4D tpoint = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * point;
	QVector4D tpointEnd = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pointEnd;
	QVector3D point1(tpoint.x() / tpoint.w(), tpoint.y() / tpoint.w(), tpoint.z() / tpoint.w());
	QVector3D point2(tpointEnd.x() / tpointEnd.w(), tpointEnd.y() / tpointEnd.w(), tpointEnd.z() / tpointEnd.w());

	Vector3d orig(point1.x(), point1.y(), point1.z());
	auto _dir = point2 - point1;
	Vector3d dir(_dir.x(), _dir.y(), _dir.z());
	Ray ray(orig, -dir);
	time_t time2 = clock();
	auto isec = bvh_tree->getIntersection(ray);
	//std::cout << "There exists " << bvh_tree->intersection_count << " intersections.\n";
	int time_bvh = clock() - time2;

	MeshKernel::iGameFaceHandle tempFaceHandle(isec.index);
	QVector3D tempIntersection(isec.pos[0], isec.pos[1], isec.pos[2]);

	if (tempFaceHandle != -1) {
		auto& tempFace = hexMesh->faces(tempFaceHandle);
		double min_distance = numeric_limits<double>::max();
		int nearest_handle(-1);
		auto ccnt = hexMesh->CellSize();
		MeshKernel::iGameVertex pos;

		switch (pickType) {
		case PickType::PickPoint:

			for (auto& vh : tempFace.getVertexHandle()) {
				auto& v = hexMesh->vertices(vh);
				QVector3D point0(v.x(), v.y(), v.z());
				double dist = point2LineDis(point1, point2, point0);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = vh;
				}
			}
			if (std::find(selected_v_indices.begin(), selected_v_indices.end(), nearest_handle) == selected_v_indices.end()) {
				selected_v_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected vertex handle is %d\n", nearest_handle);
			pos = hexMesh->vertices(MeshKernel::iGameVertexHandle(nearest_handle));
			this->pick_info = "VH " + std::to_string(nearest_handle) + ": (" + std::to_string(pos.x()) + ", " + std::to_string(pos.y()) + ", " +
				std::to_string(pos.z()) + ")";
			std::cout << pick_info << std::endl;
#endif // DEBUG


			break;

		case PickType::PickEdge:
			for (auto& eh : tempFace.getEdgeHandle()) {
				auto& e = hexMesh->edges(eh);
				auto& v1 = hexMesh->vertices(e.vh1());
				auto& v2 = hexMesh->vertices(e.vh2());
				QVector3D ev1(v1.x(), v1.y(), v1.z());
				QVector3D ev2(v2.x(), v2.y(), v2.z());
				double dist = point2LineDis(ev1, ev2, tempIntersection);
				if (dist < min_distance) {
					min_distance = dist;
					nearest_handle = eh;
				}
			}
			if (std::find(selected_e_indices.begin(), selected_e_indices.end(), nearest_handle) == selected_e_indices.end()) {
				selected_e_indices.push_back(nearest_handle);
				// added for Top Opt
				/*auto e = hexMesh->edges(MeshKernel::iGameEdgeHandle(nearest_handle));
				pickedge.push_back({(int)e.vh1(), (int)e.vh2()});*/
			}
#ifdef DEBUG
			printf("The selected edge handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickFace:
			nearest_handle = tempFaceHandle;
			if (std::find(selected_f_indices.begin(), selected_f_indices.end(), nearest_handle) == selected_f_indices.end()) {
				selected_f_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected face handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		case PickType::PickCell:
			for (auto& cp : hexMesh->allcells()) {
				if (hexMesh->isOnBoundary(cp.first)) {
					auto cell_fhs = cp.second.getFaceHandle();
					if (std::find(cell_fhs.begin(), cell_fhs.end(), tempFaceHandle) != cell_fhs.end()) {
						nearest_handle = cp.first;
						break;
					}
				}
			}
			if (std::find(selected_c_indices.begin(), selected_c_indices.end(), nearest_handle) == selected_c_indices.end()) {
				selected_c_indices.push_back(nearest_handle);
			}
#ifdef DEBUG
			printf("The selected cell handle is %d\n", nearest_handle);
#endif // DEBUG
			break;

		default:
			break;
		}
	}
	else {
#ifdef DEBUG
		printf("No intersection!!!\n");
#endif // DEBUG
	}
}



/****************************鼠标绘制事件******************************************/
void ModelDraw::mousePressEvent(QMouseEvent* _event) {
	/* 【通用】
	* 左键： 设定初始移动点供移动事件使用
	* 右键： 设定初始移动点供移动事件使用
	*/
	Qt::MouseButtons btns = _event->buttons();
	lastpoint2 = _event->pos();
	lastpointok = MapToSphere(lastpoint2, lastpoint3);
	mousemode = _event->button();						// 可以不用这个实现(可改)
	switch (editmodel) {

	case DEFULTEDIT: case DEFORMATION: {
		lastpoint2 = _event->pos();
		makeCurrent();
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		mousemode = _event->button();
		update();
		break;
	}
	case PICKITEM: {// 为框选记录鼠标点击的点
		if (_event->button() == Qt::LeftButton) {
			lastpoint2 = _event->pos();
			double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickEvent(winX, winY);
		}
		break;
	}
	case COORDINATE: {
		/*
		* 左键： 1.cell pick  2.坐标轴功能 分裂
		* 右键：
		*/
		if (btns == Qt::LeftButton) {
			if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
				switch (local_coordinater_axis->get_click_type()) {
					// 分裂
				case CoordinateAxis::ClickType::X_CUBE_CLICK:
					solutionSkel->spiltX();
					drawNormalSkel();
					break;
				case CoordinateAxis::ClickType::Y_CUBE_CLICK:
					solutionSkel->spiltY();
					drawNormalSkel();
					break;
				case CoordinateAxis::ClickType::Z_CUBE_CLICK:
					solutionSkel->spiltZ();
					drawNormalSkel();
					break;
				default:
					break;
				}
			}
			if (local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
				pickEvent(winX, winY);
				update();
			}
		}
		break;
	}
	case SET_CELL_SKEL: {
		/*
		* 左键： 1.cell pick(之前没选中，左点选中：选点) 2.cell add(之前选中过，左点到空白：增加点) 3.cell insert(之前选中过，左点选中另一个：增加点)
		* 右键： 1.cell not pick(之前选中过，右键点到空白处：取消选点)  2.cell del(之前点选中过，且再次右键选中：删除它)
		*/
		if (btns == Qt::LeftButton) {
			double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
			BMesh::Node* pre = solutionSkel->selectNode;
			initializePick();
			if (pre) {
				if (pickEvent_SKEL(winX, winY)) //cell insert
				{
					solutionSkel->rinsertNode(pre);
					selected_c_indices.clear();
				}
				else { // cell add，以之前点选的深度作为深度
					double p[3];
					double temp[3] = { pre->position.x, pre->position.y, pre->position.z };
					pick_point_with_depth(winX, winY, p, temp);
					solutionSkel->raddNode(p[0], p[1], p[2], 0, 0, 0, radius * 0.05, radius * 0.05, radius * 0.05);
				}
			}
			else {
				pickEvent_SKEL(winX, winY);
			}
		}
		else {
			double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内
			BMesh::Node* pre = solutionSkel->selectNode;
			initializePick();
			if (pre) {
				if (pickEvent_SKEL(winX, winY)) {	// 之前点选中过，且再次右键选中，删除它
					if (pre->getPosition() == solutionSkel->selectNode->getPosition())
						solutionSkel->rdelNode();
				}
				initializePick();
			}
		}
		break;
	}
	case CON_CELL_SKEL: {
		/*
		* 左键： 1.cell pick(点到空白处要取消选点)  2.在已有选中后仍选中，连接两点或取消连接两点
		* 右键：
		*/
		if (btns == Qt::LeftButton) {
			double winX = (double)(_event->x() * 2.f - width()) / (width());	// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());	// 窗口坐标，范围在[-1,1]之内

			BMesh::Node* pre = solutionSkel->selectNode;
			initializePick();
			if (pickEvent_SKEL(winX, winY) && pre) {
				solutionSkel->smartConnectNodes(pre, solutionSkel->selectNode);
			}
			update();
		}
		break;
	}
					  /*
					  case PICK_CELL_SKEL: {
						  if (_event->button() == Qt::LeftButton) {
							  mousemode = Qt::LeftButton;
							  double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
							  double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
							  pickEvent(winX, winY);
						  }
						  else {
							  SetEditModel(ModelDraw::DEFULTEDIT);
							  setPickType(1);
						  }
						  break;
					  }

					  case MOVE_CELL_SKEL: {
						  if (_event->button() == Qt::LeftButton) {
							  mousemode = Qt::LeftButton;
							  double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
							  double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
							  pickEvent(winX, winY);
						  }
						  else {
							  SetEditModel(ModelDraw::DEFULTEDIT);
							  setPickType(1);
						  }
						  break;
					  }

					  case PICK_CELLS_SKEL: {
						  if (_event->button() == Qt::LeftButton) {
							  double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
							  double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
							  pickEvent(winX, winY);
						  }
						  else {
							  SetEditModel(ModelDraw::DEFULTEDIT);
							  setPickType(1);
						  }
						  break;
					  }

					  case SET_CELL_SKEL: {
						  if (_event->button() == Qt::LeftButton) {
							  double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
							  double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
							  pickEvent(winX, winY);
						  }
						  else {
							  SetEditModel(ModelDraw::DEFULTEDIT);
							  setPickType(1);
						  }
						  break;
					  }

					  case ROTATE_CELL_SKEL: {
						  if (_event->button() == Qt::LeftButton) {
							  double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
							  double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
							  pickEvent(winX, winY);
						  }
						  else {
							  SetEditModel(ModelDraw::DEFULTEDIT);
							  setPickType(0);
						  }
						  break;
					  }
					  */
	case BSPLINESOLIDEDIT:
	{
		//std::cout << "mouse press" << std::endl;
		lastpoint2 = _event->pos();
		makeCurrent();
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		mousemode = _event->button();
		update();
		if (mousemode == Qt::LeftButton) {
			if (ifchooseface)
			{
				double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
				pickPoint_chooseface(winX, winY);
				break;
			}
		}
		if (mousemode == Qt::RightButton) {//鼠标右键点击选择移动的控制点
			controlLineOrMesh = 2;
			Is_Highlight.clear();//清除之前的高亮
			double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickPoint(winX, winY);
		}
	}
	break;

	}


}

void ModelDraw::Rotation(const QPoint& p) {
	Vector3d  newPoint3D;
	bool newPoint_hitSphere = MapToSphere(p, newPoint3D);
	if (newPoint_hitSphere) {
		Vector3d axis = lastpoint3.Cross(newPoint3D);// corss product
		if (axis.Length() < 1e-7) {
			axis = Vector3d(1, 0, 0);
		}
		else {
			axis.Normalize();
		}
		// find the amount of rotation
		Vector3d d = lastpoint3 - newPoint3D;
		double t = 0.5 * d.Length() / trackballradius;
		if (t < -1.0) t = -1.0;
		else if (t > 1.0) t = 1.0;
		double phi = 2.0 * asin(t);
		double angle = phi * 180.0 / M_PI;

		Rotate(axis, angle);
	}
}

void ModelDraw::mouseMoveEvent(QMouseEvent* _event) {

	switch (editmodel) {
	case DEFULTEDIT: {
		QPoint newPoint2D = _event->pos();
		makeCurrent();
		if (lastpointok) {
			switch (mousemode) {
			case Qt::LeftButton:
				if (drawmodule == SPARSEOPTIMIZATION) return;
				ModelDraw::Rotation(newPoint2D);
				break;
			case Qt::RightButton:
				Translation(newPoint2D);
				break;
			default:
				break;
			}
		}
		lastpoint2 = newPoint2D;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
		break;
	}
	case PICKITEM: {

		if (_event->buttons() == Qt::LeftButton) {

			mouse_move_flag = true;
			// TODO: 绘制框选区域
			QVector2D point0 = QVector2D((double)(lastpoint2.x() * 2 - width()) / (width()), -(double)(lastpoint2.y() * 2 - height()) / (height()));
			QVector2D point1 = QVector2D((double)(_event->x() * 2 - width()) / (width()), -(double)(_event->y() * 2 - height()) / (height()));
			double x_min = std::min(point0.x(), point1.x()), x_max = std::max(point0.x(), point1.x());
			double y_min = std::min(point0.y(), point1.y()), y_max = std::max(point0.y(), point1.y());
			// std::cout << "ret: (" << x_min << ", " << y_min << ") --> (" << x_max << ", " << y_max << ")\n";
			// 初始化为设备坐标
			std::vector<QVector4D> ret = {
				QVector4D(x_min, y_min, 0.9, 1.0), QVector4D(x_min, y_max, 0.9, 1.0), 
				QVector4D(x_max, y_min, 0.9, 1.0), QVector4D(x_max, y_max, 0.9, 1.0)
			};
			// 转化为世界坐标
			for (auto& pos : ret) {
				pos = modelmatrix.inverted() * viewmatrix.inverted() * projectionmatrix.inverted() * pos;
				pos /= pos.w();
			}
			std::vector<int> draw_indices = {
				0, 1, 0, 2, 3, 1, 3, 2
			};
			selected_rectangle.clear();
			for (auto idx : draw_indices) {
				auto& pos4d = ret[idx];
				QVector3D pos3d(pos4d.x(), pos4d.y(), pos4d.z());
				selected_rectangle.push_back(pos3d);
			}
			update();
		}
		break;
	}
	case COORDINATE: {
		/*
		* 移动：坐标轴PICK
		* 移动左键： 整体旋转   坐标轴触发事件 : 移动, 旋转，拉伸
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();
		switch (btns) {
		case Qt::NoButton: {
			// coordinater PICK
			local_coordinater_axis->PickEvent(*this, newWin.x(), newWin.y());
			break;
		}
		case Qt::LeftButton: {
			switch (local_coordinater_axis->get_click_type()) {
			case CoordinateAxis::ClickType::NO_CLICK: {
				// rotation
				if (lastpointok) {
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
				break;
			}
			case CoordinateAxis::ClickType::X_AXIS_CLICK:
			case CoordinateAxis::ClickType::Y_AXIS_CLICK:
			case CoordinateAxis::ClickType::Z_AXIS_CLICK: {
				if (lastpointok) {
					double x, y, z;
					local_coordinater_axis->get_orign(x, y, z);
					BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(x, y, z));
					// 移动
					if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::TRANSFROM) {
						QVector3D dir = local_coordinater_axis->MoveEvent(*this, oldWin, newWin);
						solutionSkel->rXyzCubeByMouse(x + dir.x(), y + dir.y(), z + dir.z(), rlt->root->getLength(), rlt->root->getWidth(), rlt->root->getHeight());
					}
					// 拉伸
					else if (local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
						QVector3D dir = local_coordinater_axis->FlexibleEvent(*this, oldWin, newWin);
						solutionSkel->rXyzCubeByMouse(x, y, z, rlt->root->getLength() + dir.y(), rlt->root->getWidth() + dir.x(), rlt->root->getHeight() + dir.z());
					}
					lastpoint2 = newPoint2D;
				}
				break;
			}
			case CoordinateAxis::ClickType::X_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Y_CIRCLE_CLICK:
			case CoordinateAxis::ClickType::Z_CIRCLE_CLICK: {
				// 旋转
				if (lastpointok) {
					QVector4D quat = local_coordinater_axis->RotateEvent(*this, oldWin, newWin);
					solutionSkel->rRotateCubeByMouse(V3f(quat[0], quat[1], quat[2]), quat[3]);
					lastpoint2 = newPoint2D;
				}
				break;
			}
			default:
				break;
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok && local_coordinater_axis->get_click_type() == CoordinateAxis::ClickType::NO_CLICK) {
				// translation
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}
	case PICK_CELLS_SKEL: {
		if (_event->button() == Qt::LeftButton) {
			double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickEvent(winX, winY);
		}
		else {
			SetEditModel(ModelDraw::DEFULTEDIT);
			setPickType(1);
		}
		break;
	}
	case SET_CELL_SKEL: {
		/*
		* 移动：
		* 移动左键： （未选中）整体旋转   （选中）复制生成
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();
		switch (btns) {
		case Qt::LeftButton: {
			if (solutionSkel->selectNode) {
				BMesh::Node* sNode = solutionSkel->selectNode;
				double temp[3] = { sNode->position.x, sNode->position.y, sNode->position.z };
				if (copyTime) {
					// 移动复制品
					double p[3];
					pick_point_with_depth(newWin.x(), newWin.y(), p, temp);
					solutionSkel->rXyzCubeByMouse(p[0], p[1], p[2], sNode->length, sNode->width, sNode->height);
				}
				else {
					// 生成一个复制品
					solutionSkel->rcopyNode();
					selected_c_indices.clear();
					selected_c_indices.push_back(findNearestCell(temp));
					copyTime = true;
				}
			}
			else {
				if (lastpointok) {
					Rotation(newPoint2D);
					lastpoint2 = newPoint2D;
					lastpointok = MapToSphere(lastpoint2, lastpoint3);
				}
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok) {
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}
	case CON_CELL_SKEL: {
		/* 【基础版】
		* 移动：
		* 移动左键： 整体旋转
		* 移动右键： 整体平移
		*/
		QPoint newPoint2D = _event->pos();
		QVector2D oldWin = QVector2D((float)(lastpoint2.x() * 2.f - width()) / (width()), -(float)(lastpoint2.y() * 2.f - height()) / (height()));
		QVector2D newWin = QVector2D((float)(newPoint2D.x() * 2.f - width()) / (width()), -(float)(newPoint2D.y() * 2.f - height()) / (height()));
		Qt::MouseButtons btns = _event->buttons();
		switch (btns) {
		case Qt::LeftButton: {
			if (lastpointok) {
				Rotation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		case Qt::RightButton: {
			if (lastpointok) {
				Translation(newPoint2D);
				lastpoint2 = newPoint2D;
				lastpointok = MapToSphere(lastpoint2, lastpoint3);
			}
			break;
		}
		default:
			break;
		}
		break;
	}

	case DEFORMATION: {// 实时渲染变形结果
		QPoint newPoint2D = _event->pos();
		if (mousemode == Qt::LeftButton) {
			if (selected_v_indices.empty()) {
				QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Select No Vertex! "));
				return;
			}
			UpdatePickMoveMatrix(newPoint2D);
			std::vector<std::pair<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex>> move_para;
			for (auto& v_idx : selected_v_indices) {
				MeshKernel::iGameVertexHandle vh(v_idx);
				auto& v = surface_mesh == nullptr ? tetMesh->vertices(vh) : surface_mesh->vertices(vh);
				QVector4D current_pos(v.x(), v.y(), v.z(), 1.f);
				QVector4D target_pos = pickmovematrix * current_pos;
				MeshKernel::iGameVertex new_v(target_pos[0], target_pos[1], target_pos[2]);
				move_para.emplace_back(vh, new_v);
				std::cout << "Move Information: (" << v.x() << ", " << v.y() << ", " << v.z() << ") move to (" <<
					target_pos[0] << ", " << target_pos[1] << ", " << target_pos[2] << ")" << std::endl;
			}
			if (surface_mesh) tri_deformation->Execute(fixed_indices, move_para);// 曲面网格变形
			else if (tetMesh) {// 四面体网格变形
				std::vector<std::pair<MeshKernel::iGameVertexHandle, MeshKernel::iGameVertex>> fix_para;
				for (auto& fix_idx : fixed_indices) {
					MeshKernel::iGameVertexHandle fix_vh(fix_idx);
					auto& v = tetMesh->vertices(fix_vh);
					fix_para.emplace_back(fix_vh, v);
				}
				tet_deformation->Execute(fix_para, move_para);
			}
			updateDrawData();
		}
		else if (mousemode == Qt::RightButton) {
			Translation(newPoint2D);
		}
		lastpoint2 = newPoint2D;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
		break;
	}
	case BSPLINESOLIDEDIT:
	{
		QPoint newPoint2D = _event->pos();
		makeCurrent();
		if (lastpointok) {
			switch (mousemode) {
			case Qt::LeftButton:
				Rotation(newPoint2D);
				break;
			case Qt::RightButton:
			{
				//std::cout << "mouse move" << std::endl;
				double winX = (double)(_event->x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
				double winY = -(double)(_event->y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
				double p[3];
				pick_point_with_depth(winX, winY, p, pickedItemPosition);
				CPoint old_pick_point(pickedItemPosition[0], pickedItemPosition[1], pickedItemPosition[2]);//保存move过程中上一帧的顶点
				vector<int> tmp_handle_vector;
				tmp_handle_vector = PointBelongToSolid[old_pick_point];//暂存这个顶点关联的体
				CPoint last_pick_point = old_pick_point;//保存move后的位置
				for (int u = 0; u < tmp_handle_vector.size(); u++) {//遍历所有片
					vector<vector<vector<CPoint>>> tmp = bsplinesolidmesh[u].getControlPoints();
					bool flag = 1;
					//三重循环找点
					for (int i = 0; i < tmp.size(); i++) {
						if (flag == 0)
							break;
						for (int j = 0; j < tmp[i].size(); j++) {
							if (flag == 0)
								break;
							for (int k = 0; k < tmp[i][j].size(); k++) {
								//std::cout << "i" << "j" << 'k' << std::endl;
								if (flag == 0)
									break;
								if (tmp[i][j][k] == old_pick_point) {//体中找到该点，赋新值
									tmp[i][j][k] = CPoint(p[0], p[1], p[2]);
									last_pick_point = tmp[i][j][k];//暂存这一步移动的结果
									Is_Highlight[tmp[i][j][k]] = 1;//新的点高亮
									AllControlPoint.insert(tmp[i][j][k]);
									PointBelongToSolid[tmp[i][j][k]] = tmp_handle_vector;//新的点对应的handle_vector与原来点相同
									bsplinesolidmesh[u].setControlPoints(tmp);

									flag = 0;//在该片中找到了该点，可以跳出三重循环
								}
							}
						}
					}

				}
				//删除上一帧顶点的所有关联
				Is_Highlight.erase(CPoint(pickedItemPosition[0], pickedItemPosition[1], pickedItemPosition[2]));
				AllControlPoint.erase(CPoint(pickedItemPosition[0], pickedItemPosition[1], pickedItemPosition[2]));
				PointBelongToSolid.erase(CPoint(pickedItemPosition[0], pickedItemPosition[1], pickedItemPosition[2]));
				//上一帧顶点位置赋新值
				pickedItemPosition[0] = last_pick_point.m_dX;
				pickedItemPosition[1] = last_pick_point.m_dY;
				pickedItemPosition[2] = last_pick_point.m_dZ;
				if (isRealtime)//实时渲染选项
					bsplineSolidToDiscreteXJ();

				bsplinesolid_drawControlMesh();
				update();
			}
			break;
			default:
				break;
			}
		}
		lastpoint2 = newPoint2D;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
	}
	break;

	default:
		break;
	}
}

void ModelDraw::translationCell(QPoint p) {
	double winX = (double)(p.x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
	double winY = -(double)(p.y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
	double pLocal[3] = { 0 };
	// gluUnProject(winX, winY, depth, modelview_matrix_, projection_matrix_, viewport, &pLocal[0], &pLocal[1], &pLocal[2]);
	double item[3]{ px, py, pz };
	pick_point_with_depth(winX, winY, pLocal, item);
	double relativePosition[3] = { 0 };
	relativePosition[0] = pLocal[0] - px;
	relativePosition[1] = pLocal[1] - py;
	relativePosition[2] = pLocal[2] - pz;
	px = pLocal[0]; py = pLocal[1]; pz = pLocal[2];
	//SingalFather::getInstance()->sdeltaEditCell(relativePosition[0], relativePosition[1], relativePosition[2], -1);
}



void ModelDraw::mouseReleaseEvent(QMouseEvent* _event) {
	//assert(mousemode < N_MOUSE_MODES);
	switch (editmodel) {
	case PICKITEM: {// 框选功能 Beta
		selected_rectangle.clear();
		QPoint currentpoint2 = _event->pos();
		int x_min = std::min(lastpoint2.x(), currentpoint2.x()), x_max = std::max(lastpoint2.x(), currentpoint2.x());
		int y_min = std::min(lastpoint2.y(), currentpoint2.y()), y_max = std::max(lastpoint2.y(), currentpoint2.y());
		int red, green, blue, sum_color_threshold = 75, signal_color_threshold = 30;// threshold 越大, 采样点越少, 速度越快, 但可能出现精度问题
		QScreen* screen = QApplication::primaryScreen();
		QPixmap pixmap = screen->grabWindow(this->winId(), x_min, y_min, x_max - x_min + 1, y_max - y_min + 1);
		QImage image = pixmap.toImage();
		image.save("Image/captureImg.jpg", "jpg"); // 保存截图
		int w = image.width(), h = image.height();
		vector<vector<bool>> visited(w, vector<bool>(h, false));
		vector<QPoint> dirs = {
			QPoint(-1, 0), QPoint(1, 0), QPoint(0, -1), QPoint(0, 1),
			QPoint(-1, -1), QPoint(-1, 1), QPoint(1, -1), QPoint(1, 1)
		};
		vector<QPoint> centroids;
		for (int i = 0; i < w; ++i) {
			for (int j = 0; j < h; ++j) {
				if (visited[i][j]) continue;
				queue<pair<int, int>> que;
				que.emplace(i, j);
				visited[i][j] = true;
				QPoint sum(i, j);
				int cnt = 1;
				while (!que.empty()) {
					auto pospair = que.front();
					que.pop();
					int x = pospair.first, y = pospair.second;
					QPoint pos(x, y);
					const QColor& color = image.pixelColor(pos);
					red = color.red(), green = color.green(), blue = color.blue();
					for (auto& move : dirs) {
						QPoint adj = pos + move;
						int _x = adj.x(), _y = adj.y();
						if (_x < 0 || _x >= w || _y < 0 || _y >= h || visited[_x][_y]) continue;
						const QColor& _c = image.pixelColor(adj);
						int diff_r = abs(_c.red() - red);
						int diff_g = abs(_c.green() - green);
						int diff_b = abs(_c.blue() - blue);
						int diff_signal = std::max(diff_r, std::max(diff_g, diff_b));
						int diff_sum = diff_r + diff_g + diff_b;
						if (diff_signal < signal_color_threshold && diff_sum < sum_color_threshold) {
							que.emplace(_x, _y);
							visited[_x][_y] = true;
							sum += adj;
							cnt++;
						}
					}
				}
				sum /= cnt;
				sum.setX(sum.x() + x_min);
				sum.setY(sum.y() + y_min);
				centroids.push_back(sum);
			}
		}

		time_t time1 = clock();
		for (auto& pos : centroids) {
			double winX = (double)(pos.x() * 2.f - width()) / (width());// 窗口坐标，范围在[-1,1]之内
			double winY = -(double)(pos.y() * 2.f - height()) / (height());// 窗口坐标，范围在[-1,1]之内
			pickEvent(winX, winY);
		}
		int pick_cost = clock() - time1;
		// cout << "Sample points size = " << centroids.size() << ", cost " << pick_cost << "ms.\n";
		lastpoint2 = currentpoint2;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
		break;
	}
	case SET_CELL_SKEL: {
		/*
		* 释放左键： （选中）将  复制时间与选中 重置
		* 释放右键：
		*/
		if (copyTime) {
			copyTime = false;
			initializePick();
		}
		break;
	}
	case BSPLINESOLIDEDIT: {
		QPoint newPoint2D = _event->pos();
		makeCurrent();
		if (lastpointok) {
			switch (mousemode) {
			case Qt::LeftButton:
				Rotation(newPoint2D);
				break;
			case Qt::RightButton:
			{
				std::cout << "mouse release" << std::endl;
				bsplineSolidToDiscreteXJ();
				bsplinesolid_drawControlMesh();

			}
			break;
			default:
				break;
			}
		}
		lastpoint2 = newPoint2D;
		lastpointok = MapToSphere(lastpoint2, lastpoint3);
		update();
	}
	default:
		break;
	}
	mouse_move_flag = false;
	lastpointok = false;
	mousemode = Qt::NoButton;
}

void ModelDraw::wheelEvent(QWheelEvent* _event) {
	// 滚轮拉伸
	if (editmodel == COORDINATE && local_coordinater_axis->get_mode_type() == CoordinateAxis::ModeType::EDIT) {
		double x, y, z;
		local_coordinater_axis->get_orign(x, y, z);
		BMesh::Node* rlt = solutionSkel->selectWhichNode(V3f(x, y, z));
		double d = -_event->delta() / 120.0 * 0.01 * radius;
		solutionSkel->rXyzCubeByMouse(x, y, z, rlt->root->getLength() + d, rlt->root->getWidth() + d, rlt->root->getHeight() + d);
	}
	else {
		// 0.1代表缩放速度
		double d = -_event->delta() / 120.0 * 0.1 * radius;
		Translate(Vector3d(0.0, 0.0, d));
		local_coordinater_axis->set_len(local_coordinater_axis->get_len() + d * 0.1);
	}
	update();
}

void ModelDraw::keyPressEvent(QKeyEvent* _event) {
	switch (_event->key()) {
	case Qt::Key_Escape:
		qApp->quit();
		break;
	case Qt::Key_Delete:
	case Qt::Key_D:
		deleteSelectedAll();
		break;
	case Qt::Key_A:
		add_new_cell();
		break;
	case Qt::Key_Control:
		ctlPressed = true;
		break;
	}
	_event->ignore();
}

void ModelDraw::keyReleaseEvent(QKeyEvent* _event) {
	if (_event->key() == Qt::Key_Control) {
		ctlPressed = false;
	}
	_event->ignore();
}

void ModelDraw::SetEditModel(const EditModel& editmodel) {
	this->editmodel = editmodel;
}


/**************************** 渲染函数 ******************************************/
void ModelDraw::EnableLighting(bool b) {
	isEnableLighting = b;
	update();
}

void ModelDraw::LoadRotation(void) {
	LoadCopyModelViewMatrix();
	update();
}

void ModelDraw::ResetView(void) {
	ResetModelviewMatrix();
	viewCenter();
	update();
}

void ModelDraw::getPatchStatus(bool status) {
	this->PatchChange = status;
}

void ModelDraw::viewCenter() {
	update();
}

void ModelDraw::ModelClear() {

	// 清除绘制数据
	this->vertices_data.clear();
	vector<GLfloat>().swap(this->vertices_data);
	tri_cnt = 0;
	quad_cnt = 0;

	// 清除模型数据
	delete this->hexMesh;
	this->hexMesh = nullptr;
	delete this->surface_mesh;
	this->surface_mesh = nullptr;
	delete this->tetMesh;
	this->tetMesh = nullptr;
	delete this->skeletal_mesh;
	this->skeletal_mesh = nullptr;

	// 清除稀疏优化数据
	delete this->pipeline;
	this->pipeline = nullptr;
	delete this->poly_mesh;// OpenMesh
	this->poly_mesh = nullptr;
	bdy.clear();
	axis.clear();


	// 清除点选数据
	this->selected_v_indices.clear();
	this->selected_e_indices.clear();
	this->selected_f_indices.clear();
	this->selected_c_indices.clear();
	vector<int>().swap(this->selected_v_indices);
	vector<int>().swap(this->selected_e_indices);
	vector<int>().swap(this->selected_f_indices);
	vector<int>().swap(this->selected_c_indices);
	mouse_move_flag = false;

	SetDrawModule(DrawModule::EMPTY);// 加上此行可以避免部分 Bug

}

void ModelDraw::genColors(size_t sz) {
	colors.resize(sz);
	srand((unsigned)time(NULL));
	for (int i = 0; i < sz; ++i) {
		double R = ((std::rand() % 9999) + 1) / 10000.f;
		double G = ((std::rand() % 9999) + 1) / 10000.f;
		double B = ((std::rand() % 9999) + 1) / 10000.f;
		colors[i] = Eigen::Vector3d(R, G, B);
	}
}

void ModelDraw::initColorbar() {
	std::vector<Vector3d> tmp_colors = {
		Vector3d(6, 48, 254), Vector3d(10, 97, 255), Vector3d(12, 147, 255), Vector3d(17, 195, 255), Vector3d(18, 246, 255),
		Vector3d(19, 245, 205), Vector3d(18, 245, 153), Vector3d(17, 246, 103), Vector3d(16, 246, 50), Vector3d(15, 247, 1),
		Vector3d(63, 249, 0), Vector3d(108, 250, 0), Vector3d(156, 253, 0), Vector3d(205, 255, 0), Vector3d(253, 254, 1),
		Vector3d(249, 207, 0), Vector3d(245, 158, 0), Vector3d(243, 108, 0), Vector3d(239, 59, 0), Vector3d(237, 9, 0)
	};
	int sz = tmp_colors.size();
	double step = 1.f / (sz - 1), value = 0.f;
	colorbar.clear();
	colorbar.reserve(sz);
	for (int i = 0; i < sz; ++i) {
		tmp_colors[i] /= 255.f;
		colorbar.emplace_back(value, tmp_colors[i]);
		value += step;
	}

}

Vector3d ModelDraw::getColor(double value) {// value more larger color more red, more smaller more blue.

	if (value < 0.f) value = 0.f;
	if (value > 1.f) value = 1.f;
	Vector3d color(0, 0, 0);
	//if (value < 0.5) {
	//	color[2] = 1 - value * 2;
	//	color[1] = value * 2;
	//}
	//else {
	//	color[0] = (value - 0.5) * 2;
	//	color[1] = (1 - value) * 2;
	//}
	
	for (int i = 1; i < colorbar.size(); ++i) {
		if (value <= colorbar[i].first || i == colorbar.size() - 1) {
			double pre_value = colorbar[i - 1].first;
			double next_value = colorbar[i].first;
			Vector3d pre_color = colorbar[i - 1].second;
			Vector3d next_color = colorbar[i].second;
			color = pre_color * (next_value - value) + next_color * (value - pre_value);
			color /= (next_value - pre_value);
			break;
		}
	}
	return color;

}

void ModelDraw::loadShader() {

	// 相关数据传入sahder
	shaderProgram.setUniformValue("model", modelmatrix);	// 把模型矩阵传入shader
	shaderProgram.setUniformValue("view", viewmatrix); // 视角矩阵
	shaderProgram.setUniformValue("projection", projectionmatrix); // 投影矩阵
	shaderProgram.setUniformValue("rotate", rotatematrix); // 旋转矩阵（法向量专用）
	shaderProgram.setUniformValue("viewpos", 0.0f, 0.0f, 1.0f);  // 视角方向
	shaderProgram.setUniformValue("light.position", 0.0f, 0.0f, 1.0f);  // 光照方向
	shaderProgram.setUniformValue("light.ambient", 0.4f, 0.4f, 0.4f);	// 环境光照系数
	shaderProgram.setUniformValue("light.diffuse", 0.5f, 0.5f, 0.5f);	// 漫反射系数
	shaderProgram.setUniformValue("light.specular", 0.5f, 0.5f, 0.5f);	// 镜面反射系数
	shaderProgram.setUniformValue("shininess", 32.0f);

	// 释放旧数据的内存
	core->glDeleteVertexArrays(1, &VAO);
	core->glDeleteBuffers(1, &VBO);
	// 1.新建顶点数组对象
	core->glGenVertexArrays(1, &VAO);
	core->glGenBuffers(1, &VBO);
	// 2. 绑定顶点数组对象
	core->glBindVertexArray(VAO);
	// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
	core->glBindBuffer(GL_ARRAY_BUFFER, VBO);
	core->glBufferData(GL_ARRAY_BUFFER, vertices_data.size() * sizeof(GLfloat), vertices_data.data(), GL_STATIC_DRAW);
	// 4. 设定顶点属性指针
	// 顶点坐标传入VBO，location = 0
	core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)0);
	core->glEnableVertexAttribArray(0);
	// 顶点 法向/颜色 传入VBO，location=1
	core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 6 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
	core->glEnableVertexAttribArray(1);
	core->glBindVertexArray(0);
	core->glBindBuffer(GL_ARRAY_BUFFER, 0);

//``````````````````````````````````````````````````````````````````````

//// 释放旧数据的内存
//    core->glDeleteVertexArrays(1, &VAO);
//    core->glDeleteBuffers(1, &VBO);
//
//// 1. 新建顶点数组对象
//    core->glGenVertexArrays(1, &VAO);
//    core->glGenBuffers(1, &VBO);
//
//// 2. 绑定顶点数组对象
//    core->glBindVertexArray(VAO);
//
//// 3. 把我们的顶点数组复制到一个顶点缓冲中，供OpenGL使用
//    core->glBindBuffer(GL_ARRAY_BUFFER, VBO);
//    core->glBufferData(GL_ARRAY_BUFFER, vertices_data.size() * sizeof(GLfloat), vertices_data.data(), GL_STATIC_DRAW);
//
//// 4. 设定顶点属性指针
//// 顶点坐标传入VBO，location = 0
//    core->glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (void*)0);
//    core->glEnableVertexAttribArray(0);
//
//// 顶点法向量传入VBO，location = 1
//    core->glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (void*)(3 * sizeof(GLfloat)));
//    core->glEnableVertexAttribArray(1);
//
//// 顶点颜色传入VBO，location = 2
//    core->glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, 9 * sizeof(GLfloat), (void*)(6 * sizeof(GLfloat)));
//    core->glEnableVertexAttribArray(2);
//
//// 解绑 VAO
//    core->glBindVertexArray(0);
//    core->glBindBuffer(GL_ARRAY_BUFFER, 0);


}

void ModelDraw::updateDrawData_HexMesh() {
	if (!hexMesh) return;
	vertices_data.clear();
    hexMesh->generateAllBoundaryVertexNormals();

	auto& cells = hexMesh->allcells();
    auto& edges = hexMesh->alledges();

	// 更新绘制数据
	unordered_map<int, bool> isDrawed;
	for (auto& cp : cells) {
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		auto fhs = cp.second.getFaceHandle();
		for (auto& fh : fhs) {
			if (isDrawed.count(fh)) continue;
			isDrawed[fh] = true;
			auto vhs = hexMesh->faces(fh).getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = hexMesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(0.349);
				vertices_data.push_back(0.6);
				vertices_data.push_back(0.349);
			}
		}
	}
	quad_cnt = vertices_data.size() / 24;
	vertices_data.insert(vertices_data.end(), compared_vertices_data.begin(), compared_vertices_data.end());
	used_last_vcnt = compared_vertices_data.size() / 6;

	update();

}

void ModelDraw::updateDrawData_HexMesh_Analysis() {
	if (!hexMesh) return;

	vertices_data.clear();

	auto& cells = hexMesh->allcells();
	std::cout << "cells.size(): " << cells.size() << std::endl;
	if (analysis_mode < 3) {// 等几何分析
		unordered_map<int, double> vex_val;
		if (analysis_flag && !analysis_data[0].empty() && !analysis_data[1].empty() && !analysis_data[2].empty()) {
			auto& data = analysis_data[analysis_mode];
			double min_val(Double_MAX), max_val(Double_MIN);
			//for (auto& cp : cells) {
			for (int i = 0; i < hexMesh->csize(); ++i) {
				//std::cout << "i: " << std::endl;
				MeshKernel::iGameCellHandle ch(i);
				auto& c = hexMesh->cells(ch);
				
				if (!hexMesh->isOnBoundary(ch)) continue;
				//std::cout << "isOnBoundary ok " << std::endl;
				auto vhs = c.getVertexHandle();
				//std::cout << "iGame ch = " << ch << ", vhs: ";
				//for (auto& vh : vhs) {
				for (int j = 0; j < 8; ++j) {
					MeshKernel::iGameVertexHandle vh(vhs[j]);
					//std::cout << vh << " ";
					if (!vex_val.count(vh.idx())) {
						vex_val[vh] = data[vh];
						min_val = std::min(min_val, vex_val[vh]);
						max_val = std::max(max_val, vex_val[vh]);
						//std::cout << "vex val: " << vex_val[vh] << std::endl;
					}
				}
				//std::cout << std::endl;
			}
			double range_val = max_val - min_val;
			range_val = max_val - min_val;
			std::cout << "max_val = " << max_val << ", min_val = " << min_val << std::endl;
			for (auto& vvp : vex_val) {
				vvp.second = (vvp.second - min_val) / (max_val - min_val);
			}
		}

		//std::cout << " 更新绘制数据" << std::endl;
		// 更新绘制数据

		unordered_map<int, bool> isDrawed;
		for (auto& cp : cells) {
			if (!hexMesh->isOnBoundary(cp.first)) continue;
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				if (isDrawed.count(fh)) continue;
				isDrawed[fh] = true;
				auto vhs = hexMesh->faces(fh).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					Vector3d color = getColor(vex_val[vh]);
					//std::cout << v.x() << " " << v.y() << " " << v.z() << std::endl;
					//std::cout << color[0] << " " << color[1] << " " << color[2] << std::endl;
					vertices_data.push_back(v.x());
					vertices_data.push_back(v.y());
					vertices_data.push_back(v.z());
					vertices_data.push_back(color[0]);
					vertices_data.push_back(color[1]);
					vertices_data.push_back(color[2]);
				}
			}
		}

	}
	else {//  Scaled Jacobian 体质量分析

		vector<double> scaled_jacobian;
		Analysis_HexMesh app(*hexMesh);
		app.calculateCellQuality_Jacobian(scaled_jacobian);
		std::cout << " scaled_jacobian size: " << scaled_jacobian.size() << ", cells size: " << cells.size() << std::endl;
		int idx = -1;
		unordered_map<int, bool> isDrawed;
		for (auto& cp : cells) {
			idx++;
			if (!hexMesh->isOnBoundary(cp.first)) continue;
			Vector3d color = getColor(1 - scaled_jacobian[idx]);
			auto fhs = cp.second.getFaceHandle();
			for (auto& fh : fhs) {
				if (isDrawed.count(fh)) continue;
				isDrawed[fh] = true;
				auto vhs = hexMesh->faces(fh).getVertexHandle();
				for (auto& vh : vhs) {
					auto& v = hexMesh->vertices(vh);
					vertices_data.push_back(v.x());
					vertices_data.push_back(v.y());
					vertices_data.push_back(v.z());
					vertices_data.push_back(color[0]);
					vertices_data.push_back(color[1]);
					vertices_data.push_back(color[2]);
				}
			}
		}

	}
	std::cout << " ModelDraw 4778" << std::endl;
	quad_cnt = vertices_data.size() / 24;
	vertices_data.insert(vertices_data.end(), compared_vertices_data.begin(), compared_vertices_data.end());
	used_last_vcnt = compared_vertices_data.size() / 6;
	update();

}

void ModelDraw::updateDrawData_SurfaceMesh() {
	if (!surface_mesh) return;
	vertices_data.clear();
	vertices_data.reserve(surface_mesh->FaceSize() * 18);
	surface_mesh->genAllVerticesNormal();
	//std::cout << "Generate all vertices' normals success." << std::endl;
	std::vector<GLfloat> quad_data;// 保存四边形的数据
	// 更新绘制数据
	for (auto& fp : surface_mesh->allfaces()) {
		auto vhs = fp.second.getVertexHandle();
		if (vhs.size() == 3) {
			for (auto& vh : vhs) {
				auto& v = surface_mesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
                vertices_data.push_back(v.getNormalX());
                vertices_data.push_back(v.getNormalY());
                vertices_data.push_back(v.getNormalZ());
			}
		}
		else if (vhs.size() == 4) {
			for (auto& vh : vhs) {
				auto& v = surface_mesh->vertices(vh);
				quad_data.push_back(v.x());
				quad_data.push_back(v.y());
				quad_data.push_back(v.z());
				quad_data.push_back(v.getNormalX());
				quad_data.push_back(v.getNormalY());
				quad_data.push_back(v.getNormalZ());
			}
		}
		else {
			std::cerr << "Wrong: exist face that is neither a triangle nor a quadrangle. vertices = " << vhs.size() << std::endl;
		}

	}

	tri_cnt = vertices_data.size() / 18;
	quad_cnt = quad_data.size() / 24;
	//std::cout << "The model has " << tri_cnt << " triangles and " << quad_cnt << " quadrangles." << std::endl;

	vertices_data.insert(vertices_data.end(), quad_data.begin(), quad_data.end());

	update();

}

void ModelDraw::updateDrawData_SurfaceMesh_Analysis() {
	if (!surface_mesh || quad_cnt != 0) return;
	vertices_data.clear();
	vertices_data.reserve(surface_mesh->FaceSize() * 18);
	Analysis_TriMesh app(*surface_mesh);
	if (analysis_mode == 1) {// angle
		vector<double> face_q;
		app.calculateTriangleQuality(face_q);
		std::cout << "quality size: " << face_q.size() << ", faces size: " << surface_mesh->fsize() << std::endl;
		// 更新绘制数据
		int fcnt = 0;
		for (auto& fp : surface_mesh->allfaces()) {
			auto vhs = fp.second.getVertexHandle();
			Vector3d color = getColor(face_q[fcnt++]);
			for (auto& vh : vhs) {
				auto& v = surface_mesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}
		}
	}
	else {// curvature
		vector<double> curvatures;
		app.calculateGaussCurvature(curvatures);
		int vcnt = surface_mesh->vsize();
		std::cout << "curvatures size: " << curvatures.size() << ", vertices size: " << vcnt << std::endl;
		typedef std::pair<int, double> vc;
		std::vector<vc> vec;
		vec.reserve(vcnt);
		for (int i = 0; i < vcnt; ++i) {
			vec.emplace_back(i, curvatures[i]);
		}
		std::sort(vec.begin(), vec.end(), [&](vc& v1, vc& v2) {
			return v1.second < v2.second;
			});
		for (int i = 0; i < vcnt; ++i) {
			curvatures[vec[i].first] = i * 1.0f / vcnt;
			//std::cout << "curvature: " << vec[i].second << std::endl;
		}


		for (auto& fp : surface_mesh->allfaces()) {
			auto vhs = fp.second.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = surface_mesh->vertices(vh);
				Vector3d color = getColor(curvatures[vh]);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(color[0]);
				vertices_data.push_back(color[1]);
				vertices_data.push_back(color[2]);
			}
		}

	}

	update();

}

void ModelDraw::updateDrawData_TriMesh() {
	vertices_data.clear();
	if (!tri_mesh) return;
	if (!tri_mesh->has_face_normals()) {
		tri_mesh->request_face_normals();
	}
	if (!tri_mesh->has_vertex_normals()) {
		tri_mesh->request_vertex_normals();
	}
	tri_mesh->update_face_normals();
	tri_mesh->update_vertex_normals();
	for (auto f_it = tri_mesh->faces_sbegin(); f_it != tri_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto V = tri_mesh->point(*fv_it);
			auto pos = V.data();
			vertices_data.push_back(pos[0]);
			vertices_data.push_back(pos[1]);
			vertices_data.push_back(pos[2]);
			auto vec = tri_mesh->normal(*fv_it).data();
			vertices_data.push_back(vec[0]);
			vertices_data.push_back(vec[1]);
			vertices_data.push_back(vec[2]);
		}

	}

	update();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_1() {
	vertices_data.clear();

	auto propsSke = OpenMesh::FProp<int>(*poly_mesh, "ske_flag");
	auto propsBdy = OpenMesh::FProp<int>(*poly_mesh, "bdy_flag");

	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		int s_flag = propsSke[f_it], b_flag = propsBdy[f_it];
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			if (s_flag == 1) {
				vertices_data.push_back(0.4745);
				vertices_data.push_back(0.6902);
				vertices_data.push_back(0.3608);
			}
			else if (b_flag == -1) {
				vertices_data.push_back(0.2667);
				vertices_data.push_back(0.001f);
				vertices_data.push_back(0.3294);
			}
			else {
				vertices_data.push_back(0.1882);
				vertices_data.push_back(0.4039);
				vertices_data.push_back(0.4353);
			}

		}

	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_2() {
	vertices_data.clear();

	auto propsBdy = OpenMesh::FProp<int>(*poly_mesh, "ske_sp_flag");
	int max_value = -9999999, min_value = 9999999;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		int value = propsBdy[f_it];
		min_value = std::min(min_value, value);
		max_value = std::max(max_value, value);
	}
	double range_value = max_value - min_value;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		int value = propsBdy[f_it];
		Vector3d color = getColor((value - min_value) * 1.f / range_value);
		//OpenMesh::Vec3f N = poly_mesh->normal(*f_it).normalized();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			vertices_data.push_back(color[0]);
			vertices_data.push_back(color[1]);
			vertices_data.push_back(color[2]);
		}

	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::updateDrawData_PolyMesh_OPT_V() {
	vertices_data.clear();

	auto propsBdy = OpenMesh::VProp<double>(*poly_mesh, "mmc_p_val");
	int cnt_val = 0;
	double max_value = -9999999, min_value = 9999999;
	for (auto v_it = poly_mesh->vertices_begin(); v_it != poly_mesh->vertices_end(); ++v_it) {
		double value = propsBdy[v_it];
		if (value < 0.f) value = 0.f;
		min_value = std::min(min_value, value);
		max_value = std::max(max_value, value);
		//if (cnt_val++ % 50 == 0) std::cout << "vertex value: " << value << std::endl;
	}
	//std::cout << "\n\n--------------------\n\n" << std::endl;
	double range_value = max_value - min_value;
	for (auto f_it = poly_mesh->faces_begin(); f_it != poly_mesh->faces_end(); ++f_it) {
		auto fvs = f_it->vertices();
		//OpenMesh::Vec3f N = poly_mesh->normal(*f_it).normalized();
		for (auto fv_it = fvs.begin(); fv_it != fvs.end(); ++fv_it) {
			double value = propsBdy[fv_it];
			//if (cnt_val++ % 50 == 0) std::cout << "vertex value: " << value << std::endl;
			double k = (value - min_value) * 1.f / range_value;
			// 将 0 ~ 0.1 拉长到 0 ~ 0.5, 0.1 ~ 1 压缩到 0.5 ~ 1
			if (k <= 0.1) k *= 5;
			else {
				k = 0.5 + 0.5 * (k - 0.1) / 0.9;
			}
			Vector3d color = getColor(k);
			auto V = poly_mesh->point(fv_it);
			auto P = V.data();
			vertices_data.push_back(P[0]);
			vertices_data.push_back(P[1]);
			vertices_data.push_back(P[2]);
			vertices_data.push_back(color[0]);
			vertices_data.push_back(color[1]);
			vertices_data.push_back(color[2]);
		}

	}

	for (int i = 0; i < bdy.size(); ++i) {
		auto& line = bdy[i];
		auto& color = colors[i % colors.size()];
		for (int j = 1; j < line.size(); ++j) {
			auto& p1 = line[j - 1], & p2 = line[j];
			vertices_data.push_back(p1.x()), vertices_data.push_back(p1.y()), vertices_data.push_back(0.f);
			vertices_data.push_back(color.x()), vertices_data.push_back(color.y()), vertices_data.push_back(color.z());
			vertices_data.push_back(p2.x()), vertices_data.push_back(p2.y()), vertices_data.push_back(0.f);
			vertices_data.push_back(color.x()), vertices_data.push_back(color.y()), vertices_data.push_back(color.z());
		}
	}

#define Handle(Class) opencascade::handle<Class>
	update();

}

void ModelDraw::updateDrawData_TetMesh() {
	if (!tetMesh) return;
	vertices_data.clear();
	vertices_data.reserve(tetMesh->FaceSize() * 18);
	// 更新绘制数据
	tetMesh->genAllVerticesNormal();
	for (auto& fp : tetMesh->allfaces()) {
		if (tetMesh->isOnBoundary(fp.first)) {// 只渲染表面
			auto vhs = fp.second.getVertexHandle();
			for (auto& vh : vhs) {
				auto& v = tetMesh->vertices(vh);
				vertices_data.push_back(v.x());
				vertices_data.push_back(v.y());
				vertices_data.push_back(v.z());
				vertices_data.push_back(v.getNormalX());
				vertices_data.push_back(v.getNormalY());
				vertices_data.push_back(v.getNormalZ());
			}
		}
	}
	update();

}


void ModelDraw::updateDrawData_Surface() {
	vertices_data.clear();

	for (int i = 0; i < data_surfaces.size(); ++i) {
		auto& singleSurface = data_surfaces[i];
		auto& color = colors[i % colors.size()];
		for (auto& row : singleSurface) {
			for (auto& col : row) {
				vertices_data.push_back(col.getX());
				vertices_data.push_back(col.getY());
				vertices_data.push_back(col.getZ());
//                vertices_data.push_back(0.349);
//                vertices_data.push_back(0.6);
//                vertices_data.push_back(0.349);
				vertices_data.push_back(color.x());
				vertices_data.push_back(color.y());
				vertices_data.push_back(color.z());
			}
		}
	}
	update();


//    vertices_data.clear();
//// colors[0] = { 227.0 / 255, 121.0 / 255, 77.0 / 255 };
//    for (int i = 0; i < data_surfaces.size(); ++i) {
//        auto& singleSurface = data_surfaces[i];
//        auto& color = colors[i % colors.size()];
//        auto& normals = surface_normals[i];
//        for (int j = 0; j < singleSurface.size(); j++) {
//            for (int k = 0; k < singleSurface[j].size(); ++k) {
//                auto& v = singleSurface[j][k];
//                auto& normal = normals[j][k];
////                Vex center_c = getCellCenter(ch);
////                Vex center_f = getFaceCenter(fp.first);
////                Vec dir_out = (center_f - center_c).normalized();
//                if(v.getX() * normal.getX() + v.getY() * normal.getY() + v.getZ() * normal.getZ() < 0){
//                    normal.setXYZ(normal.getX() * -1.0, normal.getY() * -1.0, normal.getZ() * -1.0) ;
//                }
//                vertices_data.push_back(v.getX());
//                vertices_data.push_back(v.getY());
//                vertices_data.push_back(v.getZ());
//                vertices_data.push_back(normal.getX());
//                vertices_data.push_back(normal.getY());
//                vertices_data.push_back(normal.getZ());
//                vertices_data.push_back(color[0]);
//                vertices_data.push_back(color[1]);
//                vertices_data.push_back(color[2]);
////                vertices_data.push_back(0.37);
////                vertices_data.push_back(0.6);
////                vertices_data.push_back(0.37);
//            }
//        }
//    }
//    update();

}

void ModelDraw::updateDrawData_Surface_CutPlane() {

	vertices_data.clear();

	double draw_min_x = cutX * (bbMax[0] - bbMin[0]) + bbMin[0];
	double draw_min_y = cutY * (bbMax[1] - bbMin[1]) + bbMin[1];
	double draw_min_z = cutZ * (bbMax[2] - bbMin[2]) + bbMin[2];
	printf("draw min: (%.2f, %.2f, %.2f)\n", draw_min_x, draw_min_y, draw_min_z);

	for (int i = 0; i < data_surfaces.size(); ++i) {
		auto& singleSurface = data_surfaces[i];
		auto& color = colors[i % colors.size()];
		for (int j = 0; j < singleSurface.size(); ++j) {
			vector<CPoint> _drawVertices(singleSurface[j].size());
			bool _drawFlag = true;
			for (int k = 0; k < singleSurface[j].size(); ++k) {
				double x = singleSurface[j][k].getX();
				double y = singleSurface[j][k].getY();
				double z = singleSurface[j][k].getZ();
				if ((x < draw_min_x - 1E-1F) || (y < draw_min_y - 1E-1F) || (z < draw_min_z - 1E-1F)) {
					_drawFlag = false;
					break;
				}
				else {
					_drawVertices[k] = singleSurface[j][k];
				}
			}
			if (!_drawFlag) continue;
			for (auto& v : _drawVertices) {
				vertices_data.push_back(v.getX());
				vertices_data.push_back(v.getY());
				vertices_data.push_back(v.getZ());
				vertices_data.push_back(color.x());
				vertices_data.push_back(color.y());
				vertices_data.push_back(color.z());
			}
		}

	}
	update();

}

void ModelDraw::updateDrawData_Surface_IGES_Color(int dirID) {
    vertices_data.clear();
    for (int i = 0; i < data_surfaces.size(); ++i) {
        auto& singleSurface = data_surfaces[i];
        auto& singleSurfaceColor = data_surfaces_color[dirID][i];
        auto& normal = surface_normals[i];
        for (int i = 0; i < singleSurface.size(); i++) {
            for (int k = 2; k < singleSurface[i].size(); ++k) {
                vector<CPoint> tri_vhs = { singleSurface[i][0], singleSurface[i][k - 1], singleSurface[i][k] };
                vector<Vector3d> tri_color = { singleSurfaceColor[i][0], singleSurfaceColor[i][k - 1], singleSurfaceColor[i][k] };
                for (int id = 0; id < 3; ++ id) {
                    vertices_data.push_back(tri_vhs[id].getX());
                    vertices_data.push_back(tri_vhs[id].getY());
                    vertices_data.push_back(tri_vhs[id].getZ());
                    vertices_data.push_back(tri_color[id](0));
                    vertices_data.push_back(tri_color[id](1));
                    vertices_data.push_back(tri_color[id].operator()(2));
                }
            }
        }
    }
    update();
}

void ModelDraw::updateDrawData(bool need_update_bvhtree) {
    cout << "update draw data" << endl;

	vertices_data.clear();
	quad_cnt = tri_cnt = line_cnt = translucence_tricnt = 0;

	if (cutFlag) {
		updateDrawData_CutPlane();
		return;
	}
	switch (drawmodule) {

	case DrawModule::EMPTY:
		return;

	case DrawModule::SPARSEOPTIMIZATION:
		break;

	case DrawModule::SKELMESH:
		updateDrawData_SkelMesh();
        updateBVH();
		break;

	case DrawModule::TRIMESH:
		updateDrawData_TriMesh();
		break;

	case DrawModule::SURFACEMESH:
		if (analysis_flag) updateDrawData_SurfaceMesh_Analysis();
		else updateDrawData_SurfaceMesh();
		//updateBVH();
		break;

	case DrawModule::TETMESH:
		updateDrawData_TetMesh();
		updateBVH();
		break;

	case DrawModule::HEXMESH:
		if (analysis_flag) updateDrawData_HexMesh_Analysis();
		else updateDrawData_HexMesh();
		updateBVH();
		break;
		/*Temp-lee*/
	case DrawModule::SKEL:
		updateDrawData_HexMesh();
		break;
	case DrawModule::BSPLINESURFACE:
		updateDrawData_Surface();
		//updateBVH();
		break;

	default:
		break;
	}
	//if (need_update_bvhtree) updateBVH();
}




void ModelDraw::updateDrawData_CutPlane() {
	switch (drawmodule) {

	case DrawModule::EMPTY:
		return;
	case DrawModule::SURFACEMESH:
		updateDrawData_SurfaceMesh_CutPlane();
		break;

	case DrawModule::TETMESH:
		updateDrawData_TetMesh_CutPlane();
		break;

	case DrawModule::HEXMESH:
		updateDrawData_HexMesh_CutPlane();
		break;

	case DrawModule::BSPLINESURFACE:
		updateDrawData_Surface_CutPlane();
		break;

	default:
		break;
	}
}

void ModelDraw::update_fixed_indices() {
	fixed_indices = selected_v_indices;
	selected_v_indices.clear();
	std::cout << "fixed: " << fixed_indices.size() << "vertices\n";
	if (surface_mesh) {
		for (auto& fix_idx : fixed_indices) {
			auto& v = surface_mesh->vertices(MeshKernel::iGameVertexHandle(fix_idx));
			printf("vex %d fixed at: (%.2f, %.2f, %.2f)\n", fix_idx, v.x(), v.y(), v.z());
		}
	}
	else if (tetMesh) {
		for (auto& fix_idx : fixed_indices) {
			auto& v = tetMesh->vertices(MeshKernel::iGameVertexHandle(fix_idx));
			printf("vex %d fixed at: (%.2f, %.2f, %.2f)\n", fix_idx, v.x(), v.y(), v.z());
		}
	}

}

std::string ModelDraw::add_moved_vertex(double x, double y, double z) {
	//assert(selected_v_indices.size() == 1);
	int idx = selected_v_indices.back();
	if (!surface_mesh) {
		std::cout << "surface mesh is nullptr!\n";
		return "surface mesh is nullptr!";
	}
	selected_v_indices.clear();
	auto v_ori = surface_mesh->vertices(MeshKernel::iGameVertexHandle(idx));
	iGameVertex v_tar(x, y, z);
	moved_vertices.emplace_back(MeshKernel::iGameVertexHandle(idx), v_tar);
	string str = "vex " + std::to_string(idx) + ": (" + to_string(v_ori.x()) + ", " + to_string(v_ori.y()) + ", " + to_string(v_ori.z()) +
		") move to (" + to_string(v_tar.x()) + ", " + to_string(v_tar.y()) + ", " + to_string(v_tar.z()) + ")";
	std::cout << str << std::endl;
	return str;
}

bool ModelDraw::BuildSkeletonToFile(std::string filename, uint32_t width, uint32_t height) {
	ModelClear();
	int dot_idx = filename.find_last_of('.');
	string fileformat = filename.substr(++dot_idx, filename.length());
	//std::cout << "file format: " << fileformat << std::endl;
	std::cout << "filename: " << filename << std::endl;
	std::cout << "width = " << width << ", height = " << height << std::endl;

	MeshKernel::IO io;
	assert(fileformat == "obj");
	time_t t1 = clock();
	//poly_mesh = pipeline->BuildSkeletonToFile(filename, 722, 242);
	pipeline = new Pipeline();
	poly_mesh = new polyMesh();
	*poly_mesh = pipeline->BuildSkeletonToFile(filename, width, height);
	std::cout << "Build Skeleton Cost: " << int(clock() - t1) << " ms" << std::endl;

	SetDrawModule(SPARSEOPTIMIZATION);
	updateMeshModel();
	updateDrawData_PolyMesh_OPT_1();
	update();
	return true;

}

void ModelDraw::SparseOptimization(double max_error) {

	time_t t1 = clock();
	*poly_mesh = pipeline->SparseOptimization(1, max_error);
	std::cout << "Sparse Optimization Cost: " << int(clock() - t1) << " ms" << std::endl;
	updateDrawData_PolyMesh_OPT_2();
	update();

}

void ModelDraw::GenerateSymAxis() {
	time_t t1 = clock();
	pipeline->GenerateSymAxis(axis);
	std::cout << "Generate SymAxis Cost: " << int(clock() - t1) << " ms" << std::endl;
	update();
}

void ModelDraw::ShapeApproximation(uint32_t max_iter) {
	time_t t1 = clock();
	*poly_mesh = pipeline->ShapeApproximation(max_iter);

	pipeline->ExtractMMCBoundary(bdy);
	axis.clear();
	std::cout << "Shpae Approximation Cost: " << int(clock() - t1) << " ms" << std::endl;
	updateDrawData_PolyMesh_OPT_V();
	update();
}
void ModelDraw::Read_Bdy(std::string& filename) {
	IGA::ReadBdy(filename, hexbdy, IGA::protected_cell);
}

void ModelDraw::savePickItem(string fileName)
{
    std::cout <<"WriteBdy\n";
    std::cout <<IGA::protected_cell.size() <<"\n";
	IGA::WriteBdy(fileName, hexbdy, IGA::protected_cell);
}
bool ModelDraw::loadPickFile(const std::string& filename)
{
	hexbdy.clear();
	IGA::ReadBdy(filename, hexbdy, topopt_pro);

    // todo:
    auto find_close_item  = [&](int type, std::vector<std::vector<double>>& points)->int{
        int ans = -1;
        double min_err = std::numeric_limits<double>::max();
        iGameVertex pp;
        int cnt = 0;
        for(auto& vec_d : points){
            pp += (iGameVertex(vec_d[0], vec_d[1], vec_d[2]));
            cnt++;
        }
        auto p_mid = pp / cnt;

        if(type == 1){
            for(const auto& [vh, v] : hexMesh->allvertices()){
                double err = v.dist(p_mid);
                if(err < min_err) {
                    min_err = err;
                    ans = vh.idx();
                }
            }
        }
        if(type == 2){
            for(const auto& [eh, e] : hexMesh->alledges()){
                iGameVertex edge_mid_vertex = (hexMesh->vertices(e.vh1()) + hexMesh->vertices(e.vh2())) / 2;
                double err = edge_mid_vertex.dist(p_mid);
                if(err < min_err) {
                    min_err = err;
                    ans = eh.idx();
                }
            }
        }
        if(type == 3){
            for(const auto& [fh, f] : hexMesh->allfaces()){
                iGameVertex face_mid_vertex;
                for(int id = 0; id < 4; ++id){
                    face_mid_vertex += hexMesh->vertices(f.vh(id));
                }
                face_mid_vertex /= 4;
                double err = face_mid_vertex.dist(p_mid);
                if(err < min_err) {
                    min_err = err;
                    ans = fh.idx();
                }
            }
        }
        if(type == 4){
            for(const auto& [ch, c] : hexMesh->allcells()){
                iGameVertex cell_mid_vertex;
                for(int id = 0; id < 8; ++id){
                    cell_mid_vertex += hexMesh->vertices(c.vh(id));
                }
                cell_mid_vertex /= 4;
                double err = cell_mid_vertex.dist(p_mid);
                if(err < min_err) {
                    min_err = err;
                    ans = ch.idx();
                }
            }
        }
        assert(ans > 0);
        return ans;
    };
    for(auto& temp_bdy : hexbdy){
        if(temp_bdy.element_type == IGA::HexBdyBase::Point){
            auto idx = find_close_item(1, temp_bdy.element_points);
            if(temp_bdy.bdy_type == IGA::BdyBase::Neumann) force_v_idx.push_back(idx);
            else fix_v_idx.push_back(idx);
        }else if(temp_bdy.element_type == IGA::HexBdyBase::Line){
            auto idx = find_close_item(2, temp_bdy.element_points);
            if(temp_bdy.bdy_type == IGA::BdyBase::Neumann) force_e_idx.push_back(idx);
            else fix_e_idx.push_back(idx);
        }else if(temp_bdy.element_type == IGA::HexBdyBase::Face){
            auto idx = find_close_item(3, temp_bdy.element_points);
            if(temp_bdy.bdy_type == IGA::BdyBase::Neumann) force_f_idx.push_back(idx);
            else fix_f_idx.push_back(idx);
        }else if(temp_bdy.element_type == IGA::HexBdyBase::Cell){
            auto idx = find_close_item(4, temp_bdy.element_points);
            if(temp_bdy.bdy_type == IGA::BdyBase::Neumann) force_c_idx.push_back(idx);
            else fix_c_idx.push_back(idx);
        }
    }
	//cout << hexbdy.size() << "\n";
	//cout << IGA::protected_cell << "\n";
	return false;
}

//最后全部的边界信息的存储




// bdy_type: 0, fixed; 1, force.
void ModelDraw::addPickData(vector<double> bdy_value, int bdy_type)
{
	auto a = IGA::BdyBase::Dirichlet;
	if (bdy_type == 1) {
		//force
		a = IGA::BdyBase::Neumann;
	}
	//对此时的所以已经点选的索引加东西
	for (auto tempInt : selected_v_indices) {
        auto vh = iGameVertexHandle(tempInt);
		iGameVertex v = hexMesh->vertices(vh);
		IGA::HexBdyBase temp;
		temp.element_type = IGA::HexBdyBase::Point;
		temp.bdy_type = a;
		IGA::HexBdyBase::VertexPos tempPos{ v.x(),v.y(),v.z() };
		temp.element_points = { tempPos };
		temp.bdy_value = bdy_value;
		hexbdy.push_back(temp);
	}
	for (auto tempInt : selected_e_indices) {
		auto eh = iGameEdgeHandle(tempInt);
		auto& e = hexMesh->edges(eh);
		auto& v1 = hexMesh->vertices(e.vh1());
		auto& v2 = hexMesh->vertices(e.vh2());
		IGA::HexBdyBase::ElementPoints tempPos{ {v1.x(),v1.y(),v1.z()},{v2.x(),v2.y(),v2.z()} };
		IGA::HexBdyBase temp;
		temp.element_type = IGA::HexBdyBase::Line;
		temp.bdy_type = a;
		temp.element_points = tempPos;
		temp.bdy_value = bdy_value;
		hexbdy.push_back(temp);

	}
	for (auto tempInt : selected_f_indices) {
		auto fh = iGameFaceHandle(tempInt);
		auto& f = hexMesh->faces(fh);
		auto vec = f.getVertexHandle();
		IGA::HexBdyBase temp;
		IGA::HexBdyBase::ElementPoints tempPos;
		for (auto& i : vec)
		{
			auto& v = hexMesh->vertices(i);
			IGA::HexBdyBase::VertexPos tempV{ v.x(),v.y(),v.z() };
			tempPos.push_back(tempV);
		}
		temp.element_points = tempPos;
		temp.element_type = IGA::HexBdyBase::Face;
		temp.bdy_type = a;
		temp.bdy_value = bdy_value;
		hexbdy.push_back(temp);
	}
	for (auto tempInt : selected_c_indices)
	{
		auto ch = iGameCellHandle(tempInt);
		auto& c = hexMesh->cells(ch);
		auto vec = c.getVertexHandle();
		IGA::HexBdyBase temp;
		IGA::HexBdyBase::ElementPoints tempPos;
		for (auto& i : vec)
		{
			auto& v = hexMesh->vertices(i);
			IGA::HexBdyBase::VertexPos tempV{ v.x(),v.y(),v.z() };
			tempPos.push_back(tempV);
		}
		temp.element_points = tempPos;
		temp.element_type = IGA::HexBdyBase::Cell;
		temp.bdy_type = a;
		temp.bdy_value = bdy_value;
		hexbdy.push_back(temp);
	}
    auto append_std_vec = [](std::vector<int>& a, std::vector<int>& b){
        a.insert(a.end(), b.begin(), b.end());
    };
    if(bdy_type == 0){
        append_std_vec(fix_v_idx, selected_v_indices);
        append_std_vec(fix_e_idx, selected_e_indices);
        append_std_vec(fix_f_idx, selected_f_indices);
        append_std_vec(fix_c_idx, selected_c_indices);
    }else{
        append_std_vec(force_v_idx, selected_v_indices);
        append_std_vec(force_e_idx, selected_e_indices);
        append_std_vec(force_f_idx, selected_f_indices);
        append_std_vec(force_c_idx, selected_c_indices);
    }
	selected_v_indices.clear();
	selected_e_indices.clear();
	selected_f_indices.clear();
	selected_c_indices.clear();

}

void ModelDraw::TobeProtected()
{
	cout << "selected_c_indices    " << selected_c_indices << "\n";
	///读出当前点选体的索引，
	for (int id : selected_c_indices)
	{
		IGA::protected_cell.push_back(id);
	}
    cout <<IGA::protected_cell<<"\n";
	selected_c_indices.clear();
}

void ModelDraw::updatePickPara() {
	assert(hexMesh);
	pickedge.clear();
	pickface.clear();
	pickcell.clear();

	for (auto& e_idx : selected_e_indices) {
		MeshKernel::iGameEdgeHandle eh(e_idx);
		auto e = hexMesh->edges(eh);
		pickedge.push_back({ (int)e.vh1(), (int)e.vh2() });
	}

	for (auto& f_idx : selected_f_indices) {
		MeshKernel::iGameFaceHandle fh(f_idx);
		auto face = hexMesh->faces(fh);
		std::vector<int> _tmp;
		for (auto v_idx : face.getVertexHandle()) {
			_tmp.push_back(int(v_idx));
		}
		pickface.push_back(_tmp);
	}
	for (auto& c_idx : selected_c_indices) {
		MeshKernel::iGameCellHandle ch(c_idx);
		auto cell = hexMesh->cells(ch);
		std::vector<int> _tmp;
		for (auto v_idx : cell.getVertexHandle()) {
			_tmp.push_back(int(v_idx));
		}
		pickcell.push_back(_tmp);
	}

}

bool ModelDraw::updateAllHandles() {
	switch (drawmodule)
	{
	case ViewOpenGL::EMPTY:
		break;
	case ViewOpenGL::SURFACEMESH:
		if (this->surface_mesh) {
			this->surface_mesh->updateAllHandles();
			removeSelectedAll();
			updateDrawData();
			return true;
		}
		break;
	case ViewOpenGL::HEXMESH:
		if (this->hexMesh) {
			updateAllHandles_kepp_constraint(*hexMesh, force_id, spc_id);
			removeSelectedAll();
			updateDrawData();
			return true;
		}
		break;
	case ViewOpenGL::SKELMESH:
		if (this->skeletal_mesh) {
			this->skeletal_mesh->updateAllHandles();
			removeSelectedAll();
			updateDrawData();
			return true;
		}
	case ViewOpenGL::TETMESH:
		break;
	case ViewOpenGL::HYBRIDMESH:
		break;
	case ViewOpenGL::BSPLINECURVE:
		break;
	case ViewOpenGL::BSPLINESURFACE:
		break;
	case ViewOpenGL::SPARSEOPTIMIZATION:
		break;
	case ViewOpenGL::TRIMESH:
		break;
	default:
		break;
	}
	QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Support iGame-Model Only! "));
	return false;
}

void ModelDraw::deformation_arap() {

	if (drawmodule == ViewOpenGL::SURFACEMESH) {
		if (surface_mesh) {
			Deformation app(surface_mesh);
			//Deformation_ARAP app(*surface_mesh);
			app.Execute(fixed_indices, moved_vertices);
			updateDrawData_SurfaceMesh();
			update();
		}
	}

}

void ModelDraw::parameterization(int mode) {
	if (surface_mesh) {
		Parameterization app(*surface_mesh, mode);
		app.Execute();
		updateMeshModel();
		updateDrawData();
		update();
	}
}

void ModelDraw::simplification_qem(int target_vertices_size, bool powered_flag) {
	if (surface_mesh) {
		if (target_vertices_size < surface_mesh->VertexSize()) {
			if (powered_flag) {
				Simplification_Powered_QEM app(*surface_mesh, target_vertices_size);
				app.Execute();
			}
			else {
				Simplification_QEM app(*surface_mesh, target_vertices_size);
				app.Execute();
			}
			updateDrawData_SurfaceMesh();
			update();
		}
	}
}

void ModelDraw::subdivision_loop(int iter) {

	if (surface_mesh) {
		Subdivision_Loop app(*surface_mesh, iter);
		app.Execute();

		updateDrawData_SurfaceMesh();
		update();
	}

}

void ModelDraw::denoising_bnf(int iter) {
	if (surface_mesh) {
		for (int i = 0; i < iter; ++i) {
			Denoising_BNF app(*surface_mesh);
			app.Execute();
		}
		updateDrawData_SurfaceMesh();
		update();
	}
}

void ModelDraw::remeshing_uniform() {
	if (surface_mesh) {
		Remeshing app(*surface_mesh, true);
		app.Execute();
		updateDrawData_SurfaceMesh();
		update();
	}
}

void ModelDraw::remeshing_adaptive() {
	if (surface_mesh) {
		Remeshing app(*surface_mesh, false);
		app.Execute();
		updateDrawData_SurfaceMesh();
		update();
	}
}

std::string ModelDraw::convertMeshToString() {
	if (!this->surface_mesh) return "";
	MeshKernel::IO io;
	return io.WriteOffString(*surface_mesh);
}

bool ModelDraw::loadTetMeshFromStr(const std::string& data) {

	//clock_t time1 = clock();
	//delete this->tetMesh;
	//this->tetMesh = nullptr;
	//MeshKernel::IO io;
	//tetMesh = new MeshKernel::TetMesh(io.ReadMeshFileFromStr(data));


	//int sides_num;
	//puts("success transfer");

	////   updateMeshModel();
	//clock_t time2 = clock();
	//std::cout << time2 - time1 << " ms" << std::endl;
	//SetDrawModule(TETMESH);
	//updateMeshModel();
	//updateDrawData_TetMesh();
	//update();

	return true;
}

void ModelDraw::test_bspline() {

	vector<CBSplineSurface> patch;
	for (int i = 0; i < bsplineSurfaces.size(); ++i) {
		patch.push_back(bsplineSurfaces[i]);
	}

	loadIGESSurface(patch);

	std::cout << "load bspline surface from patch success." << std::endl;
}

void ModelDraw::loadSkeletalMesh(std::string filename) {

	ModelClear();

	skeletal_mesh = new SkeletalMesh();
	skeletal_mesh->readFile(filename);
	if (scale_to_uint_flag) skeletal_mesh->scale2Uint();
	SetDrawModule(SKELMESH);
	if (resetViewFlag) updateMeshModel();

	updateDrawData();
	update();

}

void ModelDraw::saveSkeletalMesh(std::string filename) {

	if (skeletal_mesh) {
		if (scale_to_uint_flag) skeletal_mesh->scale2Origin();
		skeletal_mesh->writeFile(filename);
	}

}

void ModelDraw::updateDrawData_SkelMesh() {
	vertices_data.clear();
	if (!skeletal_mesh) return;

	v_cnt = skeletal_mesh->vsize();
	e_cnt = skeletal_mesh->esize();
	tri_cnt = skeletal_mesh->fsize();
	vertices_data.reserve(6 * (2 * e_cnt + 3 * tri_cnt));


	for (auto ep : skeletal_mesh->alledges()) {
		auto& v1 = skeletal_mesh->vertices(ep.second.vh1());
		vertices_data.push_back(v1.x());
		vertices_data.push_back(v1.y());
		vertices_data.push_back(v1.z());
        vertices_data.push_back(0);
        vertices_data.push_back(0);
        vertices_data.push_back(0);
		auto& v2 = skeletal_mesh->vertices(ep.second.vh2());
		vertices_data.push_back(v2.x());
		vertices_data.push_back(v2.y());
		vertices_data.push_back(v2.z());
        vertices_data.push_back(0);
        vertices_data.push_back(0);
        vertices_data.push_back(0);
	}

	for (auto fp : skeletal_mesh->allfaces()) {
		for (auto& vh : fp.second.getVertexHandle()) {
			auto& v = skeletal_mesh->vertices(vh);
			vertices_data.push_back(v.x());
			vertices_data.push_back(v.y());
			vertices_data.push_back(v.z());
			vertices_data.push_back(0.7725);
			vertices_data.push_back(0.7725);
			vertices_data.push_back(0);
		}
	}

	if (!show_tangent_plane) return;
	quad_cnt = 0;
	for (auto& vp : skeletal_mesh->allvertices()) {
		auto& vex = vp.second;
		quad_cnt += 6;
		Eigen::Vector3d pos(vex.x(), vex.y(), vex.z());
		double radius = skeletal_mesh->radius[vp.first];
		vector<vector<int>> dirs = {
			{-1, -1, -1}, {-1, -1, 1}, {1,-1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, 1, 1}, {1,1, 1}, {1, 1, -1}
		};
		vector<Eigen::Vector3d> v_plane;
		for (auto& dir : dirs) {
			Eigen::Vector3d move(dir[0], dir[1], dir[2]);
			Eigen::Vector3d v_p = pos + radius * move;
			v_plane.push_back(v_p);
		}
		vector<vector<int>> quad_faces = {
			{0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 5, 6}
		};
		for (auto& face : quad_faces) {
			for (auto& vid : face) {
				auto& vpos = v_plane[vid];
				vertices_data.push_back(vpos.x());
				vertices_data.push_back(vpos.y());
				vertices_data.push_back(vpos.z());
				vertices_data.push_back(0.7725);
				vertices_data.push_back(0.7725);
				vertices_data.push_back(0.7725);
			}
		}

	}

//    vertices_data.clear();
//    if (!skeletal_mesh) return;
//
//    v_cnt = skeletal_mesh->vsize();
//    e_cnt = skeletal_mesh->esize();
//    tri_cnt = skeletal_mesh->fsize();
//    vertices_data.reserve(6 * (2 * e_cnt + 3 * tri_cnt));
//
//    for (auto ep : skeletal_mesh->alledges()) {
//        auto& v1 = skeletal_mesh->vertices(ep.second.vh1());
//        vertices_data.push_back(v1.x());
//        vertices_data.push_back(v1.y());
//        vertices_data.push_back(v1.z());
//        vertices_data.push_back(0);
//        vertices_data.push_back(0);
//        vertices_data.push_back(0);
//
//        auto& v2 = skeletal_mesh->vertices(ep.second.vh2());
//        vertices_data.push_back(v2.x());
//        vertices_data.push_back(v2.y());
//        vertices_data.push_back(v2.z());
//        vertices_data.push_back(0);
//        vertices_data.push_back(0);
//        vertices_data.push_back(0);
//    }
//
//    for (auto fp : skeletal_mesh->allfaces()) {
////        vector<iGameVertexHandle> face_vertices = fp.second.getVertexHandle();
////        Eigen::Vector3d normal(0, 0, 0);
////        if (face_vertices.size() == 3) {
////            auto& v1 = skeletal_mesh->vertices(face_vertices[0]);
////            auto& v2 = skeletal_mesh->vertices(face_vertices[1]);
////            auto& v3 = skeletal_mesh->vertices(face_vertices[2]);
////            Eigen::Vector3d vec1(v2.x() - v1.x(), v2.y() - v1.y(), v2.z() - v1.z());
////            Eigen::Vector3d vec2(v3.x() - v1.x(), v3.y() - v1.y(), v3.z() - v1.z());
////            normal = vec1.cross(vec2).normalized();
////        }
//
//        for (auto& vh : face_vertices) {
//            auto& v = skeletal_mesh->vertices(vh);
//            vertices_data.push_back(v.x());
//            vertices_data.push_back(v.y());
//            vertices_data.push_back(v.z());
//            vertices_data.push_back(normal.x());
//            vertices_data.push_back(normal.y());
//            vertices_data.push_back(normal.z());
//            vertices_data.push_back(0.7725);
//            vertices_data.push_back(0.7725);
//            vertices_data.push_back(0);
//        }
//    }
//
//    if (!show_tangent_plane) return;
//    quad_cnt = 0;
//    for (auto& vp : skeletal_mesh->allvertices()) {
//        auto& vex = vp.second;
//        quad_cnt += 6;
//        Eigen::Vector3d pos(vex.x(), vex.y(), vex.z());
//        double radius = skeletal_mesh->radius[vp.first];
//        vector<vector<int>> dirs = {
//                {-1, -1, -1}, {-1, -1, 1}, {1,-1, 1}, {1, -1, -1}, {-1, 1, -1}, {-1, 1, 1}, {1,1, 1}, {1, 1, -1}
//        };
//        vector<Eigen::Vector3d> v_plane;
//        for (auto& dir : dirs) {
//            Eigen::Vector3d move(dir[0], dir[1], dir[2]);
//            Eigen::Vector3d v_p = pos + radius * move;
//            v_plane.push_back(v_p);
//        }
//        vector<vector<int>> quad_faces = {
//                {0, 1, 2, 3}, {4, 5, 6, 7}, {0, 1, 5, 4}, {3, 2, 6, 7}, {0, 3, 7, 4}, {1, 2, 5, 6}
//        };
//        for (auto& face : quad_faces) {
//            for (auto& vid : face) {
//                auto& vpos = v_plane[vid];
//                vertices_data.push_back(vpos.x());
//                vertices_data.push_back(vpos.y());
//                vertices_data.push_back(vpos.z());
//                vertices_data.push_back(0.7725);
//                vertices_data.push_back(0.7725);
//                vertices_data.push_back(0.7725);
//            }
//        }
//    }


}

void ModelDraw::drawSkelMesh() {

	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
	shaderProgram.setUniformValue("colormode", -1);

	glEnable(GL_POLYGON_OFFSET_FILL);
	glPolygonOffset(1.0, 1.0);
	glLineWidth(meshline_width);
	core->glDrawArrays(GL_LINES, 0, e_cnt * 2);
	core->glDrawArrays(GL_TRIANGLES, e_cnt * 2, tri_cnt * 3);

	if (show_tangent_plane) {
		// 请在绘制其他物体后再绘制半透明物体
		glDepthRange(0.01, 2.0);
		shaderProgram.setUniformValue("colormode", 518);// 半透明
		core->glDrawArrays(GL_QUADS, e_cnt * 2 + tri_cnt * 3, quad_cnt * 4);
	}

}

void ModelDraw::updateBVH() {

	if (editmodel == EditModel::DEFORMATION) return;

	bvh_tree->destoryBVH_Tree();

	switch (drawmodule)
	{
	case ViewOpenGL::EMPTY:
		break;
	case ViewOpenGL::SURFACEMESH:
		bvh_tree->buildBVH_Tree(*surface_mesh);
		break;
	case ViewOpenGL::HEXMESH:
		bvh_tree->buildBVH_Tree(*hexMesh);
		break;
	case ViewOpenGL::TETMESH:
		bvh_tree->buildBVH_Tree(*tetMesh);
		break;
	case ViewOpenGL::HYBRIDMESH:
		break;
	case ViewOpenGL::BSPLINECURVE:
		break;
	case ViewOpenGL::BSPLINESURFACE:
		break;
	case ViewOpenGL::SKELMESH:
        bvh_tree->buildBVH_Tree(*skeletal_mesh);
		break;
	case ViewOpenGL::SKEL:
		break;
	case ViewOpenGL::SPARSEOPTIMIZATION:
		break;
	case ViewOpenGL::TRIMESH:
		break;
	default:
		break;
	}
	//cout << "build BVH success.\n";
}

int ModelDraw::delete_island_HexMesh() {
	int island_count = -1, volume_count = -1;
	if (!hexMesh) return island_count;
	unordered_map<int, bool> visited;
	vector<vector<int>> islands;
	for (auto& cp : hexMesh->allcells()) {
		int ch = cp.first;
		if (visited[ch]) continue;
		queue<int> que;
		que.push(ch);
		visited[ch] = true;
		vector<int> indices = { ch };
		while (!que.empty()) {
			MeshKernel::iGameCellHandle cur_ch(que.front());
			que.pop();
			for (auto& adjch : hexMesh->NeighborCh(cur_ch)) {
				if (visited.count(adjch)) continue;
				que.push(adjch);
				visited[adjch] = true;
				indices.push_back(adjch);
			}
		}
		islands.push_back(indices);
	}

	vector<pair<int, int>> idx2cnt;
	for (auto& island : islands) {
		idx2cnt.emplace_back(island.size(), idx2cnt.size());
	}
	sort(idx2cnt.begin(), idx2cnt.end());
	island_count = idx2cnt.size() - 1;
	for (int i = 0; i + 1 < idx2cnt.size(); ++i) {
		int idx = idx2cnt[i].second;
		volume_count += idx2cnt[i].first;
		for (auto& ch : islands[idx]) {
			MeshKernel::iGameCellHandle _ch(ch);
			hexMesh->DeleteCell(_ch);
		}
	}
	cout << "Delete isolated islands success. delete " << island_count << " islands and " << volume_count << " cells.\n";

	updateAllHandles();
	updateDrawData();
	return island_count;
}

int ModelDraw::delete_redundancy_cell_HexMesh() {

	if (!hexMesh) return -1;

	std::vector<int> redundancy_chs;
	for (auto& cp : hexMesh->allcells()) {
		if (hexMesh->NeighborCh(cp.first).size() < 3) {
			redundancy_chs.push_back(cp.first);
		}
	}
	for (auto& _ch : redundancy_chs) {
		hexMesh->DeleteCell(MeshKernel::iGameCellHandle(_ch));
	}

	cout << "Delete redundancy cells success. delete " << redundancy_chs.size() << " cells.\n";

	updateAllHandles();
	updateDrawData();
	return redundancy_chs.size();
}

void ModelDraw::read_constraint_force(string str) {
	read_constraint_data(str, force_id);
}

void ModelDraw::read_constraint_spc(string str) {
	read_constraint_data(str, spc_id);
}

void ModelDraw::write_constraint_force(string str) {
	write_constraint_data(str, force_id);
}

void ModelDraw::write_constraint_spc(string str) {
	write_constraint_data(str, spc_id);
}

void ModelDraw::remove_all_constraints() {
	force_id.clear();
	spc_id.clear();
}

int ModelDraw::getVerticesSize() {
	switch (drawmodule)
	{
	case ViewOpenGL::EMPTY:
		return 0;
		break;
	case ViewOpenGL::SURFACEMESH:
		return surface_mesh->VertexSize();
		break;
	case ViewOpenGL::HEXMESH:
		return hexMesh->vsize();
		break;
	case ViewOpenGL::TETMESH:
		return tetMesh->vsize();
		break;
	case ViewOpenGL::HYBRIDMESH:
		break;
	case ViewOpenGL::BSPLINECURVE:
		break;
	case ViewOpenGL::BSPLINESURFACE:
		break;
	case ViewOpenGL::SKELMESH:
		return skeletal_mesh->vsize();
		break;
	case ViewOpenGL::SKEL:
		break;
	case ViewOpenGL::SPARSEOPTIMIZATION:
		break;
	case ViewOpenGL::TRIMESH:
		return tri_mesh->n_vertices();
		break;
	default:
		break;
	}
	return -1;
}

void ModelDraw::simplification_qmat(int target_vsize) {

	// 若不放大可能会因为精度导致结果错误
	for (auto& vp : skeletal_mesh->allvertices()) {
		auto& v = skeletal_mesh->vertices(vp.first);
		v *= 1000;
		skeletal_mesh->radius[vp.first] *= 1000;
	}

	Q_MAT app(*skeletal_mesh);
	app.simplification(target_vsize);

	for (auto& vp : skeletal_mesh->allvertices()) {
		auto& v = skeletal_mesh->vertices(vp.first);
		v *= 0.001;
		skeletal_mesh->radius[vp.first] *= 0.001;
	}

	//updateMeshModel();
	updateDrawData();

}

void ModelDraw::subdivision_skl(int iter) {
	if (skeletal_mesh) {
		Subdivision_Skl app(skeletal_mesh);
		app.subdivision(iter);
	}
	updateDrawData();
}

void ModelDraw::compared_bind_vertices() {
	if (hexMesh == nullptr) return;

	compared_to_current.clear();
	int cmp_vcnt = compared_vertices_data.size() / 6;
	for (int i = 0; i < cmp_vcnt; ++i) {
		int min_vh = -1, idx = i * 6;
		double min_dis = 999999999;
		MeshKernel::iGameVertex cmp_vex(compared_vertices_data[idx], compared_vertices_data[idx + 1], compared_vertices_data[idx + 2]);
		for (auto& vp : hexMesh->allvertices()) {
			double dis = (cmp_vex - vp.second).norm();
			if (dis < min_dis) {
				min_dis = dis;
				min_vh = vp.first;
				if (min_dis < 1E-6F) break;
			}
		}
		compared_to_current[i] = min_vh;
	}

}

double ModelDraw::getVolume() {
	double res = 0;

	if (surface_mesh) res = MeshMath::get_volume_surface_mesh(*surface_mesh);
	else if (hexMesh) res = MeshMath::get_volume_hexahedral_mesh(*hexMesh);

	return res;
}

void ModelDraw::test_action(double a) {

	if (surface_mesh) {
		MeshKernel::SurfaceMesh quadmesh;
		QuadMesher app(*surface_mesh, {});
		app.meshing(quadmesh);
		app.get_streamlines(streamlines, singularity_vertices);
	}


	updateDrawData();

}

void ModelDraw::enable_realtime_deformation() {

	if (surface_mesh == nullptr && tetMesh == nullptr) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Mesh is empty! "));
		return;
	}
	delete tri_deformation;
	delete tet_deformation;
	if (surface_mesh) tri_deformation = new Deformation(surface_mesh);
	else if (tetMesh) tet_deformation = new TetDeformation(tetMesh);
	if (tri_deformation == nullptr) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" Deformation init fail! "));
		return;
	}
	SetEditModel(EditModel::DEFORMATION);

}

void ModelDraw::calc_cross_field() {
	if (!surface_mesh || !surface_mesh->isTriangleMesh()) return;
	std::vector<Vec3> points;
	std::vector<ID3> triangles;
	std::vector<ID2> feature_edges;
	points.reserve(surface_mesh->vsize());
	triangles.reserve(surface_mesh->fsize());

	for (int i = 0; i < surface_mesh->vsize(); ++i) {
		VH vh(i);
		auto& v = surface_mesh->vertices(vh);
		points.push_back(Vec3{ v.x(), v.y(), v.z() });
	}

	for (int i = 0; i < surface_mesh->fsize(); ++i) {
		FH fh(i);
		auto& f = surface_mesh->faces(fh);
		const auto& vhs = f.getVertexHandle();
		triangles.push_back(ID3{ (ID)vhs[0], (ID)vhs[1], (ID)vhs[2] });
	}

	std::cout << "Initialize mesh success. #V = " << points.size() <<
		", #E = " << feature_edges.size() << ", #F = " << triangles.size() << std::endl;

	std::vector<ID3> singularities;// 奇异点的 id
	std::vector<std::array<double, 9>> global_triangle_dir;

	MeshMath::calc_cross_field(points, triangles, feature_edges, singularities, global_triangle_dir);

	int sv_cnt = 0, se_cnt = 0, sf_cnt = 0;

	for (auto& singular_id : singularities) {// 所有的奇异点
		if (singular_id[0] != NO_ID) {
			selected_v_indices.push_back((int)singular_id[0]);
			sv_cnt++;
		}
		else if (singular_id[1] != NO_ID) {
			int fidx = singular_id[1] / 3;
			int eidx = singular_id[1] % 3;
			FH fh(fidx);
			selected_e_indices.push_back(surface_mesh->faces(fh).eh(1));
			se_cnt++;
		}
		else if (singular_id[2] != NO_ID) {
			selected_f_indices.push_back((int)singular_id[2]);
			sf_cnt++;
		}
	}

	/*for (int fid = 0; fid < triangles.size(); ++fid) {// show result
		std::cout << "#T " << fid << ": ";
		for (int vid = 0; vid < 3; ++vid) {
			std::cout << "(" << global_triangle_dir[fid][vid * 3] << ", "
				<< global_triangle_dir[fid][vid * 3 + 1] << ", "
				<< global_triangle_dir[fid][vid * 3 + 2] << "), ";
		}
		std::cout << std::endl;
	}*/

	std::cout << "Singular Vertices: vcnt = " << sv_cnt << ", ecnt = " << se_cnt << ", fcnt = " << sf_cnt << std::endl;

	updateDrawData();

}



void ModelDraw::getBsplineSolidFromParam(int& udegree, int& vdegree, int& wdegree, vector<vector<vector<CPoint>>>& controlPoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis) {
    Poster_Needed::CBSplineSolid* bsplines = new Poster_Needed::CBSplineSolid(udegree, vdegree, wdegree, controlPoints, uknots, vknots, wknots, umultis, vmultis, wmultis);
	bsplinesolid = *bsplines;
	bsplinesolid.sethandle(bsplinesolidmesh.size());
	bsplinesolidmesh.push_back(bsplinesolid);
	auto bp = bsplinesolidmesh.back();
	SaveIndex(bp.getControlPoints(), bp.gethandle());

}

int ModelDraw::getBsplineSolidFromParamXJ(int& udegree, int& vdegree, int& wdegree, vector<vector<vector<CPoint>>>& controlPoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis) {
    Poster_Needed::CBSplineSolid* bsplines = new Poster_Needed::CBSplineSolid(udegree, vdegree, wdegree, controlPoints, uknots, vknots, wknots, umultis, vmultis, wmultis);
	bsplinesolidmesh.push_back(*bsplines);
	auto& bp = bsplinesolidmesh.back();
	bp.sethandle(bsplinesolidmesh.size() - 1);//设置solidhandle为其下标

	bp.m_UKnots_nonredundant.clear();
	for (int i = 0; i < uknots.size(); i++)
	{
		if (i + 1 < uknots.size() && uknots[i] != uknots[i + 1])
			bp.m_UKnots_nonredundant.push_back(uknots[i]);
		else if (i == uknots.size() - 1)
			bp.m_UKnots_nonredundant.push_back(uknots[i]);
	}
	bp.m_VKnots_nonredundant.clear();
	for (int i = 0; i < vknots.size(); i++)
	{
		if (i + 1 < vknots.size() && vknots[i] != vknots[i + 1])
			bp.m_VKnots_nonredundant.push_back(vknots[i]);
		else if (i == vknots.size() - 1)
			bp.m_VKnots_nonredundant.push_back(vknots[i]);
	}
	bp.m_WKnots_nonredundant.clear();
	for (int i = 0; i < wknots.size(); i++)
	{
		if (i + 1 < wknots.size() && wknots[i] != wknots[i + 1])
			bp.m_WKnots_nonredundant.push_back(wknots[i]);
		else if (i == wknots.size() - 1)
			bp.m_WKnots_nonredundant.push_back(wknots[i]);
	}
	return bp.gethandle();//
}
void ModelDraw::SaveIndex(vector<vector<vector<CPoint>>> controlmesh, int handle) {
	for (int i = 0; i < controlmesh.size(); i++) {
		for (int j = 0; j < controlmesh[i].size(); j++) {
			for (int k = 0; k < controlmesh[i][j].size(); k++) {
				AllControlPoint.insert(controlmesh[i][j][k]);//移动顶点后，删除元顶点，加入新顶点
				vector<int> tmp;
				if (PointBelongToSolid[controlmesh[i][j][k]].empty()) {
					tmp.push_back(handle);
					PointBelongToSolid[controlmesh[i][j][k]] = tmp;
				}
				else {
					tmp = PointBelongToSolid[controlmesh[i][j][k]];
					tmp.push_back(handle);
					PointBelongToSolid[controlmesh[i][j][k]] = tmp;
				}
			}
		}
	}
}
void ModelDraw::trans11(std::vector<double>& modelMatrix_, QMatrix4x4& modelmatrix)
{
	modelMatrix_[0] = modelmatrix(0, 0);
	modelMatrix_[1] = modelmatrix(1, 0);
	modelMatrix_[2] = modelmatrix(2, 0);
	modelMatrix_[3] = modelmatrix(3, 0);
	modelMatrix_[4] = modelmatrix(0, 1);
	modelMatrix_[5] = modelmatrix(1, 1);
	modelMatrix_[6] = modelmatrix(2, 1);
	modelMatrix_[7] = modelmatrix(3, 1);
	modelMatrix_[8] = modelmatrix(0, 2);
	modelMatrix_[9] = modelmatrix(1, 2);
	modelMatrix_[10] = modelmatrix(2, 2);
	modelMatrix_[11] = modelmatrix(3, 2);
	modelMatrix_[12] = modelmatrix(0, 3);
	modelMatrix_[13] = modelmatrix(1, 3);
	modelMatrix_[14] = modelmatrix(2, 3);
	modelMatrix_[15] = modelmatrix(3, 3);
}



void ModelDraw::bsplinesolid_drawControlMesh(void) {//..............changed

	glDepthRange(0.01, 1.0);
	glDepthFunc(GL_LEQUAL);
	//glColor4f(1.0, 0.0, 0.0, 1.0);

	shaderProgram.setUniformValue("colormode", -1);
	for (auto& bp : bsplinesolidmesh) {
		controlMesh = bp.getControlPoints();
		if (bp.gethandle() == 0) {
			double sum = 0;
			for (int i = 0; i < controlMesh[0][0].size() - 1; i++)
			{
				sum += (controlMesh[0][0][i + 1] - controlMesh[0][0][i]).norm();
			}
			ControlPointSize = sum / (controlMesh[0][0].size() - 1) / 10;//选择一个较为合适的顶点大小
			ControlLineWidth = 50 * ControlPointSize;
		}
		for (int i = 0; i < controlMesh.size(); i++) {
			for (int j = 0; j < controlMesh[i].size(); j++) {
				for (int k = 0; k < controlMesh[i][j].size(); k++) {
					if (Is_Highlight.count(controlMesh[i][j][k]) && Is_Highlight[controlMesh[i][j][k]] == 1) {

						shaderProgram.setUniformValue("colormode", 232);
						drawControlPoint(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ(), ControlPointSize, 16, 16);
						shaderProgram.setUniformValue("colormode", -1);
					}
					else
					{
						shaderProgram.setUniformValue("colormode", 234);
						//glColor4f(0.0, 1.0, 0.0, 1.0);
						drawControlPoint(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ(), ControlPointSize, 16, 16);
					}


				}
			}
		}
		glLineWidth(ControlLineWidth);
		glDepthRange(0.01, 1.0);
		glDepthFunc(GL_LEQUAL);
		shaderProgram.setUniformValue("colormode", 530);
		//glColor4f(1.0, 0.0, 0.0, 1.0);
		glPolygonMode(GL_FRONT_AND_BACK, GL_LINE);
		glBegin(GL_LINES);
		for (int i = 0; i < controlMesh.size(); i++) {
			for (int j = 0; j < controlMesh[i].size(); j++) {
				for (int k = 0; k < controlMesh[i][j].size(); k++) {
					if (k + 1 < controlMesh[i][j].size()) {
						glVertex3d(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ());
						glVertex3d(controlMesh[i][j][k + 1].getX(), controlMesh[i][j][k + 1].getY(), controlMesh[i][j][k + 1].getZ());
					}
					if (j + 1 < controlMesh[i].size()) {
						glVertex3d(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ());
						glVertex3d(controlMesh[i][j + 1][k].getX(), controlMesh[i][j + 1][k].getY(), controlMesh[i][j + 1][k].getZ());
					}
					if (i + 1 < controlMesh.size()) {
						glVertex3d(controlMesh[i][j][k].getX(), controlMesh[i][j][k].getY(), controlMesh[i][j][k].getZ());
						glVertex3d(controlMesh[i + 1][j][k].getX(), controlMesh[i + 1][j][k].getY(), controlMesh[i + 1][j][k].getZ());
					}
				}
			}
		}
		glEnd();
		//编辑模式下，使用填充样式
		glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
	}
}

void ModelDraw::moveSolidControlPoint(double winX, double winY) {
	//将新的位置点point设置为选中控制点的新坐标
	double p[3];
	pick_point_with_depth(winX, winY, p, pickedItemPosition);
	vector<vector<vector<CPoint>>> tmp = bsplinesolid.getControlPoints();
	tmp[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setX(p[0]);
	tmp[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setY(p[1]);
	tmp[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]].setZ(p[2]);
	Is_Highlight[tmp[surfacepointIndex[0]][surfacepointIndex[1]][surfacepointIndex[2]]] = 1;//changed

	this->controlMesh = tmp;
	bsplinesolid.setControlPoints(tmp);
	bsplineSolidToDiscreteXJ();
	bsplinesolid_drawControlMesh();
	update();
}



//void ModelDraw::pickPointDepth(int x, int y, double p[3], double& d) {// 用深度信息来 pick 点
//	GLint viewport[4];
//	glGetIntegerv(GL_VIEWPORT, viewport);
//	GLdouble winX = double(x);
//	GLdouble winY = double(viewport[3] - y);
//	std::vector<double>  modelMatrix_(16, 0.0);
//	std::vector<double>  projectionmatrix_(16, 0.0);
//	auto xx = modelmatrix * viewmatrix;
//	trans11(modelMatrix_, xx);
//	trans11(projectionmatrix_, projectionmatrix);
//	gluUnProject(winX, winY, d, &modelMatrix_[0], &projectionmatrix_[0], viewport, &p[0], &p[1], &p[2]);
//}
//
//
//void ModelDraw::getPointDepth(int x, int y, double point[3], double& dz) {
//	GLint viewport[4];
//	glGetIntegerv(GL_VIEWPORT, viewport);
//	GLdouble winX = double(x);
//	GLdouble winY = double(viewport[3] - y);
//	GLfloat winZ = 0.0;
//	glReadPixels((int)winX, (int)winY, 1, 1, GL_DEPTH_COMPONENT, GL_FLOAT, &winZ);
//	dz = (GLdouble)winZ;
//	std::vector<double> modelMatrix_(16, 0.0);
//	std::vector<double> projectionmatrix_(16, 0.0);
//	trans11(modelMatrix_, viewmatrix);
//	trans11(projectionmatrix_, projectionmatrix);
//	std::vector<double> model = modelMatrix_;
//	gluUnProject(winX, winY, (GLdouble)winZ, &modelMatrix_[0], &projectionmatrix_[0], viewport, &point[0], &point[1], &point[2]);
//
//}


void ModelDraw::aotu_createBSS(int ucnum, int vcnum, int wcnum, int udegree, int vdegree, int wdegree, vector<vector<vector<CPoint>>>& controlpoints, vector<vector<CPoint>>& singlepatch, vector<CPoint>& singlepoints, vector<double>& uknots, vector<double>& vknots, vector<double>& wknots, vector<int>& umultis, vector<int>& vmultis, vector<int>& wmultis)
{
	CPoint du(0.0, 0.0, 0.0);//初始控制点坐标
	for (int i = 0; i < ucnum; i++) {
		CPoint tmp1 = du;
		for (int j = 0; j < vcnum; j++) {
			CPoint tmp2 = du;
			for (int k = 0; k < wcnum; k++) {
				singlepoints.push_back(du);
				du = CPoint(du.getX() + 0.5, du.getY(), du.getZ());
			}
			singlepatch.push_back(singlepoints);
			singlepoints.clear();
			du = CPoint(tmp2.getX(), tmp2.getY() + 0.5, tmp2.getZ());
		}
		controlpoints.push_back(singlepatch);
		singlepatch.clear();
		du = CPoint(tmp1.getX(), tmp1.getY(), tmp1.getZ() + 0.5);
	}
	//不含重复度的节点个数
	int uknum = ucnum + udegree + 1 - (2 * udegree);
	int vknum = vcnum + vdegree + 1 - (2 * vdegree);
	int wknum = wcnum + wdegree + 1 - (2 * wdegree);
	//设置u方向均匀节点（范围都是[0,1])
	double uknot = 0.0;
	uknots.push_back(0.0);
	for (int i = 0; i < uknum - 2; i++) {
		uknot = uknot + 1.0 / (uknum - 1);
		uknots.push_back(uknot);
	}
	uknots.push_back(1.0);

	//设置v方向均匀节点
	double vknot = 0.0;
	vknots.push_back(0.0);
	for (int i = 0; i < vknum - 2; i++) {
		vknot = vknot + 1.0 / (vknum - 1);
		vknots.push_back(vknot);
	}
	vknots.push_back(1.0);

	//设置w方向均匀节点
	double wknot = 0.0;
	wknots.push_back(0.0);
	for (int i = 0; i < wknum - 2; i++) {
		wknot = wknot + 1.0 / (wknum - 1);
		wknots.push_back(wknot);
	}
	wknots.push_back(1.0);

	//设置重复度,使曲面为闭合曲面（端点处控制点在曲面上）
	umultis.push_back(udegree + 1);
	for (int i = 0; i < uknum - 2; i++) {
		umultis.push_back(1);
	}
	umultis.push_back(udegree + 1);

	vmultis.push_back(vdegree + 1);
	for (int i = 0; i < vknum - 2; i++) {
		vmultis.push_back(1);
	}
	vmultis.push_back(vdegree + 1);

	wmultis.push_back(wdegree + 1);
	for (int i = 0; i < wknum - 2; i++) {
		wmultis.push_back(1);
	}
	wmultis.push_back(wdegree + 1);
}

//void ModelDraw::bsplineSolidToDiscrete() {
//	time_t beg = clock();
//	data_surfaces.clear();
//	controlMesh.clear();
//	double umin, umax, vmin, vmax, wmin, wmax;
//	vector<CPoint> singlePatch;// 每个面片的离散点
//	vector<vector<CPoint>> singlelSurface;// 每个曲面的离散点
//	vector<vector<vector<CPoint>>> allSurface;// 所有曲面集合的离散点
//
//	//vector<CVector>singlePatch_normal;
//	//vector<vector<CVector>>singlelSurface_normal;
//	//vector<vector<vector<CVector>>>allSurface_normal;
//	//CVector normal;//保存法向量
//
//	int totalVerticesSize = 0;
//	bsplinesolid.getUVMinMax(umin, umax, vmin, vmax, wmin, wmax);
//	int unumber, vnumber, wnumber;
//	unumber = upatchnum * (bsplinesolid.getSizeUControlPoints());
//	vnumber = vpatchnum * (bsplinesolid.getSizeVControlPoints());
//	wnumber = wpatchnum * (bsplinesolid.getSizeWControlPoints());
//	unumber = unumber < 100 ? unumber : 100;
//	vnumber = vnumber < 100 ? vnumber : 100;
//	wnumber = wnumber < 100 ? wnumber : 100;
//	//std::cout << "u_num: " << unumber << "; v_num: " << vnumber << std::endl;
//	//unumber = 5, vnumber = 5;// for faster rendering
//	unumber *= 0.2, vnumber *= 0.2f; wnumber *= 0.2f;
//	double sizeu = (umax - umin) / unumber - 1e-9;// 步距
//	double sizev = (vmax - vmin) / vnumber - 1e-9;
//	double sizew = (wmax - wmin) / wnumber - 1e-9;
//	vertices_data.clear();//new
//	CPoint tmp;
//	for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u = u + sizeu) {
//		for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v = v + sizev) {
//			for (double w = wmin + 1e-9; w < wmax - sizew - 1e-9; w = w + sizew) {
//
//				tmp = bsplinesolid.getPointAtParam(u, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				//顶点坐标
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				//颜色
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//				tmp = bsplinesolid.getPointAtParam(u, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v + sizev, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//				tmp = bsplinesolid.getPointAtParam(u, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u + sizeu, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				tmp = bsplinesolid.getPointAtParam(u, v, w + sizew);
//				vertices_data.push_back(tmp.getX());
//				vertices_data.push_back(tmp.getY());
//				vertices_data.push_back(tmp.getZ());
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//				vertices_data.push_back(0.7725);
//
//			}
//		}
//
//	}
//	totalVerticesSize = vertices_data.size() / 6;
//	vector<vector<vector<CPoint>>> singleSolidControlPoints;
//	singleSolidControlPoints = bsplinesolid.getControlPoints();
//	this->controlMesh = singleSolidControlPoints;
//	time_t end = clock();
//	printf("bsplineSurfaceToDiscrete cost: %zdms, vertices size: %zd, \n", int(end - beg), totalVerticesSize);
//
//}


void ModelDraw::bsplineSolidSimulation()
{
	if (bsplinesolidmesh.size() == 0) {
		std::cout << "please input BSplineSolid\n" << "\n";
		return;
	}
	double timeSpan = 1;
	double timeStep = 0.05;
	double density = 33.2;
	double gravitationalAcc = 9.8;
	std::vector<boundary_condition> m_dirichletSides, m_neumannSides;
	CVector traction(0, 110, 0);
	CVector traction1(0, -110, 0);
	m_neumannSides.push_back(boundary_condition(2, boundary::front, condition_type::neumann, traction1));
	m_neumannSides.push_back(boundary_condition(8, boundary::front, condition_type::neumann, traction1));
	m_neumannSides.push_back(boundary_condition(11, boundary::front, condition_type::neumann, traction));
	m_neumannSides.push_back(boundary_condition(5, boundary::front, condition_type::neumann, traction));

	m_dirichletSides.push_back(boundary_condition(11, boundary::north, condition_type::dirichlet, 0, 0));
	m_dirichletSides.push_back(boundary_condition(11, boundary::north, condition_type::dirichlet, 0, 2));
	m_dirichletSides.push_back(boundary_condition(5, boundary::north, condition_type::dirichlet, 0, 0));
	m_dirichletSides.push_back(boundary_condition(5, boundary::north, condition_type::dirichlet, 0, 2));
	m_dirichletSides.push_back(boundary_condition(2, boundary::north, condition_type::dirichlet, 0, 0));
	m_dirichletSides.push_back(boundary_condition(2, boundary::north, condition_type::dirichlet, 0, 2));
	m_dirichletSides.push_back(boundary_condition(8, boundary::north, condition_type::dirichlet, 0, 0));
	m_dirichletSides.push_back(boundary_condition(8, boundary::north, condition_type::dirichlet, 0, 2));
	stiffAssembler stiff(bsplinesolidmesh, m_dirichletSides, m_neumannSides);

	massAssembler mass(bsplinesolidmesh, density, stiff.system());


	timeIntegrator timeAssemble(stiff, mass);
	timeAssemble.setInitialDisplacement(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
	timeAssemble.setInitialVelocity(Eigen::MatrixXd::Zero(stiff.numDofs(), 1));
	timeAssemble.initialize();
	Eigen::MatrixXd displacement;
	for (int t = 0; t < (int)timeSpan / timeStep; ++t) {
		std::cout <<"step:" << t << "\n";
		timeAssemble.makeTimeStep(timeStep);
		stiff.constructCoefs(timeAssemble.displacementVector(), displacement);

		auto a = bsplinesolidmesh[8].coefs();
		int begin = 0;
		for (size_t i = 0; i < bsplinesolidmesh.size(); ++i) {
			int num = bsplinesolidmesh[i].getSizeControlPoints();
			bsplinesolidmesh[i].updateCoefs(displacement.block(begin, 0, num, 3));
			begin += num;
		}

		bsplineSolidToDiscreteXJ();
		bsplinesolid_drawControlMesh();
		update();
		 
		auto b = bsplinesolidmesh[8].coefs();

		displacement *= -1;
		begin = 0;
		for (size_t i = 0; i < bsplinesolidmesh.size(); ++i) {
			int num = bsplinesolidmesh[i].getSizeControlPoints();
			bsplinesolidmesh[i].updateCoefs(displacement.block(begin, 0, num, 3));
			begin += num;
		}
	}
}

void ModelDraw::bsplineSolidToDiscreteXJ() {
	time_t beg = clock();


	double umin, umax, vmin, vmax, wmin, wmax;
	int totalVerticesSize = 0;
	vertices_data.clear();//new
	singlepatch_vertics_number.clear();//new

	for (auto& bd : bsplinesolidmesh) {
		//std::cout << "1time" << std::endl;
		bd.getUVMinMax(umin, umax, vmin, vmax, wmin, wmax);
		int unumber, vnumber, wnumber;
		unumber = upatchnum * (bd.getSizeUControlPoints());
		vnumber = vpatchnum * (bd.getSizeVControlPoints());
		wnumber = wpatchnum * (bd.getSizeWControlPoints());
		unumber = unumber < 100 ? unumber : 100;
		vnumber = vnumber < 100 ? vnumber : 100;
		wnumber = wnumber < 100 ? wnumber : 100;
		//std::cout << "u_num: " << unumber << "; v_num: " << vnumber << std::endl;
		//unumber = 5, vnumber = 5;// for faster rendering
		unumber *= 0.2, vnumber *= 0.2f; wnumber *= 0.2f;

		double sizeu = (umax - umin) / unumber - 1e-9;// 步距
		double sizev = (vmax - vmin) / vnumber - 1e-9;
		double sizew = (wmax - wmin) / wnumber - 1e-9;

		//std::cout << "size=" << vertices_data.size() << std::endl;
		CPoint tmp;
		int count = 0;
		if (isShowinside)
		{

			for (double u = umin + 1e-9; u <= umax - sizeu - 1e-9; u = u + sizeu) {
				for (double v = vmin + 1e-9; v < vmax - sizev - 1e-9; v = v + sizev) {
					for (double w = wmin + 1e-9; w < wmax - sizew - 1e-9; w = w + sizew) {

						tmp = bd.getPointAtParam(u, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
						//顶点坐标
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        //法向量信息
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						//颜色
						tmp = bd.getPointAtParam(u + sizeu, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

						tmp = bd.getPointAtParam(u, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

						tmp = bd.getPointAtParam(u, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

						tmp = bd.getPointAtParam(u, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

						tmp = bd.getPointAtParam(u, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v, w + sizew);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);

					}
				}

			}
		}

		else
		{
			bd.U0_size = 0;
			bd.U1_size = 0;
			bd.V0_size = 0;
			bd.V1_size = 0;
			bd.W0_size = 0;
			bd.W1_size = 0;
			for (double u = umin; u <= umax; u = u + umax - 1e-9) {//暂时只在U0U1加入法向，看看对不对
				vector<CPoint>* tmpvector;
				vector<vector<CPoint>>* normal_tmpvector;
				int* num;
				if (u == umin) {
					tmpvector = &bd.U_0_discrete_point;
					normal_tmpvector = &bd.U_0_point_to_calcuNormal;
					num = &bd.U0_size;
				}
				else
				{
					tmpvector = &bd.U_1_discrete_point;
					normal_tmpvector = &bd.U_1_point_to_calcuNormal;
					num = &bd.U1_size;
				}
				vector<CPoint> normal_singleLine;
				for (double v = vmin; v <= vmax - sizev; v = v + sizev) {
					for (double w = wmin; w <= wmax - sizew; w = w + sizew) {
						tmp = bd.getPointAtParam(u, v, w);
						normal_singleLine.push_back(tmp);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
						//顶点坐标
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        //法向量信息
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						//颜色
						tmp = bd.getPointAtParam(u, v + sizev, w);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w + sizew);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v, w + sizew);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						count += 24;
						*num += 24;
					}
					normal_singleLine.push_back(tmp);
					normal_tmpvector->push_back(normal_singleLine);
					normal_singleLine.clear();
				}
				for (double w = wmin; w <= wmax - sizew; w = w + sizew) {
					tmp = bd.getPointAtParam(u, vmax, w);
					normal_singleLine.push_back(tmp);
				}
				tmp = bd.getPointAtParam(u, vmax, wmax);
				normal_singleLine.push_back(tmp);
				normal_tmpvector->push_back(normal_singleLine);
			}

			for (double v = vmin; v <= vmax; v = v + vmax - 1e-9) {
				vector<CPoint>* tmpvector;
				vector<vector<CPoint>>* normal_tmpvector;
				int* num;
				if (v == vmin) {
					tmpvector = &bd.V_0_discrete_point;
					normal_tmpvector = &bd.V_0_point_to_calcuNormal;
					num = &bd.V0_size;
				}
				else
				{
					tmpvector = &bd.V_1_discrete_point;
					normal_tmpvector = &bd.V_1_point_to_calcuNormal;
					num = &bd.V1_size;
				}
				vector<CPoint> normal_singleLine;
				for (double u = umin; u <= umax - sizeu; u = u + sizeu) {

					for (double w = wmin; w <= wmax - sizew; w = w + sizew) {
						count += 24;
						*num += 24;
						tmp = bd.getPointAtParam(u, v, w);
						normal_singleLine.push_back(tmp);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
						//顶点坐标
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        //法向量信息
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						//颜色
						tmp = bd.getPointAtParam(u + sizeu, v, w);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v, w + sizew);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v, w + sizew);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
					}
					normal_singleLine.push_back(tmp);
					normal_tmpvector->push_back(normal_singleLine);
					normal_singleLine.clear();
				}
				for (double w = wmin; w <= wmax - sizew; w = w + sizew) {
					tmp = bd.getPointAtParam(umax, v, w);
					normal_singleLine.push_back(tmp);
				}
				tmp = bd.getPointAtParam(umax, v, wmax);
				normal_singleLine.push_back(tmp);
				normal_tmpvector->push_back(normal_singleLine);


			}

			for (double w = wmin; w <= wmax; w = w + wmax - 1e-9) {//
				vector<CPoint>* tmpvector;
				vector<vector<CPoint>>* normal_tmpvector;
				int* num;
				if (w == wmin) {
					tmpvector = &bd.W_0_discrete_point;
					normal_tmpvector = &bd.W_0_point_to_calcuNormal;
					num = &bd.W0_size;
				}
				else
				{
					tmpvector = &bd.W_1_discrete_point;
					normal_tmpvector = &bd.W_1_point_to_calcuNormal;
					num = &bd.W1_size;
				}
				vector<CPoint> normal_singleLine;
				for (double u = umin; u <= umax - sizeu; u = u + sizeu) {

					for (double v = vmin; v <= vmax - sizev; v = v + sizev) {//
						count += 24;
						*num += 24;
						tmp = bd.getPointAtParam(u, v, w);
						normal_singleLine.push_back(tmp);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
						//顶点坐标
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        //法向量信息
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						//颜色
						tmp = bd.getPointAtParam(u + sizeu, v, w);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u + sizeu, v + sizev, w);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						tmp = bd.getPointAtParam(u, v + sizev, w);
						tmpvector->push_back(tmp);
						vertices_data.push_back(tmp.getX());
						vertices_data.push_back(tmp.getY());
						vertices_data.push_back(tmp.getZ());
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
//                        vertices_data.push_back(0);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
						vertices_data.push_back(0.7725);
					}
					normal_singleLine.push_back(tmp);
					normal_tmpvector->push_back(normal_singleLine);
					normal_singleLine.clear();
				}
				for (double v = vmin; v <= vmax - sizev; v = v + sizev) {
					tmp = bd.getPointAtParam(umax, v, w);
					normal_singleLine.push_back(tmp);
				}
				tmp = bd.getPointAtParam(umax, vmax, w);
				normal_singleLine.push_back(tmp);
				normal_tmpvector->push_back(normal_singleLine);

			}
		}
		bd.calcuNormal();
		totalVerticesSize = vertices_data.size() / 6;
		singlepatch_vertics_number.push_back(count);
		time_t end = clock();
		printf("single bsplineSolidToDiscrete cost: %zdms, vertices size: %zd, \n", int(end - beg), totalVerticesSize);
	}
}
void ModelDraw::getControlMeshStatus(bool status) {
	this->isControlMesh = status;
}
void ModelDraw::showHexMesh(bool status) {
	this->isHexMesh = status;
}

void ModelDraw::set_isRealtime(bool status)//是否实时渲染
{
	this->isRealtime = status;
}
void ModelDraw::set_isShowinside(bool status)//是否渲染内部
{
	this->isShowinside = status;
}
void ModelDraw::set_chooseface(bool status)//是否进入选面模式
{
	this->ifchooseface = status;
}

void ModelDraw::set_F(double x, double y, double z) {
	this->My_Force.block_id = bp_index;
	this->My_Force.face_id = face_handle;
	this->My_Force.force.x() = x;
	this->My_Force.force.y() = y;
	this->My_Force.force.z() = z;
}
void ModelDraw::set_V(double x, double y, double z, int dir) {
	this->My_Offset.dir = dir;
	this->My_Offset.block_id = bp_index;
	this->My_Force.face_id = face_handle;
	this->My_Offset.force.x() = x;
	this->My_Offset.force.y() = y;
	this->My_Offset.force.z() = z;
}
//void ModelDraw::draw_BSplineSolid_HexMesh() {
//	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
//	//shaderProgram.setUniformValue("colormode", 606);
//	shaderProgram.setUniformValue("colormode", -1);
//	core->glDrawArrays(GL_QUADS, 0, vertices_data.size() / 6);
//	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
//		// 绘制线框
//		glEnable(GL_POLYGON_OFFSET_FILL);
//		glPolygonOffset(1.0, 1.0);
//		if (editmodel == EditModel::PICKITEM && pickType == PickType::PickEdge) {
//			glLineWidth(meshline_width * 2.5f);
//		}
//		else {
//			glLineWidth(meshline_width*0.5);
//		}
//		shaderProgram.setUniformValue("colormode", 530);
//		for (int i = 0; i < vertices_data.size() / 24; ++i) {
//			core->glDrawArrays(GL_LINE_LOOP, i * 4, 4);
//		}
//	}
//	update();
//}
void ModelDraw::draw_BSplineSolid_HexMesh() {
	core->glBindVertexArray(VAO);	// 前面已经记录了VBO和EBO现在只需一步调用即可
	//shaderProgram.setUniformValue("colormode", 606);
	//shaderProgram.setUniformValue("colormode", -1);
	int count = 0;
	for (int i = 0; i < singlepatch_vertics_number.size(); i++)
	{
		int num = 231 + i % 4;
		if (chooseface_index == i) {
			if (face_handle == 1) {
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count, bsplinesolidmesh[i].U0_size / 6);
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6, (singlepatch_vertics_number[i] - bsplinesolidmesh[i].U0_size) / 6);
			}
			else if (face_handle == 2) {
				
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count, bsplinesolidmesh[i].U0_size / 6);
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6, bsplinesolidmesh[i].U1_size / 6);
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6, (singlepatch_vertics_number[i] - bsplinesolidmesh[i].U0_size - bsplinesolidmesh[i].U1_size) / 6);

			}
			else if (face_handle == 3) {
				
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count, (bsplinesolidmesh[i].U0_size + bsplinesolidmesh[i].U1_size) / 6);
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6, bsplinesolidmesh[i].V0_size / 6);
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6, (singlepatch_vertics_number[i] - bsplinesolidmesh[i].U0_size - bsplinesolidmesh[i].U1_size - bsplinesolidmesh[i].V0_size) / 6);

			}
			else if (face_handle == 4) {
				
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count, (bsplinesolidmesh[i].U0_size + bsplinesolidmesh[i].U1_size + bsplinesolidmesh[i].V0_size) / 6);
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6, bsplinesolidmesh[i].V1_size / 6);
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6 + bsplinesolidmesh[i].V1_size / 6, (singlepatch_vertics_number[i] - bsplinesolidmesh[i].U0_size - bsplinesolidmesh[i].U1_size - bsplinesolidmesh[i].V0_size - bsplinesolidmesh[i].V1_size) / 6);

			}
			else if (face_handle == 5) {
				
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count, (bsplinesolidmesh[i].U0_size + bsplinesolidmesh[i].U1_size + bsplinesolidmesh[i].V0_size + bsplinesolidmesh[i].V1_size) / 6);
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6 + bsplinesolidmesh[i].V1_size / 6, bsplinesolidmesh[i].W0_size / 6);
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6 + bsplinesolidmesh[i].V1_size / 6 + bsplinesolidmesh[i].W0_size / 6, (singlepatch_vertics_number[i] - bsplinesolidmesh[i].U0_size - bsplinesolidmesh[i].U1_size - bsplinesolidmesh[i].V0_size - bsplinesolidmesh[i].V1_size - bsplinesolidmesh[i].W0_size) / 6);

			}
			else if (face_handle == 6) {
				
				shaderProgram.setUniformValue("colormode", num);
				core->glDrawArrays(GL_QUADS, count, (bsplinesolidmesh[i].U0_size + bsplinesolidmesh[i].U1_size + bsplinesolidmesh[i].V0_size + bsplinesolidmesh[i].V1_size + bsplinesolidmesh[i].W0_size) / 6);
				shaderProgram.setUniformValue("colormode", 606);
				core->glDrawArrays(GL_QUADS, count + bsplinesolidmesh[i].U0_size / 6 + bsplinesolidmesh[i].U1_size / 6 + bsplinesolidmesh[i].V0_size / 6 + bsplinesolidmesh[i].V1_size / 6 + bsplinesolidmesh[i].W0_size / 6, bsplinesolidmesh[i].W1_size / 6);

			}
		}


		else
		{
			shaderProgram.setUniformValue("colormode", num);
			core->glDrawArrays(GL_QUADS, count, singlepatch_vertics_number[i] / 6);
		}

		count += singlepatch_vertics_number[i] / 6;
	}

//	if (drawMeshLine && drawmode != DrawMode::POINTS && drawmode != DrawMode::BOUNDLINE) {
//		// 绘制线框
//		glEnable(GL_POLYGON_OFFSET_FILL);
//		glPolygonOffset(1.0, 1.0);
//		if (editmodel == EditModel::PICKITEM && pickType == PickType::PickEdge) {
//			glLineWidth(meshline_width * 2.5f);
//		}
//		else {
//			glLineWidth(meshline_width * 0.5);
//		}
//		shaderProgram.setUniformValue("colormode", 530);
//		for (int i = 0; i < vertices_data.size() / 24; ++i) {
//			core->glDrawArrays(GL_LINE_LOOP, i * 4, 4);
//		}
//	}
	update();
}
// 建议修改：将HexMesh存在solutionSkel里比较好
int ModelDraw::findNearestCell(double aim[3]) {
	iGameVertex aimV = iGameVertex(aim[0], aim[1], aim[2]);
	double minDistance = 999;
	iGameCellHandle record;
	for (auto& cp : hexMesh->allcells()) {		// this->hexMesh->boundary_face
		if (!hexMesh->isOnBoundary(cp.first)) continue;
		iGameVertex ccenter = hexMesh->getCellCenter(cp.first);
		double dist = ccenter.dist(aimV);
		if (dist < minDistance) {
			dist = minDistance;
			record = cp.first;
		}
	}
	return record;
}

// 延时函数
void ModelDraw::Sleep(int msec) {
	QTime dieTime = QTime::currentTime().addMSecs(msec);
	while (QTime::currentTime() < dieTime)
		QCoreApplication::processEvents(QEventLoop::AllEvents, 100);
}

int ModelDraw::loadMesh(std::string filename) {
	MeshFileIO io;
	int d, c;
	io.read_mesh_file(filename, d, c);
	std::vector<Vex> _vertices;
	std::vector<std::vector<VH>> _vhs;
	io.get_mesh_data(_vertices, _vhs);

	ModelClear();
	if (d == 2) {
		if (c == 3) {// 三角形网格
			surface_mesh = new MeshKernel::SurfaceMesh(_vertices, _vhs);
			std::cout << "Initialize success. #V = " << surface_mesh->vsize() << ", #F = " << surface_mesh->fsize() << std::endl;
			SetDrawModule(SURFACEMESH);
		}
	}
	else if (d == 3) {
		if (c == 6) {// 六面体网格
			hexMesh = new MeshKernel::VolumeMesh(_vertices, _vhs);
			std::cout << "Initialize success. #V = " << hexMesh->vsize() << ", #C = " << hexMesh->csize() << std::endl;
			SetDrawModule(HEXMESH);
		}
	}
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	update();

	return 0;
}

void ModelDraw::subdivisionNormal() {
	HexV3fMesh* temp = convent();
	BMesh::MLCASubdivision* mlcaSubdivision = new BMesh::MLCASubdivision();
	mlcaSubdivision->Spliting(*temp);
	free(mlcaSubdivision);
	FileHandler fileOut("crose_meshs.mesh");
	fileOut.writeFile("crose_meshs.mesh", temp);
	// 读取文件
	loadMesh("crose_meshs.mesh");
}

void ModelDraw::edgeSubdivision(int divnum) {
	if (selected_e_indices.size() != 1) return;
	static EdgeSubdivision eageSubdivision = EdgeSubdivision();
	for (int i = 0; i < selected_e_indices.size(); i++) {
		eageSubdivision.Spliting(hexMesh, selected_e_indices[i], divnum);
	}
	//hexMesh->updateAllHandles();
	selected_e_indices.clear();
	if (resetViewFlag) updateMeshModel();
	updateDrawData();
	update();
}

HexV3fMesh* ModelDraw::convent() {
	HexV3fMesh* _hexMesh = new HexV3fMesh();

	/* iGame_hexMesh to XIE_hexMesh */
	{
		uint64_t vcnt = hexMesh->VertexSize(), ccnt = hexMesh->CellSize();
		for (uint64_t i = 0; i < vcnt; ++i) {
			MeshKernel::iGameVertexHandle vh(i);
			auto& v = hexMesh->vertices(vh);
			V3f _v(v.x(), v.y(), v.z());
			_hexMesh->add_vertex(_v);
		}

		for (uint64_t i = 0; i < ccnt; ++i) {
			MeshKernel::iGameCellHandle ch(i);
			if (hexMesh->isValid(ch)) {
				auto& c = hexMesh->cells(ch);
				auto vhs = c.getVertexHandle();
				vector<VertexHandle> _vhs;
				for (auto& vh : vhs) {
					_vhs.push_back(VertexHandle(vh));
				}
				_hexMesh->add_cell(_vhs);
			}
		}

		/*std::cout << "iGameHexMesh convert to OldHexMesh success. OldHexMesh: #V = " << _hexMesh->vertices_.size() <<
		 ", #C = " << _hexMesh->cells_.size() << std::endl;*/
	}


	/*FileHandler file(filename);
	assert(file.read2HexMeshObj(_hexMesh, true, true, true));*/

	if (_hexMesh->vertices_.size() == 0) {
		QMessageBox::information(this, QString::fromLocal8Bit(" Warning "), QString::fromLocal8Bit(" HexMesh(XIE) is empty! "));
		return nullptr;
	}
	return _hexMesh;
}

void ModelDraw::generate_coverage_axis() {
    if (surface_mesh) {
        CoverageAxis app(*surface_mesh);
        app.compute_coverage_axis();

        skeletal_mesh = new SkeletalMesh();
        app.get_medial_mesh(*skeletal_mesh);
        SetDrawModule(SKELMESH);
        updateDrawData();
        //update_SurfaceMedialAxis();
    }

}

void ModelDraw::hex_meshing_test() {

//    hexMesh = new MeshKernel::VolumeMesh();
//
//    // 1. 生成面片上的六面体
//    QMessageBox::StandardButton ret1 = QMessageBox::question(this, " Question ", " Is non-manifold ? ", QMessageBox::Ok | QMessageBox::Cancel);
//    auto start_app = std::chrono::high_resolution_clock::now();
//    if (ret1 == QMessageBox::Cancel) {
//        HexMeshing_FaceSkeleton app0(*skeletal_mesh, *hexMesh);// 只适合于流形的
//        app0.hexmeshing();
//        updateDrawData();
//    } else {
//        HexMeshing_FaceSkeleton_Nonmanifold app0(*skeletal_mesh, *hexMesh);
//        app0.hexmeshing();
//        updateDrawData();
//    }
//
//    VH algorithm_insert(-1);
//    // 2. 在曲线分支节点上放置六面体
//    HexMeshing_Branch_Processor app1(*skeletal_mesh, *hexMesh);
//    app1.optimize(algorithm_insert, true);// 只放置, 不细分
//    updateDrawData();
//    auto point1 = std::chrono::high_resolution_clock::now();
//    QMessageBox::StandardButton ret = QMessageBox::question(this, " Question ", " Continue sweep hex ? ", QMessageBox::Ok | QMessageBox::Cancel);
//    if (ret == QMessageBox::Cancel) return;
//    auto point2 = std::chrono::high_resolution_clock::now();
//
//    // 3. 确定一致性细分并输出扫掠所需要的 源面, 目标面 以及 扫掠路径
//    HexMeshing_ConformingTessellation app2(*skeletal_mesh, *hexMesh);
//    app2.conforming_split(algorithm_insert);
//
//    auto end_app = std::chrono::high_resolution_clock::now();
//
//    std::chrono::duration<double> duration = end_app - start_app;
//    std::chrono::duration<double> duration_req = point2 - point1;
//
//    std::cout << "[HexMeshing]: cost " << (duration.count() - duration_req.count()) << "s\n";
//
//    updateDrawData();
//
//    /*line_cnt = 0;
//    unordered_set<EH> drawed;
//    for (auto& ep : app2.hexeh_pairs) {
//        auto& p1 = ep.first;
//        auto& p2 = ep.second;
//        vertices_data.push_back(p1.x()); vertices_data.push_back(p1.y()); vertices_data.push_back(p1.z());
//        vertices_data.push_back(0); vertices_data.push_back(0); vertices_data.push_back(0);
//        vertices_data.push_back(p2.x()); vertices_data.push_back(p2.y()); vertices_data.push_back(p2.z());
//        vertices_data.push_back(0); vertices_data.push_back(0); vertices_data.push_back(0);
//        line_cnt++;
//    }
//    elements_data.push_back(DrawElement(DrawElementType::T_LINE, -1, line_cnt, false, zFar_line));
//    std::cout << "The size of eh_pairs is " << app2.hexeh_pairs.size() << ".\n";*/
//
//    /*line_cnt = 0;
//    for (auto& vp : app2.matched_vex) {
//        auto& p1 = vp.first;
//        auto& p2 = vp.second;
//        vertices_data.push_back(p1.x()); vertices_data.push_back(p1.y()); vertices_data.push_back(p1.z());
//        vertices_data.push_back(1); vertices_data.push_back(0); vertices_data.push_back(0);
//        vertices_data.push_back(p2.x()); vertices_data.push_back(p2.y()); vertices_data.push_back(p2.z());
//        vertices_data.push_back(1); vertices_data.push_back(0); vertices_data.push_back(0);
//        line_cnt++;
//    }
//    elements_data.push_back(DrawElement(DrawElementType::T_LINE, -1, line_cnt, false, zFar_line));*/
//
//    /*int i = 0;
//    for (auto& ffhs : app2.rayid_to_sub_fhs) {
//        if (i % 2 == 0) {
//            for (auto& fhs : ffhs) {
//                for (auto& fh : fhs) {
//                    selected_f_indices.emplace_back(fh);
//                }
//            }
//        }
//        i++;
//    }*/

}

void ModelDraw::test_det_faces() {

//    if (!skeletal_mesh) return;
//
//    // 光滑分支结点
//    //BoundarySmoother::smoothing_branch(*skeletal_mesh);// 参数不好设置
//    // 抽取出三角形面片区域
//    vector<vector<VH>> domains;
//    vector<MeshKernel::SurfaceMesh> trimeshs;
//    DomainProcessor app;
//    //app.detect_faces(*skeletal_mesh, domains);// 把所有在域边界的VH组织成一串VH
//    //app.get_input_domains(*skeletal_mesh, domains, expanded_domains);// 偏移顶点, 得到最终的三维边界
//    ////app.preject_to_2d(*skeletal_mesh, domains, expanded_domains, trimeshs);
//    //app.remeshing_3d(*skeletal_mesh, domains, expanded_domains, trimeshs);
//    app.expand_boundary(*skeletal_mesh, expanded_domains);
//    app.decomposition(*skeletal_mesh);
//
//    updateDrawData();

}