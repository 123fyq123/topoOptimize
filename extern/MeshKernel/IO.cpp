#pragma once
#include "IO.h"

MeshKernel::SurfaceMesh MeshKernel::IO::ReadObjFile(const std::string& _InputFile, int& sides_num) {
    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> faces;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<double>> uvs;
    std::unordered_map<int, int> V2N;// vertex to normal
    std::unordered_map<int, int> V2T;// vertex to uv
    std::string line;

    std::cout << "Reading " << _InputFile << " File" << std::endl;
    // std::cout << inputfile.good() << std::endl;
    while (inputfile) {
        line.clear();
        getline(inputfile, line);
        if (line[0] == '#') {
            continue;//
        }
        std::stringstream linestream;
        linestream.str(line);

        std::string flag;
        linestream >> flag;
        if (flag == "v") {
            double x, y, z;
            linestream >> x >> y >> z;
            vertices.push_back(iGameVertex(x, y, z));
        } else if (flag == "f") {
            std::vector<std::string> vex;
            std::string tmp;
            while (linestream >> tmp) vex.push_back(tmp);
            auto n = vex.size();
            sides_num = n;
            std::vector<iGameVertexHandle> face(n);
            for (size_t i = 0; i < n; i++) {
                size_t idx = 0;
                while (idx < vex[i].length() && std::isdigit(vex[i][idx])) idx++;
                int vh = std::stoi(vex[i].substr(0, idx)) - 1;//
                face[i] = (iGameVertexHandle)(vh);
            }
            faces.push_back(face);
        } else if (flag == "vt") {
            double u, v;
            linestream >> u >> v;
            uvs.push_back({ u, v });
        } else if (flag == "vn") {
            double x, y, z;
            linestream >> x >> y >> z;
            normals.push_back({ x, y, z });
        }
    }
    //printf("read file success, fcnt: %d, vcnt: %d, vtcnt: %d, vncnt: %d\n", faces.size(), vertices.size(), uvs.size(), normals.size());
    if (!normals.empty()) {
        int ncnt = normals.size();
        for (int i = 0; i < vertices.size(); ++i) {
            int nidx = V2N[i];
            //if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug ��
            assert(nidx >= 0 && nidx < ncnt);
            vertices[i].setNormal(normals[nidx]);
        }
    }

    auto mesh = SurfaceMesh(vertices, faces);
    inputfile.close();
    return mesh;
}


std::string MeshKernel::IO::WriteOffString(const SurfaceMesh& _mesh) {
    std::stringstream ret;
    ReOrderiGameVertexHandle(_mesh);
    ret << "OFF" << std::endl;
    ret << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << _mesh.alledges().size() << std::endl;
    for (iGameVertexHandle vh : reorderedvh_) {
        iGameVertex v(_mesh.vertices(vh));
        ret << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }
    auto allf = _mesh.allfaces();
    for (auto f : allf) {
        ret << f.second.size();
        for (int i = 0; i < f.second.size(); ++i) {
            ret << " " << newvh_[f.second.vh(i)];
        }
        ret << std::endl;
    }
    return std::string(ret.str());
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadOffFile(const std::string& _InputFile, int& sides_num) {
    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> faces;
    std::vector<std::vector<double>> normals;
    std::vector<std::vector<double>> uvs;
    std::unordered_map<int, int> V2N;// vertex to normal
    std::unordered_map<int, int> V2T;// vertex to uv
    std::string line;
    int v_size, f_size, e_size;

    std::cout << "Reading " << _InputFile << " File" << std::endl;
    line.clear();
    getline(inputfile, line);
    std::stringstream linestream;

    if (line == "OFF") {
        line.clear();
        getline(inputfile, line);
        linestream.str(line);
        linestream >> v_size >> f_size >> e_size;
    }
    for (int i = 0; i < v_size; i++) {
        line.clear();
        getline(inputfile, line);
        std::stringstream linestream;
        linestream.str(line);
        double x, y, z;
        linestream >> x >> y >> z;
        vertices.push_back(iGameVertex(x, y, z));
        //std::cout << x << " " << y << " " << z << std::endl;
    }
    for (int i = 0; i < f_size; i++) {
        line.clear();
        getline(inputfile, line);
        std::stringstream linestream;
        linestream.str(line);
        int v;
        linestream >> v;
        //std::cout << v << " " << x << " " << y << " " << z << std::endl;
        std::vector<iGameVertexHandle> face(v);
        sides_num = v;
        for (int i = 0; i < v; i++) {
            int temp;
            linestream >> temp;
            face[i] = (iGameVertexHandle)(temp);
        }
        faces.push_back(face);
    }
    printf("read file success, fcnt: %d, vcnt: %d, vtcnt: %d, vncnt: %d\n", faces.size(), vertices.size(), uvs.size(), normals.size());
    if (!normals.empty()) {
        int ncnt = normals.size();
        for (int i = 0; i < vertices.size(); ++i) {
            int nidx = V2N[i];
            //if (nidx < 0 || nidx >= ncnt) printf("error: nidx = %d\n", nidx);// debug ��
            assert(nidx >= 0 && nidx < ncnt);
            vertices[i].setNormal(normals[nidx]);
        }
    }

    auto mesh = SurfaceMesh(vertices, faces);
    inputfile.close();
    return mesh;
}

MeshKernel::SurfaceMesh MeshKernel::IO::ReadStlFile(const std::string& _InputFile, int& sides_num) {
    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> faces;
    std::vector<iGameVertexHandle> face(3);
    std::vector<std::vector<double>> normals;
    std::string line;
    size_t vex_num = 0;  // 追踪当前facet中的顶点数量

    struct Vertex {
        double x, y, z;
        size_t id;  // 添加一个成员变量用于存储编号

        struct hash {   // 自定义哈希函数
            size_t operator()(const Vertex& v) const {
                return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
            }
        };
        struct equal {  // 自定义相等函数
            bool operator()(const Vertex& v1, const Vertex& v2) const {
                return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
            }
        };
    };
    std::unordered_set<Vertex, Vertex::hash, Vertex::equal> vertexSet;

    std::cout << "Reading " << _InputFile << " File" << std::endl;
    while (getline(inputfile, line)) {
        std::istringstream iss(line);
        std::string token;
        iss >> token;

        if (token == "facet" && (iss >> token) && token == "normal") {
            double x, y, z;
            iss >> x >> y >> z;
            std::vector<double> normal = {x, y, z};
            normals.push_back(normal);

            face.clear();
            face.resize(3);
        } else if (token == "vertex") {
            double x, y, z;
            iss >> x >> y >> z;
            iGameVertex v(x, y, z);
            Vertex vex{x, y, z, vertexSet.size()};  // 分配编号
            auto result = vertexSet.insert(vex);
            if (result.second) {
                vertices.push_back(v);
                face[vex_num] =(iGameVertexHandle) (vertexSet.size() - 1);
            } else {
                Vertex vex{x, y, z, 0};
                face[vex_num] =(iGameVertexHandle) (vertexSet.find(vex)->id);
            }
            vex_num++;
        } else if (token == "endfacet") {
            faces.push_back(face);
            vex_num = 0;
        }
    }
//    std::cout << "Vertices:" << std::endl;
//    for (iGameVertex& v : vertices) {
//        std::cout << "(" << v.x() << ", " << v.y() << ", " << v.z() << ") ";
//    }
//
//    std::cout << "\nFaces:" << std::endl;
//    for (const auto& f : faces) {
//        for (const auto& vh : f) {
//            std::cout << vh << " ";
//        }
//        std::cout << std::endl;
//    }

    auto mesh = SurfaceMesh(vertices, faces);
    inputfile.close();
    return mesh;
}

//MeshKernel::SurfaceMesh MeshKernel::IO::ReadStlFile(const std::string& _InputFile, int& sides_num) {
//
//	char const* file_path = _InputFile.c_str();
//
//	size_t file_size;
//
//	FILE* file = fopen(file_path, "rb");
//	if (!file) {
//		fprintf(stderr, "file %s does not exist\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//
//	fseek(file, 0, SEEK_END);
//	file_size = ftell(file);
//	rewind(file);
//
//	char* file_contents = (char*)malloc(file_size);
//	if (fread(file_contents, 1, file_size, file) != file_size) {
//		fprintf(stderr, "I/O error while reading %s\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//	fclose(file);
//
//	size_t buf_size;
//	if (one_stl_buf_size(&buf_size, file_contents, file_size, ONE_STL_NVVV)) {
//		fprintf(stderr, "file \"%s\" is not valid stl data\n", file_path);
//		exit(EXIT_FAILURE);
//	}
//
//	float* buf = (float*)malloc(buf_size);
//	size_t trig_count = one_stl_parse(buf, file_contents, ONE_STL_NVVV);
//
//	printf("triangle count: %ld\n", trig_count);
//
//	typedef struct {
//		float normal[3];
//		float v0[3];
//		float v1[3];
//		float v2[3];
//	} trig_nvvv_s;
//	trig_nvvv_s* triangles = (trig_nvvv_s*)buf;
//	float EPS = 0;
//	//for (size_t i = 0; i < trig_count; ++i){
//	//    EPS=std::max(EPS,(triangles + i)->v0[0]);
//	//    EPS=std::max(EPS,(triangles + i)->v0[1]);
//	//    EPS=std::max(EPS,(triangles + i)->v0[2]);
//	//}
//	//EPS/=100000;
//	//EPS = 0;
//
//	struct PointEps {
//		PointEps() {}
//		PointEps(float a, float b, float c, float EPS) {
//			v0[0] = a; v0[1] = b;
//			v0[2] = c;
//			this->EPS = EPS;
//		}
//		float v0[3];
//		float EPS;
//		bool operator < (const PointEps& other) const {
//			for (int i = 0; i < 3; i++) {
//				if (abs(v0[i] - other.v0[i]) > EPS) {
//					return  v0[i] + EPS < other.v0[i];
//				}
//			}
//			return false;
//		}
//	};
//
//	std::map<PointEps, int>mp;
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	std::function<int(float, float, float)> getHandleId = [&](float x, float y, float z) {
//		if (mp.count(PointEps(x, y, z, EPS))) {
//			return mp[PointEps(x, y, z, EPS)];
//		}
//		int cnt = mp.size();
//		mp[PointEps(x, y, z, EPS)] = cnt;
//		vertices.push_back(MeshKernel::iGameVertex{ x,y,z });
//		return cnt;
//	};
//	std::vector<std::vector<float> >normal;
//	for (size_t i = 0; i < trig_count; ++i) {
//		faces.push_back({ (iGameVertexHandle)getHandleId((triangles + i)->v0[0],(triangles + i)->v0[1],(triangles + i)->v0[2]),
//						(iGameVertexHandle)getHandleId((triangles + i)->v1[0],(triangles + i)->v1[1],(triangles + i)->v1[2]),
//						(iGameVertexHandle)getHandleId((triangles + i)->v2[0],(triangles + i)->v2[1],(triangles + i)->v2[2])
//			});
//		normal.push_back(std::vector<float>{ (triangles + i)->normal[0], (triangles + i)->normal[1], (triangles + i)->normal[2] });
//	}
//	std::cout << vertices.size() << " " << faces.size() << std::endl;
//	auto mesh = SurfaceMesh(vertices, faces);
//	std::cout << trig_count << " " << faces.size() << std::endl;
//	//for(int i=0;i<faces.size();i++){
//	//    mesh.faces(iGameFcaeHandle(i)).setNormal(normal[i][0], normal[i][1], normal[i][2]);
//	//}
//	free(file_contents);
//	free(buf);
//	return mesh;
//}
MeshKernel::SurfaceMesh MeshKernel::IO::ReadIgesFile(const std::string& _InputFile, int& sides_num) {
    // 创建IGES文件读取器
    IGESControl_Reader aReader_Iges;
    // 设置IGES文件路径
    aReader_Iges.ReadFile(_InputFile.c_str());

    IFSelect_PrintCount mode = IFSelect_ItemsByEntity;
    Interface_Static::SetCVal("xstep.cascade.unit", "M");
    aReader_Iges.PrintCheckLoad(Standard_False, mode);

    Handle_IGESData_IGESModel model = aReader_Iges.IGESModel();

    Standard_Integer modelNum = aReader_Iges.NbShapes();

    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> faces;
    Standard_Integer Vertex_size = 0, Face_size = 0;

    struct Vertex {
        double x, y, z;
        size_t id;

        struct hash {
            size_t operator()(const Vertex& v) const {
                return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
            }
        };
        struct equal {
            bool operator()(const Vertex& v1, const Vertex& v2) const {
                return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
            }
        };
    };
    std::unordered_set<Vertex, Vertex::hash, Vertex::equal> vertexSet;
    std::vector<iGameVertex> vertexSum;

    std::cout << "Reading " << _InputFile << " File" << std::endl;

    for (Standard_Integer i = 1; i <= modelNum; i++) {
        TopoDS_Shape cur = aReader_Iges.Shape(i);
        Bnd_Box bounds;
        BRepBndLib::Add(cur, bounds);
        bounds.SetGap(0.1);
        Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
        bounds.Get(xMin, yMin, zMin, xMax, yMax, zMax);
        Standard_Real deflection = ((xMax - xMin) + (yMax - yMin) + (zMax - zMin)) / 300.0 * 0.5;
        BRepMesh_IncrementalMesh(cur, deflection);
        for (TopExp_Explorer expIter(cur, TopAbs_FACE); expIter.More(); expIter.Next()) {
            const TopoDS_Shape& curFace = expIter.Current();
            TopLoc_Location aLoc;
            Handle_Poly_Triangulation mesh = BRep_Tool::Triangulation(TopoDS::Face(curFace), aLoc);
            if (!mesh.IsNull()) {
                Vertex_size += mesh->NbNodes();
                Face_size += mesh->NbTriangles();
            }
        }
    }

    Handle(Poly_Triangulation) tempMesh = new Poly_Triangulation(Vertex_size, Face_size, Standard_False);
    Standard_Integer nodeOffset = 0;
    Standard_Integer triangleOffset = 0;

    for (Standard_Integer i = 1; i <= modelNum; i++) {
        TopoDS_Shape cur = aReader_Iges.Shape(i);
        for (TopExp_Explorer expIter(cur, TopAbs_FACE); expIter.More(); expIter.Next()) {
            const TopoDS_Shape& curFace = expIter.Current();
            TopLoc_Location aLoc;
            Handle_Poly_Triangulation mesh = BRep_Tool::Triangulation(TopoDS::Face(curFace), aLoc);
            gp_Trsf aTrsf = aLoc.Transformation();
            for (int j = 1; j <= mesh->NbNodes(); j++) {
                Standard_Integer id = (Standard_Integer)j;
                auto tmp = mesh->Node(id);
                tmp.Transform(aTrsf);
                tempMesh->SetNode(id + nodeOffset, tmp);
            }
            for (int j = 1; j <= mesh->NbTriangles(); j++) {
                Poly_Triangle tempTri = mesh->Triangle(j);  //.hxx unlocked
                Standard_Integer N1, N2, N3;
                tempTri.Get(N1, N2, N3);
                if (curFace.Orientation() == TopAbs_REVERSED) {
                    Standard_Integer tmp = N1;
                    N1 = N2;
                    N2 = tmp;
                }
                N1 += nodeOffset;
                N2 += nodeOffset;
                N3 += nodeOffset;
                tempTri.Set(N1, N2, N3);
                tempMesh->SetTriangle((Standard_Integer)j + triangleOffset, tempTri);
            }
            nodeOffset += mesh->NbNodes();
            triangleOffset += mesh->NbTriangles();
        }
    }

    for (int i = 1; i <= tempMesh->NbNodes(); i++) {
        Standard_Integer id = (Standard_Integer)i;
        auto tmp = tempMesh->Node(id);
        //std::cout << tmp.X() << " " << tmp.Y() << " " << tmp.Z() << std::endl;
        iGameVertex v(tmp.X(), tmp.Y(), tmp.Z());
        Vertex vex{tmp.X(), tmp.Y(), tmp.Z(), vertexSet.size()};
        auto result = vertexSet.insert(vex);
        if (result.second) {
            vertices.push_back(v);
        }
        vertexSum.push_back(v);
    }
    for (int i = 1; i <= tempMesh->NbTriangles(); i++)
    {
        Poly_Triangle tempTri = tempMesh->Triangle(i);  //.hxx unlocked
        Standard_Integer N1, N2, N3;
        tempTri.Get(N1, N2, N3);
        N1 = vertexSet.find({vertexSum[N1-1].x(), vertexSum[N1-1].y(), vertexSum[N1-1].z(), 0})->id;
        N2 = vertexSet.find({vertexSum[N2-1].x(), vertexSum[N2-1].y(), vertexSum[N2-1].z(), 0})->id;
        N3 = vertexSet.find({vertexSum[N3-1].x(), vertexSum[N3-1].y(), vertexSum[N3-1].z(), 0})->id;
        //std::cout << N1 << " " << N2 << " " << N3 << std::endl;
        std::vector<iGameVertexHandle> face = {(iGameVertexHandle)N1, (iGameVertexHandle)N2, (iGameVertexHandle)N3};
        faces.push_back(face);
    }

    auto mesh = SurfaceMesh(vertices, faces);
    return mesh;
}


MeshKernel::SurfaceMesh MeshKernel::IO::ReadStepFile(const std::string& _InputFile, int& sides_num) {
    STEPControl_Reader aReader_Step;
    //Standard_CString tmp = (Standard_CString)filename_in;
    aReader_Step.ReadFile(_InputFile.c_str());
    //检查文件加载状态
    aReader_Step.PrintCheckLoad(Standard_False, IFSelect_ItemsByEntity);
    ////加载step文件
    Standard_Integer NbRoots = aReader_Step.NbRootsForTransfer();// 获取可转移根的数量
    Standard_Integer num = aReader_Step.TransferRoots();//翻译所有可转换的根，并返回//成功翻译的次数

    Standard_Integer modelNum = aReader_Step.NbShapes();
    //std::cout << "读取到 ：" << modelNum << "个零件" << std::endl;
    // read Face
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> faces;
    Standard_Integer Vertex_size = 0, Face_size = 0;

    struct Vertex {
        double x, y, z;
        size_t id;  // 添加一个成员变量用于存储编号

        struct hash {   // 自定义哈希函数
            size_t operator()(const Vertex& v) const {
                return std::hash<float>()(v.x) ^ std::hash<float>()(v.y) ^ std::hash<float>()(v.z);
            }
        };
        struct equal {  // 自定义相等函数
            bool operator()(const Vertex& v1, const Vertex& v2) const {
                return v1.x == v2.x && v1.y == v2.y && v1.z == v2.z;
            }
        };
    };
    std::unordered_set<Vertex, Vertex::hash, Vertex::equal> vertexSet;
    std::vector<iGameVertex> vertexSum;

    std::cout << "Reading " << _InputFile << " File" << std::endl;

    for (Standard_Integer i = 1; i <= modelNum; i++) {
        //printf("第%d个部件\n", i);
        TopoDS_Shape cur = aReader_Step.Shape(i);
        Bnd_Box bounds;
        BRepBndLib::Add(cur, bounds);
        bounds.SetGap(0.1);
        Standard_Real xMin, yMin, zMin, xMax, yMax, zMax;
        bounds.Get(xMin, yMin, zMin, xMax, yMax, zMax);
        Standard_Real deflection = ((xMax - xMin) + (yMax - yMin) + (zMax - zMin)) / 300.0 * 0.5;
        BRepMesh_IncrementalMesh(cur, deflection);
        for (TopExp_Explorer expIter(cur, TopAbs_FACE); expIter.More(); expIter.Next()) {
            const TopoDS_Shape& curFace = expIter.Current();
            TopLoc_Location aLoc;
            Handle_Poly_Triangulation mesh = BRep_Tool::Triangulation(TopoDS::Face(curFace), aLoc);
            if (!mesh.IsNull()) {
                Vertex_size += mesh->NbNodes();
                Face_size += mesh->NbTriangles();
            }
        }
    }

    //printf("Vertex_size = %d , Face_size = %d\n", Vertex_size, Face_size);
    Handle(Poly_Triangulation) tempMesh = new Poly_Triangulation(Vertex_size, Face_size, Standard_False);
    Standard_Integer nodeOffset = 0;
    Standard_Integer triangleOffset = 0;

    for (Standard_Integer i = 1; i <= modelNum; i++) {
        //printf("第%d个部件\n", i);
        TopoDS_Shape cur = aReader_Step.Shape(i);
        for (TopExp_Explorer expIter(cur, TopAbs_FACE); expIter.More(); expIter.Next()) {
            const TopoDS_Shape& curFace = expIter.Current();
            TopLoc_Location aLoc;
            Handle_Poly_Triangulation mesh = BRep_Tool::Triangulation(TopoDS::Face(curFace), aLoc);
            gp_Trsf aTrsf = aLoc.Transformation();
            for (int j = 1; j <= mesh->NbNodes(); j++) {
                Standard_Integer id = (Standard_Integer)j;
                auto tmp = mesh->Node(id);
                tmp.Transform(aTrsf);
                tempMesh->SetNode(id + nodeOffset, tmp);
            }
            for (int j = 1; j <= mesh->NbTriangles(); j++) {
                Poly_Triangle tempTri = mesh->Triangle(j);  //.hxx unlocked
                Standard_Integer N1, N2, N3;
                tempTri.Get(N1, N2, N3);
                if (curFace.Orientation() == TopAbs_REVERSED) {
                    Standard_Integer tmp = N1;
                    N1 = N2;
                    N2 = tmp;
                }
                N1 += nodeOffset;
                N2 += nodeOffset;
                N3 += nodeOffset;
                tempTri.Set(N1, N2, N3);
                tempMesh->SetTriangle((Standard_Integer)j + triangleOffset, tempTri);
            }
            nodeOffset += mesh->NbNodes();
            triangleOffset += mesh->NbTriangles();
        }
    }

    for (int i = 1; i <= tempMesh->NbNodes(); i++) {
        Standard_Integer id = (Standard_Integer)i;
        auto tmp = tempMesh->Node(id);
        //std::cout << tmp.X() << " " << tmp.Y() << " " << tmp.Z() << std::endl;
        iGameVertex v(tmp.X(), tmp.Y(), tmp.Z());
        Vertex vex{tmp.X(), tmp.Y(), tmp.Z(), vertexSet.size()};  // 分配编号
        auto result = vertexSet.insert(vex);
        if (result.second) {
            vertices.push_back(v);
        }
        vertexSum.push_back(v);
    }
    for (int i = 1; i <= tempMesh->NbTriangles(); i++)
    {
        Poly_Triangle tempTri = tempMesh->Triangle(i);  //.hxx unlocked
        Standard_Integer N1, N2, N3;
        tempTri.Get(N1, N2, N3);
        N1 = vertexSet.find({vertexSum[N1-1].x(), vertexSum[N1-1].y(), vertexSum[N1-1].z(), 0})->id;
        N2 = vertexSet.find({vertexSum[N2-1].x(), vertexSum[N2-1].y(), vertexSum[N2-1].z(), 0})->id;
        N3 = vertexSet.find({vertexSum[N3-1].x(), vertexSum[N3-1].y(), vertexSum[N3-1].z(), 0})->id;
        //std::cout << N1 << " " << N2 << " " << N3 << std::endl;
        std::vector<iGameVertexHandle> face = {(iGameVertexHandle)N1, (iGameVertexHandle)N2, (iGameVertexHandle)N3};
        //从零开始编码
        faces.push_back(face);
    }

    auto mesh = SurfaceMesh(vertices, faces);
    return mesh;
}

bool MeshKernel::IO::WriteObjFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    ReOrderiGameVertexHandle(_mesh);
    for (iGameVertexHandle vh : reorderedvh_) {
        iGameVertex v(_mesh.vertices(vh));
        outputfile << "v " << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }
    auto allf = _mesh.allfaces();
    for (auto f : allf) {
        outputfile << "f";
        for (int i = 0; i < f.second.size(); ++i) {
            outputfile << " " << newvh_[f.second.vh(i)] + 1;
        }
        outputfile << std::endl;
    }
    outputfile.close();
    return true;
}

bool MeshKernel::IO::WriteObjFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& outputFile) {
    std::ofstream outputfile(outputFile, std::ios::out);
    if (!outputfile.is_open()) {
        std::cerr << "Failed to open output file: " << outputFile << std::endl;
        return false;
    }

    if (surfacesDraw.empty()) {
        std::cerr << "Error: surfacesDraw is empty!" << std::endl;
        return false;
    }

    struct VertexHash {
        size_t operator()(const Eigen::Vector3d& v) const {
            return std::hash<double>()(v.x()) ^ std::hash<double>()(v.y()) ^ std::hash<double>()(v.z());
        }
    };

    struct VertexEqual {
        bool operator()(const Eigen::Vector3d& v1, const Eigen::Vector3d& v2) const {
            const double epsilon = 1e-9;
            bool flag = std::fabs(v1.x() - v2.x()) <= epsilon &&
                        std::fabs(v1.y() - v2.y()) <= epsilon &&
                        std::fabs(v1.z() - v2.z()) <= epsilon;
            return flag;
        }
    };

    std::vector<Eigen::Vector3d> vertices;
    std::vector<std::vector<int>> faces;
    std::unordered_map<Eigen::Vector3d, int, VertexHash> vertexMap;

    auto findVertex = [&](const Eigen::Vector3d& vertex) -> int {
        for (const auto& pair : vertexMap) {
            if (VertexEqual()(pair.first, vertex)) {
                return pair.second;
            }
        }
        return -1;
    };

    auto addVertex = [&](const Eigen::Vector3d& vertex) -> int {
        int vertexIndex = findVertex(vertex);
        if (vertexIndex == -1) {
            vertexIndex = vertices.size() + 1;
            vertexMap[vertex] = vertexIndex;
            vertices.push_back(vertex);
        }
        return vertexIndex;
    };

    for (size_t i = 0; i < surfacesDraw.size(); ++i) {
        double umin, umax, vmin, vmax;
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);

        int unumber = 4; // u方向控制点数量
        int vnumber = 4; // v方向控制点数量

        double sizeu = (umax - umin) / (unumber - 1); // u方向步长
        double sizev = (vmax - vmin) / (vnumber - 1); // v方向步长

        for (int u_index = 0; u_index < unumber; ++u_index) {
            for (int v_index = 0; v_index < vnumber; ++v_index) {
                CPoint pnt = surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + v_index * sizev);
                Eigen::Vector3d vertex(pnt.getX(), pnt.getY(), pnt.getZ());

                int vertexIndex = findVertex(vertex);
                if (vertexIndex == -1) {
                    vertexIndex = vertices.size() + 1;
                    vertexMap[vertex] = vertexIndex;
                    vertices.push_back(vertex);
                }
            }
        }

        for (int u_index = 0; u_index < unumber - 1; ++u_index) {
            for (int v_index = 0; v_index < vnumber - 1; ++v_index) {
                Eigen::Vector3d v0 = Eigen::Vector3d(
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + v_index * sizev).getX(),
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + v_index * sizev).getY(),
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + v_index * sizev).getZ()
                );
                Eigen::Vector3d v1 = Eigen::Vector3d(
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + v_index * sizev).getX(),
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + v_index * sizev).getY(),
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + v_index * sizev).getZ()
                );
                Eigen::Vector3d v2 = Eigen::Vector3d(
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + (v_index + 1) * sizev).getX(),
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + (v_index + 1) * sizev).getY(),
                        surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + (v_index + 1) * sizev).getZ()
                );
                Eigen::Vector3d v3 = Eigen::Vector3d(
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + (v_index + 1) * sizev).getX(),
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + (v_index + 1) * sizev).getY(),
                        surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + (v_index + 1) * sizev).getZ()
                );

                faces.push_back({addVertex(v0), addVertex(v1), addVertex(v2)});
                faces.push_back({addVertex(v0), addVertex(v2), addVertex(v3)});
            }
        }
    }

    for (const auto& vertex : vertices) {
        outputfile << "v " << vertex.x() << " " << vertex.y() << " " << vertex.z() << std::endl;
    }

    for (const auto& face : faces) {
        outputfile << "f";
        for (const auto& idx : face) {
            outputfile << " " << idx;
        }
        outputfile << std::endl;
    }

    outputfile.close();
    return true;
}

bool MeshKernel::IO::WriteOffFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    ReOrderiGameVertexHandle(_mesh);
    outputfile << "OFF" << std::endl;
    outputfile << _mesh.allvertices().size() << " " << _mesh.allfaces().size() << " " << _mesh.alledges().size() << std::endl;
    for (iGameVertexHandle vh : reorderedvh_) {
        iGameVertex v(_mesh.vertices(vh));
        outputfile << v.x() << " " << v.y() << " " << v.z() << std::endl;
    }
    auto allf = _mesh.allfaces();
    for (auto f : allf) {
        outputfile << f.second.size();
        for (int i = 0; i < f.second.size(); ++i) {
            outputfile << " " << newvh_[f.second.vh(i)];
        }
        outputfile << std::endl;
    }
    outputfile.close();
    return true;
}

bool MeshKernel::IO::WriteStlFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    ReOrderiGameVertexHandle(_mesh);
    int front = _OutputFile.find_last_of('/');
    std::string output = _OutputFile.substr(++front, _OutputFile.length()).c_str();

    auto allf = _mesh.allfaces();
    std::vector<Eigen::Vector3d> Normals;//每个面的法向量
    Normals.resize(allf.size(), Eigen::Vector3d::Zero());

    outputfile << "solid " << output << std::endl;
    for (auto f : allf) {
        Eigen::Vector3d Normals;//面的法向量
        auto vhs = f.second.getVertexHandle();//f.second代表的就是allfaces<map>中第二个参数，也就是face的实体
        assert(vhs.size() > 2);
        Eigen::Vector3d pos1(_mesh.vertices(vhs[0]).x(), _mesh.vertices(vhs[0]).y(), _mesh.vertices(vhs[0]).z());
        Eigen::Vector3d pos2(_mesh.vertices(vhs[1]).x(), _mesh.vertices(vhs[1]).y(), _mesh.vertices(vhs[1]).z());
        Eigen::Vector3d pos3(_mesh.vertices(vhs[2]).x(), _mesh.vertices(vhs[2]).y(), _mesh.vertices(vhs[2]).z());
        Normals = (pos3 - pos2).cross(pos1 - pos2);//cross叉乘
        Normals.normalize();

        outputfile << "  facet normal " << Normals.x() << " " << Normals.y() << " " << Normals.z() << std::endl;
        outputfile << "    outer loop" << std::endl;
        outputfile << "      vertex " << pos1.x() << " " << pos1.y() << " " << pos1.z() << std::endl;
        outputfile << "      vertex " << pos2.x() << " " << pos2.y() << " " << pos2.z() << std::endl;
        outputfile << "      vertex " << pos3.x() << " " << pos3.y() << " " << pos3.z() << std::endl;
        outputfile << "    endloop" << std::endl;
        outputfile << "  endfacet" << std::endl;
    }
    outputfile << "endsolid " << output << std::endl;

    outputfile.close();
    return true;
}

bool MeshKernel::IO::WriteStlFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& outputFile) {
    std::ofstream outputfile(outputFile, std::ios::out);
    if (!outputfile.is_open()) {
        std::cerr << "Failed to open output file: " << outputFile << std::endl;
        return false;
    }

    if (surfacesDraw.empty()) {
        std::cerr << "Error: surfacesDraw is empty!" << std::endl;
        return false;
    }

    int front = outputFile.find_last_of('/');
    std::string output = outputFile.substr(++front, outputFile.length() - front);

    outputfile << "solid " << output << std::endl;

    for (size_t i = 0; i < surfacesDraw.size(); i++) {
        double umin, umax, vmin, vmax;
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);

        int unumber = 4; // u direction control points
        int vnumber = 4; // v direction control points

        double sizeu = (umax - umin) / (unumber - 1); // step size in u direction
        double sizev = (vmax - vmin) / (vnumber - 1); // step size in v direction

        for (int u_index = 0; u_index < unumber - 1; ++u_index) {
            for (int v_index = 0; v_index < vnumber - 1; ++v_index) {
                std::vector<Eigen::Vector3d> vertices(4);
                CPoint pnt0 = surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + v_index * sizev);
                CPoint pnt1 = surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + v_index * sizev);
                CPoint pnt2 = surfacesDraw[i].getPointAtParam(umin + (u_index + 1) * sizeu, vmin + (v_index + 1) * sizev);
                CPoint pnt3 = surfacesDraw[i].getPointAtParam(umin + u_index * sizeu, vmin + (v_index + 1) * sizev);
                vertices[0] = Eigen::Vector3d(pnt0.getX(), pnt0.getY(), pnt0.getZ());
                vertices[1] = Eigen::Vector3d(pnt1.getX(), pnt1.getY(), pnt1.getZ());
                vertices[2] = Eigen::Vector3d(pnt2.getX(), pnt2.getY(), pnt2.getZ());
                vertices[3] = Eigen::Vector3d(pnt3.getX(), pnt3.getY(), pnt3.getZ());

                Eigen::Vector3d normal = (vertices[2] - vertices[1]).cross(vertices[0] - vertices[1]);
                normal.normalize();

                // First triangle
                outputfile << "  facet normal " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
                outputfile << "    outer loop" << std::endl;
                outputfile << "      vertex " << vertices[0].x() << " " << vertices[0].y() << " " << vertices[0].z() << std::endl;
                outputfile << "      vertex " << vertices[1].x() << " " << vertices[1].y() << " " << vertices[1].z() << std::endl;
                outputfile << "      vertex " << vertices[2].x() << " " << vertices[2].y() << " " << vertices[2].z() << std::endl;
                outputfile << "    endloop" << std::endl;
                outputfile << "  endfacet" << std::endl;

                // Second triangle
                outputfile << "  facet normal " << normal.x() << " " << normal.y() << " " << normal.z() << std::endl;
                outputfile << "    outer loop" << std::endl;
                outputfile << "      vertex " << vertices[0].x() << " " << vertices[0].y() << " " << vertices[0].z() << std::endl;
                outputfile << "      vertex " << vertices[2].x() << " " << vertices[2].y() << " " << vertices[2].z() << std::endl;
                outputfile << "      vertex " << vertices[3].x() << " " << vertices[3].y() << " " << vertices[3].z() << std::endl;
                outputfile << "    endloop" << std::endl;
                outputfile << "  endfacet" << std::endl;
            }
        }
    }

    outputfile << "endsolid " << output << std::endl;
    outputfile.close();
    return true;
}

//bool MeshKernel::IO::WriteStlFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile) {
//    std::ofstream outputfile(_OutputFile, std::ios::out);
//    if (!outputfile.is_open()) {
//        std::cerr << "Cannot open output file!" << std::endl;
//        return false;
//    }
//
//    double umin, umax, vmin, vmax;
//    vector<CPoint> singlePatch; // 每个面片的离散点
//    vector<vector<CPoint>> singlelSurface; // 每个曲面的离散点
//    vector<vector<vector<CPoint>>> allSurface; // 所有曲面集合的离散点
//
//    vector<CVector> singlePatch_normal;
//    vector<vector<CVector>> singlelSurface_normal;
//    vector<vector<vector<CVector>>> allSurface_normal;
//    CVector normal; // 保存法向量
//
//    STEPControl_Writer writer;
//
//    for (size_t i = 0; i < surfacesDraw.size(); i++) {
//        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
//
//        int unumber = 4; // u方向的控制顶点数
//        int vnumber = 4; // v方向的控制顶点数
//
//        double sizeu = (umax - umin) / (unumber - 1); // 步距
//        double sizev = (vmax - vmin) / (vnumber - 1);
//
//        // 构建样条表面点阵
//        TColgp_Array2OfPnt points(1, unumber, 1, vnumber);
//        for (int u_index = 0; u_index < unumber; ++u_index) {
//            for (int v_index = 0; v_index < vnumber; ++v_index) {
//                double u = umin + u_index * sizeu;
//                double v = vmin + v_index * sizev;
//                auto pnt = surfacesDraw[i].getPointAtParam(u, v);
//                points.SetValue(u_index + 1, v_index + 1, gp_Pnt(pnt.getX(), pnt.getY(), pnt.getZ()));
//            }
//        }
//
//        // 创建样条表面
//        GeomAPI_PointsToBSplineSurface splineBuilder(points, 3, 3, GeomAbs_C2);
//        Handle(Geom_BSplineSurface) splineSurface = splineBuilder.Surface();
//
//        if (splineSurface.IsNull()) {
//            std::cerr << "Failed to create BSpline surface!" << std::endl;
//            return false;
//        }
//
//        // 写入 STL 文件
//        outputfile << "solid " << _OutputFile << std::endl;
//
//        for (int u_index = 1; u_index < unumber; ++u_index) {
//            for (int v_index = 1; v_index < vnumber; ++v_index) {
//                gp_Pnt p1 = splineSurface->Value(u_index, v_index);
//                gp_Pnt p2 = splineSurface->Value(u_index + 1, v_index);
//                gp_Pnt p3 = splineSurface->Value(u_index + 1, v_index + 1);
//                gp_Pnt p4 = splineSurface->Value(u_index, v_index + 1);
//
//                Eigen::Vector3d pos1(p1.X(), p1.Y(), p1.Z());
//                Eigen::Vector3d pos2(p2.X(), p2.Y(), p2.Z());
//                Eigen::Vector3d pos3(p3.X(), p3.Y(), p3.Z());
//                Eigen::Vector3d pos4(p4.X(), p4.Y(), p4.Z());
//
//                Eigen::Vector3d normal1 = (pos2 - pos1).cross(pos3 - pos1);
//                Eigen::Vector3d normal2 = (pos4 - pos3).cross(pos1 - pos3);
//
//                outputfile << "  facet normal " << normal1.x() << " " << normal1.y() << " " << normal1.z() << std::endl;
//                outputfile << "    outer loop" << std::endl;
//                outputfile << "      vertex " << pos1.x() << " " << pos1.y() << " " << pos1.z() << std::endl;
//                outputfile << "      vertex " << pos2.x() << " " << pos2.y() << " " << pos2.z() << std::endl;
//                outputfile << "      vertex " << pos3.x() << " " << pos3.y() << " " << pos3.z() << std::endl;
//                outputfile << "    endloop" << std::endl;
//                outputfile << "  endfacet" << std::endl;
//
//                outputfile << "  facet normal " << normal2.x() << " " << normal2.y() << " " << normal2.z() << std::endl;
//                outputfile << "    outer loop" << std::endl;
//                outputfile << "      vertex " << pos3.x() << " " << pos3.y() << " " << pos3.z() << std::endl;
//                outputfile << "      vertex " << pos4.x() << " " << pos4.y() << " " << pos4.z() << std::endl;
//                outputfile << "      vertex " << pos1.x() << " " << pos1.y() << " " << pos1.z() << std::endl;
//                outputfile << "    endloop" << std::endl;
//                outputfile << "  endfacet" << std::endl;
//            }
//        }
//
//        outputfile << "endsolid " << _OutputFile << std::endl;
//    }
//
//    outputfile.close();
//    return true;
//}

bool MeshKernel::IO::WriteStepFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    ReOrderiGameVertexHandle(_mesh);
    int front = _OutputFile.find_last_of('/');
    std::string output = _OutputFile.substr(++front, _OutputFile.length()).c_str();

    auto allf = _mesh.allfaces();
    auto allv = _mesh.allvertices();
    std::vector<Eigen::Vector3d> Normals;//每个面的法向量
    Normals.resize(allf.size(), Eigen::Vector3d::Zero());

    STEPControl_Writer writer;// 创建STEP文件写入器
    BRep_Builder builder; // 创建顶点
    const Standard_Real tol = 1.0e-6;//公差

    std::vector<TopoDS_Vertex> vertexList;
    for (const auto& vertex : allv) {
        auto vs = vertex.second;
        TopoDS_Vertex v;
        builder.MakeVertex(v, gp_Pnt(vs.x(), vs.y(), vs.z()), tol);
        std::cout << vs.x() << " " << vs.y() << " " << vs.z() << std::endl;
        writer.Transfer(v, STEPControl_AsIs);
        vertexList.push_back(v);
    }
    for (const auto& face : allf) {
        BRepBuilderAPI_MakePolygon polygonBuilder;
        auto vhs = face.second.getVertexHandle();
        for (auto vertex : vhs) {
            polygonBuilder.Add(vertexList[vertex]); // 使用之前创建的顶点
            std::cout << vertex << std::endl;
        }
        if (polygonBuilder.IsDone()) {
            TopoDS_Wire wire = polygonBuilder.Wire();
            if (!wire.IsNull()) {
                BRepBuilderAPI_MakeFace makeFace(wire, true);
                writer.Transfer(makeFace.Face(), STEPControl_AsIs); // 添加线到STEP控制器
            }
        }
    }

    const char* filename_cstr = _OutputFile.c_str();
    Standard_Boolean status = writer.Write(filename_cstr);// 写入STEP文件
    if (!status) {
        std::cout << "Failed to write STEP file!" << std::endl;
    } else {
        std::cout << "STEP file successfully written!" << std::endl;
    }

    outputfile.close();
    return true;
}

//bool MeshKernel::IO::WriteStepFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
//    std::ofstream outputfile(_OutputFile, std::ios::out);
//    if (!outputfile.is_open()) {
//        std::cout << "Cannot open output file!" << std::endl;
//        return false;
//    }
//
//    ReOrderiGameVertexHandle(_mesh);
//
//    std::string output = _OutputFile.substr(_OutputFile.find_last_of('/') + 1);
//
//    auto allf = _mesh.allfaces();
//    auto allv = _mesh.allvertices();
//
//    STEPControl_Writer writer;
//    BRep_Builder builder;
//    const Standard_Real tol = 1.0e-6;
//
//    std::vector<TopoDS_Vertex> vertexList(allv.size());
//    for (const auto& vertex : allv) {
//        const auto& vs = vertex.second;
//        TopoDS_Vertex v;
//        builder.MakeVertex(v, gp_Pnt(vs.x(), vs.y(), vs.z()), tol);
//        writer.Transfer(v, STEPControl_AsIs);
//        vertexList[vertex.first] = v;
//    }
//
//    for (const auto& face : allf) {
//        BRepBuilderAPI_MakePolygon polygonBuilder;
//        for (const auto& vertexHandle : face.second.getVertexHandle()) {
//            polygonBuilder.Add(vertexList[vertexHandle]);
//        }
//        if (polygonBuilder.IsDone()) {
//            TopoDS_Wire wire = polygonBuilder.Wire();
//            if (!wire.IsNull()) {
//                BRepBuilderAPI_MakeFace makeFace(wire, true);
//                writer.Transfer(makeFace.Face(), STEPControl_AsIs);
//            }
//        }
//    }
//
//    Standard_Boolean status = writer.Write(_OutputFile.c_str());
//    if (!status) {
//        std::cout << "Failed to write STEP file!" << std::endl;
//        outputfile.close();
//        return false;
//    } else {
//        std::cout << "STEP file successfully written!" << std::endl;
//    }
//
//    outputfile.close();
//    return true;
//}

bool MeshKernel::IO::WriteStepFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    if (!outputfile.is_open()) {
        std::cerr << "Cannot open output file!" << std::endl;
        return false;
    }

    if (surfacesDraw.empty()) {
        std::cerr << "Error: surfacesDraw is empty!" << std::endl;
        return false;
    }

    double umin, umax, vmin, vmax;
    vector<CPoint> singlePatch; // 每个面片的离散点
    vector<vector<CPoint>> singlelSurface; // 每个曲面的离散点
    vector<vector<vector<CPoint>>> allSurface; // 所有曲面集合的离散点

    vector<CVector> singlePatch_normal;
    vector<vector<CVector>> singlelSurface_normal;
    vector<vector<vector<CVector>>> allSurface_normal;
    CVector normal; // 保存法向量

    STEPControl_Writer writer;

    for (size_t i = 0; i < surfacesDraw.size(); i++) {
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
        int unumber = 4; // u方向的控制顶点数
        int vnumber = 4; // v方向的控制顶点数

        double sizeu = (umax - umin) / (unumber - 1); // 步距
        double sizev = (vmax - vmin) / (vnumber - 1);

        // 构建样条表面点阵
        TColgp_Array2OfPnt points(1, unumber, 1, vnumber);
        for (int u_index = 0; u_index < unumber; ++u_index) {
            for (int v_index = 0; v_index < vnumber; ++v_index) {
                double u = umin + u_index * sizeu;
                double v = vmin + v_index * sizev;
                auto pnt = surfacesDraw[i].getPointAtParam(u, v);
                // 检查 pnt 是否有效
//                if (pnt == nullptr) {
//                    std::cerr << "Error: getPointAtParam returned null!" << std::endl;
//                    return false;
//                }
                points.SetValue(u_index + 1, v_index + 1, gp_Pnt(pnt.getX(), pnt.getY(), pnt.getZ()));
            }
        }

        // 创建样条表面
        GeomAPI_PointsToBSplineSurface splineBuilder(points, 3, 3, GeomAbs_C2); // 注意维度
        Handle(Geom_BSplineSurface) splineSurface = splineBuilder.Surface();

        if (splineSurface.IsNull()) {
            std::cerr << "Failed to create BSpline surface!" << std::endl;
            return false;
        }

        // 创建面
        TopoDS_Face face = BRepBuilderAPI_MakeFace(splineSurface, 1.0e-6);
        if (face.IsNull()) {
            std::cerr << "Failed to create face from BSpline surface!" << std::endl;
            return false;
        }

        writer.Transfer(face, STEPControl_AsIs);
    }

    // 写入STEP文件
    Standard_Boolean status = writer.Write(_OutputFile.c_str());
    if (!status) {
        std::cerr << "Failed to write STEP file!" << std::endl;
        outputfile.close();
        return false;
    } else {
        std::cout << "STEP file successfully written!" << std::endl;
    }

    outputfile.close();
    return true;
}

bool MeshKernel::IO::WriteIgesFileWithSpline(std::vector<CBSplineSurface>& surfacesDraw, const std::string& _OutputFile) {
    std::ofstream outputfile(_OutputFile, std::ios::out);
    if (!outputfile.is_open()) {
        std::cerr << "Cannot open output file!" << std::endl;
        return false;
    }

    if (surfacesDraw.empty()) {
        std::cerr << "Error: surfacesDraw is empty!" << std::endl;
        return false;
    }

    double umin, umax, vmin, vmax;

    IGESControl_Writer writer;
    BRep_Builder builder;

    for (size_t i = 0; i < surfacesDraw.size(); i++) {
        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
        int unumber = 4; // u方向的控制顶点数
        int vnumber = 4; // v方向的控制顶点数

        double sizeu = (umax - umin) / (unumber - 1); // 步距
        double sizev = (vmax - vmin) / (vnumber - 1);
        // 构建样条表面点阵
        TColgp_Array2OfPnt points(1, unumber, 1, vnumber);
        for (int u_index = 0; u_index < unumber; ++u_index) {
            for (int v_index = 0; v_index < vnumber; ++v_index) {
                double u = umin + u_index * sizeu;
                double v = vmin + v_index * sizev;
                auto pnt = surfacesDraw[i].getPointAtParam(u, v);
                points.SetValue(u_index + 1, v_index + 1, gp_Pnt(pnt.getX(), pnt.getY(), pnt.getZ()));
            }
        }

        // 创建样条表面
        GeomAPI_PointsToBSplineSurface splineBuilder(points, 3, 3, GeomAbs_C2);
        Handle(Geom_BSplineSurface) splineSurface = splineBuilder.Surface();

        if (splineSurface.IsNull()) {
            std::cerr << "Failed to create BSpline surface!" << std::endl;
            return false;
        }

        // 创建面
        TopoDS_Face face = BRepBuilderAPI_MakeFace(splineSurface, 1.0e-6);
        if (face.IsNull()) {
            std::cerr << "Failed to create face from BSpline surface!" << std::endl;
            return false;
        }

        writer.AddShape(face);
    }

    // 写入IGES文件
    if (writer.Write(_OutputFile.c_str()) != IFSelect_RetDone) {
        std::cerr << "Failed to write IGES file!" << std::endl;
        outputfile.close();
        return false;
    } else {
        std::cout << "IGES file successfully written!" << std::endl;
    }

    outputfile.close();
    return true;
}

//bool MeshKernel::IO::WriteStepFileWithSpline(vector<CBSplineSurface>& surfacesDraw, int upatchnum, int vpatchnum, const std::string& _OutputFile) {
//    std::ofstream outputfile(_OutputFile, std::ios::out);
//    if (!outputfile.is_open()) {
//        std::cerr << "Cannot open output file!" << std::endl;
//        return false;
//    }
//
//    if (surfacesDraw.empty()) {
//        std::cerr << "Error: surfacesDraw is empty!" << std::endl;
//        return false;
//    }
//
//    // 检查 upatchnum 和 vpatchnum 是否为正数
//    if (upatchnum <= 0 || vpatchnum <= 0) {
//        std::cerr << "Error: upatchnum and vpatchnum must be positive integers!" << std::endl;
//        return false;
//    }
//
//    double umin, umax, vmin, vmax;
//    vector<CPoint> singlePatch; // 每个面片的离散点
//    vector<vector<CPoint>> singlelSurface; // 每个曲面的离散点
//    vector<vector<vector<CPoint>>> allSurface; // 所有曲面集合的离散点
//
//    vector<CVector> singlePatch_normal;
//    vector<vector<CVector>> singlelSurface_normal;
//    vector<vector<vector<CVector>>> allSurface_normal;
//    CVector normal; // 保存法向量
//
//    int totalVerticesSize = 0;
//
//    STEPControl_Writer writer;
//    BRep_Builder builder;
//
//    for (int i = 0; i < surfacesDraw.size(); i++) {
//        surfacesDraw[i].getUVMinMax(umin, umax, vmin, vmax);
//        int unumber, vnumber;
//        unumber = upatchnum * (surfacesDraw[i].getSizeUControlPoints());
//        vnumber = vpatchnum * (surfacesDraw[i].getSizeVControlPoints());
//        unumber = unumber < 100 ? unumber : 100;
//        vnumber = vnumber < 100 ? vnumber : 100;
//        // For faster rendering, reduce the number of patches
//        unumber *= 0.2;
//        vnumber *= 0.2f;
//        double sizeu = (umax - umin) / unumber - 1e-9; // 步距
//        double sizev = (vmax - vmin) / vnumber - 1e-9;
//
//        // 构建样条表面点阵
//        TColgp_Array2OfPnt points(1, unumber, 1, vnumber);
//        for (int i = 0; i < unumber; ++i) {
//            for (int j = 0; j < vnumber; ++j) {
//                double u = umin + i * sizeu;
//                double v = vmin + j * sizev;
//                auto pnt = surfacesDraw[i].getPointAtParam(u, v);
//                points.SetValue(i + 1, j + 1, gp_Pnt(pnt.getX(), pnt.getY(), pnt.getZ()));
//            }
//        }
//
//        // 创建样条表面
//        GeomAPI_PointsToBSplineSurface splineBuilder(points, 4, 4, GeomAbs_C2);
//        Handle(Geom_BSplineSurface) splineSurface = splineBuilder.Surface();
//
//        if (splineSurface.IsNull()) {
//            std::cerr << "Failed to create BSpline surface!" << std::endl;
//            return false;
//        }
//
//        // 创建面
//        TopoDS_Face face = BRepBuilderAPI_MakeFace(splineSurface, 1.0e-6);
//        if (face.IsNull()) {
//            std::cerr << "Failed to create face from BSpline surface!" << std::endl;
//            return false;
//        }
//
//        writer.Transfer(face, STEPControl_AsIs);
//    }
//
//    // 写入STEP文件
//    Standard_Boolean status = writer.Write(_OutputFile.c_str());
//    if (!status) {
//        std::cerr << "Failed to write STEP file!" << std::endl;
//        outputfile.close();
//        return false;
//    } else {
//        std::cout << "STEP file successfully written!" << std::endl;
//    }
//
//    outputfile.close();
//    return true;
//}

std::vector<std::string> MeshKernel::IO::SplitFileName(const std::string& fileName)
{
    // JFR DO NOT CHANGE TO std::vector<std::string> s(3), it segfaults while
    // destructor si called
    std::vector<std::string> s;
    s.resize(3);
    if (fileName.size()) {
        // returns [path, baseName, extension]
        int idot = (int)fileName.find_last_of('.');
        int islash = (int)fileName.find_last_of("/\\");
        if (idot == (int)std::string::npos) idot = -1;
        if (islash == (int)std::string::npos) islash = -1;
        if (idot > 0) s[2] = fileName.substr(idot);
        if (islash > 0) s[0] = fileName.substr(0, islash + 1);
        s[1] =
                fileName.substr(s[0].size(), fileName.size() - s[0].size() - s[2].size());
    }
    return s;
}

bool MeshKernel::IO::WriteVtkFile(const SurfaceMesh& _mesh, const std::string& _OutputFile) {
    // 还不支持约束线的输出

    std::vector<std::array<double, 3>>verts;
    std::vector<std::vector<uint32_t>> faces;
    std::vector<std::array<int, 2>> constrains;

    for (int i = 0; i < _mesh.allvertices().size(); i++) {
        std::array<double, 3> vert;
        vert[0] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).x();
        vert[1] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).y();
        vert[2] = _mesh.vertices(MeshKernel::iGameVertexHandle(i)).z();
        verts.push_back(vert);
    }

    for (auto fp : _mesh.allfaces()) {
        auto fh = fp.first;
        auto face = fp.second;
        std::vector<uint32_t> oneface;
        for (int i = 0; i < face.getVertexHandle().size(); i++) {
            oneface.push_back(face.getVertexHandle()[i]);
        }
        faces.push_back(oneface);
    }

    int numVertices = verts.size();

    FILE* fp = fopen(_OutputFile.c_str(), "w");

    std::vector<std::string> s = SplitFileName(_OutputFile);
    s.resize(3);

    if (!fp) {
        fprintf(stdout, "Unable to open file '%s'\n", _OutputFile.c_str());
        return false;
    }

    fprintf(fp, "# vtk DataFile Version 2.0\n");
    fprintf(fp, "%s, Created by labb \n", s[1].c_str());
    fprintf(fp, "ASCII\n");
    fprintf(fp, "DATASET UNSTRUCTURED_GRID\n");

    // write mesh vertices
    fprintf(fp, "POINTS %d double\n", numVertices);
    for (std::size_t i = 0; i < verts.size(); i++)
        fprintf(fp, "%.16g %.16g %.16g\n", verts[i][0], verts[i][1], verts[i][2]);
    fprintf(fp, "\n");

    int numElements = faces.size();
    //int totalNumInt = numElements * 5;
    int totalNumInt = 0;
    for (std::size_t i = 0; i < faces.size(); i++) {
        if (faces[i].size() == 3) totalNumInt += 4;
        else if (faces[i].size() == 4) totalNumInt += 5;
    }

    for (std::size_t i = 0; i < constrains.size(); ++i) {
        ++numElements;
        totalNumInt += 3;
    }

    // print vertex indices in ascii or binary
    fprintf(fp, "CELLS %d %d\n", numElements, totalNumInt);
    for (std::size_t i = 0; i < constrains.size(); ++i) {
        fprintf(fp, "%d", 2);
        fprintf(fp, " %ld", constrains[i][0]);
        fprintf(fp, " %ld", constrains[i][1]);
        fprintf(fp, "\n");
    }

    for (std::size_t i = 0; i < faces.size(); i++) {
        fprintf(fp, "%d", faces[i].size());
        for (int j = 0; j < faces[i].size(); j++)
            fprintf(fp, " %ld", faces[i][j]);
        fprintf(fp, "\n");
    }
    fprintf(fp, "\n");

    // print element types in ascii or binary
    fprintf(fp, "CELL_TYPES %d\n", numElements);
    for (std::size_t i = 0; i < numElements - faces.size(); i++) {
        fprintf(fp, "%d\n", 3);
    }
    for (std::size_t i = 0; i < faces.size(); i++) {
        if (faces[i].size() == 3) fprintf(fp, "%d\n", 5);
        else if (faces[i].size() == 4) fprintf(fp, "%d\n", 9);
    }
    fclose(fp);
    return true;
}

void MeshKernel::IO::ReOrderiGameVertexHandle(const SurfaceMesh& _mesh) {
    auto allv = _mesh.allvertices();
    int idx = 0;
    for (auto v : allv) {
        reorderedvh_.push_back(v.first);
        newvh_[v.first] = idx++;
    }
}


//MeshKernel::TetMesh MeshKernel::IO::ReadMeshFileFromStr(const std::string& data) {
//	std::stringstream ss(data);
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> surface_faces;
//	std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	std::string str;
//	enum State
//	{
//		USELESS = 1, VERTEX, FACE, TET, EDGE
//	}state = USELESS;
//	while (std::getline(ss, str)) {
//		int len = str.length();
//		std::vector<std::string>info;
//		std::string s;
//		for (int i = 0; i < len; i++) {
//			if (str[i] == '#')break;
//			if (str[i] == ' ') {
//				if (s.length() > 0)
//					info.push_back(s);
//				s = "";
//			}             else {
//				s.push_back(str[i]);
//			}
//		}
//		if (s.length() > 0)
//			info.push_back(s);
//		if (info.size() == 0)continue;
//		if (info[0] == "MeshVersionFormatted") {
//			state = USELESS;
//		}         else if (info[0] == "Dimension") {
//			std::getline(ss, str);
//			state = USELESS;
//		}         else if (info[0] == "Vertices") {
//			std::getline(ss, str);
//			state = VERTEX;
//		}         else if (info[0] == "Tetrahedra") {
//			std::getline(ss, str);
//			state = TET;
//		}         else if (info[0] == "Triangles") {
//			std::getline(ss, str);
//			state = FACE;
//		}         else if (info[0] == "iGameEdges") {
//			std::getline(ss, str);
//			state = EDGE;
//		}         else if (info[0] == "End") {
//			state = USELESS;
//		}         else {
//			if (state == USELESS) {
//				continue;
//			}             else if (state == VERTEX) {
//				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
//					, std::stod(info[1]), std::stod(info[2])));
//			}             else if (state == TET) {
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				cells.push_back(faces);
//				faces.clear();
//			}             else if (state == FACE) {
//				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//			}
//		}
//	}
//	return TetMesh(vertices, cells, surface_faces);
//}
//MeshKernel::TetMesh MeshKernel::IO::ReadTetMeshFile(const std::string& _InputFile) {
//	FILE* fp = fopen(_InputFile.c_str(), "r");
//	char str[100];
//	enum State
//	{
//		USELESS = 1, VERTEX, FACE, TET, EDGE
//	}state = USELESS;
//	std::vector<iGameVertex> vertices;
//	std::vector<std::vector<iGameVertexHandle>> surface_faces;
//	std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
//	std::vector<std::vector<iGameVertexHandle>> faces;
//	while (fscanf(fp, "%[^\n]\n", str) != EOF) {
//		int len = strlen(str);
//		std::vector<std::string>info;
//		std::string s;
//		for (int i = 0; i < len; i++) {
//			if (str[i] == '#')break;
//			if (str[i] == ' ') {
//				if (s.length() > 0)
//					info.push_back(s);
//				s = "";
//			} 			else {
//				s.push_back(str[i]);
//			}
//		}
//		if (s.length() > 0)
//			info.push_back(s);
//		if (info.size() == 0)continue;
//		if (info[0] == "MeshVersionFormatted") {
//			state = USELESS;
//		} 		else if (info[0] == "Dimension") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = USELESS;
//		} 		else if (info[0] == "Vertices") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = VERTEX;
//		} 		else if (info[0] == "Tetrahedra") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = TET;
//		} 		else if (info[0] == "Triangles") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = FACE;
//		} 		else if (info[0] == "iGameEdges") {
//			fscanf(fp, "%[^\n]\n", str);
//			state = EDGE;
//		} 		else if (info[0] == "End") {
//			state = USELESS;
//		} 		else {
//			if (state == USELESS) {
//				continue;
//			} 			else if (state == VERTEX) {
//				vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
//					, std::stod(info[1]), std::stod(info[2])));
//			} 			else if (state == TET) {
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//				});
//				faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[3]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//				cells.push_back(faces);
//				faces.clear();
//			} 			else if (state == FACE) {
//				surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
//					(iGameVertexHandle)(std::stoi(info[0]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[1]) - 1)
//						, (iGameVertexHandle)(std::stoi(info[2]) - 1)
//				});
//			}
//		}
//	}
//	std::cout << "read file success" << std::endl;
//	return TetMesh(vertices, cells, surface_faces);
//}

void MeshKernel::IO::iGameWriteVolumeFile(const VolumeMesh& _mesh, const std::string& _filename) {


    std::ofstream off(_filename.c_str(), std::ios::out);

    if (!off.good()) {
        std::cerr << "Error: Could not open file " << _filename << " for writing!" << std::endl;
        off.close();
        return;
    }

    // Write header
    off << "MeshVersionFormatted 1" << std::endl;
    off << "Dimension 3" << std::endl;
    uint64_t n_vertices(_mesh.vsize());
    off << "Vertices" << std::endl;
    off << n_vertices << std::endl;

    // write vertices
    for (uint64_t v_it = 0; v_it < n_vertices; ++v_it) {
        auto& v = _mesh.vertices((MeshKernel::iGameVertexHandle)v_it);
        off <<  v.x() << " "
            <<  v.y() << " "
            <<  v.z() << " "
            << "-1" << std::endl;
    }
    uint64_t n_faces(_mesh.fsize());
    /*off << "Quadrilaterals" << std::endl;
    off << n_faces << std::endl;

    for (auto& fp : _mesh.allfaces()) {
        auto vhs = fp.second.getVertexHandle();
        for (auto& vh : vhs) {
            off << (int)vh + 1 << " ";
        }
        off << "1" << std::endl;
    }*/

    uint64_t n_cells(_mesh.csize());
    off << "Hexahedra" << std::endl;
    off << n_cells << std::endl;

    for (auto& cp : _mesh.allcells()) {
        auto vhs = cp.second.getVertexHandle();
        for (auto& vh : vhs) {
            off << (int)vh + 1 << " ";
        }
        off << "1" << std::endl;
    }

    off << "End" << std::endl;
    off.close();

}

MeshKernel::VolumeMesh MeshKernel::IO::iGameReadVolumeFile_TopOptMeshGeneration(const std::string& _InputFile) {
    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> cells;
    std::string line;
    std::cout << "ok" << std::endl;
    VolumeMesh mesh;
    int cnt = 0;
    std::cout << cnt << std::endl;
    while (inputfile) {
        std::cout << cnt++ << std::endl;
        line.clear();
        getline(inputfile, line);
        std::stringstream linestream;
        linestream.str(line);
        std::vector<iGameVertexHandle> vhs;
        for (int i = 0; i < 8; ++i) {
            double x, y, z;
            linestream >> x >> y >> z;
            iGameVertex v(x, y, z);
            auto vh = mesh.AddVertex(iGameVertex(v.x(), v.y(), v.z()));
            vhs.push_back(vh);
        }
        mesh.AddCell(vhs);
    }

    return mesh;
}

MeshKernel::VolumeMesh MeshKernel::IO::iGameReadVolumeFile(const std::string& _InputFile) {

    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> cells;
    std::string line;
    //std::cout << "Is Reading : " << _InputFile << " File." << std::endl;

    std::vector<int> invalid_vhs;

    while (inputfile) {
        line.clear();
        getline(inputfile, line);
        //std::cout << line << "\n";
        if (line[0] == '#') continue;

        else if (line[0] == 'V') {

            std::istringstream isstream(line);
            size_t num = 0;
            std::string Vertices;
            isstream >> Vertices;
            assert(Vertices == "Vertices");

            if(!isstream.eof())
            {
                isstream >> num;
                if(num == 0)
                {
                    line.clear();
                    getline(inputfile, line);
                    num = stoi(line);
                }
            }
            else
            {
                line.clear();
                getline(inputfile, line);
                num = stoi(line);
            }


            for (int i = 0; i < num; ++i) {

                line.clear();
                getline(inputfile, line);
                std::stringstream linestream;
                linestream.str(line);
                double x, y, z;
                int tag = -1;
                linestream >> x >> y >> z >> tag;
                iGameVertex vv(x, y, z);
                //printf("v: %.4f, %.4f, %.4f\n", x, y, z);
                vertices.push_back(vv);
                if (tag == 0)
                    invalid_vhs.push_back(vertices.size() - 1);
            }
        }

        else if (line[0] == 'H') {

            std::istringstream isstream(line);

            size_t num = 0;
            std::string Hexahedra;
            isstream >> Hexahedra;
            assert(Hexahedra == "Hexahedra");
            if(!isstream.eof())
            {
                isstream >> num;
                if(num == 0)
                {
                    line.clear();
                    getline(inputfile, line);
                    num = stoi(line);
                }
            }
            else
            {
                line.clear();
                getline(inputfile, line);
                num = stoi(line);
            }


            std::stringstream linestream;

            for (int i = 0; i < num; i++) {

                line.clear();
                getline(inputfile, line);
                linestream.str(line);
                std::vector<iGameVertexHandle> cellH;
                int n = 0;
                for (auto& c : line) if (c == ' ') n++;
                while (n--) {
                    int t;
                    linestream >> t;
                    cellH.push_back(iGameVertexHandle(t - 1));// ������ʱ����Ҫ - 1
                }
                cells.push_back(cellH);
            }
        }
    }
    inputfile.close();

    auto volume_mesh = VolumeMesh(vertices, cells);

    return volume_mesh;
}

MeshKernel::TetMesh MeshKernel::IO::iGameReadTetMeshFile(const std::string& _InputFile) {

    std::ifstream inputfile(_InputFile, std::ios::in);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> cells;
    std::string line;
    std::cout << "Is Reading : " << _InputFile << " File." << std::endl;

    std::vector<int> invalid_vhs;

    while (inputfile) {
        line.clear();
        getline(inputfile, line);
        //std::cout << line << "\n";
        if (line[0] == '#') continue;

        else if (line[0] == 'V') {
            line.clear();
            getline(inputfile, line);
            int num = stoi(line);
            std::cout << "!: " << num << std::endl;
            for (int i = 0; i < num; ++i) {

                line.clear();
                getline(inputfile, line);
                std::stringstream linestream;
                linestream.str(line);
                double x, y, z;
                int tag = -1;
                linestream >> x >> y >> z >> tag;
                iGameVertex vv(x, y, z);
                //printf("v: %.4f, %.4f, %.4f\n", x, y, z);
                vertices.push_back(vv);
            }
        }

        else if (line[0] == 'T') {

            line.clear();
            getline(inputfile, line);
            std::stringstream linestream;
            int num = stoi(line);
            std::cout << "!: " << num << std::endl;
            for (int i = 0; i < num; i++) {
                line.clear();
                getline(inputfile, line);
                linestream.str(line);
                std::vector<iGameVertexHandle> cellH;
                int n = 0;
                //for (auto& c : line) if (c == ' ') n++;
                n = 4;
                while (n--) {
                    int t;
                    linestream >> t;
                    cellH.push_back(iGameVertexHandle(t - 1));// ������ʱ����Ҫ - 1
                }
                cells.push_back(cellH);
            }
        }
    }
    inputfile.close();
    auto volume_mesh = TetMesh(vertices, cells);
    std::cout << "the number of vertices is " << volume_mesh.VertexSize() << std::endl;
    std::cout << "the number of cells is " << volume_mesh.CellSize() << std::endl;

    return volume_mesh;

}

MeshKernel::TetMesh MeshKernel::IO::ReadMeshFileFromStr(const std::string& data) {
    std::stringstream ss(data);
    std::vector<iGameVertex> vertices;
    std::vector<std::vector<iGameVertexHandle>> surface_faces;
    std::vector<std::vector<std::vector<iGameVertexHandle>>> cells;
    std::vector<std::vector<iGameVertexHandle>> faces;
    std::string str;
    enum State
    {
        USELESS = 1, VERTEX, FACE, TET, EDGE
    }state = USELESS;
    while (std::getline(ss, str)) {
        int len = str.length();
        std::vector<std::string>info;
        std::string s;
        for (int i = 0; i < len; i++) {
            if (str[i] == '#')break;
            if (str[i] == ' ') {
                if (s.length() > 0)
                    info.push_back(s);
                s = "";
            }             else {
                s.push_back(str[i]);
            }
        }
        if (s.length() > 0)
            info.push_back(s);
        if (info.size() == 0)continue;
        if (info[0] == "MeshVersionFormatted") {
            state = USELESS;
        }         else if (info[0] == "Dimension") {
            std::getline(ss, str);
            state = USELESS;
        }         else if (info[0] == "Vertices") {
            std::getline(ss, str);
            state = VERTEX;
        }         else if (info[0] == "Tetrahedra") {
            std::getline(ss, str);
            state = TET;
        }         else if (info[0] == "Triangles") {
            std::getline(ss, str);
            state = FACE;
        }         else if (info[0] == "Edges") {
            std::getline(ss, str);
            state = EDGE;
        }         else if (info[0] == "End") {
            state = USELESS;
        }         else {
            if (state == USELESS) {
                continue;
            }             else if (state == VERTEX) {
                vertices.push_back(MeshKernel::iGameVertex(std::stod(info[0])
                        , std::stod(info[1]), std::stod(info[2])));
            }             else if (state == TET) {
                faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
                        (iGameVertexHandle)(std::stoi(info[0]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[1]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[2]) - 1)
                });
                faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
                        (iGameVertexHandle)(std::stoi(info[0]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[2]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[3]) - 1)
                });
                faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
                        (iGameVertexHandle)(std::stoi(info[0]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[3]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[1]) - 1)
                });
                faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
                        (iGameVertexHandle)(std::stoi(info[1]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[3]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[2]) - 1)
                });
                cells.push_back(faces);
                faces.clear();
            }             else if (state == FACE) {
                surface_faces.push_back(std::vector<MeshKernel::iGameVertexHandle>{
                        (iGameVertexHandle)(std::stoi(info[0]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[1]) - 1)
                        , (iGameVertexHandle)(std::stoi(info[2]) - 1)
                });
            }
        }
    }
    return TetMesh(vertices, cells, surface_faces);
}

