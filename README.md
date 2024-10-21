# **LabViewer** 

<p align="center">
    <img src="./Apps/Poster/Resources/iGameLogo.png" width="80%"><br>
</p>

## **Contents**
- [**About**](#about)
- [**Compilation**](#compilation)
  * [**Dependencies**](#dependencies)
- [**Organization**](#organization)
- [**Replicability**](#replicability)
- [**Bibtex**](#bibtex)

## **About**
iGame实验室CMake实例项目

## **Compilation**

Windows和Linux平台可以成功运行

### **Dependencies**

内置依赖
- [OpenMesh](https://www.graphics.rwth-aachen.de:9000/OpenMesh/OpenMesh) to verify the applications against reference CPU implementation
- [RapidJson](https://github.com/Tencent/rapidjson) to report the results in JSON file(s)
- [spdlog](https://github.com/gabime/spdlog) for logging
- freeglut
- nlopt
- opencascade(occt)
- tetgen
- CGAL
- boost
- Eigen
- glm

需要额外安装的依赖(安装完成后记得重启)
- Qt (>=5.12, <6): 安装链接=https://download.qt.io/archive/qt/5.12/5.12.3/qt-opensource-windows-x86-5.12.3.exe
- intel-mkl (<=2023.3.0): 现在不能支持mkl2024版本，请使用ReferenceDLL/mkl-2023.2.0/ 下的安装包，安装mkl2023版本
- cuda (>= 11.8): 安装链接=https://developer.nvidia.com/cuda-downloads


Note: 
1. need to add `-DQT_CMAKE_PATH="your-qt-file-path"` in the cmake command
to find QT. ` "your-qt-file-path" ` can be like `"C:\\Qt\\Qt5.12.3\\5.12.3\\msvc2017_64\\lib\\cmake"`.

编译操作

```bash
> git clone (LabViewer-url)
> cd LabViewer
> mkdir build 
> cd build 
> cmake ../
```
cmake选项
```cmake
#暂时适配ninja生成器
-G Ninja 
#指定mkl interface
-DMKL_INTERFACE_FULL:STRING=intel_lp64
#指定mkl 链接方式，对于使用Pardiso求解器，是必须的 
-DMKL_LINK:STRING=static
#指定qt的cmake目录
-DQT_CMAKE_PATH="Qt\\14\\5.14.2\\msvc2017_64\\lib\\cmake"
#指定mkl的cmake目录，暂时不支持mkl2024版本
-DMKL_DIR="OneAPI/mkl/2023.x.y/lib/cmake/mkl"
#关闭shared libs生成
-DBUILD_SHARED_LIBS:BOOL=OFF
```
cmake选项 exmaple
```cmake
-G Ninja
-DMKL_INTERFACE_FULL:STRING=intel_lp64
-DMKL_LINK:STRING=static
-DMKL_DIR="xx/OneAPI/mkl/202x.x.y/lib/cmake/mkl"
-DQT_CMAKE_PATH=xx\\Qt\\xx\\x.x.x\\msvc20xx_64\\lib\\cmake 
-DBUILD_SHARED_LIBS:BOOL=OFF
-DCPACK_BINARY_TBZ2:BOOL=OFF
-DCPACK_BINARY_TGZ:BOOL=OFF
-DFREEGLUT_BUILD_DEMOS:BOOL=OFF
-DFREEGLUT_BUILD_SHARED_LIBS:BOOL=ON
-DFREEGLUT_PRINT_ERRORS:BOOL=OFF
-DFREEGLUT_PRINT_WARNINGS:BOOL=OFF
```

手动选择cuda arch
30系显卡，可以只保留 "86"
40系显卡，可以只保留 "89"
如果为了保证exe可以在多平台使用，需要保留全部架构
```cmake
#should set for specific gpu structure
# 3090 - "86"
# 4090 - "89"
set(CMAKE_CUDA_ARCHITECTURES "70" "80" "86" "89")
```

## **Organization**

