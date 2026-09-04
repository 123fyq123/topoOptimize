# AGENTS.md

## 目的

本文件用于给后续进入本仓库的代码代理提供“仓库级长期说明”。

重点不是重复通用编码规范，而是说明：

1. 这个项目里 `MMC` 和 `SIMP` 分别指什么。
2. 后续凡是涉及拓扑优化任务，应该先看哪些文件。
3. 编译、验证、定位问题时有哪些仓库特有注意事项。

## 先看什么

如果后续任务涉及 `MMC`、`SIMP`、顺序 `MMC+SIMP`、材料插值、拓扑优化显示或求解器装配，优先阅读：

1. `AGENTS.md`
2. `docs/mmc-simp-framework-notes.md`
3. `Apps/Poster/source/Topopt_GPU_Version.cpp`
4. `Apps/SolverGPU/Cuda_TopOpt_Pipeline.h`
5. `Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu`
6. `Apps/SolverGPU/GP_Spline_v0.h`
7. `Apps/SolverGPU/GP_IGA_kernel.cuh`
8. `Apps/SolverGPU/MMASolver.h`

如果任务专门涉及顺序多材料 `MMC+SIMP`，还要额外看：

1. `Apps/SolverGPU/Cuda_TopOpt_MMC_SIMP_Sequential.h`
2. `Apps/SolverGPU/Cuda_TopOpt_MMC_SIMP_Sequential.cu`
3. `Apps/SolverGPU/Cuda_TopOpt_MultiMaterial_SIMP.h`
4. `Apps/SolverGPU/Cuda_TopOpt_MultiMaterial_SIMP.cu`

## 算法索引

后续如果任务中提到下面这些词，默认按下面方式理解。

### `SIMP`

默认指当前项目中的密度法主线，而不是泛指所有拓扑优化。

主参考：

- `Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu`
- `Cuda_TopOpt_Pipeline(...)`

### `3D MMC`

默认指当前项目 `Apps/SolverGPU` 中的 3D 显式 MMC 拓扑优化主线。

主参考：

- `Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu`
- `Cuda_TopOpt_MMC_Node(...)`
- `Apps/SolverGPU/MMASolver.h`

### `顺序 MMC+SIMP`

默认指当前项目中新加的多材料顺序链路：

1. 先跑 `MMC`
2. 再在 `MMC` 有效设计域内跑材料选择型 `SIMP`

主参考：

- `Apps/SolverGPU/Cuda_TopOpt_MMC_SIMP_Sequential.cu`
- `Apps/SolverGPU/Cuda_TopOpt_MultiMaterial_SIMP.cu`

## UI 与入口

用户从界面点击按钮后，和拓扑优化最相关的入口主要在：

- `Apps/Poster/source/iGamePoster.cpp`
- `Apps/Poster/source/Topopt_GPU_Version.cpp`

其中：

- `Topopt_GPU_Version.cpp` 负责网格转换、Bezier 导出、边界条件整理、调用 `SolverGPU`
- 真正的数值迭代主逻辑主要在 `Apps/SolverGPU`

## 编译注意事项

本项目是明显偏 Windows / MSVC / CUDA 的工程。

重要注意事项：

1. 普通 PowerShell 环境下不一定具备完整的 MSVC / SDK / CUDA 编译环境。
2. 真正做编译验证时，优先使用 Visual Studio Developer Command Prompt。
3. 当前仓库常用构建目录是：
   - `out/build/x64-Debug`

常用编译命令模式：

```powershell
cmd /c "call ""D:\apps\VS\Community\Common7\Tools\VsDevCmd.bat"" -arch=x64 && cmake --build out\build\x64-Debug --config Debug"
```

## 验证注意事项

后续如果构建失败，要先区分是哪一类问题：

1. 环境问题
   - 没进 `VsDevCmd`
   - SDK / MSVC / CUDA 路径问题
2. 源码问题
   - 真实编译错误
   - 接口不匹配
   - 类型或语法错误
3. 进程占用问题
   - 尤其是 `iGamePoster.exe` 正在运行时，常见 `LNK1168`

这个区分在本仓库里非常重要，因为很多“编译失败”并不等于代码本身有错。

## 修改建议

后续改动时，默认遵守下面这些仓库级原则：

1. 优先做局部、小范围修改，不做无关重构。
2. 生成目录如 `out/build` 默认不手改。
3. 如果工作区是脏的，不要随意覆盖用户原有修改。
4. 涉及 `MMC / SIMP` 时，先确认当前讨论的是哪条算法分支，再动手。
5. 如果只是为了可视化或导出，不要随便改变主求解逻辑，优先新加函数或新字段。
6. 新增或修改实现代码时，必须为关键算法步骤、公式到代码的映射、状态更新及非直观分支添加简洁中文注释；显而易见的逐行赋值不需要注释。

## 额外说明

仓库中的 `docs/mmc-simp-framework-notes.md` 是当前 `MMC / SIMP` 的详细中文说明。

后续如果任务涉及下面内容，默认先回看该文档：

- `SIMP` 推导
- `MMC` 推导
- `MMC` 与 `SIMP` 的区别
- 顺序 `MMC+SIMP` 的含义
- 刚度矩阵、灵敏度、材料插值的物理解释
