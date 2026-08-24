## Why

项目当前已经支持独立的 3D MMC 优化和顺序执行的 MMC 后接多材料 SIMP 流程，但还不能基于同一次有限元状态同步更新 MMC 几何与材料分布。需要通过新增函数和独立模块实现同步优化，使各单元当前材料模量能够影响 MMC 几何灵敏度，同时 MMC 占有率能够影响材料灵敏度；现有 MMC、独立 SIMP 和顺序 MMC+SIMP 的源码流程均不得修改。

## What Changes

- 在新头文件和 CUDA 源文件中新增公开函数 `Cuda_TopOpt_MMC_MultiMaterial_Simultaneous(...)`，独立持有求解入口、迭代状态和收敛控制。
- 在新流程中复用现有 3D MMC 的组件参数化、TDF 导数、多组件并集聚合、Heaviside 投影和几何体积灵敏度逻辑。
- 复用现有双材料插值公式，并通过逐单元等效刚度系数与 MMC 占有率耦合。
- 每次同步迭代只求解一次共享平衡方程，再由同一位移场和单元能量场计算 MMC 与材料两组梯度。
- 将材料设计变量固定在完整网格上，避免 MMC 移动时因硬阈值筛选导致设计变量维数发生变化。
- 为新流程增加独立的收敛信息、诊断、解析灵敏度输出和结果导出。
- 将新模式入口放在现有 MMC 优化界面中，并通过独立分发分支只调用新增函数。
- 不修改独立 MMC、独立 SIMP 和顺序 MMC+SIMP 求解函数的函数体、控制流程、参数语义与原有入口映射。

## Capabilities

### New Capabilities

- `mmc-multimaterial-simultaneous-optimization`：使用共享有限元求解和耦合解析灵敏度，同步优化显式 3D MMC 几何与隐式双材料分布。

### Modified Capabilities

无。

## Impact

- 在 `Apps/SolverGPU` 下增加新的头文件、CUDA 实现和公开函数；`Cuda_TopOpt_Pipeline.*`、`Cuda_TopOpt_MMC_SIMP_Sequential.*` 与 `Cuda_TopOpt_MultiMaterial_SIMP.*` 属于不得修改的算法兼容性边界。
- 复用现有 IGA 有限元装配、MMC 设计变量约定、材料插值参数和 MMA 基础设施。
- 在现有 MMC 优化界面中增加“MMC-多材料同步优化”入口，并在 `Apps/Poster/source/Topopt_GPU_Version.cpp` 等界面/分发文件中进行最小接线；现有 MMC 和顺序优化分支仍调用各自原函数。
- 新增 MMC 占有率、材料变量、插值材料模量和等效刚度系数等输出字段。
- 预计不增加新的第三方依赖。
