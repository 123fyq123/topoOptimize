# MMC / SIMP 框架说明

## 说明范围

这份文档只保留当前项目里和 `MMC`、`SIMP` 直接相关的内容，重点回答下面几个问题：

1. 当前项目里的 `SIMP` 是什么实现。
2. 当前项目里的 `MMC` 是什么实现。
3. 两者在数学上分别优化什么变量。
4. 目标函数、刚度矩阵、灵敏度是怎样一步步推出来的。
5. 当前项目里如果要继续改 `MMC` 或 `SIMP`，应该优先看哪些文件。

本文默认讨论的是当前项目中的两条主线：

- `SIMP / 密度法`：`Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu`
- `3D MMC`：`Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu` 中的 `Cuda_TopOpt_MMC_Node(...)`

本文不再展开 `2D skeleton MMC`、UI 其他按钮、离散密度法等旁支内容。

## 一、SIMP 部分

### 1.1 当前项目里的 SIMP 是什么

当前项目里的 `SIMP` 本质上是经典密度法拓扑优化：

- 设计变量是每个单元的密度 `rho_e`
- `rho_e` 取值范围在 `[rho_min, 1]`
- 单元越接近 `1`，说明材料越“实”
- 单元越接近 `rho_min`，说明材料越“空”

在当前工程实现里，这条主线对应：

- [Cuda_TopOpt_Pipeline.cu](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu)
- 入口通常是 `Cuda_TopOpt_Pipeline(...)`

从代码结构上看，它是：

1. 用当前密度场装配刚度矩阵
2. 解平衡方程 `KU=F`
3. 计算柔顺度与灵敏度
4. 做灵敏度滤波
5. 用 `OC` 风格更新密度

### 1.2 SIMP 的基本公式

经典单材料 SIMP 的材料插值写成：

```text
E_e(rho_e) = E_min + rho_e^p (E_0 - E_min)
```

其中：

- `E_0` 是实材料弹性模量
- `E_min` 是避免刚度矩阵奇异的极小模量
- `p` 是惩罚因子，常取 `3`

如果把单位弹性模量下的单元刚度矩阵记为 `k0`，则第 `e` 个单元刚度矩阵为：

```text
k_e = E_e(rho_e) * k0
```

在很多实现里，如果 `E_min` 很小，常直接写成近似形式：

```text
k_e ≈ rho_e^p * k0
```

当前项目里单材料密度法的核心思想就是这一类写法。

### 1.3 从总势能到平衡方程

线弹性结构的离散总势能可写为：

```text
Pi(U) = 1/2 U^T K U - F^T U
```

其中：

- `U` 是全局位移向量
- `K` 是全局刚度矩阵
- `F` 是外载荷向量

对 `U` 取驻值：

```text
dPi/dU = K U - F = 0
```

于是得到平衡方程：

```text
K U = F
```

这就是程序里每轮优化都要解的主方程。

### 1.4 从单元刚度到全局刚度

设共有 `N` 个单元，则全局刚度矩阵可写为：

```text
K(rho) = Σ_e k_e(rho_e)
```

若采用 SIMP 插值，则：

```text
K(rho) = Σ_e rho_e^p k0_e
```

这里的 `k0_e` 表示第 `e` 个单元在单位弹性模量下的刚度矩阵。

程序实现时通常不是显式地写求和公式，而是：

1. 先为每个单元算出缩放系数
2. 用这个缩放系数去乘基础单元刚度
3. 再装配到全局稀疏矩阵中

### 1.5 柔顺度目标函数

当前项目里默认采用柔顺度最小化。柔顺度写为：

```text
C = F^T U
```

在满足 `KU=F` 时，也可写为：

```text
C = U^T K U
```

进一步展开为单元形式：

```text
C = Σ_e U_e^T k_e U_e
```

如果把单元应变能记为：

```text
Ce = U_e^T k0_e U_e
```

则在单材料 SIMP 下有：

```text
C = Σ_e rho_e^p Ce
```

### 1.6 体积分数约束

最常见的约束是体积分数约束：

```text
V(rho) / V0 <= volfrac
```

离散后通常写为：

```text
Σ_e rho_e v_e <= volfrac * Σ_e v_e
```

若各单元体积相同，则可简化为：

```text
1/N Σ_e rho_e <= volfrac
```

也就是单元密度平均值不超过目标体积分数。

### 1.7 柔顺度对密度的灵敏度推导

先写出目标函数：

```text
C(rho) = F^T U(rho)
```

由于 `U` 依赖于 `rho`，严格求导要用伴随法。对于柔顺度最小化问题，最终可得到经典结果：

```text
dC/drho_e = - U^T (dK/drho_e) U
```

因为：

```text
K = Σ_e rho_e^p k0_e
```

所以：

```text
dK/drho_e = p rho_e^(p-1) k0_e
```

代入得：

```text
dC/drho_e = - p rho_e^(p-1) U_e^T k0_e U_e
```

若定义：

```text
Ce = U_e^T k0_e U_e
```

则：

```text
dC/drho_e = - p rho_e^(p-1) Ce
```

这就是 SIMP 最核心的灵敏度公式。

当前项目里你看到类似：

```text
dc = -3 * density^2 * Ce
```

就是在 `p=3` 情况下这个公式的直接实现。

### 1.8 滤波与 OC 更新

为了抑制棋盘格和网格依赖，当前项目对灵敏度做邻域滤波。形式上可以理解为：

```text
dC_filtered(i) =
Σ_j w_ij rho_j dC_j
/
(rho_i Σ_j w_ij)
```

其中：

- `w_ij` 是邻接权重
- `neighbor_info` 保存了局部邻域

优化更新通常采用 `OC` 风格。其典型形式为：

```text
rho_e^new =
clip(
rho_e * sqrt(-dc_e / lambda),
rho_e - move,
rho_e + move
)
```

然后通过二分法求拉格朗日乘子 `lambda`，使体积约束满足。

### 1.9 当前项目里 SIMP 相关最重要的代码位置

- [Cuda_TopOpt_Pipeline.cu](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu)
- [GP_Spline_v0.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/GP_Spline_v0.h)
- [GP_IGA_kernel.cuh](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/GP_IGA_kernel.cuh)

理解顺序建议是：

1. 先看 `Cuda_TopOpt_Pipeline(...)` 的迭代主循环
2. 再看求解器接口如何根据密度装配刚度
3. 最后看灵敏度与 `OC` 更新

## 二、MMC 部分

### 2.1 当前项目里的 MMC 是什么

当前项目里的 `MMC` 是显式参数化拓扑优化方法。

和 SIMP 最大的不同是：

- `SIMP` 直接以单元密度为设计变量
- `MMC` 以几何构件参数为设计变量

也就是说，`MMC` 不是“每个单元一个变量”，而是“每个构件若干几何参数”。

当前工程中这条主线对应：

- [Cuda_TopOpt_Pipeline.cu](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu)
- 入口是 `Cuda_TopOpt_MMC_Node(...)`

### 2.2 MMC 的设计变量

当前项目里的 3D MMC 构件可概括为一个参数化实体，每个构件通常包含 9 个设计变量：

1. 中心位置
   - `x0`
   - `y0`
   - `z0`
2. 三个方向的尺度参数
   - `l1`
   - `l2`
   - `l3`
3. 三个姿态角
   - `alpha`
   - `beta`
   - `gamma`

如果一共有 `m` 个 MMC 构件，则总设计变量维数大致为：

```text
n = 9m
```

### 2.3 MMC 初始组件摆放参数

当前代码里的初始组件摆放主要由下面几组参数控制：

- `xyzRange`
  - 设计域范围，格式是 `{xmin, xmax, ymin, ymax, zmin, zmax}`
- `xInt, yInt, zInt`
  - 初始组件中心点的“半间距参数”
  - 实际中心点步长分别是 `2*xInt`、`2*yInt`、`2*zInt`
- `BaseComponent = {l1, l2, l3, alpha, beta_amp, gamma_amp}`
  - `l1, l2, l3`：组件在局部三轴方向上的尺度参数
  - `alpha`：绕 `x` 轴的初始旋转角
  - `beta_amp`：绕 `y` 轴的初始旋转角幅值
  - `gamma_amp`：绕 `z` 轴的初始旋转角幅值
- `type`
  - 控制每个采样点放几个组件、以及姿态符号组合方式

这里要特别注意：

- `alpha / beta / gamma` 不是“和某个轴的夹角”
- 它们是“绕某个轴旋转多少角度”
- 当前代码里对应关系是：
  - `alpha`：绕 `x` 轴
  - `beta`：绕 `y` 轴
  - `gamma`：绕 `z` 轴

### 2.4 局部拓扑描述函数

每个 MMC 构件都会定义一个局部拓扑描述函数 `phi_i(x)`。

它的几何含义可以理解为：

- `phi_i(x) > 0`：点 `x` 在第 `i` 个构件内部
- `phi_i(x) = 0`：点 `x` 在第 `i` 个构件边界
- `phi_i(x) < 0`：点 `x` 在第 `i` 个构件外部

在当前项目里，这部分核心逻辑对应：

- `calc_Phi_Node(...)`

### 2.5 多构件合成：从局部场到全局场

理想的集合并操作是：

```text
Phi(x) = max_i phi_i(x)
```

但 `max` 不可导，因此工程里通常采用 `KS` 聚合：

```text
Phi(x) = 1/lambda * log(Σ_i exp(lambda * phi_i(x)))
```

当 `lambda` 足够大时：

```text
Phi(x) ≈ max_i phi_i(x)
```

### 2.6 从全局场到材料指示：Heaviside 投影

得到平滑全局场 `Phi(x)` 后，再通过 Heaviside 投影得到更接近 `0-1` 的材料指示：

```text
H(Phi)
```

当前项目里的意义是：

- `H` 接近 `1` 表示实体
- `H` 接近 `alpha` 或接近 `0` 表示空域

### 2.7 从顶点场到单元密度

MMC 经过 Heaviside 投影后，先得到的是顶点值，而有限元装配需要的是单元级量。

对于八节点六面体单元，最简单的做法是：

```text
rho_e^M = 1/8 Σ_(a=1)^8 H_a
```

这个 `rho_e^M` 就是 MMC 阶段输出的单元拓扑密度。

### 2.8 MMC 下的刚度矩阵

如果采用 MMC 单独优化，当前项目更接近下面这种线性缩放逻辑：

```text
k_e = rho_e^M * k0_e
```

因此全局刚度矩阵可写为：

```text
K = Σ_e rho_e^M k0_e
```

再由：

```text
K U = F
```

求位移。

### 2.9 MMC 的目标函数与体积约束

MMC 单独优化时，目标函数仍然通常取柔顺度最小：

```text
C = F^T U
```

体积约束则写成：

```text
V / V0 <= target_rate
```

离散后，可理解为：

```text
1/N Σ_e rho_e^M <= target_rate
```

### 2.10 MMC 灵敏度的链式求导

若把某个 MMC 设计变量记为 `x_j`，则柔顺度对它的导数要经过：

```text
x_j
-> phi_i
-> Phi
-> H
-> rho_e^M
-> K
-> C
```

所以：

```text
dC/dx_j
= Σ_e (dC/dK) (dK/drho_e^M) (drho_e^M/dH) (dH/dPhi) (dPhi/dphi_i) (dphi_i/dx_j)
```

更常见的写法是用伴随法把前半段压缩，得到：

```text
dC/dx_j
= - U^T (dK/dx_j) U
```

而：

```text
dK/dx_j
= Σ_e k0_e * d(rho_e^M)/dx_j
```

再进一步展开：

```text
d(rho_e^M)/dx_j
= (drho_e^M/dH)
 * (dH/dPhi)
 * (dPhi/dphi_i)
 * (dphi_i/dx_j)
```

### 2.11 MMA 更新

SIMP 常用 `OC` 更新，而 MMC 因为设计变量是几何参数，且上下界、尺度、耦合更复杂，所以更适合 `MMA`。

当前项目中，这部分核心是：

- [MMASolver.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/MMASolver.h)

MMC 主循环里会把全部构件参数展平为一个大向量，再交给 `MMA` 更新。

### 2.12 当前项目里 MMC 相关最重要的代码位置

- [Cuda_TopOpt_Pipeline.cu](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu)
- [MMASolver.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/MMASolver.h)
- [GP_Spline_v0.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/GP_Spline_v0.h)

理解顺序建议是：

1. 先看 `Cuda_TopOpt_MMC_Node(...)` 主循环
2. 再看 `calc_Phi_Node(...)` 怎样计算场函数和导数
3. 再看 `KS + Heaviside + 单元平均`
4. 最后看 `MMA` 如何更新构件参数

## 三、MMC 与 SIMP 的本质区别

### 3.1 设计变量不同

`SIMP`：

```text
设计变量 = 每个单元的密度 rho_e
```

`MMC`：

```text
设计变量 = 每个构件的几何参数 x_j
```

### 3.2 中间变量不同

`SIMP` 直接得到单元密度。

`MMC` 则需要经过：

```text
构件参数 -> 场函数 -> KS 聚合 -> Heaviside -> 单元密度
```

### 3.3 优化器不同

`SIMP`：

- 更适合 `OC`
- 变量多但形式简单

`MMC`：

- 更适合 `MMA`
- 变量少但几何非线性更强

### 3.4 几何表达方式不同

`SIMP`：

- 几何是隐含在密度场里的

`MMC`：

- 几何由显式构件直接表达

## 四、当前项目里继续改 MMC / SIMP 时建议先看的文件

如果后续继续修改这两部分，建议优先回看：

1. [AGENTS.md](D:/labCode/LabViewer-Vsub/AGENTS.md)
2. [mmc-simp-framework-notes.md](D:/labCode/LabViewer-Vsub/docs/mmc-simp-framework-notes.md)
3. [Cuda_TopOpt_Pipeline.cu](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.cu)
4. [Cuda_TopOpt_Pipeline.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/Cuda_TopOpt_Pipeline.h)
5. [GP_Spline_v0.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/GP_Spline_v0.h)
6. [GP_IGA_kernel.cuh](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/GP_IGA_kernel.cuh)
7. [MMASolver.h](D:/labCode/LabViewer-Vsub/Apps/SolverGPU/MMASolver.h)

如果只改 `SIMP`，先看：

- `Cuda_TopOpt_Pipeline(...)`
- `Solver_With_Density(...)`
- 灵敏度滤波与 `OC` 更新

如果只改 `MMC`，先看：

- `Cuda_TopOpt_MMC_Node(...)`
- `calc_Phi_Node(...)`
- `KS` 聚合、Heaviside 投影、`MMA` 更新

## 五、简短结论

当前项目里的 `SIMP` 和 `MMC` 可以用一句话概括：

- `SIMP`：直接优化单元密度，靠惩罚插值和 `OC` 更新做柔顺度最小化。
- `MMC`：直接优化几何构件参数，靠场函数、KS 聚合、Heaviside 投影和 `MMA` 更新做柔顺度最小化。

二者最终都会落到有限元方程：

```text
K U = F
```

区别主要不在求解器，而在“设计变量如何映射到刚度矩阵”这一步。
