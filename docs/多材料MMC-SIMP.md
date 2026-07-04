

# 《基于 MMC-SIMP 顺序拓扑优化方法的多材料结构设计》论文解读
> **原文标题**：Design of the multi-material structure using an MMC-SIMP sequential topology optimization method  
**作者**：Zhao Li, Hongyu Xu*, Shuai Zhang*, Jintao Cui, Xiaofeng Liu（河南科技大学）  
**期刊**：PLOS ONE, 2025, 20(5): e0321100  
**DOI**：[https://doi.org/10.1371/journal.pone.0321100](https://doi.org/10.1371/journal.pone.0321100)  
**时间线**：收稿 2025-01-22 / 接收 2025-03-01 / 发表 2025-05-09
>

---

## 一、研究背景与动机
### 1.1 多材料结构的意义
相较于单材料结构，多材料结构能够更合理地利用不同性能的材料共同承担载荷，满足多样化设计需求，并进一步提升结构轻量化水平。随着增材制造（3D 打印）的快速发展，多材料结构设计成为研究热点。

### 1.2 单材料 vs 多材料拓扑优化
| 类型 | 核心问题 |
| --- | --- |
| 单材料拓扑优化 | 只判断设计域中各点**有/无**材料 |
| 多材料拓扑优化 | 判断设计域中各点是否有材料、以及**是哪种**材料 |


多材料拓扑优化的核心在于**材料的参数化**，关键是构造"材料插值机制"来表征材料属性指示函数：

+ SIMP 方法 → 密度函数
+ 水平集方法 → 水平集函数
+ MMC 方法 → 组件描述函数

### 1.3 两种优化思路
+ **并行优化（Concurrent）**：结构拓扑与材料选择同时进行，实时交换拓扑信息与材料属性信息。能自然获得更优的多材料结构，但**过程复杂、单次迭代计算规模大**。
+ **顺序优化（Sequential）**：先确定结构拓扑，再在已知拓扑内优化材料分布。**易于实现、收敛可控**，但拓扑由单材料确定，最终结构未必是多材料下的全局最优。

> 本文采用**顺序优化**思路，并通过实验论证：不同杨氏模量得到的 MMC 拓扑非常相似，说明材料差异对拓扑形状影响有限，从而佐证了顺序优化思路的可行性。
>

---

## 二、核心思想：MMC-SIMP 两步顺序法
```plain
┌─────────────┐        组件映射单元         ┌─────────────┐
│  Step 1     │  components mapping        │  Step 2     │
│  MMC 方法   │ ─────  elements  ────────► │  SIMP 方法  │
│ 确定结构拓扑 │  (代理材料模型识别单元)     │ 优化材料布局 │
└─────────────┘                            └─────────────┘
   显式边界                                    内部多材料分配
（力传递路径）                              （基体 E₁ / 增强 E₂）
```

**关键创新点**：通过"组件映射单元"建立 MMC 与 SIMP 的数据衔接。映射单元内部的材料布局变化**不会影响外部结构拓扑**，从而兼得：

+ MMC 的**显式、清晰的外部拓扑边界**
+ SIMP 的**高自由度内部多材料布局**
+ 实现简单、易于工程落地

---

## 三、方法详解
### 3.1 第一步：MMC 确定结构拓扑
MMC（Moving Morphable Component，移动变形组件法）用一组可移动、可变形的几何组件来描述结构。

**组件描述函数**（采用二次变厚度组件）：

$ \phi(x,y) = \left(\frac{x'}{L}\right)^p + \left(\frac{y'}{f(x')}\right)^p - 1 $

其中局部坐标通过旋转变换得到：

$ \begin{Bmatrix} x' \\ y' \end{Bmatrix} =
\begin{bmatrix} \cos\theta & \sin\theta \\ -\sin\theta & \cos\theta \end{bmatrix}
\begin{Bmatrix} x-x_0 \\ y-y_0 \end{Bmatrix} $

厚度变化函数（二次变厚度）：

$ f(x') = \frac{t_1 + t_3 - 2t_2}{2L^2}(x')^2 + \frac{t_3 - t_1}{2L}x' + t_2 $

+ $ L $：组件半长；$ \theta $：组件旋转角；$ (x_0,y_0) $：组件中心
+ $ t_1, t_2, t_3 $：组件三处控制厚度；$ p $ 一般取 6（较大偶数）

### 3.2 衔接：代理材料模型与组件映射单元
MMC 优化完成后，用**代理材料模型（Surrogate Material Model）**识别哪些单元被组件覆盖：

$ \rho^M_e = \left(\frac{\sum_{i=1}^{4} H(\phi^e_i)}{4}\right)^q $

+ $ H(\cdot) $：Heaviside 函数；$ \phi^e_i $ 为单元 $ e $ 四个节点处的拓扑描述函数值
+ $ q $：惩罚因子（通常 $ q=2 $）

正则化的 Heaviside 函数（避免刚度矩阵奇异）：

$ H(x)=
\begin{cases}
1, & x>\varepsilon \\[4pt]
\dfrac{3(1-\alpha)}{4}\left(\dfrac{x}{\varepsilon}-\dfrac{x^3}{3\varepsilon^3}\right)+\dfrac{1+\alpha}{2}, & -\varepsilon \le x \le \varepsilon \\[6pt]
\alpha, & x<-\varepsilon
\end{cases} $

+ $ \varepsilon $：控制光滑过渡宽度的正则化参数
+ $ \alpha $：模拟空洞刚度的小正数，避免整体刚度矩阵奇异

> **判定准则**：单元密度 $ \rho^M_e > \alpha^q $ 的单元被识别为**组件映射单元**，构成 SIMP 的优化设计域。
>

### 3.3 第二步：SIMP 进行材料选择（多材料布局）
> 这是论文重点分析的部分。SIMP 在此**不做传统的"有/无材料"优化**，而是在 MMC 已确定的映射单元内做**两相材料分配**。
>

**两相材料插值模型（无空洞）**：

$ E = \rho^p_e E_1 + (1-\rho^p_e)E_2 $

+ $ \rho_e $：人工密度（设计变量）；$ p $：惩罚幂（通常 $ p \ge 3 $）
+ $ E_1 $：基体材料（base material）
+ $ E_2 $：增强材料（reinforcement material）
+ 插值不含 void，仅在两种实体材料间切换；增强材料布局由 $ \rho_e $ 决定

**单元刚度矩阵**（MMC 代理密度固定，SIMP 密度可变）：

$ k_e = \rho^M_e \cdot \big[\rho^p_e E_1 + (1-\rho^p_e)E_2\big] \cdot k_0 $

### 3.4 优化问题建模
$ \begin{aligned}
\min_{\rho} \quad & C = U^T K U \\
\text{s.t.} \quad & KU = F \\
& V(D) = V_1 + V_2 \le \bar{V} \quad (\text{由 MMC 结果确定}) \\
& 0 < \rho_{\min} \le \rho_e \le 1
\end{aligned} $

引入**材料体积比** $ f = V_2 / V_1 $（增强材料与基体材料体积之比），算例取 $ f = [0.2, 0.5, 0.8, 1.0, 1.5, 2.0] $。

**灵敏度分析**（仅对映射单元计算，其余单元置 0）：

$ \frac{\partial k_e}{\partial \rho_e} = \rho^M_e \cdot \big[p \rho^{p-1}_e (E_1 - E_2)\big] \cdot k_0 $

**数值实现**：

+ 采用 Sigmund 灵敏度过滤（滤波半径 $ r_{\min}=2 $）抑制棋盘格现象
+ 由于是两相材料问题（仅一个设计变量），使用**优化准则法（OC）**更新设计变量

---

## 3.5 关键公式推导
> 本节对优化建模中的核心公式给出逐步推导，便于理解灵敏度与 OC 更新的来源。
>

### 3.5.1 单元刚度对设计变量的导数
由两相材料插值模型与单元刚度矩阵：

$ E(\rho_e) = \rho^p_e E_1 + (1-\rho^p_e)E_2, \qquad
k_e = \rho^M_e \cdot E(\rho_e)\cdot k_0 $

其中 $ \rho^M_e $ 为 MMC 代理密度（顺序优化中**固定不变**），$ k_0 $ 为单位弹性模量下的单元刚度矩阵。对设计变量 $ \rho_e $ 求导：

$ \frac{\partial E}{\partial \rho_e}
= p\rho^{p-1}_e E_1 - p\rho^{p-1}_e E_2
= p\rho^{p-1}_e (E_1 - E_2) $

$ \boxed{\;\frac{\partial k_e}{\partial \rho_e}
= \rho^M_e \cdot p\rho^{p-1}_e (E_1 - E_2)\cdot k_0\;} $

> 说明：因为 $ \rho^M_e $ 固定，求导只作用于 SIMP 插值项。对非映射单元 $ \rho^M_e \le \alpha^q $，其贡献近似为零，故灵敏度只在映射单元上计算。
>

### 3.5.2 柔度（目标函数）灵敏度推导
目标函数为结构柔度：

$ C = U^T K U = F^T U, \qquad K = \sum_e k_e $

由平衡方程 $ KU=F $ 且外载 $ F $ 与设计变量无关（$ \partial F/\partial \rho_e = 0 $），两边对 $ \rho_e $ 求导：

$ \frac{\partial K}{\partial \rho_e}U + K\frac{\partial U}{\partial \rho_e} = 0
\;\Longrightarrow\;
\frac{\partial U}{\partial \rho_e} = -K^{-1}\frac{\partial K}{\partial \rho_e}U $

对 $ C = U^T K U $ 求全导数：

$ \frac{\partial C}{\partial \rho_e}
= \frac{\partial U^T}{\partial \rho_e}KU
+ U^T\frac{\partial K}{\partial \rho_e}U
+ U^T K\frac{\partial U}{\partial \rho_e} $

利用 $ K $ 对称（$ K^T=K $）及 $ KU=F $，第一项与第三项满足：

$ \frac{\partial U^T}{\partial \rho_e}KU
= U^T K\frac{\partial U}{\partial \rho_e} $

代入 $ \dfrac{\partial U}{\partial \rho_e} = -K^{-1}\dfrac{\partial K}{\partial \rho_e}U $：

$ 2\,U^T K\frac{\partial U}{\partial \rho_e}
= -2\,U^T K K^{-1}\frac{\partial K}{\partial \rho_e}U
= -2\,U^T\frac{\partial K}{\partial \rho_e}U $

三项相加：

$ \frac{\partial C}{\partial \rho_e}
= U^T\frac{\partial K}{\partial \rho_e}U - 2\,U^T\frac{\partial K}{\partial \rho_e}U
= -\,U^T\frac{\partial K}{\partial \rho_e}U $

由于 $ \dfrac{\partial K}{\partial \rho_e}=\dfrac{\partial k_e}{\partial \rho_e} $（仅单元 $ e $ 相关），用单元位移 $ u_e $ 表达：

$ \boxed{\;\frac{\partial C}{\partial \rho_e}
= -\,u_e^T\frac{\partial k_e}{\partial \rho_e}u_e
= -\,\rho^M_e\, p\rho^{p-1}_e (E_1 - E_2)\,u_e^T k_0\, u_e\;} $

> 自伴随特性：柔度问题的灵敏度无需额外求解伴随方程，直接用位移结果即可。当 $ E_1<E_2 $（增强材料更硬）时 $ E_1-E_2<0 $，灵敏度为正，OC 会推动 $ \rho_e $ 减小 → 该处倾向选用增强材料 $ E_2 $。
>

### 3.5.3 体积约束灵敏度
体积约束 $ V(\rho)=\sum_e v_e \rho_e \le \bar V $，对设计变量求导：

$ \frac{\partial V}{\partial \rho_e} = v_e $

其中 $ v_e $ 为单元体积（均匀网格下为常数）。材料体积比 $ f=V_2/V_1 $ 用于确定增强材料的总用量约束 $ \bar V $。

### 3.5.4 灵敏度过滤（Sigmund 滤波）
为抑制棋盘格与网格依赖，对灵敏度做加权平均：

$ \widehat{\frac{\partial C}{\partial \rho_e}}
= \frac{1}{\rho_e \sum_{i\in N_e} H_{ei}}
\sum_{i\in N_e} H_{ei}\,\rho_i\,\frac{\partial C}{\partial \rho_i} $

权重因子由滤波半径 $ r_{\min} $ 决定（取邻域 $ N_e $ 内单元中心距离 $ \text{dist}(e,i)\le r_{\min} $）：

$ H_{ei} = \max\big(0,\; r_{\min} - \text{dist}(e,i)\big) $

### 3.5.5 OC（优化准则）更新公式推导
构造拉格朗日函数（柔度最小 + 体积约束）：

$ \mathcal{L} = C + \lambda\Big(\sum_e v_e\rho_e - \bar V\Big) $

KKT 最优条件 $ \partial \mathcal{L}/\partial \rho_e = 0 $ 给出：

$ \frac{\partial C}{\partial \rho_e} + \lambda v_e = 0
\;\Longrightarrow\;
B_e \equiv \frac{-\,\partial C/\partial \rho_e}{\lambda v_e} = 1 $

定义更新因子 $ B_e $，按其偏离 1 的程度迭代（启发式 OC 格式，$ \eta $ 为阻尼系数，通常 $ \eta=0.5 $，$ m $ 为移动限）：

$ \boxed{\;
\rho_e^{\text{new}} =
\begin{cases}
\max(\rho_{\min},\,\rho_e - m), & \rho_e B_e^{\eta} \le \max(\rho_{\min},\,\rho_e - m) \\[4pt]
\rho_e B_e^{\eta}, & \text{otherwise} \\[4pt]
\min(1,\,\rho_e + m), & \rho_e B_e^{\eta} \ge \min(1,\,\rho_e + m)
\end{cases}
\;} $

拉格朗日乘子 $ \lambda $ 通过**二分法**调节，使更新后满足体积约束 $ \sum_e v_e\rho_e^{\text{new}} = \bar V $。

> 由于本方法在两相材料下只有一个密度变量 $ \rho_e $，OC 法即可高效收敛，无需 MMA 等更复杂的优化器。
>

---

## 四、数值算例
论文用三个经典算例验证方法：

| 算例 | 设计域尺寸 | 网格分辨率 | 备注 |
| --- | --- | --- | --- |
| MBB 梁 | 3 × 1 | 90 × 30 | 二次变厚度组件 |
| 短梁（short beam） | 2 × 1 | 80 × 40 | 二次变厚度组件 |
| 薄壁短梁 | 2 × 1 | 120 × 60 | 单孔薄壁组件（缩放因子 a=b=0.5，来自参考文献[32]） |


在最终优化结构中：**黑色区域为基体材料 **$ E_1 $** 的布局，绿色区域为增强材料 **$ E_2 $** 的布局。**

### 4.1 优化效果分析
**(1) 柔度下降率（Flexibility Reduction Rate）**

$ \eta = \frac{C(f_0) - C(f_i)}{C(f_0)} \times 100\% $

其中 $ C(f_0) $ 为 $ f=0 $（单材料）时的柔度，$ C(f_i) $ 为不同材料比下的柔度。

+ 随着 $ f $ 增大，结构柔度持续下降（即刚度提升），但柔度下降率的**增速逐渐放缓**
+ 结论：增强材料越多刚度提升越显著，但应在**合适位置**适量增加增强材料

参考单材料柔度值：MBB 梁 $ c=230.147 $；短梁 $ c=75.4630 $；薄壁短梁 $ c=85.575 $（均为 $ f=0 $）。

**(2) 计算效率与精度**

+ 影响效率的两方面：设计变量规模（灵敏度分析+变量更新）、有限元分析（刚度矩阵+平衡方程求解）
+ 顺序法总耗时 ≈ MMC 耗时 + SIMP 耗时；代理材料模型识别映射单元**不增加额外计算成本**
+ 两相材料下 SIMP 仅一个变量，可用 OC 法高效更新
+ 整体精度取决于 MMC 与 SIMP 各自的固有精度

**(3) 与其他方法对比**  
论文将 MBB 多材料结果与四种主流方法对比：水平集法[33]、交替活跃相算法[34]、新映射插值函数法[35]、材料场级数展开（MFSE）模型[36]。本方法的优势在于：**外部拓扑边界清晰显式、内部多材料布局分明、棋盘格与灰度单元等数值问题更少**。

---

## 五、结论与展望
1. **核心贡献**：提出 MMC-SIMP 顺序拓扑优化方法，通过"组件映射单元"衔接两种方法，兼得显式外部拓扑与高自由度内部多材料布局，且易于实现、提供了 MATLAB 关键代码。
2. **顺序优化的可行性**：实验表明不同杨氏模量得到的 MMC 拓扑高度相似，说明材料差异对拓扑形状影响有限，从而支撑了"先定拓扑、再分材料"的顺序思路。
3. **工程价值**：由于整体有效结构已由 MMC 确定，SIMP 得到的局部增强结构具有更强的可调性，便于对不连续、孤立、难加工区域进行后处理而不影响整体性能。适用于"已有结构件局部替换高强轻质材料以提升刚度、实现轻量化"等场景。
4. **局限**：顺序优化的拓扑通常由单材料确定，最终结构**未必是多材料下的全局最优**（这是顺序法相对并行法的固有代价）。本文仅研究两相材料，可进一步拓展到多相材料与三维结构。

---

## 六、总结
> **MMC 先画出"骨架"（显式力传递路径），SIMP 再在骨架内部"填肉、选材"（决定哪里用基体材料、哪里用增强材料），两步之间用"组件映射单元"无缝衔接，互不干扰。**
>

---

## 七、三维扩展：MMC 的三维公式推导
> 本节专门推导 MMC 部分从二维到三维的扩展。**SIMP 第二步（两相材料插值、柔度灵敏度、OC 更新）的公式与二维完全一致**，仅 $ u_e $ 的自由度由 8（4 节点 × 2）变为 24（8 节点 × 3），$ k_0 $ 改为三维实体单元刚度，故不再重复，下面只推导 MMC 几何描述与其灵敏度。
>

### 7.1 三维组件描述函数（TDF）
二维组件用平面内的变厚度条带描述：

$ \phi(x,y) = \left(\frac{x'}{L}\right)^p + \left(\frac{y'}{f(x')}\right)^p - 1 $

三维组件升级为**沿轴线 **$ x' $** 延伸、在 **$ y' $**-**$ z' $** 截面内变截面的实体**（柱状/锥状组件）：

$ \boxed{\;\phi(x,y,z) = \left(\frac{x'}{L}\right)^p
+ \left(\frac{y'}{f_y(x')}\right)^p
+ \left(\frac{z'}{f_z(x')}\right)^p - 1\;} $

+ $ \phi>0 $：点在组件外部；$ \phi=0 $：组件表面；$ \phi<0 $：组件内部。
+ $ p $ 取较大偶数（如 6），使组件接近"带圆角的长方体/椭球柱"。
+ 若 $ f_y=f_z=r(x') $，组件退化为**回转体（圆截面变半径柱）**，可减少设计变量。

### 7.2 三维局部坐标变换（两个旋转角）
二维只需一个旋转角 $ \theta $；三维需**两个角度**确定组件轴线方向：方位角 $ \varphi $（绕 $ z $ 轴）与仰角 $ \psi $（绕中间 $ y $ 轴）。坐标变换：

$ \begin{Bmatrix} x' \\ y' \\ z' \end{Bmatrix}
= R\,
\begin{Bmatrix} x-x_0 \\ y-y_0 \\ z-z_0 \end{Bmatrix},
\qquad R = R_y(\psi)\,R_z(\varphi) $

$ R_z(\varphi)=\begin{bmatrix}
\cos\varphi & \sin\varphi & 0\\
-\sin\varphi & \cos\varphi & 0\\
0 & 0 & 1
\end{bmatrix},
\qquad
R_y(\psi)=\begin{bmatrix}
\cos\psi & 0 & -\sin\psi\\
0 & 1 & 0\\
\sin\psi & 0 & \cos\psi
\end{bmatrix} $

展开后 $ R $ 的各分量 $ R_{ij} $ 是 $ \varphi,\psi $ 的三角函数组合，后续灵敏度需要 $ \partial R/\partial\varphi $、$ \partial R/\partial\psi $。

### 7.3 三维变截面函数
沿 $ y' $、$ z' $ 两个方向各用一条二次变厚度曲线：

$ f_y(x') = \frac{t^y_1+t^y_3-2t^y_2}{2L^2}(x')^2 + \frac{t^y_3-t^y_1}{2L}\,x' + t^y_2 $

$ f_z(x') = \frac{t^z_1+t^z_3-2t^z_2}{2L^2}(x')^2 + \frac{t^z_3-t^z_1}{2L}\,x' + t^z_2 $

+ $ t^y_1,t^y_2,t^y_3 $：$ y' $ 方向首/中/尾三处半厚；$ t^z_1,t^z_2,t^z_3 $：$ z' $ 方向对应半厚。

### 7.4 三维组件设计变量
| 维度 | 设计变量向量 $ \mathbf{d} $ | 数量 |
| --- | --- | --- |
| 2D | $ (x_0,y_0,L,t_1,t_2,t_3,\theta) $ | 7 |
| 3D | $ (x_0,y_0,z_0,L,t^y_1,t^y_2,t^y_3,t^z_1,t^z_2,t^z_3,\varphi,\psi) $ | 12 |
| 3D 回转体 | $ (x_0,y_0,z_0,L,r_1,r_2,r_3,\varphi,\psi) $ | 9 |


> 单组件变量由 7 增至 12（或 9）。$ N $ 个组件总变量数为 $ 12N $，三维优化的变量规模因此显著上升。
>

### 7.5 全局拓扑描述函数（多组件并集）
整个结构由 $ N $ 个组件取并集，用 $ \max $ 或其光滑近似（KS 聚合函数）：

$ \phi^s(\mathbf{x}) = \max\big(\phi^1,\phi^2,\dots,\phi^N\big)
\;\approx\;
\frac{1}{\mu}\ln\!\left(\sum_{k=1}^{N} e^{\mu\,\phi^k}\right) $

$ \mu $ 为聚合参数（越大越逼近真实 $ \max $）。三维与二维**形式完全相同**，仅 $ \phi^k $ 改为三维描述函数。

### 7.6 三维代理材料模型（8 节点）
二维四边形单元取 4 节点，三维六面体单元取 8 节点：

$ \boxed{\;\rho^M_e = \left(\frac{\sum_{i=1}^{8} H(\phi^{s}(\mathbf{x}^e_i))}{8}\right)^q\;} $

正则化 Heaviside $ H(\cdot) $、判定准则 $ \rho^M_e>\alpha^q $ 与二维完全一致。

### 7.7 MMC 阶段几何灵敏度（链式法则推导）
MMC 第一步本身是单材料最小柔度优化，设计变量是组件几何参数 $ \mathbf{d} $。柔度对几何变量 $ d_j $ 的灵敏度需经**多层链式法则**：

$ \frac{\partial C}{\partial d_j}
= \sum_e \frac{\partial C}{\partial \rho^M_e}
\sum_{i=1}^{8}\frac{\partial \rho^M_e}{\partial H_i}
\frac{\partial H_i}{\partial \phi^s_i}
\frac{\partial \phi^s_i}{\partial d_j} $

**第 1 层——柔度对代理密度**（与二维同构，单材料惩罚 $ k_e=(\rho^M_e)^{p_M}E_0 k_0 $）：

$ \frac{\partial C}{\partial \rho^M_e}
= -\,u_e^T\frac{\partial k_e}{\partial \rho^M_e}u_e
= -\,p_M(\rho^M_e)^{p_M-1}E_0\,u_e^T k_0\,u_e $

**第 2 层——代理密度对节点 Heaviside 值**：

$ \frac{\partial \rho^M_e}{\partial H_i}
= q\left(\frac{\sum_{i} H_i}{8}\right)^{q-1}\cdot\frac{1}{8} $

**第 3 层——Heaviside 对描述函数**（正则化区间 $ -\varepsilon\le\phi\le\varepsilon $ 内非零）：

$ \frac{\partial H}{\partial \phi}
= \frac{3(1-\alpha)}{4\varepsilon}\left(1-\frac{\phi^2}{\varepsilon^2}\right),
\qquad |\phi|\le\varepsilon;\quad \text{否则为 } 0 $

**第 4 层——描述函数对几何变量**（三维相对二维新增工作量所在）。令 $ u=\dfrac{x'}{L},\;v=\dfrac{y'}{f_y},\;w=\dfrac{z'}{f_z} $，则 $ \phi=u^p+v^p+w^p-1 $，由复合求导：

$ \frac{\partial \phi}{\partial d_j}
= p\,u^{p-1}\frac{\partial u}{\partial d_j}
+ p\,v^{p-1}\frac{\partial v}{\partial d_j}
+ p\,w^{p-1}\frac{\partial w}{\partial d_j} $

下面给出几个关键分量（以中心坐标和旋转角为例）。

**(a) 对中心坐标 **$ x_0,y_0,z_0 $：由 $ \{x',y',z'\}^T = R\,\{x-x_0,\,y-y_0,\,z-z_0\}^T $，得

$ \frac{\partial x'}{\partial x_0}=-R_{11},\quad
\frac{\partial y'}{\partial x_0}=-R_{21},\quad
\frac{\partial z'}{\partial x_0}=-R_{31}\;\;(\text{对 }y_0,z_0\text{ 类推取 }R\text{ 第 2、3 列}) $

于是（以 $ x_0 $ 为例，注意 $ f_y,f_z $ 仅依赖 $ x' $）：

$ \frac{\partial u}{\partial x_0}=\frac{1}{L}(-R_{11}),\quad
\frac{\partial v}{\partial x_0}=\frac{-R_{21}}{f_y}-\frac{y'}{f_y^2}f_y'(-R_{11}),\quad
\frac{\partial w}{\partial x_0}=\frac{-R_{31}}{f_z}-\frac{z'}{f_z^2}f_z'(-R_{11}) $

其中 $ f_y'=\dfrac{df_y}{dx'} $。

**(b) 对旋转角 **$ \varphi,\psi $：

$ \frac{\partial}{\partial \varphi}
\begin{Bmatrix} x' \\ y' \\ z' \end{Bmatrix}
= \frac{\partial R}{\partial \varphi}
\begin{Bmatrix} x-x_0 \\ y-y_0 \\ z-z_0 \end{Bmatrix},
\qquad
\frac{\partial}{\partial \psi}
\begin{Bmatrix} x' \\ y' \\ z' \end{Bmatrix}
= \frac{\partial R}{\partial \psi}
\begin{Bmatrix} x-x_0 \\ y-y_0 \\ z-z_0 \end{Bmatrix} $

由 $ R=R_y(\psi)R_z(\varphi) $：

$ \frac{\partial R}{\partial \varphi}=R_y(\psi)\,R_z'(\varphi),
\qquad
\frac{\partial R}{\partial \psi}=R_y'(\psi)\,R_z(\varphi) $

$ R_z'(\varphi)=\begin{bmatrix}
-\sin\varphi & \cos\varphi & 0\\
-\cos\varphi & -\sin\varphi & 0\\
0&0&0
\end{bmatrix},
\quad
R_y'(\psi)=\begin{bmatrix}
-\sin\psi & 0 & -\cos\psi\\
0&0&0\\
\cos\psi & 0 & -\sin\psi
\end{bmatrix} $

得到 $ \partial x'/\partial\varphi $ 等后，再按 $ u,v,w $ 的链式关系代入即可。

**(c) 对长度/厚度 **$ L,t^y_i,t^z_i $：直接对 $ u=x'/L $、$ f_y,f_z $ 求偏导，例如

$ \frac{\partial u}{\partial L}=-\frac{x'}{L^2},\qquad
\frac{\partial v}{\partial t^y_2}=-\frac{y'}{f_y^2}\frac{\partial f_y}{\partial t^y_2},\qquad
\frac{\partial f_y}{\partial t^y_2}=-\frac{(x')^2}{L^2}+1 $

（其余 $ t^y_1,t^y_3,t^z_i $ 按 $ f_y,f_z $ 的二次式系数逐项求导即可。）

### 7.8 二维 → 三维差异小结
| 项目 | 2D | 3D | 改动量 |
| --- | --- | --- | --- |
| 描述函数 $ \phi $ | $ u^p+v^p-1 $ | $ u^p+v^p+w^p-1 $ | 加一项 $ w^p $ |
| 坐标变换 | 1 个角 $ \theta $ | 2 个角 $ \varphi,\psi $，$ R=R_yR_z $ | 旋转矩阵升维 |
| 变截面 | 1 条 $ f(x') $ | 2 条 $ f_y,f_z $ | 翻倍 |
| 单组件变量 | 7 | 12（回转体 9） | 增 5 |
| 代理模型 | 4 节点 | 8 节点 | 节点数变化 |
| 全局并集 | $ \max $/KS | 形式不变 | 无 |
| 几何灵敏度 | $ \partial\phi/\partial\theta $ 等 | 多 $ z_0 $、第二角 $ \psi $、$ z $ 向厚度偏导 | 链式法则项增多 |
| SIMP 第二步 | — | **完全不变** | 0 |


> **结论**：三维化的"方法层"改动集中在 §7.1–7.4 的几何描述与 §7.7(b)(c) 的新增偏导（多一个旋转角与一组 $ z $ 方向厚度变量）；链式法则的整体框架与二维一致，可直接照搬。真正的成本仍在三维有限元求解（参见前述讨论）。3D MMC 已有成熟文献（如郭旭团队的三维 MMC/MMV 工作）可供参考实现。
>
