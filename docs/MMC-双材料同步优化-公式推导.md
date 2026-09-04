# MMC–双材料同步优化：全链公式与灵敏度推导

> 对应实现：`Apps/SolverGPU/Cuda_TopOpt_MMC_MultiMaterial_Simultaneous.cu`
> 本文完整写出同步优化从设计变量到目标/约束灵敏度的全部公式，自成体系，不引用其他流程的定义。

---

## 1. 问题设定

### 1.1 设计变量

同步优化同时持有两块设计变量：

**几何块（MMC 组件）**：$M$ 个组件，每个组件 9 个物理变量

$$
d_i=\bigl(x_0^i,\;y_0^i,\;z_0^i,\;\ell_1^i,\;\ell_2^i,\;\ell_3^i,\;\alpha^i,\;\beta^i,\;\gamma^i\bigr)\in\mathbb R^{9},\qquad i=1,\dots,M
$$

含义依次为：组件中心平移 $(x_0,y_0,z_0)$、半轴长 $(\ell_1,\ell_2,\ell_3)$、姿态角 $(\alpha,\beta,\gamma)$。变量以物理量纲直接进入优化器，不做归一化换元。

**材料块（双材料插值变量）**：固定网格上每个单元一个

$$
\eta_e\in[\eta_{\min},\,1],\qquad e=1,\dots,N
$$

$\eta_e=1$ 表示纯材料 1（软，$E_1$），$\eta_e=0$ 表示纯材料 2（硬，$E_2$）。材料变量始终固定在完整网格上，MMC 组件移动时不按阈值增删变量。

### 1.2 优化问题

$$
\begin{aligned}
\min_{d,\;\eta}\quad & C(d,\eta)=F^{\mathsf T}U \\
\text{s.t.}\quad & K(d,\eta)\,U=F \\
& g_{\mathrm{vol}}(d)\le 0 \qquad(\text{MMC 几何体积约束，由 MMA 处理})\\
& \text{材料比例由 OC 二分满足（见 \S8，不进 MMA)}
\end{aligned}
$$

两条变量链通过**单元等效刚度系数** $a_e$ 耦合：

$$
a_e=m_e(d)\cdot E_e(\eta_e),\qquad K_e=1000\cdot a_e\,K_e^{0}
$$

其中 $m_e$ 是单元的 MMC 占有率（\S2 推导），$E_e$ 是单元插值模量（\S3 推导），$K_e^0$ 是单位杨氏模量下的单元刚度矩阵，系数 $1000$ 为刚度装配内核的固定缩放（见 \S4.2）。

### 1.3 链式总览

```
几何块 d_ir ─► 组件TDF φ_c(v) ─► KS并集 Φ(v) ─► Heaviside H(Φ) ─► 占有率 m_e ─┐
                                                                             ├─► a_e = m_e·E_e ─► K_e=1000·a_e·K_e⁰ ─► K U = F ─► U, q_e ─► C
材料块 η_e ──────────────────────────────► 插值模量 E_e(η_e) ────────────────┘
```

每轮迭代只求解**一次**共享平衡方程，同一个 $U$ 与单元能量 $q_e$ 同时用于几何块和材料块两组梯度。

### 1.4 参数表（默认值）

| 符号 | 代码名 | 默认值 | 含义 |
| --- | --- | --- | --- |
| $n$ | `kTdfPower` | 6 | TDF 超椭球幂次 |
| $p$ | `material_penalty` | 3 | 材料插值惩罚 |
| $\lambda$ | `ks_lambda` | 100 | KS 并集聚合参数 |
| $\alpha_H$ | `heaviside_alpha` | $10^{-3}$ | Heaviside 空域下限（占有率地板） |
| $\varepsilon_H$ | `heaviside_epsilon` | 0.20 | Heaviside 过渡带半宽 |
| $E_1,E_2$ | `E1,E2` | 1, 100 | 软/硬材料杨氏模量 |
| $f$ | `material_ratio_f` | 1.5 | 目标体积比 $V_2=f\,V_1$ |
| $V^*$ | `mmc_target_rate` | 0.3 | MMC 几何体积上限（体积分数） |
| $\eta_{\min}$ | `eta_min` | $10^{-3}$ | 材料变量下界 |
| $\eta^{0}$ | `initial_eta` | $1/(1+f)=0.4$ | 材料变量统一初值 |
| $\nu$ | `Poissons` | 0.3 | 泊松比（两种材料共用） |
| $\kappa$ | `kTiny` | $10^{-12}$ | 坐标正则项，防 $\ pow(0,\cdot)$ |

---

## 2. 几何链：组件变量 → 单元占有率 $m_e$

### 2.1 组件局部坐标与旋转矩阵

对任一网格顶点 $X=(X,Y,Z)^{\mathsf T}$，先平移到组件中心：

$$
\tilde X = X-X_0^{(i)}=(X-x_0,\;Y-y_0,\;Z-z_0)^{\mathsf T}
$$

再旋转到组件局部系：

$$
\begin{pmatrix}x\\y\\z\end{pmatrix}=R(\alpha,\beta,\gamma)\,\tilde X,
\qquad
R=\begin{pmatrix}
c_\beta c_\gamma & c_\beta s_\gamma & -s_\beta\\
s_\alpha s_\beta c_\gamma-c_\alpha s_\gamma & s_\alpha s_\beta s_\gamma+c_\alpha c_\gamma & s_\alpha c_\beta\\
c_\alpha s_\beta c_\gamma+s_\alpha s_\gamma & c_\alpha s_\beta s_\gamma-s_\alpha c_\gamma & c_\alpha c_\beta
\end{pmatrix}
$$

其中 $s_\alpha=\sin\alpha$、$c_\alpha=\cos\alpha$ 等。可以验证 $R=R_x(-\alpha)R_y(-\beta)R_z(-\gamma)$。实现中再取

$$
x_1=x+\kappa,\quad y_1=y+\kappa,\quad z_1=z+\kappa,\qquad
\hat\ell_j=\ell_j+\varepsilon_{\mathrm{mach}}
$$

作为数值正则（$\kappa=10^{-12}$），导数同样对 $x_1,y_1,z_1,\hat\ell_j$ 计算。

$R$ 对三个姿态角的偏导矩阵（用于 \S2.3 的角度导数，逐项对上表求导即可验证）：

$$
R_\alpha=\frac{\partial R}{\partial \alpha}=
\begin{pmatrix}
0&0&0\\
R_{31}&R_{32}&R_{33}\\
-R_{21}&-R_{22}&-R_{23}
\end{pmatrix},\qquad
R_\beta=\frac{\partial R}{\partial \beta}=
\begin{pmatrix}
-s_\beta c_\gamma&-s_\beta s_\gamma&-c_\beta\\
s_\alpha c_\beta c_\gamma&s_\alpha c_\beta s_\gamma&-s_\alpha s_\beta\\
c_\alpha c_\beta c_\gamma&c_\alpha c_\beta s_\gamma&-c_\alpha s_\beta
\end{pmatrix}
$$

$$
R_\gamma=\frac{\partial R}{\partial \gamma}=
\begin{pmatrix}
-c_\beta s_\gamma&c_\beta c_\gamma&0\\
-s_\alpha s_\beta s_\gamma-c_\alpha c_\gamma&s_\alpha s_\beta c_\gamma-c_\alpha s_\gamma&0\\
-c_\alpha s_\beta s_\gamma+s_\alpha c_\gamma&c_\alpha s_\beta c_\gamma+s_\alpha s_\gamma&0
\end{pmatrix}
$$

（$R_{ij}$ 表示上面 $R$ 的第 $i$ 行第 $j$ 列元素。）

### 2.2 单组件 TDF（超椭球带符号距离函数）

组件 $i$ 在顶点处的拓扑描述函数：

$$
T=\Bigl(\frac{x_1}{\ell_1}\Bigr)^{n}+\Bigl(\frac{y_1}{\ell_2}\Bigr)^{n}+\Bigl(\frac{z_1}{\ell_3}\Bigr)^{n},\qquad
\boxed{\;\varphi_i=1-T^{1/n}\;},\qquad n=6
$$

- $\varphi_i>0 \Leftrightarrow T<1$：顶点在超椭球**内部**；
- $\varphi_i=0$：在边界上；$\varphi_i<0$：在外部。

### 2.3 TDF 对 9 个组件变量的解析导数

由 $\varphi=1-T^{1/n}$ 得公共因子 $T^{1/n-1}$。先算坐标方向的三个中间量：

$$
A_x\triangleq T^{\frac1n-1}\Bigl(\frac{x_1}{\ell_1}\Bigr)^{n-1}\frac{1}{\ell_1},\qquad
A_y\triangleq T^{\frac1n-1}\Bigl(\frac{y_1}{\ell_2}\Bigr)^{n-1}\frac{1}{\ell_2},\qquad
A_z\triangleq T^{\frac1n-1}\Bigl(\frac{z_1}{\ell_3}\Bigr)^{n-1}\frac{1}{\ell_3}
$$

**对局部坐标的导数**（后续链式用）：

$$
\frac{\partial\varphi}{\partial x_1}
=-\frac1n T^{\frac1n-1}\cdot\frac{\partial T}{\partial x_1}
=-\frac1n T^{\frac1n-1}\cdot n\Bigl(\frac{x_1}{\ell_1}\Bigr)^{n-1}\frac{1}{\ell_1}
=-A_x
$$

同理 $\partial\varphi/\partial y_1=-A_y$、$\partial\varphi/\partial z_1=-A_z$。

**对半轴长的导数**（$\partial T/\partial \ell_1=n(x_1/\ell_1)^{n-1}\cdot(-x_1/\ell_1^2)$）：

$$
\frac{\partial\varphi}{\partial\ell_1}
=-\frac1n T^{\frac1n-1}\cdot n\Bigl(\frac{x_1}{\ell_1}\Bigr)^{n-1}\Bigl(-\frac{x_1}{\ell_1^{2}}\Bigr)
=+A_x\,\frac{x_1}{\ell_1},\qquad
\frac{\partial\varphi}{\partial\ell_2}=A_y\,\frac{y_1}{\ell_2},\qquad
\frac{\partial\varphi}{\partial\ell_3}=A_z\,\frac{z_1}{\ell_3}
$$

**对平移量的导数**：角度不变时 $x,y,z$ 是 $\tilde X$ 的线性函数，$\partial(x,y,z)^{\mathsf T}/\partial x_0=-R$ 的第一列，于是

$$
\frac{\partial\varphi}{\partial x_0}=A_xR_{11}+A_yR_{21}+A_zR_{31},\qquad
\frac{\partial\varphi}{\partial y_0}=A_xR_{12}+A_yR_{22}+A_zR_{32},\qquad
\frac{\partial\varphi}{\partial z_0}=A_xR_{13}+A_yR_{23}+A_zR_{33}
$$

**对姿态角的导数**：中心不变时 $\partial(x,y,z)^{\mathsf T}/\partial\alpha=R_\alpha\tilde X\triangleq(x_\alpha,y_\alpha,z_\alpha)^{\mathsf T}$（$\beta,\gamma$ 同理换 $R_\beta,R_\gamma$），于是

$$
\frac{\partial\varphi}{\partial\alpha}=-(A_xx_\alpha+A_yy_\alpha+A_zz_\alpha),\qquad
\frac{\partial\varphi}{\partial\beta}=-(A_xx_\beta+A_yy_\beta+A_zz_\beta),\qquad
\frac{\partial\varphi}{\partial\gamma}=-(A_xx_\gamma+A_yy_\gamma+A_zz_\gamma)
$$

以上 9 个导数在实现中与 $\varphi$ 同步解析计算（不做数值差分）。

### 2.4 多组件并集：KS 聚合

结构是各组件的**并集**，需要平滑近似 $\max_c\varphi_c$。对每个顶点 $v$：

$$
\Phi(v)=\frac1\lambda\ln\sum_{c=1}^{M}e^{\lambda\varphi_c(v)},\qquad
\Phi(v)\leftarrow\max\bigl(\Phi(v),\,-10^{3}\bigr)
$$

KS 满足 $\max_c\varphi_c\le\Phi\le\max_c\varphi_c+\ln M/\lambda$（$\lambda=100$ 时偏差上界 $\ln M/100$）。

KS 对各组件 TDF 的链式权：

$$
\boxed{\;w_{v,c}\triangleq\frac{\partial\Phi(v)}{\partial\varphi_c(v)}
=\frac{e^{\lambda\varphi_c(v)}}{\sum_{j=1}^{M}e^{\lambda\varphi_j(v)}+10^{-9}}\;}
$$

分母的 $+10^{-9}$ 是数值稳定项：它使 $w_{v,c}$ 之和略小于 1，解析导数因此是严格 KS 导数的轻微近似；该偏差由有限差分日志量化（\S5.4）。

### 2.5 正则化 Heaviside

把全局 TDF 映射到 $[\alpha_H,\,1]$，过渡带 $[-\varepsilon_H,\varepsilon_H]$ 内用三次多项式保证 $C^1$：

$$
H(\Phi)=
\begin{cases}
1, & \Phi>\varepsilon_H\\[2pt]
\dfrac{3(1-\alpha_H)}{4}\Bigl(\dfrac{\Phi}{\varepsilon_H}-\dfrac{\Phi^{3}}{3\varepsilon_H^{3}}\Bigr)+\dfrac{1+\alpha_H}{2}, & |\Phi|\le\varepsilon_H\\[8pt]
\alpha_H, & \Phi<-\varepsilon_H
\end{cases}
$$

$$
H'(\Phi)=
\begin{cases}
0, & |\Phi|>\varepsilon_H\\[2pt]
\dfrac{3(1-\alpha_H)}{4\varepsilon_H}\Bigl(1-\dfrac{\Phi^{2}}{\varepsilon_H^{2}}\Bigr), & |\Phi|\le\varepsilon_H
\end{cases}
$$

可验证 $H(\varepsilon_H)=1$、$H(-\varepsilon_H)=\alpha_H$、$H'(\pm\varepsilon_H)=0$。$\alpha_H=10^{-3}$ 是空域占有率地板（保证刚度矩阵非奇异）；**$H'$ 只在过渡带内非零**，因此几何灵敏度只在组件边界带上有支撑。

### 2.6 单元占有率

六面体单元 $e$ 的占有率取其 8 个顶点 Heaviside 值的算术平均：

$$
\boxed{\;m_e=\frac18\sum_{v\in e}H\bigl(\Phi(v)\bigr)\;},\qquad m_e\in[\alpha_H,\,1]
$$

### 2.7 占有率对几何变量的全链导数

综合 2.3–2.6 四级链式：

$$
\boxed{\;
\frac{\partial m_e}{\partial d_{ir}}
=\frac18\sum_{v\in e}H'\bigl(\Phi(v)\bigr)\cdot w_{v,i}\cdot\frac{\partial\varphi_i}{\partial d_{ir}}(v)
\;}
$$

只有拥有变量 $d_{ir}$ 的组件 $i$ 出现在求和中（其他组件的 TDF 不含 $d_{ir}$，其 KS 权乘以零导数自动消去）。

---

## 3. 材料链：$\eta_e$ → $E_e$ 与耦合刚度

### 3.1 双材料插值

$$
\boxed{\;E_e=\eta_e^{\,p}E_1+\bigl(1-\eta_e^{\,p}\bigr)E_2\;},\qquad p=3
$$

导数：

$$
\frac{dE_e}{d\eta_e}=p\,\eta_e^{\,p-1}\bigl(E_1-E_2\bigr)=3\eta_e^{2}(E_1-E_2)
$$

由于 $E_2>E_1$，有 $dE_e/d\eta_e<0$：**$\eta_e$ 越大材料越软**。这一符号约定贯穿后续所有灵敏度与 OC 分支方向。

### 3.2 耦合点与单元刚度

$$
a_e=m_e\,E_e,\qquad
K_e=1000\cdot a_e\,K_e^{0}=1000\cdot m_e\,E_e\,K_e^{0}
$$

$K_e^{0}\in\mathbb R^{192\times192}$ 为单位杨氏模量、共同泊松比 $\nu$ 下的单元刚度矩阵（每个六面体是 64 个控制点 × 3 自由度的三次 Bezier 体单元），**只在初始化时计算一次**；每轮仅用当前 $a_e$ 缩放装配。MMC 占有率与材料模量的耦合**只**发生在 $a_e$ 这一个点上。

---

## 4. 结构分析与柔顺度

### 4.1 共享平衡方程

每轮迭代对当前联合状态求解一次：

$$
K(d,\eta)\,U=F,\qquad K=\underset{e}{\text{assemble}}\;1000\,a_e\,K_e^{0}
$$

$F$ 由 Neumann 载荷点组装，Dirichlet 约束在系统层面消元；$F$ 不依赖设计变量。

### 4.2 柔顺度及其双口径

单元能量采用单位模量矩阵：

$$
q_e\triangleq U_e^{\mathsf T}K_e^{0}\,U_e
$$

则

$$
C=F^{\mathsf T}U=U^{\mathsf T}KU=1000\sum_{e=1}^{N}a_e\,q_e
=1000\sum_{e=1}^{N}m_e\,E_e\,q_e
$$

实现以**归一化口径**为目标与梯度基准：

$$
C_n\triangleq\sum_{e=1}^{N}a_e\,q_e=\frac{C}{1000}
$$

两个口径每轮做相对差校验（$F^{\mathsf T}U/1000$ 与 $\sum a_eq_e$，容差 $5\times10^{-3}$），超差即抛错——它是"装配的 $a_e$"与"求和用的 $a_e$"一致性的直接检查。

### 4.3 伴随法求 $dC/dx$（完整推导）

对任一设计变量 $x\in\{d_{ir},\eta_e\}$，平衡方程两边取微分（$F$ 与 $x$ 无关）：

$$
dK\,U+K\,dU=0\quad\Longrightarrow\quad dU=-K^{-1}\,dK\,U
$$

$$
dC=F^{\mathsf T}dU=U^{\mathsf T}dU=-U^{\mathsf T}\,dK\,U
$$

（柔顺度是自伴随的，不需要另解伴随方程。）而 $K_e$ 只通过标量 $a_e$ 依赖设计：

$$
dK=1000\sum_{e}K_e^{0}\,da_e
\quad\Longrightarrow\quad
dC=-1000\sum_{e}U_e^{\mathsf T}K_e^{0}U_e\;da_e=-1000\sum_{e}q_e\,da_e
$$

对耦合点 $a_e=m_eE_e$ 全微分：

$$
da_e=E_e\;dm_e+m_e\;\frac{dE_e}{d\eta_e}\;d\eta_e
$$

代入即得两个偏导数：

$$
\boxed{\;
\frac{\partial C}{\partial m_e}=-1000\,E_e\,q_e,
\qquad
\frac{\partial C}{\partial \eta_e}=-1000\,m_e\,q_e\,\frac{dE_e}{d\eta_e}=-1000\,m_e\,q_e\cdot3\eta_e^{2}(E_1-E_2)
\;}
$$

（归一化口径 $C_n$ 去掉因子 $1000$，下同。）

---

## 5. 两组设计变量的完整灵敏度

### 5.1 几何块：$\partial C_n/\partial d_{ir}$

把 \S4.3 的 $\partial C/\partial m_e$ 与 \S2.7 的 $\partial m_e/\partial d_{ir}$ 链起来：

$$
\boxed{\;
\frac{\partial C_n}{\partial d_{ir}}
=-\sum_{e=1}^{N}E_e\,q_e\cdot\frac18\sum_{v\in e}H'(\Phi_v)\,w_{v,i}\,\frac{\partial\varphi_i}{\partial d_{ir}}(v)
\;}
$$

实现采用**顶点累加形式**：先把单元权重摊到 8 个顶点

$$
\hat w_v=\frac18\sum_{e\ni v}\bigl(-E_e\,q_e\bigr),
$$

再对每个几何变量做一次顶点求和：

$$
\frac{\partial C_n}{\partial d_{ir}}
=\sum_{v}\hat w_v\;H'(\Phi_v)\;w_{v,i}\;\frac{\partial\varphi_i}{\partial d_{ir}}(v)
$$

> **关键点**：$E_e$ 必须在单元贡献累加到顶点**之前**逐单元乘进权重。同步状态下各单元 $E_e$ 不同，若先按 $-q_e$ 累加再乘单一模量，等于隐式假设了均匀材料场，结果错误。

### 5.2 材料块：$\partial C_n/\partial\eta_e$

$$
\boxed{\;
\frac{\partial C_n}{\partial\eta_e}=-m_e\,q_e\cdot3\eta_e^{2}\bigl(E_1-E_2\bigr)\;>0\quad(E_2>E_1)
\;}
$$

符号恒正：增大 $\eta_e$（趋向软材料 1）柔顺度上升；顺序版复刻的 OC 候选 $\eta_e\sqrt{g_e/\lambda}$ 直接消费该正梯度（\S8）。

同一个 $U$、同一个 $q_e$ 同时喂给 5.1 和 5.2——这就是"同步"在计算结构上的含义。

### 5.3 均匀材料退化一致性校验

当 $E_e\equiv\bar E$（全场相同）时，5.1 必须退化为

$$
\frac{\partial C_n}{\partial d_{ir}}\Big|_{E\equiv\bar E}
=\bar E\sum_{v}\hat w_v^{0}H'(\Phi_v)w_{v,i}\frac{\partial\varphi_i}{\partial d_{ir}}(v),
\qquad
\hat w_v^{0}=\frac18\sum_{e\ni v}(-q_e)
$$

实现检测到全场模量近似一致时逐轮计算该参考值并校验相对误差（容差 $10^{-9}$），验证 "$E_e$ 逐单元前乘" 的链式实现没有错位。

### 5.4 有限差分数值验证

对抽样变量用中心差分 $(C_n(x^+)-C_n(x^-))/2h$ 重解全状态并与解析值比对；差分通过/警告阈值分别为 $10^{-3}$、$10^{-2}$ 相对容差。输出 `analytical_sensitivities.tsv`、`all_finite_difference_sensitivities.tsv` 等日志用于离线核对（含 KS 稳定项 $10^{-9}$ 造成的理论偏差量化）。

---

## 6. 约束及其灵敏度

以下均用归一化单元体积 $v_e=1/N$（均匀网格）。

### 6.1 MMC 几何体积约束（进 MMA）

$$
g_{\mathrm{vol}}=\frac1N\sum_{e}m_e-V^*\le 0,\qquad V^*=0.3
$$

对几何变量的导数：以单元权重 $1/N$ 走 \S2.7 同一条反传链，

$$
\frac{\partial g_{\mathrm{vol}}}{\partial d_{ir}}
=\sum_{e}\frac1N\cdot\frac{\partial m_e}{\partial d_{ir}}
=\sum_{v}\Bigl(\frac18\sum_{e\ni v}\frac1N\Bigr)H'(\Phi_v)\,w_{v,i}\,\frac{\partial\varphi_i}{\partial d_{ir}}(v)
$$

它不依赖 $\eta$（$g_{\mathrm{vol}}$ 只含 $m_e$）。

### 6.2 材料比例约束（当前实现中仅作诊断）

占有率加权的双材料资源量：

$$
V_1=\frac1N\sum_{e}m_e\,\eta_e,\qquad
V_2=\frac1N\sum_{e}m_e\,(1-\eta_e),\qquad
R=V_2-f\,V_1
$$

**对材料变量的导数**（$m_e$ 视为常数）：

$$
\frac{\partial R}{\partial\eta_e}
=\frac{m_e}{N}\bigl(-1-f\bigr)
=-\frac{(1+f)\,m_e}{N}
$$

**对几何变量的导数**（$\eta_e$ 视为常数）：以单元权重 $\bigl[(1-\eta_e)-f\eta_e\bigr]/N$ 反传，

$$
\frac{\partial R}{\partial d_{ir}}
=\sum_{e}\frac{(1-\eta_e)-f\,\eta_e}{N}\cdot\frac{\partial m_e}{\partial d_{ir}}
$$

即使 $\eta$ 固定，组件移动也会改变 $V_1,V_2$——两组导数都在每轮诊断中输出。

初始 $\eta^{0}=1/(1+f)$ 的意义：若全场 $m_e$ 相同，则 $V_2/V_1=(1-\eta^0)/\eta^0=f$，即均匀初值精确满足比例约束（第 0 轮诊断 $R=0$）。

### 6.3 当前实现中约束的归属

- MMC-MMA 的约束集合**只有** $g_{\mathrm{vol}}$；$R$ 及其两块梯度只写入诊断日志，不进任何优化器。
- 材料比例由 OC 二分的目标量满足；二分目标为活动单元（$m_e>0.3$）上未加权的 $\sum_e\eta_e=N_{\mathrm{act}}/(1+f)$，与顺序版口径一致（\S8）。由于活动单元内 $m_e$ 并非严格为 1，该目标只**近似**对应物理比例 $V_2=fV_1$。

---

## 7. 每轮执行顺序

1. 由当前 $d$ 计算 $\varphi_c\to\Phi\to H\to m_e$ 及全部几何链式因子（\S2）；
2. 由当前 $\eta$ 计算 $E_e$、$a_e=m_eE_e$（\S3）；
3. 装配 $K=1000\sum a_eK_e^0$，求解**一次** $KU=F$（\S4.1）；
4. 由同一个 $U$ 计算 $q_e$、$C_n$、两组目标梯度（\S5）、约束值与约束梯度（\S6）；
5. 用同一当前状态分别更新：MMA 更新一次 $d$，过滤+OC 更新一次 $\eta$（\S8）；
6. 接受两块更新，记录诊断，进入下一轮。固定运行 `max_iterations` 轮，变化量与残差只作诊断、不触发早停。

循环结束后用最终设计变量**再重算一次**完整状态用于导出（不能把更新前状态当最终结果）。

---

## 8. 优化器更新（简述）

- **几何块（MMA）**：持久化 MMA 求解器，九个物理变量连同完整物理上下界直接进入；约束只带 $g_{\mathrm{vol}}$；目标/约束梯度为 \S5.1、\S6.1 的解析值，不做归一化或缩放。
- **材料块（过滤 + OC，逐项复刻顺序版）**：变量集合为当前 $m_e>0.3$ 的活动单元（每轮随几何重建），活动单元外的 $\eta_e$ 保持不变。滤波输入中活动单元外的 $\rho$ 与 $dc$ 置零，解析梯度在活动单元上过 Sigmund 灵敏度过滤

  $$
  \widehat{g}_i=\frac{\sum_j w_{ij}\,\eta_j\,g_j}{\max(\eta_i,\eta_{\min})\sum_j w_{ij}},\qquad g_j=\partial C_n/\partial\eta_j
  $$

  （$w_{ij}$ 为固定邻域权重；过滤结果只作为 OC 输入，原始解析梯度单独保存诊断。）随后 OC 二分更新

  $$
  \eta_e^{\text{new}}=\mathrm{clamp}_{[\eta_{\min},1]}\Bigl(\eta_e\pm\Delta_{\max},\;\;\eta_e\sqrt{\widehat{g}_e/\lambda}\Bigr),\qquad \Delta_{\max}=0.2
  $$

  与顺序版 `x.*sqrt(-dc./lmid)` 代数相同（顺序版以 $dc=-g$ 的负值约定传入同一量）。候选随 $\lambda$ 增大而减小，二分在 $\sum$ 超标时抬高 $\lambda$；二分资源目标为活动单元上的未加权 $\sum_e\eta_e=N_{\mathrm{act}}/(1+f)$。

  > **设计备注**：曾实现并验证过一个替代口径——资源目标按占有率加权 $\sum_e m_e\eta_e=\frac{1}{1+f}\sum_e m_e$、候选取 $\eta_e\sqrt{\lambda/g_e}$（此时最优性条件 $g_e=\mu$ 要求高灵敏度单元减小 $\eta_e$）。该变体使物理比例 $V_2/V_1$ 精确收敛到 $f$（残差约 $3\times10^{-4}$）且软材料落在低应变能单元；离线重放对照显示顺序版方向将软材料推向高应变能单元（软/硬单元平均能量比 $1.98$ 对 $0.40$）。按当前决定，为与顺序版逐项可对照，实现回退为完全复刻顺序版 OC；两种口径的对照数据保留在 openspec 设计文档中。

---

## 9. 数值一致性校验汇总

| 校验 | 公式 | 容差 | 失败行为 |
| --- | --- | --- | --- |
| 柔顺度双口径 | $\|C_n-F^{\mathsf T}U/1000\|/C_n$ | $5\times10^{-3}$ | 抛错终止 |
| 均匀材料退化 | \S5.3 相对误差 | $10^{-9}$ | 抛错终止 |
| 有限差分 | 解析 vs 中心差分 | $10^{-3}$/$10^{-2}$ | 记录日志 |
| 全链有限值 | 所有场/梯度 `isfinite` | — | 抛错终止 |

---

## 10. 实现对照表

| 公式 | 实现位置（函数） |
| --- | --- |
| $H(\Phi),H'(\Phi)$ | `Heaviside` / `HeavisideDerivative` |
| $\varphi_i$ 及 9 导数（\S2.2–2.3） | `calc_Phi_Node` |
| $\Phi$、$w_{v,c}$（\S2.4） | `ComputeGlobalTDF` |
| $m_e$（\S2.6） | `EvaluateGeometry` |
| $E_e,dE_e/d\eta_e$（\S3.1） | `InterpolateYoungsModulus` / `InterpolateYoungsModulusDerivative` |
| 单轮全部场、求解、两组梯度、约束（\S4–6） | `EvaluateJointState` |
| 通用几何反传（\S2.7、5.1、6.1、6.2 的 $d$ 部分） | `ComputeMmcDerivativeFromElementWeights`（单元权重分别为 $-E_eq_e$、$1/N$、$[(1-\eta_e)-f\eta_e]/N$） |
| $K_e^0$ 一次性构造 | 主函数 `Get_Local_Matrix(..., 1.0, ν)` |
| 灵敏度过滤、OC（\S8） | `FilterSensitivities` / `OCUpdate` |
| 有限差分验证（\S5.4） | `RunFiniteDifferenceValidation` / `RunAllFiniteDifferenceValidation` |
