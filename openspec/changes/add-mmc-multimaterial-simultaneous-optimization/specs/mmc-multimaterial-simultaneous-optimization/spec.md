## Purpose

提供一条独立的拓扑优化模式，基于同一次结构响应同步更新显式 3D MMC 几何与双材料单元分布，同时保持所有现有优化模式不变。

## ADDED Requirements

### Requirement: 独立的同步优化模式
系统 SHALL 通过新增公开函数 `Cuda_TopOpt_MMC_MultiMaterial_Simultaneous(...)` 及独立求解模块实现 MMC-多材料同步优化。实现 MUST NOT 修改独立 MMC、独立 SIMP 或顺序 MMC+SIMP 的数值算法、控制流程、参数语义、调用关系、结果或可用性。为比较同步版与独立 MMC 的公共计算链，独立 MMC MAY 增加只读诊断日志，但日志值 MUST NOT 反馈到求解、更新或收敛判定。

#### Scenario: 选择新模式
- **WHEN** 用户使用有效的几何、边界条件、MMC 参数和双材料参数选择 MMC-多材料同步优化
- **THEN** 系统运行新的同步优化流程，而不是先完整执行独立 MMC 或顺序 MMC+SIMP 流程

#### Scenario: 新模式加入后运行现有模式
- **WHEN** 用户选择任意已有优化模式
- **THEN** 系统仍调用该模式原有函数并沿用原求解路径和语义，不进入新增同步函数

### Requirement: MMC 优化界面中的独立入口
系统 SHALL 在现有 MMC 优化界面中提供“MMC-多材料同步优化”独立选项。新增选项 SHALL 收集同步模式所需的 MMC 与双材料参数，并且 MUST NOT 替换或改变原 MMC 优化选项。

#### Scenario: 打开 MMC 优化界面
- **WHEN** 用户进入当前用于启动 3D MMC 优化的界面
- **THEN** 标准 MMC 选项和 MMC-多材料同步优化选项在同一界面中可见

#### Scenario: 选择标准 MMC
- **WHEN** 用户在该界面选择原有标准 MMC 优化
- **THEN** 系统调用原 MMC 入口，参数、默认值和求解流程保持不变

#### Scenario: 选择同步优化
- **WHEN** 用户在该界面选择 MMC-多材料同步优化并提供有效的双材料参数
- **THEN** 系统只调用新增公开函数 `Cuda_TopOpt_MMC_MultiMaterial_Simultaneous(...)`

### Requirement: 每轮共享结构分析
每次同步迭代中，系统 SHALL 计算当前 MMC 单元占有率 `m_e`、当前插值单元模量 `E_e` 和等效单元刚度系数 `a_e = m_e E_e`。系统 SHALL 先由当前 `a_e` 场装配并求解一次平衡方程，再计算两组设计变量。

#### Scenario: 计算一次同步迭代
- **WHEN** 当前 MMC 组件变量或材料变量与上一轮不同
- **THEN** 系统在更新任意变量块之前重新计算 `m_e`、`E_e`、`a_e`、位移和单元能量

#### Scenario: MMC 外部单元
- **WHEN** 一个单元只有 MMC 替代材料占有率
- **THEN** 该单元的等效刚度贡献由该占有率及其当前插值材料模量共同缩放，而不是被当作完全占有的材料单元

### Requirement: 固定网格的双材料设计变量
系统 SHALL 在整个同步优化过程中，为固定设计网格的每个单元保留一个材料设计变量。系统 MUST NOT 根据移动的 MMC 占有率设置硬阈值来增加或删除材料变量。

#### Scenario: MMC 边界跨越单元
- **WHEN** MMC 更新使某些单元进入或离开占有区域
- **THEN** 材料设计向量的维数及单元到变量的映射保持不变

### Requirement: MMC 物理变量与完整边界
系统 SHALL 在联合 MMA 中直接使用与独立 MMC 相同的九个物理控制变量，并 SHALL 将其完整物理上下界传给优化器。系统 MUST NOT 在同步流程外层将 MMC 变量归一化为单位区间，也 MUST NOT 用当前值附近的局部移动区间替代完整物理上下界。目标和约束对 MMC 块的输入 SHALL 是相对于物理控制变量的解析导数。

#### Scenario: 联合更新 MMC 变量
- **WHEN** 同步流程构造联合 MMA 设计向量与上下界
- **THEN** MMC 块依次包含每个组件的 `x0,y0,z0,l1,l2,l3,alpha,beta,gamma` 物理值及其完整物理上下界
- **AND** 日志中的 MMC 最大变化量表示这些物理变量的实际最大变化，而不是归一化变化

### Requirement: 耦合的 MMC 目标灵敏度
系统 SHALL 复用从组件控制变量 `d_ir` 到单元占有率 `m_e` 的既有 3D MMC 导数。对于同步目标，每个单元的 MMC 导数贡献 SHALL 在单元累加内部使用该单元当前的插值材料模量加权：

`dC/dd_ir = -sum_e(E_e q_e dm_e/dd_ir)`。

系统 MUST NOT 使用单个全局材料模量去乘已经完成单元汇总的 MMC 梯度，以近似上述表达式。

#### Scenario: 单元材料模量不同
- **WHEN** 两个占有单元的材料变量不同，因此具有不同的 `E_e`
- **THEN** 它们对同一 MMC 控制变量梯度的贡献分别使用各自的单元模量加权

#### Scenario: 材料场暂时均匀
- **WHEN** 某次迭代中所有单元具有相同材料变量
- **THEN** 耦合 MMC 梯度与均匀材料模量缩放下的既有 MMC 占有率梯度等价

### Requirement: 耦合的材料目标灵敏度
初始能力不增加材料密度滤波或投影，系统 SHALL 使用现有插值公式：

`E_e = eta_e^p E1 + (1 - eta_e^p) E2`

并 SHALL 根据当前 MMC 占有率和共享单元能量计算材料灵敏度：

`dC/deta_e = -m_e q_e dE_e/deta_e`。

#### Scenario: 单元位于 MMC 结构内部
- **WHEN** 一个单元具有完全 MMC 占有率
- **THEN** 系统根据该单元的插值导数和当前单元能量计算材料梯度

#### Scenario: 单元位于 MMC 结构外部
- **WHEN** 一个单元只有 MMC 替代材料占有率
- **THEN** 其材料梯度幅值按该占有率相应衰减

### Requirement: 几何与材料资源约束
系统 SHALL 根据 MMC 占有率场施加 MMC 几何体积约束，并 SHALL 使用按占有率加权的材料度量维持请求的双材料资源比例。几何体积对 MMC 控制变量的灵敏度 SHALL 复用既有 MMC 体积灵敏度计算。

#### Scenario: 材料变量固定时 MMC 几何变化
- **WHEN** MMC 更新改变占有体积
- **THEN** 几何体积约束和按占有率加权的材料资源度量均反映新的占有率场

#### Scenario: MMC 几何固定时材料变量变化
- **WHEN** 材料更新改变软材料与硬材料的分布
- **THEN** 材料资源度量发生变化，而 MMC 几何体积度量保持不变

### Requirement: 固定外层迭代与诊断
系统 SHALL 在每次迭代中报告目标值、MMC 几何体积、材料资源度量、MMC 变量最大变化量和材料变量最大变化量。除输入、有限元求解或结果导出失败外，系统 SHALL 固定执行请求的 `max_iterations` 次联合更新。变量变化量与活动约束残差 SHALL 只用于诊断，不得触发提前终止。

#### Scenario: 变量变化量很小
- **WHEN** MMC 与材料变量变化均低于诊断容差
- **THEN** 同步优化仍继续运行，直到完成 `max_iterations` 次联合更新

#### Scenario: 完成请求的迭代次数
- **WHEN** 同步优化已经接受 `max_iterations` 次联合更新
- **THEN** 系统成功终止并基于最终联合设计重新求解和导出结果字段

### Requirement: 结果与解析灵敏度输出
系统 SHALL 将 MMC 占有率、材料设计变量、插值材料模量和等效刚度系数作为相互独立的逐单元输出。系统 SHALL 支持输出选定 MMC 和材料变量的解析灵敏度，并且该诊断过程 SHALL 只读取当前已接受状态。

#### Scenario: 导出完成的结果
- **WHEN** 同步优化终止
- **THEN** 导出结果包含可分别识别的 `m_e`、`eta_e`、`E_e` 和 `a_e` 字段

#### Scenario: 查看组件内的二值材料分布
- **WHEN** 用户在 ParaView 中查看同步优化 VTK 的材料二值显示字段
- **THEN** 系统按与顺序流程一致的 MMC 活动阈值只将组件内单元分类为材料 1 或材料 2
- **AND** 组件外单元使用独立类别及浅灰色，不能被误判为任一种组件内材料
- **AND** 系统额外提供仅包含 MMC 活动单元的材料 VTK，并在该网格上按 `eta_e>=0.5` 输出与顺序 MMC-SIMP 同语义的 `material_binary` 和红蓝直接颜色

#### Scenario: 查看解析梯度
- **WHEN** 为选定变量启用灵敏度验证
- **THEN** 系统按变量块和变量索引报告当前已接受状态下的解析灵敏度及相关单元贡献摘要，并且不推进优化器历史

#### Scenario: 对照同步版与独立 MMC
- **WHEN** 用户分别运行同步优化和独立 MMC
- **THEN** 两条路径分别输出同列结构的 `mmc_debug_trace.tsv`，其中包含固定迭代状态、步长、梯度关系以及 MMC 设计变量、全局 TDF、Heaviside、单元占有率和目标/体积梯度的统计量与固定索引样本
- **AND** 独立 MMC 的日志插桩不改变任何设计变量、灵敏度、MMA 输入、更新顺序或退出条件
