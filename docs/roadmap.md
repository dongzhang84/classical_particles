# 经典带电质点的运动：研究路线图

> 版本：v2.1 · 2026-10-23 · 草稿
> 状态：第 0 步（文献核查）已完成三轮，进入第 1 步。所有"结果"仍为待检验的猜想。
> 文献核查的过程记录见 [`docs/literature-review.md`](literature-review.md)（下称"调研记录"）。旧版本见第 9 节版本记录。

---

## 0. 这件事要做什么

**出发点**：只用经典理论，即麦克斯韦电动力学和爱因斯坦广义相对论，把"一个质量为 m、电荷为 q 的质点会怎样运动"这个问题推到底，再推到两个质点。看理论自己会不会给出意想不到的结论。

**方法论**：不去修补方程，而是把方程当真。爱因斯坦当年就是这样做的：他没有修改麦克斯韦方程，而是认真对待它和力学之间的矛盾，矛盾自己推出了相对论。

### 0.1 为什么用经典理论，而不用量子电动力学

不是因为量子电动力学（QED）完全不能和引力放在一起。QED 可以放在一个给定的弯曲时空上（弯曲时空量子场论）；引力本身也可以作为低能有效场论量子化，在远低于普朗克能标时给出可靠的修正（Donoghue 1994）。

真正的障碍在于本项目关心的区域。ADM 的结果表明，引力"治愈"点电荷自能发散，发生在
$$\frac{\sqrt G\, q}{c^2} \approx \sqrt\alpha\,\ell_{\rm Pl} \approx 1.4\times10^{-34}\ {\rm cm}$$
的尺度上。在这里，引力自能与电磁自能同量级，时空弯曲必须**非微扰地、动态地**反作用于电荷本身。而在这个尺度上：

- 弯曲时空量子场论不够：它的时空背景是固定的，没有反作用。
- 引力有效场论失效：它只在远低于普朗克尺度时成立。
- 完整的量子引力理论还没有。

所以，能**同时、自洽地**包含"电荷的自场"和"这个自场产生的完整时空弯曲"的理论，目前只有经典的 Maxwell–Einstein 理论。本项目的问题因此是经典理论内部的问题：**Maxwell 电动力学加 Einstein 引力，作为一个封闭的经典理论，能否自洽地描述一个带电质点？** Dirac (1938)、ADM (1960)、Newman 一系（第 1.2 节）都是在这个框架下提问的。

**代价，必须说清楚**：对真实的电子，量子效应在康普顿波长 $\hbar/mc$ 处就已介入，它比经典电子半径大 $1/\alpha \approx 137$ 倍；非相对论 QED 中也没有 runaway（Moniz & Sharp 1974）。因此本项目的结论首先是关于经典理论自身是否自洽的结论，用到真实电子上时要谨慎表述。这与第 2 节、第 3 步的适用范围问题是同一件事。

**在 $10^{-34}$ cm 这个尺度上，本该用什么理论？** 对真实电子，三套理论各缺一块：

| 理论 | 有引力吗 | 在 $\sqrt\alpha\,\ell_{\rm Pl}$ 处适用于电子吗 |
|---|---|---|
| 经典 Maxwell–Einstein | 有，而且是非微扰的 | 不适用：电子在康普顿波长 $\sim4\times10^{-11}$ cm 处就已经是量子的 |
| QED | 没有（或只有固定背景） | 适用；QED 中自能只是对数发散，不需要引力来救 |
| QED + 量子引力 | 有 | 理论不存在 |

**经典处理何时自洽：一个与引力无关的独立判据。** 对电荷为 $Q$、质量为 $m$ 的物体，
$$\frac{\text{经典半径}}{\text{康普顿波长}} = \frac{Q^2/mc^2}{\hbar/mc} = \frac{Q^2}{\hbar c}.$$
这个比值与质量无关。

- **电荷为一个基本电荷 $e$**：比值恒为 $\alpha = 1/137$。只带一个电荷的东西，经典结构永远藏在量子尺度里面。ADM 极限下的那个物体也不例外：质量 $q/\sqrt G = \sqrt\alpha\,m_{\rm Pl}$，康普顿波长 $\ell_{\rm Pl}/\sqrt\alpha$，经典尺寸 $\sqrt\alpha\,\ell_{\rm Pl}$，二者之比仍是 137。
- **电荷为 $Ne$**：比值为 $N^2\alpha$。$N\gtrsim12$ 时，经典半径已大于康普顿波长，经典处理自洽。

所以从量子这一侧也能看出：**本项目的经典结论天然适用于带大量电荷的宏观物体（如带电天体），而不适用于电子。** 这与第 2 节从引力一侧得到的结论一致，两条线互相印证。

### 0.2 想回答的两个核心问题

1. **狄拉克的选择题**：狄拉克 1938 年推导电子运动方程时承认，能动量守恒只约束了一部分，粒子的"机械动量" $P_{\rm mech}$（也就是 $B^\mu$）该怎么写，是他挑了最简单的一种。其他写法是不同的物理，还是同一个物理的不同写法？
2. **引力能否消除"自加速"**：狄拉克方程（Lorentz–Dirac 方程，下称 LD）有一个荒谬的解：电子不受力也会自己加速到接近光速（runaway）。把广义相对论加进来，这个解还在不在？如果在某些情况下消失，消失的条件是什么？

### 0.3 目标结果

**最初的设想**是一句话："困扰物理学家近百年的'电子自加速'悖论，在加入爱因斯坦引力之后消失了。"

三轮文献核查和数量级估算之后，这句话**对电子基本不可能成立**（第 2 节）。现行的目标表述是：

> 在经典 Maxwell–Einstein 理论中，带电物体是否出现自加速，由一个可以写下来的判据决定。对电荷质量比 $\xi = q/(\sqrt G\, m) \le 1$ 的物体，引力保证自加速不出现：它们无法被压到经典半径以内而不形成视界。对 $\xi>1$ 的物体，在一个确定的尺寸范围内，引力是否消除自加速取决于 GR 中正确的判据是什么；这是本项目要判决的问题。对电子这样 $\xi\sim10^{21}$ 的对象，引力无能为力。

这仍是目标，不是结论。$\xi\le1$ 的部分几乎是平凡的（第 3.1 节）；不平凡的部分在 $\xi>1$ 一侧，依赖第 2 步的判决（猜想 A）。

---

## 1. 已知的事实

做之前先把地图画清楚。本节只列核实过的结论；出处和核实程度见调研记录。

### 1.1 一个带电质点，只有电动力学

| 现象 | 说明 | 来源 |
|---|---|---|
| 自能发散 | 点电荷的场能量无穷大 | 经典 |
| 4/3 问题 | 由场的能量和动量算出的电磁质量之比为 4/3 而非 1 | Abraham, Lorentz；Poincaré 应力 |
| runaway 与预加速 | LD 方程的解可以不受力自加速，或在力到来之前开始运动 | Dirac 1938 |
| 经典电子半径 | $r_e = e^2/mc^2$，比康普顿波长小 $1/\alpha = 137$ 倍：经典理论在病灶出现前已被量子力学接管 | 经典 |
| 有尺寸电荷的 runaway 判据 | 粒子尺寸**大于经典半径**（静电自能等于质量时的半径）时，理论自洽，无 runaway、无非因果；运动方程是微分–差分方程。平直时空中这等价于裸质量为正 | Levine, Moniz & Sharp 1977；Moniz & Sharp 1977 |
| 降阶 | LD 的物理解落在一个"临界流形"上，Landau–Lifshitz（LL）方程是它的一阶近似 | Spohn 2000；Gralla, Harte & Wald 2009 |
| 推导的歧义 | 在不同的面上积分能动量守恒，得到 LD 方程或一个变静质量、无 Schott 项的方程；分布论的"自然"正则化唯一选出 LD；场张量中不可能出现高阶导数项，因为经典电磁学没有长度常数 | Rowe 1975、1978 |
| 量子理论中没有 runaway | 非相对论 QED 中点电子自质量为零、无 runaway | Moniz & Sharp 1974 |

### 1.2 一个带电质点，加上广义相对论

| 现象 | 说明 | 来源 |
|---|---|---|
| **引力治好了自能发散** | 带电球壳在 Einstein–Maxwell 理论中，半径 → 0 时总质量趋于 $|q|/\sqrt{G}$，有限，且与裸质量无关 | Arnowitt, Deser & Misner 1960 |
| 经典理论算出普朗克尺度 | 对电子，$q/\sqrt G \approx 1.9\times10^{-6}$ g $=\sqrt\alpha\, m_{\rm Pl}$；对应长度 $\sqrt G\, q/c^2 \approx 1.4\times10^{-34}$ cm $=\sqrt\alpha\,\ell_{\rm Pl}$。全程没有 ħ | 由 ADM 结果推出 |
| 电子是裸奇点 | 电子的 $\xi = q/(\sqrt G\, m) \approx 2\times10^{21}$，Reissner–Nordström 解在此情形下没有视界 | 经典 GR |
| **runaway 在完整 GR 中仍在** | Robinson–Trautman Einstein–Maxwell 时空的演化方程有 runaway 解；RN 附近线性化，ℓ=1 增长率**恰为** $3M/2Q^2 = 1/\tau_0$，与 $Q/M$ 无关 | Newman & Posadas 1969；Lun & Chow 1994 |
| 视界是选解条件 | Lun–Chow **假设**过去表观视界存在，以此挑选物理解（原文："a global condition selecting the physically meaningful solutions"），推测它排除 runaway，但没有证明 | Lun & Chow 1994 |
| 该模型不适定 | 带电 RT 方程所有 ℓ≥1 模式不稳定，增长率随 ℓ 无界，初值问题不适定 | Kozameh, Kreiss & Reula 2008 |
| **公开问题** | 从 Einstein–Maxwell 渐近场推出慢速的 Abraham–Lorentz 方程；引力能否稳定运动，原文："to our understanding, the answer is unknown" | Kozameh, Newman 等 2009；Adamo, Kozameh & Newman 2009/2012；Newman 2011 |
| 超极端 RN 的不稳定性 | 按曲率条件选边界，$Q>M$ 的 RN 对**所有 ℓ≥2** 不稳定；ℓ=1 不在其中 | Dotti, Gleiser & Pullin 2007；Dotti & Gleiser 2010 |
| 超极端 RN 的"沉默"奇点 | 按能量有限选边界（Friedrichs 扩张），ℓ≥2 没有增长模式；Dotti–Gleiser 的模式能量发散。ℓ=1 被当作纯规范剔除，电荷固定静止，没有运动方程 | Delucchi 2026 |
| 带电壳层的径向坍缩 | 带电薄壳能坍缩成裸奇点，**当且仅当**壳的物质能量密度为负 | Boulware 1973 |
| 弯曲时空自力的"尾巴" | 自力依赖整段过去历史，因为弯曲时空中波不只沿光锥传播 | DeWitt & Brehme 1960；Hobbs 1968 |
| 黑洞推开静止电荷 | 静止在 Schwarzschild 黑洞外的电荷受到排斥性自力 | Smith & Will 1980 |
| 电磁与引力自力耦合 | 小带电物体在电真空背景中的耦合自力，一阶已有完整方程；未讨论 runaway | Zimmerman & Poisson 2014；Linz, Friedman & Wiseman 2014 |

**空白**（三轮检索未找到）：
- 带电有尺寸物体在完整 Einstein–Maxwell 理论中的**平移**运动方程；
- 把 Levine–Moniz–Sharp 判据推广到含引力；
- 超极端 RN 的 ℓ=1 扇区（让点电荷能运动）。

### 1.3 两个质点

| 现象 | 说明 | 来源 |
|---|---|---|
| 经典氢原子塌缩 | 电子绕核辐射，约 $10^{-11}$ s 落入 | 经典 |
| 时间对称电动力学的稳定圆轨道 | 一半推迟、一半超前的相互作用下，两电荷可永久绕转而不辐射 | Schild 1963 |
| 两体问题是延迟微分方程 | 需要整段历史，而非初始条件 | 经典；TheoPhysEval 题 1e |
| 静态多体平衡 | 每个粒子满足 $q = \sqrt G\, m$ 时，电斥力与引力精确抵消，任意摆放均可静止 | Majumdar 1947；Papapetrou 1947 |
| 极端黑洞的慢速散射 | 从静态平衡出发让两个极端带电黑洞运动 | Ferrell & Eardley 1987 |

### 1.4 一条贯穿的线索：选解条件

三处文献里，runaway 都没有被方程本身消除，而是被一个**额外的选解条件**去掉：

| 场合 | 选解条件 |
|---|---|
| 平直时空 LD 方程 | Dirac 的渐近条件（加速度在远未来趋于零） |
| Robinson–Trautman Einstein–Maxwell | 要求过去表观视界存在（Lun–Chow） |
| 超极端 RN 裸奇点 | 奇点处的边界条件：能量有限（Delucchi）vs 曲率条件（Dotti–Gleiser） |

本项目的一个副产品，可能就是回答：这些条件是同一个物理原则的不同写法，还是各自独立的人为选择。有尺寸的模型（第 2 步）不需要任何选解条件，正好可以作为对照。

---

## 2. 一个必须先面对的数量级问题

**在电子的参数下，引力弱到几乎不存在。**

$$\frac{G m_e^2}{e^2} \approx 2.4\times10^{-43}$$

如果把引力当成小修正逐阶加入，它不可能改变 runaway 这种量级的行为。

ADM 的惊喜之所以出现，是因为它是**非微扰**的：电荷被压缩到 $\sqrt G\, q/c^2 \approx 10^{-34}$ cm 时，引力自能与电磁自能同量级。但在那个极限下，总质量变成 $q/\sqrt G$，比电子重 $10^{21}$ 倍。

**静态估算：对电子，引力不起作用。** 在 ADM 球壳模型里把总质量固定为电子质量 $m=m_e$，反解裸质量（$c=1$）：

$$m_0 = m_e - \frac{q^2}{2R} + \frac{G m_e^2}{2R}$$

引力项与电磁项之比为 $Gm_e^2/q^2\sim10^{-43}$。所以 $m_0<0$ 的条件仍是 $R\lesssim r_e/2$，引力对这个临界半径的相对修正只有 $10^{-43}$。ADM 的"点极限下总质量有限"靠的是总质量本身趋于 $|q|/\sqrt G$（即 $\xi\to1$）；只要质量钉在 $m_e$，这条路就不通。电子要么半径 $\sim r_e$，要么裸质量为负。

无论第 3 节猜想 A 用哪种判据（裸质量或尺寸），对电子的结论都一样，两者只差 $10^{-43}$。**真正开放的问题在 $\xi>1$ 一侧、$\xi$ 不太大的区域**（第 3.1 节）。

---

## 3. 核心猜想

### 3.1 猜想 A：GR 中的 runaway 判据

平直时空中，Levine–Moniz–Sharp 的判据可以有两种等价的说法：尺寸大于经典半径，或裸质量为正。**加入引力后，这两种说法不再等价。** ADM 球壳中（$G=c=1$）：

| 判据 | runaway 的条件 |
|---|---|
| 裸质量判据 | $m_0 < 0 \iff R < R_{m_0} = \dfrac{q^2 - m^2}{2m}$ |
| 尺寸判据（LMS 原文的字面推广） | $R < R_{\rm cl} = \dfrac{q^2}{2m}$ |

**$\xi \le 1$：两种判据一致，答案几乎是平凡的。** 静态球壳必须在自己的视界外面，
$$R > r_+ = m + \sqrt{m^2 - q^2} \;\ge\; m \;\ge\; \frac{q^2}{m} \;>\; R_{\rm cl} \;\ge\; R_{m_0}.$$
任何没有塌成黑洞的球壳，都自动比经典半径大，裸质量也自动为正。两种判据都说没有 runaway。换句话说，**$\xi\le1$ 的物体无法被压到经典半径以内而不形成视界。** 这正是 Lun–Chow 视界条件背后的物理，只是这里不需要假设视界，它是自动的。（前提是 LMS 判据的某种形式在 GR 中成立。）

**$\xi > 1$：存在一条两种判据预言相反的带。** 没有视界，球壳可以有任意半径。在
$$R_{m_0} = \frac{q^2 - m^2}{2m} \;<\; R \;<\; \frac{q^2}{2m} = R_{\rm cl}$$
这条带里，裸质量判据说没有 runaway，尺寸判据说有。带宽为 $m/2$（即 $Gm/2c^2$），与引力半径同量级：$\xi\to1^+$ 时它占据 $0<R<m/2$ 的全部；$\xi\gg1$ 时它相对于 $R_{\rm cl}$ 的比例为 $1/\xi^2$，对电子只有 $10^{-43}$。

> **猜想 A**
> 在 Einstein–Maxwell 理论中，正确的判据是**裸质量判据**：一个有尺寸的带电物体出现 runaway，当且仅当它的裸质量为负。因此 $\xi\le1$ 的物体没有 runaway；$\xi>1$ 的物体只在 $R<R_{m_0}$ 时有，而在上述判决带中没有。

**支持的证据**：Boulware 1973 在完整 GR 里得到了这个判据的**径向**版本：带电薄壳坍缩成裸奇点，当且仅当壳的物质能量密度为负。但径向运动没有偶极辐射（Birkhoff 定理），而 runaway 是偶极（平移）现象，所以这只是旁证。

**判决方法**：第 2 步的 ADM 球壳**平移**动力学。取 $\xi$ 略大于 1、$R$ 落在判决带内，检查 runaway 根是否存在。

**已撤回的版本**（演变过程见调研记录第 1.4、2.7 节）：
- "引力从根上消除 runaway"：Robinson–Trautman 模型中 runaway 原样存在。
- "视界消除 runaway"：Lun–Chow 的视界是假设出来的选解条件，消除 runaway 也没有证明。

### 3.2 猜想 B：狄拉克选择题的答案

> **猜想 B**
> Eliezer 型的 $B^\mu$ 推广（在 $P_{\rm mech}$ 中加入 $P_2\ddot u$、$P_4 u^{(4)}$ 等项），降阶后在 $\tau_0$ 一阶全部等价于 LL；从二阶开始出现的差别，对应粒子的内部结构（如电极化率），而不是同一个点电荷的不同运动定律。任何有尺寸的电荷取点极限，都会选出 $P_{2i}=0$，也就是 LD。

一维非相对论下的初步计算支持猜想 B：

$$m\dot v + P_2\dddot v = F + m\tau_0\ddot v
\;\Rightarrow\;
m\dot v = F + \tau_0\dot F + \left(\tau_0^2 - \frac{P_2}{m}\right)\ddot F + \cdots$$

$P_2$ 从二阶开始进入，形式上等价于一个 $\propto E^2$ 的极化率项。相对论情形尚未验证。

**文献中的部分支持**：
- Rowe 1978 的量纲论证：经典电磁学没有长度常数，纯电磁的场张量不能含 $e^2/m$，所以不可能产生 Eliezer 型的高阶项。换句话说，**这类项只能来自粒子本身，也就是内部结构**。
- Gal'tsov & Spirin 2006：Schott 项来自束缚场动量，不能作为附加项随手加进机械动量。

**尚未覆盖**：Rowe 的 $p=mv$ 是假设，不是推出来的；他没有讨论点极限、降阶或与 LL 的比较。猜想 B 的核心论断仍然要自己证明。

**风险**：它可能是世界线有效场论圈内的常识（第 6 节）。定位为独立短文，不作为主结论。

### 3.3 备注：自加速与宇宙监督

第一轮检索时提出过一个问题：runaway 只在没有视界时存在，没有视界就是裸奇点，那么**自加速和宇宙监督是不是同一件事的两面？**

按现有证据**不成立**：
- Robinson–Trautman 的 runaway 增长率与 $Q/M$ 无关，次极端也有；
- 已知的超极端 RN 不稳定性在 ℓ≥2，ℓ=1 没有，增长率形式也与 $3M/2Q^2$ 不同。

仍有一个弱版本成立，但几乎是平凡的：$\xi\le1$ 的物体无法被压到经典半径以内而不形成视界（第 3.1 节），所以在这一侧，"有视界"与"没有 runaway"自动同时成立。它来自 $Q^2/M$ 与 $M$ 的比较，不需要更深的原因。超极端一侧是否有更多内容，第 2′ 步会顺带给出部分答案。

---

## 4. 路线图

### 总览

| 步骤 | 内容 | 状态 / 决策点 |
|---|---|---|
| 第 0 步 | 文献核查 | **已完成三轮**：核心问题未见有人回答 |
| 第 1 步 | 平直时空：把狄拉克方程彻底吃透 | 能否复现 Levine–Moniz–Sharp；猜想 B 是否成立 |
| 第 2 步 | 加入引力：ADM 带电球壳 | **哪个判据正确：$\xi>1$ 判决带中有没有 runaway（最关键）** |
| 第 2′ 步 | 并行支线：超极端 RN 的 ℓ=1 扇区 | runaway 是否出现；是否取决于自伴扩张的选择 |
| 第 3 步 | 确定适用范围 | 对电子：静态部分已有答案（不成立）；只需动力学确认 |
| 第 4 步 | 独立验证（从第 2 步起并行） | 有没有专家签字 |
| 第 5 步 | 两个质点（第二阶段） | — |

---

### 第 0 步：文献核查（已完成三轮）

过程见调研记录。结论：

1. **核心问题是一个被明确提出、没人回答的公开问题。** Kozameh–Newman 等从 Einstein–Maxwell 渐近场推出了慢速 Abraham–Lorentz 方程，原话说引力能否稳定运动"the answer is unknown"。
2. **runaway 在完整 GR 中仍然存在**（Robinson–Trautman，增长率与 LD 完全相同）。去掉它的只有人为的选解条件（视界），而且没有证明。该模型本身不适定，不适合用来回答问题。
3. **三个空白**：带电有尺寸物体在 GR 中的平移运动方程；Levine–Moniz–Sharp 判据在 GR 中的推广；超极端 RN 的 ℓ=1 扇区。
4. **猜想 B 的风险仍在**：没找到"Eliezer 推广等价于内部结构"的明确陈述，但 Rowe 1978 已经很接近。

**尚未完成**（不阻塞第 1 步）：
- Boulware 1973、Levine–Moniz–Sharp 1977 的全文：都在付费墙后，需要机构访问。
- ℓ=1 工具文献：Zerilli 1974、Burko 1995/1999、Kodama–Ishibashi、Moncrief。第 2′ 步启动前读。
- Boulware 1973 的引用网络：APS 显示被引 116 次，查有没有人做过平移运动。

---

### 第 1 步：平直时空，把狄拉克方程彻底吃透

这是地基。第 2 步的所有结果，在 $G \to 0$ 时都必须回到这里。它也是检验 AI 推导可靠性的标尺。

**1a. 两种方法独立重推 LD**

- Dirac 1938 的管面积分法。
- 推迟坐标法（Poisson 的处理）。
- 要求：两种方法给出同一个方程；动量通量中的发散项正比于 $a^\mu$，能被质量重整化吸收。
- 顺带复现 Rowe 1975 的歧义：换成光锥端面，得到变质量方程。

**1b. 精确刻画 runaway**

- 一维情形求出特征根，确认 runaway 来自 $e^{\tau/\tau_0}$ 那一支。
- 写出非 runaway 解的积分形式 $m\dot v(t)=\int_0^\infty e^{-s}F(t+\tau_0 s)\,ds$，明确预加速的来源。
- 降阶级数 $\sum \tau_0^n F^{(n)}$：在什么条件下收敛，什么条件下只是渐近级数。

**1c. 有尺寸电荷模型**

- 带电球壳的自力（微分–差分方程），复现 Levine–Moniz–Sharp：尺寸大于经典半径时无 runaway。
- 明确写出平直时空中"尺寸判据"与"裸质量判据"的等价性，以及等价性在哪一步用到了 $E_{\rm self} = q^2/2R$。**这一步决定第 2 步怎么区分两种判据。**
- 写成方便推广到 GR 的形式：第 2 步直接搬进 Einstein–Maxwell。

**1d. 狄拉克的选择题（猜想 B）**

- 写出 Eliezer 型 $P_{\rm mech}$ 的一般形式（2006 年论文式 (1.3.3)）。
- 相对论情形下逐阶降阶到 $\tau_0^2$，与"LD + 极化率世界线作用量"逐项比较。
- 判断哪些 $P_{2i}$ 可以通过重新定义世界线消掉，哪些不能。
- 用 1c 的有尺寸模型检查：点极限下 $P_{2i}$ 是否自动为零。
- 把 Rowe 1978 的倾角依赖论证推广到粒子张量：$\int(vw+wv)\delta$ 型的世界线项是否必然要求偶极型结构项（调研记录第 3.1 节）。

**1e. 练手题：Robinson–Trautman 的 ℓ=1 增长率**

- 在 RN 背景下写出 RT–Maxwell 方程的线性化，推出 Lun–Chow 的特征值
  $$\lambda_{\ell,\pm} = \frac{3M}{8Q^2}\,\ell(\ell+1)\left[1 \pm \sqrt{1 + \frac{4Q^2(\ell+2)(\ell-1)}{9M^2}}\right]$$
  和 $\lambda_{1,+}=3M/2Q^2$。
- 答案已知，用来检验 AI 推导流程（两次独立运行 + SymPy）。

**产出**

- 一套经过两种方法、两次独立 AI 运行和 SymPy 校验的推导。
- 可复用的代码：LD、LL、降阶级数、有尺寸电荷模型的数值求解器。
- 猜想 B 的结论（无论成立与否），本身可以写成一篇短文。

**决策**

- 无法复现 Levine–Moniz–Sharp → 先查清原因，不进入第 2 步。
- 猜想 B 成立 → 记录，作为独立结果。

---

### 第 2 步：加入引力，ADM 带电球壳

**为什么是这个模型**：它是有尺寸的，不需要任何选解条件（第 1.4 节）；它直接回答 Kozameh–Newman 留下的公开问题；它能检验 Levine–Moniz–Sharp 判据在 GR 中的形式，这是 Newman 的渐近方法做不到的。

**模型**：ADM 的带电薄球壳，放入完全的 Einstein–Maxwell 理论。壳内平直，壳外 Reissner–Nordström，壳面满足 Israel 结合条件。

**2a. 复现 ADM 静态结果**

- 总质量满足 $m = m_0 + \dfrac{q^2}{2R} - \dfrac{G m^2}{2R}$（$c=1$），解出
  $m = \dfrac{R}{G}\left[-1+\sqrt{1+\dfrac{2G}{R}\left(m_0+\dfrac{q^2}{2R}\right)}\right]$。
- 验证 $R\to0$ 时 $m\to |q|/\sqrt G$，与 $m_0$ 无关。
- 写出两个临界半径 $R_{m_0}$、$R_{\rm cl}$ 和视界半径 $r_+$，在 $(\xi, R)$ 平面上画出：视界内（不可达）、两种判据都说有 runaway、判决带、两种判据都说没有 runaway 四个区域。
- 注意 $R\to0$ 时 $\xi\to 1$：点极限自动落在极端条件上。

**2b. 让球壳动起来**

- **径向**（热身）：球对称，无辐射，Birkhoff 定理保证壳外仍是 RN。复现 Boulware 1973 的结论，并检查是否有指数增长的模式。注意：径向运动没有偶极辐射，**判决不了 runaway**。
- **平移**：施加一个小的外力，计算球壳中心的运动。这一步破坏球对称，需要 ℓ=1 微扰（壳内平直、壳外 RN 的 ℓ=1 扰动，用结合条件拼接）或数值相对论。
- 两种情况都要和第 1 步 $G\to0$ 的结果对接。

**2c. 核心问题**

1. 在判决带（$\xi>1$，$R_{m_0}<R<R_{\rm cl}$）中，平移运动方程有没有 runaway 特征根？
2. $\xi\to1^+$ 时，判决带占据 $0<R<m/2$，这里引力是强场；runaway 根随 $R$ 怎样变化？
3. 检查 $\xi\le1$、$R>r_+$ 时确实没有 runaway（第 3.1 节的"平凡"结论）。

**决策（整个项目最关键的一个）**

- 裸质量判据成立 → 猜想 A 成立：引力把 runaway 的边界从 $R_{\rm cl}$ 移到 $R_{m_0}$，移动量恰为 $Gm/2c^2$。进入第 3 步。
- 尺寸判据成立 → 引力不移动 runaway 的边界；$\xi\le1$ 时没有 runaway 只是因为视界挡住了。也是一个干净的结果。
- 两者都不对，或结果依赖模型细节（球壳 vs 实心球 vs 尘埃）→ 问题本身不是普适的，需要重新表述。

**技术风险**

- 平移运动需要壳内外 ℓ=1 扰动的拼接，或数值相对论，工作量会明显增加。
- 点极限与极端条件重合，接近 $q=\sqrt G m$ 时几何可能退化，需要特别小心。

---

### 第 2′ 步：并行支线，超极端 Reissner–Nordström 的 ℓ=1 扇区

**问题**：把超极端 RN（$Q>M$）当作一个能运动的点电荷，在 ℓ=1 扇区推出它的世界线方程。看会不会出现 $3M/2Q^2$ 的 runaway，以及 Delucchi 的能量有限条件（Friedrichs 扩张）能否把它排除。

**为什么值得做**

- 这是点粒子模型，不需要内部结构，和第 2 步互补。
- Delucchi 2026 的论证说 ℓ=1 是纯规范，这在 Einstein–Maxwell 中不成立：电磁偶极自由度是真实的（调研记录第 3.2 节）。
- 静态背景和 ℓ≥2 的自伴扩张分析已有现成结果可以直接借用。
- 如果 runaway 只在非 Friedrichs 扩张下出现，那"runaway 是否存在"就等价于"选哪个自伴扩张"。这与 Dirac 的渐近条件同构（第 1.4 节），本身就是一个结果。

**步骤与工作量估计**

| 子步骤 | 内容 | 工作量 |
|---|---|---|
| 2′a | 轴向 ℓ=1：电磁一支的势 $V = f(2/r^2 + 4Q^2/r^4)$（已手算，未与文献对照），套用 Hardy/Friedrichs 检查谱 | 约一天 |
| 2′b | 极性 ℓ=1 真空主方程（$Q>M$）：改造 Burko 1995 或 Kodama–Ishibashi，检查与轴向是否等谱 | 几天 |
| 2′c | 加入世界线：把 ℓ=1 规范矢量 $\xi(t)$ 当作集体坐标，由 ℓ=1 约束或无穷远处 ADM 动量平衡读出它的方程 | 几周 |
| 2′d | 求各自伴扩张下的谱或散射矩阵极点，检验 $\kappa = 3M/2Q^2$ | 视 2′c 而定 |

**决策**

- runaway 在 Friedrichs 扩张下也出现 → 点电荷模型中引力不能消除 runaway，与 Delucchi 的"沉默奇点"图景冲突，值得写。
- runaway 只在非 Friedrichs 扩张下出现 → "选解条件"线索的一个具体实例。
- 2′b 发现极性 ℓ=1 没有闭合形式，或 2′c 无法定义世界线 → 降为备注，不影响主线。

---

### 第 3 步：确定适用范围

**目标**：回答第 2 节的数量级问题。

**静态部分已有答案**：第 2 节的估算表明，在 ADM 球壳模型里固定 $m=m_e$ 时，引力对 runaway 临界半径的修正只有 $10^{-43}$ 量级，无论用哪种判据。**结论主语预定为"带电天体"，而不是电子。**

动力学部分（第 2 步之后）只需确认没有静态估算看不到的效应：

- 用 $\xi = q/(\sqrt G\, m)$ 给出 runaway 消失的条件。
- 区分三种情况：$\xi < 1$（亚极端，可以形成黑洞）、$\xi \approx 1$（极端）、$\xi \gg 1$（超极端，电子 $\xi\sim10^{21}$）。

**结论主语**

- 第 2 步裸质量判据成立 → "引力把带电物体自加速的临界尺寸缩小了一个引力半径；对接近极端的带电天体这是决定性的，对电子毫无影响"。
- 第 2 步尺寸判据成立 → "引力本身不改变自加速的判据，只是视界让 $\xi\le1$ 的物体够不着病态区"。
- 量子一侧的独立判据（第 0.1 节）同样指向带电天体：只有电荷数 $N\gtrsim12$ 时经典处理才自洽。

---

### 第 4 步：独立验证（从第 2 步起并行）

**原则**：没有领域内专家独立核对，任何结果都不对外讲。

**做法**

- 第 2 步一出现苗头，就把推导和代码整理成可独立运行的形式。
- 联系做 Reissner–Nordström、自力、带电壳层的研究者，请他们核对关键一步。
- 候选方向：
  - Newman、Kozameh 一系（提出公开问题的人）；
  - Delucchi；Dotti、Gleiser（第 2′ 步）；
  - Gralla、Harte、Wald（严格的点粒子极限）；
  - Capra 自力社区；
  - TheoPhysEval 中修正电动力学题目的出题人（通讯作者 Suvrat Raju）。

**产出**：至少一位专家书面确认关键结果。

---

### 第 5 步：两个质点（第二阶段）

在第 2、3 步有结论之后再展开。

**5a. 平直时空基准**

- 带推迟相互作用的两体问题（延迟微分方程，需要历史条件）。
- 经典氢原子塌缩的完整数值复现，作为代码校验。
- 对照 Schild 的时间对称圆轨道。

**5b. 加入引力**

- 从 Majumdar–Papapetrou 静态平衡（$q=\sqrt G m$）出发，让两个物体动起来。
- 参考 Ferrell–Eardley 的极端黑洞慢速散射。
- 问题：偏离极端条件时，辐射损失与引力束缚如何竞争？有没有经典图景下的新行为？

---

## 5. AI 的使用规范

**分工**

- 人负责：提出问题、在每个决策点做判断、找人验证、决定对外讲什么。
- AI 负责：推导、符号计算、数值实现、文献检索与整理。

**验证规则（每个关键结果必须同时满足）**

1. 至少两种独立方法推出同一结果。
2. AI 两次独立运行（互不看对方输出）结果一致。
3. 解析结果配一个数值检验。
4. 有已知极限的（$G\to0$、非相对论、球对称），必须回到已知结果。

**文献规则**（第 0 步的教训）

- 摘要不等于结论。第一轮基于摘要的两条判断（"视界消除 runaway""runaway ⇔ 超极端"），读全文后都被推翻。
- 每条引用标注核实程度：【全文】【摘要】【二手】。
- 转述别人的结论时，与原文措辞对照。例如 Levine–Moniz–Sharp 原文说的是"尺寸"，不是"裸质量"。

**记录**

- 所有推导、代码、AI 会话保存在仓库中，按步骤分目录。
- 每一步结束写一页小结：做了什么、结论、是否通过决策点。
- 这份记录本身就是"AI 如何参与物理研究"的一手材料。

---

## 6. 风险清单

| 风险 | 可能性 | 影响 | 应对 |
|---|---|---|---|
| 核心问题已有人做过 | 低（三轮未找到；Kozameh–Newman 等明确称答案未知） | 高 | 第 1、2 步期间持续跟踪新文献；查 Boulware 1973 的引用网络 |
| 引力对电子无效 | 已基本确定（第 2 节） | 中 | 结论主语预定为带电天体 |
| 尺寸判据才是对的（猜想 A 不成立） | 中 | 中 | 也是干净的结果，可以写 |
| $\xi\le1$ 的结论被认为平凡 | 高 | 中 | 主动说明它的平凡性（第 3.1 节）；把卖点放在 $\xi>1$ 判决带 |
| 平移动力学需要数值相对论 | 中 | 中 | 先做径向，再做 ℓ=1 微扰拼接，最后才考虑数值相对论 |
| 第 2′ 步的极性 ℓ=1 没有可用的主方程 | 中 | 低 | 支线，不阻塞主线 |
| 猜想 B 是 EFT 圈内常识 | 中 | 低 | 作为独立短文，不作为主结论 |
| AI 推导出现难以察觉的错误 | 中 | 高 | 第 5 节的验证规则；第 1e 步先用已知答案练手 |
| AI 文献检索误读原文 | 中 | 高 | 第 5 节的文献规则；关键文献人工核对 |
| 找不到专家愿意核对 | 中 | 高 | 从第 2 步起提前联系，不等结果完整 |

---

## 7. 对外发布的条件

只有同时满足以下三条才对外发布：

1. 第 2 步得到明确结果，第 3 步明确了适用范围。
2. 至少一位领域专家书面确认。
3. 结论能用一句外行听得懂的话讲清楚，不出现任何需要解释的术语。

新闻稿模板待写（放在 `docs/` 下）。模板中的"独立验证"和"局限"两段是硬条件，空着不发。

---

## 8. 参考文献

> 标注：**【已核】** 出版信息经 Crossref、INSPIRE 或 arXiv 核实；**【待核】** 尚未核实。内容的核实程度（全文 / 摘要 / 二手）见调研记录。

### 8.1 平直时空：LD 方程及其推导

- P. A. M. Dirac, *Classical theory of radiating electrons*, *Proc. R. Soc. A* 167, 148–169 (1938).【已核】
- C. J. Eliezer, *The interaction of electrons and an electromagnetic field*, *Rev. Mod. Phys.* 19, 147–184 (1947).【已核】
- A. Schild, *Electromagnetic two-body problem*, *Phys. Rev.* 131, 2762–2766 (1963).【已核】
- E. J. Moniz, D. H. Sharp, *Absence of runaways and divergent self-mass in nonrelativistic quantum electrodynamics*, *Phys. Rev. D* 10, 1133–1136 (1974).【已核】
- E. G. P. Rowe, *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975).【已核，本地有】
- R. Tabensky, *Electrodynamics and the electron equation of motion*, *Phys. Rev. D* 13, 267–273 (1976).【已核】
- H. Levine, E. J. Moniz, D. H. Sharp, *Motion of extended charges in classical electrodynamics*, *Am. J. Phys.* 45, 75–78 (1977).【已核】
- E. J. Moniz, D. H. Sharp, *Radiation reaction in nonrelativistic quantum electrodynamics*, *Phys. Rev. D* 15, 2850–2865 (1977).【已核】
- E. G. P. Rowe, *Structure of the energy tensor in the classical electrodynamics of point particles*, *Phys. Rev. D* 18, 3639 (1978).【已核，本地有】
- E. G. P. Rowe, *Angular momentum radiated by classical point particles*, *Phys. Rev. D* 20, 580 (1979).【已核，本地有】
- H. Spohn, *The critical manifold of the Lorentz–Dirac equation*, *Europhys. Lett.* 50, 287–292 (2000), arXiv:physics/9911027.【已核】
- H. Spohn, *Dynamics of Charged Particles and Their Radiation Field* (Cambridge University Press, 2004).【待核】
- D. V. Gal'tsov, P. Spirin, *Radiation reaction reexamined: bound momentum and Schott term*, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121.【已核】
- A. D. Yaghjian, *Relativistic Dynamics of a Charged Sphere: Updating the Lorentz–Abraham Model* (Springer; 1992, 2006, 2022).【已核】
- S. E. Gralla, A. I. Harte, R. M. Wald, *Rigorous derivation of electromagnetic self-force*, *Phys. Rev. D* 80, 024031 (2009), arXiv:0905.2391.【已核】
- O. Birnholtz, arXiv:1410.5871.【待核】
- B. P. Kosyakov, *Phys. Rep.* (2019), arXiv:1812.03290.【待核】
- W. D. Goldberger, I. Z. Rothstein, *Phys. Rev. D* 73, 104029 (2006).【待核】
- A. Ramesh Chandra *et al.*, *TheoPhysEval: Assessing AI agents at autonomous theoretical physics problem solving*, arXiv:2610.16733 (2026)（题 1d、1e）.【已核】

### 8.2 引力：自能、自力与渐近方法

- B. S. DeWitt, R. W. Brehme, *Radiation damping in a gravitational field*, *Ann. Phys.* 9, 220–259 (1960).【已核】
- R. Arnowitt, S. Deser, C. W. Misner, *Gravitational-electromagnetic coupling and the classical self-energy problem*, *Phys. Rev.* 120, 313–320 (1960).【已核】
- J. M. Hobbs, *A vierbein formalism of radiation damping*, *Ann. Phys.* 47, 141–165 (1968).【已核】
- A. G. Smith, C. M. Will, *Force on a static charge outside a Schwarzschild black hole*, *Phys. Rev. D* 22, 1276–1284 (1980).【已核】
- J. C. Herrera, *Equation of motion in classical electrodynamics. II*, *Phys. Rev. D* 21, 384 (1980).【已核，本地有】
- J. F. Donoghue, *General relativity as an effective field theory: the leading quantum corrections*, *Phys. Rev. D* 50, 3874–3888 (1994), arXiv:gr-qc/9405057.【已核】
- E. Poisson, A. Pound, I. Vega, *The motion of point particles in curved spacetime*, *Living Rev. Relativ.* 14, 7 (2011), arXiv:1102.0529.【已核】
- P. Zimmerman, E. Poisson, *Gravitational self-force in nonvacuum spacetimes*, *Phys. Rev. D* 90, 084030 (2014), arXiv:1406.5111.【已核】
- T. M. Linz, J. L. Friedman, A. G. Wiseman, *Self force on an accelerated particle*, *Phys. Rev. D* 90, 024064 (2014), arXiv:1404.7039.【已核】
- T. M. Linz, J. L. Friedman, A. G. Wiseman, *Phys. Rev. D* 90, 084031 (2014), arXiv:1406.5112.【已核】
- C. N. Kozameh, E. T. Newman, J. G. Ortega, G. Silva-Ortigoza, *Electrodynamic radiation reaction and general relativity*, *Gen. Rel. Grav.* 41, 1139 (2009), arXiv:0803.3921.【已核】
- T. M. Adamo, C. N. Kozameh, E. T. Newman, *Living Rev. Relativ.* 12, 6 (2009)；更新版 15, 1 (2012), arXiv:0906.2155.【已核】
- E. T. Newman, *Newton's 2nd law, radiation reaction & type II Einstein–Maxwell fields*, *Class. Quantum Grav.* 28, 245003 (2011), arXiv:1109.4106.【已核】
- E. T. Newman, *Classical mechanics via general relativity and Maxwell's theory: a bit of magic*, *Gen. Rel. Grav.* 50, 111 (2018).【已核】

### 8.3 Robinson–Trautman Einstein–Maxwell

- E. T. Newman, R. Posadas, *Motion and structure of singularities in general relativity*, *Phys. Rev.* 187, 1784–1791 (1969)；II, *J. Math. Phys.* 12, 2319–2327 (1971).【已核】
- A. W. C. Lun, E. W. M. Chow, *The role of the apparent horizon in the evolution of Robinson–Trautman Einstein–Maxwell space-times*, arXiv:gr-qc/9409024 (1994)；收入 *Confronting the Infinite*, World Scientific (1995), pp. 279–292.【arXiv 已核，文集信息为二手】
- C. N. Kozameh, H.-O. Kreiss, O. Reula, *Class. Quantum Grav.* 25, 025004 (2008), arXiv:0708.1933.【已核】
- P. T. Chruściel, arXiv:2112.02126.【待核】
- T. Tahamtan, *Phys. Rev. D* 103, 064052 (2021), arXiv:2010.01689.【已核】
- T. Tahamtan, A. Flores-Alfonso, O. Svitek, *Phys. Rev. D* 108, 124076 (2023), arXiv:2311.03110.【已核】

### 8.4 超极端 Reissner–Nordström 与点电荷时空

- G. Dotti, R. J. Gleiser, J. Pullin, *Instability of charged and rotating naked singularities*, *Phys. Lett. B* 644, 289 (2007), arXiv:gr-qc/0607052.【已核】
- G. Dotti, R. J. Gleiser, *Gravitational instability of the inner static region of a Reissner–Nordström black hole*, *Class. Quantum Grav.* 27, 185007 (2010), arXiv:1001.0152.【已核】
- Delucchi, *Unitary and finite self-energy of a massive classical point charge as a naked point singularity*, *Class. Quantum Grav.* 43, 195008 (2026), arXiv:2602.19015.【已核】
- A. S. Tahvildar-Zadeh, *On the static spacetime of a single point charge*, *Rev. Math. Phys.* 23, 309–346 (2011), arXiv:1012.1400.【已核】
- A. Burtscher, M. K.-H. Kiessling, A. S. Tahvildar-Zadeh, *Weak second Bianchi identity for static, spherically symmetric spacetimes with timelike singularities*, *Class. Quantum Grav.* 38, 185001 (2021), arXiv:1901.00813.【已核】
- 第 2′ 步要读的工具文献（均【待核】）：F. J. Zerilli, *Phys. Rev. D* 9, 860 (1974)；L. M. Burko, arXiv:gr-qc/9507012、gr-qc/9902002；H. Kodama, A. Ishibashi（ℓ=0、1 低多极处理）；V. Moncrief (1974–75)；D. Bini, A. Geralico, R. Ruffini, arXiv:gr-qc/0609041、gr-qc/0608139；arXiv:2010.12984 附录 C。

### 8.5 带电壳层

- W. Israel, *Singular hypersurfaces and thin shells in general relativity*, *Nuovo Cim. B* 44, 1–14 (1966)；erratum 48, 463 (1967).【已核】
- V. de la Cruz, W. Israel, *Gravitational bounce*, *Nuovo Cim. A* 51, 744–760 (1967).【已核】
- K. Kuchař, *Charged shells in general relativity and their gravitational collapse*, *Czech. J. Phys. B* 18, 435–463 (1968).【已核】
- D. G. Boulware, *Naked singularities, thin shells, and the Reissner–Nordström metric*, *Phys. Rev. D* 8, 2363–2368 (1973).【已核】
- G. Dillon, *The gravitational self-energy of a spherical shell*, *Eur. Phys. J. Plus* 128, 72 (2013), arXiv:1303.2706.【已核】
- *Electromagnetic Radiation Reaction in General Relativity*，博士论文（1988），可能为 N. O'Donnell, Cardiff University.【待核】

### 8.6 两个质点

- S. D. Majumdar, *Phys. Rev.* 72, 390 (1947)；A. Papapetrou, *Proc. R. Irish Acad. A* 51, 191 (1947).【待核】
- R. C. Ferrell, D. M. Eardley, *Slow-motion scattering and coalescence of maximally charged black holes*, *Phys. Rev. Lett.* 59, 1617–1620 (1987).【已核】

---

## 9. 版本记录

旧版本保存在 git 历史中，可用 `git show <commit>:docs/roadmap.md` 查看。

| 版本 | commit | 主要内容 |
|---|---|---|
| v1 | `3f3549e` | 初版路线图；加入"对电子，引力修正只有 $10^{-43}$"的静态估算 |
| v1 + 第一轮 | `722f1a7` | 附加第一轮文献核查（摘要级）；猜想 A′：视界或尺寸消除 runaway |
| v1 + 第二轮 | `081b942` | 附加第二轮（读全文）；视界版本撤回；猜想 A″：裸质量判据；新增超极端 RN ℓ=1 方向；改正 arXiv 编号 |
| v1 + QED 说明 | `2ef7ba8` | 第 0 节加入"为什么用经典理论，而不用 QED" |
| v2 | `db4ed75` | 正文整体重写为现行版；文献核查过程移到 `docs/literature-review.md`；猜想 A 改为区分裸质量判据与尺寸判据（Levine–Moniz–Sharp 原文是尺寸判据）；猜想 B 加入 Rowe 1978；第 2′ 步成为正式支线；新增"选解条件"线索；参考文献合并去重并标注核实状态 |
| **v2.1** | 本次 | 第 0.1 节加入量子一侧的独立判据（经典半径/康普顿波长 $=Q^2/\hbar c$）；改正第 3.1 节：$\xi\le1$ 时两种判据一致且结论几乎平凡，真正的判决带在 $\xi>1$、$R_{m_0}<R<R_{\rm cl}$；指出径向运动判决不了 runaway；相应修改第 0.3、2、3.3 节和第 2、3 步 |
