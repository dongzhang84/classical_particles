# 经典带电物体在引力下的运动：研究路线图

> 版本：v3 · 2026-10-23 · 草稿
> 状态：第 0 步（文献核查）已完成四轮，下一步是第 1 步。所有"结果"仍为待检验的猜想。
> 文献核查的过程记录见 [`docs/literature-review.md`](literature-review.md)（下称"调研记录"）。旧版本见第 9 节。

---

## 0. 先把话说清楚

### 0.1 这个项目研究什么、不研究什么

**研究对象：经典带电物体。** 也就是带有大量电荷的宏观物体，比如带电小球、带电薄壳、带电星体。我们问的是：在麦克斯韦电动力学加爱因斯坦广义相对论这套**经典**理论里，这样一个质量为 $m$、电荷为 $q$ 的物体会怎样运动；再推到两个物体。

**不研究：电子、质子以及任何单个基本粒子。** 不是因为它们不重要，而是因为经典理论管不到它们（第 0.2 节）。

**两个核心问题**

1. **狄拉克的选择题。** 狄拉克 1938 年推导带电粒子运动方程时承认，能动量守恒只约束了一部分，"机械动量" $P_{\rm mech}$（即他的 $B^\mu$）怎么写，是他挑了最简单的一种。其他写法是不同的物理，还是同一个物理的不同写法？
2. **引力会不会改变"自加速"出现的条件。** 平直时空里，一个有尺寸的带电物体，只要比自己的"经典半径"大，就不会出现自加速（runaway）。加上引力以后，这个条件会不会变？变成什么？

**一句话目标**（目标，不是结论）：

> 在经典 Maxwell–Einstein 理论中，写出带电物体不出现自加速的精确条件，并算出引力把这个条件改变了多少。

### 0.2 为什么不研究电子和质子

**对任何带一个基本电荷 $e$ 的粒子，量子效应总是比经典的病先出现。** 原因是一个与质量无关的比值：

$$\frac{\text{经典半径}}{\text{康普顿波长}} = \frac{e^2/mc^2}{\hbar/mc} = \frac{e^2}{\hbar c} = \alpha \approx \frac{1}{137}.$$

经典半径是经典电磁自能开始出问题的尺度；康普顿波长是量子力学接管的尺度。前者永远比后者小 137 倍。所以把电子（或质子、μ 子……）往小里压，还没压到经典的病发作的地方，它就已经是量子的了。

**以电子为例，各个尺度从大到小：**

| 尺度 | 大小 | 在这里发生什么 |
|---|---|---|
| 康普顿波长 $\hbar/m_ec$ | $4\times10^{-11}$ cm | **量子力学接管**。比这更小，电子不能当经典粒子看 |
| 经典电子半径 $r_e = e^2/m_ec^2$ | $3\times10^{-13}$ cm | 经典电磁自能等于质量；经典理论的自加速问题出现在这里 |
| 实验上看不到电子结构的上限 | $\lesssim10^{-18}$ cm 量级 | 电子比 $r_e$ 小得多，"经典小球"模型被实验排除 |
| ADM 尺度 $\sqrt G\,e/c^2 = \sqrt\alpha\,\ell_{\rm Pl}$ | $1.4\times10^{-34}$ cm | **引力自能**才与电磁自能同量级 |
| 电子的史瓦西半径 $2Gm_e/c^2$ | $10^{-55}$ cm | 比普朗克长度还小，没有物理意义 |

**由此得到三点：**

1. **电子的自加速问题，是被量子力学解决的，不是被引力。** 在量子电动力学（QED）里没有 runaway（Moniz & Sharp 1974）。
2. **"引力能不能治好电子的自加速"在经典理论里不是一个能问的问题。** 引力要到 $10^{-34}$ cm 才变强，那里需要量子引力，而量子引力理论还不存在。
3. **电子不是黑洞。** 它的电荷远大于质量（$\xi = e/\sqrt G m_e \approx 2\times10^{21}$），经典 GR 里根本形成不了视界。

**质子也一样。** 质子的经典半径 $\sim1.5\times10^{-16}$ cm，康普顿波长 $\sim2\times10^{-14}$ cm，比值同样是 137。（质子另有真实尺寸 $\sim0.84\times10^{-13}$ cm，比它的经典半径大得多，所以即使当经典物体看也没有自加速。）

**什么样的物体经典理论管得到？** 电荷为 $Ne$ 时，上面的比值变成 $N^2\alpha$。**电荷数 $N\gtrsim12$ 时，经典半径就大于康普顿波长，经典处理自洽。** 任何宏观带电物体都远远满足这个条件。

**$10^{-34}$ cm 以下的问题有没有意义？** 有，而且学术界一直在讨论，但都是在量子引力的框架里，不在本项目范围内。相关方向（均【待核】，凭记忆列出）：引力能否作为量子场论的"普适正规化子"（DeWitt 1964；Isham, Salam & Strathdee 1971–72，这是 ADM 想法的量子版本）；QED 加引力时耦合常数的跑动（Robinson & Wilczek 2006 及其后的争论）；弱引力猜想（Arkani-Hamed 等 2006，已核）。

### 0.3 为什么用经典理论

对第 0.1 节的研究对象（宏观带电物体），经典 Maxwell–Einstein 理论就是正确的理论：量子效应可以忽略（$N\gg12$），而引力必须完整、非微扰地保留，因为我们关心的是物体被压到接近自己引力半径时的行为。

QED 可以放在给定的弯曲时空上（弯曲时空量子场论），引力也可以作为低能有效场论量子化（Donoghue 1994）。但能同时、自洽地包含"物体的自场"和"这个自场产生的完整时空弯曲"的理论，目前只有经典 Maxwell–Einstein 理论。

**方法论**：不去修补方程，而是把方程当真。爱因斯坦当年就是这样做的：他没有修改麦克斯韦方程，而是认真对待它和力学之间的矛盾，矛盾自己推出了相对论。

### 0.4 这个判据什么时候真正起作用

自加速的判据只在物体被压到与它的经典半径、引力半径同量级时才起作用。对日常物体，这永远碰不到：

- 一个 1 克的小球，只要多带约 $5\times10^5$ 个电子，它的电荷质量比 $\xi = q/(\sqrt G\,m)$ 就等于 1。所以 $\xi\sim1$ 的宏观物体很容易有。
- 但对 1 克的物体，经典半径和引力半径都在 $10^{-28}$ cm 量级，而真实小球有 1 cm。

**所以这个问题真正相关的，是接近黑洞那样致密的带电天体**：带电星体的坍缩、带电黑洞的形成、带电壳层的动力学。本项目的结论首先是关于经典 Maxwell–Einstein 理论自身的结论，其次是关于致密带电天体的结论。

---

## 1. 已知与未知

本节只列核实过的结论；出处和核实程度见调研记录。

### 1.1 学术界已经解决了什么

| 问题 | 现状 | 怎么解决的 | 来源 |
|---|---|---|---|
| 平直时空中，有尺寸带电物体有没有自加速 | **已解决（精确结果）** | 只用麦克斯韦方程和洛伦兹力，精确求解带电球壳的运动（微分–差分方程）。**尺寸大于经典半径，即裸质量为正时，没有自加速，也没有预加速** | Levine, Moniz & Sharp 1977；Moniz & Sharp 1977；Medina 2006；Steane 2015；D'Andrea, Luty & Verhaaren 2019 |
| LD 方程和 LL 方程是什么 | **已解决（严格近似）** | 当物体远小于运动变化的尺度时，精确方程可以展开；一阶就是 LL 方程。LD 方程的 runaway 是展开写错了形式带进来的 | Spohn 2000；Gralla, Harte & Wald 2009；Poisson, Pound & Vega 2011 |
| 电子等基本粒子的自加速 | **交给量子理论** | QED 中点电子自质量为零，没有 runaway | Moniz & Sharp 1974 |
| 引力对电子的影响 | **可以忽略** | 引力修正是 $Gm_e^2/e^2\sim10^{-43}$ | Steane 2015 专门讨论过并因此舍弃引力 |

"裸质量"指扣掉电磁自能之后物体本身剩下的质量。裸质量为负意味着惯性为负，推它它往回走，这在物理上没有意义。所以"裸质量为正"是一个物理要求。Steane 2015 的说法是："没有任何已知的、裸质量为正的物理对象，其运动方程有 runaway 解。"

### 1.2 加上广义相对论：已知的事实

| 现象 | 说明 | 来源 |
|---|---|---|
| 引力治好了自能发散 | 带电球壳在 Einstein–Maxwell 理论中，半径 → 0 时总质量趋于 $|q|/\sqrt G$，有限，且与裸质量无关 | Arnowitt, Deser & Misner 1960 |
| runaway 在完整 GR 中仍然存在 | Robinson–Trautman Einstein–Maxwell 时空：RN 附近线性化，偶极模式的增长率**恰为** $3M/2Q^2 = 1/\tau_0$，与 $Q/M$ 无关 | Newman & Posadas 1969；Lun & Chow 1994 |
| 视界只是一个选解条件 | Lun–Chow 假设过去表观视界存在，以此挑选物理解，但没有证明它能排除 runaway；而且该模型的初值问题不适定 | Lun & Chow 1994；Kozameh, Kreiss & Reula 2008 |
| GR 中的点电荷：黑洞或负裸质量 | 裸质量为正的点电荷会被视界包住；没有视界的点电荷只能取负裸质量 | Kiessling & Tahvildar-Zadeh 2019（引 Geroch–Traschen）；Ponce de León 2004 |
| 带电壳层的径向坍缩 | 带电薄壳能坍缩成裸奇点，**当且仅当**壳的物质能量密度为负 | Boulware 1973 |
| 超极端 RN 的稳定性 | 曲率边界条件下 $\ell\ge2$ 不稳定（Dotti–Gleiser）；能量有限边界条件下没有增长模式（Delucchi）。两者都不讨论运动方程 | Dotti, Gleiser & Pullin 2007；Delucchi 2026 |
| 弯曲时空的自力 | 自力依赖整段过去历史；静止在黑洞外的电荷受排斥性自力；电磁与引力自力的耦合一阶已有完整方程，都不讨论 runaway | DeWitt & Brehme 1960；Smith & Will 1980；Zimmerman & Poisson 2014；Linz, Friedman & Wiseman 2014 |
| 自力与宇宙监督 | 带电粒子落入近极端黑洞，自力阻止它把黑洞"充成"裸奇点；轨迹事先给定，不讨论 runaway | Zimmerman, Vega, Poisson & Haas 2013；Sorce & Wald 2017 |

### 1.3 还没有人回答的问题

四轮检索（覆盖到 2026 年）都没有找到：

1. **有尺寸带电物体的自加速判据，加上引力后怎么变。** 连 $G$ 一阶（后牛顿近似）的计算都没有人做过。Steane 2015 只做了能量论证，没有动力学。
2. **ADM 带电球壳的动力学**（尤其是平移运动）。
3. **Newman 一系明确提出的问题**："引力能否稳定运动……据我们所知，答案未知"（Kozameh, Newman 等 2009）。

**为什么没人做？** 不是没人想到，而是主流的关注点在电子上，而对电子引力只有 $10^{-43}$，可以忽略（第 1.1 节）。对致密带电物体，这个问题仍然开放。

### 1.4 狄拉克选择题的现状

没有人完整写下第 3.2 节的猜想 B，但几块零件已经发表：

- **Steane 2015**：LL 方程与 Eliezer/Ford–O'Connell 方程"只在下一阶不同"，差别取决于电荷分布的细节。
- **Galley, Leibovich & Rothstein 2012**（对 Forgács 等人评论的回复）：世界线上的 $\int a^2\,d\tau$ 项可以通过平移世界线消掉，并入偶极项；它的系数正比于 $e^2R$，点极限下为零。这一项正是 Eliezer 型的 $P_2$ 项，但他们没有这样说。
- **Rowe 1975、1978**：推导 LD 方程的歧义（Schott 项还是变质量）；量纲论证：经典电磁学没有长度常数，场张量里不可能出现高阶导数项。
- **Sokolov 2009**：重新提起狄拉克的"任意选择"，但没有把它当作结构问题。
- **Kiessling 2012、2019**：在改造的电动力学（Bopp–Podolsky）中，守恒律能确定受力；粒子动量 $p=m\gamma v$ 是假设的。

### 1.5 一条贯穿的线索：选解条件

在点粒子的处理里，runaway 从来不是被方程本身消除的，而是被一个额外的条件去掉：

| 场合 | 选解条件 |
|---|---|
| 平直时空 LD 方程 | Dirac 的渐近条件（加速度在远未来趋于零），或降阶 |
| Robinson–Trautman Einstein–Maxwell | 要求过去表观视界存在 |
| 超极端 RN 裸奇点 | 奇点处的边界条件：能量有限 vs 曲率条件 |

**有尺寸的物体不需要任何选解条件**：裸质量为正就够了。这是本项目坚持用有尺寸模型的原因。

### 1.6 两个物体

| 现象 | 说明 | 来源 |
|---|---|---|
| 时间对称电动力学的稳定圆轨道 | 一半推迟、一半超前的相互作用下，两电荷可永久绕转而不辐射 | Schild 1963 |
| 两体问题是延迟微分方程 | 需要整段历史，而非初始条件 | 经典；TheoPhysEval 题 1e |
| 静态多体平衡 | 每个物体满足 $q = \sqrt G\, m$ 时，电斥力与引力精确抵消，任意摆放均可静止 | Majumdar 1947；Papapetrou 1947 |
| 极端黑洞的慢速散射 | 从静态平衡出发让两个极端带电黑洞运动 | Ferrell & Eardley 1987 |

---

## 2. 核心计算：两个临界半径

用 ADM 带电薄球壳（半径 $R$，总质量 $m$，电荷 $q$，取 $G=c=1$）。总质量满足

$$m = m_0 + \frac{q^2}{2R} - \frac{m^2}{2R},$$

其中 $m_0$ 是裸质量，$q^2/2R$ 是电磁自能，$-m^2/2R$ 是引力束缚能。

平直时空里，"尺寸大于经典半径"和"裸质量为正"是同一个条件。**加上引力后，它们分开了：**

| 判据 | 出现自加速的条件 |
|---|---|
| 裸质量判据 | $m_0 < 0 \iff R < R_{m_0} = \dfrac{q^2 - m^2}{2m}$ |
| 尺寸判据（Levine–Moniz–Sharp 原文的字面推广） | $R < R_{\rm cl} = \dfrac{q^2}{2m}$ |

两者之差恰为 $R_{\rm cl} - R_{m_0} = m/2$（即 $Gm/2c^2$），是 $G$ 的**一阶**量。

**按电荷质量比 $\xi = q/(\sqrt G\,m)$ 分三种情况：**

- **$\xi \le 1$：两种判据一致，都说没有自加速，而且答案几乎是平凡的。** 静态球壳必须在自己的视界外面：
  $$R > r_+ = m + \sqrt{m^2 - q^2} \;\ge\; m \;\ge\; \frac{q^2}{m} \;>\; R_{\rm cl} \;\ge\; R_{m_0}.$$
  任何没有塌成黑洞的球壳，都自动比经典半径大。**这样的物体无法被压到经典半径以内而不形成视界。**
- **$\xi > 1$，且 $\xi$ 不太大：两种判据在一条带里预言相反。** 在
  $$R_{m_0} < R < R_{\rm cl}$$
  这条带里，裸质量判据说没有自加速，尺寸判据说有。带宽 $m/2$ 与引力半径同量级。**这是本项目要判决的地方。**
- **$\xi \gg 1$：引力可以忽略。** 判决带相对于 $R_{\rm cl}$ 的比例是 $1/\xi^2$。对电子这样 $\xi\sim10^{21}$ 的对象，只有 $10^{-43}$（但电子本来就不在研究范围内，见第 0.2 节）。

---

## 3. 核心猜想

### 3.1 猜想 A：引力把自加速的临界半径缩小了 $Gm/2c^2$

> **猜想 A**
> 在 Einstein–Maxwell 理论中，正确的判据是**裸质量判据**：一个有尺寸的带电物体出现自加速，当且仅当它的裸质量为负。等价地说，引力把自加速的临界半径从 $q^2/2mc^2$ 缩小到 $(q^2 - Gm^2)/2mc^2$。

**支持的证据**
- 平直时空中，现代文献都用裸质量来表述这个判据（Steane 2015；D'Andrea, Luty & Verhaaren 2019；Medina 2006）。
- Boulware 1973 在完整 GR 里得到了它的**径向**版本：带电薄壳坍缩成裸奇点，当且仅当物质能量密度为负。但径向运动没有偶极辐射（Birkhoff 定理），而自加速是偶极（平移）现象，所以这只是旁证。

**判决方法**：第 2 步。先在弱场下算到 $G$ 的一阶（第 2-0 步），再到强场确认。

**如果猜想 A 不成立**（尺寸判据才对）：引力不改变自加速的判据；$\xi\le1$ 时没有自加速，只是因为视界挡住了。这也是一个干净的结果。

### 3.2 猜想 B：狄拉克选择题的答案

> **猜想 B**
> Eliezer 型的推广（在 $P_{\rm mech}$ 中加入 $P_2\ddot u$、$P_4 u^{(4)}$ 等项），降阶后在 $\tau_0$ 的一阶全部等价于 LL 方程；从二阶开始出现的差别，对应物体的内部结构（如电极化率），而不是同一个点电荷的不同运动定律。任何有尺寸的电荷取点极限，都会选出 $P_{2i}=0$，也就是 LD 方程。

一维非相对论下的初步计算支持猜想 B：

$$m\dot v + P_2\dddot v = F + m\tau_0\ddot v
\;\Rightarrow\;
m\dot v = F + \tau_0\dot F + \left(\tau_0^2 - \frac{P_2}{m}\right)\ddot F + \cdots$$

$P_2$ 从二阶开始进入，形式上等价于一个 $\propto E^2$ 的极化率项。

**我们要做的**：第 1.4 节的零件已经存在；我们的工作是把它们拼成一个完整的陈述（Eliezer 的 $P_{2i}$ 歧义正是世界线重定义的冗余），推广到相对论情形和所有阶，并用有尺寸模型证明点极限选出 $P_{2i}=0$。

**定位**：独立短文，新意中等。不作为主结论。

---

## 4. 路线图

### 4.1 总览

| 步骤 | 内容 | 状态 |
|---|---|---|
| 第 0 步 | 文献核查 | ✅ 已完成四轮 |
| 第 1 步 | 平直时空：复现 Levine–Moniz–Sharp；狄拉克选择题（猜想 B） | ⏭ 下一步 |
| 第 2 步 | 加入引力：ADM 带电球壳，判决两种判据（猜想 A） | 主线，最关键 |
| 第 3 步 | 独立验证（从第 2 步起并行） | |
| 第 4 步 | 两个带电物体加引力（第二阶段） | |
| 附：数学支线 | 超极端 RN 点电荷的 $\ell=1$ 扇区 | 可选，不影响主线 |

---

### 第 0 步：文献核查（已完成四轮）

过程见调研记录。结论已写进第 1 节。

**尚未完成**（不阻塞第 1 步）：
- Boulware 1973、Levine–Moniz–Sharp 1977 的全文（在付费墙后，需要机构访问）。
- 第 0.2 节"量子引力方向"的几篇文献（DeWitt 1964 等）的核对。

---

### 第 1 步：平直时空

这是地基。第 2 步的所有结果，在 $G \to 0$ 时都必须回到这里。它也是检验 AI 推导可靠性的标尺。

**1a. 练手：检验推导流程**

- 推出 Lun–Chow 的特征值 $\lambda_{1,+} = 3M/2Q^2$（Robinson–Trautman Einstein–Maxwell 在 RN 附近线性化）。
- 答案已知，用来跑通"两次独立推导 + SymPy 校验"的流程。

**1b. 两种方法独立重推 LD 方程**

- Dirac 1938 的管面积分法；推迟坐标法（Poisson 的处理）。
- 顺带复现 Rowe 1975 的歧义：换成光锥端面，得到变质量方程。

**1c. 复现 Levine–Moniz–Sharp（关键路径）**

- 带电球壳的自力，写出微分–差分方程。
- 数值上确认：尺寸大于经典半径时没有自加速，小于时有。
- 写明平直时空中"尺寸判据"与"裸质量判据"的等价性，以及等价性在哪一步用到了 $E_{\rm self} = q^2/2R$。**这一步决定第 2 步怎么区分两种判据。**
- 写成方便加入引力的形式。

**1d. 狄拉克的选择题（猜想 B）**

- 写出 Eliezer 型 $P_{\rm mech}$ 的一般形式。
- 相对论情形下逐阶降阶到 $\tau_0^2$，与"LD + 极化率世界线作用量"逐项比较。
- 判断哪些 $P_{2i}$ 能通过重新定义世界线消掉，与 Galley–Leibovich–Rothstein 2012 的 $\int a^2$ 结果对照。
- 用 1c 的球壳模型检查：点极限下 $P_{2i}$ 是否自动为零。

**产出**：一套经过两种方法、两次独立 AI 运行和 SymPy 校验的推导；可复用的数值代码；猜想 B 的短文初稿。

**决策**：无法复现 Levine–Moniz–Sharp → 先查清原因，不进入第 2 步。

---

### 第 2 步：加入引力，判决两种判据

**模型**：ADM 带电薄球壳，放入完整的 Einstein–Maxwell 理论。壳内平直，壳外 Reissner–Nordström，壳面满足 Israel 结合条件。

**为什么是这个模型**：它有尺寸，不需要任何选解条件（第 1.5 节）；它直接回答 Newman 一系留下的公开问题；它能检验 Levine–Moniz–Sharp 判据在 GR 中的形式。

**2-0. 弱场判决（先做，几周）**

- 两个临界半径之差 $Gm/2c^2$ 是 $G$ 的一阶量（第 2 节），所以**不需要强场就能判决**。
- 把 1c 的球壳模型算到 $G$ 的一阶（后牛顿近似）：引力束缚能、电磁场能量的引力效应、壳的应力。
- 问：自加速的临界半径移动了 $Gm/2c^2$（裸质量判据），还是没有移动（尺寸判据）？
- 在 $\xi \gg 1$ 时，壳的半径 $\sim q^2/m \gg Gm$，引力确实是弱场，这个近似是可控的。
- 四轮检索都没有找到这个计算。

**2a. 复现 ADM 静态结果**

- 解出 $m(R, m_0, q)$，验证 $R\to0$ 时 $m\to |q|/\sqrt G$。
- 在 $(\xi, R)$ 平面上画出四个区域：视界内（不可达）、两种判据都说有自加速、判决带、两种判据都说没有。

**2b. 强场动力学**

- **径向（热身）**：复现 Boulware 1973；检查有没有指数增长的模式。径向运动判决不了自加速。
- **平移**：壳内外 $\ell=1$ 扰动用结合条件拼接，或用数值相对论。取 $\xi$ 略大于 1、$R$ 落在判决带内，检查自加速的特征根。
- 与 2-0 的弱场结果和第 1 步 $G\to0$ 的结果对接。

**决策（整个项目最关键的一个）**

- 裸质量判据成立 → 猜想 A 成立：引力把自加速的临界半径缩小了 $Gm/2c^2$。
- 尺寸判据成立 → 引力不改变判据，也是干净的结果。
- 都不对，或结果依赖模型细节（球壳 vs 实心球 vs 尘埃）→ 问题不是普适的，需要重新表述。

**技术风险**：平移需要 $\ell=1$ 拼接或数值相对论；$\xi\to1$ 时几何可能退化。2-0 的意义正在于先绕开这两点拿到结论。

---

### 第 3 步：独立验证（从第 2 步起并行）

**原则**：没有领域内专家独立核对，任何结果都不对外讲。

- 第 2 步一出现苗头，就把推导和代码整理成可独立运行的形式。
- 候选专家：Newman、Kozameh 一系（提出公开问题的人）；Steane（最接近的平直时空工作）；Gralla、Harte、Wald；Kiessling、Tahvildar-Zadeh；Capra 自力社区；TheoPhysEval 修正电动力学题目的出题人（通讯作者 Suvrat Raju）。

**产出**：至少一位专家书面确认关键结果。

---

### 第 4 步：两个带电物体（第二阶段）

在第 2 步有结论之后再展开。

- **平直时空基准**：带推迟相互作用的两体问题（延迟微分方程）；对照 Schild 的时间对称圆轨道。
- **加入引力**：从 Majumdar–Papapetrou 静态平衡（$q=\sqrt G m$）出发，让两个物体动起来；参考 Ferrell–Eardley 的极端黑洞慢速散射。问题：偏离极端条件时，辐射损失与引力束缚如何竞争？

---

### 附：数学支线，超极端 RN 点电荷的 $\ell=1$ 扇区

**性质**：经典点电荷不对应任何真实粒子（第 0.2 节），所以这只是"经典理论自身"的数学问题，**不是主线**。

**问题**：Delucchi 2026 把超极端 RN 当作点电荷的自场，但把 $\ell=1$ 当作纯规范剔除了。这在 Einstein–Maxwell 中不成立：电磁偶极自由度是真实的。让点电荷能运动，在 $\ell=1$ 扇区推出世界线方程，看自加速是否出现，以及是否取决于奇点处自伴扩张的选择。

**工作量**：轴向 $\ell=1$ 约一天（势 $V = f(2/r^2 + 4Q^2/r^4)$，已手算，未与文献对照）；极性 $\ell=1$ 几天；加入世界线几周。

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

- 摘要不等于结论。第一轮基于摘要的两条判断，读全文后都被推翻。
- 每条引用标注核实程度：【全文】【摘要】【二手】。
- 转述别人的结论时，与原文措辞对照。
- 凭记忆给出的文献，一律标【待核】。

**记录**

- 所有推导、代码保存在仓库中，按步骤分目录。
- 每一步结束写一页小结：做了什么、结论、是否通过决策点。

---

## 6. 风险清单

| 风险 | 可能性 | 影响 | 应对 |
|---|---|---|---|
| 核心问题已有人做过 | 低（四轮未找到） | 高 | 持续跟踪新文献 |
| 尺寸判据才是对的（猜想 A 不成立） | 中 | 中 | 也是干净的结果，可以写 |
| 结论被认为只对致密天体有意义 | 高 | 中 | 第 0.4 节已主动说明；卖点放在"经典理论自身"的判据 |
| 后牛顿近似下物体的应力、刚性处理不当 | 中 | 高 | 用薄壳（结构最简单）；与 ADM 精确静态解对照 |
| 强场平移需要数值相对论 | 中 | 中 | 先做 2-0 弱场判决 |
| 猜想 B 是 EFT 圈内常识 | 中 | 低 | 作为独立短文 |
| AI 推导或文献检索出错 | 中 | 高 | 第 5 节的规则；第 1a 步先用已知答案练手 |
| 找不到专家愿意核对 | 中 | 高 | 从第 2 步起提前联系 |

---

## 7. 对外发布的条件

只有同时满足以下三条才对外发布：

1. 第 2 步得到明确结果。
2. 至少一位领域专家书面确认。
3. 结论能用一句外行听得懂的话讲清楚，并且**明确说清楚它针对的是什么物体**（第 0.1 节）。

---

## 8. 参考文献

> 标注：**【已核】** 出版信息经 Crossref、INSPIRE 或 arXiv 核实；**【待核】** 尚未核实。内容的核实程度见调研记录。

### 8.1 平直时空：运动方程及其推导

- P. A. M. Dirac, *Classical theory of radiating electrons*, *Proc. R. Soc. A* 167, 148–169 (1938).【已核】
- C. J. Eliezer, *The interaction of electrons and an electromagnetic field*, *Rev. Mod. Phys.* 19, 147–184 (1947).【已核】
- E. J. Moniz, D. H. Sharp, *Absence of runaways and divergent self-mass in nonrelativistic quantum electrodynamics*, *Phys. Rev. D* 10, 1133–1136 (1974).【已核】
- E. G. P. Rowe, *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975).【已核，本地有】
- R. Tabensky, *Electrodynamics and the electron equation of motion*, *Phys. Rev. D* 13, 267–273 (1976).【已核】
- H. Levine, E. J. Moniz, D. H. Sharp, *Motion of extended charges in classical electrodynamics*, *Am. J. Phys.* 45, 75–78 (1977).【已核】
- E. J. Moniz, D. H. Sharp, *Radiation reaction in nonrelativistic quantum electrodynamics*, *Phys. Rev. D* 15, 2850–2865 (1977).【已核】
- E. G. P. Rowe, *Structure of the energy tensor in the classical electrodynamics of point particles*, *Phys. Rev. D* 18, 3639 (1978).【已核，本地有】
- G. W. Ford, R. F. O'Connell, *Phys. Lett. A* 157, 217 (1991).【已核】
- H. Spohn, *The critical manifold of the Lorentz–Dirac equation*, *Europhys. Lett.* 50, 287–292 (2000), arXiv:physics/9911027.【已核】
- R. Medina, *Radiation reaction of a classical quasi-rigid extended particle*, *J. Phys. A* 39, 3801 (2006), arXiv:physics/0508031.【已核】
- D. V. Gal'tsov, P. Spirin, *Radiation reaction reexamined: bound momentum and Schott term*, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121.【已核】
- A. D. Yaghjian, *Relativistic Dynamics of a Charged Sphere* (Springer; 1992, 2006, 2022).【已核】
- I. V. Sokolov, *JETP* 109, 207 (2009), arXiv:0906.1150.【已核】
- S. E. Gralla, A. I. Harte, R. M. Wald, *Rigorous derivation of electromagnetic self-force*, *Phys. Rev. D* 80, 024031 (2009), arXiv:0905.2391.【已核】
- C. R. Galley, A. K. Leibovich, I. Z. Rothstein, *Finite size corrections to the radiation reaction force in classical electrodynamics*, *Phys. Rev. Lett.* 105, 094802 (2010), arXiv:1005.2617.【已核】
- P. Forgács, T. Herpay, P. Kovács, *Phys. Rev. Lett.* 109, 029501 (2012), arXiv:1202.6289；C. R. Galley, A. K. Leibovich, I. Z. Rothstein, Reply, *Phys. Rev. Lett.* 109, 029502 (2012), arXiv:1206.4773.【已核】
- M. K.-H. Kiessling, arXiv:1107.2286 (2012)，收入 *Quantum Field Theory and Gravity* (Springer), pp. 299–335.【已核】
- A. M. Steane, *Am. J. Phys.* 83, 256 (2015), arXiv:1402.1106.【已核】
- C. Burton, A. Noble, *Aspects of electromagnetic radiation reaction in strong fields*, *Contemp. Phys.* 55, 110 (2014), arXiv:1409.7707.【已核】
- M. K.-H. Kiessling, *Force on a point charge source of the classical electromagnetic field*, *Phys. Rev. D* 100, 065012 (2019), arXiv:1907.11239.【已核】
- B. P. Kosyakov, *Phys. Rep.* 812, 1 (2019), arXiv:1812.03290.【已核】
- M. D'Andrea, M. A. Luty, C. B. Verhaaren, arXiv:1909.00960 (2019).【已核】
- O. Birnholtz, arXiv:1410.5871.【待核】
- W. D. Goldberger, I. Z. Rothstein, *Phys. Rev. D* 73, 104029 (2006), arXiv:hep-th/0409156.【已核】
- A. Di Piazza, C. Müller, K. Z. Hatsagortsyan, C. H. Keitel, *Rev. Mod. Phys.* 84, 1177 (2012).【已核】
- A. Gonoskov, T. Blackburn, M. Marklund, S. S. Bulanov, *Rev. Mod. Phys.* 94, 045001 (2022).【已核】
- A. Ramesh Chandra *et al.*, *TheoPhysEval*, arXiv:2610.16733 (2026).【已核】

### 8.2 引力：自能、自力与渐近方法

- B. S. DeWitt, R. W. Brehme, *Radiation damping in a gravitational field*, *Ann. Phys.* 9, 220–259 (1960).【已核】
- R. Arnowitt, S. Deser, C. W. Misner, *Gravitational-electromagnetic coupling and the classical self-energy problem*, *Phys. Rev.* 120, 313–320 (1960).【已核】
- J. M. Hobbs, *A vierbein formalism of radiation damping*, *Ann. Phys.* 47, 141–165 (1968).【已核】
- A. G. Smith, C. M. Will, *Force on a static charge outside a Schwarzschild black hole*, *Phys. Rev. D* 22, 1276–1284 (1980).【已核】
- J. F. Donoghue, *Phys. Rev. D* 50, 3874–3888 (1994), arXiv:gr-qc/9405057.【已核】
- C. N. Kozameh, E. T. Newman, J. G. Ortega, G. Silva-Ortigoza, *Gen. Rel. Grav.* 41, 1139 (2009), arXiv:0803.3921.【已核】
- T. M. Adamo, C. N. Kozameh, E. T. Newman, *Living Rev. Relativ.* 12, 6 (2009)；15, 1 (2012), arXiv:0906.2155.【已核】
- E. T. Newman, *Class. Quantum Grav.* 28, 245003 (2011), arXiv:1109.4106.【已核】
- E. Poisson, A. Pound, I. Vega, *Living Rev. Relativ.* 14, 7 (2011), arXiv:1102.0529.【已核】
- C. N. Kozameh, J. I. Quiroga, arXiv:1202.5972 (2012)；J. I. Quiroga, *Gen. Rel. Grav.* 49, 20 (2017), arXiv:1711.11384.【已核】
- P. Zimmerman, E. Poisson, *Phys. Rev. D* 90, 084030 (2014), arXiv:1406.5111.【已核】
- T. M. Linz, J. L. Friedman, A. G. Wiseman, *Phys. Rev. D* 90, 024064 (2014), arXiv:1404.7039；*Phys. Rev. D* 90, 084031 (2014), arXiv:1406.5112.【已核】
- E. T. Newman, *Gen. Rel. Grav.* 50, 111 (2018).【已核】
- M. K.-H. Kiessling, A. S. Tahvildar-Zadeh, *Int. J. Mod. Phys. D* 28, 1930017 (2019), arXiv:1906.00308.【已核】

### 8.3 Robinson–Trautman Einstein–Maxwell

- E. T. Newman, R. Posadas, *Phys. Rev.* 187, 1784–1791 (1969)；*J. Math. Phys.* 12, 2319–2327 (1971).【已核】
- A. W. C. Lun, E. W. M. Chow, arXiv:gr-qc/9409024 (1994).【arXiv 已核，文集信息为二手】
- C. N. Kozameh, H.-O. Kreiss, O. Reula, *Class. Quantum Grav.* 25, 025004 (2008), arXiv:0708.1933.【已核】

### 8.4 超极端 Reissner–Nordström、点电荷时空与宇宙监督

- G. Dotti, R. J. Gleiser, J. Pullin, *Phys. Lett. B* 644, 289 (2007), arXiv:gr-qc/0607052.【已核】
- G. Dotti, R. J. Gleiser, *Class. Quantum Grav.* 27, 185007 (2010), arXiv:1001.0152.【已核】
- Delucchi, *Class. Quantum Grav.* 43, 195008 (2026), arXiv:2602.19015.【已核】
- A. S. Tahvildar-Zadeh, *Rev. Math. Phys.* 23, 309–346 (2011), arXiv:1012.1400.【已核】
- A. Burtscher, M. K.-H. Kiessling, A. S. Tahvildar-Zadeh, *Class. Quantum Grav.* 38, 185001 (2021), arXiv:1901.00813.【已核】
- H. Ponce de León, *Gen. Rel. Grav.* 36, 1451 (2004), arXiv:gr-qc/0310117.【已核】
- P. Zimmerman, I. Vega, E. Poisson, R. Haas, *Self-force as a cosmic censor*, *Phys. Rev. D* 87, 041501 (2013), arXiv:1211.3889.【已核】
- J. Sorce, R. M. Wald, *Phys. Rev. D* 96, 104014 (2017), arXiv:1707.05862.【已核】
- 数学支线的工具文献：F. J. Zerilli, *Phys. Rev. D* 9, 860 (1974)；L. M. Burko, arXiv:gr-qc/9507012、gr-qc/9902002.【待核】

### 8.5 带电壳层与带电星体

- W. Israel, *Nuovo Cim. B* 44, 1–14 (1966)；erratum 48, 463 (1967).【已核】
- V. de la Cruz, W. Israel, *Gravitational bounce*, *Nuovo Cim. A* 51, 744–760 (1967).【已核】
- K. Kuchař, *Czech. J. Phys. B* 18, 435–463 (1968).【已核】
- D. G. Boulware, *Naked singularities, thin shells, and the Reissner–Nordström metric*, *Phys. Rev. D* 8, 2363–2368 (1973).【已核】
- D. Anninos, T. Rothman, *Phys. Rev. D* 65, 024003 (2002), arXiv:gr-qc/0108082.【已核】
- G. Dillon, *Eur. Phys. J. Plus* 128, 72 (2013), arXiv:1303.2706.【已核】
- R. Casadio, *Int. J. Mod. Phys. A* 28, 1350088 (2013), arXiv:1303.1274.【已核】

### 8.6 量子引力方向（第 0.2 节，不在本项目范围内）

- N. Arkani-Hamed, L. Motl, A. Nicolis, C. Vafa, *JHEP* 0706, 060 (2007), arXiv:hep-th/0601001.【已核】
- D. Harlow, B. Heidenreich, M. Reece, T. Rudelius, *Rev. Mod. Phys.* 95, 035003 (2023), arXiv:2201.08380.【已核】
- B. S. DeWitt, *Gravity: a universal regulator?*, *Phys. Rev. Lett.* 13, 114 (1964).【待核】
- C. J. Isham, A. Salam, J. Strathdee, *Phys. Rev. D* 3, 1805 (1971).【待核】
- S. P. Robinson, F. Wilczek, *Phys. Rev. Lett.* 96, 231601 (2006).【待核】

### 8.7 两个物体

- A. Schild, *Electromagnetic two-body problem*, *Phys. Rev.* 131, 2762–2766 (1963).【已核】
- S. D. Majumdar, *Phys. Rev.* 72, 390 (1947)；A. Papapetrou, *Proc. R. Irish Acad. A* 51, 191 (1947).【待核】
- R. C. Ferrell, D. M. Eardley, *Phys. Rev. Lett.* 59, 1617–1620 (1987).【已核】

---

## 9. 版本记录

旧版本保存在 git 历史中，可用 `git show <commit>:docs/roadmap.md` 查看。

| 版本 | commit | 主要内容 |
|---|---|---|
| v1 | `3f3549e` | 初版路线图，目标是"引力消除电子自加速" |
| v1 + 第一、二轮 | `722f1a7`、`081b942` | 附加文献核查；视界版本撤回 |
| v2 | `db4ed75` | 整体重写；文献核查移到独立文件 |
| v2.1 | `3187f4a` | 判决带改到 $\xi>1$；加入量子一侧的判据 |
| **v3** | 本次 | **研究对象明确为经典宏观带电物体，不研究基本粒子**（第 0.1、0.2 节，含尺度表）；加入第四轮检索：学术界已解决什么、还剩什么（第 1.1、1.3 节）；第 2 步加入弱场判决（2-0）；超极端 RN 点电荷降为数学支线；删去"宇宙监督"备注、撤回猜想的历史和原"第 3 步：适用范围"（已并入第 0、2 节） |
