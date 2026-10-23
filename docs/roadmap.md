# 经典带电质点的运动：研究路线图

> 版本：2026-10-23 · 草稿
> 状态：规划阶段，所有"结果"均为待检验的猜想

---

## 0. 这件事要做什么

**出发点**：只用经典理论，即麦克斯韦电动力学和爱因斯坦广义相对论，把"一个质量为 m、电荷为 q 的质点会怎样运动"这个问题推到底，再推到两个质点。看理论自己会不会给出意想不到的结论。

**方法论**：不去修补方程，而是把方程当真。爱因斯坦当年就是这样做的：他没有修改麦克斯韦方程，而是认真对待它和力学之间的矛盾，矛盾自己推出了相对论。

**想回答的两个核心问题**

1. **狄拉克的选择题**：狄拉克 1938 年推导电子运动方程时承认，能动量守恒只约束了一部分，粒子的"机械动量" $P_{\rm mech}$（也就是 $B^\mu$）该怎么写，是他挑了最简单的一种。其他写法是不同的物理，还是同一个物理的不同写法？
2. **引力能否消除"自加速"**：狄拉克方程（Lorentz–Dirac 方程，下称 LD）有一个荒谬的解：电子不受力也会自己加速到接近光速。把广义相对论加进来，这个解还在不在？

**目标结果（如果成功）**，一句话版本：

> 困扰物理学家近百年的"电子自加速"悖论，在加入爱因斯坦引力之后消失了。

这句话是目标，不是结论。它可能成立，可能只对一部分对象成立，也可能不成立。路线图的设计就是要尽早知道是哪一种。

---

## 1. 已知的事实：经典理论已经给过哪些"惊喜"

做之前先把地图画清楚，避免重复别人的工作。

### 1.1 一个带电质点，只有电动力学

| 现象 | 说明 | 来源 |
|---|---|---|
| 自能发散 | 点电荷的场能量无穷大 | 经典 |
| 4/3 问题 | 由场的能量和动量算出的电磁质量之比为 4/3 而非 1 | Abraham, Lorentz；Poincaré 应力 |
| runaway 与预加速 | LD 方程的解可以不受力自加速，或在力到来之前开始运动 | Dirac 1938 |
| 经典电子半径 | $r_e = e^2/mc^2$，比康普顿波长小 $1/\alpha = 137$ 倍：经典理论在病灶出现前已被量子力学接管 | 经典 |
| 有尺寸电荷的 runaway 条件 | 带电球半径大于经典电子半径（裸质量为正）时没有 runaway；半径更小（裸质量为负）时出现 | Levine, Moniz & Sharp, *Am. J. Phys.* 45, 75 (1977)；Moniz & Sharp, *Phys. Rev. D* 15, 2850 (1977)【待核】 |
| 降阶 | LD 的物理解落在一个"临界流形"上，Landau–Lifshitz（LL）方程是它的一阶近似 | Spohn 2000；Gralla–Harte–Wald 2009 |

### 1.2 一个带电质点，加上广义相对论

| 现象 | 说明 | 来源 |
|---|---|---|
| **引力治好了自能发散** | 带电球壳在 Einstein–Maxwell 理论中，半径 → 0 时总质量趋于 $m \to |q|/\sqrt{G}$，与裸质量无关 | Arnowitt, Deser & Misner, *Phys. Rev.* 120, 313 (1960) |
| 经典理论算出普朗克尺度 | 对电子，$q/\sqrt G \approx 1.9\times10^{-6}$ g $=\sqrt\alpha\, m_{\rm Pl}$；对应长度 $\sqrt G\, q/c^2 \approx 1.4\times10^{-34}$ cm $=\sqrt\alpha\,\ell_{\rm Pl}$。全程没有 ħ | 由 ADM 结果推出 |
| 电子是裸奇点 | 电子的 $q/(\sqrt G\, m) \approx 2\times10^{21}$，Reissner–Nordström 解在此情形下没有视界 | 经典 GR |
| 弯曲时空自力的"尾巴" | 自力依赖整段过去历史，因为弯曲时空中波不只沿光锥传播 | DeWitt & Brehme 1960；Hobbs 1968；DeWitt & DeWitt 1964 |
| 黑洞推开静止电荷 | 静止在 Schwarzschild 黑洞外的电荷受到排斥性自力 | Smith & Will, *PRD* 22, 1276 (1980) |
| 自由下落电荷是否辐射 | 等效原理与"加速电荷辐射"之间的悖论，答案依赖观测者 | Rohrlich；Boulware 1980 |
| 电磁与引力自力耦合 | 粒子的场扰动度规，度规扰动又反过来作用于场，一阶已有完整方程 | Zimmerman & Poisson, *PRD* 90, 084030 (2014), arXiv:1406.5111；Linz, Friedman & Wiseman, *PRD* 90, 084031 (2014), arXiv:1406.5112 |

### 1.3 两个质点

| 现象 | 说明 | 来源 |
|---|---|---|
| 经典氢原子塌缩 | 电子绕核辐射，约 $10^{-11}$ s 落入 | 经典 |
| 时间对称电动力学的稳定圆轨道 | 一半推迟、一半超前的相互作用下，两电荷可永久绕转而不辐射 | Schild, *Phys. Rev.* 131, 2762 (1963) |
| 两体问题是延迟微分方程 | 需要整段历史，而非初始条件 | 经典；TheoPhysEval 题 1e |
| 静态多体平衡 | 每个粒子满足 $q = \sqrt G\, m$ 时，电斥力与引力精确抵消，任意摆放均可静止 | Majumdar 1947；Papapetrou 1947 |
| 极端黑洞的慢速散射 | 从静态平衡出发让两个极端带电黑洞运动 | Ferrell & Eardley, *PRL* 59, 1617 (1987) |

---

## 2. 一个必须先面对的数量级问题

**在电子的参数下，引力弱到几乎不存在。**

$$\frac{G m_e^2}{e^2} \approx 2.4\times10^{-43}$$

如果把引力当成小修正逐阶加入，它不可能改变 runaway 这种量级的行为。

ADM 的惊喜之所以出现，是因为它是**非微扰**的：电荷被压缩到 $\sqrt G\, q/c^2 \approx 10^{-34}$ cm 时，引力自能与电磁自能同量级。但在那个极限下，总质量变成 $q/\sqrt G$，比电子重 $10^{21}$ 倍。

**风险**：就算"引力消除自加速"成立，它也可能只适用于质量与电荷相当的对象（接近极端的带电黑洞），而不适用于电子。这会直接改变最终结论的表述，所以第 3 步专门用来判断这一点。

**静态估算已经指向否定的答案。** 在 ADM 球壳模型里把总质量固定为电子质量 $m=m_e$，反解裸质量：

$$m_0 = m_e - \frac{q^2}{2R} + \frac{G m_e^2}{2R}$$

引力项与电磁项之比为 $Gm_e^2/q^2\sim10^{-43}$，所以 $m_0<0$ 的条件仍是 $R\lesssim r_e/2$，引力对这个临界半径的相对修正只有 $10^{-43}$。ADM 的"点极限下总质量有限"靠的是总质量本身趋于 $|q|/\sqrt G$（即 $\xi\to1$）；只要质量钉在 $m_e$，这条路就不通。电子要么半径 $\sim r_e$，要么裸质量为负，正好回到 Levine–Moniz–Sharp 的 runaway 条件。

因此，至少在球壳模型中，猜想 A 对电子不成立；值得做的问题是 $\xi\sim1$ 附近 runaway 特征根的行为。

---

## 3. 核心猜想

> **猜想 A（引力消除 runaway）**
> 平直时空中，点电荷的 runaway 根源于点极限要求裸质量为负（Levine–Moniz–Sharp 的结果指向这一点）。在 Einstein–Maxwell 理论中，ADM 表明点极限下总质量有限，裸质量不必为负。因此引力可能从根上消除 runaway。
>
> *第 0 步第一轮检索后需要改写，修订版猜想 A′ 见第 9.4 节；第二轮检索后再次修订，现行版本为第 10.4 节的猜想 A″。*

> **猜想 B（狄拉克选择题的答案）**
> Eliezer 型的 $B^\mu$ 推广（在 $P_{\rm mech}$ 中加入 $P_2\ddot u$、$P_4 u^{(4)}$ 等项），降阶后在 $\tau_0$ 一阶全部等价于 LL；从二阶开始出现的差别，对应粒子的内部结构（如电极化率），而不是同一个点电荷的不同运动定律。任何有尺寸的电荷取点极限，都会选出 $P_{2i}=0$，也就是 LD。

一维非相对论下的初步计算支持猜想 B：

$$m\dot v + P_2\dddot v = F + m\tau_0\ddot v
\;\Rightarrow\;
m\dot v = F + \tau_0\dot F + \left(\tau_0^2 - \frac{P_2}{m}\right)\ddot F + \cdots$$

$P_2$ 从二阶开始进入，形式上等价于一个 $\propto E^2$ 的极化率项。相对论情形尚未验证。

---

## 4. 路线图

### 总览

| 步骤 | 内容 | 决策点 |
|---|---|---|
| 第 0 步 | 文献核查 | 有没有人做过 |
| 第 1 步 | 平直时空：把狄拉克方程彻底吃透 | 能否复现 Levine–Moniz–Sharp；猜想 B 是否成立 |
| 第 2 步 | 加入引力：带电球壳模型 | **runaway 是否消失（最关键）** |
| 第 3 步 | 确定适用范围（静态部分可提前到第 1 步之前，见第 2 节） | 对电子是否成立 |
| 第 4 步 | 独立验证（从第 2 步起并行） | 有没有专家签字 |
| 第 5 步 | 两个质点（第二阶段） | — |

---

### 第 0 步：文献核查

**目标**：确认核心问题"在完全的 Einstein–Maxwell 理论中，带电质点是否还有 runaway"有没有人系统做过。

**要查的线索**

- ADM 1960 及其后续：有没有人从静态推广到动力学？
- Levine–Moniz–Sharp 1977：有尺寸电荷 runaway 条件的原始论证，以及后来的推广（Yaghjian 的 *Relativistic Dynamics of a Charged Sphere*）。
- 带电尘埃球和带电壳层的引力坍缩（Bonnor 等）。
- 弯曲时空自力综述：Poisson, Pound & Vega, *Living Reviews in Relativity* 14, 7 (2011)。
- 电磁与引力耦合自力（arXiv:1406.5111、1406.5112）：它们是否讨论过 runaway？（第二轮：都没有，见第 10 节）
- 修正场论下的良定性：Kiessling 等关于 Bopp–Landé–Thomas–Podolsky 理论的工作。
- 猜想 B：世界线有效场论（Goldberger–Rothstein 一路）中，高阶导数项与极化率的对应是否已有人写下。

**产出**

- 一份文献表，每篇标注：做了什么、没做什么、与本项目的关系。
- 一张"空白清单"：哪些问题确实没人做。

**决策**

- 已有人完整做过 → 换问题，或改为在其基础上推进。
- 只做了一部分 → 锁定空着的那一块。
- 完全空白 → 进入第 1 步。

**进展**：第一轮检索结果见第 9 节，第二轮（逐篇读全文）见第 10 节。第二轮结论：核心问题未见有人回答，且被 Newman 一系明确列为公开问题；进入第 1 步。

---

### 第 1 步：平直时空，把狄拉克方程彻底吃透

这是地基。第 2 步的所有结果，在 $G \to 0$ 时都必须回到这里。它也是检验 AI 推导可靠性的标尺。

**1a. 两种方法独立重推 LD**

- Dirac 1938 的管面积分法。
- 推迟坐标法（Poisson 的处理）。
- 要求：两种方法给出同一个方程；动量通量中的发散项正比于 $a^\mu$，能被质量重整化吸收。

**1b. 精确刻画 runaway**

- 一维情形求出特征根，确认 runaway 来自 $e^{\tau/\tau_0}$ 那一支。
- 写出非 runaway 解的积分形式 $m\dot v(t)=\int_0^\infty e^{-s}F(t+\tau_0 s)\,ds$，明确预加速的来源。
- 降阶级数 $\sum \tau_0^n F^{(n)}$：在什么条件下收敛，什么条件下只是渐近级数。

**1c. 有尺寸电荷模型**

- 带电球壳（或均匀带电球）的自力，复现 Levine–Moniz–Sharp 的结论：裸质量为正时无 runaway，为负时有。
- 找出临界半径，确认它就是经典电子半径的量级。
- **这一步的模型在第 2 步会直接搬进广义相对论**，所以要写成方便推广的形式。

**1d. 狄拉克的选择题（猜想 B）**

- 写出 Eliezer 型 $P_{\rm mech}$ 的一般形式（2006 年论文式 (1.3.3)）。
- 相对论情形下逐阶降阶到 $\tau_0^2$，与"LD + 极化率世界线作用量"逐项比较。
- 判断哪些 $P_{2i}$ 可以通过重新定义世界线消掉，哪些不能。
- 用 1c 的有尺寸模型检查：点极限下 $P_{2i}$ 是否自动为零。

**产出**

- 一套经过两种方法、两次独立 AI 运行和 SymPy 校验的推导。
- 可复用的代码：LD、LL、降阶级数、有尺寸电荷模型的数值求解器。
- 猜想 B 的结论（无论成立与否），本身可以写成一篇短文。

**决策**

- 无法复现 Levine–Moniz–Sharp → 先查清原因，不进入第 2 步。
- 猜想 B 成立 → 记录，作为独立结果。

---

### 第 2 步：加入引力，带电球壳模型

**模型**：ADM 的带电薄球壳，放入完全的 Einstein–Maxwell 理论。壳内平直，壳外 Reissner–Nordström，壳面满足 Israel 结合条件。

**2a. 复现 ADM 静态结果**

- 总质量满足 $m = m_0 + \dfrac{q^2}{2R} - \dfrac{G m^2}{2R}$（$c=1$），解出
  $m = \dfrac{R}{G}\left[-1+\sqrt{1+\dfrac{2G}{R}\left(m_0+\dfrac{q^2}{2R}\right)}\right]$。
- 验证 $R\to0$ 时 $m\to |q|/\sqrt G$，与 $m_0$ 无关。
- 注意此时 $q/(\sqrt G m)\to 1$：点极限自动落在极端条件上。这一点要特别留意，它可能是整个问题的关键。

**2b. 让球壳动起来**

- 先做最简单的动力学：径向振荡（球对称，无辐射，Birkhoff 定理保证壳外仍是 RN）。检查是否有指数增长的模式。
- 再做平移：施加一个小的外力，计算球壳中心的运动。这一步破坏球对称，需要微扰处理或数值相对论。
- 两种情况都要和第 1 步 $G\to0$ 的结果对接。

**2c. 核心问题**

随着球壳半径趋于零，运动方程里还有没有 runaway 那个特征根？

**决策（整个项目最关键的一个）**

- 按第 2 节的静态估算，2b 应直接瞄准 $\xi\sim1$ 附近，那里几何会退化，也是唯一可能出现新行为的区域。
- runaway 消失 → 猜想 A 成立，进入第 3 步。
- runaway 仍在 → 也是一个结果："连引力也救不了经典电子"。可以写成文章，但不是新闻；转向第 5 步或回到猜想 B。
- 结果依赖模型细节（球壳 vs 实心球 vs 尘埃） → 说明问题本身不是普适的，需要重新表述。

**技术风险**

- 非球对称的动力学可能需要数值相对论，工作量会明显增加。
- 点极限与极端条件重合，接近 $q=\sqrt G m$ 时几何可能退化，需要特别小心。

---

### 第 3 步：确定适用范围

**目标**：回答第 2 节的数量级问题。

**静态部分可先在纸上完成**（不依赖 2b）：第 2 节的反解已表明，在球壳模型里固定 $m=m_e$ 时，引力对 runaway 临界半径的修正只有 $10^{-43}$ 量级。第 1 步之前先把这一估算写成一页小结，据此决定第 2 步的投入，并预先把结论主语定为"近极端带电天体"。动力学部分（2b 之后）只需确认没有静态估算看不到的效应。

- 用比值 $\xi = q/(\sqrt G\, m)$ 给出 runaway 消失的条件。
- 区分三种情况：$\xi \le 1$（亚极端，可以形成黑洞）、$\xi \approx 1$（极端）、$\xi \gg 1$（超极端，电子 $\xi\sim10^{21}$）。
- 回答：电子落在哪个区间？结论对它是否成立？

**这一步决定最终结论的主语**

- 对电子成立 → "电子自加速悖论被引力消除"。
- 只对 $\xi\lesssim1$ 成立 → "带电黑洞/带电天体的自加速悖论被消除"，吸引力降一档，但仍是干净的结果。
- 都不成立 → 回到第 2 步的负面结论。

---

### 第 4 步：独立验证（从第 2 步起并行）

**原则**：没有领域内专家独立核对，任何结果都不对外讲。

**做法**

- 第 2 步一出现苗头，就把推导和代码整理成可独立运行的形式。
- 联系做 Reissner–Nordström、自力、带电壳层的研究者，请他们核对关键一步。
- 候选方向：Capra 自力社区；ADM 结果的后续研究者；TheoPhysEval 中修正电动力学题目的出题人（可通过通讯作者 Suvrat Raju 联系）。

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

**记录**

- 所有推导、代码、AI 会话保存在仓库中，按步骤分目录。
- 每一步结束写一页小结：做了什么、结论、是否通过决策点。
- 这份记录本身就是"AI 如何参与物理研究"的一手材料。

---

## 6. 风险清单

| 风险 | 可能性 | 影响 | 应对 |
|---|---|---|---|
| 核心问题已有人做过 | 低（第二轮未找到；Kozameh–Newman 等明确称答案未知） | 高 | 第 0 步已做两轮；第 1、2 步期间持续跟踪新文献 |
| 引力只对 $\xi\lesssim1$ 有效，对电子无效 | 很高（第 2 节静态估算已指向此结论） | 中 | 第 3 步静态部分提前完成；结论主语改为带电天体 |
| 非球对称动力学需要数值相对论 | 中 | 中 | 先做球对称振荡，再决定是否投入 |
| 猜想 B 是 EFT 圈内常识 | 中 | 低 | 作为独立短文，不作为主结论 |
| AI 推导出现难以察觉的错误 | 中 | 高 | §5 AI 使用规范中的验证规则 |
| 找不到专家愿意核对 | 中 | 高 | 从第 2 步起提前联系，不等结果完整 |

---

## 7. 对外发布的条件

只有同时满足以下三条才对外发布：

1. 第 2 步得到正面结果，第 3 步明确了适用范围。
2. 至少一位领域专家书面确认。
3. 结论能用一句外行听得懂的话讲清楚，不出现任何需要解释的术语。

新闻稿模板待写（放在 `docs/` 下）。模板中的"独立验证"和"局限"两段是硬条件，空着不发。

---

## 8. 主要参考文献（待第 0 步逐条核对）

- P. A. M. Dirac, *Proc. R. Soc. A* 167, 148 (1938).
- C. J. Eliezer, *Rev. Mod. Phys.* 19, 147 (1947).
- R. Arnowitt, S. Deser, C. W. Misner, *Phys. Rev.* 120, 313 (1960).
- B. S. DeWitt, R. W. Brehme, *Ann. Phys.* 9, 220 (1960).
- J. M. Hobbs, *Ann. Phys.* 47, 141 (1968).
- A. Schild, *Phys. Rev.* 131, 2762 (1963).
- H. Levine, E. J. Moniz, D. H. Sharp, *Am. J. Phys.* 45, 75 (1977).
- E. J. Moniz, D. H. Sharp, *Phys. Rev. D* 15, 2850 (1977).
- A. G. Smith, C. M. Will, *Phys. Rev. D* 22, 1276 (1980).
- R. C. Ferrell, D. M. Eardley, *Phys. Rev. Lett.* 59, 1617 (1987).
- H. Spohn, *Europhys. Lett.* 50, 287 (2000)；*Dynamics of Charged Particles and Their Radiation Field* (2004).
- A. D. Yaghjian, *Relativistic Dynamics of a Charged Sphere* (Springer).
- S. E. Gralla, A. I. Harte, R. M. Wald, *Phys. Rev. D* 80, 024031 (2009).
- E. Poisson, A. Pound, I. Vega, *Living Rev. Relativ.* 14, 7 (2011).
- P. Zimmerman, E. Poisson, *Phys. Rev. D* 90, 084030 (2014), arXiv:1406.5111.
- T. M. Linz, J. L. Friedman, A. G. Wiseman, *Phys. Rev. D* 90, 084031 (2014), arXiv:1406.5112.
- T. M. Linz, J. L. Friedman, A. G. Wiseman, *Self force on an accelerated particle*, *Phys. Rev. D* 90, 024064 (2014), arXiv:1404.7039.
- TheoPhysEval, arXiv:2610.16733（题 1d、1e：修正电动力学中的辐射阻尼与二体问题）.

---

## 9. 第 0 步文献核查：第一轮结果

> 2026-10-23 · 第一轮检索
> 说明：以下内容来自摘要和检索工具对原文的摘读，**还没有逐篇通读原文**。标【待读】的条目是关键文献，结论依赖它们，需要读全文确认。文献编号、作者和期刊信息也都还没有逐条核对。
> 9.1、9.2 是检索记录；9.3 是对检索结果的评注和交叉检查；9.4、9.5 是由此对路线图做的调整。
> **第二轮读了全文后，本节有几条结论被推翻或降级，见第 10 节。本节保留原样，作为记录。**

### 9.1 核心问题：引力能否消除自加速

**初步结论：有人做过一部分，答案倾向于"引力本身不能消除，视界可以"。**

#### (1) Lun & Chow，arXiv:gr-qc/9409024（1994）【待读】

*The role of the apparent horizon in the evolution of Robinson–Trautman Einstein–Maxwell space-times*

- **研究对象**：Robinson–Trautman Einstein–Maxwell 时空，即完整的非线性广义相对论加电磁场。它描述一个带电、不旋转、同时辐射引力波和电磁波的紧致源。
- **结果**：这类时空的演化方程**存在 runaway 解**。在 Reissner–Nordström 附近线性化，偶极（ℓ = 1）模式的增长率为
  $$\lambda = \frac{3M}{2Q^2}$$
  这**正好等于** LD 方程中自由电荷 runaway 的指数 $1/\tau_0 = 3m/2e^2$（取 $G=c=1$）。
- **视界条件**：如果要求存在过去的表观视界，runaway 被排除。视界存在时，Bondi 质量满足 $M_B^2 \ge Q^2$，即次极端条件。
- **作者猜想**：带视界的这类时空会稳定地演化到 Reissner–Nordström。

**含义**：在完整的 Einstein–Maxwell 理论里，狄拉克方程的 runaway 原封不动地出现了，增长率一模一样。去掉它的不是引力本身，而是视界。视界只在 $Q \le M$ 时存在，而电子的 $Q/M \approx 2\times10^{21}$，远超这个条件。

#### (2) Newman、Kozameh、Adamo 等的系列工作

- Adamo, Newman, Kozameh, *Living Rev. Relativ.*（arXiv:0906.2155）
- Newman, *Gen. Rel. Grav.* 50, 111 (2018)
- 只用 Einstein 方程和 Maxwell 方程，从无穷远处零质量渐近场推出 Abraham–Lorentz–Dirac 辐射反作用力，以及牛顿第二定律、Dirac 旋磁比等经典结果。
- 【待读】他们如何处理 runaway。

#### (3) 超极端 Reissner–Nordström 的稳定性

- **Dotti, Gleiser 等，arXiv:gr-qc/0607052**：超极端（$Q > M$）Reissner–Nordström 裸奇点**线性不稳定**，存在无穷多个精确的不稳定解。【待读：增长率是否等于 $3M/2Q^2$？】
- **Delucchi，*Class. Quantum Grav.*（2026，arXiv:2602.19015）**：把超极端 Reissner–Nordström 当作点电荷的自场。在适当的边界条件下，线性扰动的演化是幺正的，没有束缚态，奇点"动力学上沉默"。未讨论运动方程或 runaway。
- 两者结论不同，差别应该来自奇点处边界条件的选取。

#### (4) ADM 带电球壳的动力学

- 没有找到把 ADM 1960 的静态结果推广到运动方程的工作。arXiv:1303.2706 只做静态。
- 这一块目前看是空白，但需要再查一轮（带电壳层坍缩、Boulware 1973 等）。

### 9.2 狄拉克的选择题（$B^\mu$ 的不唯一性）

**初步结论：有人讨论过歧义本身，但没找到"Eliezer 型推广等价于粒子内部结构"的明确陈述。**

| 文献 | 要点 |
|---|---|
| Kosyakov, *Phys. Rep.* (2019), arXiv:1812.03290 | 认为取 $p^\mu = mv^\mu$ "除了类比力学之外没有依据" |
| Tabensky, *Phys. Rev. D* 13, 267 (1976) | Lorentz 力项与截断方式无关；辐射反作用项依赖截断方式和推导约定 |
| arXiv:hep-th/0405121 | Schott 项不能作为附加项随意加进机械动量 |
| Birnholtz, arXiv:1410.5871 | 在有效场论框架下处理 LD 的初始条件 |
| *Phys. Rev. D* 12, 1576 (1975) | *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*【待读，本轮未能打开】 |
| 世界线有效场论（Goldberger–Rothstein 一路） | 高阶导数项与极化率的对应：未找到直接针对 Dirac/Eliezer 的文章，可能是圈内默认的常识 |

### 9.3 评注：对第一轮结果的交叉检查

#### 9.3.1 关键数字还没有核实

全篇最关键的一个数，是 $\lambda = 3M/2Q^2$ 与 $1/\tau_0$ 正好相等，目前它只来自摘要。量纲是对的：$G=c=1$ 时 $\tau_0 = 2e^2/3m$，$1/\tau_0 = 3m/2e^2$，两边形式一致。但"正好相等"必须读全文确认。

更便宜、也更符合 §5 AI 使用规范中验证规则的做法，是**自己推一遍**：在 RN 背景下写出 RT–Maxwell 方程 ℓ=1 部分的线性化，求特征根。这一步工作量不大，也可以顺带检验 AI 推导的可靠性。

#### 9.3.2 视界：动力学结果，还是选择条件？

Lun–Chow 的条件是**要求存在过去的表观视界**。要读清楚这个条件的地位：

- 如果视界是演化的结果（给定一般初始数据，有视界的解自动没有 runaway 模式），那确实是"视界消除了 runaway"。
- 如果视界是人为加的边界条件，那它和 Dirac 用渐近条件剔除 runaway 是同一类操作：只是换了一种选解的方式。

另外，线性化增长率 $3M/2Q^2$ 本身与 $Q/M$ 的取值无关，次极端时同样存在这个模式。所以问题在于：次极端情形下，这个模式是被方程本身排除的，还是被视界条件排除的。

"自加速和宇宙监督是同一件事"这个说法是强还是弱，取决于这一点。

#### 9.3.3 Q ≤ M 这个巧合可能很平凡，而 ADM 球壳恰好能检验它

**(a) 与 Levine–Moniz–Sharp 的关系。** LMS 的判据是：物体尺寸小于经典半径 $Q^2/M$ 时出现 runaway。有视界时，物体的有效尺寸不小于 $r_+ \ge M$；而
$$M \ge \frac{Q^2}{M} \iff Q \le M .$$
所以 Lun–Chow 的视界条件，可能只是把 LMS 判据里的"尺寸"换成了视界半径。

**(b) ADM 球壳不需要视界，也得到同一条边界。** 由
$$m_0 = m - \frac{q^2}{2R} + \frac{G m^2}{2R}$$
（$c=1$），裸质量为负的条件是
$$R < \frac{q^2 - G m^2}{2m}.$$

- $q \le \sqrt G\, m$（$\xi \le 1$）时，右边非正，**任何半径的球壳裸质量都为正**，按 LMS 判据没有 runaway。
- $q > \sqrt G\, m$（$\xi > 1$）时，才存在会 runaway 的小半径；临界半径在 $\xi \to 1$ 时连续地趋于零。
- 对电子，临界半径 $\approx r_e/2$，引力修正为 $10^{-43}$ 量级（第 2 节）。

**(c) 结论。** 两个模型从不同机制得到同一条边界 $\xi = 1$：

| 模型 | 去除 runaway 的机制 | 边界 |
|---|---|---|
| Robinson–Trautman（Lun–Chow） | 视界 | $Q \le M$ |
| ADM 球壳（若 LMS 判据在 GR 中成立） | 引力束缚能使裸质量为正 | $q \le \sqrt G\, m$ |

因此，ADM 动力学**不是**换个模型把 RT 的结论重复一遍，而是**区分这两种机制的判决性计算**。前提是 LMS 判据在广义相对论中仍然成立，这正是第 2b 步要检验的。

#### 9.3.4 与 Dotti–Gleiser 的比较要具体

- 要比的是 **ℓ=1** 模式的增长率，不是任意不稳定模式。
- 裸奇点处的边界条件不唯一。Delucchi 选的边界条件下没有增长模式，说明增长率本身依赖边界条件。
- 这与 Dirac 的情形结构相同：runaway 在不在，取决于加什么条件（Dirac 的渐近条件 ↔ 奇点处的边界条件）。这个结构对应本身，可能比"增长率是否相等"更有意思。

#### 9.3.5 狄拉克选择题

"可能是 EFT 圈内常识"的判断合理，与第 6 节风险表中的定位（可能性中、影响低）一致。猜想 B 维持原定位：作为独立短文，不作为主结论。

### 9.4 对路线图的影响

**1. 猜想 A 需要改写。** "引力消除 runaway"在 Robinson–Trautman 模型里不成立：引力把 runaway 原样保留了下来。修订为：

> **猜想 A′（runaway 的尺寸判据）**
> runaway 出现，当且仅当物体的尺寸小于它的经典半径。在 Einstein–Maxwell 理论中，"尺寸"可以是物质半径（ADM 球壳），也可以是视界半径（Robinson–Trautman）；两种情况下，引力的作用都是把 runaway 存在与否的边界移到 $\xi = q/(\sqrt G\, m) = 1$。

**2. 一个可能更强的问题：自加速与宇宙监督。**

> runaway 只在没有视界时存在；没有视界，就是超极端，就是裸奇点。
> 那么，**狄拉克方程的自加速和宇宙监督假设，是不是同一件事的两面？**

可检验的一步：超极端 Reissner–Nordström 的 ℓ=1 不稳定模式（Dotti–Gleiser），增长率是否就是 $3M/2Q^2$？如果是，"电子的自加速"就是"裸奇点的不稳定性"。

需要查清 Lun–Chow 或后人是否已经明确写过这个联系。同时要注意 9.3.3：如果 $\xi = 1$ 这条边界只是 LMS 判据的改写，这个联系就可能是平凡的。

**3. ADM 动力学重新定位。** 它仍然是空白，而且按 9.3.3，它的意义比原先设想的更大：它是区分"视界机制"与"裸质量机制"的判决性计算。第 2 步保留，2b 的目标改为：

- 检验 LMS 判据在 GR 中是否成立；
- 重点研究 $\xi \to 1^+$ 附近 runaway 根的行为（临界半径趋于零，几何趋于退化）。

### 9.5 下一轮要做的事

按优先级排列：

1. **读 Lun & Chow 全文**：确认增长率 $3M/2Q^2$；确认视界条件是演化结果还是边界条件（9.3.2）；查它的引用，看谁后来接着做了。
2. **自己推 ℓ=1 增长率**：RN 背景下 RT–Maxwell 方程的线性化，两次独立推导加 SymPy 校验。
3. **读 Dotti–Gleiser**：比较 ℓ=1 不稳定模式的增长率与 $3M/2Q^2$；记录所用的奇点边界条件，与 Delucchi 对比。
4. **查 Newman–Kozameh 如何处理 runaway。**
5. **补读 *PRD* 12, 1576 (1975)。**
6. **再查一轮 ADM 带电球壳的动力学**（带电壳层坍缩、Boulware 1973、Israel 结合条件下的壳层运动）。
7. **核对第 8 节和 9.6 节全部文献的作者、卷号、页码。**

### 9.6 本轮新增参考文献（均待核）

**Einstein–Maxwell 中的 runaway 与视界**
- Lun & Chow, *The role of the apparent horizon in the evolution of Robinson–Trautman Einstein–Maxwell space-times*, arXiv:gr-qc/9409024 (1994).
- T. M. Adamo, E. T. Newman, C. N. Kozameh, *Living Rev. Relativ.*, arXiv:0906.2155.
- E. T. Newman, *Gen. Rel. Grav.* 50, 111 (2018).

**超极端 Reissner–Nordström**
- G. Dotti, R. J. Gleiser 等, arXiv:gr-qc/0607052.
- Delucchi, *Class. Quantum Grav.* (2026), arXiv:2602.19015.

**带电球壳**
- arXiv:1303.2706（ADM 带电球壳，静态）.
- D. G. Boulware (1973)，带电壳层【待查具体文献】.

**$B^\mu$ 的不唯一性与 LD 方程的推导**
- B. P. Kosyakov, *Phys. Rep.* (2019), arXiv:1812.03290.
- R. Tabensky, *Phys. Rev. D* 13, 267 (1976).
- *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975).
- arXiv:hep-th/0405121.
- O. Birnholtz, arXiv:1410.5871.

**世界线有效场论**
- W. D. Goldberger, I. Z. Rothstein, *Phys. Rev. D* 73, 104029 (2006).

---

## 10. 第 0 步文献核查：第二轮结果

> 2026-10-23 · 第二轮检索
> 方法：三路并行检索（arXiv 全文、INSPIRE 引用网络、本地 `references/`），关键原文逐篇通读。每条结论标注核实程度：**【全文】** 读过原文；**【摘要】** 只看过摘要；**【二手】** 来自他人引用。

### 10.0 总的判断

1. **目前没有任何文献表明引力或视界能消除自加速。** 第一轮"引力不能，视界可以"的说法，后半句没有依据。
2. **我们的核心问题是一个被人明确提出、又没人回答的公开问题。** Kozameh–Newman 等人原话："引力是否能稳定运动，是一个开放问题……据我们所知，答案未知。"
3. **没有人推出过带电有尺寸物体在完整 Einstein–Maxwell 理论里的平移运动方程，也没有人把 Levine–Moniz–Sharp 判据推广到含引力。** ADM 球壳动力学仍是空白。
4. 新发现一个具体空白：**超极端 Reissner–Nordström 的 ℓ=1 扇区**，也就是让作为背景的点电荷能运动。

### 10.1 Lun & Chow 1994：逐条核实【全文】

**文献**：A. W. C. Lun, E. W. M. Chow (Monash University), *The role of the apparent horizon in the evolution of Robinson–Trautman Einstein–Maxwell space-times*, arXiv:gr-qc/9409024；报告号 MU-DM-App-94/16。INSPIRE 上无期刊信息；据 Chruściel (arXiv:2112.02126) 的引用，收入会议文集 *Confronting the Infinite* (Adelaide 1994), World Scientific (1995), pp. 279–292【二手】。

| 第一轮说法 | 核实结果 |
|---|---|
| ℓ=1 增长率为 $3M/2Q^2$ | ✅ 成立。由线性化方程推出（原文式 (5)–(7)），不是直接声称 |
| 等于 LD 的 runaway 指数 | ✅ 论文自己指出，并把这一联系归功于 Newman & Posadas, *Phys. Rev.* 187, 1784 (1969) |
| 要求过去表观视界即排除 runaway | ❌ **降级**：视界是**假设**出来的，而且排除 runaway 只是推测，没有证明 |
| 视界存在时 $M_B^2\ge Q^2$ | ✅ 成立（原文 Prop. 2） |
| 只有超极端才有 runaway | ❌ 原文没有这个结论；线性化增长率与 $Q/M$ 无关 |

**完整的特征值**（原文式 (7)）：
$$\lambda_{\ell,\pm} = \frac{3M_0}{8Q_0^2}\,\ell(\ell+1)\left[1 \pm \sqrt{1 + \frac{4Q_0^2(\ell+2)(\ell-1)}{9M_0^2}}\right]$$
ℓ=1 时根号内的耦合项为零，$\lambda_{1,+}=3M_0/2Q_0^2$，$\lambda_{1,-}=0$。原文："all modes with ℓ ≥ 1 of the linearised equations are unstable. The critical point of the system is a saddle."

**视界的地位**（回答 9.3.2 的问题）。原文 §1："we assume the existence of a past apparent horizon"；§5："The requirement for the existence of an apparent horizon amounts to the imposition of a global condition selecting the physically meaningful solutions."

结论：**这是选解条件，和 Dirac 用渐近条件剔除 runaway 是同一类操作。** 视界并没有从动力学上消除 runaway。

**排除 runaway 没有证明。** 他们证明的是：在有视界的假设下 $M_B\ge0$、$M_B^2\ge Q_0^2$，视界面积单调，并构造了一个 Lyapunov 泛函 $L = 16\pi M_B^2 + A_T$。逐点条件 $M - Q_0^2/\Re \ge 0$ "would be sufficient to rule out the exponentially growing modes"，但他们只证明了它在球面上的积分非负。原文："We have not proved the existence of solutions"。

另外，arXiv 元数据里的摘要说"证明了"偶极模式被消除，而 PDF 正文里的摘要没有这句话，正文也只说"details in preparation"。INSPIRE 上没找到后续论文。

**更严重的问题：模型本身不适定。** $\lambda_{\ell,+}$ 随 ℓ 无界增长。Kozameh, Kreiss & Reula, *Class. Quantum Grav.* 25, 025004 (2008), arXiv:0708.1933 证明：带电时 Robinson–Trautman–Maxwell 方程不稳定，"the growth rate increases without bound as the wave number tends to infinity"，初值问题不适定【全文，检索确认该文未提 LD 或 runaway】。Chruściel (arXiv:2112.02126) 也指出两个方程一个对未来抛物、一个对过去抛物。

**含义**：Robinson–Trautman 模型的病不只是偶极 runaway，整个初值问题都有问题。**它不适合用来回答"引力能否消除自加速"。**

**引用网络**（INSPIRE 列出 6 篇）：都是讨论非线性电动力学下 Robinson–Trautman 解的适定性（Tahamtan 2021；Tahamtan, Flores-Alfonso & Svitek 2023；Kokoška & Ortaggio 2021、2022；Nozawa 2024），**没有一篇把它和宇宙监督或裸奇点联系起来**。

### 10.2 Newman–Kozameh 一系：runaway 怎么处理【全文】

**他们得到的不是精确的 LD 方程**，而是慢速、二阶微扰下的类 Abraham–Lorentz 方程。方程写在复世界线（H 空间）上，不在物理时空里。他们**没有降阶，也没有加选解条件**；质量是 Bondi 质量，不需要质量重整化。

- **Kozameh, Newman, Ortega, Silva-Ortigoza, *Gen. Rel. Grav.* 41, 1139 (2009), arXiv:0803.3921**（最关键）
  - 式 (31)：$M v'_R + v' M' - \frac{2Q^2}{3c^3} v''_R + R = F$
  - 结论节：$M'v$ 项"has the correct sign to try to stabilize the motion. Whether or not it does stabilize is an open question……Unfortunately, this is a difficult question and, to our understanding, the answer is unknown."
- **Adamo, Kozameh, Newman, *Living Rev. Relativ.* 12, 6 (2009)；更新版 15, 1 (2012), arXiv:0906.2155**
  - §6.4，式 (6.53)–(6.56)：只有当复电荷中心和复质心两条世界线重合时，$\frac23 q^2\ddot{\vec\xi}_R$ 项才成为辐射反作用力。
  - 原文："The problem of the runaway solutions, though not solved here, is converted to the stability of the Einstein–Maxwell equations with the 'coinciding' condition on the two world lines."
  - §8.5：如果引力不能稳定运动，那么"the Einstein–Maxwell equations are unstable"；另一种可能是电荷中心有自己的动力学。
- **Newman, *Class. Quantum Grav.* 28, 245003 (2011), arXiv:1109.4106**：在 RN 上做到四极。§4："Could the momentum recoil force damp out the run-away behavior? Unlikely?"
- **Newman, *Gen. Rel. Grav.* 50, 111 (2018)**【摘要】：综述，摘要未涉及 runaway。

**含义**：这正是我们的核心问题，而这个领域里最接近它的人明确说答案未知。他们的方法是微扰的、渐近的，没有给出有尺寸物体的内部结构，所以也没法检验 Levine–Moniz–Sharp 判据。**ADM 球壳正好补上这一块。**

### 10.3 超极端 Reissner–Nordström：不稳定性与 ℓ=1【全文】

**文献**
- G. Dotti, R. J. Gleiser, J. Pullin, *Instability of charged and rotating naked singularities*, *Phys. Lett. B* 644, 289 (2007), arXiv:gr-qc/0607052。
- G. Dotti, R. J. Gleiser, *Gravitational instability of the inner static region of a Reissner–Nordström black hole*, *Class. Quantum Grav.* 27, 185007 (2010), arXiv:1001.0152。末节说明方法可直接推广到超极端情形，补全了 2007 年文章的细节。
- Delucchi, *Unitary and finite self-energy of a massive classical point charge as a naked point singularity*, *Class. Quantum Grav.* 43, 195008 (2026), arXiv:2602.19015（PDF 页眉年份误印为 2025）。

**Dotti–Gleiser–Pullin**
- 扇区：Regge–Wheeler 规范下的极性扰动，引力与电磁耦合后化为两个 Zerilli 型场。形式体系只对 ℓ = 2, 3, … 写出。
- 不稳定模式是 Chandrasekhar 的代数特殊模式，增长率为
  $$k = \frac{(\ell-1)\ell(\ell+1)(\ell+2)}{2\left[\sqrt{9M^2+4Q^2(\ell-1)(\ell+2)} - 3M\right]}$$
  原文结论："divergent perturbations … for all ℓ ≥ 2"。
- **ℓ=1 不在其中**：分子 $(\ell-1)$ 为零，$k=0$。
- **与 $3M/2Q^2$ 不符**：ℓ=2、$M\ll Q$ 时 $k\to 3/|Q|$，函数形式不同（这是我们自己的代数，原文没有比较）。
- 奇点处边界条件：要求扰动后的 Kretschmann 标量发散不快于背景的 $56Q^4/r^8$（原文式 (6)），称之为"不使解平凡化的最严格边界条件"。2010 年文章指出该问题属于 limit circle 情形，不同自伴扩张给出不同演化，"physical input must then dictate" 选择。

**Delucchi 2026**
- 不另加边界条件：用 Einstein–Maxwell 的能量形式（Hardy 不等式）取 Friedrichs 扩张。
- 奇点附近两支解为 $x^{1/3}$ 和 $x^{2/3}$，Friedrichs 扩张选 $x^{2/3}$。分解为 $H = A^*A \ge 0$，没有本征值、没有零能共振，所以没有增长模式。
- 附录 D 专门讨论 Dotti–Gleiser：他们的模式属于 $x^{1/3}$ 一支，能量像 $\int x^{-4/3}$ 一样发散。原文："no Dotti–Gleiser instability exists in the action-selected finite-T-energy dynamics, although the profiles remain formal solutions"。
- **ℓ=1 被明确剔除**（原文 §2）："For ℓ = 1, the vector (axial) sector is pure gauge … The scalar (polar) sector describes dipole perturbations … that can be absorbed into a redefinition of the center-of-mass worldline"。也就是说，背景电荷被固定为静止，**没有运动方程**。

**含义**
1. 两篇的分歧只在奇点处边界条件：一个按曲率不变量选，一个按能量有限选。这与 Dirac 的情形结构相同（第 9.3.4 节）。
2. "裸奇点不稳定 = 电子自加速"（第 9.4 节第 2 条）**按现有文献不成立**：已知的不稳定性在 ℓ≥2，ℓ=1 没有，增长率形式也不同。
3. **ℓ=1 扇区正是运动方程所在的地方，两篇都没有做。** 这是一个具体、干净的空白：把超极端 RN 当作能运动的点电荷，在 ℓ=1 扇区推出它的世界线方程，看会不会出现 $3M/2Q^2$ 的 runaway，以及 Delucchi 的能量有限条件会不会把它排除。
4. 相关但不直接的工作：Tahvildar-Zadeh, *Rev. Math. Phys.* 23, 309 (2011), arXiv:1012.1400（单个点电荷的静态时空）【摘要】；Burtscher, Kiessling & Tahvildar-Zadeh, *Class. Quantum Grav.* 38, 185001 (2021), arXiv:1901.00813（类时奇点的弱 Bianchi 恒等式，涉及点电荷运动方程）【摘要】。没有找到推导超极端 RN"电子"运动方程的论文。

### 10.4 ADM 球壳、带电壳层与 Levine–Moniz–Sharp

**引用核实**
- H. Levine, E. J. Moniz, D. H. Sharp, *Motion of extended charges in classical electrodynamics*, *Am. J. Phys.* 45, 75 (1977)：标题、作者确认；"裸质量为负才有 runaway"的表述未从原文确认【摘要缺】。
- E. J. Moniz, D. H. Sharp, *Radiation reaction in nonrelativistic quantum electrodynamics*, *Phys. Rev. D* 15, 2850 (1977)：摘要确认经典理论"shows no runaway solutions or preacceleration when the electron size exceeds the classical electron radius"【摘要】。
- 姊妹篇：Moniz & Sharp, *Absence of runaways and divergent self-mass in nonrelativistic QED*, *Phys. Rev. D* 10, 1133 (1974)。

**广义相对论中的带电壳层**：全部是径向坍缩、反弹或静态能量，**没有一篇处理平移运动或 runaway**。

| 文献 | 内容 | 核实 |
|---|---|---|
| ADM, *Phys. Rev.* 120, 313 (1960) | 静态自能，总质量有限且只依赖电荷 | 【摘要】 |
| Boulware, *Phys. Rev. D* 8, 2363 (1973) | 带电壳层径向演化；坍缩成裸奇点**当且仅当**壳的物质能量密度为负 | 【摘要】 |
| Israel, *Nuovo Cim. B* 44, 1 (1966) | 薄壳结合条件 | 【摘要】 |
| de la Cruz & Israel, *Nuovo Cim. A* 51, 744 (1967) | RN 中带电尘埃壳的径向反弹 | 【二手】 |
| Kuchař, *Czech. J. Phys. B* 18, 435 (1968) | 带压强和电荷的壳层径向运动方程 | 【二手】 |
| Dillon, *Eur. Phys. J. Plus* 128, 72 (2013), arXiv:1303.2706 | 牛顿引力下的球壳模型，复现 ADM 并给出另一个解 | 【全文】 |
| 1988 年博士论文 *Electromagnetic Radiation Reaction in General Relativity*（ADS 1988PhDT........27O） | 外引力场中的刚性带电球体，不含自引力 | 【二手，作者未查到】 |

**Boulware 1973 与猜想 A′ 的关系。** Boulware 的结论是 Levine–Moniz–Sharp 判据（裸质量为负 ⇔ 病态）的**径向版本**，而且是在完整 GR 里得到的。这是猜想 A′ 目前最直接的旁证。需要读全文。

**其他确认**
- Gralla, Harte & Wald, *PRD* 80, 024031 (2009), arXiv:0905.2391：平直时空，小物体标度极限下严格推出降阶的 LD 方程；不含引力【全文】。
- Zimmerman & Poisson (arXiv:1406.5111)、Linz, Friedman & Wiseman (arXiv:1406.5112)：电磁–引力耦合自力，都没有讨论 runaway 或降阶【全文检索】。
- INSPIRE 上 Levine–Moniz–Sharp 有 37 篇引用，ADM 有 136 篇，按引力、辐射反作用、runaway 筛选后，**没有找到把 Levine–Moniz–Sharp 推广到含引力的工作**。

### 10.5 狄拉克选择题

- **E. G. P. Rowe, *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975)**【摘要和第一页，本地 `references/LD equation/1975_Rowe.pdf`】
  - 歧义：用 $\int T^{\mu\nu}dS_\mu = 0$，在一个面上积分得到 LD 方程；在另一个面上积分得到一个**变静质量**的方程，质量变化恰好抵消 Larmor 辐射，且**没有 Schott 项**。两种情况都需要形式相同的无穷质量重整化，Rowe 认为歧义正来自这里。
  - 解决：用分布理论写 $\partial_\mu T^{\mu\nu}=0$，在"很自然的选择"下得到处处有限的表述，LD 方程唯一地出现。
  - 与猜想 B 直接相关：Rowe 说的"自然选择"，可能就对应 $P_{2i}=0$。**需要细读全文**，本地还有 Rowe 1978、1979 两篇。
- **R. Tabensky, *Electrodynamics and the electron equation of motion*, *Phys. Rev. D* 13, 267 (1976)**【摘要】：Lorentz 力与截断无关，辐射反作用力依赖截断。
- **D. V. Gal'tsov, P. Spirin, *Radiation reaction reexamined: bound momentum and Schott term*, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121**【摘要】：Schott 项来自束缚场动量，"can not be introduced by hand as an additional term in the mechanical momentum"。

### 10.6 本地 `references/` 的情况

本地两个文件夹（`LD equation/`、`Motion of Particles/`）主要是平直时空的 LD 方程文献和 Dirac 的原始论文，与引力直接相关的很少。
- J. C. Herrera, *Phys. Rev. D* 21, 384 (1980)：用等效原理把一个**修改过的**运动方程推广到引力场，不是完整的 Einstein–Maxwell 计算。
- P. A. M. Dirac, *An extensible model of the electron*, *Proc. R. Soc. A* 268, 57 (1962)：不含引力。
- 引力自力方面有 Mino、Quinn、Detweiler、Barack、Poisson 等人的文献，是第 1、2 步推导时的工具书。
- 不少文件是扫描件，没有文字层。

### 10.7 对路线图的影响

**1. 猜想 A 第二次修订。** 视界版本撤掉（Robinson–Trautman 模型里视界是假设出来的，而且模型本身不适定），保留裸质量版本：

> **猜想 A″（裸质量判据在 GR 中成立）**
> 一个有尺寸的带电物体出现 runaway，当且仅当它的裸质量为负。在 Einstein–Maxwell 理论中，引力束缚能把裸质量为负的临界半径改为 $R_c = (q^2 - Gm^2)/2m$（ADM 球壳，第 9.3.3 节）。因此 $\xi \le 1$ 的物体对任何尺寸都没有 runaway，$\xi > 1$ 的物体只在 $R < R_c$ 时有。对电子，$R_c \approx r_e/2$，引力修正为 $10^{-43}$。
>
> 旁证：Boulware 1973 在完整 GR 里得到了这个判据的径向版本。

**2. "自加速与宇宙监督"降为备注。** 已知的裸奇点不稳定性在 ℓ≥2，与偶极 runaway 对不上（第 10.3 节）；Robinson–Trautman 的 runaway 与 $Q/M$ 无关。目前没有证据支持这条联系。

**3. 第 2 步（ADM 球壳）地位上升。** 它直接回答 Kozameh–Newman 明确留下的公开问题，而且能检验 Levine–Moniz–Sharp 判据在 GR 中是否成立，这是 Newman 的渐近方法做不到的。

**4. 新增候选方向：超极端 RN 的 ℓ=1 扇区**（可作为第 2 步的并行支线或第 2′ 步）。
- 问题：让 Delucchi 的背景点电荷能运动，在 ℓ=1 扇区推出世界线方程。
- 要看的：是否出现 $3M/2Q^2$ 的 runaway；Delucchi 的能量有限条件能否把它排除；排除机制与 Dirac 渐近条件的关系。
- 优点：点粒子模型，不需要内部结构；有现成的静态背景和边界条件分析可以直接用。

**5. Robinson–Trautman 不再作为主要依据**，只作为"runaway 在完整 GR 中仍然存在"的一个例证。

**6. 猜想 B 新增要读的文献**：Rowe 1975 全文，以及本地的 Rowe 1978、1979。

### 10.8 下一轮要做的事

按优先级排列：

1. **读 Rowe 1975 全文**（本地有），弄清他的"自然选择"与 $P_{2i}=0$ 的关系。为第 1d 步做准备。
2. **读 Boulware 1973 全文**：确认"裸奇点 ⇔ 负能量密度"的准确表述和推导，评估能否推广到平移运动。
3. **读 Levine–Moniz–Sharp 1977 全文**：确认判据原文表述；这是第 1c 步要复现的结果。
4. **读 Delucchi 2026 的 §2 和附录 D**：弄清 ℓ=1 被"吸收进世界线"的具体做法，评估第 4 条新方向的工作量。
5. **自己推 ℓ=1 增长率**（第 9.5 节第 2 条）：重要性下降，因为 Lun–Chow 的推导已核实；改为第 1 步的练手题。
6. 查 1988 年那篇博士论文的作者和内容。
7. 补核所有【摘要】【二手】条目。

### 10.9 本轮新增或修正的参考文献

**修正**
- arXiv:1404.7039 是 T. M. Linz, J. L. Friedman, A. G. Wiseman, *Self force on an accelerated particle*, *Phys. Rev. D* 90, 024064 (2014)，不是 Zimmerman–Poisson（第 1.2、8 节已改正）。
- P. Zimmerman, E. Poisson, *Gravitational self-force in nonvacuum spacetimes*, *Phys. Rev. D* 90, 084030 (2014), arXiv:1406.5111。
- Dotti–Gleiser 2007 的作者还有 J. Pullin，发表在 *Phys. Lett. B* 644, 289。
- *PRD* 12, 1576 (1975) 的作者是 E. G. P. Rowe。

**Robinson–Trautman Einstein–Maxwell**
- A. W. C. Lun, E. W. M. Chow, arXiv:gr-qc/9409024 (1994)；*Confronting the Infinite*, World Scientific (1995), pp. 279–292.
- E. T. Newman, J. Posadas, *Phys. Rev.* 187, 1784 (1969)；*J. Math. Phys.* 12, 2319 (1971).
- C. N. Kozameh, H.-O. Kreiss, O. Reula, *Class. Quantum Grav.* 25, 025004 (2008), arXiv:0708.1933.
- P. T. Chruściel, arXiv:2112.02126.
- T. Tahamtan, *Phys. Rev. D* 103, 064052 (2021), arXiv:2010.01689.
- T. Tahamtan, A. Flores-Alfonso, O. Svitek, *Phys. Rev. D* 108, 124076 (2023), arXiv:2311.03110.

**Newman–Kozameh**
- C. N. Kozameh, E. T. Newman, J. G. Ortega, G. Silva-Ortigoza, *Gen. Rel. Grav.* 41, 1139 (2009), arXiv:0803.3921.
- T. M. Adamo, C. N. Kozameh, E. T. Newman, *Living Rev. Relativ.* 12, 6 (2009)；15, 1 (2012), arXiv:0906.2155.
- E. T. Newman, *Class. Quantum Grav.* 28, 245003 (2011), arXiv:1109.4106.
- E. T. Newman, *Gen. Rel. Grav.* 50, 111 (2018).

**超极端 Reissner–Nordström 与点电荷**
- G. Dotti, R. J. Gleiser, J. Pullin, *Phys. Lett. B* 644, 289 (2007), arXiv:gr-qc/0607052.
- G. Dotti, R. J. Gleiser, *Class. Quantum Grav.* 27, 185007 (2010), arXiv:1001.0152.
- Delucchi, *Class. Quantum Grav.* 43, 195008 (2026), arXiv:2602.19015.
- S. Tahvildar-Zadeh, *Rev. Math. Phys.* 23, 309 (2011), arXiv:1012.1400.
- A. Burtscher, M. K.-H. Kiessling, A. S. Tahvildar-Zadeh, *Class. Quantum Grav.* 38, 185001 (2021), arXiv:1901.00813.

**带电壳层**
- D. G. Boulware, *Naked singularities, thin shells, and the Reissner–Nordström metric*, *Phys. Rev. D* 8, 2363 (1973).
- W. Israel, *Nuovo Cim. B* 44, 1 (1966)；erratum 48, 463 (1967).
- V. de la Cruz, W. Israel, *Nuovo Cim. A* 51, 744 (1967).
- K. Kuchař, *Czech. J. Phys. B* 18, 435 (1968).
- G. Dillon, *Eur. Phys. J. Plus* 128, 72 (2013), arXiv:1303.2706.

**有尺寸电荷与 runaway**
- E. J. Moniz, D. H. Sharp, *Phys. Rev. D* 10, 1133 (1974).
- E. J. Moniz, D. H. Sharp, *Phys. Rev. D* 15, 2850 (1977).

**狄拉克选择题**
- E. G. P. Rowe, *Phys. Rev. D* 12, 1576 (1975).
- R. Tabensky, *Phys. Rev. D* 13, 267 (1976).
- D. V. Gal'tsov, P. Spirin, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121.
