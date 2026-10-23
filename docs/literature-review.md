# 第 0 步文献核查记录

> 这是路线图（`docs/roadmap.md`）第 0 步的过程记录，按检索轮次排列，保留每一轮当时的判断，**后一轮可能推翻前一轮**。现行结论以路线图为准。
>
> 核实程度标注：**【全文】** 读过原文；**【摘要】** 只看过摘要；**【二手】** 来自他人引用；**【待读】** 尚未读。
>
> 编号说明：本文件内"第 1.x、2.x、3.x 节"指本文件的小节。其余"第 N 节""第 N 步"指当时路线图的编号（第一、二轮写于 commit `722f1a7`、`081b942` 的路线图，见 git 历史）；路线图 v2 的步骤编号与之基本一致（第 0–5 步），另新增第 2′ 步。

## 1. 第一轮（摘要级检索）

> 2026-10-23 · 第一轮检索
> 说明：以下内容来自摘要和检索工具对原文的摘读，**还没有逐篇通读原文**。标【待读】的条目是关键文献，结论依赖它们，需要读全文确认。文献编号、作者和期刊信息也都还没有逐条核对。
> 1.1、1.2 是检索记录；1.3 是对检索结果的评注和交叉检查；1.4、1.5 是由此对路线图做的调整。
> **第二轮读了全文后，本节有几条结论被推翻或降级，见第 2 节。本节保留原样，作为记录。**

### 1.1 核心问题：引力能否消除自加速

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

### 1.2 狄拉克的选择题（$B^\mu$ 的不唯一性）

**初步结论：有人讨论过歧义本身，但没找到"Eliezer 型推广等价于粒子内部结构"的明确陈述。**

| 文献 | 要点 |
|---|---|
| Kosyakov, *Phys. Rep.* (2019), arXiv:1812.03290 | 认为取 $p^\mu = mv^\mu$ "除了类比力学之外没有依据" |
| Tabensky, *Phys. Rev. D* 13, 267 (1976) | Lorentz 力项与截断方式无关；辐射反作用项依赖截断方式和推导约定 |
| arXiv:hep-th/0405121 | Schott 项不能作为附加项随意加进机械动量 |
| Birnholtz, arXiv:1410.5871 | 在有效场论框架下处理 LD 的初始条件 |
| *Phys. Rev. D* 12, 1576 (1975) | *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*【待读，本轮未能打开】 |
| 世界线有效场论（Goldberger–Rothstein 一路） | 高阶导数项与极化率的对应：未找到直接针对 Dirac/Eliezer 的文章，可能是圈内默认的常识 |

### 1.3 评注：对第一轮结果的交叉检查

#### 1.3.1 关键数字还没有核实

全篇最关键的一个数，是 $\lambda = 3M/2Q^2$ 与 $1/\tau_0$ 正好相等，目前它只来自摘要。量纲是对的：$G=c=1$ 时 $\tau_0 = 2e^2/3m$，$1/\tau_0 = 3m/2e^2$，两边形式一致。但"正好相等"必须读全文确认。

更便宜、也更符合 §5 AI 使用规范中验证规则的做法，是**自己推一遍**：在 RN 背景下写出 RT–Maxwell 方程 ℓ=1 部分的线性化，求特征根。这一步工作量不大，也可以顺带检验 AI 推导的可靠性。

#### 1.3.2 视界：动力学结果，还是选择条件？

Lun–Chow 的条件是**要求存在过去的表观视界**。要读清楚这个条件的地位：

- 如果视界是演化的结果（给定一般初始数据，有视界的解自动没有 runaway 模式），那确实是"视界消除了 runaway"。
- 如果视界是人为加的边界条件，那它和 Dirac 用渐近条件剔除 runaway 是同一类操作：只是换了一种选解的方式。

另外，线性化增长率 $3M/2Q^2$ 本身与 $Q/M$ 的取值无关，次极端时同样存在这个模式。所以问题在于：次极端情形下，这个模式是被方程本身排除的，还是被视界条件排除的。

"自加速和宇宙监督是同一件事"这个说法是强还是弱，取决于这一点。

#### 1.3.3 Q ≤ M 这个巧合可能很平凡，而 ADM 球壳恰好能检验它

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

#### 1.3.4 与 Dotti–Gleiser 的比较要具体

- 要比的是 **ℓ=1** 模式的增长率，不是任意不稳定模式。
- 裸奇点处的边界条件不唯一。Delucchi 选的边界条件下没有增长模式，说明增长率本身依赖边界条件。
- 这与 Dirac 的情形结构相同：runaway 在不在，取决于加什么条件（Dirac 的渐近条件 ↔ 奇点处的边界条件）。这个结构对应本身，可能比"增长率是否相等"更有意思。

#### 1.3.5 狄拉克选择题

"可能是 EFT 圈内常识"的判断合理，与第 6 节风险表中的定位（可能性中、影响低）一致。猜想 B 维持原定位：作为独立短文，不作为主结论。

### 1.4 对路线图的影响

**1. 猜想 A 需要改写。** "引力消除 runaway"在 Robinson–Trautman 模型里不成立：引力把 runaway 原样保留了下来。修订为：

> **猜想 A′（runaway 的尺寸判据）**
> runaway 出现，当且仅当物体的尺寸小于它的经典半径。在 Einstein–Maxwell 理论中，"尺寸"可以是物质半径（ADM 球壳），也可以是视界半径（Robinson–Trautman）；两种情况下，引力的作用都是把 runaway 存在与否的边界移到 $\xi = q/(\sqrt G\, m) = 1$。

**2. 一个可能更强的问题：自加速与宇宙监督。**

> runaway 只在没有视界时存在；没有视界，就是超极端，就是裸奇点。
> 那么，**狄拉克方程的自加速和宇宙监督假设，是不是同一件事的两面？**

可检验的一步：超极端 Reissner–Nordström 的 ℓ=1 不稳定模式（Dotti–Gleiser），增长率是否就是 $3M/2Q^2$？如果是，"电子的自加速"就是"裸奇点的不稳定性"。

需要查清 Lun–Chow 或后人是否已经明确写过这个联系。同时要注意 1.3.3：如果 $\xi = 1$ 这条边界只是 LMS 判据的改写，这个联系就可能是平凡的。

**3. ADM 动力学重新定位。** 它仍然是空白，而且按 1.3.3，它的意义比原先设想的更大：它是区分"视界机制"与"裸质量机制"的判决性计算。第 2 步保留，2b 的目标改为：

- 检验 LMS 判据在 GR 中是否成立；
- 重点研究 $\xi \to 1^+$ 附近 runaway 根的行为（临界半径趋于零，几何趋于退化）。

### 1.5 下一轮要做的事

按优先级排列：

1. **读 Lun & Chow 全文**：确认增长率 $3M/2Q^2$；确认视界条件是演化结果还是边界条件（1.3.2）；查它的引用，看谁后来接着做了。
2. **自己推 ℓ=1 增长率**：RN 背景下 RT–Maxwell 方程的线性化，两次独立推导加 SymPy 校验。
3. **读 Dotti–Gleiser**：比较 ℓ=1 不稳定模式的增长率与 $3M/2Q^2$；记录所用的奇点边界条件，与 Delucchi 对比。
4. **查 Newman–Kozameh 如何处理 runaway。**
5. **补读 *PRD* 12, 1576 (1975)。**
6. **再查一轮 ADM 带电球壳的动力学**（带电壳层坍缩、Boulware 1973、Israel 结合条件下的壳层运动）。
7. **核对第 8 节和 1.6 节全部文献的作者、卷号、页码。**

### 1.6 本轮新增参考文献（均待核）

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

## 2. 第二轮（逐篇读全文）

> 2026-10-23 · 第二轮检索
> 方法：三路并行检索（arXiv 全文、INSPIRE 引用网络、本地 `references/`），关键原文逐篇通读。每条结论标注核实程度：**【全文】** 读过原文；**【摘要】** 只看过摘要；**【二手】** 来自他人引用。

### 2.0 总的判断

1. **目前没有任何文献表明引力或视界能消除自加速。** 第一轮"引力不能，视界可以"的说法，后半句没有依据。
2. **我们的核心问题是一个被人明确提出、又没人回答的公开问题。** Kozameh–Newman 等人原话："引力是否能稳定运动，是一个开放问题……据我们所知，答案未知。"
3. **没有人推出过带电有尺寸物体在完整 Einstein–Maxwell 理论里的平移运动方程，也没有人把 Levine–Moniz–Sharp 判据推广到含引力。** ADM 球壳动力学仍是空白。
4. 新发现一个具体空白：**超极端 Reissner–Nordström 的 ℓ=1 扇区**，也就是让作为背景的点电荷能运动。

### 2.1 Lun & Chow 1994：逐条核实【全文】

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

**视界的地位**（回答 1.3.2 的问题）。原文 §1："we assume the existence of a past apparent horizon"；§5："The requirement for the existence of an apparent horizon amounts to the imposition of a global condition selecting the physically meaningful solutions."

结论：**这是选解条件，和 Dirac 用渐近条件剔除 runaway 是同一类操作。** 视界并没有从动力学上消除 runaway。

**排除 runaway 没有证明。** 他们证明的是：在有视界的假设下 $M_B\ge0$、$M_B^2\ge Q_0^2$，视界面积单调，并构造了一个 Lyapunov 泛函 $L = 16\pi M_B^2 + A_T$。逐点条件 $M - Q_0^2/\Re \ge 0$ "would be sufficient to rule out the exponentially growing modes"，但他们只证明了它在球面上的积分非负。原文："We have not proved the existence of solutions"。

另外，arXiv 元数据里的摘要说"证明了"偶极模式被消除，而 PDF 正文里的摘要没有这句话，正文也只说"details in preparation"。INSPIRE 上没找到后续论文。

**更严重的问题：模型本身不适定。** $\lambda_{\ell,+}$ 随 ℓ 无界增长。Kozameh, Kreiss & Reula, *Class. Quantum Grav.* 25, 025004 (2008), arXiv:0708.1933 证明：带电时 Robinson–Trautman–Maxwell 方程不稳定，"the growth rate increases without bound as the wave number tends to infinity"，初值问题不适定【全文，检索确认该文未提 LD 或 runaway】。Chruściel (arXiv:2112.02126) 也指出两个方程一个对未来抛物、一个对过去抛物。

**含义**：Robinson–Trautman 模型的病不只是偶极 runaway，整个初值问题都有问题。**它不适合用来回答"引力能否消除自加速"。**

**引用网络**（INSPIRE 列出 6 篇）：都是讨论非线性电动力学下 Robinson–Trautman 解的适定性（Tahamtan 2021；Tahamtan, Flores-Alfonso & Svitek 2023；Kokoška & Ortaggio 2021、2022；Nozawa 2024），**没有一篇把它和宇宙监督或裸奇点联系起来**。

### 2.2 Newman–Kozameh 一系：runaway 怎么处理【全文】

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

### 2.3 超极端 Reissner–Nordström：不稳定性与 ℓ=1【全文】

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
1. 两篇的分歧只在奇点处边界条件：一个按曲率不变量选，一个按能量有限选。这与 Dirac 的情形结构相同（第 1.3.4 节）。
2. "裸奇点不稳定 = 电子自加速"（第 1.4 节第 2 条）**按现有文献不成立**：已知的不稳定性在 ℓ≥2，ℓ=1 没有，增长率形式也不同。
3. **ℓ=1 扇区正是运动方程所在的地方，两篇都没有做。** 这是一个具体、干净的空白：把超极端 RN 当作能运动的点电荷，在 ℓ=1 扇区推出它的世界线方程，看会不会出现 $3M/2Q^2$ 的 runaway，以及 Delucchi 的能量有限条件会不会把它排除。
4. 相关但不直接的工作：Tahvildar-Zadeh, *Rev. Math. Phys.* 23, 309 (2011), arXiv:1012.1400（单个点电荷的静态时空）【摘要】；Burtscher, Kiessling & Tahvildar-Zadeh, *Class. Quantum Grav.* 38, 185001 (2021), arXiv:1901.00813（类时奇点的弱 Bianchi 恒等式，涉及点电荷运动方程）【摘要】。没有找到推导超极端 RN"电子"运动方程的论文。

### 2.4 ADM 球壳、带电壳层与 Levine–Moniz–Sharp

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

### 2.5 狄拉克选择题

- **E. G. P. Rowe, *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975)**【摘要和第一页，本地 `references/LD equation/1975_Rowe.pdf`】
  - 歧义：用 $\int T^{\mu\nu}dS_\mu = 0$，在一个面上积分得到 LD 方程；在另一个面上积分得到一个**变静质量**的方程，质量变化恰好抵消 Larmor 辐射，且**没有 Schott 项**。两种情况都需要形式相同的无穷质量重整化，Rowe 认为歧义正来自这里。
  - 解决：用分布理论写 $\partial_\mu T^{\mu\nu}=0$，在"很自然的选择"下得到处处有限的表述，LD 方程唯一地出现。
  - 与猜想 B 直接相关：Rowe 说的"自然选择"，可能就对应 $P_{2i}=0$。**需要细读全文**，本地还有 Rowe 1978、1979 两篇。
- **R. Tabensky, *Electrodynamics and the electron equation of motion*, *Phys. Rev. D* 13, 267 (1976)**【摘要】：Lorentz 力与截断无关，辐射反作用力依赖截断。
- **D. V. Gal'tsov, P. Spirin, *Radiation reaction reexamined: bound momentum and Schott term*, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121**【摘要】：Schott 项来自束缚场动量，"can not be introduced by hand as an additional term in the mechanical momentum"。

### 2.6 本地 `references/` 的情况

本地两个文件夹（`LD equation/`、`Motion of Particles/`）主要是平直时空的 LD 方程文献和 Dirac 的原始论文，与引力直接相关的很少。
- J. C. Herrera, *Phys. Rev. D* 21, 384 (1980)：用等效原理把一个**修改过的**运动方程推广到引力场，不是完整的 Einstein–Maxwell 计算。
- P. A. M. Dirac, *An extensible model of the electron*, *Proc. R. Soc. A* 268, 57 (1962)：不含引力。
- 引力自力方面有 Mino、Quinn、Detweiler、Barack、Poisson 等人的文献，是第 1、2 步推导时的工具书。
- 不少文件是扫描件，没有文字层。

### 2.7 对路线图的影响

**1. 猜想 A 第二次修订。** 视界版本撤掉（Robinson–Trautman 模型里视界是假设出来的，而且模型本身不适定），保留裸质量版本：

> **猜想 A″（裸质量判据在 GR 中成立）**
> 一个有尺寸的带电物体出现 runaway，当且仅当它的裸质量为负。在 Einstein–Maxwell 理论中，引力束缚能把裸质量为负的临界半径改为 $R_c = (q^2 - Gm^2)/2m$（ADM 球壳，第 1.3.3 节）。因此 $\xi \le 1$ 的物体对任何尺寸都没有 runaway，$\xi > 1$ 的物体只在 $R < R_c$ 时有。对电子，$R_c \approx r_e/2$，引力修正为 $10^{-43}$。
>
> 旁证：Boulware 1973 在完整 GR 里得到了这个判据的径向版本。

**2. "自加速与宇宙监督"降为备注。** 已知的裸奇点不稳定性在 ℓ≥2，与偶极 runaway 对不上（第 2.3 节）；Robinson–Trautman 的 runaway 与 $Q/M$ 无关。目前没有证据支持这条联系。

**3. 第 2 步（ADM 球壳）地位上升。** 它直接回答 Kozameh–Newman 明确留下的公开问题，而且能检验 Levine–Moniz–Sharp 判据在 GR 中是否成立，这是 Newman 的渐近方法做不到的。

**4. 新增候选方向：超极端 RN 的 ℓ=1 扇区**（可作为第 2 步的并行支线或第 2′ 步）。
- 问题：让 Delucchi 的背景点电荷能运动，在 ℓ=1 扇区推出世界线方程。
- 要看的：是否出现 $3M/2Q^2$ 的 runaway；Delucchi 的能量有限条件能否把它排除；排除机制与 Dirac 渐近条件的关系。
- 优点：点粒子模型，不需要内部结构；有现成的静态背景和边界条件分析可以直接用。

**5. Robinson–Trautman 不再作为主要依据**，只作为"runaway 在完整 GR 中仍然存在"的一个例证。

**6. 猜想 B 新增要读的文献**：Rowe 1975 全文，以及本地的 Rowe 1978、1979。

### 2.8 下一轮要做的事

按优先级排列：

1. **读 Rowe 1975 全文**（本地有），弄清他的"自然选择"与 $P_{2i}=0$ 的关系。为第 1d 步做准备。
2. **读 Boulware 1973 全文**：确认"裸奇点 ⇔ 负能量密度"的准确表述和推导，评估能否推广到平移运动。
3. **读 Levine–Moniz–Sharp 1977 全文**：确认判据原文表述；这是第 1c 步要复现的结果。
4. **读 Delucchi 2026 的 §2 和附录 D**：弄清 ℓ=1 被"吸收进世界线"的具体做法，评估第 4 条新方向的工作量。
5. **自己推 ℓ=1 增长率**（第 1.5 节第 2 条）：重要性下降，因为 Lun–Chow 的推导已核实；改为第 1 步的练手题。
6. 查 1988 年那篇博士论文的作者和内容。
7. 补核所有【摘要】【二手】条目。

### 2.9 本轮新增或修正的参考文献

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
- J. F. Donoghue, *General relativity as an effective field theory: the leading quantum corrections*, *Phys. Rev. D* 50, 3874 (1994)【待核】.
- E. J. Moniz, D. H. Sharp, *Phys. Rev. D* 15, 2850 (1977).

**狄拉克选择题**
- E. G. P. Rowe, *Phys. Rev. D* 12, 1576 (1975).
- R. Tabensky, *Phys. Rev. D* 13, 267 (1976).
- D. V. Gal'tsov, P. Spirin, *Grav. Cosmol.* 12, 1 (2006), arXiv:hep-th/0405121.

---

## 3. 第三轮（专项精读与文献补核）

> 2026-10-23 · 第三轮
> 内容：按第 2.8 节的清单，精读 Rowe 1975/1978/1979（本地）、Delucchi 2026 的 ℓ=1 处理与附录 D；用 Chrome 尝试获取 Boulware 1973 全文；用 Crossref/INSPIRE/arXiv 补核全部【摘要】【二手】条目的出版信息。

### 3.0 本轮要点

1. **Rowe 1975 的歧义不是 Eliezer 型**：是"LD + Schott 项"与"变静质量、无 Schott 项"之争。Rowe 1978 的量纲论证说明场张量里不可能出现 Eliezer 型高阶项，因为经典电磁学没有可用的长度常数。这部分支持猜想 B，但 $p=mv$ 在 Rowe 那里是假设的。
2. **Delucchi 对 ℓ=1 的剔除论证有漏洞**：Einstein–Maxwell 的 ℓ=1 中有真实的电磁偶极自由度，不是纯规范。超极端 RN 的 ℓ=1 扇区仍是空白，且有现成工具（Burko 1995、1999）可借用。
3. **Levine–Moniz–Sharp 的原文判据是"尺寸大于经典半径"，摘要中没有"裸质量"**。平直时空中二者等价，加入引力后不再等价。猜想 A 必须明确选哪个判据，这正是第 2 步要判决的。
4. **Boulware 1973 全文需要机构访问**（APS 付费墙），本轮只核实了摘要原文。
5. 文献出版信息基本全部核实；改正三处（Newman–Posadas 的作者名、Hobbs 1968 的篇目、de la Cruz–Israel 的首页）。

### 3.1 Rowe 1975、1978、1979【全文，本地扫描件】

**E. G. P. Rowe, *Resolution of an ambiguity in the derivation of the Lorentz–Dirac equation*, *Phys. Rev. D* 12, 1576 (1975)**

- **歧义**：动量守恒取 $\int dS_\mu T^{\mu\nu}=0$（式 (2)），$T = K + \Theta$。比较两种闭合面（图 1）：
  - (a) Dirac 的细管，平端面：得到 LD 方程（式 (1)），$m\dot v = eF_{\rm ext}\cdot v + \tfrac23 e^2(\dot a - a^2 v)$，质量重整化 $m = m_0 + e^2/2\eta$。
  - (b) 两端为未来光锥的有限管：得到式 (3)，$\dfrac{d}{d\tau}(mv^\mu) = eF_{\rm ext}^{\mu\nu}v_\nu - \tfrac23 e^2 a^2 v^\mu$，重整化形式相同。与 $v$ 缩并得 $m\dot v = eF_{\rm ext}\cdot v$，$\dot m = -\tfrac23 e^2 a^2$：静质量的变化恰好抵消 Larmor 辐射，没有 Schott 项。
  - Rowe 的诊断：矛盾来自质量重整化，"the finite Schott term is lost between the canceling infinities"（p. 1582）。
- **解决**：把 $\partial_\mu T^{\mu\nu}=0$ 写成分布方程。$\Theta_{\rm I+II}$（$\rho^{-4}$、$\rho^{-3}$ 奇异）需要正则化；"很自然的选择"是仿照 Gel'fand–Shilov 对 $x^{-1}$、$x^{-2}$ 的典范正则化（Hadamard 有限部分），推广到光锥坐标 $(\tau,\rho,u)$ 中的 $r^{-3}$、$r^{-4}$（式 (52)–(61)）。结果 $\partial\cdot\Theta_{\rm I+II} = -\tfrac23 e^2\dot a\,\delta$（式 (67)），由此得到 LD（式 (68)），"No mass renormalization of any sort … is needed"。
- **残留的任意性**：Rowe 承认"to any regularization … one can add a distribution concentrated on the world line"，需要"a convincing physical or mathematical principle to cut through this arbitrariness"。1975 年文中没有分类这些附加项。
- **其他**：粒子动量 $p = mv$（式 (42)，Minkowski 粒子能动张量，取实验质量）是**假设**；runaway 与预加速只在引言提及，没有分析。Schott 项属于场："There is certainly no need to introduce an acceleration energy concept, intrinsic to the particle"（p. 1587）。

**E. G. P. Rowe, *Structure of the energy tensor in the classical electrodynamics of point particles*, *Phys. Rev. D* 18, 3639 (1978)**

- 把 $\theta_{\rm ret}$ 分解为 $\theta_1$（自能，处处无散度）、$\theta_2$（Schott 张量，任意类空面上的动量为 $-\tfrac23 e^2 a$）、$\theta_3$（辐射）。
- **§VIII 唯一性论证**（p. 3651）：按量纲，可以加在世界线上的项只能由 $(v,a)$ 和 $(v,\partial)$ 组成；更高阶导数项被排除，因为"There are no physical constants, available within classical electromagnetism, of dimension L"，且"Since θ is purely electromagnetic, it cannot depend on the mass of the particle, so that in classical theory $e^2/m$ (L) is excluded"。剩下的两种附加项会让平面动量依赖于平面的倾角，从而"no equation at all"。结论："the distribution form of $\theta_{\rm ret}$ is uniquely determined"。

**E. G. P. Rowe, *Angular momentum radiated by classical point particles*, *Phys. Rev. D* 20, 580 (1979)**：把 1978 年的分解用于角动量，与动量歧义无关。

**评注（与猜想 B 的关系）**

- 高把握：Rowe 1975 的歧义不是 Eliezer 型。他的 $B^\mu$ 始终是 $mv$（式 (3) 中 $m$ 可变），从未出现 $mv + k\ddot v$。
- 中等把握：$P_{2i}=0$ 在 Rowe 那里是**假设**进粒子张量的，不是推出来的。他的唯一性论证只针对场张量。但他的量纲论证换个说法就是：**Eliezer 型的项需要额外的长度尺度 $e^2/m$，而纯电磁的场张量不能提供它**，只能来自粒子本身，也就是内部结构。这部分支持猜想 B。
- 未仔细检验的推断：1978 年的"倾角依赖"论证似乎也适用于粒子张量。形如 $\int(vw+wv)\delta$、$w_\perp\ne0$ 的世界线项，其散度含 $\delta$ 的横向导数，除非给粒子加上偶极型（结构）项，否则"no equation"。这值得在第 1d 步自己推一遍。
- Rowe 没有讨论点极限，也没有讨论降阶或与 Landau–Lifshitz 的比较。

### 3.2 Delucchi 2026：ℓ=1 与附录 D【全文】

**ℓ=1 是怎么被去掉的**（§3"Master field reduction"，p. 9；不是第 2.3 节记录的 §2）

- 全部论证只有一段："For ℓ = 1, the vector (axial) sector is pure gauge. (The axial gauge-invariant variable vanishes once the constraints are imposed.) The scalar (polar) sector describes dipole perturbations of the Maxwell field and metric that can be absorbed into a redefinition of the center-of-mass worldline and of the global Maxwell gauge."
- 没有写出规范矢量或扰动函数，只引用 Kodama–Ishibashi 的 §3、§6.1（脚注 6）。
- 全文是真空扰动，没有源；背景库仑场被解释为 $r=0$ 处的分布源。§6 定义的"fixed-background single-source sector"按定义商掉了 ℓ = 0、1 的规范与参数变化。
- §7 进一步说："for ℓ = 0, 1 the corresponding sectors are pure gauge and do not carry propagating gravitational or electromagnetic degrees of freedom"。ℓ=1 的辐射从未讨论；能量求和（式 (19)）只对 ℓ≥2。

**附录 D：与 Dotti–Gleiser 的比较**

- 顶点附近 $V_\pm \approx C_\pm x^{-2}$，$C_+=-2/9$，$C_-=4/9$，对应 $\nu_+=1/6$、$\nu_-=5/6$（式 (176)–(178)）。局部行为 $x^{1/2\mp\nu}$：极性通道 $x^{1/3}$ 或 $x^{2/3}$，轴向通道 $x^{-1/3}$ 或 $x^{4/3}$（式 (180)–(181)）。
- 局部能量 $\int_0^\epsilon x^{2\gamma-2}dx$ 仅当 $\gamma>1/2$ 时有限（式 (186)–(187)），保留 $x^{2/3}$ 和 $x^{4/3}$。
- Dotti–Gleiser 的不稳定模式 $\sim x^{1/3}$，属于更奇异的一支，能量 $\sim 3\epsilon^{-1/3}\to\infty$（式 (199)）："formal solution"，但在物理 Hilbert 空间之外。
- 文本里的两处不一致：§3 定义 $\varphi_+$ 为轴向、$\varphi_-$ 为极性，附录 D 反过来；附录 D 写"$x\sim r^{1/3}$"，按附录 A 的 $x\sim r^3/3Q^2$ 应为 $r\sim x^{1/3}$。

**能量泛函与自伴扩张**：每个模式 $E = \tfrac12\|\dot\varphi\|^2 + q_\pm[\varphi]$，$q_\pm[u] = \tfrac12\int(|u'|^2 + V_\pm|u|^2)dx$；用半直线上的锐 Hardy 不等式 $\int|u'|^2 \ge \tfrac14\int|u|^2/x^2$（关键余量 $-2/9 = -1/4 + 1/36$），在 $H^1_0$ 上闭合，取 Friedrichs 扩张，$u(0)=0$。端点是 limit circle，存在其他扩张，作者以"不是由 Einstein–Maxwell 能量本身产生的端点结构"为由排除。对 ℓ=1 未作任何说明。

**展望（§8）**：Kerr–Newman、极端极限、非线性稳定性、量子化。没有提运动方程、辐射反作用或 ℓ=1。

**评注**

- **"ℓ=1 是纯规范"在 Einstein–Maxwell 中不成立**。真空 GR 里 ℓ=1 确实只有规范和守恒量；但 Maxwell 场在两种宇称下都有真实的偶极辐射。$Q\ne0$ 时引力部分可规范掉，每种宇称仍剩一个电磁规范不变量。
- **轴向 ℓ=1 势**（我们自己的推导，未与文献对照）：在 Chandrasekhar 的 RN 轴向势 $V_i = \dfrac{\Delta}{r^5}\Big[(\mu^2+2)r - q_j\Big(1+\dfrac{q_i}{\mu^2 r}\Big)\Big]$（$q_1+q_2=6M$，$q_1q_2=-4Q^2\mu^2$）中取 $\mu^2=(\ell-1)(\ell+2)\to0$，电磁一支（$q_i=q_1\to6M$，$q_j=q_2\to0$）给出
  $$V = f(r)\left(\frac{2}{r^2}+\frac{4Q^2}{r^4}\right),\qquad f = 1-\frac{2M}{r}+\frac{Q^2}{r^2}.$$
  已手算复核。$Q=0$ 时退化为 Maxwell 偶极势 $2f/r^2$。$Q>M$ 时处处为正，顶点附近 $V\approx 4/9x^2$（$\nu=5/6$），与 Delucchi 已处理的通道相同，直接套用 Friedrichs 机制即无增长模式。
- **极性 ℓ=1 是世界线所在之处**：平移规范矢量给出"质心"模式，剩下的电偶极规范不变量是入射偶极波摇动电荷的通道（经典 Thomson 散射的结构）。推测（未检验）：$3M/2Q^2$ 的 runaway 在 Friedrichs 扩张下不出现，在非 Friedrichs 扩张下可能出现。若如此，runaway 问题就等价于"选哪个自伴扩张"，与 Dirac 的渐近条件、Lun–Chow 的视界条件同构。
- **可借用的文献**（均【二手】，待读）：
  - Zerilli 1974（*PRD* 9, 860）：RN 上耦合扰动的形式体系。
  - L. M. Burko：Zerilli/Chandrasekhar 形式在 ℓ=1 失效，推广到极性偶极（arXiv:gr-qc/9507012, 1995）和轴向偶极（arXiv:gr-qc/9902002, 1999）。黑洞区间。
  - arXiv:2010.12984 附录 C：K=0 规范下的偶宇称 ℓ=1。
  - Bini, Geralico, Ruffini（arXiv:gr-qc/0609041、gr-qc/0608139）：RN 上带 δ 源的带电粒子，含 ℓ=1，黑洞区间。
  - Moncrief 1974–75 的 Hamilton 约化（Delucchi 的参考文献 [4]）。
  - **没有找到 $Q>M$ 的 ℓ=1 分析，也没有找到极性 ℓ=1 主势的闭合形式。**
- **工作量估计**：轴向 ℓ=1 约一天；极性 ℓ=1 真空主方程（$Q>M$）几天；加入世界线源（或把 ℓ=1 规范矢量 $\xi(t)$ 当集体坐标、由 ℓ=1 约束或无穷远处 ADM 动量平衡读出其方程）几周；最后求各自伴扩张的谱或散射矩阵极点，检验 $\kappa = 3M/2Q^2$。

### 3.3 Boulware 1973 与 Levine–Moniz–Sharp 1977

**Boulware 1973**：用 Chrome 打开 APS 页面，全文需要机构登录（"Authorization Required"），未登录、未购买。摘要原文已核实：

> "The time development of a thin charged shell is studied and it is shown that it can collapse to form a naked singularity if, and only if, the matter energy density of the shell is negative."

APS 页面显示被引 116 次（Citing Articles 90），可作为下一轮的引用网络入口。

**Levine–Moniz–Sharp 1977**（*Am. J. Phys.* 45, 75–78）：摘要原文（Crossref）说理论"internally consistent (i.e., no runaways or acausality) whenever the size of the particle exceeds the classical radius (defined as the radius for which the electrostatic self-energy equals the mass of the particle)"；运动方程是微分–差分方程。**摘要中没有"裸质量"。**

**评注：两种判据在 GR 中分道扬镳**

- 平直时空：尺寸等于经典半径 ⇔ 静电自能等于质量 ⇔ 裸质量为零。两种说法等价。
- 加入引力（ADM 球壳，$c=1$）：
  - 裸质量判据：$m_0<0 \iff R < R_{m_0} = \dfrac{q^2 - Gm^2}{2m}$；
  - 尺寸判据（按 LMS 字面定义，静电自能 $q^2/2R$ 等于 $m$）：$R < R_{\rm cl} = \dfrac{q^2}{2m}$。
- $\xi\le1$ 时两者预言完全不同：裸质量判据说任何尺寸都没有 runaway；尺寸判据说 $R<q^2/2m$ 仍会 runaway。对电子两者相差 $10^{-43}$，无法区分。
- **旧版路线图第 1.1 节的表格和本文件第 1.3.3 节把 LMS 写成"裸质量为负才有 runaway"，是我们自己的转述，不是原文。** 哪个判据在 GR 中正确，是第 2 步要判决的问题。

### 3.4 文献补核结果（Crossref / INSPIRE / arXiv）

| 条目 | 状态 | 说明 |
|---|---|---|
| ADM, *Phys. Rev.* 120, 313–320 (1960) | 已核 | 标题 *Gravitational-electromagnetic coupling and the classical self-energy problem*；摘要："finite (although still quite large), and … depends only on the charge and not on the bare mechanical mass" |
| Boulware, *PRD* 8, 2363–2368 (1973) | 已核 | 摘要见 3.3 |
| Israel, *Nuovo Cim. B* 44, 1–14 (1966)；erratum 48, 463 (1967) | 已核 | |
| de la Cruz & Israel, *Nuovo Cim. A* 51, 744–760 (1967) | **改正** | 首页是 744（第二轮有一处来源写 745） |
| Kuchař, *Czech. J. Phys. B* 18, 435–463 (1968) | 已核 | |
| Levine, Moniz, Sharp, *Am. J. Phys.* 45, 75–78 (1977) | 已核 | 判据是尺寸，见 3.3 |
| Moniz & Sharp, *PRD* 10, 1133–1136 (1974) | 已核 | *Absence of runaways and divergent self-mass in nonrelativistic QED*：点电子在量子力学中自质量为零，无 runaway |
| Moniz & Sharp, *PRD* 15, 2850–2865 (1977) | 已核 | 经典理论在尺寸大于经典电子半径时无 runaway、无预加速 |
| Tabensky, *PRD* 13, 267–273 (1976) | 已核 | |
| Newman, *GRG* 50, 111 (2018) | 已核 | 标题 *Classical mechanics via general relativity and Maxwell's theory: a bit of magic* |
| Newman & Posadas, *Phys. Rev.* 187, 1784–1791 (1969)；*J. Math. Phys.* 12, 2319–2327 (1971) | **改正** | 第二作者是 **R.**（Roger）Posadas；标题 *Motion and structure of singularities in general relativity*（I、II） |
| Donoghue, *PRD* 50, 3874–3888 (1994), arXiv:gr-qc/9405057 | 已核 | |
| Tahvildar-Zadeh, *Rev. Math. Phys.* 23, 309–346 (2011) | 已核 | 一类非线性电动力学中，静态点电荷时空只有锥形奇点，场能有限且等于 ADM 质量 |
| Burtscher, Kiessling, Tahvildar-Zadeh, *CQG* 38, 185001 (2021) | 已核 | 带类时奇点的静态球对称时空中，二次缩并 Bianchi 恒等式在分布意义下成立 |
| Spohn, *Europhys. Lett.* 50, 287–292 (2000), arXiv:physics/9911027 | 已核 | *The critical manifold of the Lorentz–Dirac equation* |
| Gralla, Harte, Wald, *PRD* 80, 024031 (2009) | 已核 | |
| Smith & Will, *PRD* 22, 1276–1284 (1980) | 已核 | |
| Ferrell & Eardley, *PRL* 59, 1617–1620 (1987) | 已核 | |
| Schild, *Phys. Rev.* 131, 2762–2766 (1963) | 已核 | |
| Hobbs, *Ann. Phys.* 47, 141–165 (1968) | **改正** | 这一篇是 *A vierbein formalism of radiation damping*；同卷 166–172 页是 *Radiation damping in conformally flat universes*，勿混 |
| DeWitt & Brehme, *Ann. Phys.* 9, 220–259 (1960) | 已核 | |
| Eliezer, *Rev. Mod. Phys.* 19, 147–184 (1947) | 已核 | |
| Dirac, *Proc. R. Soc. A* 167, 148–169 (1938) | 已核 | |
| Poisson, Pound, Vega, *Living Rev. Relativ.* 14, 7 (2011), arXiv:1102.0529 | 已核 | |
| Yaghjian, *Relativistic Dynamics of a Charged Sphere* | 已核 | Springer，1992、2006、2022 三个版本 |
| 1988 年博士论文 *Electromagnetic Radiation Reaction in General Relativity* | 部分 | 可能是 Nuala O'Donnell, Cardiff University, 导师 B. F. Schutz；仅有一个二手来源，未能确认就是 ADS 1988PhDT........27O |
| TheoPhysEval, arXiv:2610.16733 | 已核 | A. Ramesh Chandra 等 35 人，*TheoPhysEval: Assessing AI agents at autonomous theoretical physics problem solving*；题 1d 为修正电动力学中的辐射阻尼，题 1e 为修正电动力学中的二体问题；S. Raju 是唯一通讯作者 |

### 3.5 对路线图的影响

这一轮之后，路线图整体重写为 v2（2026-10-23）。主要变化：

1. 猜想 A 改为在 GR 中**区分两种判据**（裸质量 vs 尺寸），以第 2 步判决。
2. 猜想 B 加入 Rowe 1978 的量纲论证作为部分支持。
3. 新增第 2′ 步：超极端 RN 的 ℓ=1 扇区。
4. "选解条件"作为贯穿的线索写进路线图：Dirac 的渐近条件、Lun–Chow 的视界条件、Delucchi 的 Friedrichs 扩张。
5. 参考文献合并去重，按本节改正。

### 3.6 仍未完成

- Boulware 1973、Levine–Moniz–Sharp 1977 的全文（需要机构访问）。
- 1988 年博士论文的确认。
- 第 3.2 节列出的 ℓ=1 工具文献（Zerilli、Burko、Kodama–Ishibashi、Moncrief）尚未读。
- Kosyakov 2019、Birnholtz 2014、Goldberger–Rothstein 2006、Chruściel 2021 的出版信息未在本轮补核。

### 3.7 补记：第 3.3 节评注的改正

第 3.3 节说"$\xi\le1$ 时两者预言完全不同"，**这是错的**。静态球壳必须在视界外，$R > r_+ \ge m \ge q^2/m > R_{\rm cl}$，所以 $\xi\le1$ 时任何可达的球壳都同时满足两种判据，两者一致地预言没有 runaway。两种判据真正分歧的区域在 $\xi>1$：$R_{m_0} < R < R_{\rm cl}$，带宽 $m/2$。路线图 v2.1 第 3.1 节已改正。


---

## 4. 第四轮（2010–2026 年的现代文献）

> 2026-10-23 · 第四轮
> 起因：前三轮主要从老文献往回追参考文献，结果偏老。本轮改为往前追引用（谁引用了 Gralla–Harte–Wald、Spohn 2000、Levine–Moniz–Sharp、ADM、Kozameh–Newman），并专查三个现代领域：数学物理（Kiessling 一系）、弱引力猜想与超极端物质、有尺寸电荷与世界线有效场论。四路并行，主要用 INSPIRE（arXiv 和 Semantic Scholar 的接口中途限流）。
> 关键原话已由主线程在原文中复核：Steane 2015 关于引力的一节和"no known physical object"一句；Galley–Leibovich–Rothstein 2012 回复中"worldline shifts"一段。

### 4.0 本轮结论

1. **不是半个世纪没人想，而是主流认为对电子已经解决。** 平直时空中：有尺寸物体裸质量为正就没有 runaway（精确结果）；LD/LL 是它的近似（Gralla–Harte–Wald 2009 严格化）；电子交给 QED。引力被舍弃，是因为对电子只有 $10^{-43}$（Steane 2015 明确这样论证）。
2. **2010–2026 年没有人问过或算过"自引力能否改变有尺寸带电体的 runaway 判据"**，包括 $G$ 一阶的后牛顿计算。$\xi>1$ 判决带也没人做过。
3. **GR 中的点电荷：黑洞或负裸质量**（Kiessling–Tahvildar-Zadeh 2019 引 Geroch–Traschen；Ponce de León 2004）。
4. **猜想 B 没有人完整写过，但零件都已发表**（Steane 2015；Galley–Leibovich–Rothstein 2012 回复；Sokolov 2009）。
5. **改正**：第三轮之后我曾凭记忆说 Kiessling 2019 论证"守恒律不能唯一确定点电荷受力"，**说反了**。见 4.2。

### 4.1 往前追引用：现代对 runaway 的共识

自 2009 年起的引用数（INSPIRE）：Gralla–Harte–Wald 148，Spohn 2000 72，Levine–Moniz–Sharp 20，ADM 1960 61，Adamo–Kozameh–Newman 165。

**共识：runaway 是点粒子近似的产物，降阶是正当的。**
- Gralla–Harte–Wald 2009【全文】："no runaway behavior occurs in the context of perturbation theory … There is no corresponding justification for using the ALD equation."
- Poisson–Pound–Vega 2011【全文】：点粒子只是有尺寸物体的理想化，在这个近似下用 $f_{\rm ext}/m$ 代替加速度是自洽的。
- Burton–Noble 2014【全文】："an extended charge is not subject to runaways … provided its radius exceeds cτ"。
- **Steane, *Am. J. Phys.* 83, 256 (2015), arXiv:1402.1106【全文，已复核】**："there is no known physical object, having a positive bare mass, whose equation of motion has runaway solutions"。§II B 用牛顿模型加入引力自能 $-Gm^2/2R$，重新得到 ADM 的 $q/m=\sqrt G$，随即指出电子的 $q/m$ 大 $2\times10^{21}$ 倍而舍弃引力。只有能量论证，没有动力学。
- 异议或细化：Ekman, Heinzl & Ilderton 2021（降阶级数发散，需要 Borel–Padé 重求和）【摘要】；Kosyakov 2019（把 runaway 当作物理）【全文】。

**其他现代工作**（均不涉及引力改变判据）：Moxon & Flanagan 2018（Gralla–Harte–Wald 推广到电荷二阶，平直时空）【全文】；Harte, Blanco & Flanagan 2025（自相互作用有尺寸物体满足有效场中的 MPD 方程，形式体系）【摘要】；Casadio 2013（ADM 球壳的静态重推，加最小长度）【摘要】；Verma 等 2026（Einstein–Maxwell 后牛顿带电双星，直接用 LL 降阶）【全文检索】。

**问过"引力能否消除 runaway"的只有 Newman 一系**：Kozameh & Quiroga 2012 声称没有 runaway，但细看是因为把零质量无穷远处的辐射数据当作给定的，不是引力的作用【全文，解读】；Quiroga 2017 断言额外的 Bondi 数据项"would give well behaved solutions"，没有证明【全文】。

**$G$ 一阶或后牛顿的 runaway 判据计算：没有找到。**

### 4.2 数学物理：Kiessling 一系

- **Kiessling, *PRD* 100, 065012 (2019), arXiv:1907.11239【全文】**：粒子动量直接假设为 $p = m_b\gamma v$（式 (12)）。受力由总动量守恒定义，只有场动量局部可积时才有意义，所以排除标准 Maxwell–Lorentz 理论，适用于 Bopp–Landé–Thomas–Podolsky（BLTP）电动力学。结果：受力只依赖 $q, v, a$，"the infamous third-order time derivative … does not show up"，运动方程是加速度的 Volterra 积分方程，局部适定。**这与我此前凭记忆的说法相反。**
- **Kiessling 2012, arXiv:1107.2286【全文】**：才是讨论歧义的那篇。自力的正则化"cannot be expected to be unique"；静态 Born–Infeld 两电荷解守恒一切，所以"simply postulating their conservation does not suffice to produce the right motions"。最后用 Hamilton–Jacobi 定律确定 $P$。
- **裸质量的符号**：McGuigan & Kiessling 2026（*IJMPA* 41, 2650044, arXiv:2506.08799）【全文】：常场问题中，截断到 $O(\kappa^3)$ 时 $m_b<0$ 出现 runaway；$O(\kappa^4)$ 时出现来回振荡、幅度增长的"over-stability"，问题开放。氢原子光谱要求 $\kappa$ 大，因而 $m_b$ 为大的负值。
- **加入 GR**：Kiessling & Tahvildar-Zadeh 2019（*IJMPD* 28, 1930017, arXiv:1906.00308）【全文】：打开引力后，只改电磁定律"is by itself not sufficient to obtain even a merely well-defined joint initial value problem"，引力耦合也要改。引 Geroch–Traschen：没有黑洞的奇点不能有正裸质量。Burtscher–Kiessling–Tahvildar-Zadeh 2021：弱 Bianchi 恒等式对裸 RN 奇点不成立，对 Hoffmann 的负裸质量 Born–Infeld 点电荷时空成立。
- **这一系没有人讨论 Eliezer 型的高阶项**，$p = m_b\gamma v$ 总是直接假设。

### 4.3 弱引力猜想与超极端物质

- **没有任何 2000–2026 年的论文研究自引力下 $\xi>1$ 带电物体的 runaway 或辐射反作用稳定性。**
- 弱引力猜想文献里的"自力"指两个相同粒子之间的静态力，不是辐射反作用。Harlow 等 2023 综述【全文】："the electron satisfies the WGC by more than 20 orders of magnitude"；脚注说电荷的经典半径小于康普顿波长，经典自能应在康普顿半径处截断（与路线图第 0.2 节一致）。全文没有 runaway、radiation reaction、Lorentz–Dirac 等词。
- "自力作为宇宙监督者"（Zimmerman, Vega, Poisson & Haas 2013）【全文】：$q\gg m$ 的粒子落入近极端 RN 黑洞，自力阻止过充；轨迹事先给定，不反馈，不讨论 runaway。Sorce & Wald 2017：二阶计算，Kerr–Newman 不能被过充或过旋【摘要】。
- 带电星体（Anninos & Rothman 2002；Arbañil, Lemos & Zanchin；Bonnor 星）：都是 $Q\le M$ 或恰好极端，静态或径向，没有辐射反作用。
- Ponce de León 2004【摘要】：GR 中的电磁质量电子模型一般需要负质量的核心。
- D'Andrea, Luty & Verhaaren 2019【全文】：平直时空中 runaway 来自负裸质量，尺寸大于经典半径即可消除。

### 4.4 有尺寸电荷与世界线有效场论（猜想 B）

- **Galley, Leibovich & Rothstein, *PRL* 105, 094802 (2010)【全文】**：带电球壳的有限尺寸修正在 $O(R^2)$；极化率算符列出但留待以后。
- **Forgács, Herpay & Kovács, *PRL* 109, 029501 (2012)** 与 Kazinski（arXiv:1103.0537）：指出漏掉了 $O(R)$ 的 $\int a^2\,d\tau$ 项，球壳系数 $-\tfrac29 e^2R$（取自 Nodvik 1964）。
- **Galley, Leibovich & Rothstein 的回复, *PRL* 109, 029502 (2012)【全文，已复核】**："the worldline of the point charge in the effective theory is not an observable. Hence, one may perform worldline shifts by amounts of order R without changing the predictions"；该项可由 $x^\mu \to x^\mu - (C/m)a^\mu$ 消去，并入偶极算符。**评注**：对 $C\int a^2$ 变分得到正比于 $\dot a$ 的动量修正，即 Eliezer/Bhabha 型的 $P_2$ 项；系数 $\propto e^2R$，点极限下为零。原文没有这样说。
- Steane 2015【全文】：LL 与 Eliezer/Ford–O'Connell 方程"only differ at the next higher order in … the width of the charge distribution"。
- Sokolov 2009（*JETP* 109, 207）【摘要】：明确重提狄拉克的"arbitrary choice"，选了另一个与 $u$ 不共线的四动量，但没有当作结构问题。
- Medina 2006（*J. Phys. A* 39, 3801，**不是** *Am. J. Phys.*）【摘要】：需要 $m_0 > m_e/3$，给出半径下限；满足时"never run away"。
- Goldberger & Rothstein 2006【全文】（引力情形）：Ricci 世界线算符"can be removed by a redefinition of the metric tensor"。
- Jakobsen 2024（*PRL* 132, 151601）【全文】：世界线 EFT 中的自旋、磁化、电极化修正。
- 激光等离子体综述（Di Piazza 等 2012；Blackburn 2020；Gonoskov 等 2022）【全文】：都认为 LAD 有病态，模拟几乎只用 LL 方程。
- **结论**：猜想 B 没有被完整写下；"Eliezer 的 $P_{2i}$ 歧义正是世界线重定义的冗余"以及"点极限选出 $P_{2i}=0$"都没有找到。

### 4.5 讨论中澄清的问题（2026-10-23）

本轮汇报后，讨论中澄清了几件事，已写进路线图 v3 第 0 节：

- 经典半径与康普顿波长之比 $=e^2/\hbar c=\alpha$，与质量无关。**任何单个基本粒子（电子、质子、μ 子……）都在经典的病发作之前进入量子区域**，所以本项目不研究基本粒子。
- 电子的自加速被量子力学解决，不是被引力；引力要到 $10^{-34}$ cm 才变强，那里属于量子引力，经典理论回答不了。电子不是黑洞。
- 电荷数 $N\gtrsim12$ 时经典处理自洽。本项目的对象是经典宏观带电物体；判据真正起作用的是接近黑洞那样致密的带电天体。
- $10^{-34}$ cm 以下的问题在量子引力框架里有人讨论（引力作为普适正规化子：DeWitt 1964；Isham–Salam–Strathdee 1971–72；QED 加引力的跑动：Robinson–Wilczek 2006；弱引力猜想）。前三项凭记忆，待核。

### 4.6 仍未完成

- Boulware 1973、Levine–Moniz–Sharp 1977 全文（付费墙）。
- 4.5 节中凭记忆列出的量子引力方向文献。
- Ekman, Heinzl & Ilderton 2021 关于降阶级数发散的全文。
