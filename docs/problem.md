# 带电黑洞在两种辐射反作用下的运动：问题陈述与阅读计划

> 2026-10-30 写定 · 2026-11-01 补进度 · 2026-11-18 按新颖性核查修改 P1–P3
> 这是项目现在的主线。[`roadmap.md`](roadmap.md)（v4）的带电双星后牛顿计算是另一条线，已暂停，作为备选。它是质量相当的 PN 展开，这里是小物体展开，两者只在"极端质量比 + 弱场"一角重叠，那条线不是这个问题的特例。

---

## 一、问题

### 1.1 出发点

经典电磁学对带电物质是不完备的：点电荷的自场在自身位置发散，Lorentz–Dirac 方程要靠质量重整化和一个约定才写得出来；有尺寸的带电体又需要一种电磁学之外的力（Poincaré 应力）才能维持。

加入广义相对论后，有两处根本变化：

1. **运动不再是额外假设。**Einstein 方程本身蕴含能动量守恒，物体怎么运动由场方程决定（Einstein–Infeld–Hoffmann 纲领）。
2. **理论自己提供了带电"粒子"。**Reissner–Nordström 黑洞是真空 Einstein–Maxwell 方程的精确解，不需要任何物质模型。代价是只对 $\xi = q/(\sqrt G\, m) \le 1$ 成立：$f(r)=1-2M/r+Q^2/r^2$ 的根 $r_\pm = M\pm\sqrt{M^2-Q^2}$ 只在 $Q\le M$ 时存在，$\xi>1$ 时没有视界，只剩裸奇点。

因此，**只用 Einstein–Maxwell 场方程、不加任何物质就能描述的带电孤立物体，是 $\xi\le1$ 的黑洞。**本项目问：对这类物体，经典理论能否自洽地给出它的运动。

这类物体还有一个性质：**两种辐射反作用天然同一量级。**小黑洞的电磁自力 $\propto q^2=\xi^2m^2$，引力自力 $\propto m^2$，两者之比 $\xi^2$，$\xi\sim1$ 时完全相当，不需要任何调节（roadmap v4 里要靠偶极–四极交叉区才能做到）。

### 1.2 设定

- **理论**：真空 Einstein–Maxwell 方程，不加任何物质。
- **物体**：一个小的带电黑洞，质量 $\lambda m$、电荷 $\lambda q$，$\xi = q/m \le 1$ 固定（$G=c=1$）。
- **背景**：一个 Einstein–Maxwell 电真空解 $(g_0, F_0)$，例如平直时空加外电磁场，或一个大的 Reissner–Nordström 黑洞。
- **方法框架**：取一族精确解 $(g(\lambda), F(\lambda))$，$\lambda\to0$ 时小黑洞缩成一点（Gralla–Wald 的"普通极限"和"标度极限"）。极限下得到一条世界线 $z(\tau)$，它的运动方程按 $\lambda$ 逐阶修正。

### 1.3 核心问题

> **对带电黑洞（$\xi \le 1$），能否只从 Einstein–Maxwell 场方程出发，严格推出一个适定、没有 runaway 的运动方程，同时包含两种辐射反作用？**

拆成三层：

**P1（一阶，严格推导）**
在上述设定下，世界线零阶满足外场中的洛伦兹力方程（Gralla 2010、2013 已证）。要证明的是一阶偏离满足的方程，它同时包含电磁自力、引力自力，以及两者之间的耦合项（$\propto em$ 的局部项和交叉尾巴项）。
目标方程已知：Zimmerman–Poisson 2014 式 (8.16)（arXiv 版编号），已降阶。要回答：严格推导是否给出同一个方程？黑洞近区（RN）的结构是否带来额外的项？
定位：这是 Gralla 2010 的 Einstein–Maxwell 框架从零阶到一阶的推广，也是 Gralla–Wald 2008（纯引力）和 Gralla–Harte–Wald 2009（平直时空、不含自身引力的电磁）的合并。Linz–Friedman–Wiseman 2014 的结论明说这一步 "has not yet been done"。

**P2（自洽演化的长时间有效性）**
P1 推出的是偏离矢量方程：线性、二阶，本身适定，没有 runaway；降阶方程和未降阶的三阶方程在 $O(\lambda)$ 都与它一致。所以"降阶还是不降阶"在一阶上不是独立的问题，旧版 P2 只是 P1 的推论。
真正的问题在长时间。偏离矢量方程只在固定的有限时间内有效（Gralla–Wald 2008 第 VII 节），旋近要靠自洽方程。Gralla–Harte–Wald 2009 给出自洽方程的三条判据，其中第 (3) 条（自洽解在一阶描述应当成立的时间内与一阶解接近）只写了 "we believe"，没有证明。
问：**对耦合的 Einstein–Maxwell 情形，降阶的自洽方程能否在 $\sim1/\lambda$ 的时间上一致逼近真实解？**纯引力情形这也还没有完全解决。

**P3（近极端）**
Gralla–Wald 2008 中，近区背景稳态、渐近平直不是假设，而是由标度极限（假设 (ii)）和一致性条件（假设 (iii)）推出的（原文："derived here rather than assumed"）。GW 还指出（脚注 7），存在不稳定性时，满足这些假设的一族解是否存在 "far from clear"。$\xi\to1$ 时有几个迹象表明 (ii)、(iii) 可能对 $\xi$ 不一致：表面引力趋于零，近区弛豫（准正则模衰减）时间发散，出现 AdS$_2$ 近视界喉，极端视界有 Aretakis 不稳定性。问：
- $\lambda\to0$ 和 $\xi\to1$ 两个极限能否交换？
- 一阶展开的余项对 $\xi\in[0,1]$ 是否一致有界？如果不一致，失效发生在什么标度（例如 $\lambda$ 与 $1-\xi^2$ 的某个幂相当时）？
- 在极端点，"运动方程"这个概念是否还站得住？

（旧版用 $c\tau_0/r_+\to2/3$ 作为动机，是错的：$\tau_0$ 和 $r_+$ 都正比于 $\lambda$，比值对任何 $\xi$ 都是 $O(1)$，不是展开的小参数。那是把经典电子的结构问题错搬到了黑洞上。）

### 1.4 不在范围内

- $\xi>1$ 的物体（需要物质模型或量子理论）。
- 量子效应。
- 天体物理应用和波形建模（Placidi、Zhang 等组的方向）。

### 1.5 已知与空白

| | 状态 |
|---|---|
| Einstein–Maxwell 初值问题 | 适定（经典定理） |
| 纯引力自力的严格推导 | 已有：Gralla–Wald 2008 |
| 纯电磁自力的严格推导 | 已有：Gralla–Harte–Wald 2009（平直时空，物体是物质，不含自身引力） |
| 一般场论框架下的运动 | Gralla 2010、2013：含 Einstein–Maxwell 和自身引力，允许黑洞，背景弯曲；**只到零阶洛伦兹力**（"I do not treat self-force corrections"） |
| 非线性有效场框架 | Harte–Blanco–Flanagan 2025：固定度规上的标量场；引力只有草图（"left for later work"），没有 Einstein–Maxwell 和黑洞 |
| 电磁 + 引力耦合自力 | 2014 年写出（ZP 式 (8.16)），但**把物体设为点粒子**：ZP 的洛伦兹力方程 "on a formal level" 假定、正规化沿用 Detweiler–Whiting 处方；LFW 明说重整化是 "assume"。不是严格推导 |
| 从无穷远渐近场推运动 | Newman、Kozameh 一系；runaway 仍在，所用模型不适定，"答案未知" |
| 近极端 | 静态平衡解（Majumdar–Papapetrou）、慢速散射（Ferrell–Eardley）已知；辐射反作用和有效方程的适定性没有人做 |

"耦合情形没有人严格做过"：2026-11-18 对照 ZP、LFW、Gralla 2010、Gralla 2013、Harte–Blanco–Flanagan 2025 的原文确认。**2026-11-24 补查**：INSPIRE 上检索引用 GW 2008、GHW 2009 的论文，发现漏查的 Futamase–Hogan–Itoh 2008（PRD 78, 104014, arXiv:0811.4020）。他们从 Einstein–Maxwell 真空方程推出了外场中小带电黑洞的运动方程，但依赖对展开形式的具体假设，尾巴项含未定函数，看不到引力自力项，背景场为零时尾巴消失，也没有与 ZP 对照。它是最接近的前人工作，已写入论文引言；"空缺"的表述随之收窄为：没有人在一般假设下、连同引力自力一起推出，并与 ZP 对照。LFW 结论原话："What has not yet been done is to extend to the coupled fields we consider here the underlying justification based on variants of matched asymptotic expansion"。详见第四节。

---

## 二、阅读计划

按顺序读。每篇列出读什么、带着什么问题读。

### 第一阶段：方法的核心（P1、P2 的基础）

**1. Gralla & Wald 2008**，*A rigorous derivation of gravitational self-force* — https://arxiv.org/abs/0806.3293
全文精读，这是整个项目的方法模板。
- 一族解要满足的两个假设（普通极限、标度极限）是怎么表述的？为什么需要两个？
- 世界线是怎么从一族精确解中定义出来的？
- 推出 MiSaTaQuWa 方程用到了哪些步骤？哪一步依赖于"真空"？
- 他们怎样看待"自洽"演化和降阶？

**2. Gralla, Harte & Wald 2009**，*Rigorous derivation of electromagnetic self-force* — https://arxiv.org/abs/0905.2391
全文精读，与第 1 篇对照。
- 物体自身的引力场是否被包含？背景是平直还是弯曲？
- 为什么降阶后的方程有依据、未降阶的没有？这正是 P2 要推广的论证。
- 和第 1 篇相比，电磁情形多了哪些技术困难（电荷守恒、偶极矩、自旋）？

**2b. Gralla 2010**，*Motion of small bodies in classical field theory* — https://arxiv.org/abs/1002.5045
P1 的零阶起点：Einstein–Maxwell（含自身引力、允许黑洞）的三条假设，以及零阶洛伦兹力的推导。读假设的表述和 Einstein–Maxwell 一节。配套：Gralla 2013（https://arxiv.org/abs/1303.0269 ）给出 ADM 质量和 Gauss 电荷的面积分公式。

**3. Poisson, Pound & Vega 2011**，*The motion of point particles in curved spacetime*（Living Review）— https://arxiv.org/abs/1102.0529
工具书，按需查：Hadamard 展开、Detweiler–Whiting 分解、Fermi 坐标、匹配渐近展开。

### 第二阶段：已有的耦合结果（P1 要对照的对象）

**4. Zimmerman & Poisson 2014**，*Gravitational self-force in nonvacuum spacetimes* — https://arxiv.org/abs/1406.5111
- 物体模型是什么？哪些是假设，哪些是推导？
- 最终运动方程（arXiv 版式 (8.16)；旧记录写作 (119)，可能是期刊版编号，待核）里耦合项长什么样？是否降阶？

**5. Linz, Friedman & Wiseman 2014**，*Combined gravitational and electromagnetic self-force on charged particles in electrovac spacetimes* — https://arxiv.org/abs/1406.5112 （配套：https://arxiv.org/abs/1404.7039 ）
- 他们明确"假设"了什么？
- "重整化中电磁场和引力场不混合"是怎么得到的？严格推导中是否仍然成立？

**6. Pound 2012**，二阶引力自力（选读，P1 完成后）— https://arxiv.org/abs/1201.5089

### 第三阶段："答案未知"的来源（短读）

**7. Kozameh, Newman, Ortega & Silva-Ortigoza 2009** — https://arxiv.org/abs/0803.3921
只读引言和结论。开放问题原话怎么写？"位置"怎么定义，和 Gralla–Wald 有什么不同？

**8. Lun & Chow 1994** — https://arxiv.org/abs/gr-qc/9409024
ℓ=1 增长率 $3M/2Q^2$ 是怎么出来的？视界条件是假设还是推论？

**9. Kozameh, Kreiss & Reula 2008** — https://arxiv.org/abs/0708.1933
只读摘要和主要定理：Robinson–Trautman Einstein–Maxwell 为什么不适定。

### 第四阶段：近极端（P3）

**10. Giorgi 2019**，RN 线性稳定性（全次极端范围）— https://arxiv.org/abs/1910.05630
只读引言。"稳定"指什么？引力和电磁扰动在 RN 上怎样耦合？

**11. 极端 RN 的视界不稳定性**
- Aretakis 2011，极端 RN 标量扰动 I — https://arxiv.org/abs/1110.2007
- 极端 RN 引力与电磁扰动的不稳定性 — https://arxiv.org/abs/2211.09182
- 极端 RN 的非线性稳定性（2024，选读）— https://arxiv.org/abs/2410.16234 。如果极端黑洞在非线性意义下稳定，P3 的出发点要调整。
- 问题：这种不稳定性会不会影响小黑洞在极端背景附近、或小黑洞本身接近极端时的运动？

**12. Ferrell & Eardley 1987**，*Slow-motion scattering and coalescence of maximally charged black holes* — https://doi.org/10.1103/PhysRevLett.59.1617
辐射在这里是怎么处理的？

**13. Zimmerman, Vega, Poisson & Haas 2013**，*Self-force as a cosmic censor* — https://arxiv.org/abs/1211.3889
自力怎样阻止黑洞被充到超极端？

### 背景（选读）

- Einstein, Infeld & Hoffmann 1938 — https://edition-open-sources.org/sources/10/17/index.html
- Geroch & Traschen 1987，为什么广义相对论里不能放点粒子 — https://doi.org/10.1103/PhysRevD.36.1017

### 之后的综述（2026-11-01 补充）

| 综述 | 链接 | 用处 |
|---|---|---|
| Harte 2015，*Motion in classical field theories and the foundations of the self-force problem* | https://arxiv.org/abs/1405.5077 | **下一篇要读。**有尺寸物体的运动，标量、电磁、引力用同一套框架。决定 P1 走 Gralla–Wald 还是 Harte 的路线 |
| Harte, Blanco & Flanagan 2025 | https://arxiv.org/abs/2504.11912 | 看"生成泛函"是什么，以及对引力情形说了什么 |
| Pound 2015，自力入门讲义 | https://arxiv.org/abs/1506.06245 | Gralla–Wald 和 Pound 两种严格方法的区别 |
| Barack & Pound 2019 | https://arxiv.org/abs/1805.10385 | 2011 年之后领域的全景 |
| Pound & Wardell 2021 | https://arxiv.org/abs/2101.04592 | 最全，含二阶自力，当工具书用 |

Ferrell–Eardley 和 Geroch–Traschen 只有期刊版本，需要机构访问。

---

## 三、检查点

- [ ] **读完第 1、2 篇**：用自己的话写出"一族解"的两个假设，以及它们在 Einstein–Maxwell 里需要怎样改。这一步决定 P1 能不能做。
- [ ] **读完第 4、5 篇**：列出 2014 年两篇里**所有的假设**，逐条判断严格推导时是否需要、能否去掉。这张清单就是 P1 的工作清单。
- [ ] **读完第 10–13 篇**：P3 是一个物理问题，还是一个数学上已经有答案的问题？

---

## 四、进度（截至 2026-11-18）

**已读过一遍：**第 1–5 篇，第 7 篇。检查点都还没写。

几条结论：

- **第 3 篇**：第 18 节（电荷）是固定背景上的 DeWitt–Brehme–Hobbs 方程，明确做了降阶；第 19 节（点质量）靠放宽 Lorenz 规范，正是 Gralla–Wald 批评的那一步；第五部分（Pound）只处理引力。耦合情形的空白在这篇里同样存在。
- **第 4 篇**：式 (119)（arXiv 版 (8.16)）就是 P1 **要严格证明的目标**。它已降阶；G→0 时回到 Landau–Lifshitz 加电磁尾巴。
- **第 5 篇**：重整化这一步电磁和引力不混合，耦合只出现在正则场里。但重整化是"假设"成立的，不是推出来的。
- **第 7 篇**：得到的是未降阶的三阶方程。其中的 $\ddot v$ 项可以用降阶换掉，但这只是一阶上的推论（见 P2），不单独算一个问题。
- **第 8、9 篇**：结论已经知道（runaway 还在，模型不适定），不用读。
- **第 10–13 篇**：等做 P3 时再读。

**2026-11-18 新颖性核查**（对照原文，引文已逐条核对）：

- **ZP 2014**：点粒子 + δ 源，Detweiler–Whiting 正规化，Lorenz 规范直接施加；给出完整降阶方程 (8.16)；不声称严格，没有把物体设为小黑洞。
- **LFW 2014**：点粒子；重整化在摘要中写明是 "assume"；没有给出带尾巴项的显式方程；结论明说耦合场的严格论证 "has not yet been done"，并说有信心 "no surprise will emerge"。
- **Gralla 2010、2013**：覆盖 Einstein–Maxwell、自身引力、黑洞、弯曲背景，但只到零阶洛伦兹力；近区稳态同样是推论。
- **Harte–Blanco–Flanagan 2025**：不覆盖。
- **GW 2008**：近区稳态是推论（第 III 节），不是假设，见 P3。**GHW 2009**：平直时空，物质物体，不含自身引力；摘要说降阶方程有依据、未降阶的没有，但依据是三条判据，第 (3) 条只是 "we believe"，见 P2。

**结论**：P1 的空白是真的，仍是主体。但目标方程已知，P1 的新意在于推导本身和黑洞近区可能带来的项；物理上的意外更可能来自 P2、P3。

**下一步：**先写第一个检查点（以 Gralla 2010 的 Einstein–Maxwell 假设为起点），再读 Harte 2015。
