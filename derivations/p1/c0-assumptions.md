# C0（检查点 1）：一族解的假设——Einstein–Maxwell，带电黑洞

> 草稿 · 2026-11-19
> 来源：Gralla 2010 第 II、IV 节【全文】；Gralla–Wald 2008（下称 GW）第 II–III 节【全文】；Gralla–Harte–Wald 2009（下称 GHW）第 II 节【全文】。
> 原文在 `references/Self-force rigorous/`。

这一份回答 problem.md 第三节的第一个检查点：用自己的话写出"一族解"的假设，以及它们在 Einstein–Maxwell 里需要怎样改。

---

## 0. 约定

- 号差 $(-,+,+,+)$，$G=c=1$，Gauss 单位。
- $F=dA$，$4\pi T_{ab}=F_{ac}F_b{}^c-\tfrac14 g_{ab}F_{cd}F^{cd}$。
- 场方程记为
  $$E_{ab}\equiv G_{ab}-8\pi T_{ab}=0,\qquad \mathcal M^a\equiv\nabla_bF^{ab}=0 .$$
  有源时 $\mathcal M^a=4\pi j^a$。检验：Coulomb 场 $A=-(q/r)\,dt$ 给出 $\mathcal M^t=4\pi q\,\delta^3$。洛伦兹力写作 $m\,u^b\nabla_bu^a=q\,F^a{}_bu^b$，同号电荷相斥。
- 小黑洞的质量 $\lambda M$、电荷 $\lambda q$，$\xi=q/M$ 固定。

## 1. 标度

Einstein–Maxwell 方程在 $(g,A)\to(c^2g,\,cA)$（$c$ 为常数）下不变：
- $G_{ab}$ 在常数共形变换下不变；
- $T_{ab}$（指标全在下）中 $F_{ac}F_{bd}g^{cd}$ 和 $g_{ab}F^2$ 都按 $c^0$ 变；
- $\nabla_bF^{ab}$ 整体乘以 $c^{-1}$，所以方程不变。

这正是 Gralla 2010 选 $\bar A=\lambda^{-1}A$ 的理由（原文："This is also the scaling of $A_a$ that leaves the Einstein-Maxwell theory invariant"）。

对每个 $t_0$，取 $\bar t=(t-t_0)/\lambda$、$\bar x^i=x^i/\lambda$，定义
$$\bar g_{ab}=\lambda^{-2}g_{ab},\qquad \bar A_a=\lambda^{-1}A_a .$$
在带横线的坐标里，分量满足（Gralla 2010 式 (2) 及其矢量类比）
$$\bar g_{\bar\mu\bar\nu}(\lambda;t_0;\bar t,\bar x^i)=g_{\mu\nu}(\lambda;\,t_0+\lambda\bar t,\,\lambda\bar x^i),\qquad
\bar A_{\bar\mu}(\lambda;t_0;\bar t,\bar x^i)=A_{\mu}(\lambda;\,t_0+\lambda\bar t,\,\lambda\bar x^i).$$
也就是直接把 $t=t_0+\lambda\bar t$、$x^i=\lambda\bar x^i$ 代进去。

## 2. 三条假设

照 Gralla 2010 第 IV 节，把 GW 的度规换成一对 $(g_{\mu\nu},A_\mu)$。要求对每个 $\lambda\ge0$ 存在坐标 $(t,x^i)$ 和一个电磁规范，使下面三条成立。记 $r=|x|$，$\alpha=r$，$\beta=\lambda/r$。

- **(i) 普通极限。**在 $r\ge\lambda\bar R$ 上，$(g_{\mu\nu},A_\mu)$ 满足 $E_{ab}=0$、$\mathcal M^a=0$，并且对 $(t,x^i,\lambda)$ 光滑。$\lambda=0$ 时，$(g,A)$ 在包括 $r=0$ 在内处处光滑，$\gamma$（$r=0$）是类时曲线。
- **(ii) 标度极限。**在 $\bar r\ge\bar R$ 上，$(\bar g_{\bar\mu\bar\nu},\bar A_{\bar\mu})$ 对 $(\lambda,t_0,\bar t,\bar x^i)$ 光滑。
- **(iii) 一致性。**在固定的 $(t,\theta,\phi)$ 下，$(g_{\mu\nu},A_\mu)$ 作为 $(\alpha,\beta)$ 的函数在 $(0,0)$ 处联合光滑。GW 后来补注，实际需要的是对全部变量联合光滑。

对 $r<\lambda\bar R$ **不做任何假设**：那里可以是黑洞，也可以是物质（GW："could equally well be 'filled in' with ordinary matter or a black hole"）。

### 和已有文献相比改了什么

| | GW 2008 | GHW 2009 | 这里 |
|---|---|---|---|
| 场 | 只有 $g$，真空 | 只有 $F$，平直时空 | $(g,A)$，耦合 |
| 物体 | 任意（含黑洞） | 物质：$J^a$、$T_{ab}$ 紧支 | 任意（含黑洞） |
| 对场的假设 | 三条极限假设 | $F=F_{\rm ret}+F_{\rm ext}$，外场随 $\lambda$ 光滑 | 三条极限假设（Gralla 2010 的写法） |
| 物体自身引力 | 有 | 没有 | 有 |

对 GHW 只需要改一处：不再假设推迟解加外场的分解，这个分解要在一阶远区里推出来（见 C1 第 6 节）。

### 需要自己补的地方

1. **电磁规范。**$A_\mu$ 的分量依赖规范，所以假设应读作"存在一个规范，使 (i)–(iii) 成立"。允许的规范变换 $A\to A+d\psi$ 中，$\psi(\lambda,x)$ 应满足和 GW 式 (21) 同类的光滑条件。Gralla 2010 没有讨论这一点【已查全文，没有找到】，要自己写清楚。直接用 $F_{\mu\nu}$ 表述可以避开规范问题，但标度极限里 $\bar A$ 的常数部分 $b^{[A]}(t_0)$ 只是纯规范，用 $F$ 表述反而更干净。**倾向：(i)–(iii) 改用 $F_{\mu\nu}$ 表述**，待定。
2. **背景外场。**$\lambda=0$ 的背景 $(g_0,A_0)$ 是电真空解，$F_0\neq0$。它在 $\gamma$ 上光滑，在近区里怎样出现，见下一节。

## 3. 由假设推出的结果（照搬 Gralla 2010 式 (3)–(8)、(45)–(47)）

由 (iii)，在 $\beta$ 上展开到一阶：
$$g_{\mu\nu}=b_{\mu\nu}(t)+O(r)+\lambda\Big[\frac{c_{\mu\nu}(t,\theta,\phi)}{r}+O(1)\Big]+O(\lambda^2),$$
$A_\mu$ 同理，系数记为 $b^{[A]}_\mu$、$c^{[A]}_\mu$。由此得到：

1. **远区背景** $(g_0,A_0)$：一个在 $\gamma$ 上光滑的电真空解。
2. **一阶远区扰动** $(h,a)\equiv\partial_\lambda(g,A)|_{\lambda=0}$：在 $\gamma$ 附近最多像 $1/r$ 那样奇异。
3. **近区背景**
   $$\bar g^{(0)}_{\bar\mu\bar\nu}=b_{\mu\nu}(t_0)+\frac{c_{\mu\nu}(t_0,\theta,\phi)}{\bar r}+O(\bar r^{-2}),\qquad
   \bar A^{(0)}_{\bar\mu}=b^{[A]}_\mu(t_0)+\frac{c^{[A]}_\mu(t_0,\theta,\phi)}{\bar r}+O(\bar r^{-2}).$$
   它不依赖 $\bar t$（稳态），在 $\bar r\to\infty$ 时趋于常数（渐近平直），并且满足 Einstein–Maxwell 方程。**稳态是推出来的，不是假设**（Gralla 2010："Stationarity follows from smoothness of $g_{\mu\nu}(\lambda)$ in $t$"；GW："derived here rather than assumed"）。常数 $b^{[A]}(t_0)$ 是纯规范。

### 外场在近区里出现在一阶

$\bar F_{\bar\mu\bar\nu}=\lambda F_{\mu\nu}$，所以背景外场 $F_0(\gamma)$ 在近区里是 $O(\lambda)$ 的均匀场。背景外场的梯度、背景曲率在近区里是 $O(\lambda^2)$。这和 GW 里潮汐场出现在近区二阶是一样的结构。

物理上的推论：外场在黑洞上感应出的偶极，大小 $\sim(\text{极化率}\sim(\lambda M)^3)\times F_0$。它产生的力比洛伦兹力高 $\lambda^2$ 阶，不进入 $O(\lambda)$ 运动。

## 4. 黑洞的附加假设 (H)

GW 框架本身不要求物体是黑洞。如果要说"这是一个带电黑洞"，需要额外加一条：

- **(H)** 近区背景 $(\bar g^{(0)},\bar A^{(0)})$ 是一个 Kerr–Newman 解 $(M,q,a)$，并且严格次极端：$a^2+q^2<M^2$。极端情形留给 P3。

**依据。**如果近区背景含有正则的非简并视界，稳态电真空黑洞的唯一性定理会给出 Kerr–Newman【待核：具体定理（Mazur 1982、Bunting 1983 一系）及其条件，例如解析性、视界连通性】。(H) 要作为假设写出来，不能说是推论。

**(H) 在 $O(\lambda)$ 运动里固定了什么。**
- 自旋 $S=\lambda^2 Ma$ 和偶极矩都是 $O(\lambda^2)$，它们产生的力相对洛伦兹力是 $O(\lambda)$，所以**会进入一阶**。GW 式 (111) 的第一项就是 Papapetrou 自旋力；GHW 式 (134)–(136) 有电磁偶极 $Q_{ab}$ 和自旋。
- 对一般物体，这些多极矩是自由参数。对 Kerr–Newman，它们由 $(M,q,a)$ 固定：相对质心的电偶极为零，磁偶极 $\mu=qa$，即旋磁比 $g=2$。
- 这是"黑洞结构带来什么"的第一个具体回答：在偶极这一层，黑洞就是 $Q_{ab}$ 由 $g=2$ 固定的 GHW 物体。
- 和 ZP（无自旋点粒子）对照时取 $a=0$，即 Reissner–Nordström，所有偶极项都消失。

## 5. 为 P3 留的记号

- (i)–(iii) 是对**固定的** $\xi$ 提的。常数 $\bar R$，以及 (ii)、(iii) 中"光滑"背后隐含的导数界，都可能依赖 $\xi$。
- P3 问的就是：当 $\xi\to1$ 时，这些量能否对 $\xi$ 一致选取。
- GW 允许一阶近区扰动随 $\bar t$ 线性增长、二阶随 $\bar t$ 二次增长（GW 第 VI 节）。近极端时近区的弛豫时间 $\sim1/\kappa$ 发散，这里可能是失一致的地方。这一点到 P3 再展开。

## 6. 待核清单

- [ ] Kerr–Newman 唯一性定理的确切表述和条件（第 4 节）。
- [ ] 电磁规范自由度的光滑条件，以及是否改用 $F_{\mu\nu}$ 表述（第 2 节）。
- [ ] Gralla 2010 附录 A 的分布论证（$1/r$ 奇异的场经二阶线性算子作用后只给出 $\delta^3$，不出现 $\delta$ 的导数）在耦合方程组上是否照搬成立。C1 用到这一点。
