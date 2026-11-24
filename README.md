# classical_particles

三篇关于经典电荷运动的旧稿（2006、2010、2011），二十年后对它们的回顾，以及由此展开的研究。

## 当前主线（2026-11 起）

**带电小黑洞在电磁与引力辐射反作用同时存在时的一阶运动**：用 Gralla–Wald 的"一族精确解"方法，从 Einstein–Maxwell 场方程推导运动方程，与 Zimmerman–Poisson 2014 的方程对照。

| 文件 | 说明 |
|---|---|
| [`paper/paper.pdf`](paper/paper.pdf) | 论文草稿（中文，xelatex）。主结果（一阶运动方程）尚未证明，首页列出全部待办 |
| [`derivations/p1/calc.pdf`](derivations/p1/calc.pdf) | 验证计算工作文档：分四级逐步推进，记录每一步的中间式子和检验 |
| [`derivations/p1/`](derivations/p1/) | 推导草稿（C0、C1）和计算代码 |
| [`docs/problem.md`](docs/problem.md) | 问题陈述与新颖性核查 |
| [`CHANGELOG.md`](CHANGELOG.md) | 每一步做了什么、得到什么、还欠什么 |

进度：验证计算第 1 级（复现 Gralla–Wald 2008）、第 2 级（物体带电，背景不带电：运动方程恰好多出尾巴电场的洛伦兹力）已通过；第 3 级（物体不带电，背景带电）进行中。

## 早期内容

## 内容

### 原始文稿

| 文件 | 说明 |
|---|---|
| [`original_papers/2006paper.pdf`](original_papers/2006paper.pdf) | 《论电荷和质点的运动》，南京大学天文学系本科毕业论文，2006 年 6 月，79 页 |
| [`original_papers/2010paper.pdf`](original_papers/2010paper.pdf) | *How Much Energy Electrons Can Emit? Radiation Theory Based on the Lorentz-Dirac Equation and other Modified Equations*，Ohio State，18 页，**未完成** |
| [`original_papers/2011paper.pdf`](original_papers/2011paper.pdf) | *Radiation Reaction in Ultra-High Magnetic Fields and Highest Energy Cosmic Rays with Their Neutrinos*，4 页，**未完成** |

2010 那篇用两个算例判别 LD 及其修正方程：一维谐振束缚电子的**辐射谱轮廓**因方程而异；强磁场中超相对论电子的同步辐射在 `Bγ > 10¹⁶ G` 时功率显著偏离经典的 `P ∝ v²γ²`。

2011 那篇处理超高磁场中的辐射反作用与超高能宇宙线。值得注意的是它**已经改用降阶后的 Landau–Lifshitz 形式**（把 `ȧ` 换成 `Ḟ_ext/m`），而不再是 2006 年那条修改 LD 的路线；结论是 `bω₀ > 1` 时辐射功率与 `γ` 无关，趋于常数，与同步辐射的 `γ²` 和曲率辐射的 `γ⁴` 都不同。

### 回顾与重写

| 文件 | 说明 |
|---|---|
| [`REVIEW_2026.md`](REVIEW_2026.md) | 2026 年的回顾：原文梳理 + 领域二十年进展 + 可能的后续方向 |
| [`notes/notes.pdf`](notes/notes.pdf) | 正在重写的笔记，AASTeX preprint 格式。把 2006 年那篇的经典电子部分按一条线索重排，并与它当年所依据的文献对照 |

## 原文讲了什么

经典点电荷的运动由 Lorentz–Dirac (LD) 方程描述，有质量粒子由 MiSaTa(QuWa) 方程描述。两者都有问题：LD 的解存在自加速（runaway）与预加速，且未计入电磁场能动量张量引起的时空弯曲。

论文试图**修改方程**来消除这些病态，走了四步：

1. **分离奇异项法** —— 证明「能动量守恒」与「用 ½(A_ret − A_adv) 算 Lorentz 力」的等价定理，然后把奇异项放宽成双参数族 `A = λA_ret + κA_adv`。得到的 `κ > 0` 分支表现为「自减速」，论文称之为 **"Aristotle 运动"**。
2. **Eliezer 方程系列** —— 把 Dirac 待定四矢量 `B^μ` 推广成含任意阶导数的级数，得到一族更高阶的运动方程，数值扫描参数空间。
3. **广义相对论修正** —— 提出「有效度规」假设，用正则部分的势构造能动量张量，迭代求解 Einstein 方程。
4. **实验判别方案** —— 设计用定角辐射功率 `dP/dΩ` 的时间演化来拟合模型参数。

## 二十年后

`REVIEW_2026.md` 的三个主要结论：

- **领域走了相反的路**。主流答案不是修改 LD 方程，而是**降阶**（Landau–Lifshitz / Ford–O'Connell）：LD 只是有效方程，其物理解全部落在慢流形上。论文第 3 节提高导数阶数的做法，在原理上会加剧问题而非解决它。

- **实验从「不可能」变成「已完成」**。论文判定判别实验不可行，依据的特征时标算错了 10 个量级（漏了 1/4πε₀；真值 ~10⁻²³ s 而非 10⁻³³ s）。2018 年起的激光—电子束对撞实验正是在做这件事，2026 年已有 5σ 结果。

- **强磁场方向的直觉数值上是对的，但被 QED 抢先了 137 倍**。经典辐射反作用变成非微扰的临界场 `B_cl = m²c⁴/e³ = B_Schwinger/α = 6×10¹⁵ G`，恰好落在磁星量级；但 Schwinger 场比它低 1/α，Landau 量子化与强场 QED 早已接管。而论文里那个 **"Aristotle 运动"**，正是今天脉冲星磁层建模中的 **Aristotelian electrodynamics**。

## 后续方向

`REVIEW_2026.md` 第三部分列了 7 个方向，按可独立完成度排序。最推荐的两个：

- **(B, γ) 参数相图** —— 把 τ₀ω_B = 1、B = B_S、χ = 1、γ_RRL 等所有边界画在一张图上，标出脉冲星、磁星、FRB 源与实验室激光，一次性定位「当年那个问题在参数空间的哪一块」。
- **Aristotelian electrodynamics 的量子推广** —— χ ≫ 1 时那个「场 → 速度」的代数关系该是什么形式。这是唯一一条能把当年的 κ 族直接推到当代前沿的路。

## 说明

`REVIEW_2026.md` 中标注的两处原文问题：单位错（已复算确认）与附录 A.3 的 T[F_S] 积分（与 Detweiler–Whiting / Harte–Pound 定理冲突，**待独立核对**）。文中所有数值均可用附录中的脚本复现。
