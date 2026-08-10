# 《论电荷和质点的运动》二十年回顾与后续方向

> **原文**：张栋，《论电荷和质点的运动》，南京大学天文学系本科毕业论文，指导教师罗新炼副教授，2006 年 6 月（79 页，`original_paper.pdf`）
> **本文档**：2026 年 8 月整理。内容为原文梳理 + 领域二十年进展回顾 + 可能的后续方向。
> **原始研究动机（论文正文未明写）**：考察强磁场（中子星／磁星 10¹⁵–10¹⁶ G）下经典电荷运动是否会给出新结果，以及强场下的加速机制是否与其他情形不同。

---

## 目录

- [第一部分　二十年前这篇文章做了什么](#第一部分二十年前这篇文章做了什么)
- [第二部分　二十年后的回顾、展望、什么时候用经典](#第二部分二十年后的回顾展望什么时候用经典)
- [第三部分　接下来可以做什么](#第三部分接下来可以做什么)
- [附录　数值复算](#附录数值复算)
- [参考文献](#参考文献)

---

# 第一部分　二十年前这篇文章做了什么

## 1.1 出发点

经典电动力学的三大支柱（守恒律、Maxwell 方程、Lorentz 力）并不包含电荷自身的运动方程。点粒子模型下自场在电荷处发散，Lorentz 力形式失效。历史上三条路线：

- **Poincaré (1906)**：广延观 + 额外的束缚能动量张量（Poincaré 应力）——为解决单一矛盾而假设一种力，应者寥寥；
- **Lorentz (1909)**：均匀带电球模型，得到 Lorentz 方程 $m\dot{v} - \tfrac{2}{3}e^2\ddot{v} + \cdots = F$——问题是 4/3 因子、非协变、以及自加速崩溃；
- **Dirac (1938)**：粒子观 + 能动量守恒（管形积分），得到 Lorentz–Dirac (LD) 方程

$$ma^\alpha = F^\alpha_{\rm ext} + \tfrac{2}{3}e^2\left(\delta^\alpha{}_\beta + u^\alpha u_\beta\right)\dot a^\beta$$

论文认定 LD 方程有三个问题：(i) Dirac 待定四矢量 $B^\mu$ 的选择只受 $v\cdot B=0$ 约束，形式有任意性；(ii) 两种推导都只算了管面积分、没算盖面（cap）积分；(iii) **未消除自加速崩溃**——Dirac 本人把初条件放到 $\tau\to\pm\infty$ 来回避，代价是"未卜先知"的预加速，违反因果律。

论文另加了第四条自己的观察：电子自场的能动量张量必然弯曲时空，而世界线上 $T_{\rm em}\sim O(r^{-4})$ 发散，Einstein 方程在那里失效。

论文在引言中给出三条"合乎物理图像"的判据：(1) 可由常规微分方程 + 初值完全解出；(2) 稳定态受小扰动后回归稳定，而非无限逼近光速或在 $\pm c$ 之间振荡；(3) 符合能量守恒与因果律。

## 1.2 第一部分：经典电子运动论

### §1 意义与历史背景
综述上述三条路线，给出 LD 方程的一维化简。用双曲参数 $\dot t=\cosh q,\ \dot x=\sinh q$（$q$ 为快度）化 LD 为

$$\dot q - \tfrac{2}{3}\ddot q = F^{01}$$

自由电子解 $q=C_1e^{3\tau/2}+C_2$：$C_1\neq0$ 即 runaway，$C_1=0$ 则需人为限定初条件。

### §2 分离奇异项方法（本文第一个原创点）

**等价定理**：由能动量守恒得到的运动，等价于以 $\tfrac12(A_{\rm ret}-A_{\rm adv})$ 作矢势按 Lorentz 力公式算出的运动。论文给出了一个比张宗燧／Poisson 更直接的初等证明（推迟势 Taylor 展开逐项相消）。

**关键一步**：论文计算 $A_S=\tfrac12(A_{\rm ret}+A_{\rm adv})$ 产生的能动量张量对电子的影响，得到

$$\frac{dP^\alpha}{d\tau}=\frac{a^\alpha}{2r}+\frac{2}{3}\dot a^\alpha \qquad (1.2.17)$$

并由第二项非零断言：**"这表示 $\tfrac12(A_{\rm ret}+A_{\rm adv})$ 对于电子的运动，还是有影响的"**。既然标准的奇异项选择并非"无影响"，论文就把它放宽成双参数族

$$A^\alpha = \lambda A^\alpha_{\rm ret} + \kappa A^\alpha_{\rm adv}$$

分两支讨论：

**(a) $\lambda+\kappa=0$**（有效部分满足齐次 Maxwell 方程）。一维自由电子的 $a$–$v$ 相图有解析解

$$a=\left[-\frac{3}{8\kappa}\ln\frac{1+v}{1-v}+C\right](1-v^2)^{3/2}$$

$\kappa=-1/2$ 即回到 LD（此时系数为 $+3/4$）。最终方程为

$$ma^\alpha = F^\alpha_{\rm ext} - \frac{4\kappa}{3}e^2\left(\dot a^\alpha - a^2u^\alpha\right) \qquad (1.5.1)$$

**这实质上就是 LD 方程把辐射反作用系数 $2e^2/3$ 换成自由参数 $-4\kappa e^2/3$。**
- $\kappa<0$：凸函数，仍然自加速；
- $\kappa>0$：凹函数，与 $v$ 轴必有交点，表现为**自减速**，电子受扰后回到某个匀速状态。论文称之为 **"Aristotle 运动"**（因为像 Aristotle 描述的运动图像：力决定速度而非加速度）。

**(b) $1-\lambda=-\kappa$**（去掉的奇异项中 ret 与 adv 地位相同）。以 $\kappa$ 为参量数值求解，发现 $\kappa=-1$ 是分水岭：$\kappa>-1$ 仍自加速；$\kappa<-1$ 出现**临界速度**——初速低于 $v_c$ 则永远低于 $v_c$，高于 $v_c$ 则稳定自减速。数值：

| $\kappa$ | $-1.5$ | $-2$ | $-3$ | $-500$ | $-5000$ |
|---|---|---|---|---|---|
| $v_c$ | 0.3473 | 0.2261 | 0.1341 | $6.67\times10^{-3}$ | $6.67\times10^{-4}$ |

论文注意到这支模型推广到二维以上会遇到 $\mathbf{n}\cdot\mathbf{v}$ 在电子自身位置不确定的问题，建议试 $\mathbf{n}\cdot\mathbf{v}\to\tfrac13 v$ 的平均，但没有做。

### §3 Eliezer 方程系列（第二个原创点）

Dirac 的 $B_\mu$ 只需满足 $v\cdot B=0$。论文把它推广成

$$B_\mu = P_0v_\mu + P_1\dot v_\mu + P_2\ddot v_\mu + \cdots + P_iv^{(i)}_\mu+\cdots$$

由 $v^\mu v_\mu=-1$ 逐阶微分得到一串恒等式（系数呈杨辉三角形式，附录 A.4 给了通项），推出一般化的运动方程 (1.3.3)。特例：

- $P_i=0$ → LD 方程；
- 只有 $P_2$ → 论文称的"Eliezer 方程"，一维化简为三阶非线性 ODE

  $$\dddot q - B\ddot q + A\dot q - \tfrac12\dot q^3 = 0,\qquad A=\frac{3m}{2k},\ B=\frac{e^2}{k}$$

  数值分析结论：**无论 $(A,B)$ 同正还是同负，电子总是自行加到光速**；$A>0,B>0$ 且 $A\gg B$ 时出现更糟的行为——速度在 $+c$ 与 $-c$ 之间大幅高频振荡（Fig 1-3-1 ~ 1-3-4）。
- $P_2,P_4$ 非零 → 五阶方程 (1.3.11)。扫描了 22 组 $(A,B,C)$，结论表里 (\*) 表示 $\pm c$ 振荡、(\*\*) 表示单调加到光速、"待定"表示既没到光速也不回归匀速、而是**在有限速度区间内周期性脉动**（"很像一个振子，只是振子本身一直向前运动"）。

论文自己承认："自然界难道真会去选择更为复杂的 $B_\mu$ 吗？倘若真如此，我们真是无法想象自然界作这种选择的根源。"

### §4 考虑广义相对论效应（第三个原创点）

指出世界线上 $T^{\alpha\beta}_{\rm em}=T_{\rm rad}+T_{\rm bnd}$，其中 $T_{\rm rad}\sim O(r^{-2})$、$T_{\rm bnd}\sim O(r^{-4})$，故 $G^{\alpha\beta}\to\infty$。

**处理方案（论文称为"有效度规"假设）**：类比 §2 从 $A_{\rm ret}$ 中分离奇异项保留 $A_{(R)}$，这里也从 $T_{\rm em}$ 中分离出"有效"部分 $T_{{\rm em}(g)}$——具体做法是**用有效势 $A_{(R)}$ 去构造能动量张量**，理由是电子等价地处在一个 $A_{(R)}$ 的外场空间中。

然后联立三个方程迭代：

$$\Box A^\alpha_{(R)}-R^\alpha{}_\beta A^\beta_{(R)}=-4\pi j^\alpha,\qquad G^{\alpha\beta}=\kappa T^{\alpha\beta}_{{\rm em}(g)},\qquad f^\alpha=eF^{(R)\alpha}{}_\beta u^\beta+F^\alpha_{\rm ext}$$

从 Minkowski 度规出发做一级近似 → 算 $T_{{\rm em}(g)}$ → 得 $G^{\alpha\beta}$ → 修正度规 → 新的运动方程 → 二级近似 → 迭代。一维化简下算出

$$R^\alpha{}_\beta = \frac{2\kappa^2 k}{9\pi}\dot q^2\begin{pmatrix}0&1\\-1&0\end{pmatrix},\qquad k=8\pi G$$

再代入推广的 DeWitt–Brehme–Hobbs 方程 (1.4.19)（含 tail 积分项）求解，得出修正量级 $A\sim e^2/c^6$，即使 $e^2\sim c^3$ 也只有 $\sim c^{-3}$。

### §5 联系实验观测（第四个原创点）

设计判别实验：静止电子受一个电脉冲，测其后的运动。对 $\lambda+\kappa=0$ 模型解出

$$\frac{dt}{d\tau}=\cosh\left(-\frac{E}{m}e^{-\varepsilon\tau}+\frac{E}{m}\right),\quad \varepsilon=\frac{3m}{4\kappa e^2}$$

弛豫时标 $\varepsilon^{-1}$。观测量：定角辐射功率

$$\frac{dP(\tau)}{d\Omega}=\frac{e^2}{4\pi}\frac{(E\varepsilon/m)^2e^{-2\varepsilon\tau}\sin^2\theta/\cosh^6}{(1-\cos\theta\tanh)^5}$$

由 $dP(t)/d\Omega$ 的时间演化拟合 $\varepsilon$ 即 $\kappa$。并指出三种模型的观测特征本质不同：分离奇异项法给出**衰减**的辐射，Eliezer 高阶系列给出**周期性**辐射（电子受脉冲后永远在辐射，"完全符合能量动量守恒，这也正是很奇特的地方"）。

**论文自己的可行性判断**：算出 $\varepsilon^{-1}\sim e^2/mc^3$，电子为 $1.044\times10^{-33}$ s、质子为 $5.684\times10^{-37}$ s，判定"在目前的实验条件下……可以认为这种时间尺度是瞬间完成的，根本无法判断出 Dirac 的解释是否可能"。于是提出需要寻找"大电荷、小质量"的假想粒子。**——这个判断是错的，见 §2.1。**

## 1.3 第二部分：质点的运动

### §1 MiSaTa 方程介绍
完全平行于第一部分：有质量粒子的自引力场使世界线上度规奇异，测地线方程失效。两条路线——Mino–Sasaki–Tanaka 的能动量守恒管形法，和 Poisson 的奇异项／辐射项分解法。论文详细复述了 MiSaTa 三步推导（Hadamard 形式 Green 函数 → 引力辐射能动量张量 → 管形能动量守恒）和 Poisson 的 Detweiler–Whiting 型分解，得到

$$\frac{Du^\mu}{d\tau}=-\tfrac12(g^{\mu\nu}+u^\mu u^\nu)(2h^{\rm tail}_{\nu\lambda;\rho}-h^{\rm tail}_{\lambda\rho;\nu})u^\lambda u^\rho$$

### §2 对两种推导方法的评述（本部分核心）
- **对 Poisson**：摄动展开的前提是 $h_{\alpha\beta}\ll g_{\alpha\beta}$，而解在 $r\to0$ 时 $\propto 1/r\to\infty$，**完全不满足展开条件**。Einstein 方程的非线性使得分离辐射部分比电动力学困难得多。Poisson 本人清楚这点，强调渐进展开，但那要求把质量物体当成小黑洞。论文认为这个不足"暂时无法改善"。
- **对 Mi–Sa–Ta**：(1) $T^{\mu\nu}$ 在世界线上仍发散，用 Gauss 定理的条件不具备——论文的补救是再套一个更窄的管形 $V_{\rm tube}$ 并令其趋零，把那一项当边界条件；(2) cap 面积分 $\int_0^r(r')^{-2}dr'$、$\int_0^r(r')^{-1}dr'$ 都发散，Mi–Sa–Ta 忽略了第三项、第二项只是估计，由此得出的**质量演化 $\dot m(\tau)$ "具有很大的推测成分"**——论文主张质量应像 Dirac 那样退回成外加参数，而不是理论输出。

然后把第一部分的 $B_\mu$ 参数族思路搬过来，得到比 MiSaTa 多一项的方程 (2.2.12) 及其参数化推广 (2.2.14)。

### §3 背景为平直空间
MiSaTa 在 Minkowski 背景下 $V_{\alpha\beta\gamma\delta}\equiv0$，给出 $m\ddot z^\alpha=0$（匀速直线运动）。论文用 Fourier 分解 + 引力平面波能动量张量 $\langle t^{\alpha\beta}\rangle=\frac{k^\alpha k^\beta}{16\pi G}(e^{\lambda\rho}e_{\lambda\rho}-\tfrac12|e^\lambda{}_\lambda|^2)$ 重算，得到

$$\frac{dP^\alpha_{\rm gra}}{d\tau}=\frac{G}{\pi}Su^\alpha \neq 0$$

**结论：平直背景下有质量粒子也未必匀速直线运动**——与 MiSaTa 不一致。

### §4 补充
$t^{\mu\nu}$ 的高阶项（只有偶数阶起作用，仍带 $k_\alpha k_\beta$ 因子）；宇宙学常数（线性化方程形式不变，但底空间度规改变，Minkowski 不再是真空解）；以及最后设想把 LD 改进型与 MiSaTa 改进型相加，得到"既带电又有质量的电子"的综合方程。

## 1.4 作者当时的自我评价

> "本文是对于解决经典粒子运动问题的一个初步的尝试……模型中所带的参数不能由理论本身得到，而要求用实验的方法解决……目前的工作，只能作为算学上的一种尝试。"
> "我个人认为本文的亮点是我们遇到无穷大的物理量所可采取的思路和解决方法，也许对其他问题也有着启发作用。"

---

# 第二部分　二十年后的回顾、展望、什么时候用经典

## 2.1 原文中发现的具体问题

### ① §5 特征时标算错了 10 个量级 —— 直接推翻了论文自己"实验不可行"的结论 ⚠️

论文写 $e^2/mc^3|_{\text{电子}}=1.044\times10^{-33}$ s、质子 $5.684\times10^{-37}$ s。复算发现，这是把 SI 制的库仑数直接平方、却**漏掉了 $1/(4\pi\varepsilon_0)$**（而论文开头明确声明电动力学用高斯制）：

```
SI 裸算   e²/(mc³) = 1.0459e-33 s   ← 与论文的 1.044e-33 吻合
正确(高斯) ke²/(mc³) = 9.400e-24 s
τ₀ = 2e²/3mc³        = 6.266e-24 s
质子: 裸算 5.696e-37 s  vs  正确 5.119e-27 s
```

差一个因子 $1/(4\pi\varepsilon_0)=8.99\times10^9$。**真实的经典辐射反作用时标是 $10^{-23}\sim10^{-24}$ s，不是 $10^{-33}$ s。**

后果：(a) §5 "实验无法判断"的结论不成立；(b) §5 由此提出的"需要寻找大电荷小质量的假想粒子"整段失去必要性；(c) 这个 $\tau_0$ 是后面所有强场判据（$\tau_0\omega_B$）的输入，见 §2.6。

### ② §3 走错了方向：提高导数阶数必然增加 runaway 分支 ⚠️

论文把 $B^\mu$ 扩展到含 $\ddot v,\ v^{(4)},\dots$，得到三阶、五阶 ODE，然后发现所有参数下解都病态。**这不是巧合**：LD 的病根就是三阶项引入了一个正指数特征根（$e^{3\tau/2}$ 那支）；再加高阶导数只会再添更多不稳定分支（Ostrogradsky 型）。

领域后来的答案恰好是**反方向**——**降阶（reduction of order）**：把 $\dot a$ 用 $\dot F_{\rm ext}/m$ 替换，得到 **Landau–Lifshitz (LL) 方程**，二阶、无 runaway、无预加速。Spohn (2000) 用 Fenichel 慢流形理论证明 LD 的物理解全部落在一个"临界流形"上，LL 是它的一阶近似；Gralla–Harte–Wald (2009) 从推广的点粒子极限严格导出。

所以论文纠结的"自然界到底选哪个 $B^\mu$"，现代的回答是：**这些选择在有效场论的每一阶上互为场重定义，物理上等价**——问题被溶解了，不是被解决了。

**另注**：论文的"Eliezer 方程"与文献通行的 **Eliezer / Ford–O'Connell (EFO) 方程不是一回事**。EFO 正是降阶得到的那个（$\dot a\to\dot F/m$），恰恰无 runaway，且今天是实验判别的主要候选之一——它预言存在**加速度上限**，在极强场下与 LL 分道扬镳。论文取的是 Eliezer 1947 综述中另一个 $B^\mu$ 拟设。

### ③ §2 那个 "$\tfrac12(A_{\rm ret}+A_{\rm adv})$ 有影响" 的结论，是整个第一部分的立论基础，但很可能站不住 ⚠️

论文由 $A_S$ 构造 $T^{\alpha\beta}(S)$，积分得 $dP^\alpha/d\tau = a^\alpha/2r + \tfrac23\dot a^\alpha$，由此推断奇异项并非无影响，才有理由放宽成 $\lambda,\kappa$ 族。

**问题在于**：「奇异场自身应力张量的动量流」$\neq$「奇异场施加的力」。总的 $T[F_{\rm ret}]=T[F_S]+T[F_R]+2\times$交叉项，只取 $T[F_S]$ 一项算不出前者。而 **Detweiler–Whiting (2003)** 明确断言、**Harte 与 Pound 后来严格证明**：奇异场只做质量重整化，**不施加力**；全部自力来自正则场 $F_R$。

所以那个多出来的 $\tfrac23\dot a$ 要么是漏了交叉项，要么是附录 A.3 的代数问题。**这是最值得优先核对的一条。**

不过要说明：即使 ③ 成立，$\lambda+\kappa=0$ 那支模型本身（把 LD 的 $2/3$ 系数当自由参数）作为**唯象修正**仍有独立意义——见第三部分方向 D。

### ④ $\kappa>0$ 分支的能量收支没有检查

$\lambda+\kappa=0$ 模型就是 LD 把 $2e^2/3$ 换成 $-4\kappa e^2/3$。$\kappa>0$ 相当于**把辐射项符号翻过来**——runaway 是没了，但辐射阻尼变成了辐射增益，这直接抵触论文自己在引言里列的第三条判据（能量守恒）。论文全篇没做这个检查。需要补一个显式的能量收支计算。

### ⑤ 文献遗漏与笔误（小）

- **MiSaTa 应为 MiSaTaQuWa**：Quinn & Wald (PRD 56, 3381, 1997) 独立得到同一方程，论文完全没引；
- 参考文献【14】"Mion, Y." 应为 **Mino, Y.**；
- 【11】Eliezer, Rev. Mod. Phys. 19 (1947) 是综述文章，与今天所称的 Eliezer 方程（Proc. R. Soc. A 194, 543, 1948；及 Ford & O'Connell 1991）需区分。

## 2.2 经典电子运动论：领域给出的答案

论文的整个框架建立在"LD 方程必须被修改"这个前提上。二十年后，主流答案是**不修方程，而是正确理解它**：

1. **LD 是有效方程，不是基本方程**。它只在"临界流形／慢流形"上有物理意义。Spohn 用 Fenichel 定理证明 Abraham–Lorentz 模型存在一个正规双曲的慢流形，其上无 runaway，且一阶就是 LL 方程。
2. **降阶 = 唯一物理内容**。LL 方程是 LD 在 $\tau_0$ 的一阶展开；EFO 是另一种降阶方案。二者在 $\tau_0$ 的高阶上不同，但都无病态。
3. **严格推导已完成**。Gralla–Harte–Wald (PRD 80, 024031, 2009) 用推广的点粒子极限严格导出了电磁自力，避开了 Dirac 的管形／$B^\mu$ 任意性。
4. **$B^\mu$ 任意性的现代解释**：等价于在有效作用量中加正比于运动方程的项，即场重定义，物理不变。
5. 近年仍有活跃工作：transseries 结构与降阶（2021）、ALD 病态的变分处理（2025）、LD 动力学的后牛顿框架（2025）。

**对本文的意义**：§2、§3 试图通过修改方程来消除病态，这条路在二十年后被证明是不必要的，而且 §3 的方向（升阶）在原理上就会加剧问题。但 §2 的**参数化**思路可以改换目的——不是为了"修好"方程，而是为了**给辐射反作用系数一个可被实验约束的自由度**（见第三部分方向 D）。

## 2.3 实验：从"不可能"到"已完成"

这是二十年里最大的变化，也正是 §5 的题目。

| 年份 | 进展 |
|---|---|
| 2018 | Cole et al. (PRX 8, 011020) 与 Poder et al. (PRX 8, 031004)：$a_0>10$ 激光与 LWFA 电子束（>500 MeV）对撞，**首次观测到辐射反作用**，能损与硬光子信号关联，与量子描述一致 |
| 2024 | 多拍瓦激光全光学非线性 Compton 散射（Nature Photonics） |
| 2024–25 | 量子辐射反作用的解析解（PRA 109, 022234）；LL + 随机扩散与 Monte Carlo QED 的一致性建立 |
| 2026 | **Nature Communications：5σ 观测辐射反作用，给出偏向量子模型、不利于经典描述的定量证据**（全光学对撞位形） |
| 2026 | 经典辐射主导区的精密定量测量（arXiv:2601.10728）；EPJ Plus 综述《Radiation reaction: status and recent developments》 |
| 在建 | LUXE @ DESY（强场 QED 专用实验） |

**对本文的意义**：§5 设计的"给静止粒子一个脉冲、测定角辐射功率随时间演化来拟合模型参数"，在概念上正是今天的激光—电子束对撞实验在做的事。把 §5 的 $\kappa$ 换成可测的辐射反作用系数标度因子，是一个**今天就能做**的题目（方向 D）。

## 2.4 引力自力（第二部分）：从争议到工业化

论文对 MiSaTa 的两条批评（摄动在世界线失效、cap 积分发散、$m(\tau)$ 可疑），后来基本都被解决：

- **Detweiler–Whiting 分解 (2003)**：奇异／正则场的干净分离，奇异场球对称因而不施力，只重整化质量；这条已被 Harte 与 Pound **严格证明**（论文当年只能引用为"Detweiler-Whiting 原则"，并称之为"实用性的原则"）。
- **Gralla–Wald (2008)、Gralla–Harte–Wald (2009)、Pound**：把 MiSaTaQuWa 建立在严格的点粒子极限（一族有限尺寸物体的极限）之上；论文担心的 cap 积分发散问题在这个框架里不再出现。
- **计算技术**：mode-sum 正规化 → puncture / effective source 方案 → 自洽演化 vs. Gralla–Wald 展开。
- **二阶自力**：Schwarzschild 准圆轨道已完成（Pound–Wardell–Warburton–Miller），2025 年推广到**偏心轨道**（PRD 112, 064048）；二阶自力波形在质量比 1:10 时与完全非线性数值相对论的相位差已到**亚弧度**量级；正在为 LISA 的 EMRI 波形库供货；2026 年开始用数值相对论方法攻 Kerr 的一般二阶计算。

**对本文的意义**：第二部分的技术批评已被历史解决，硬啃没有边际收益。真正有价值的是做一次"清算"——把当年提出的每条反对意见与后来的解决方式逐条对照（方向 E）。

## 2.5 §4 的综合问题已被做完

论文 §4 那个"电子既带电又有质量、自身电磁场还要弯曲时空、需要迭代求解"的问题，已有严格版本：

- **Zimmerman & Poisson, PRD 90, 084031 (2014)**，《Combined gravitational and electromagnetic self-force on charged particles in electrovac spacetimes》。核心结果：尽管 Einstein–Maxwell 把扰动的电磁场与引力场耦合起来，**重整化中发生了一个漂亮的相消**——重整化质量只需分别减去纯电荷的电磁贡献与纯质量的引力贡献。这正面回答了论文"世界线上度规无穷大怎么办"的担忧。
- **2026 年 PRD**：给出 DeWitt–Brehme–Hobbs 方程在共传播电磁波 + 引力波背景下的**首个精确解**——即论文式 (1.4.19) 的精确解，并指出引力波的存在会定性改变电磁辐射反作用效应。

论文的"有效度规"假设，实质上就是 Detweiler–Whiting 正则场——只是 DW 在 2003 年才发表，论文（2006）只能从 Poisson 综述里间接引到，并把它当作一个"目前只能用实用性原则代替"的假设。**这个直觉是对的，方向也是对的，只是当时缺少工具。**

## 2.6 强磁场方向的定量检验

这是原始研究动机所在，单独检验。

### (a) 好消息：磁星场强正是经典辐射反作用崩溃的标度

经典 RR 的展开参数是 $\tau_0\omega_B$（辐射反作用时标 × 回旋频率）。令其为 1：

$$\tau_0\omega_B=\frac{2e^3B}{3m^2c^4}=1 \;\Longrightarrow\; B\approx 9\times10^{15}\ {\rm G}$$

$$\boxed{B_{\rm cl}=\frac{m^2c^4}{e^3}=6.05\times10^{15}\ {\rm G}=\frac{B_{\rm Schwinger}}{\alpha}}$$

数值表（$\tau_0=6.266\times10^{-24}$ s，$B_S=4.414\times10^{13}$ G）：

| $B$ (G) | $b=B/B_S$ | $\tau_0\omega_B$ | 第一 Landau 能级 $mc^2(\sqrt{1+2b}-1)$ | 经典冷却时间 $t_{\rm cool}(\gamma{=}10)$ | $\chi=\gamma B/B_S$ |
|---|---|---|---|---|---|
| $10^{12}$ | 0.023 | $1.1\times10^{-4}$ | 0.011 MeV | $5.2\times10^{-17}$ s | 0.23 |
| $10^{13}$ | 0.23 | $1.1\times10^{-3}$ | 0.105 MeV | — | 2.3 |
| $4.41\times10^{13}$ | **1.00** | $4.9\times10^{-3}$ | 0.374 MeV | — | 10 |
| $10^{14}$ | 2.27 | $1.1\times10^{-2}$ | 0.691 MeV | $5.2\times10^{-21}$ s | 23 |
| $10^{15}$ | 22.7 | 0.110 | **2.97 MeV** | **$5.2\times10^{-23}$ s $\approx 8\tau_0$** | 227 |
| $6\times10^{15}$ | 136 | 0.66 | 7.93 MeV | — | — |
| $10^{16}$ | 227 | **1.10** | 10.4 MeV | — | 2266 |

**在 $10^{15}$–$10^{16}$ G，LD 方程里那个"病态"的 $\tau_0$ 尺度和真实动力学尺度合并了。** 在 $B=10^{15}$ G、$\gamma=10$ 时，经典同步冷却时间只有约 $8\tau_0$——runaway 时标与物理时标已经分不开。

**当年"强场下经典电子运动论的病态会变成 $O(1)$ 效应"这个直觉，数量上完全正确**，而且特征场强 $B_S/\alpha$ 恰好落在磁星量级。这不是巧合，是量纲分析的必然。

### (b) 坏消息：QED 提前 $1/\alpha=137$ 倍到达

Schwinger 场 $B_S=4.41\times10^{13}$ G 比 $B_{\rm cl}$ 低整整 137 倍，而磁星表面场已超过 $B_S$ 一到两个量级。所以那个"经典 RR 变成非微扰"的窗口，**整个埋在 QED 区里面，埋得多深恰好是 $1/\alpha$**：

- **Landau 量子化**：$B=10^{15}$ G 时第一激发 Landau 能级在 2.97 MeV。电子被钉在基态 Landau 能级，横向运动冻结 → **一维沿磁力线运动**。经典回旋运动在磁星磁层里根本不存在，"修正回旋运动的 LD 方程"是无的放矢。
- **$\chi=\gamma B/B_S$**：$\gamma=10$ 就已经 $\chi\approx200$，$\gamma=10^3$ 时 $\chi\approx2\times10^4$。深度量子区，经典辐射公式连零阶都不对。
- 单光子对产生、光子分裂、真空双折射（IXPE 2022 在 4U 0142+61 上已测到）全部开启。

**结论**：磁星是"经典点电荷运动论"这条路的终点，不是它的用武之地。这个结论本身干净、可量化，且解释了为什么当年那条"改造经典 LD 方程"的路在磁星语境下注定走不通——不是不重要，是被提前 137 倍到达的 QED 覆盖掉了。

### (c) 顺带结掉 §4 在磁星中的地位

磁能密度的等效质量密度：

| $B$ (G) | $u_B=B^2/8\pi$ (erg/cm³) | $\rho_{\rm equiv}$ (g/cm³) | / 核密度 |
|---|---|---|---|
| $10^{15}$ | $3.98\times10^{28}$ | $4.4\times10^{7}$ | $1.6\times10^{-7}$ |
| $10^{16}$ | $3.98\times10^{30}$ | $4.4\times10^{9}$ | $1.6\times10^{-5}$ |
| $10^{18}$ | $3.98\times10^{34}$ | $4.4\times10^{13}$ | 0.16 |

即：**"自场弯曲时空"在磁星里是一个良定义的小修正**（且磁化中子星的结构与形变早有成熟文献），不是当年那个发散的点粒子问题。只有假想的内部场 $\sim10^{18}$ G 才达到核密度的量级。

## 2.7 "Aristotle 运动" = Aristotelian Electrodynamics

这是整份回顾里最值得注意的一条。

强辐射阻尼下的加速，现在有标准名字：**radiation-reaction-limited regime**，其极限形式在脉冲星磁层建模里就叫 **Aristotelian electrodynamics (AE)**，定义是：

> *"由于辐射过阻尼，局部电磁场决定电荷的速度、而不是加速度。"*

而论文 §2 在讨论 $\kappa>0$ 那支解时写的是：

> *"这里电子运动的表现性质很像 Aristotle 所提到的物质运动现象，不妨称之为'Aristotle 运动'。"*

**同一个物理**：辐射阻尼强到把牛顿二阶动力学压成一阶的场—速度代数关系。Gruzinov 在 2012–13 正式提出 AE，此后 Contopoulos (2016)、Pétri (2020)、Cerutti 等把它做成脉冲星磁层耗散建模的标准工具（三维耗散磁层、斜转子、乃至 AR Sco 建模）。

**也就是说：当年被当作"$\kappa$ 参数一个奇怪分支"的东西，恰好是强场天体物理里的正确极限。** 只不过正确的到达路径不是"修改 LD 方程"，而是"在标准 LL 方程的辐射主导极限下取过阻尼"。

### 强场下加速机制确实不同（回答原始问题 2）

1. **横向动量在 $10^{-21}$–$10^{-23}$ s 内被辐射殆尽** → 粒子瞬间落到最低 Landau 能级，一维沿磁力线运动；
2. **唯一有效的加速是平行电场 $E_\parallel$**（gap / 磁层 twist），并被**曲率辐射反作用**封顶：$\gamma_{\rm RRL}=(3E_\parallel\rho^2/2e)^{1/4}\sim10^7$；
3. **磁星里主导能损不是同步辐射而是共振逆 Compton 散射 (RICS)**（$\gamma\sim10$–$10^3$ 时）——本质是回旋吸收 + 中间 Landau 态再发射的一阶 QED 过程；
4. RICS 光子按偏振分别走单光子对产生和光子分裂，触发多代 pair cascade（2025 ApJ）；
5. 这一整套正是磁星硬 X 射线尾、IXPE 偏振、以及 FRB 模型的核心。

**已知的开放缺口**：AE 并非处处成立——磁层里存在辐射阻尼不够强的局部区域，现在的做法是 force-free 与 AE 拼接；且 PIC 模型因等离子体频率与自转频率的巨大尺度分离，尚无法用真实脉冲星参数自洽建模。

---

## 2.8 什么时候用经典，什么时候必须用量子

这一节把判据讲清楚，因为它决定了后面所有工作的适用边界。

### 一句话判据

看每次辐射出去的那一个光子，拿走了电子多少能量。拿走一点点，用经典。拿走一大块，必须用量子。

### 为什么是这个

经典描述里电子的能量是连续、平滑地漏出去的，像水从桶底慢慢渗。真实情况是电子一个光子一个光子地往外扔能量。

如果电子有一百块钱，每个光子只花掉一分钱，那"连续漏水"这个近似非常好，一分一分地花和连续地流没有区别。如果一个光子一下子花掉五十块，那就完全不是漏水了，电子的能量是被一大口一大口随机咬掉的，而且咬完之后运动方向都变了。这时候还用"平滑漏水"的方程去算，答案会错得离谱。

量子非线性参数 $\chi$ 量的正是"一个光子大概占电子总能量的多少"：

$$\chi=\frac{\gamma B}{B_S},\qquad B_S=\frac{m^2c^3}{e\hbar}=4.41\times10^{13}\ {\rm G}$$

- $\chi\ll1$：一分一分地花，**用经典**
- $\chi\gtrsim1$：一口咬掉一大半，**必须用量子**

两个因素让 $\chi$ 变大：电子跑得越快（$\gamma$ 越大）辐射出的光子越硬；磁场越强（$B$ 越大）电子被掰弯得越急、辐射越硬。

### 另一件独立的事：轨道本身被量子化

上面说的是辐射。另有一件事关系到运动本身。

磁场里的电子走螺旋线，场越强圈越小。强到一定程度，圈小到量子力学不允许它连续变化，电子只能待在分立的 Landau 能级上。门槛同样是 $B_S$。**超过 $B_S$，横向打转被冻住，电子只能沿磁力线一维地走**，这时候连"回旋运动"这个图像都不存在了。

### 三条判据合起来

|判据|临界值|触发之后什么失效|
|---|---|---|
|$\chi=\gamma B/B_S$|$\chi\gtrsim1$|光子反冲不可忽略，经典辐射公式连零阶都不对|
|$b=B/B_S$|$4.41\times10^{13}$ G|Landau 量子化，横向运动冻结|
|$\tau_0\omega_B$|$B_{\rm cl}=6.05\times10^{15}$ G|经典辐射反作用变成非微扰，LD 的病态时标与物理时标合并|

### 关键：这三条有先后，而且顺序是致命的

$$B_{\rm cl}=\frac{B_S}{\alpha}$$

经典理论自己崩掉的场强，比量子接管的场强**高 137 倍**。这不是巧合，$\alpha=e^2/\hbar c$ 就是这两个尺度的比。

后果是：**经典理论从来没有机会自己崩掉。**你还没等到经典方程出问题，量子力学早就把场子接管了。

### 落到具体天体上

|地方|$B$|经典还能用吗|
|---|---|---|
|普通脉冲星表面|$10^{12}$ G|勉强，但电子稍快一点（$\gamma\gtrsim5$）就出界|
|Schwinger 场|$4.41\times10^{13}$ G|**这是分界线**|
|磁星表面|$10^{15}$ G|完全不行，螺旋运动都没有了，$\gamma=10$ 时 $\chi\approx227$|
|经典自己崩掉处|$6\times10^{15}$ G|永远轮不到它|

### 经典的东西今天在哪还活着

三个地方：

1. **实验室。**激光对撞实验刻意在 $\chi\sim0.1$ 到 1 这个过渡区工作，那里经典与量子的预言开始分开，正好用来判别。这也是辐射反作用系数唯一能被卡住的地方。
2. **脉冲星磁层的低场区。**远离表面处 $B$ 降下来，经典处理重新可用，Aristotelian electrodynamics 就在这一区。
3. **作为量子结果的极限检查。**强场 QED 的公式在 $\chi\to0$ 时必须回到经典，这是标准的一致性检验。

## 2.9 两代展望的对照

### 当年这篇论文的展望

论文结尾自己写的是：

> "本文是对于解决经典粒子运动问题的一个初步的尝试……模型中所带的参数不能由理论本身得到，而要求用实验的方法解决……目前的工作，只能作为算学上的一种尝试。"

拆开看，当年押了三件事：

1. **参数要靠实验定。**这个判断是对的，而且是全文最有生命力的一句。
2. **修改方程能治好病态。**这个方向错了，领域后来走的是降阶而非升阶。
3. **判别实验不可行。**这个结论建立在算错 10 个量级的特征时标上，是错的，而实验在 2018 年就做出来了。

三押一中，中的那一条恰好是唯一还能往下走的那条。

### 领域现在的展望

两支的驱动力完全不同。

**经典电子运动论这一支，驱动力是实验设施。**形式理论的争论已经关闭，现在的问题是拿数据去卡。在跑的和在建的有 SLAC FACET-II 的 E320（10 GeV 电子束对撞 10 TW 激光，目标是在电子静止系达到 Schwinger 强度）和 DESY 的 LUXE（强场 QED 专用，2026 年开始取数）。它们的展望是把 $\chi\sim1$ 这个过渡区扫干净，判别经典与量子描述，并给出定量精度。

**引力自力这一支，驱动力是 LISA 的截止日期。**EMRI 波形库需要二阶自力。一阶在 Kerr 上的一般轨道 2018 年做完，二阶目前只做到 Schwarzschild 的圆轨道与偏心轨道，**Kerr 的一般二阶还没做出来**，这是当前的主攻方向，2026 年的新思路是借数值相对论的方法去做（arXiv:2606.04998）。

### 两支的可进入性差别很大

|  |经典电子运动论|引力自力|
|---|---|---|
|形式理论|已关闭|一阶关闭，二阶在攻|
|活跃前沿|实验对撞、强场 QED|Kerr 二阶、LISA 波形库|
|社区规模|中等，分散在激光等离子体与 QED 两边|小而紧密|
|工具链|公开工具多，可复算|自建代码，半公开|
|外部人能不能插进去|**部分可以**|**基本不行**|

引力自力那一支已经工业化，它需要的是能写代码接上现有工具链的人，不是提想法的人。经典这一支还有缝，原因是它现在分成互不引用的两拨：做激光对撞的实验组在测辐射反作用系数偏不偏离标准值，做磁层建模的天体组在用 Aristotelian electrodynamics，而 AE 恰恰就是强阻尼极限下的同一个方程。**那个自由系数横跨这两边，而没有人把两边的数放在一起。**

---

# 第三部分　接下来可以做什么

结论先写在前面：**只有一条主线值得做，就是把辐射反作用系数变成一个可被数据约束的参数。**其余都是它的配件或者补充。

## 3.1 主线：辐射反作用系数的实验约束

### 问题的形式

把 LD 方程的辐射反作用系数放开成一个自由参数：

$$ma^\alpha = F^\alpha_{\rm ext} + \zeta\,\tau_0\,m\left(\dot a^\alpha - a^2u^\alpha\right),\qquad \zeta=1\ \text{为标准值}$$

论文 §2 那个 $\kappa$ 族正是 $\zeta=-2\kappa$。

要问的问题是一句话：**现有数据把 $\zeta$ 卡在什么区间？**

这个问题有明确的参数、有现成的数据、而且据检索没人系统做过。今天的实验判的都是"哪个方程形式对"（经典对量子、LL 对 LAD），而 LAD、Mo–Papas、LL、EFO 在经典适用区里的 Lorentz 因子相对差异小于 $10^{-6}$，**形式之争在那个区间几乎测不出来，系数之争才是能测的那个**。

### 第 0 步：先算敏感度（半天，决定这条线成不成立）

方程不能直接用，三阶形式带 runaway，第一步必须降阶，把 $\dot a$ 换成 $\dot F_{\rm ext}/m$，得到带自由 $\zeta$ 的 Landau–Lifshitz 形式。这一步是纯代数。

降阶之后在均匀强磁场里有闭式解，同步冷却率整体正比于 $\zeta$。于是每一个"加速功率等于辐射损失"的平衡条件都变成 $\zeta$ 的函数：

- $\gamma_{\rm RRL}$（辐射反作用限制的 Lorentz 因子）
- 同步 burn-off 封顶（理想 MHD 下约 160 MeV）
- 曲率辐射的上限

**要算的就是这三条各自对 $\zeta$ 的幂次依赖。**

这一步决定后面还有没有戏：

- 幂次浅（比如 $\zeta^{1/4}$），$\zeta$ 差一倍才移动 20%，天体观测卡不住它，这条线只能走实验室对撞
- 幂次陡（$\zeta^{1}$ 或更陡），磁星与 FRB 的观测能直接给出 $\zeta$ 的区间，那就是一篇有实质结论的论文

**这一步必须排在所有事情前面**，包括那张相图。图是把已知边界画出来，这一步是判断这些边界会不会动。

### 第 1 步：修两个 bug（周末，地基）

1. **单位错**（§2.1 ①）：已确认并复算，需在任何后续使用 $\tau_0$ 的地方改正。
2. **附录 A.3 的 $T[F_S]$ 积分**（§2.1 ③）：用符号计算在 Retarded 坐标下独立重算，判定 $\tfrac23\dot a$ 项是否应当出现。这一条决定 $\lambda$-$\kappa$ 族的立论是否成立，也就决定了 $\zeta$ 这个参数**能不能说出自你自己的工作**。作为唯象参数它照样能用，但来源要改写。
3. 顺带核对 §3 "所有 $(A,B)$ 都 runaway" 是真结论还是数值假象。

### 第 2 步：$(B,\gamma)$ 参数相图（一两周，论文的核心图）

在对数 $(B,\gamma)$ 平面上画出全部边界：$\tau_0\omega_B=1$、$B=B_S$、$\chi=1$、$t_{\rm cool}=\tau_0$、单光子对产生阈值、$\gamma_{\rm RRL}$ 曲线、经典适用区的楔形。再把普通脉冲星、磁星表面与磁层、FRB 源、以及实验室激光实验标上去。

**这张图单独拿出来只是一个 finding，撑不起一篇论文**，但作为主线论文的核心图它是必需的：它一眼说明"经典能用的那块楔形在哪，离磁星有多远"，以及第 0 步算出的那些边界会往哪个方向移。

公式在 §2.6 与 §2.8 已全部列出，核心代码二三十行。

### 第 3 步：拿数据卡 $\zeta$

1. 用降阶形式正演（避免 runaway 污染拟合）；
2. 在 2018 年 Gemini 与 2026 年那个 5σ 对撞位形下，算电子能谱与光子谱对 $\zeta$ 的依赖；
3. 拿公开数据反推 $\zeta$ 的置信区间；
4. 把论文 §5 强调的观测量 $dP(t)/d\Omega$ 翻译成全光学对撞里的角分辨 Compton 谱。

先读 EPJ Plus 2026 那篇《Radiation reaction: status and recent developments》，弄清实验组现在怎么参数化偏离、精度到什么量级。一天的事，但能省掉重复劳动。

### 这条线为什么不是从零开始

经典电子运动论现在分成互不引用的两拨人：做激光对撞的实验组，和做磁层建模的天体组。$\zeta$ 正好横跨两边。这不是自己开一条路，是站在两个已有社区中间做一次翻译。

## 3.2 三个补充方向

以下都不是主线，视兴趣与时间补。

**A. 磁星 cascade / RICS 中辐射反作用处理的自洽性审计。**现在的动理学模拟用 LL 加 Monte Carlo 随机发射，但在 $b=B/B_c>1$ 时套用的仍是 $b\ll1$ 推出的发射公式。系统检查这些近似在磁星参数下的有效范围。具体、有限、有人会读，不需要新的形式理论。

**B. 证明 Eliezer 系列的参数是冗余的。**明确证明论文 (1.3.3) 那一族方程通过降阶在每一阶 $\tau_0$ 上塌缩到同一个 LL/EFO，即 $B^\mu$ 自由度等价于场重定义。这样论文最困惑的那句"自然界难道真会去选择更复杂的 $B_\mu$"就有了确定答案：**它不选，因为选什么都一样。**纯解析，中等工作量。

**C. 第二部分的清算笔记。**把 Part 2 当年提出的每条反对意见与后来的解决方式逐条对照（Gralla–Wald 点粒子极限、Harte–Pound 奇异场不施力定理、二阶自力中 $m(\tau)$ 的地位、Quinn–Wald 的独立推导）。价值主要是对作者自身的，能看清当年的物理判断力准在哪、偏在哪。

## 3.3 明确不做的

**引力自力那一支。**理由见 §2.9：它已经工业化，工具链自建、有 LISA 的硬截止日期，需要的是能接上现有代码的人。当年那些技术批评的价值在"清算笔记"，不在"参与"。

## 3.4 执行顺序

```
第 0 步  算 ζ 的幂次依赖      半天    ← 决定这条线成不成立，必须最先做
   ↓
第 1 步  修两个 bug           周末    ← 决定 ζ 能不能挂自己的名
   ↓
第 2 步  (B,γ) 相图           一两周  ← 论文的图 1
   ↓
第 3 步  拿公开数据卡 ζ       主体
   ↓
3.2 的三个补充方向（可选）
```

关键点是把顺序倒过来：**先做那个能一票否决整篇的计算，再做工作量大的部分。**

# 附录　数值复算

以下脚本复现本文所有数值（高斯制 CGS）：

```python
import math
e=4.80320471e-10; m=9.1093837015e-28; c=2.99792458e10; hbar=1.054571817e-27

tau0 = 2*e**2/(3*m*c**3)          # 6.266e-24 s
BS   = m**2*c**3/(e*hbar)         # 4.414e13 G  (Schwinger)
alpha= e**2/(hbar*c)              # 1/137.04
Bcl  = m**2*c**4/e**3             # 6.049e15 G  = BS/alpha

# 论文 §5 的单位错还原
e_SI=1.602176634e-19; m_SI=9.1093837015e-31; c_SI=2.99792458e8; k=8.9875517923e9
print(e_SI**2/(m_SI*c_SI**3))     # 1.0459e-33  <- 论文的 1.044e-33（漏了 1/4πε₀）
print(k*e_SI**2/(m_SI*c_SI**3))   # 9.400e-24   <- 正确值

# 强场判据
for B in [1e12,1e13,4.414e13,1e14,1e15,6e15,1e16,1e18]:
    wB = e*B/(m*c)                              # 非相对论回旋频率
    E1 = 0.511*(math.sqrt(1+2*B/BS)-1)          # 第一 Landau 能级 (MeV)
    print(B, B/BS, tau0*wB, E1)

# 经典同步冷却时间 (beta_perp = 1)
def t_cool(B, g):
    P = (2/3)*e**4*B**2*g**2/(m**2*c**3)
    return g*m*c**2/P

# 曲率辐射的辐射反作用极限 Lorentz 因子
def gamma_RRL(E_par, rho):        # E_par in statvolt/cm, rho in cm
    return (3*E_par*rho**2/(2*e))**0.25
```

**关键数字速查**：

| 量 | 值 |
|---|---|
| $\tau_0=2e^2/3mc^3$ | $6.266\times10^{-24}$ s |
| $B_{\rm Schwinger}=m^2c^3/e\hbar$ | $4.414\times10^{13}$ G |
| $B_{\rm cl}=m^2c^4/e^3=B_S/\alpha$ | $6.05\times10^{15}$ G |
| $\tau_0\omega_B=1$ 处 | $B\approx9\times10^{15}$ G |
| 论文 §5 的错误值 | $1.044\times10^{-33}$ s（差 $8.99\times10^{9}$） |

---

# 参考文献

## 原文引用的（1906–2004）

- Poincaré, *Rend. Circ. Mat. Palermo* **21**, 129 (1906)
- Lorentz, *The Theory of Electrons* (1909)
- Dirac, *Proc. R. Soc. London* **A167**, 148 (1938)
- Eliezer, *Rev. Mod. Phys.* **19**, 147 (1947)
- DeWitt & Brehme, *Ann. Phys. (N.Y.)* **9**, 220 (1960)
- Hobbs, *Ann. Phys. (N.Y.)* **47**, 141 (1968)
- Mino, Sasaki & Tanaka, *Phys. Rev. D* **55**, 3457 (1997)
- Landau & Lifshitz, *The Classical Theory of Fields*
- 张宗燧《电动力学及狭义相对论》
- Poisson, [gr-qc/0012057](https://arxiv.org/abs/gr-qc/0012057); *The Motion of Point Particles in Curved Spacetime* (2004)
- Weinberg《引力论和宇宙论》第十章

## 原文遗漏的同期工作

- **Quinn & Wald, *Phys. Rev. D* 56, 3381 (1997)** — 与 MiSaTa 独立、同时得到相同方程（"MiSaTaQuWa" 的后半）
- **Detweiler & Whiting, *Phys. Rev. D* 67, 024025 (2003)** — 奇异／正则分解
- Ford & O'Connell, *Phys. Lett. A* **157**, 217 (1991) — EFO 方程

## 经典自力与降阶（2000–2026）

- [Spohn, *Dynamics of charged particles and their radiation field*](https://arxiv.org/pdf/math-ph/9908024) — 慢流形 / 临界流形
- Gralla, Harte & Wald, *Phys. Rev. D* **80**, 024031 (2009) — 电磁自力的严格推导
- [Reduced-order Abraham-Lorentz-Dirac equation and the consistency of classical electromagnetism](https://arxiv.org/pdf/1402.1106)
- [Reduction of Order and Transseries Structure of Radiation Reaction](https://arxiv.org/pdf/2112.10235)
- [Radiation reaction and limiting acceleration (EFO), *Phys. Rev. D* 105, 016024](https://arxiv.org/abs/2112.04444)
- [Variational Resolution of the Abraham–Lorentz–Dirac Equation Pathologies (2025)](https://arxiv.org/html/2505.02870v2)
- [A Framework for Lorentz-Dirac Dynamics and Post-Newtonian (2025)](https://arxiv.org/pdf/2512.18637)

## 实验（2018–2026）

- [Cole et al., *Phys. Rev. X* 8, 011020 (2018)](https://arxiv.org/abs/1707.06821) — 首个辐射反作用实验证据
- [Poder et al., *Phys. Rev. X* 8, 031004 (2018)](https://journals.aps.org/prx/abstract/10.1103/PhysRevX.8.031004) — 量子辐射反作用的实验特征
- [All-optical nonlinear Compton scattering with a multi-petawatt laser, *Nature Photonics* (2024)](https://www.nature.com/articles/s41566-024-01550-8)
- [Analytical solutions for quantum radiation reaction, *Phys. Rev. A* 109, 022234 (2024)](https://arxiv.org/abs/2312.03592)
- [**Observation of quantum effects on radiation reaction in strong fields, *Nature Communications* (2026)**](https://www.nature.com/articles/s41467-025-67918-8) — 5σ
- [Towards precision quantitative measurement of radiation reaction (2026)](https://arxiv.org/abs/2601.10728)
- [Radiation reaction: status and recent developments, *EPJ Plus* (2026)](https://link.springer.com/article/10.1140/epjp/s13360-026-07436-8) — 综述
- [Input to the European Strategy for Particle Physics: Strong-Field QED (2025)](https://arxiv.org/pdf/2504.02608) — 含 LUXE

## 引力自力（2011–2026）

- [Poisson, Pound & Vega, *Living Rev. Relativity* 14, 7 (2011)](https://pmc.ncbi.nlm.nih.gov/articles/PMC5255936/) — 标准综述
- [The Capra Research Program for Modelling EMRIs](https://arxiv.org/pdf/1102.2857)
- [Zimmerman & Poisson, *Phys. Rev. D* 90, 084031 (2014)](https://journals.aps.org/prd/abstract/10.1103/PhysRevD.90.084031) — 电磁 + 引力联合自力
- [Toward second-order self-force for eccentric EMRIs in Schwarzschild, *Phys. Rev. D* 112, 064048 (2025)](https://arxiv.org/abs/2504.09640)
- [Self-force calculations with numerical relativity methods (2026)](https://arxiv.org/html/2606.04998)
- [Exact solution of the DeWitt-Brehme-Hobbs equation in copropagating EM and gravitational waves, *Phys. Rev. D* (2026)](https://journals.aps.org/prd/abstract/10.1103/z4xf-rqcq)

## 强场天体物理 / Aristotelian Electrodynamics

- [Pétri, *Three-dimensional Dissipative Pulsar Magnetospheres with Aristotelian Electrodynamics* (2020)](https://arxiv.org/pdf/1912.00335)
- [Radiative pulsar magnetospheres: oblique rotators (2022)](https://arxiv.org/pdf/2202.12712)
- [Pulsar gamma-ray emission in the radiation reaction regime, *MNRAS* 484, 5669](https://academic.oup.com/mnras/article/484/4/5669/5307092)
- [Testing convergence to Aristotelian electrodynamics (AR Sco), *MNRAS* 540, 3863 (2025)](https://academic.oup.com/mnras/article/540/4/3863/8161028)
- [Particle acceleration and radiation reaction in a strongly magnetised rotating dipole, *A&A* (2022)](https://www.aanda.org/articles/aa/full_html/2022/10/aa43634-22/aa43634-22.html)
- [Radiative back-reaction on charged particle motion in NS dipole magnetospheres (2024)](https://arxiv.org/pdf/2412.04996)
- [General-relativistic and non-ideal radiative cooling in NS magnetospheres (2026)](https://arxiv.org/pdf/2603.25949)

## 磁星 QED

- [Resonant Inverse Compton Scattering and Hard X-ray Emission in Magnetar Magnetospheres (2026)](https://arxiv.org/abs/2605.07255)
- [Pair Cascades in Magnetar Magnetospheres, *ApJ* (2025)](https://iopscience.iop.org/article/10.3847/1538-4357/adfa06)
- [Impact of Resonant Compton Scattering on Magnetar X-Ray Polarization with QED Vacuum Resonance (2026)](https://arxiv.org/html/2603.08119)
- [Synchro-curvature radiation of charged particles in strong curved magnetic fields](https://arxiv.org/pdf/1501.04994)

---

*文档生成于 2026-08-10。第二、三部分的文献状态基于当时的检索结果；具体数值均已本地复算（见附录）。第 2.1 节 ① 的单位错为已确认；③ 的 $T[F_S]$ 积分问题为待核对项。*
