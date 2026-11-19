# C1：零阶运动与一阶远区方程

> 草稿 · 2026-11-19
> 前提：C0 的假设 (i)–(iii)；这一份不用 (H)。
> 目标：(1) 从假设推出世界线满足洛伦兹力方程，复现 Gralla 2010 式 (52)、Gralla 2013 式 (224)；(2) 写出一阶远区扰动 $(h,a)$ 满足的耦合方程组，作为 C2 的起点。
> 数值检验：`check_c1_identities.py`，检验本文的式 (1)、(3)、(4)。

记号同 C0。背景记为 $(g,A,F)$；$\lambda=0$ 的量不再加下标 0；一阶扰动记为 $h_{ab}=\partial_\lambda g_{ab}|_0$、$a_a=\partial_\lambda A_a|_0$，并记 $f=da$，$h=g^{ab}h_{ab}$。

---

## 1. 两个恒等式

对**任意**度规和任意 $F=dA$，不要求满足场方程：

$$4\pi\,\nabla^aT_{ab}=-F_b{}^c\,\mathcal M_c
\quad\Longrightarrow\quad
\nabla^aE_{ab}=2F_b{}^c\,\mathcal M_c, \tag{1}$$

$$\nabla_a\mathcal M^a=0. \tag{2}$$

**推导 (1)。**
- 展开：$4\pi\nabla^aT_{ab}=(\nabla^aF_{ac})F_b{}^c+F^{ac}\nabla_aF_{bc}-\tfrac12F^{cd}\nabla_bF_{cd}$。
- 由 $dF=0$，中间一项等于 $\tfrac12F^{ac}\nabla_bF_{ac}$，与最后一项相消。
- 剩下的 $\nabla^aF_{ac}=-\mathcal M_c$ 给出 (1)。
- 第二式由 $\nabla^aG_{ab}=0$ 得到。
- 物理含义：$\mathcal M=4\pi j$ 时，(1) 就是 $\nabla^aT^{\rm EM}_{ab}=-F_{bc}j^c$。

**推导 (2)。**$\mathcal M^a$ 是反对称张量的散度，再取一次散度恒为零。

这两式是 Gralla 2010 式 (43)、(44) 在我们归一化下的写法，数值检验见 (a)。

## 2. 线性化

背景是电真空解：$E_{ab}=0$，$\mathcal M^a=0$。对 $(g+\lambda h,\,A+\lambda a)$ 线性化，得到：

$$\mathcal M^{(1)a}\equiv\delta(\nabla_bF^{ab})=\nabla_bY^{ab},\qquad
Y^{ab}=f^{ab}+\tfrac12hF^{ab}-2h^{c[a}F_c{}^{b]}. \tag{3}$$

一般背景下右边还多一项 $-\tfrac12h\,\nabla_bF^{ab}$，背景满足无源 Maxwell 方程时为零。数值检验 (b) 检验的是带这一项的一般式。

$$\delta(4\pi T_{ab})=2F_{(a}{}^cf_{b)c}-h^{cd}F_{ac}F_{bd}-\tfrac14h_{ab}F^2-\tfrac12g_{ab}\big(f_{cd}F^{cd}-h^{ce}F_c{}^dF_{ed}\big), \tag{4}$$

$$E^{(1)}_{ab}=\delta G_{ab}[h]-2\,\delta(4\pi T_{ab}). \tag{5}$$

右边的指标升降和 $\nabla$ 都用背景度规。数值检验 (c) 检验的是 (4)。

对 (1)、(2) 线性化时，凡是乘在背景 $E_{ab}$、$\mathcal M^a$ 上的项都为零，所以

$$\nabla^aE^{(1)}_{ab}=2F_b{}^c\,\mathcal M^{(1)}_c,\qquad \nabla_a\mathcal M^{(1)a}=0. \tag{6}$$

两式都是**对任意 $(h,a)$ 成立的恒等式**，不要求 $(h,a)$ 满足方程，所以对分布意义下的 $(h,a)$ 也成立。

## 3. 分布源

由 C0 第 3 节，$(h,a)$ 在 $\gamma$ 附近最多像 $1/r$ 那样奇异，在 $r>0$ 上满足 $E^{(1)}=0$、$\mathcal M^{(1)}=0$。把它们看成整个背景时空上的分布，两个算子都是二阶线性微分算子。照 Gralla 2010 附录 A 的论证，有

$$E^{(1)}_{\mu\nu}=N_{\mu\nu}(t)\,\delta^3(x^i),\qquad \mathcal M^{(1)\mu}=N^\mu(t)\,\delta^3(x^i). \tag{7}$$

右边只有 $\delta^3$，不出现 $\delta^3$ 的导数，因为场最多只有 $1/r$ 的奇异性。这一步对耦合方程组是否照搬成立，列在 C0 的待核清单里。

## 4. 在 Fermi 坐标里求解 (6)

取背景度规下沿 $\gamma$ 的 Fermi 坐标。在 $\gamma$ 上：
- 度规 $g_{\mu\nu}=\eta_{\mu\nu}$；
- 联络只有 $\Gamma^i{}_{00}=\Gamma^0{}_{0i}=a^i$ 不为零，$a^i$ 是 $\gamma$ 的加速度；
- 逆度规 $g^{\mu\nu}$ 的空间一阶导数为零，所以光滑函数乘 $\partial_i\delta^3$ 时不产生额外的 $\delta^3$ 项。

**Maxwell 部分。**
$$\nabla_a(N^a\delta^3)=\dot N^0\delta^3+N^i\partial_i\delta^3+a_iN^i\delta^3=0 .$$
$\partial_i\delta^3$ 与 $\delta^3$ 线性无关，所以 $N^i=0$、$\dot N^0=0$。归一化 $N^0\equiv4\pi q$，用 Coulomb 场检验。

**Einstein 部分。**
- 含 $\partial_i\delta^3$ 的项只有 $N_{ib}\,\partial_i\delta^3$，而 (6) 右边只有 $\delta^3$，所以 $N_{ib}=0$，即 $N_{0i}=N_{ij}=0$。
- 剩下 $S_{ab}=N_{00}\,\delta^3\,t_at_b$。直接计算得
$$(\nabla^aS_{ab})_0=-\dot N_{00}\,\delta^3,\qquad (\nabla^aS_{ab})_i=a_i\,N_{00}\,\delta^3 .$$
- (6) 右边 $2F_b{}^cN_c\,\delta^3$，其中 $N_c$ 只有 $N_0=-4\pi q$。
- $b=0$：$F_0{}^0=0$，所以 $\dot N_{00}=0$。归一化 $N_{00}\equiv8\pi M$，用线性化 Schwarzschild $h_{00}=2M/r$ 检验。
- $b=i$：$8\pi M\,a_i=2F_i{}^0(-4\pi q)=8\pi q\,F_{i0}$。

**协变写法。**

$$E^{(1)}_{ab}=8\pi M\!\int_\gamma u_au_b\,\delta_4(x,z(\tau))\,d\tau,\qquad
\mathcal M^{(1)a}=4\pi q\!\int_\gamma u^a\,\delta_4(x,z(\tau))\,d\tau, \tag{8}$$

$$M\,u^b\nabla_bu^a=q\,F^a{}_bu^b,\qquad \dot M=\dot q=0. \tag{9}$$

这就是 Gralla 2013 式 (224) 和 Gralla 2010 式 (52)，只是归一化不同。

- **$M\neq0$ 是必需的**（GW 的注记）。$M=0$ 时 (9) 退化为约束 $qF_{ab}u^b=0$。对黑洞，$M\ge|q|>0$，没有问题。
- **$M$、$q$ 的含义**：$M$ 是近区背景的 ADM 质量，$q$ 是它的 Gauss 电荷（Gralla 2013 式 (222)、(223)）。在 (H) 下，它们就是 Kerr–Newman 的参数，$\xi=q/M$。

**小结：零阶完全复现，没有新东西。**新的只有一点：由 (6) 可以直接看出，正是耦合项 $2F_b{}^c\mathcal M_c$ 把"点源守恒"变成了洛伦兹力。

## 5. 一阶远区方程组

把 (3)–(5)、(8) 合起来，$(h,a)$ 在整个背景时空上分布意义下满足：

$$\delta G_{ab}[h]-2\Big[2F_{(a}{}^cf_{b)c}-h^{cd}F_{ac}F_{bd}-\tfrac14h_{ab}F^2-\tfrac12g_{ab}\big(f_{cd}F^{cd}-h^{ce}F_c{}^dF_{ed}\big)\Big]
=8\pi M\!\int u_au_b\,\delta_4\,d\tau, \tag{10}$$

$$\nabla_b\Big[f^{ab}+\tfrac12hF^{ab}-2h^{c[a}F_c{}^{b]}\Big]=4\pi q\!\int u^a\,\delta_4\,d\tau. \tag{11}$$

**规范。**允许的变换是
$$h\to h+\mathcal L_\zeta g,\qquad a\to a+\mathcal L_\zeta A+d\psi,\qquad
\mathcal L_\zeta A=\zeta^cF_{ca}+d(\zeta^cA_c).$$
- 背景满足场方程，所以 (10)、(11) 的左边在这些变换下不变。
- **两种规范混在一起**：微分同胚 $\zeta$ 会让 $a$ 改变 $\zeta^cF_{ca}$，模去电磁规范。这正是 C2 的难点之一。
- Lorenz 规范 $\nabla^b(h_{ab}-\tfrac12g_{ab}h)=0$、$\nabla^aa_a=0$ 共 5 个条件，正好对应 $\zeta$ 和 $\psi$ 的 5 个自由度。ZP 对两种场都取了 Lorenz 规范（原文 L333："we impose a Lorenz-gauge condition on both"）。

**推迟解。**
- (10)、(11) 的分布解可以写成"推迟解 + 齐次解"。
- 由假设 (i)，$(h,a)$ 在 $\gamma$ 外光滑；推迟解在 $\gamma$ 外也光滑。所以差值在 $\gamma$ 外光滑，并且跨过 $\gamma$ 也满足齐次方程。
- 由奇异性传播定理（Lorenz 规范下主部是对角的 $\square$），这个差值应该处处光滑【推断，待核：定理对耦合方程组的适用条件】。
- 所以 $(h,a)$ 等于推迟解加上一个光滑齐次解。后者相当于入射辐射，GW 取为零（原文："A choice of retarded solution (corresponding to 'no incoming radiation') … is not necessary"），我们照做。
- **结果**：GHW 直接假设的 $F=F_{\rm ret}+F_{\rm ext}$，在这里推了出来，并且两种场同时成立。

**和 ZP 的关系。**(10)、(11) 就是 ZP 的点粒子方程组：带 $\delta$ 源，取 Lorenz 规范。区别只在来历：ZP 把它当作点粒子模型直接写下，这里是从 (i)–(iii) 推出来的。**到一阶远区为止，ZP 的出发点被证明是对的。**ZP 真正多假设的是下一步：运动方程 (8.2) 中只由正则场驱动。这一步属于 C2。

## 6. C2 要做什么

照 GW 第 V 节：
1. 写出二阶远区扰动 $(j,b)$ 在 $\gamma$ 附近 $1/r^2$ 阶的结构。特解部分来自一阶场的平方：度规里有 RN 的 $q^2/r^2$ 项和 Schwarzschild 的 $M^2$ 项，势里有相应的 $Mq$ 项。剩下的齐次部分 $C_{ab}/r^2$、$C^{[A]}_a/r^2$ 待定。
2. 用近区一阶扰动和"无质量偶极"条件定出 $\gamma$ 的一阶偏离 $Z^i$。
3. **新的问题**：现在有两个偶极，质心和电荷中心。对一般物体它们是两个独立的自由度（GHW 的 $Q_{ab}$）；在 (H) 下，Kerr–Newman 近区背景把电偶极固定为零。要检查"无质量偶极"和"无电偶极"两个条件能否在同一个规范里同时成立。
4. 规范混合：从 Lorenz 规范变到"无偶极"规范时，$\zeta$ 同时改变 $a$。
5. 目标：推出 ZP 式 (8.2)，其中推迟场换成正则场（Detweiler–Whiting），再降阶得到 (8.16)。
6. 检验：$q\to0$ 回到 GW 式 (111)；$G\to0$（不含自身引力）回到 GHW 在弯曲时空的推广，即 DeWitt–Brehme–Hobbs。

独立检验（Ernst 1976 精确解）放在 C2 之后。
