# 2 假设

> 草稿 · 2026-11-24。依据：`derivations/p1/c0-assumptions.md`；Gralla 2010 第 II、IV 节；GW 2008 第 II–III 节；GHW 2009 第 II 节。
>
> **定稿前要核对的地方：**
> 1. 2.5 节唯一性定理的条件，目前根据 Lopes Costa 2010 的摘要和 Chruściel–Lopes Costa–Heusler 2012 的综述概括，还没有逐条对照原文的定理表述。
> 2. 2.3 节对电磁规范变换的光滑性要求是我们补的，Gralla 2010 没有讨论；它是否足够用，要到一阶推导（第 4 节）里检验。
> 3. 2.5 节"四极及更高的多极矩不进入一阶运动"是量级估计，不是证明。

---

## 2.1 约定与场方程

号差取 $(-,+,+,+)$，单位取 $G=c=1$，电磁量用高斯单位。电磁场 $F=dA$，能量–动量张量为
$$4\pi T_{ab}=F_{ac}F_b{}^c-\tfrac14g_{ab}F_{cd}F^{cd}. \tag{2.1}$$
Einstein–Maxwell 方程记为
$$E_{ab}\equiv G_{ab}-8\pi T_{ab}=0,\qquad \mathcal M^a\equiv\nabla_bF^{ab}=0. \tag{2.2}$$
有源时 $\mathcal M^a=4\pi j^a$。在这一约定下，Coulomb 势 $A=-(q/r)\,dt$ 给出 $\mathcal M^t=4\pi q\,\delta^3$，电荷 $q$、质量 $m$ 的检验粒子满足 $m\,u^b\nabla_bu^a=q\,F^a{}_bu^b$，同号电荷相斥。

## 2.2 一族解与标度

我们考虑 Einstein–Maxwell 方程的一族解 $(g_{ab}(\lambda),A_a(\lambda))$，其中参数 $\lambda\ge0$ 同时控制物体的尺寸、质量和电荷：物体的质量为 $\lambda M$，电荷为 $\lambda q$，荷质比 $\xi=q/M$ 在整族解中保持不变。$\lambda\to0$ 时物体缩成一条世界线，留下的是一个不含物体的背景。

确定近区的标度要用到 Einstein–Maxwell 方程的一个对称性：对任意常数 $c$，变换
$$(g_{ab},A_a)\ \to\ (c^2g_{ab},\,cA_a) \tag{2.3}$$
把解映为解。原因是 $G_{ab}$ 在常数共形变换下不变；$T_{ab}$（指标全在下）中的 $F_{ac}F_{bd}g^{cd}$ 和 $g_{ab}F^2$ 都不随 $c$ 改变；$\nabla_bF^{ab}$ 整体乘以 $c^{-1}$，方程本身不变。Gralla（2010）正是据此选取电磁势的标度。对每个时刻 $t_0$，定义近区坐标 $\bar t=(t-t_0)/\lambda$、$\bar x^i=x^i/\lambda$，以及近区的场
$$\bar g_{ab}=\lambda^{-2}g_{ab},\qquad \bar A_a=\lambda^{-1}A_a. \tag{2.4}$$
在近区坐标里，分量之间的关系就是把 $t=t_0+\lambda\bar t$、$x^i=\lambda\bar x^i$ 直接代入：
$$\bar g_{\bar\mu\bar\nu}(\lambda;t_0;\bar t,\bar x^i)=g_{\mu\nu}(\lambda;t_0+\lambda\bar t,\lambda\bar x^i),\qquad
\bar A_{\bar\mu}(\lambda;t_0;\bar t,\bar x^i)=A_\mu(\lambda;t_0+\lambda\bar t,\lambda\bar x^i). \tag{2.5}$$
由 (2.3)，$(\bar g,\bar A)$ 仍然是 Einstein–Maxwell 方程的解。在近区里，物体的尺寸、质量和电荷都是 $O(1)$ 的量。

## 2.3 三条极限假设

我们沿用 GW 的三条假设，按 Gralla（2010）第 IV 节的写法，把其中的度规换成一对场 $(g_{\mu\nu},A_\mu)$。记 $r=|x|$。我们要求：对每个 $\lambda\ge0$，存在坐标 $(t,x^i)$ 和一个电磁规范，使下面三条成立。

**(i) 普通极限。**存在常数 $\bar R>0$，使得在 $r\ge\lambda\bar R$ 上，$(g_{\mu\nu},A_\mu)$ 满足 (2.2)，并且对 $(\lambda,t,x^i)$ 联合光滑。$\lambda=0$ 时，$(g_{\mu\nu},A_\mu)$ 在包括 $r=0$ 在内的区域上处处光滑，曲线 $\gamma$：$r=0$ 是类时的。

**(ii) 标度极限。**在 $\bar r\ge\bar R$ 上，$(\bar g_{\bar\mu\bar\nu},\bar A_{\bar\mu})$ 对 $(\lambda,t_0,\bar t,\bar x^i)$ 联合光滑，包括 $\lambda=0$。

**(iii) 一致性。**记 $\alpha=r$、$\beta=\lambda/r$。在固定的 $(t,\theta,\phi)$ 下，$(g_{\mu\nu},A_\mu)$ 作为 $(\alpha,\beta)$ 的函数，在 $(\alpha,\beta)=(0,0)$ 处联合光滑。GW 后来补注说，推导实际需要的是对全部变量联合光滑，我们按这个较强的形式使用。

对 $r<\lambda\bar R$ 的区域不做任何假设。按 GW 的说法，这一区域"could equally well be 'filled in' with ordinary matter or a black hole"。只有到 2.5 节加上附加假设 (H)，物体才被指定为黑洞。

三条假设是用势 $A_\mu$ 而不是场强 $F_{\mu\nu}$ 表述的，这带来两点需要说明。第一，$A_\mu$ 的分量依赖规范，所以假设应读作"存在一个规范，使 (i)–(iii) 成立"。这和 GW 中"存在坐标"的地位相同。允许的规范变换 $A\to A+d\psi$ 中，$\psi(\lambda,x)$ 应满足与 GW 式 (21) 中坐标变换同类的光滑性要求，使得 (i)–(iii) 在变换下保持。Gralla（2010）没有讨论这一点，这一要求是我们补上的。第二，区域 $r\ge\lambda\bar R$ 含有包住物体的二维球面，在这样的区域上，闭形式 $F$ 能写成 $dA$ 的条件是 $F$ 在该球面上的积分为零，也就是物体没有磁荷。因此用势表述的假设隐含了零磁荷，我们把它明确列为假设的一部分。采用势而不是场强，是因为一阶推导中要处理两种规范的混合：无穷小微分同胚 $\zeta$ 会使 $A_a$ 改变 $\zeta^cF_{ca}$，这一步离不开势。

## 2.4 直接推论

这一节的结果对任意满足 (i)–(iii) 的物体成立，是 GW 和 Gralla（2010）相应结果在 Einstein–Maxwell 下的直接推广。由 (iii)，把 $(g_{\mu\nu},A_\mu)$ 在 $\beta$ 上展开到一阶，得
$$g_{\mu\nu}=b_{\mu\nu}(t)+O(r)+\lambda\Big[\frac{c_{\mu\nu}(t,\theta,\phi)}{r}+O(1)\Big]+O(\lambda^2), \tag{2.6}$$
$A_\mu$ 的展开形式相同，系数记为 $b^{[A]}_\mu$、$c^{[A]}_\mu$。从这个展开读出三个结果。

第一，远区背景 $(g,A)|_{\lambda=0}$ 是一个在 $\gamma$ 上光滑的电真空解，其中的电磁场 $F$ 一般不为零，这就是物体所处的外场。第二，一阶远区扰动 $(h_{ab},a_a)\equiv\partial_\lambda(g_{ab},A_a)|_{\lambda=0}$ 在 $\gamma$ 附近最多像 $1/r$ 那样奇异。第三，近区背景
$$\bar g^{(0)}_{\bar\mu\bar\nu}=b_{\mu\nu}(t_0)+\frac{c_{\mu\nu}(t_0,\theta,\phi)}{\bar r}+O(\bar r^{-2}),\qquad
\bar A^{(0)}_{\bar\mu}=b^{[A]}_\mu(t_0)+\frac{c^{[A]}_\mu(t_0,\theta,\phi)}{\bar r}+O(\bar r^{-2}) \tag{2.7}$$
不依赖 $\bar t$，在 $\bar r\to\infty$ 时趋于常数，并且满足 Einstein–Maxwell 方程。也就是说，近区背景是一个稳态、渐近平直的电真空解。稳态性是推出来的，不是假设的：它来自 $g_{\mu\nu}(\lambda)$ 对 $t$ 的光滑性（Gralla 2010；GW 称之为"derived here rather than assumed"）。常数 $b^{[A]}_\mu(t_0)$ 是纯规范。

外场在近区里出现的阶次可以直接读出。由 (2.5)，$\bar F_{\bar\mu\bar\nu}=\lambda F_{\mu\nu}$，所以远区背景的电磁场在近区里是一个 $O(\lambda)$ 的均匀场；它的梯度和背景曲率在近区里都是 $O(\lambda^2)$。这和 GW 中潮汐场出现在近区二阶是同样的结构。一个推论是，外场在物体上感应出的偶极大小约为 $(\lambda M)^3F$，它产生的力比洛伦兹力高 $\lambda^2$ 阶，不进入一阶运动。

## 2.5 黑洞的附加假设

(i)–(iii) 不要求物体是黑洞。为了让研究对象是带电黑洞，我们加一条假设：

**(H)** 近区背景 $(\bar g^{(0)},\bar A^{(0)})$ 是一个 Kerr–Newman 解，参数为质量 $M$、电荷 $q$、比角动量 $\bar a$，并且严格次极端：$\bar a^2+q^2<M^2$。

这条假设的动机来自电真空黑洞的唯一性定理。Lopes Costa（2010）证明，在稳态、解析、$\mathcal I^+$ 正则、视界连通且非简并的电真空时空中，黑洞只能是 Kerr–Newman 解；这一结果建立在 Mazur（1982）和 Bunting（1983）的工作之上，综述见 Chruściel、Lopes Costa 和 Heusler（2012）。2.4 节已经推出近区背景是稳态、渐近平直的电真空解，但解析性、视界的存在和正则性、视界的连通性和非简并性，都不能从 (i)–(iii) 推出。因此 (H) 必须作为假设写出，不能说是推论。视界连通的条件不能省：Majumdar–Papapetrou 解给出了多个极端黑洞静态共存的例子。严格次极端的要求对应于视界非简并；极端情形留给后续的工作。

(H) 在一阶运动中固定了什么，可以用量级来看。近区里黑洞的自旋是 $S=\lambda^2M\bar a$，电偶极和磁偶极也都是 $O(\lambda^2)$。它们产生的力（Papapetrou 自旋力，以及偶极在外场梯度中受的力）相对于洛伦兹力是 $O(\lambda)$，所以进入一阶运动。GW 式 (111) 的第一项就是 Papapetrou 自旋力，GHW 式 (134)–(136) 中出现的是电磁偶极和自旋。对一般物体，这些偶极矩是自由参数；对 Kerr–Newman 解，它们由 $(M,q,\bar a)$ 完全固定：相对于质心的电偶极为零，磁偶极为 $\mu=q\bar a$，即旋磁比 $g=2$（Carter 1968）。所以在偶极这一层，带电黑洞就是偶极矩由 $g=2$ 固定的 GHW 物体。四极及更高的多极矩在近区里至少是 $O(\lambda^3)$，按同样的估计，它们的力相对于洛伦兹力至少是 $O(\lambda^2)$，不进入一阶运动。黑洞的结构是否还会通过近区的一阶响应带来别的项，要在第 4 节的推导中检验。与 ZP 的无自旋点粒子对照时取 $\bar a=0$，近区背景是 Reissner–Nordström 解，所有偶极项都消失。

## 2.6 与已有框架的比较

和 GW 相比，这里的场从真空中的度规换成了耦合的 $(g,A)$，物体同样可以是任意的，包括黑洞，物体的自身引力同样保留。和 GHW 相比，GHW 处理的是平直时空中的物质物体，不含自身引力，并且直接假设电磁场是物体电流的推迟解加上一个随 $\lambda$ 光滑的外场；这里不做这一假设，推迟解的结构要在一阶远区里推出来（第 3 节）。和 Gralla（2010）相比，三条假设的写法相同，我们补充了三点：电磁规范变换的光滑性要求，零磁荷条件，以及指定物体为黑洞的附加假设 (H)。Gralla（2010）只推到零阶，本文要把同一组假设用到一阶。

最后要指出，(i)–(iii) 是对固定的 $\xi$ 提出的。常数 $\bar R$，以及 (ii)、(iii) 中"光滑"背后隐含的导数界，都可能依赖 $\xi$。$\xi\to1$ 时它们能否一致地选取，决定了推导在近极端情形是否仍然成立，这一点在第 6 节讨论。
