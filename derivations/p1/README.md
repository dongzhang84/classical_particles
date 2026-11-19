# P1 推导：带电黑洞一阶运动方程的严格推导

主线问题见 [`docs/problem.md`](../../docs/problem.md) 的 P1。走 Gralla–Wald 路线：以 Gralla 2010 的 Einstein–Maxwell 框架为起点，从零阶推广到一阶。

## 分工

ZP 2014 已经给出耦合推迟 Green 函数的 Hadamard 展开、奇异场和尾巴项，这些直接沿用。我们要推的只有一步：**从场方程推出"驱动运动的只是正则场"**。ZP 在这一步写的是 "on a formal level"，LFW 写的是 "assume"。

## 文件

| 文件 | 内容 | 状态 |
|---|---|---|
| [`c0-assumptions.md`](c0-assumptions.md) | 检查点 1：三条假设在 Einstein–Maxwell 下的写法；黑洞附加假设 (H)；为 P3 留的记号 | 草稿 |
| [`c1-zeroth-order.md`](c1-zeroth-order.md) | 零阶洛伦兹力（复现 Gralla 2010、2013）；一阶远区耦合方程组；推迟解的来历 | 草稿 |
| [`check_c1_identities.py`](check_c1_identities.py) | C1 式 (1)、(3)、(4) 的数值检验（任意度规、精确有理数） | 三项通过 |
| C2 | 一阶运动：两个偶极条件、规范混合，目标是 ZP (8.16) | 未开始 |
| C3 | 独立检验：Ernst 1976 精确解按 $m$ 展开 | 未开始，先查是否已有人做过 |

## 约定

号差 $(-,+,+,+)$，$G=c=1$，Gauss 单位。$4\pi T_{ab}=F_{ac}F_b{}^c-\tfrac14g_{ab}F^2$，$\nabla_bF^{ab}=4\pi j^a$，洛伦兹力为 $m\,u^b\nabla_bu^a=qF^a{}_bu^b$。

## 运行检验

```sh
python3 derivations/p1/check_c1_identities.py   # 约 10 秒，三行都应为 0
```
