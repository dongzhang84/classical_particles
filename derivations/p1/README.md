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
| [`calc.tex`](calc.tex) → [`calc.pdf`](calc.pdf) | C2 验证计算工作文档：一阶运动分四级推进（见下表），每级记录输入、中间式子、得到的关系、代码、是否通过 | 第 1、2、3 级完成 |
| [`code/`](code/) | sympy 代码。`series.py` 是共用工具（线性化 Ricci、按阶次取系数、ℓ = 1 拟合）；每级一个脚本，另有检验脚本 | 见下 |
| C3 | 独立检验：Ernst 1976 精确解按 $m$ 展开 | 未开始，先查是否已有人做过 |

## C2：一阶运动的四级计算

| 级 | 情形 | 对照的答案 | 状态 |
|---|---|---|---|
| 1 | q = 0，F = 0 | GW 式 (106)–(108) | 通过（2026-11-24）。顺带发现 GW 式 (106) 的 R^A、R^F 两行有误，不影响运动方程 |
| 2 | q ≠ 0，F = 0 | MiSaTaQuWa + DeWitt–Brehme 尾巴 | 通过（2026-11-24）。运动方程恰好多出 q f^tail_{i0}，4 个规范模式全部通过 |
| 3 | q = 0，F ≠ 0 | ZP 取 q → 0：(5/3)M²FF 局部项 + 耦合尾巴 | 完成（2026-11-25）。决定运动的块 I 与第 1 级相同，运动由正则场梯度决定；代入 ZP 的正则场局部项得 (5/3)M²(E×B)。块 II 的自洽性和局部项的独立推导未做 |
| 4 | 完整 | ZP 式 (8.2) 的线性化 | 待做 |

## 约定

号差 $(-,+,+,+)$，$G=c=1$，Gauss 单位。$4\pi T_{ab}=F_{ac}F_b{}^c-\tfrac14g_{ab}F^2$，$\nabla_bF^{ab}=4\pi j^a$，洛伦兹力为 $m\,u^b\nabla_bu^a=qF^a{}_bu^b$。

## 运行检验

```sh
python3 derivations/p1/check_c1_identities.py   # 约 10 秒，三行都应为 0

# C2 代码要用 sympy，装在项目内的 .venv（已 gitignore）
python3 -m venv .venv && .venv/bin/pip install sympy
cd derivations/p1/code
../../../.venv/bin/python level1.py       # 第 1 级，约 2.5 分钟
../../../.venv/bin/python check_dc.py     # 检验 (a)：不截断的完整 Ricci 计算
../../../.venv/bin/python gauge_test.py   # 检验 (b)：三个纯规范扰动
../../../.venv/bin/python level2.py       # 第 2 级，约 19 分钟，输出 level2_matrix.txt
../../../.venv/bin/python gauge_test2.py  # 第 2 级的四个规范检验
../../../.venv/bin/python level3_detail.py  # 第 3 级：近区一阶场、匹配、T^(2) 的结构，约 3 秒
```

工作文档用 xelatex 编译：`cd derivations/p1 && latexmk calc`。
