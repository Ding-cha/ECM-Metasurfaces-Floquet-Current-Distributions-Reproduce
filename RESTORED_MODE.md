# Restored 模式说明

本文档说明 `run_ecm.py --mode restored` 的计算流程、关键公式和关键假设。

`restored` 模式的目标不是简单插值作者图，而是在现有论文信息不完整的情况下，尽量恢复一个能自然产生 Fig. 15 主要特征的解析 ECM 候选。

## 运行命令

在 Anaconda Prompt / cmd 中运行：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --anchors
```

如果省略 `--mode restored`，默认也是 restored 模式：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --anchors
```

输出文件：

```text
outputs/s_parameters.png
outputs/s_parameters.csv
outputs/run_config.json
```

## 通用参数输入

当前分支已经把 restored 模式改成通用 symmetric-cross 单元求解入口。默认仍使用论文 Fig. 13 / Table II 参数；如果需要计算任意一组结构参数，可以直接在命令行输入：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --lx-mm 24 --ly-mm 18 --w-mm 4 --px-mm 34 --py-mm 30 --h-mm 1.0 --epsilon-r 2.5 --freq-start-ghz 4 --freq-stop-ghz 10
```

也可以使用 JSON 文件：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --params-json examples\symmetric_cross_case.json
```

JSON 可以只写几何参数，也可以使用 `run_config.json` 这种包含 `geometry` 字段的文件。支持的参数名包括：

```text
lx_mm / l_x / lx
ly_mm / l_y / ly
w_mm / width / w
px_mm / p_x / px
py_mm / p_y / py
h_mm / h
epsilon_r / eps_r
```

程序会检查基本几何约束：

```text
所有参数必须为正数
l_x <= P_x
l_y <= P_y
w <= min(l_x, l_y)
```

如果同时给出 JSON 和命令行参数，命令行参数优先。

## restored 自动基板策略

`restored` 模式现在包含一个自动基板策略：

```text
--restored-profile auto
```

程序用下面的量判断是否属于电厚基板：

```text
tau = h_mm * sqrt(epsilon_r)
```

当：

```text
tau < 4.0
```

使用原先为论文 Fig. 15 调整出的薄基板设置：

```text
mode_order = 9
high_order_substrate = False
width_profile = "cosine"
```

当：

```text
tau >= 4.0
```

使用更适合厚基板或高介电常数基板的设置：

```text
mode_order = 3
high_order_substrate = True
width_profile = "uniform"
epsilon_eff_scale = clamp(1.03 - 0.04 * h_mm, 0.88, 0.97)
arm_length_correction_mm = 0.5
```

原因是厚基板 / 高介电常数基板中，高阶 evanescent Floquet 模态看到的介质层输入导纳不再能简单近似为空气侧导纳；继续使用薄基板 restored 假设会把透射陷波和反射谷推到错误频率。对于 `h = 3 mm, epsilon_r = 4.4` 的 HFSS 对照，完整介质负载会明显更接近仿真中的 4.6 GHz 透射陷波和 5.6-5.8 GHz 反射谷。

厚基板 profile 还默认使用有效介电常数：

```text
epsilon_eff = 1 + epsilon_eff_scale * (epsilon_r - 1)
```

其中厚基板 profile 默认使用随厚度变化的经验值：

```text
epsilon_eff_scale = clamp(1.03 - 0.04 * h_mm, 0.88, 0.97)
```

这个经验修正表示表面电流附近的场并不是完全填充在介质内部。若直接用 bulk `epsilon_r = 4.4` 加载所有模态，程序会把透射陷波算到偏低频；使用有效介电常数后，陷波位置更接近 HFSS 图中的 4.5-4.6 GHz。

`arm_length_correction_mm = 0.5` 只用于电流傅里叶积分：实际输入的几何长度、周期和基板参数不变，但用于电流积分的 `l_x`、`l_y` 会各缩短 0.5 mm。这个修正来自 `h = 2 mm, epsilon_r = 4.4, l_x = l_y = 20..25 mm` 的 HFSS 扫参：不修正时，长度越大，程序的 S21 陷波越偏低频；加入 0.5 mm 有效长度修正后，陷波频率能落回 HFSS 的 0.1 GHz 采样网格附近。

若要关闭这些经验修正，可显式设置：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --restored-profile substrate --epsilon-eff-scale 1 --arm-length-correction-mm 0 --params-json examples\symmetric_cross_case.json
```

也可以手动指定：

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --restored-profile substrate --params-json examples\symmetric_cross_case.json
```

## 几何参数

使用论文 Fig. 13 / Table II 的 symmetric cross MTS unit cell 参数：

| 参数 | 数值 |
| --- | ---: |
| `l_x` | 20 mm |
| `l_y` | 20 mm |
| `w` | 3 mm |
| `P_x` | 30 mm |
| `P_y` | 30 mm |
| `h` | 1.524 mm |
| `epsilon_r` | 3.02 |

代码位置：

```text
src/symmetric_cross_mts/model.py
```

## 电流模型

论文 Eq. (20) 给出 x 极化法向入射下的解析电流分布。坐标单位为 mm。

水平臂电流：

```text
I_1x(x) = 0.95 cos(0.12 x) - 0.2 cos(0.64 x)
```

垂直臂电流：

```text
I_2y(y) = 0.5 I_1x(0) |(|y| - L_y/2) / (L_y/2)|
```

当前代码中垂直电流带有上下方向符号，用于模拟 Fig. 14 中上下臂电流方向相反的情况。

对非 `l_x = 20 mm` 的通用几何，代码把水平臂坐标映射到归一化坐标：

```text
x_ref = x * 20 mm / l_x
I_1x(x; l_x) = 0.95 cos(0.12 x_ref) - 0.2 cos(0.64 x_ref)
```

这个处理保留了论文给出的电流形状，但把电流峰谷位置按实际臂长缩放。它是通用化求解的关键假设之一，因为论文只公开了 `l_x = 20 mm` 情况下的电流拟合式。

对应代码：

```text
src/symmetric_cross_mts/model.py
src/symmetric_cross_mts/fourier.py
```

## Floquet 波数

周期为：

```text
a = P_x
b = P_y
```

Floquet 横向波数：

```text
k_xmn = k0 sin(theta) cos(phi) + 2 m pi / a
k_ymn = k0 sin(theta) sin(phi) + 2 n pi / b
```

法向入射时：

```text
theta = 0
```

所以：

```text
k_xmn = 2 m pi / a
k_ymn = 2 n pi / b
```

纵向波数：

```text
k_zmn^+ = sqrt(k0^2 - k_xmn^2 - k_ymn^2)
k_zmn^- = sqrt(k0^2 epsilon_r - k_xmn^2 - k_ymn^2)
```

`restored` 模式中使用：

```text
kz_branch = "decay"
```

也就是对 evanescent 模态采用衰减分支。这个选择对谐振频率和高阶模态电抗非常敏感。

## 表面电流傅里叶变换

对电流分布做二维傅里叶变换：

```text
J_tilde_i(k_xmn, k_ymn)
= 1/(4 pi^2) ∫∫ J_i(x, y) exp(j k_xmn x) exp(j k_ymn y) dx dy
```

其中 `i = x, y`。

代码中利用 cross 的条带结构，把二维积分拆成一维数值积分和宽度方向解析积分。

对应代码：

```text
src/symmetric_cross_mts/fourier.py
```

## TE/TM 投影

论文 Eq. (9) 中的矢量电流投影为：

```text
FT_mn^TE = J_tilde_x k_ymn - J_tilde_y k_xmn
FT_mn^TM = J_tilde_x k_xmn + J_tilde_y k_ymn
```

原始 `raw` 模式直接采用论文排版形式的 TTR 公式。但探索发现，这样会使金属结构响应过弱，结果接近裸介质板。

`restored` 模式采用单位极化方向归一化后的投影形式：

```text
T_mn^TE = C (FT_mn^TE)^2 / (k_xmn^2 + k_ymn^2)
T_mn^TM = C (FT_mn^TM)^2 / (k_xmn^2 + k_ymn^2)
```

其中：

```text
C = (4 pi^2)^2 / (a b)
```

这个设置在代码中对应：

```text
ttr_variant = "polarized"
```

## 模态导纳

空气侧导纳：

```text
Y_mn^TE+ = k_zmn^+ / (omega mu0)
Y_mn^TM+ = omega epsilon0 / k_zmn^+
```

介质侧导纳：

```text
Y_mn^TE- = k_zmn^- / (omega mu0)
Y_mn^TM- = omega epsilon0 epsilon_r / k_zmn^-
```

基模 `m = 0, n = 0` 的介质层输入导纳仍按照论文 Eq. (8) 的传输线形式计算。

## 高阶模态介质负载假设

这是 restored 模式最关键的恢复假设。

`raw` 模式对所有 Floquet 模态都使用介质层传输线输入导纳：

```text
B_mn = dielectric slab input admittance
```

但这样会把主谐振推到约 4.95 GHz，与作者 Fig. 15 的约 5.8 GHz 不符。

探索发现，如果只对基模保留介质层负载，而对高阶 evanescent 模态使用空气侧导纳近似：

```text
for m,n != 0:
    B_mn^TE ≈ Y_mn^TE+
    B_mn^TM ≈ Y_mn^TM+
```

主谐振会回到 5.8 GHz 附近，并自然产生作者图中的反射峰和透射陷波。

代码中对应：

```text
high_order_substrate = False
```

这不是论文明确写出的公式，而是根据 Fig. 15 反推得到的恢复性假设。

## 等效电路求和

对所有截断范围内的 Floquet 模态求和：

```text
Z_AB = Σ T_mn^TE / y_mn^TE + Σ T_mn^TM / y_mn^TM
```

其中激励模态 `TE00` 从高阶求和中排除，用于端口输入阻抗计算。

`restored` 模式参数：

```text
mode_order = 9
fourier_samples = 501
```

也就是：

```text
m, n = -9, ..., +9
```

## S11 计算

按照论文 Eq. (12)-(13)：

```text
Z_in,00^TE = [ (Z_AB / T_00^TE)^(-1) + B_00^TE ]^(-1)
```

反射系数：

```text
S_11^TE00 = (Z_in,00^TE - Z_00^TE+) / (Z_in,00^TE + Z_00^TE+)
```

其中：

```text
Z_00^TE+ = 1 / Y_00^TE+
```

## S21 计算

当前代码使用论文 Eq. (14) 的传输线电压分配形式：

```text
S_21^TE00 =
2 Z_in,00^TE / (Z_in,00^TE + Z_00^TE+)
* (1 + Gamma) / (exp(j beta_00 h) + Gamma exp(-j beta_00 h))
```

其中：

```text
Gamma = (Z_load - Z_c) / (Z_load + Z_c)
```

这部分仍是近似实现，但 restored 模式已经能给出与 Fig. 15 接近的主趋势。

## 当前 restored 模式结果

典型输出：

```text
5.80 GHz: S11 ≈ -0.03 dB, S21 ≈ -21.46 dB
5.85 GHz: S11 ≈ -0.0006 dB, S21 ≈ -38.33 dB
7.80 GHz: S11 ≈ -18.77 dB, S21 ≈ -0.058 dB
8.27 GHz: S11 ≈ -23.01 dB, S21 ≈ -0.022 dB
```

这些结果能复现作者图中的主要特征：

- 5.8 GHz 附近 `S11` 接近 0 dB；
- 5.8 GHz 附近 `S21` 出现深陷波；
- 8 GHz 附近 `S11` 出现反射谷；
- 高频段 `S21` 回到接近 0 dB。

## 与其他模式的区别

### raw

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode raw
```

直接按论文公式形式实现。结果差异很大，主要用于诊断。

### calibrated

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode calibrated --anchors
```

使用 Fig. 15 手工读数锚点插值，最接近论文图，但不是独立的 ECM 计算。

### restored

```cmd
set PYTHONPATH=src && D:\anaconda3\python.exe -m symmetric_cross_mts.run_ecm --mode restored --anchors
```

当前最接近作者结果的解析恢复候选。它不直接插值作者图，但包含恢复性假设。

## 限制

restored 模式仍不能保证完全等同作者原始 ECM，因为论文没有公开以下细节：

- `I0` 的具体归一化；
- 法向入射 `TE00` 极限的完整推导；
- 作者程序中高阶 evanescent 模态的精确分支和负载处理；
- 金属厚度、损耗、介质损耗等仿真细节；
- 作者用于 Fig. 15 的模态截断和数值积分设置。

因此 restored 模式应理解为：

```text
基于论文公式和 Fig. 15 反推得到的最接近解析 ECM 候选
```

而不是作者代码的唯一确定实现。
