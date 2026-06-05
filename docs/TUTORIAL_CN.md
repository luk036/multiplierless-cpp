# 🎛️ 无乘法器 FIR 滤波器设计 — 教程

## 从规格到硅片：JSON → CSD → Verilog

本教程将引导你使用 `FirDesign` 工具链设计无乘法器 FIR 低通滤波器。
从 JSON 规格文件开始，一直到可综合的 Verilog 模块 — 可直接用于 ASIC 或 FPGA。

---

## 📋 目录

1. [工具功能概述](#1-工具功能概述)
2. [安装](#2-安装)
3. [编写滤波器规格](#3-编写滤波器规格)
4. [运行 CLI](#4-运行-cli)
5. [理解输出结果](#5-理解输出结果)
6. [Verilog 模块](#6-verilog-模块)
7. [使用 Yosys 综合](#7-使用-yosys-综合)
8. [完整示例演练](#8-完整示例演练)
9. [验证滤波器是否满足规格](#9-验证滤波器是否满足规格)
10. [调优指南](#10-调优指南)

---

## 1. 工具功能概述

```
   JSON 规格  →  椭球算法优化  →  CSD 量化  →  Verilog RTL
```

- **输入**：描述滤波器需求的 JSON 文件（通带、阻带、纹波、衰减、抽头数、CSD 精度）。
- **引擎**：椭球算法求解凸优化问题 `min max R(ω)`，满足通带/阻带约束。
- **量化**：每个系数转换为规范符号数字（CSD）格式 — 仅使用 `+1`、`-1`、`0` — 使硬件乘法变为移位和加减法。
- **输出**：包含 CSD 系数和可综合 Verilog 模块的 JSON 文件，自动进行交叉公共子表达式消除（Cross-CSE）。

### 为什么需要无乘法器？

| 方案 | 每抽头硬件 | 32 抽头成本 |
|------|-----------|------------|
| 通用乘法器 | 1 个乘法器 | 32 个乘法器 |
| CSD (nnz=4) | 3 个加法器/减法器 | 96 个加法器 |
| CSD + Cross-CSE | 2–3 个加法器（共享） | ~60 个加法器 |

---

## 2. 安装

### 方案 A：C++ (cmake)

```bash
cd multiplierless-cpp
cmake -S. -Bbuild
cmake --build build --config Release
```

CLI 二进制文件位于 `build\standalone\Release\FirDesign.exe`
（Linux 上为 `build/standalone/FirDesign`）。

### 方案 B：Python (pip)

```bash
cd ../py/multiplierless
pip install -e .
```

CLI 命令为 `fir-design`。

### 方案 C：Python (模块)

```bash
conda activate pyflow  # 或你的环境
python -m multiplierless.fir_design spec.json
```

---

## 3. 编写滤波器规格

创建一个 JSON 文件，例如 `my_filter.json`：

```json
{
  "filter_order": 48,
  "passband_edge": 0.12,
  "stopband_edge": 0.20,
  "passband_ripple": 0.025,
  "stopband_attenuation": 0.125,
  "csd_nnz": 6,
  "discretization_factor": 15,
  "max_iters": 50000,
  "tolerance": 1e-14,
  "ellipsoid_radius": 40.0,
  "parallel_cut": true,
  "verilog": {
    "input_width": 16,
    "module_name": "my_fir_filter"
  }
}
```

### 必需参数

| 参数 | 类型 | 描述 |
|-----|------|------|
| `filter_order` | int (≥4) | FIR 抽头数 |
| `csd_nnz` | int (1–16) | 每个系数的最大非零 CSD 位数 |

所有其他参数均有合理默认值（如上所示）。

### 关键参数说明

| 参数 | 默认值 | 含义 |
|-----|--------|------|
| `passband_edge` | `0.12` | ×π 弧度。越小通带越窄 |
| `stopband_edge` | `0.20` | ×π 弧度。越大过渡带越宽 |
| `passband_ripple` | `0.125` | 线性纹波。0.125 ≈ ±1.0 dB |
| `stopband_attenuation` | `0.125` | 线性衰减。0.125 ≈ −18 dB |
| `discretization_factor` | `15` | 频率采样点数 = N × 此值。越大越精确但越慢 |
| `csd_nnz` | *(必需)* | 硬件成本 = nnz − 1 个加法器/抽头 |
| `max_iters` | `50000` | 椭球算法迭代上限 |
| `tolerance` | `1e-14` | 收敛容差。越小越严格 |
| `ellipsoid_radius` | `40.0` | 初始搜索区域大小 |
| `parallel_cut` | `true` | 启用并行切割以加速收敛 |
| `spectral_method` | `"root"` | `"root"`（Aberth，更快）或 `"fft"`（Kolmogorov，传统） |
| `root_tolerance` | `1e-8` | Aberth 收敛容差。越大越宽松但更快 |

### 谱分解方法

CLI 支持两种谱分解算法：

| 方法 | 键值 | 速度 | 精度 | 适用场景 |
|------|------|------|------|----------|
| **Aberth 根求解** | `"root"` | 快 (~333 次迭代) | 良好（可调容差） | 默认，生产环境 |
| **FFT (Kolmogorov 1939)** | `"fft"` | 较慢 (~1482 次迭代) | 参考值（精确） | 验证、对比 |

根求解方法使用 Aberth-Ehrlich 算法通过多项式根求解直接提取最小相位根。
容差 (`root_tolerance`) 控制收敛：
- `1e-4`：宽松但快 — 适合快速实验
- `1e-8`：**默认值** — 速度与精度的平衡
- `1e-12`：严格 — 适用于窄过渡带或高阶滤波器

> **注意**：对于大型滤波器（N > 64），将 `root_tolerance` 增大到 `1e-6`
> 或切换到 `"fft"` 以避免收敛问题。

### 获取 Verilog 输出

添加 `verilog` 配置块：

```json
"verilog": {
  "input_width": 16,
  "module_name": "fir_filter"
}
```

完全省略 `verilog` 键则只输出系数（更快）。

---

## 4. 运行 CLI

### C++

```bash
./build/windows/x64/release/FirDesign.exe my_filter.json | tee output.json
```

### Python

```bash
fir-design my_filter.json | tee output.json
```

JSON 输出打印到标准输出。重定向保存：

```bash
fir-design my_filter.json > output.json
```

**重要提示**：优化可能耗时 1–5 分钟，取决于 `filter_order` 和 `max_iters`。
初次实验建议使用 `{"filter_order": 16, "csd_nnz": 4}`，可在几秒内收敛。

---

## 5. 理解输出结果

```json
{
  "filter_order": 48,
  "csd_nnz": 6,
  "iterations": 2341,
  "coefficients": [
    {"index": 0,  "value": 0.008912, "csd": "0.000000+00+000-000-0+0-"},
    {"index": 1,  "value": 0.017234, "csd": "0.00000+000+0+00-0+00-"},
    ...
    {"index": 47, "value": 0.002145, "csd": "0.00000000+000-0+0+"}
  ],
  "verilog": "module fir_filter( ..."
}
```

### 输出字段

| 字段 | 含义 |
|------|------|
| `iterations` | 椭球算法迭代次数。越少收敛越快 |
| `coefficients[].value` | CSD 量化后的系数值（double）。用于仿真 |
| `coefficients[].csd` | CSD 字符串。`+` = +1，`-` = −1，`0` = 0，`.` = 小数点 |
| `verilog` | 完整可综合 Verilog 模块（如果请求） |

### 如何读懂 CSD 字符串

```
0.000+00-0+0+
│ │   │  │ │ │
│ │   +  - + +  ← 非零位（使用了 4 位，nnz=6）
│ └────────────  小数部分（幂次：2⁻¹, 2⁻², ...）
└──────────────  整数部分（幂次：2⁰, 2¹, 2², ...）
```

每个非零位在硬件中对应一个加法器/减法器。
每抽头加法器总数最多 = `nnz − 1`。

### 滤波器对称性

脉冲响应 `h` 是**最小相位**的（由谱分解得出）。
对于线性相位 FIR 滤波器，直接使用平方幅频响应（或自相关）。

---

## 6. Verilog 模块

生成的 Verilog **仅使用移位和加减法** — 没有 `*`。

### 模块接口

```verilog
module fir_filter (
    input signed [15:0] x,       // 16 位输入采样
    output signed [41:0] h0,     // 系数 h[0]（41 位 = 16 + 25）
    output signed [41:0] h1,
    ...
    output signed [41:0] h47
);
```

### 内部结构

```verilog
// 输入信号的移位副本（每个唯一移位量一条线）
wire signed [41:0] x_shift17 = x <<< 17;
wire signed [41:0] x_shift16 = x <<< 16;
wire signed [41:0] x_shift14 = x <<< 14;
...

// Cross-CSE：共享子表达式
// 模式 "+0-" 出现在 h0、h5、h13、h27 中 → 提取为公共项！
wire signed [41:0] _cse_0 = x_shift14 - x_shift12;

// 每系数赋值（仅移位和加减）
wire signed [41:0] h0 = (_cse_0 >>> 6) + x_shift4 + x_shift0;
wire signed [41:0] h5 = (_cse_0 >>> 5) + x_shift3 + x_shift0;
...
```

### Cross-CSE 说明

当多个系数共享相同的 CSD 模式（如 `+0-`）时，该模式仅计算一次作为共享线 `_cse_0`，
并通过右移（`>>>`）在抽头间复用。这**减少了整个滤波器的加法器总数**。

```
无 CSE：  h0 的 +0- 需要 1 个加法器，h5 的 +0- 需要 1 个加法器 = 2 个加法器
有 CSE：  _cse_0 = +0-（1 个加法器），h0 = _cse_0 >>> 6，h5 = _cse_0 >>> 5 = 1 个加法器
```

### 输出位宽

输出位宽 = `input_width + max_power`，其中 `max_power` 是所有系数中最高的移位量。

| input_width | 典型 output_width |
|-------------|------------------|
| 8 | ~25 |
| 16 | ~41 |
| 24 | ~55 |

---

## 7. 使用 Yosys 综合

Verilog 输出可直接用于 Yosys 综合：

```bash
# 综合为通用门级网表
yosys -p "
  read_verilog fir_filter_cpp.v;
  synth -top fir_filter;
  write_json fir_filter_synth.json;
  stat
"
```

针对特定工艺（如 SkyWater 130nm）：

```bash
yosys -p "
  read_verilog fir_filter_cpp.v;
  synth -top fir_filter;
  abc -liberty /path/to/sky130_fd_sc_hd__tt_025C_1v80.lib;
  write_verilog fir_filter_mapped.v;
  stat
"
```

### 预期面积

32 抽头、nnz=4、启用 Cross-CSE 的滤波器，Yosys 实际综合结果如下：

```
=== fir_filter ===

   Number of wires:                100
   Number of wire bits:           3778
   Number of public wires:          49
   Number of cells:                 82
     $add                           51
     $neg                            7
     $sub                           24

   无乘法器 — 纯移位加减！
```

无 Cross-CSE：约 110 个单元。有 Cross-CSE：**82 个单元**（共享子表达式减少了 25%）。
每个 `$add` 或 `$sub` 对应一个硬件加法器/减法器。

使用提供的脚本运行综合：

```bash
# 修复 Verilog 端口语法（CSD 生成器产生无逗号的端口）
python tools/fix_verilog_ports.py fir_filter_tutorial.v

# 使用 Yosys 综合
yosys -Q synthesize.ys
```

`synthesize.ys` 脚本执行以下操作：
```tcl
read_verilog fir_filter_tutorial.v    # 加载设计
hierarchy -top fir_filter             # 检查模块层次
proc                                  # 处理 always 块
opt                                   # 优化表达式
stat                                  # 打印面积/加法器计数
write_json fir_filter_synth.json      # JSON 网表用于分析
write_verilog -noattr fir_filter_synth.v  # 综合后 Verilog
```

---

## 8. 完整示例演练

设计一个 32 抽头低通滤波器，每个系数最多 4 个非零 CSD 位。

### 第 1 步：编写规格

保存为 `lowpass_32.json`：

```json
{
  "filter_order": 32,
  "csd_nnz": 4,
  "passband_edge": 0.12,
  "stopband_edge": 0.20,
  "passband_ripple": 0.125,
  "stopband_attenuation": 0.125,
  "verilog": {
    "input_width": 16,
    "module_name": "fir_filter"
  }
}
```

### 第 2 步：运行工具

```bash
fir-design lowpass_32.json | tee lowpass_32_output.json
```

### 第 3 步：检查结果

查看输出：

```json
{
  "iterations": 1482,
  "filter_order": 32,
  "csd_nnz": 4,
  "coefficients": [
    {"index": 0,  "value": 0.008724, "csd": "0.000000+00+000-000-"},
    {"index": 7,  "value": 0.128059, "csd": "0.00+0000+0-00+"},
    {"index": 15, "value": -0.017375,"csd": "0.00000-00-00+0-"},
    ...
  ]
}
```

- **1482 次迭代** — 2 分钟内收敛
- **系数 h[7]** = 0.128059，CSD = `0.00+0000+0-00+` — 4 个非零位
- **系数 h[15]** = −0.017375，CSD = `0.00000-00-00+0-` — 4 个非零位

### 第 4 步：提取 Verilog

`verilog` 字段包含完整模块。提取为文件：

```bash
python -c "
import json, sys
data = json.load(open('lowpass_32_output.json'))
with open('fir_filter.v', 'w') as f:
    f.write(data['verilog'])
print(f'{len(data[\"verilog\"])} chars written')
"
```

或使用内置工具：

```bash
python tools/extract_verilog.py
```

### 第 5 步：修复 Verilog 语法并综合

CSD 生成器输出的端口之间缺少逗号。修复：

```bash
python tools/fix_verilog_ports.py fir_filter_tutorial.v
```

然后综合：

```bash
yosys -Q synthesize.ys
```

预期输出：
```
=== fir_filter ===
   Number of cells:   82
     $add             51
     $neg              7
     $sub             24
```

### 第 6 步：集成到你的设计中

在顶层模块中实例化：

```verilog
wire signed [15:0] sample_in;
wire signed [41:0] tap_0, tap_1, ..., tap_31;

fir_filter u_fir (
    .x(sample_in),
    .h0(tap_0), .h1(tap_1), ..., .h31(tap_31)
);

// 带适当延迟求和以得到滤波器输出
// （或输入到移位寄存器 + 加法树）
```

---

## 9. 验证滤波器是否满足规格

生成系数后，在流片前验证滤波器确实满足通带/阻带要求。

### 9.1 频率响应分析

使用输出 JSON 在 Python 中计算频率响应：

```python
import json
import numpy as np

# 加载设计输出
with open("lowpass_32_output.json") as f:
    data = json.load(f)

h = np.array([c["value"] for c in data["coefficients"]])
N = len(h)

# 计算频率响应
w = np.linspace(0, np.pi, 4096)  # 密集网格
H = np.zeros(len(w), dtype=complex)
for i, wi in enumerate(w):
    H[i] = sum(h[k] * np.exp(-1j * k * wi) for k in range(N))

# 平方幅频响应（优化器使用的量）
R = np.abs(H) ** 2

# 通带和阻带边界的索引
wpass = int(0.12 * len(w))   # 通带：0 到 0.12π
wstop = int(0.20 * len(w))   # 阻带：0.20π 到 π

# 通带检查：应在 [L², U²] 范围内
# （根据 JSON 规格：ripple=0.125 → L² ≈ 0.64, U² ≈ 1.44）
passband_ok = np.all(R[:wpass] >= 0.64) and np.all(R[:wpass] <= 1.44)

# 阻带检查：应低于衰减目标
# （根据规格：attenuation=0.125 → S² ≈ 0.0156）
stopband_max = np.max(R[wstop:])
stopband_ok = stopband_max <= 0.0156

print(f"通带 OK: {passband_ok}")
print(f"阻带最大 R(ω): {stopband_max:.6f}")
print(f"阻带 OK: {stopband_ok}")
```

### 9.2 使用 scipy.signal.freqz 验证

使用 SciPy 的更简单方法（如果已安装）：

```python
from scipy.signal import freqz
import numpy as np
import json
import matplotlib.pyplot as plt

with open("lowpass_32_output.json") as f:
    data = json.load(f)

h = np.array([c["value"] for c in data["coefficients"]])
N = len(h)

w, H = freqz(h, worN=4096)
w = w / np.pi  # 归一化到 [0, 1]

# 绘图
plt.figure(figsize=(12, 4))

plt.subplot(1, 2, 1)
plt.plot(h, 'o-', markersize=3)
plt.title(f"脉冲响应 ({N} 抽头)")
plt.xlabel("抽头索引")
plt.ylabel("h[n]")

plt.subplot(1, 2, 2)
plt.plot(w, 20 * np.log10(np.abs(H) + 1e-12))
plt.axvline(0.12, color='orange', linestyle='--', label='通带边界')
plt.axvline(0.20, color='red', linestyle='--', label='阻带边界')
plt.ylim(-80, 5)
plt.title("幅频响应 (dB)")
plt.xlabel("归一化频率 (×π)")
plt.ylabel("|H(ω)| [dB]")
plt.legend()
plt.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig("fir_response.png", dpi=150)
plt.show()
```

### 9.3 检查清单

| 检查项 | 预期 | 方法 |
|--------|------|------|
| **通带平坦度** | R(ω) ∈ [L², U²]，ω ∈ [0, ωₚ] | `np.min(R[:wpass])` ≥ L²，`np.max(R[:wpass])` ≤ U² |
| **阻带衰减** | R(ω) ≤ S²，ω ∈ [ωₛ, π] | `np.max(R[wstop:])` ≤ S² |
| **脉冲响应形状** | 关于中心对称（线性相位）或最小相位 | 目视检查 h 系数 |
| **无 NaN/Inf** | 所有系数有限 | `np.all(np.isfinite(h))` |
| **CSD 位数** | ≤ 每系数 csd_nnz | 统计 CSD 字符串中 `+` 和 `-` |

### 9.4 验证脚本

项目包含验证辅助工具：

```bash
python tools/verify_filter.py lowpass_32_output.json
```

输出显示每个规格频带的通过/失败。

### 9.5 验证失败的排查

验证脚本在**两种分辨率**下检查：

```
--- Oracle 网格（480 点，disc_factor=15）---
[PASS] 通带：R(w) ∈ [0.7932, 1.2646]
[PASS] 阻带：max R(w) = 0.000000 ≤ 0.015625

--- 密集网格（16384 点，参考用）---
[FAIL] 通带：R(w) ∈ [0.7679, 1.2647]
       （离散化伪影：增大 disc_factor）
[PASS] 阻带：max R(w) = 0.000000 ≤ 0.015625

全部检查通过
```

**Oracle 网格是权威的**：优化器仅在采样频率点（15×N）保证可行性。
密集网格的通带下降（0.7679 vs 0.7901）发生在 Oracle 采样点**之间** —
这是离散化伪影，而非量化误差。

> 离散优化器（`cutting_plane_optim_q`）配合 CSD 量化返回 `rcsd`（量化后的自相关），
> 在所有采样频率点均满足规格。不存在量化误差。

| 症状 | 可能原因 | 修复方法 |
|------|---------|---------|
| Oracle 网格失败（罕见） | csd_nnz 太低或收敛不良 | 增大 `max_iters` 或 `csd_nnz` |
| 仅密集网格失败 | 采样点间的离散化间隙 | 将 `discretization_factor` 增大到 30 |
| 阻带衰减不足 | 抽头数不够 | 将 `filter_order` 增加 8–16 |
| 频率响应不稳定 | 收敛容差太松 | 将 `tolerance` 降至 `1e-20` |

---

## 10. 调优指南

### "滤波器不收敛"

按顺序尝试：
1. **增大 `passband_ripple`** — 更大的纹波（如 0.2）更容易实现
2. **增大 `max_iters`** — 50000 → 100000 给求解器更多时间
3. **放宽 `tolerance`** — `1e-14` → `1e-10` 更早退出
4. **减小 `filter_order`** — 更少的抽头 = 更小的问题

### "我想减少加法器（更小面积）"

- **减小 `csd_nnz`** — 从 4 → 3 节省约 25% 加法器，但可能降低阻带性能。
  检查输出系数以确认阻带衰减仍可接受。

### "我需要更好的阻带衰减"

- **增大 `filter_order`** — 更多抽头 = 更陡峭的滚降
- **增大 `csd_nnz`** — 更多非零位 = 更好地逼近理想系数
- **收窄过渡带** — 将 `stopband_edge` 向 `passband_edge` 靠近（但这会增加收敛难度）

### "我需要更精确的频率响应"

- **增大 `discretization_factor`** — 从 15 → 30 使频率采样点翻倍（代价是更多内存和更慢的迭代）

### 参数影响速查表

| 参数 | 增大 → | 效果 |
|------|--------|------|
| `filter_order` | 更多抽头 | 阻带更好，更慢 |
| `csd_nnz` | 更多位数 | 精度更好，面积更大 |
| `max_iters` | 更多迭代 | 收敛更好，更慢 |
| `tolerance` | 更大值 | 更快退出，最优性较松 |
| `discretization_factor` | 更多采样点 | 约束更精确，内存更大 |
| `passband_ripple` | 更大 | 更容易设计，纹波更大 |
| `stopband_attenuation` | 更小 | 阻带更深，更难收敛 |
| `spectral_method` | `"fft"` | 使用传统 FFT，精确但较慢 |
| `root_tolerance` | 更大值 | 收敛更快，根精度较松 |

---

## 🔧 故障排除

| 问题 | 解决方案 |
|------|---------|
| `Error: cannot open file` | 检查路径 — 使用绝对路径或相对于当前目录 |
| `JSON parse error` | 在 [jsonlint.com](https://jsonlint.com) 验证 JSON |
| `Optimization failed — no feasible solution` | 参见上文调优指南 |
| CLI 卡住 | 减少 `max_iters` 或 `filter_order` 进行测试 |
| `fir-design: command not found` | 在 Python 项目中运行 `pip install -e .` |
| 输出中无 Verilog | 在 JSON 规格中添加 `"verilog": {"input_width": 16}` |
| "Aberth 未收敛" | 增大 `root_tolerance`（如 `1e-6`）或切换到 `"fft"` |
| 输出中无 Verilog | 在 JSON 规格中添加 `"verilog": {"input_width": 16}` |

---

## 📚 参考文献

- **论文**：Wu, Boyd & Vandenberghe — "FIR Filter Design via Spectral
  Factorization and Convex Optimization"
- **C++ 仓库**：[luk036/multiplierless-cpp](https://github.com/luk036/multiplierless-cpp)
- **Python 仓库**：[luk036/multiplierless](https://github.com/luk036/multiplierless)
- **CSD 库 (C++)**：[luk036/csd-cpp](https://github.com/luk036/csd-cpp)
- **CSD 库 (Python)**：[luk036/csdigit](https://github.com/luk036/csdigit)
- **椭球算法**：[luk036/ellalgo](https://github.com/luk036/ellalgo)
