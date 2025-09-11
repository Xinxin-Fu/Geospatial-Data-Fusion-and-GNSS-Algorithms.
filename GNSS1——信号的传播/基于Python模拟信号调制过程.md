# GNSS 信号调制过程模拟（Python 示例）

下面通过 Python 对 GNSS 信号的三层调制过程进行模拟，每一层都进行输出和可视化。

---

## 1. 载波层级

载波层提供高频正弦信号作为信息承载基础。  

```python
import numpy as np
import matplotlib.pyplot as plt

fs = 1000  # 采样频率 Hz
t = np.arange(0, 1, 1/fs)  # 1秒时间
fc = 50  # 载波频率 Hz
carrier = np.sin(2 * np.pi * fc * t)

plt.figure(figsize=(12, 3))
plt.plot(t[:200], carrier[:200])  # 只画前200个点便于观察
plt.title("1. 载波信号（Carrier）")
plt.xlabel("时间 [s]")
plt.ylabel("幅度")
plt.show()
```
## 2. 码层级

码层通过伪随机序列实现卫星区分和粗测距。

```python
nav_bits = np.random.choice([1, -1], size=int(len(t)/10))  # 低速比特
nav_bits_upsampled = np.repeat(nav_bits, 10)  # 上采样到与载波长度一致
nav_signal = code_signal * nav_bits_upsampled  # 叠加电文

plt.figure(figsize=(12, 3))
plt.plot(t[:200], nav_signal[:200])
plt.title("3. 导航电文叠加（Navigation Message Layer）")
plt.xlabel("时间 [s]")
plt.ylabel("幅度")
plt.show()
```
## 3. 导航电文层级（Navigation Message）

导航电文承载卫星轨道、时间和状态信息。

```python
# 假设导航电文为简单的二进制序列
nav_bits = np.random.choice([1, -1], size=int(len(t)/10))  # 低速比特
nav_bits_upsampled = np.repeat(nav_bits, 10)  # 上采样到与载波长度一致
nav_signal = code_signal * nav_bits_upsampled  # 叠加电文

plt.figure(figsize=(12, 3))
plt.plot(t[:200], nav_signal[:200])
plt.title("3. 导航电文叠加（Navigation Message Layer）")
plt.xlabel("时间 [s]")
plt.ylabel("幅度")
plt.show()
```
