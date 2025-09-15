# RTK 定位原理与应用（含公式推导）

## 1. 基本观测方程

GNSS 接收机的两类主要观测量为伪距 $P$ 和载波相位 $\Phi$：

- 伪距观测方程：
  
$$
P_r^s = \rho_r^s + c(\delta t_r - \delta t^s) + I_r^s + T_r^s + \varepsilon_P
$$

- 载波相位观测方程：
  
$$
\Phi_r^s = \rho_r^s + c(\delta t_r - \delta t^s) + I_r^s + T_r^s + \lambda N_r^s + \varepsilon_\Phi
$$

其中：  
- $\rho_r^s$：接收机 $r$ 到卫星 $s$ 的几何距离  
- $\delta t_r, \delta t^s$：接收机与卫星钟差  
- $I_r^s, T_r^s$：电离层、对流层延迟  
- $\lambda N_r^s$：整周模糊度项  
- $\varepsilon$：观测噪声  

---

## 2. 差分处理与残差推导

### 2.1 单差（Between Receivers）
对同一颗卫星 $s$，在接收机 $r_1$ 和 $r_2$ 之间做差：

$$
\Delta \Phi_{r_1,r_2}^s = (\Phi_{r_1}^s - \Phi_{r_2}^s) 
= (\rho_{r_1}^s - \rho_{r_2}^s) + c(\delta t_{r_1} - \delta t_{r_2}) + \Delta I + \Delta T + \lambda \Delta N + \Delta \varepsilon
$$

- 消除了卫星钟差 $\delta t^s$。  
- 仍然含有接收机钟差。  

---

### 2.2 双差（Between Satellites and Receivers）
在单差的基础上，对两颗卫星 $s_1, s_2$ 做差：

$$
\nabla\Delta \Phi_{r_1,r_2}^{s_1,s_2} 
= \Delta \Phi_{r_1,r_2}^{s_1} - \Delta \Phi_{r_1,r_2}^{s_2}
$$

展开后得到：

$$
\nabla\Delta \Phi_{r_1,r_2}^{s_1,s_2} 
= (\rho_{r_1}^{s_1} - \rho_{r_2}^{s_1}) - (\rho_{r_1}^{s_2} - \rho_{r_2}^{s_2}) + \lambda (\Delta N^{s_1} - \Delta N^{s_2}) + \nabla\Delta \varepsilon
$$

此时：  
- **卫星钟差**与**接收机钟差**均被完全消除。  
- 残差主要由**几何距离差**、**整数模糊度差**和**剩余误差项**组成。  

---

### 2.3 线性化与残差方程

几何距离 $\rho$ 与接收机未知坐标 $(X,Y,Z)$ 之间是非线性的，需在近似坐标 $(X_0,Y_0,Z_0)$ 处进行一阶线性化展开：

$$
\nabla\Delta \Phi \approx \mathbf{H} \cdot \Delta \mathbf{x} + \lambda \cdot \mathbf{N} + \varepsilon
$$

其中：  
- $\Delta \mathbf{x} = [\Delta X, \Delta Y, \Delta Z]^T$ 为位置改正量  
- $\mathbf{H}$ 为 **设计矩阵（特征向量矩阵）**，由观测方向余弦构成  
- $\mathbf{N}$ 为待解的整数模糊度向量  

---

## 3. 特征向量（设计矩阵 H）的构建

对第 $i$ 颗卫星，方向余弦可写为：

$$
l_i = \frac{X^s_i - X_r}{\rho_i}, \quad
m_i = \frac{Y^s_i - Y_r}{\rho_i}, \quad
n_i = \frac{Z^s_i - Z_r}{\rho_i}
$$

于是，线性化后的观测方程为：

$$
v_i = -l_i \Delta X - m_i \Delta Y - n_i \Delta Z + \lambda N_i + \varepsilon_i
$$

将多个卫星的方程组合，矩阵形式为：

$$
\mathbf{v} = \mathbf{H} \Delta \mathbf{x} + \lambda \mathbf{N} + \varepsilon
$$

其中：

$$
\mathbf{H} = \begin{bmatrix}-l_1 & -m_1 & -n_1 \\-l_2 & -m_2 & -n_2 \\\vdots & \vdots & \vdots \\-l_n & -m_n & -n_n\end{bmatrix}
$$

---

## 4. 解算与滤波

- 在 **模糊度固定**之前， $$\mathbf{N}$$ 是实数解（浮点模糊度）。  
- 通过 **LAMBDA 算法**进行整数搜索与固定。  
- 固定后的位置参数 $$\Delta \mathbf{x}$$ 可通过最小二乘或卡尔曼滤波解算：

$$
\hat{\Delta \mathbf{x}} = (\mathbf{H}^T \mathbf{P}^{-1} \mathbf{H})^{-1} \mathbf{H}^T \mathbf{P}^{-1} (\mathbf{v} - \lambda \mathbf{N})
$$

---
