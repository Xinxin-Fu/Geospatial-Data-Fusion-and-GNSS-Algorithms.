# 精密单点定位（PPP, Precise Point Positioning）

## 1. 概念  
PPP 是一种基于单站观测的精密定位方法。它直接利用 GNSS 卫星的 **精密轨道和钟差产品**（由 IGS 等机构提供），结合双频伪距和载波相位观测值，经过误差改正和参数估计，可以在单接收机条件下实现 **厘米级至分米级精度** 的定位。  

与传统的 **差分定位（如 RTK、DGPS）** 不同，PPP 不依赖于基准站或基准站网，属于 **绝对定位** 方法。  

---

## 2. 优势  
1. **不依赖基准站**：用户只需要一台接收机和外部精密产品即可定位，减少基站建设和维护成本。  
2. **全球统一基准框架**：PPP 使用 IGS 等提供的产品，定位结果直接在 ITRF（国际地球参考框架）下，具有全球统一性。  
3. **高精度**：理论上可达到 **厘米级定位精度**（静态长期观测）和 **分米级精度**（动态、实时）。  
4. **可扩展性强**：PPP 能够适配多系统（GPS、BDS、GLONASS、Galileo）和多频观测，精度与收敛速度持续提升。  

---

## 3. 劣势  
1. **收敛时间长**：传统 PPP 需要 20–60 分钟才能达到厘米级精度。虽然 **PPP-AR**（整周模糊度固定）和 **PPP-RTK**（结合区域改正信息）能缩短收敛，但仍比 RTK 慢。  
2. **依赖外部精密产品**：轨道、钟差、大气改正等精密产品由 IGS 或商业服务商提供，实时 PPP 需要稳定的数据链路。  
3. **对应用环境有要求**：在遮挡严重或多路径干扰强烈的环境下，收敛更慢、精度更差。  

---

## 4. 主要应用方向  
### （1）科学研究  
- 地壳形变监测、地震与火山活动研究  
- 精密大地测量、参考框架维护  
- 大气监测（对流层水汽、空间天气电离层）  

### （2）工程与行业  
- 海洋测量与精密导航（船舶、航空）  
- 遥感卫星的高精度定轨  
- 精密农业、工程测量  

### （3）大众服务  
- 智能交通与自动驾驶（结合 PPP-RTK）  
- 无人机导航  
- 实时高精度地图更新  

---

## 5. 总结  
PPP 的最大优势是 **不依赖基站、全球统一、成本低**，但 **收敛时间长** 是主要瓶颈。  
它特别适合全球分布式的长期精密监测和低成本高精度应用；  
而对 **实时性要求高** 的场景，通常采用 PPP-RTK 或 RTK 作为补充。  

---

# PPP 公式推导（从观测方程到残差和设计矩阵）

## 1. 伪距与载波观测方程

在卫星导航中，接收机对卫星的伪距 (code) 和载波相位 (carrier phase) 观测值可表示为：

$$
P_i^s = \rho^s + c \cdot (dt_r - dt_s) + T^s + I_i^s + d_{P,i}^s + \epsilon_{P,i}^s
$$

$$
\Phi_i^s = \rho^s + c \cdot (dt_r - dt_s) + T^s - I_i^s + \lambda_i N_i^s + d_{\Phi,i}^s + \epsilon_{\Phi,i}^s
$$

其中：

- $P_i^s, \Phi_i^s$ ：频率 $ i $ 的伪距与载波观测值  
-  $$\rho^s $$ ：几何距离  
-  $$dt_r, dt_s $$ ：接收机钟差与卫星钟差  
-  $$T^s $$ ：对流层延迟  
-  $$I_i^s $$ ：电离层延迟  
-  $$\lambda_i $$ ：载波波长  
-  $$N_i^s $$ ：整周模糊度  
-  $$d_{P,i}^s, d_{\Phi,i}^s $$ ：硬件延迟  
-  $$\epsilon $$ ：测量噪声  

---

## 2. 无电离层组合 (IF)

为消除一阶电离层延迟，采用无电离层组合：

$$
P_{IF}^s = \alpha P_1^s - \beta P_2^s
$$

$$
\Phi_{IF}^s = \alpha \Phi_1^s - \beta \Phi_2^s
$$

其中组合系数为：

$$
\alpha = \frac{f_1^2}{f_1^2 - f_2^2}, \quad
\beta = \frac{f_2^2}{f_1^2 - f_2^2}
$$

组合后的观测方程：

$$
P_{IF}^s = \rho^s + c \cdot (dt_r - dt_s) + T^s + d_{P,IF}^s + \epsilon_{P,IF}^s
$$

$$
\Phi_{IF}^s = \rho^s + c \cdot (dt_r - dt_s) + T^s + \lambda_{IF} N_{IF}^s + d_{\Phi,IF}^s + \epsilon_{\Phi,IF}^s
$$

其中一阶电离层项已完全消除。

---

## 3. 线性化

假设接收机坐标为：

$$
\mathbf{x} = (X, Y, Z)^T
$$

卫星坐标为 $$ (X^s, Y^s, Z^s) $$，几何距离为：

$$
\rho^s = \sqrt{(X^s - X)^2 + (Y^s - Y)^2 + (Z^s - Z)^2}
$$

对未知参数 $$ \mathbf{x}, dt_r $$ 在线性化展开：

$$
\delta \rho^s \approx -e_x^s \cdot \delta X - e_y^s \cdot \delta Y - e_z^s \cdot \delta Z
$$

其中方向余弦为：

$$
\mathbf{e}^s = \frac{1}{\rho^s} 
\begin{bmatrix}
X^s - X \\
Y^s - Y \\
Z^s - Z
\end{bmatrix}
$$

---

## 4. 残差方程

对 IF 组合的观测方程线性化后，残差写作：

$$
v^s = L^s - \hat{\rho}^s - c \cdot \delta dt_r - \delta T^s - \lambda_{IF} \delta N^s
$$

其中：
- $$L^s $$ ：观测值 (伪距或载波)  
- $$\hat{\rho}^s $$ ：由近似坐标计算的几何距离  
- $$v^s $$ ：残差  

---

## 5. 矩阵形式

将所有卫星的观测方程写成矩阵形式：

$$
\mathbf{v} = \mathbf{A} \cdot \mathbf{x} - \mathbf{l}
$$

其中：  

- 残差向量  

$$
\mathbf{v} =
\begin{bmatrix}
v^1 \\
v^2 \\
\vdots \\
v^n
\end{bmatrix}
$$

- 设计矩阵  

$$
\mathbf{A} =
\begin{bmatrix}
-e_x^1 & -e_y^1 & -e_z^1 & 1 & M^1 & \cdots \\
-e_x^2 & -e_y^2 & -e_z^2 & 1 & M^2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \vdots & \\
-e_x^n & -e_y^n & -e_z^n & 1 & M^n & \cdots
\end{bmatrix}
$$

其中各列对应：  
- 接收机位置改正量 $$(\delta X, \delta Y, \delta Z) $$  
- 接收机钟差 $$\delta dt_r $$  
- 对流层参数 $$\delta T $$ (通过映射函数 $$ M^s $$)  
- 模糊度参数 $$\delta N^s $$  

- 观测向量  

$$
\mathbf{l} =
\begin{bmatrix}
L^1 - \hat{\rho}^1 \\
L^2 - \hat{\rho}^2 \\
\vdots \\
L^n - \hat{\rho}^n
\end{bmatrix}
$$

---

## 6. 最终 PPP 观测方程形式

因此，PPP 的线性化观测方程可写为：

$$
\mathbf{v} = \mathbf{A} \mathbf{x} - \mathbf{l}
$$

这就是 PPP 由原始观测值推导到残差方程和设计矩阵的全过程。



