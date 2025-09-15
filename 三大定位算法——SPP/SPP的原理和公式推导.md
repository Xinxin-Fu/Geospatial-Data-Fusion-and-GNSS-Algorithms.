# 伪距单点定位（SPP）

## 原理
- **伪距（Pseudorange）**：接收机测得的卫星信号传播时间乘以光速得到的距离。由于接收机时钟与卫星时钟存在偏差，以及电离层、对流层延迟等因素，这个“距离”并不是真实的几何距离，因此称为“伪距”。  
- **定位思路**：接收机同时接收 ≥4 颗卫星的伪距观测值，建立方程组，解算接收机三维坐标（X, Y, Z）及接收机时钟偏差。

## 观测方程如下：
$$
\begin{aligned}
P_i &= \|\mathbf{r} - \mathbf{r}_i\| + c\bigl(\delta t_r - \delta t_{s,i}\bigr) + I_i + T_i + M_i + \varepsilon_i \\
    &= \rho_i + c\,\delta t_r + I_i + T_i + M_i + \varepsilon_i
\end{aligned}
$$


### 上式中：
- $$P_i$$ ：第 $$i$$ 颗卫星的伪距观测值
- $$\mathbf{r} = (x, y, z)^\top$$ ：接收机三维位置
- $$\mathbf{r}_i = (x_i, y_i, z_i)^\top$$ ：第 $$i$$ 颗卫星位置
- $$\rho_i = \|\mathbf{r} - \mathbf{r}_i\|$$ ：几何距离
- $$c$$ ：光速
- $$\delta t_r, \delta t_{s,i}$$ ：接收机钟差与第 $$i$$ 颗卫星钟差
- $$I_i, T_i$$ ：电离层与对流层延迟
- $$M_i$$ ：多路径误差
- $$\varepsilon_i$$ ：观测噪声及其他小误差
  
## 线性化（用于最小二乘迭代求解）
在近似位置 $$\hat{r} $$与钟差  $$\hat{\delta t} $$处线性化：

$$
P_i - \hat{P}_i \approx \frac{(\hat{\mathbf r} - \mathbf r_i)^\top}{\hat{\rho}_i}\,\Delta\mathbf r + c\,\Delta\delta t + \nu_i
$$

### 上式中：
- $$\Delta \mathbf r = \mathbf r - \hat{\mathbf r}$$ ：位置改正量  
- $$\Delta \delta t = \delta t_r - \hat{\delta t}$$ ：钟差改正量  
- $$\hat{\rho}_i = \|\hat{\mathbf r} - \mathbf r_i\|$$ ：近似几何距离  
- $$\nu_i$$ ：线性化残差（含未建模误差）  
---

接下来可通过 **构建残差向量**和**设计矩阵** \(H\) 来进行最小二乘迭代求解：

$$
\mathbf{v} = \mathbf{P} - \hat{\mathbf{P}} \approx H \Delta \mathbf{x} + \mathbf{\nu}
$$

其中：
- $$\mathbf{v}$$ ：观测残差向量  
- $$\mathbf{P}$$ ：观测伪距向量  
- $$\hat{\mathbf{P}}$$ ：预测伪距向量  
- $$\Delta {x} $$ ：状态改正量  
- $$H$$ ：设计矩阵，由各卫星几何关系构成  
- $$\mathbf{\nu}$$ ：未建模误差向量
---
# 伪距定位的特点和应用场景
## 特点
- **观测值简单**：只使用伪距，不需要载波相位。  
- **无需基站**：仅依赖广播星历即可（但若要求更高精度，需进行差分或PPP/RTK 等处理）。  
- **精度有限**：一般位于米级（约 5–15 m），受系统、环境、星历精度等影响。  

---

## 主要应用
1. **大众导航**（手机、车载导航、手持 GPS）。  
2. **紧急与通用服务**（紧急呼救定位、车辆调度、物流跟踪）。  
3. **作为高精度方法的初始化**（为 RTK、PPP提供初始解）
