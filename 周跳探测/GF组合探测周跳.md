# GF组合探测周跳原理与公式推导

## 1. 定义  
GNSS 双频载波相位观测值：  

$$
\Phi_i = \rho + c(\delta t_r - \delta t_s) + T - I \cdot \frac{f_1^2}{f_i^2} + \lambda_i N_i + \epsilon_i, \quad i=1,2
$$

- $\rho$：几何距离  
- $c(\delta t_r - \delta t_s)$：接收机和卫星钟差  
- $T$：对流层延迟  
- $I$：电离层延迟在频率 $f_1$ 下的量值  
- $\lambda_i$：波长  
- $N_i$：整周模糊度  
- $\epsilon_i$：测量噪声  

---

## 2. GF组合构造  
GF组合定义为：  

$$
\Phi_{\mathrm{GF}} = \Phi_1 - \Phi_2
$$

代入各项：  

$$
\Phi_{\mathrm{GF}} = (\lambda_1 N_1 - \lambda_2 N_2) + \left(-I + I \cdot \frac{f_1^2}{f_2^2}\right) + (\epsilon_1 - \epsilon_2)
$$

化简：  

$$
\Phi_{\mathrm{GF}} = \lambda_1 N_1 - \lambda_2 N_2 + I\left(\frac{f_1^2}{f_2^2} - 1\right) + \epsilon
$$

其中： $\epsilon = \epsilon_1 - \epsilon_2$。

---

## 3. 特点  
- **钟差、对流层延迟、几何距离**在GF组合中完全消除。  
- 剩下的主要是 **电离层项 + 整周模糊度项 + 噪声**。  

---

## 4. 用于周跳探测
