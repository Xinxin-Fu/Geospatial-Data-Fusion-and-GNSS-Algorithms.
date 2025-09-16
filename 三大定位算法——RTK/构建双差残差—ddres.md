# 构建双差残差
## 我们先介绍构建单插残差`zdres` 函数
### zdres 是 RTKLIB 中 基础的非差残差计算模块，输出残差和线性化信息

```c
static int zdres(
    int base, // 基准站索引（0 表示基准站/参考站）          | 输入    |
    const obsd_t *obs, // 观测数组（伪距 P / 载波 L / SNR）| 输入    |
    int n, // 观测数量（卫星数）                           | 输入    |
    const double *rs, // 卫星位置数组（ECEF 坐标）         | 输入    |
    const double *dts, // 卫星钟差                        | 输入    |
    const double *var, // 卫星观测方差                    | 输入    |
    const int *svh, // 卫星健康状态                       | 输入    |
    const nav_t *nav, // 导航电文/星历信息                 | 输入    |
    const double *rr, // 接收机坐标（ECEF）                | 输入    |
    const prcopt_t *opt, // RTK 解算配置参数               | 输入    |
    int index, //  观测类型标识（0 = 非差残差）             | 输入    |
    double *y, // 残差向量                    | 输出    |
    double *e, // 单位矢量数组（偏导数）       | 输出    |
    double *azel, // 卫星方位角/高度角         | 输出    |
    double *freq  // 卫星频率索引                  | 输出    |
)
```
---
## 其中调用最为核心的函数即是`geodist` 函数，我们将简单介绍

### 函数定义

```c
/* geometric distance ----------------------------------------------------------
* compute geometric distance and receiver-to-satellite unit vector
* args   : double *rs       I   satellilte position (ecef at transmission) (m) 卫星位置（ECEF 坐标，发射时刻）
*          double *rr       I   receiver position (ecef at reception) (m) 接收机位置（ECEF 坐标，接收时刻
*          double *e        O   line-of-sight vector (ecef) 卫星到接收机方向的单位矢量 
* return : geometric distance (m) (0>:error/no satellite position) 返回卫星到接收机的 **几何距离（米）**，包含 Sagnac 修正

* notes  : distance includes sagnac effect correction  如果卫星坐标无效（距离地心小于地球半径），返回 `-1.0` 表示错误
*-----------------------------------------------------------------------------*/
extern double geodist(const double *rs, const double *rr, double *e)
{
    double r;
    int i;
    
    if (norm(rs,3)<RE_WGS84) return -1.0;
    for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
    r=norm(e,3);
    for (i=0;i<3;i++) e[i]/=r;
    return r+OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
}
```

---

### 主要功能

`geodist` 是 RTKLIB 中用于计算 **卫星到接收机的几何距离** 和 **单位方向矢量** 的函数。其核心功能包括：

1. **几何距离计算**

   * 计算卫星与接收机的三维直线距离 `r`

2. **单位矢量计算**

   * 输出卫星指向接收机的 **单位矢量** `e`，用于构建设计矩阵 H 或计算方位角/仰角

3. **Sagnac 修正**

   * 考虑地球自转对信号传播引入的微小误差

---

### 代码解析

```c
if (norm(rs,3)<RE_WGS84) return -1.0;
```

* 判断卫星位置是否有效（距离地心必须大于地球半径）
* 若无效则返回 `-1.0`

```c
for (i=0;i<3;i++) e[i]=rs[i]-rr[i];
r=norm(e,3);
for (i=0;i<3;i++) e[i]/=r;
```

* 计算卫星到接收机的向量 `e = rs - rr`
* 计算向量长度 `r`（几何距离）
* 将向量归一化得到 **单位矢量**

```c
return r + OMGE*(rs[0]*rr[1]-rs[1]*rr[0])/CLIGHT;
```

* 对距离做 **Sagnac 修正**：考虑地球自转对信号传播路径的影响

---

# `ddres` 函数处理流程介绍

## 函数目的

`ddres` 用于计算 **双差（Double-Differenced, DD）伪距和载波相位残差**，并构建滤波器所需的残差向量 `v`、设计矩阵 `H` 和观测协方差 `R`。主要处理过程如下。

---

## 主要处理过程分块说明

### 1️⃣ 基线向量与位置转换

```c
bl=baseline(x,rtk->rb,dr);
ecef2pos(x,posu);
ecef2pos(rtk->rb,posr);
```

* 计算流动站和基站之间的基线向量
* 将 ECEF 坐标转换为大地坐标用于气象模型计算

### 2️⃣ 电离层与对流层延迟计算

```c
for (i=0;i<ns;i++) {
    if (opt->ionoopt>=IONOOPT_EST) ...
    if (opt->tropopt>=TROPOPT_EST) ...
}
```

* 根据观测卫星方位角计算**电离层映射函数**和**对流层延迟**
* 输出 `im`（电离层因子）、`tropu` 和 `tropr`（对流层延迟）、以及对流层延迟对位置的偏导数

### 3️⃣ 系统与频率循环

```c
for (m=0;m<6;m++) ...
for (f=opt->mode>PMODE_DGPS?0:nf;f<nf*2;f++) ...
```

* 对每个卫星系统（GPS/GLO/GAL/BDS/QZS/IRN）和每个频率循环
* 为每个系统选择**参考卫星**（仰角最高）

### 4️⃣ 双差残差计算

```c
v[nv]=(y[f+iu[i]*nf*2]-y[f+ir[i]*nf*2])-(y[f+iu[j]*nf*2]-y[f+ir[j]*nf*2]);
```

* 计算流动站与基站之间的**双差观测残差**
* 考虑伪距或载波观测值

### 5️⃣ 设计矩阵构建（H）

```c
for (k=0;k<3;k++) Hi[k]=-e[k+iu[i]*3]+e[k+iu[j]*3];
```

* 根据卫星单位向量 `e` 构建对位置的偏导数
* 同时考虑电离层、对流层和载波相位偏置的偏导数

### 6️⃣ 电离层、对流层与载波偏置修正

```c
if (opt->ionoopt==IONOOPT_EST) ...
if (opt->tropopt==TROPOPT_EST||opt->tropopt==TROPOPT_ESTG) ...
if (f<nf) { ... DD phase-bias term ... }
```

* 根据观测模型修正残差 `v[nv]`
* 更新设计矩阵 `H` 对应的偏导数列

### 7️⃣ 异常值检测

```c
if (opt->maxinno>0.0&&fabs(v[nv])>opt->maxinno) { ... continue; }
```

* 如果双差残差超过设定阈值，则记录异常卫星并跳过该观测

### 8️⃣ 观测误差方差与有效标记

```c
Ri[nv]=varerr(...);
Rj[nv]=varerr(...);
rtk->ssat[sat[i]-1].vsat[f]=1; ...
```

* 计算单差观测误差方差 `Ri` 和 `Rj`
* 标记卫星观测为有效，供滤波器使用

### 9️⃣ 位运算编码与计数

```c
vflg[nv++]=(sat[i]<<16)|(sat[j]<<8)|((f<nf?0:1)<<4)|(f%nf);
nb[b]++;
```

* 使用整数编码卫星编号、频率和观测类型，节省空间
* 记录每个系统的双差观测计数

### 🔟 移动基线约束与协方差矩阵

```c
if (opt->mode==PMODE_MOVEB&&constbl(...)) ...
ddcov(nb,b,Ri,Rj,nv,R);
```

* 对移动基线模式施加基线长度约束
* 构建双差观测协方差矩阵 `R`

### 1️⃣1️⃣ 释放临时内存与返回

```c
free(Ri); free(Rj); free(im);
free(tropu); free(tropr); free(dtdxu); free(dtdxr);
return nv;
```

* 释放动态分配的临时数组
* 返回有效双差观测数量 `nv`

---




