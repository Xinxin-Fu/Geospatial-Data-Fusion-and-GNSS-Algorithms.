# RTKLIB 最小二乘函数 `lsq` 解析

## 1. 源码

```c
/* --------------------------------------------------------------------------
 * Least Squares Solution
 * --------------------------------------------------------------------------
 * 功能：
 *   解最小二乘问题 Ax ≈ y，求解 x = (A^T A)^(-1) A^T y
 *
 * 输入参数：
 *   A  - 设计矩阵 (m × n)，观测方程系数           [输入]
 *   y  - 观测向量 (m × 1)                         [输入]
 *   m  - 观测数                                     [输入]
 *   n  - 待估参数数                                 [输入]
 *
 * 输出参数：
 *   x  - 最小二乘解向量 (n × 1)                   [输出]
 *   Q  - 协方差矩阵 (n × n)，可选                  [输出，可选]
 *
 * 返回值：
 *   0   - 成功
 *  -1   - 失败（观测数不足或求解失败）
 * -------------------------------------------------------------------------- */
extern int lsq(const double *A, const double *y, int m, int n, double *x,
               double *Q)
{
    double *AtA,*Aty;
    int info;

    if (m<n) return -1;

    AtA=mat(n,n); Aty=mat(n,1);
    matmul("TN",n,n,m,1.0,A,A,0.0,AtA); /* AtA = A'*A */
    matmul("TN",n,1,m,1.0,A,y,0.0,Aty); /* Aty = A'*y */

    /* solve normal equation AtA*x = Aty */
    info = solve("N", AtA, Aty, n, 1, x);
    if (info) {
        free(AtA); free(Aty);
        return -1;
    }
    if (Q) {
        matcpy(Q, AtA, n, n);
    }
    free(AtA); free(Aty);
    return 0;
}
```
## 2. 整体的调用关系（SPP）
```c
pntpos()   # 单点定位入口
   ├──> rescode()    # 构建设计矩阵 A 和观测向量 y
   ├──> lsq()        # 解 (A^T A)x = A^T y
   ├──> 更新解 x     # x = x + dx
   └──> 收敛判断      # 若未收敛则继续迭代
```
## 3. 牛顿迭代的步骤
```c
/* --------------------------------------------------------------------------
 * Single Point Positioning (SPP)
 * --------------------------------------------------------------------------
 * 功能：
 *   使用伪距观测数据进行单点定位，通过牛顿迭代和最小二乘法求解
 *
 * 输入参数：
 *   obs   - 观测数据数组（伪距、载波等）          [输入]
 *   n     - 观测数量                                [输入]
 *   nav   - 卫星导航电文（星历、钟差等）           [输入]
 *   opt   - 定位选项（解算模式、测站类型等）       [输入]
 *
 * 输出参数：
 *   sol   - 接收机解（位置、钟差等）               [输出]
 *   azel  - 卫星高度角和方位角                     [输出]
 *   ssat  - 卫星状态信息                             [输出]
 *   msg   - 错误或状态信息                           [输出]
 *
 * 返回值：
 *   状态码（SOLQ_NONE / SOLQ_SINGLE 等）
 * -------------------------------------------------------------------------- */
extern int pntpos(const obsd_t *obs, int n, const nav_t *nav,
                  const prcopt_t *opt, sol_t *sol, double *azel,
                  ssat_t *ssat, char *msg)
{
    int i, j, iter, stat;           // 循环计数和状态标志
    int vsat[MAXOBS];                // 有效卫星索引
    int ns = 0;                      // 有效观测数量
    double x[4] = {0};               /* 接收机位置和钟差 {X, Y, Z, dt} */
    double dx[4];                    // 最小二乘增量
    double *A, *y, *Q;               // 设计矩阵、观测向量、协方差矩阵
    double resp[MAXOBS], vmat[MAXOBS*MAXOBS]; // 临时残差和方差矩阵

    for (iter = 0; iter < MAXITR; iter++) {
        // 构建伪距残差和设计矩阵
        ns = rescode(obs, n, nav, x, opt, y, A, azel, vsat, resp, vmat);
        if (ns < 4) break; // 有效观测数不足，退出

        /* 最小二乘估计 */
        if (lsq(A, y, ns, 4, dx, Q)) {
            stat = SOLQ_NONE; 
            break; // 求解失败
        }

        /* 更新解 */
        for (j = 0; j < 4; j++) x[j] += dx[j];

        /* 收敛判断 */
        if (norm(dx, 4) < 1E-4) { 
            stat = SOLQ_SINGLE; 
            break;
        }
    }

    // 最终结果赋值到 sol 等输出参数
    ...
}
```




