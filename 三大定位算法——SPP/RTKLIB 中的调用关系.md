# SPP 定位在 RTKLIB 中的调用关系
# # 主要 RTKLIB 函数调用关系图如下图（这里再lixiaohang的基础上进行修改）：
![替代文字](SPP_代码调用逻辑.png)

整体的思路是（读取数据-构建方程（也就是计算残差和设计矩阵这一步也体现了线性化过程）-利用最小二乘法计算结果-对结果进行检验）
# # 我们先对查看最小二乘法的实现：
```python
/* 最小二乘估计 -----------------------------------------------------
* 通过解法方程进行最小二乘估计 (x = (A*A')^-1 * A * y)
* 参数   : double *A        I   （加权）设计矩阵的转置 (n x m)
*          double *y        I   （加权）观测值向量 (m x 1)
*          int    n,m       I   参数个数和观测数 (n <= m)
*          double *x        O   估计的参数向量 (n x 1)
*          double *Q        O   参数协方差矩阵 (n x n)
* 返回值 : 状态 (0: 正常, >0: 出错)
* 说明   : 对于加权最小二乘，将 A 和 y 替换为 A*w 和 w*y (w = W^(1/2))
*          矩阵按列优先存储 (Fortran 方式)
*-----------------------------------------------------------------------------*/
extern int lsq(const double *A, const double *y, int n, int m, double *x,
               double *Q)
{
    double *Ay;
    int info;
    
    if (m<n) return -1;
    Ay=mat(n,1);
    matmul("NN",n,1,m,1.0,A,y,0.0,Ay); /* Ay=A*y */
    matmul("NT",n,n,m,1.0,A,A,0.0,Q);  /* Q=A*A' */
    if (!(info=matinv(Q,n))) matmul("NN",n,1,n,1.0,Q,Ay,0.0,x); /* x=Q^-1*Ay */
    free(Ay);
    return info;
}
```

# # RTKLIB 里 SPP 的最小二乘实现确实非常“公式化”，几乎就是把线性代数公式直接套用。因此更应该关注的部分实际上是如何构建残差和设计矩阵，也就是上述的第二部（rescode函数），首先给出其大体的结构图
![替代文字](SPP_代码调用逻辑.png)
