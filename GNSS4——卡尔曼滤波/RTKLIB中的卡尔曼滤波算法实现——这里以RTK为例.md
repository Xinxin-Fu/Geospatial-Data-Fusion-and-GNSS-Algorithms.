# RTKLIB中的卡尔曼滤波算法
这里的filter和filter_和主要是针对“矫正”部分，P和H以及v已经在前面处理好了；
```python
/* 卡尔曼滤波 ---------------------------------------------------------------
* 卡尔曼滤波状态更新公式如下：
*
*   K = P * H * (H' * P * H + R)^-1
*   xp = x + K * v
*   Pp = (I - K * H') * P
*
* 参数说明 :
*   double *x        I   状态向量 (n x 1)
*   double *P        I   状态协方差矩阵 (n x n)
*   double *H        I   设计矩阵的转置 (n x m)
*   double *v        I   创新量（观测值 - 预测值）(m x 1)
*   double *R        I   观测误差协方差矩阵 (m x m)
*   int    n,m       I   状态数量和观测数量
*   double *xp       O   更新后的状态向量 (n x 1)
*   double *Pp       O   更新后的状态协方差矩阵 (n x n)
* 返回值 :
*   状态 (0: 正常, <0: 出错)
* 注意事项 :
*   矩阵以列优先存储（Fortran 风格）
*   如果状态 x[i]==0.0，则不更新对应的状态 x[i] 及 P[i+i*n]
*-----------------------------------------------------------------------------*/
extern int filter(double *x, double *P, const double *H, const double *v,
                  const double *R, int n, int m)
{
    double *x_,*xp_,*P_,*Pp_,*H_;
    int i,j,k,info,*ix;
    
    ix=imat(n,1); for (i=k=0;i<n;i++) if (x[i]!=0.0&&P[i+i*n]>0.0) ix[k++]=i;
    x_=mat(k,1); xp_=mat(k,1); P_=mat(k,k); Pp_=mat(k,k); H_=mat(k,m);
    for (i=0;i<k;i++) {
        x_[i]=x[ix[i]];
        for (j=0;j<k;j++) P_[i+j*k]=P[ix[i]+ix[j]*n];
        for (j=0;j<m;j++) H_[i+j*k]=H[ix[i]+j*n];
    }
    info=filter_(x_,P_,H_,v,R,k,m,xp_,Pp_);
    for (i=0;i<k;i++) {
        x[ix[i]]=xp_[i];
        for (j=0;j<k;j++) P[ix[i]+ix[j]*n]=Pp_[i+j*k];
    }
    free(ix); free(x_); free(xp_); free(P_); free(Pp_); free(H_);
    return info;

}
```

## 由于矩阵已经构造好了，下面的过程几乎完全是对公式的实现过程，这里可以对照着课本或者前面给出的公式进行比对：

```python
/* kalman filter ---------------------------------------------------------------
* kalman filter state update as follows:
*
*   K=P*H*(H'*P*H+R)^-1, xp=x+K*v, Pp=(I-K*H')*P
*
* args   : double *x        I   states vector (n x 1)
*          double *P        I   covariance matrix of states (n x n)
*          double *H        I   transpose of design matrix (n x m)
*          double *v        I   innovation (measurement - model) (m x 1)
*          double *R        I   covariance matrix of measurement error (m x m)
*          int    n,m       I   number of states and measurements
*          double *xp       O   states vector after update (n x 1)
*          double *Pp       O   covariance matrix of states after update (n x n)
* return : status (0:ok,<0:error)
* notes  : matirix stored by column-major order (fortran convention)
*          if state x[i]==0.0, not updates state x[i]/P[i+i*n]
*-----------------------------------------------------------------------------*/
static int filter_(const double *x, const double *P, const double *H,
                   const double *v, const double *R, int n, int m,
                   double *xp, double *Pp)
{
    double *F=mat(n,m),*Q=mat(m,m),*K=mat(n,m),*I=eye(n);
    int info;
    
    matcpy(Q,R,m,m);
    matcpy(xp,x,n,1);
    matmul("NN",n,m,n,1.0,P,H,0.0,F);       /* Q=H'*P*H+R */
    matmul("TN",m,m,n,1.0,H,F,1.0,Q);
    if (!(info=matinv(Q,m))) {
        matmul("NN",n,m,m,1.0,F,Q,0.0,K);   /* K=P*H*Q^-1 */
        matmul("NN",n,1,m,1.0,K,v,1.0,xp);  /* xp=x+K*v */
        matmul("NT",n,n,m,-1.0,K,H,1.0,I);  /* Pp=(I-K*H')*P */
        matmul("NN",n,n,n,1.0,I,P,0.0,Pp);
    }
    free(F); free(Q); free(K); free(I);
    return info;
}
```
## 我们这里关注下所谓预测过程是如何实现的
