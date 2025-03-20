title: 垂直扩散
---

## 3.1 垂直扩散

大气垂直扩散（vertical diffusion）是指由于湍流扩散和分子扩散作用，空气中的污染物、热量和水汽在垂直方向上的混合过程。分子扩散由气体分子的热运动引起，其较小，在大气模式中通常可以忽略。行星边界层中湍流作用较强，特别是在白天，垂直扩散在污染物的垂直浓度分布中起着至关重要的作用。

\[
\frac{\partial c}{\partial t} =
\frac{\partial}{\partial z} \left(K_z \rho\frac{\partial (\frac{c}{\rho})}{\partial z} \right)
\]

垂直方向的湍流发展强弱，受温度梯度的影响和风切变的共同影响，
因此垂直扩散系数的求解更为复杂，一般由边界层方案给出（比如YSU方案）。
该项目中不实现垂直扩散系数的计算。

## 3.2 边界条件
假设地面和垂直层顶无湍流扩散（齐次纽曼边界条件）。

\[
  \frac{\partial c}{\partial z} |_{z=0} = 0
\]
\[
  \frac{\partial c}{\partial z} |_{z=h} = 0
\]

## 3.3 数值方案

如果已知垂直扩散系数，采用数值方法求解垂直扩散方程的主要难点在于垂直层一般为非结构化网格（约 10 ～ 1000 米），这可能影响计算稳定性和精度。

由于在边界层内，垂直网格大小在 10 米谅解，如果采用显式时间积分方案会对时间步长 \(dt\) 施加严格限制，要求其极小，从而显著增加计算成本。因此，通常采用隐式时间积分方案来求解离散方程组，以提高计算稳定性并允许较大的时间步长。

与水平扩散一样，在计算浓度梯度时，采用中值差分方案计算浓度梯度（注意计算的是 w-stag 网格上的浓度梯度）。计算 w-stag 网格出的空气密度时（静力平衡且等温的情况，随高度指数衰减），采用一阶线性近似（指数函数是凹函数，一阶线性近似，增加了边界点的密度，会导致垂直扩散增加）。

\[
\rho \frac{\partial (\frac{c}{\rho})}{\partial x} \approx
   \nabla_z c_{i+1/2} =
   \frac{\Delta z_{i+1} \rho_{i} + \Delta z_{i} \rho_{i+1}}
        {\Delta z_{i+1}+\Delta z_{i}}
   \times
   \frac{\frac{c_{i+1}}{\rho_{i+1}} - \frac{c_{i}}{\rho_{i}} }
        {0.5(\Delta z_{i+1}+\Delta z_{i})}
\]

而时间差分方案选用后向欧拉。

\[
\frac{c^{t+1}_i - c^{t}_i }{\Delta t} =
\frac{
  {Kz_{i+1/2} \times \nabla c^{t+1}_{i+1/2}} -
  {Kz_{i-1/2} \times \nabla c^{t+1}_{i-1/2}}
  }{\Delta z_i}
\]

带入浓度梯度公式

\[
\frac{c^{t+1}_i - c^{t}_i }{\Delta t} =
\frac{
   Kz_{i+1/2}
   \frac{\Delta z_{i+1} \rho_{i} + \Delta z_{i} \rho_{i+1}}
        {\Delta z_{i+1}+\Delta z_{i}}
   \frac{\frac{c^{t+1}_{i+1}}{\rho_{i+1}} - \frac{c^{t+1}_{i}}{\rho_{i}}}
         {(\Delta z_{i+1}+\Delta z_{i})} -
   Kz_{i-1/2} 
   \frac{\Delta z_{i} \rho_{i-1} + \Delta z_{i-1} \rho_{i}}
        {\Delta z_{i}+\Delta z_{i-1}}
   \frac{\frac{c^{t+1}_{i}}{\rho_{i}} - \frac{c^{t+1}_{i-1}}{\rho_{i-1}}}
        {(\Delta z_{i}+\Delta z_{i-1})}
}{\Delta z_i}
\]

设
\[
  K^{\rho}_{i+1/2} = 
    2 \frac{Kz_{i+1/2} (\Delta z_{i+1} \rho_{i} + \Delta z_{i} \rho_{i+1})}
           {(\Delta z_{i}+\Delta z_{i+1})^2}
\]

则
\[
\frac{c^{t+1}_i - c^{t}_i }{\Delta t} =
\frac{
  K^{\rho}_{i+1/2} (\frac{c^{t+1}_{i+1}}{\rho_{i+1}} - \frac{c^{t+1}_{i}}{\rho_{i}}) -
  K^{\rho}_{i-1/2} (\frac{c^{t+1}_{i}}{\rho_{i}} - \frac{c^{t+1}_{i-1}}{\rho_{i-1}})
}{\Delta z_i}
\]

将后向欧拉方程改写为矩阵形式 \( A \mathbf{c}^{t+1} = \mathbf{c}^t \)，其中\( A \)为

\[
\begin{pmatrix}
b_1 & c_1 & 0 & \cdots & 0 \\
a_2 & b_2 & c_2 & \cdots & 0 \\
0 & \ddots & \ddots & \ddots & 0 \\
\vdots & \cdots & a_{n-1} & b_{n-1} & c_{n-1} \\
0 & \cdots & 0 & a_n & b_n
\end{pmatrix}
\]

\( A \) 的三对角元素表达式如下：

1. **下对角线元素（第 \( i \) 行，第 \( i-1 \) 列）**：
   \[
   a_i = - \frac{\Delta t}{\Delta z_i} \cdot 
          \frac{K^{\rho}_{i-1/2}}{\rho_{i-1}} < 0
   \]

2. **主对角线元素（第 \( i \) 行，第 \( i \) 列）**：
   \[
   b_i = 
      1 + \frac{\Delta t}{\Delta z_i} 
      \left( \frac{K^{\rho}_{i-1/2}}{\rho_i} + 
             \frac{K^{\rho}_{i+1/2}}{\rho_i} \right) > 0
   \]

3. **上对角线元素（第 \( i \) 行，第 \( i+1 \) 列）**：
   \[
   c_i = - \frac{\Delta t}{\Delta z_i} \cdot 
           \frac{K^{\rho}_{i+1/2}}{\rho_{i+1}} < 0
   \]

\( A \mathbf{c}^{t+1} = \mathbf{c}^t \) 可由 Thomas 算法计算。
如何保障 \(\mathbf{c}^{t+1} > 0 \)，且解是稳定的？

其充分条件是 \(A\) 是严格对角占优的 M-矩阵，
也就是 \(b_i > -c_i - a_i\)

   \[
      1 + \frac{\Delta t}{\Delta z_i} 
      \left( \frac{K^{\rho}_{i-1/2}}{\rho_i} + 
             \frac{K^{\rho}_{i+1/2}}{\rho_i} \right) > 
      \frac{\Delta t}{\Delta z_i} \cdot 
          \frac{K^{\rho}_{i-1/2}}{\rho_{i-1}} +
      \frac{\Delta t}{\Delta z_i} \cdot 
           \frac{K^{\rho}_{i+1/2}}{\rho_{i+1}} 
   \]

\[
\Delta t < \frac{\Delta z_i}
 { \frac{K^{\rho}_{i-1/2}}{\rho_{i-1}} + \frac{K^{\rho}_{i+1/2}}{\rho_{i+1}}
   - \frac{K^{\rho}_{i-1/2}}{\rho_i} - \frac{K^{\rho}_{i+1/2}}{\rho_i}}
\]

\[
\Delta t < \frac{\Delta z_i}{ 2 \left[
\frac{Kz_{i-1/2} (\Delta z_{i} \rho_{i-1} + \Delta z_{i-1} \rho_{i})}
           {(\Delta z_{i-1}+\Delta z_{i})^2} 
\left(\frac{1}{\rho_{i-1}} - \frac{1}{\rho_i} \right) 
+ 
\frac{Kz_{i+1/2} (\Delta z_{i+1} \rho_{i} + \Delta z_{i} \rho_{i+1})}
           {(\Delta z_{i}+\Delta z_{i+1})^2} 
\left(\frac{1}{\rho_{i+1}} - \frac{1}{\rho_i} \right) \right]}.
\]

垂直层中间的 \(\rho\) 比较接近时，\(dt\) 取任意值，差分分案无条件稳定。

## 3.4 单元测试

设空气密度和扩散系数为定值，比如 \(\rho = 1.0\) 且 \(Kz = 1.0\)，
则垂直扩散方程可以写作标准的热传导方程

\[
\frac{\partial c}{\partial t} = K_z \nabla^2 c
\]

纽曼边界条件条件为 \(\frac{\partial c}{\partial z} = 0\)，
用余弦函数构造初始条件
\[
c(z, 0) = \cos\left(\frac{\pi z}{L}\right) + 1
\]

其解析解为
\[
c(z, t) = \cos\left(\frac{\pi z}{L}\right) e^{-\frac{\pi^2 K_z}{L^2} t} + 1
\]

某个网格的平均值，在 \( z_0 \) 到 \( z_1 \) 之间的平均值，我们使用**定积分的平均值公式**：

\[
c_{\text{avg}} = \frac{1}{z_1 - z_0} \int_{z_0}^{z_1} c(z, t) \, dz.
\]

\[
c_{\text{avg}} = e^{-\frac{\pi^2 K_z}{L^2} t} \cdot \frac{L}{\pi(z_1 - z_0)} \left[ \sin\left(\frac{\pi z_1}{L}\right) - \sin\left(\frac{\pi z_0}{L}\right) \right] + 1.
\]

---

 > Hong S Y, Noh Y, Dudhia J. A new vertical diffusion package with an explicit treatment of entrainment processes[J]. Monthly weather review, 2006, 134(9): 2318-2341.
