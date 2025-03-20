title: 水平扩散
---

## 2.1 水平扩散

相较于水平输送，水平扩散的影响通常较小，在某些情况下甚至可以忽略不计。

\[
\begin{align*}
\frac{\partial c}{\partial t} &=
\frac{\partial}{\partial x} \left(K_x \rho\frac{\partial (\frac{c}{\rho})}{\partial x} \right) \\
\frac{\partial c}{\partial t} &=
\frac{\partial}{\partial y} \left(K_y \rho\frac{\partial (\frac{c}{\rho})}{\partial y} \right) \\
\end{align*}
\]

根据 Smagorinsky (1963) 公式，可求解水平扩散系数 
\[
   K_{x, y} = C_s \Delta^2 |D|
\]

其中：

- \( C_s \) 是 Smagorinsky 常数（通常取 **0.1 - 0.25**，大气应用中常用 **0.2**）。
- \( \Delta \) 是网格尺度（通常定义为网格间距的几何平均：\( \Delta = (\Delta x \Delta y )^{1/2} \)）。
- \( |\mathbf{D}| \) 是水平速度的二阶张量的模量，定义为：
  \[
  |\mathbf{D}| = [(\frac{\partial u}{\partial y} +
    \frac{\partial v}{\partial x})^2
  +(\frac{\partial u}{\partial x} -
    \frac{\partial v}{\partial y})^2]^{1/2}
  \]

可以增加一个常数，让模拟的浓度变得更加光滑（Anthes and Warner, 1978）。
\[
   K_{x, y} = K_0 + C_s \Delta^2 |D|
\]

其中\(K_0\)的值为

\[
  K_0 = 3 \times 10^{-3} \frac{\Delta x \Delta y}{\Delta t}
\]

为了便于通量计算，水平扩散系数\(K_{x, y}\) 通常是在 Arakawa C 网格上定义，对于一个网格单元，其四条边上的每一条边均需计算一个\(K\)值，以准确描述水平扩散过程。

\(K\)的计算涉及到速度的空间差分，采用中值差分方案求速度梯度。

\[
(\frac{\partial u}{\partial x})_{i+1/2, j} = 
\frac{u_{i+1/2+1, j} - u_{i+1/2-1, j}} {2 \Delta x}
\]

\[
(\frac{\partial u}{\partial y})_{i+1/2, j} = 
\frac{u_{i+1/2, j+1} - u_{i+1/2, j-1}}{2 \Delta y}
\]

在求解 v-stag 网格点上的 \(v\) 在 u-stag 网格上的梯度时，需要首先将 
\(v\)从 v-stag 网格插值到 mass 网格。为了简单，通常采用线性插值方法。

\[
v(i, j) = 
\frac{v(i, j+\frac{1}{2}) + v(i, j-\frac{1}{2})}{2}
\]

因此
\[
(\frac{\partial v}{\partial x})_{i+1/2, j} = 
\frac{
  v(i+1, j) - 
  v(i, j)
}{\Delta x}
\]

对于扩散方程，采用前向欧拉进行时间差分。
\[
\frac{\partial c}{\partial t} \approx \frac{c_{t+1} - c_{t}}{\Delta t} 
\]

而在计算浓度梯度时，采用中值差分方案计算浓度梯度（注意计算的是 u-stag 或者 v-stag 网格上的浓度梯度）。

\[
\rho\frac{\partial (\frac{c}{\rho})}{\partial x} \approx
\Delta c^{\rho}_{i+1/2} =
\frac{\rho_{i} + \rho_{i+1}}{\Delta x}
\frac
{\frac{c_{i+1}}{\rho_{i+1}} - \frac{c_{i}}{\rho_{i}} }{\Delta x}
\]

因此，最后的积分表达式为
\[
\frac{c_{t+1} - c_{t}}{\Delta t} =
\frac{
  {Kx_{i+1/2} \times \Delta c^{\rho}_{i+1/2}} -
  {Kx_{i-1/2} \times \Delta c^{\rho}_{i-1/2}}
}{\Delta x}
\]

注意扩散是从高浓度到低浓度，\(K\)始终大于0，浓度梯度决定通量的方向。

 > Smagorinsky J. General circulation experiments with the primitive equations: I. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164.
 > 
 > Anthes R A, Warner T T. Development of hydrodynamic models suitable for air pollution and other mesometerological studies[J]. Monthly Weather Review, 1978, 106(8): 1045-1078.
