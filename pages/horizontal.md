title: 水平扩散
---

水平湍流主要是由风切变引起的，在分辨率较高时，相较于水平输送，水平扩散的影响相对较小，
而在分辨率叫粗时，加入水平扩散方案
可以有效弥补粗网格尺度导致的次网格混合不足。

\[
\begin{align*}
\frac{\partial c}{\partial t} &=
\frac{\partial}{\partial x} \left(K_x \rho\frac{\partial (\frac{c}{\rho})}{\partial x} \right) \\
\frac{\partial c}{\partial t} &=
\frac{\partial}{\partial y} \left(K_y \rho\frac{\partial (\frac{c}{\rho})}{\partial y} \right) \\
\end{align*}
\]

## 2.1 水平扩散系数

由于水平湍流主要是由风切变引起的，因此水平湍流与流体的水平形变密切相关。
Smagorinsky (1963) 基于形变法（Deformation Method），给出了水平扩散系数的求解公式，
即假设次网格动能耗散主要由局地应变率张量（strain rate tensor）决定。

\[
   K_h = C_s \Delta^2 |\mathbb{D}|
\]

其中：

- \( C_s \) 为 Smagorinsky 常数（通常取 **0.1 - 0.25**，大气应用中常用 **0.2**）
- \( \Delta \) 为网格尺度（通常定义为网格间距的几何平均：\( \Delta = (\Delta x \Delta y )^{1/2} \)）

\( \mathbb{D} \) 为 应变率张量，定义为
\[
\mathbb{D} =
\begin{bmatrix}
\frac{\partial u}{\partial x} & \frac{1}{2} \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) \\
\frac{1}{2} \left( \frac{\partial u}{\partial y} + \frac{\partial v}{\partial x} \right) & \frac{\partial v}{\partial y}
\end{bmatrix}
\]

\( |\mathbb{D}| \) 是应变率张量的某种范数（两个特征值之差），
用来度量形变的整体程度，定义为:
  \[
  |\mathbb{D}| = 
  [(b + c)^2 + (a - b)^2]^{1/2} = 
  [(\frac{\partial u}{\partial y} +
    \frac{\partial v}{\partial x})^2 +
   (\frac{\partial u}{\partial x} -
    \frac{\partial v}{\partial y})^2]^{1/2}
  \]

\( |\mathbb{D}| \) 可以看作是剪切应变（形变）和拉伸应变（膨胀、压缩）的共同作用.

---

在 Smagorinsky 公式的基础上，可以增加一个常数项，让模拟的浓度变得更加光滑（Anthes and Warner, 1978）。
\[
   K_h = K_0 + C_s \Delta^2 |\mathbb{D}|
\]
其中 \(K_0\) 的值为
\[
  K_0 = 3 \times 10^{-3} \frac{\Delta x \Delta y}{\Delta t}
\]

## 2.2 边界条件
水平扩散采用严格的定值边界条件( Dirichlet 边界条件)。
\[
  c(t, b) = c_b
\]
\(c_b\) 为边界处的固定值，通过并行通信获得或者外部文件读入。

## 2.3 数值方案

为了便于通量计算，水平扩散系数\(K_h\) 通常定义在 Arakawa C 网格，
对于一个网格单元，其四条边（v-stag 和 u-stag 网格）均需计算一个\(K_h\)值。

根据 Smagorinsky 公式，计算 \(K_h\) 涉及速度梯度，需要空间差分，
一般采用中值差分方案求速度梯度。

\[
(\frac{\partial u}{\partial x})_{i+1/2, j} = 
\frac{u_{i+1/2+1, j} - u_{i+1/2-1, j}} {2 \Delta x}
\]

\[
(\frac{\partial u}{\partial y})_{i+1/2, j} = 
\frac{u_{i+1/2, j+1} - u_{i+1/2, j-1}}{2 \Delta y}
\]

在 u-stag 网格上的 \(v\)（定义在 v-stag 网格）梯度时，需要首先将 
\(v\) 从 v-stag 网格插值到 mass 网格。为了简单，通常采用线性插值方法。

\[
v(i, j) = 
\frac{v(i, j+\frac{1}{2}) + v(i, j-\frac{1}{2})}{2}
\]

因此
\[
(\frac{\partial v}{\partial x})_{i+1/2, j} = 
\frac{ v(i+1, j) - v(i, j)}{\Delta x}
\]

\[
(\frac{\partial v}{\partial y})_{i+1/2, j} = 
\frac{v(i, j+1) - v(i, j)}{\Delta y}
\]

将速度梯度带入 Smagorinsky 公式，即可求得 \(K_h\)

---

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

其中
\(
flux_{i+1/2} = {Kx_{i+1/2} \times \Delta c^{\rho}_{i+1/2}}
\) 可以等效理解为通量。

注意扩散是从高浓度到低浓度，\(K_h\) 始终大于0，浓度梯度决定通量的方向。

## 2.4 单元测试

设空气密度和扩散系数为定值，
比如 \(\rho = 1.0\) 且 \(K_h = 1.0\)，
则水平扩散方程可以写作标准的热传导方程

\[
\frac{\partial c}{\partial t} = K_h \nabla^2 c
\]

采用狄利克雷条件 \(c_b = 0\)，
用正玄函数构造初始条件
\[
c(x, 0) = \sin\left(\frac{\pi x}{L}\right) > 0
\]

其解析解为
\[
c(x, t) = \sin\left(\frac{\pi x}{L}\right) e^{-\frac{\pi^2 K_h}{L^2} t}
\]

某个网格的平均值( \( x_0 \) 到 \( x_1 \) 之间)，使用**定积分的平均值公式**：

\[
c_{\text{avg}} = \frac{1}{x_1 - x_0} \int_{x_0}^{x_1} c(x, t) \, dx.
\]

可以得到
\[
c_{\text{avg}} =
  - e^{-\frac{\pi^2 K_h}{L^2} t} \cdot \frac{L}{\pi(x_1 - x_0)} 
  \left[ \cos\left(\frac{\pi x_1}{L}\right) - 
         \cos\left(\frac{\pi x_0}{L}\right)
  \right]
\]

用该公式验证数值方案的精度。


 > Smagorinsky J. General circulation experiments with the primitive equations: I. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164.
 > 
 > Anthes R A, Warner T T. Development of hydrodynamic models suitable for air pollution and other mesometerological studies[J]. Monthly Weather Review, 1978, 106(8): 1045-1078.
