module mod_hdiff
   !! 水平扩散
   use mod_tool, only: fp, eps
   implicit none

   private
   public cal_kh_by_deformation_method, hdiff1d_by_k_theory

   contains
   
   subroutine cal_kh_by_deformation_method(dt, dx, dy, u, v, kx, ky, area)
      !! 计算次网格湍流扩散的水平扩散系数 in the Arakawa C grid。
      !! 采用形变方法计算水平扩散系数，参考 Smagorinsky (1963)。
      !! 
      !! Smagorinsky J. General circulation experiments with the primitive equations: 
      !! I. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164.
      !! 
      real(fp), intent(in) :: dt       !! 积分时间: s
      real(fp), intent(in) :: dx(:, :) !! x方向长度: m
      real(fp), intent(in) :: dy(:, :) !! y方向长度: m
      real(fp), intent(in) :: u(:, :)  !! 纬向风 u-stag: m/s
      real(fp), intent(in) :: v(:, :)  !! 经向风 v-stag: m/s

      real(fp), intent(out) :: kx(:, :) !! horizontal diffusion coefficient in x: m^2/s
      real(fp), intent(out) :: ky(:, :) !! horizontal diffusion coefficient in y: m^2/s
      real(fp), optional, intent(in) :: area !! 传入全局平均面积，避免重复计算: m^2

      ! local
      real(fp), parameter :: C_s = 1._fp/(4._fp*sqrt(2._fp)) !! Smagorinsky常数, 1/(4*sqrt(2)) = 0.18
      real(fp), parameter :: kh_min = 1.0_fp !! 最小扩散系数: m^2/s
      real(fp) :: kh_max !! 最大扩散系数: m^2/s
      real(fp) :: kh0 !! 静态扩散系数 (Anthes and Warner, 1978)：m^2/s

      real(fp) :: cell_area !! 网格面积：m^2
      real(fp) :: dudx, dudy, dvdx, dvdy !! 二维速度梯度张量（形变+旋转）的风量: s^-1

      integer :: i, j

      if (present(area)) then
         cell_area = area
         kh_max = max(area/(32._fp*dt), 1.e5_fp)
      end if
      do j = 2, size(u, 2)-1
         do i = 2, size(u, 1)-1
            if (.not. present(area)) then
               cell_area = dx(i, j)*dy(i, j)
               kh_max = max(cell_area/(32._fp*dt), 1.e5_fp)
            end if
            kh0 = 0.003_fp*cell_area/dt !! 可以让评估更加平滑
            ! Kx at u-stag cell
            dudx = ( u(i+1, j) - u(i-1, j) )/(2*dx(i, j))
            dudy = ( u(i, j+1) - u(i, j-1) )/(2*dy(i, j))
            dvdx = ((v(i+1, j) - v(i, j)) + (v(i+1, j-1) - v(i, j-1)))/(2*dx(i, j))
            dvdy = ((v(i, j) - v(i, j-1)) + (v(i+1, j) - v(i+1, j-1)))/(2*dy(i, j))
            kx(i, j) = kh0 + C_s*cell_area*sqrt((dudy+dvdx)**2 + (dudx-dvdy)**2)

            ! Ky at v-stag cell
            dudx = ((u(i, j) - u(i-1, j)) + (u(i, j+1) - u(i-1, j+1)))/(2*dx(i, j))
            dudy = ((u(i, j+1) - u(i, j)) + (u(i-1, j+1) - u(i-1, j)))/(2*dy(i, j))
            dvdx = (v(i+1, j) - v(i-1, j))/(2*dx(i, j))
            dvdy = (v(i, j+1) - v(i, j-1))/(2*dy(i, j))
            ky(i, j) = kh0 + C_s*cell_area*sqrt((dudy+dvdx)**2 + (dudx-dvdy)**2)

            ! 确保数值计算稳定性
            kx(i, j) = max(kh_min, min(kx(i, j), kh_max))
            ky(i, j) = max(kh_min, min(ky(i, j), kh_max))
         end do
      end do
   end subroutine cal_kh_by_deformation_method

   subroutine hdiff1d_by_k_theory(dt, dx, kh, rho, c, dc, volume)
      !! 采用前向欧拉法求解 1D 次网格湍流扩散（K-theory），最外 1 圈为边界，不做更新。
      !! 假设密度的扰动相对于平均密度很小，可以忽略。
      !! 采用体积比计算扩散，确保污染物浓度分布与密度分布的一致性。
      real(fp), intent(in) :: dt     !! 积分时间: s
      real(fp), intent(in) :: dx(:)  !! 网格长度: m
      real(fp), intent(in) :: kh(:)  !! 扩散系数: m^2/s
      real(fp), intent(in) :: rho(:) !! 密度: kg/m^3

      real(fp), intent(inout) :: c(:)  !! 浓度: ug/m^3
      real(fp), intent(out)   :: dc(:) !! 浓度变化: ug/m^3
      real(fp), optional, intent(in) :: volume(:) !! 网格体积: m^3

      ! local
      real(fp) :: u      !! 虚拟速度，扩散等效为速度: m/s
      real(fp) :: fR, fL !! 某个网格右侧和左侧通量: ug/m^2/s
      real(fp) :: flux(size(c)) !! 所有网格的右侧通量: ug/m^2/s
      integer :: i

      dc = 0.0_fp
      ! 计算通量
      do i = 1, size(c) - 1 ! flux(i) > 为流出入
         u = kh(i)*2._fp/(dx(i) + dx(i+1))
         flux(i) = (rho(i) + rho(i+1))/2.0 * u * (c(i+1)/rho(i+1) - c(i)/rho(i)) 
      end do
      ! 时间积分, 前向欧拉
      do i = 2, size(c) - 1
         fR = flux(i)
         fL = flux(i-1)
         if (present(volume)) then ! 进行校正体积校正，保证质量守恒
            if (fL < 0) fL = fL*volume(i-1)/volume(i)
            if (fR > 0) fR = fR*volume(i+1)/volume(i)
         end if
         dc(i) = dt/dx(i) * (fR - fL)
         c(i)  = c(i) + dc(i)
         if(c(i) < eps) c(i) = 0. ! 数值稳定性
      end do

   end subroutine hdiff1d_by_k_theory

end module mod_hdiff
