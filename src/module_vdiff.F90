module mod_vdiff
   !! 垂直扩散
   use mod_tool, only: fp, eps

   implicit none

   private
   public vdiff_by_k_theory, thomas_solver

   contains

   subroutine vdiff_by_k_theory(dt, kz, dz, rho, conc)
      !! 采用后向欧拉法求解垂直方向的网格湍流扩散（K-theory）。
      !! 用 Thomas 算法求解后向时间差分格式的垂直扩散方程。
      !! 假设底层和顶层边界通量均为 0。
      real(fp), intent(in) :: dt    !! 积分时间: s
      real(fp), intent(in) :: kz(:) !! 垂直扩散系数：m^2/s
      real(fp), intent(in) :: dz(:) !! 网格垂直层宽度: m

      real(fp), intent(in) :: rho(:) !! 空气密度: kg/m^3
      real(fp), intent(inout) :: conc(:) !! 浓度: ug/m^3

      real(fp) :: aa(size(conc)) !! 下对角线
      real(fp) :: bb(size(conc)) !! 主对角线
      real(fp) :: cc(size(conc)) !! 上对角线
      real(fp) :: kz_rho(size(conc)) !! 辅助变量
      integer :: k, nz

      nz = size(dz)

      do k = 1, nz-1
         kz_rho(k) = kz(k) * (rho(k) + rho(k+1))/(dz(k) + dz(k+1))
      end do
      ! 后向欧拉
      ! Lower boundary condition
      aa(1) = 0.
      bb(1) = 1. + dt/dz(1) * kz_rho(1)/rho(1)
      ! Upper boundary condition
      bb(nz) = 1. + dt/dz(nz) * kz_rho(nz-1)/rho(nz)
      cc(nz) = 0.

      ! middle
      do k = 2, nz
         aa(k) = - dt/dz(k) * kz_rho(k-1)/rho(k-1)
      end do

      do k = 2, nz-1
         bb(k) = 1. + dt/dz(k) * (kz_rho(k-1)/rho(k) + kz_rho(k)/rho(k))
      end do

      do k = 1, nz-1
         cc(k) = - dt/dz(k) * kz_rho(k)/rho(k+1)
      end do

      ! A c_{t+1} = c_t
      call thomas_solver(nz, aa, bb, cc, conc)

      where(conc < eps) conc = 0.

   end subroutine vdiff_by_k_theory

   subroutine thomas_solver(n, a, b, c, x)
      !! Thomas算法求解三对角矩阵线性方程组 Ax = d，
      !! 其中，A = [a, b, c]
      implicit none
      integer, intent(in) :: n     !! 矩阵规模
      real(fp), intent(in) :: a(n) !! 下对角线
      real(fp), intent(in) :: b(n) !! 主对角线
      real(fp), intent(in) :: c(n) !! 上对角线
      real(fp), intent(inout) :: x(n) !! 求解得到的解向量

      real(fp) :: b_(n) !! 临时数组用于存储修正后的对角线
      integer :: i
      ! 前向消去, 消去下对角线(a)
      b_(1) = b(1)
      x(1) = x(1)/b_(1) ! 最开始的时候 d = x
      do i = 2, n
         b_(i) = b(i) - a(i)/b_(i-1) * c(i-1) ! c(i-1) ==> b(i)
         x(i) = (x(i) - a(i)*x(i-1))/b_(i) ! x(i-1) ==> x(i)
      end do

      ! 回带求解
      x(n) = x(n)
      do i = n-1, 1, -1
         x(i) = x(i) - c(i) * x(i+1) / (b_(i))
      end do
   end subroutine thomas_solver

end module mod_vdiff
