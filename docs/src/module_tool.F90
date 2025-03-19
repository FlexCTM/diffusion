module mod_tool
   !! 工具集
   implicit none

   private

   public cal_cfl_time_step, fp, eps

#ifdef REAL_KIND
   integer, parameter :: fp = REAL_KIND
#else
   integer, parameter :: fp = 8
   ! use, intrinsic :: iso_fortran_env, only: fp => real64
#endif

   real(fp), parameter :: eps = epsilon(1.0_fp) !! 极小值

contains
   
   subroutine cal_cfl_time_step(dx, kx, dt, dt_, nt)
      !! 计算满足 CFL(Courant–Friedrichs–Lewy) 条件的时间积分步长
      real(fp), intent(in)  :: dx(:) !! 网格分辨率: m
      real(fp), intent(in)  :: kx(:) !! 扩散系数: m2/s
      real(fp), intent(in)  :: dt    !! 全局迭代步长: s
      real(fp), intent(out) :: dt_   !! 局部迭代步长: s
      integer, intent(out)  :: nt    !! 迭代多少次

      ! local variable
      real(fp), parameter :: factor = 0.8_fp !! 满足CFL条件的程度 < 1.0

      ! 计算 CFL 约束的最小时间步长
      dt_ = factor * minval(dx**2/kx)
   
      ! 确保 dt_ 能整除 dt，并计算迭代次数
      if (dt_ < dt) then
         nt = int(dt/dt_) + 1
         dt_ = dt/nt
      else
         dt_ = dt
         nt  = 1
      end if
   end subroutine cal_cfl_time_step

end module mod_tool
