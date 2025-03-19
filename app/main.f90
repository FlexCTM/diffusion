program main
  use mod_tool, only: fp
  use diffusion, only: vdiff_by_k_theory

  implicit none

  integer, parameter :: nz = 5

  integer :: i
  real(fp) :: dt  !! 积分时间: s
  real(fp) :: kz(nz) !! 扩散系数: m2/s
  real(fp) :: dz(nz) !! 垂直网格宽度: m

  real(fp) :: rho(nz) !! 密度
  real(fp) :: conc(nz) !! 浓度

  dt = 10.
  kz = 0.1
  dz = 1.0
  rho = 1.0
  conc = [10, 8, 5, 10, 3]

  write(*, *) '========== vdiff ==========='

  write(*, '(10F10.3)') conc
  do i = 1, 20
    call vdiff_by_k_theory(dt, kz, dz, rho, conc)
    write(*, '(10F10.3)') conc
  end do

end program main
