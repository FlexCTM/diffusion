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

  dt = 100.
  kz = 100.
  dz = [50, 70, 100, 200, 300 ]
  ! rho = [10., 1., 0.1, 0.01, 0.001] ! 稳定
  rho = [10., 9., 8., 7., 6. ] ! 稳定

  conc = [10, 8, 5, 10, 3]

  write(*, *) '========== vdiff ==========='

  write(*, '(10F10.3)') conc, sum(conc)
  do i = 1, 20
    call vdiff_by_k_theory(dt, kz, dz, rho, conc)
    write(*, '(10F10.3)') conc, sum(conc)
  end do

end program main
