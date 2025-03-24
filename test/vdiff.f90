module test_vdiff
    use mod_tool, only: fp
    use mod_vdiff, only: vdiff_by_k_theory
  
    use testdrive, only : error_type, unittest_type, new_unittest, check
    implicit none
    private
  
    public :: collect_case

  contains
  
    subroutine test_conservation(error)
        type(error_type), allocatable, intent(out) :: error
    
        integer, parameter :: nz = 5

        integer :: i
        real(fp) :: dt  !! 积分时间: s
        real(fp) :: kz(nz) !! 扩散系数: m2/s
        real(fp) :: dz(nz) !! 垂直网格宽度: m
      
        real(fp) :: rho(nz) !! 密度
        real(fp) :: conc(nz) !! 浓度
        real(fp) :: expected !! 浓度
      
        dt = 900.
        kz = [20., 10., 5., 1., 0.1]*20. ! 快速收敛
        dz = [50, 70, 100, 200, 300 ]
        rho = [10., 9., 8., 7., 6. ] ! 稳定
      
        conc = [10, 8, 5, 6, 3]
        expected = sum(conc*dz)
        
        ! 守恒性检查
        write(*, '(10F10.3)') conc, expected
        do i = 1, 20
          call vdiff_by_k_theory(dt, dz, kz, rho, conc)
          write(*, '(10F10.3)') conc, sum(conc*dz)
          call check(error, abs(sum(conc*dz) - expected) < expected*0.0001, "no conservation")
        end do
        
        ! 一致性检查: 扩散到平衡态，检查才有意义
        write(*,*), "     conc,      rho"
        do i = 1, nz-1
          write(*, '(10F10.3)')  conc(i)/conc(i+1), rho(i)/rho(i+1)
          call check(error, abs(conc(i)/conc(i+1) - rho(i)/rho(i+1)) < 0.01, "no consistency")
        end do

      end subroutine test_conservation

      subroutine test_accuracy(error)
        type(error_type), allocatable, intent(out) :: error

        integer, parameter :: nz = 10
        real(fp), parameter :: pi = 3.141592653589793_fp

        integer :: i, k, nt
        real(fp) :: L   !! 总高度: m
        real(fp) :: dt  !! 积分时间: s
        real(fp) :: kz(nz) !! 扩散系数: m2/s
        real(fp) :: dz(nz) !! 垂直网格宽度: m
      
        real(fp) :: rho(nz) !! 密度
        real(fp) :: conc(nz), conc_exact(nz) !! 浓度

        real(fp) :: z0, z1, t
      
        nt = 15
        dt = 5.0_fp
        kz = 1.0_fp
        dz = 1.0_fp
        rho = 1.0_fp
        L = dz(1)*nz

        z1 = 0.0_fp
        ! 初始化浓度场: 注意是平均值
        do k = 1, nz
          z0 = z1
          z1 = z0 + dz(k)
          ! conc(k) = cos(pi * (k-1+0.5) * dz(1) / L) + 1
          conc(k) = (L / (pi * (z1 - z0))) * (sin(pi*z1/L) - sin(pi*z0/L)) + 1.0_fp
        end do
  
        do i = 1, nt
          ! 计算数值解
          call vdiff_by_k_theory(dt, dz, kz, rho, conc)
          ! 计算解析解
          t = i * dt
          z1 = 0.0_fp
          do k = 1, nz
            z0 = z1
            z1 = z0 + dz(k)
            conc_exact(k) = exp(- (pi**2 * Kz(1) / L**2) * t) * &
            (L / (pi * (z1 - z0))) * (sin(pi * z1 / L) - sin(pi * z0 / L)) + 1.0_fp
          end do
          write(*, '(15F9.3)') conc, sum(conc*dz)
          write(*, '(15F9.3)') conc_exact, sum(conc_exact*dz)
          ! 数值解在求解浓度差时，做的时线性假设，因此存在较大误差
          call check(error, all(abs( conc - conc_exact) < conc_exact*0.5), "no accuracy")
        end do

        ! 输出误差
        write(*,*), "numerical,  exact,    error"
        do k = 1, nz
          write(*, '(15F9.3)') conc(k), conc_exact(k), (conc(k) - conc_exact(k))/conc_exact(k)
        end do
        ! 整理偏差不要查过 5%
        call check(error, abs( sum(conc) - sum(conc_exact)) < sum(conc_exact)*0.05, "no accuracy")

      end subroutine test_accuracy

    !> Collect all exported unit tests
    subroutine collect_case(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [new_unittest("test vdiff conservation", test_conservation), &
                     new_unittest("test vdiff accuracy", test_accuracy)]
    end subroutine collect_case

  end module test_vdiff
  
  program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite
    use test_vdiff, only : collect_case
    implicit none
    integer :: stat
  
    stat = 0
    call run_testsuite(collect_case, error_unit, stat)
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if
  
  end program tester
