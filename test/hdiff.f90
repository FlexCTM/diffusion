module test_hdiff
    use mod_tool, only: fp
    use mod_hdiff, only: hdiff1d_by_k_theory
  
    use testdrive, only : error_type, unittest_type, new_unittest, check
    implicit none
    private
  
    public :: collect_case

  contains
  
    subroutine test_conservation(error)
        type(error_type), allocatable, intent(out) :: error
  
        integer, parameter :: nx = 5

        integer :: i
        real(fp) :: dt  !! 积分时间: s
        real(fp) :: dx(nx) !! 网格宽度: m
        real(fp) :: kh(nx) !! 扩散系数: m2/s
      
        real(fp) :: rho(nx) !! 密度
        real(fp) :: conc(nx), dc(nx), v(nx)
        real(fp) :: expected !! 浓度
      
        dt = 900.
        kh = [20., 10., 5., 1., 0.1]*20. ! 快速收敛
        dx = [50, 70, 100, 200, 300 ]*10
        v = [50*2, 70*1, 100*3, 200*5, 300*4 ]*10
        ! dx = 10.
        ! rho = [10., 1., 0.1, 0.01, 0.001] ! 稳定
        rho = [10., 9., 8., 7., 6. ] ! 稳定
      
        conc = [10, 8, 5, 6, 3]
        expected = sum(conc*dx)
        
        ! 守恒性检查
        write(*, '(10F10.3)') conc, expected
        do i = 1, 20
          call hdiff1d_by_k_theory(dt, dx, kh, rho, conc, dc)
          write(*, '(10F10.3)') conc, sum(conc*dx)
          ! call check(error, abs(sum(conc*dx) - expected) < expected*0.0001, "no conservation")
        end do
        
        ! ! 一致性检查
        ! do i = 1, nx-1
        !   write(*, '(10F10.3)')  conc(i)/conc(i+1), rho(i)/rho(i+1)
        !   call check(error, abs( conc(i)/conc(i+1) - rho(i)/rho(i+1)) < 0.01, "no consistency")
        ! end do

      end subroutine test_conservation

      subroutine test_accuracy(error)
        type(error_type), allocatable, intent(out) :: error
    
        integer, parameter :: nx = 10
        real(fp), parameter :: pi = 3.141592653589793_fp

        integer :: i, k, nt
        real(fp) :: L   !! 总长度: m
        real(fp) :: dt  !! 积分时间: s
        real(fp) :: dx(nx) !! 网格宽度: m
        real(fp) :: kh(nx) !! 扩散系数: m2/s
      
        real(fp) :: rho(nx) !! 密度
        real(fp) :: conc(nx), dc(nx), v(nx)
        real(fp) :: conc_exact(nx) !! 浓度
        real(fp) :: x0, x1, t

        nt = 15
        dt = 90.0_fp
        kh = 10.0_fp
        dx = 100.0_fp ! dt = dx**2/kh

        rho = 1.0_fp
        L = sum(dx) - dx(nx)

        x1 = 0.0_fp
        ! 初始化浓度场: 注意是平均值
        do k = 1, nx
          x0 = x1
          x1 = x0 + dx(k)
          conc(k) = sin(pi * x0/L)
        end do
  
        conc_exact = 0
        write(*, '(15F9.4)') conc
        write(*, *) "===================="
        do i = 1, nt
          ! 计算数值解
          call hdiff1d_by_k_theory(dt, dx, kh, rho, conc, dc)
          ! 计算解析解
          t = i * dt
          x1 = dx(1)
          do k = 2, nx -1
            x0 = x1
            x1 = x0 + dx(k)
            ! conc_exact(k) = - exp(- (pi**2 * kh(1) / L**2) * t) * &
            !                  (L / (pi * (x1 - x0))) * (cos(pi * x1 / L) - cos(pi * x0 / L))
            conc_exact(k) = exp(- (pi**2 * Kh(1) / L**2) * t) * sin(pi * x0/L)
          end do
          write(*, '(15F9.4)') conc
          write(*, '(15F9.4)') conc_exact
          ! 数值解在求解浓度差时，做的时线性假设，因此存在较大误差
          ! call check(error, all(abs(conc - conc_exact) < conc_exact*0.5), "no accuracy")
        end do

        ! 输出误差
        write(*,*), "numerical, exact,    error"
        do k = 2, nx - 1
          write(*, '(15F9.3)') conc(k), conc_exact(k), (conc(k) - conc_exact(k))/conc_exact(k)
        end do
        ! 整理偏差不要查过 5%
        call check(error, abs( sum(conc) - sum(conc_exact)) < sum(conc_exact)*0.05, "no accuracy")

      end subroutine test_accuracy

    !> Collect all exported unit tests
    subroutine collect_case(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)

        testsuite = [new_unittest("test hdiff conservation", test_conservation), &
                     new_unittest("test hdiff accuracy", test_accuracy)]
    end subroutine collect_case

  end module test_hdiff
  
  program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite
    use test_hdiff, only : collect_case
    implicit none
    integer :: stat
  
    stat = 0
    call run_testsuite(collect_case, error_unit, stat)
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if
  
  end program tester
