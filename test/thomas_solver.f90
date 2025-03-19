module test_case
    use mod_tool, only: fp
    use mod_vdiff, only: thomas_solver
  
    use testdrive, only : error_type, unittest_type, new_unittest, check
    implicit none
    private
  
    public :: collect_case

  contains
  
    function diag_tridiag(n, a, b, c) result(AA)
        integer, intent(in) :: n
        real(fp), intent(in) :: a(n), b(n), c(n)
        real(fp) :: AA(n, n)
        integer :: i
        AA = 0.0_fp
        do i = 1, n
          AA(i, i) = b(i)
          if (i > 1) AA(i, i-1) = a(i)
          if (i < n) AA(i, i+1) = c(i)
        end do
    end function diag_tridiag

    !> Check substitution of a single line
    subroutine test_thomas1(error)
      !> Error handling
      type(error_type), allocatable, intent(out) :: error
      integer, parameter :: n = 4
      real(fp) :: a(n), b(n), c(n), x(n), d(n), x_expected(n)
    
      ! 定义三对角矩阵 A
      a = [0.0_fp,  -1.0_fp, -1.0_fp, -1.0_fp]  ! 下对角线 (a(1) 不使用)
      b = [2.0_fp,   2.0_fp,  2.0_fp,  2.0_fp]  ! 主对角线
      c = [-1.0_fp, -1.0_fp, -1.0_fp,  0.0_fp]  ! 上对角线 (c(n) 不使用)

      x_expected = [1.0_8, 2.0_8, 3.0_8, 4.0_8]
      
      d = matmul(diag_tridiag(n, a, b, c), x_expected) ! d = Ax

      x = d
      call thomas_solver(n, a, b, c, x)
  
      call check(error, all(abs(x - x_expected) < 1.0e-6_8), "Test Case 1 (n=4): Solution matches expected values")

    end subroutine test_thomas1

    subroutine test_thomas2(error)
        type(error_type), allocatable, intent(out) :: error
        integer, parameter :: n = 8
        real(fp) :: a(n), b(n), c(n), x(n), d(n), x_expected(n)
    
        a = [0.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8]
        b = [3.0_8, 3.0_8, 3.0_8, 3.0_8, 3.0_8, 3.0_8, 3.0_8, 3.0_8]
        c = [-1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, -1.0_8, 0.0_8]
    
        x_expected = [1.0_8, 2.0_8, 3.0_8, 4.0_8, 5.0_8, 6.0_8, 7.0_8, 8.0_8]
        d = matmul(diag_tridiag(n, a, b, c), x_expected)
        x = d
        call thomas_solver(n, a, b, c, x)
        call check(error, all(abs(x - x_expected) < 1.0e-6_8), "Test Case 2 (n=8): Solution matches expected values")
      end subroutine test_thomas2

    !> Collect all exported unit tests
    subroutine collect_case(testsuite)
        !> Collection of tests
        type(unittest_type), allocatable, intent(out) :: testsuite(:)
    
        testsuite = [new_unittest("thomas solver", test_thomas1), &
                     new_unittest("thomas solver", test_thomas2)]
    end subroutine collect_case

  end module test_case
  
  program tester
    use, intrinsic :: iso_fortran_env, only : error_unit
    use testdrive, only : run_testsuite
    use test_case, only : collect_case
    implicit none
    integer :: stat
  
    stat = 0
    call run_testsuite(collect_case, error_unit, stat)
  
    if (stat > 0) then
      write(error_unit, '(i0, 1x, a)') stat, "test(s) failed!"
      error stop
    end if
  
  end program tester
