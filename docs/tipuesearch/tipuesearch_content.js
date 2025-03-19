var tipuesearch = {"pages":[{"title":" diffusion ","text":"diffusion 扩散模块 Developer Info Linhong Xiao","tags":"home","loc":"index.html"},{"title":"cal_cfl_time_step – diffusion","text":"public  subroutine cal_cfl_time_step(dx, kx, dt, dt_, nt) 计算满足 CFL(Courant–Friedrichs–Lewy) 条件的时间积分步长 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dx (:) 网格分辨率: m real(kind=fp), intent(in) :: kx (:) 扩散系数: m2/s real(kind=fp), intent(in) :: dt 全局迭代步长: s real(kind=fp), intent(out) :: dt_ 局部迭代步长: s integer, intent(out) :: nt 迭代多少次","tags":"","loc":"proc/cal_cfl_time_step.html"},{"title":"vdiff_by_k_theory – diffusion","text":"public  subroutine vdiff_by_k_theory(dt, kz, dz, rho, conc) 采用后向欧拉法求解垂直方向的网格湍流扩散（K-theory）。\n用 Thomas 算法求解后向时间差分格式的垂直扩散方程。\n假设底层和顶层边界通量均为 0。 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: kz (:) 垂直扩散系数：m&#94;2/s real(kind=fp), intent(in) :: dz (:) 网格垂直层宽度: m real(kind=fp), intent(in) :: rho (:) 空气密度: kg/m&#94;3 real(kind=fp), intent(inout) :: conc (:) 浓度: ug/m&#94;3 Calls proc~~vdiff_by_k_theory~~CallsGraph proc~vdiff_by_k_theory vdiff_by_k_theory proc~thomas_solver thomas_solver proc~vdiff_by_k_theory->proc~thomas_solver Help Graph Key Nodes of different colours represent the following: Graph Key Subroutine Subroutine Function Function Interface Interface Type Bound Procedure Type Bound Procedure Unknown Procedure Type Unknown Procedure Type Program Program This Page's Entity This Page's Entity Solid arrows point from a procedure to one which it calls. Dashed \narrows point from an interface to procedures which implement that interface.\nThis could include the module procedures in a generic interface or the\nimplementation in a submodule of an interface in a parent module.","tags":"","loc":"proc/vdiff_by_k_theory.html"},{"title":"thomas_solver – diffusion","text":"public  subroutine thomas_solver(n, a, b, c, x) Thomas算法求解三对角矩阵线性方程组 Ax = d，\n其中，A = [a, b, c] Arguments Type Intent Optional Attributes Name integer, intent(in) :: n 矩阵规模 real(kind=fp), intent(in) :: a (n) 下对角线 real(kind=fp), intent(in) :: b (n) 主对角线 real(kind=fp), intent(in) :: c (n) 上对角线 real(kind=fp), intent(inout) :: x (n) 求解得到的解向量 Called by proc~~thomas_solver~~CalledByGraph proc~thomas_solver thomas_solver proc~vdiff_by_k_theory vdiff_by_k_theory proc~vdiff_by_k_theory->proc~thomas_solver Help Graph Key Nodes of different colours represent the following: Graph Key Subroutine Subroutine Function Function Interface Interface Type Bound Procedure Type Bound Procedure Unknown Procedure Type Unknown Procedure Type Program Program This Page's Entity This Page's Entity Solid arrows point from a procedure to one which it calls. Dashed \narrows point from an interface to procedures which implement that interface.\nThis could include the module procedures in a generic interface or the\nimplementation in a submodule of an interface in a parent module.","tags":"","loc":"proc/thomas_solver.html"},{"title":"cal_kh_by_deformation_method – diffusion","text":"public  subroutine cal_kh_by_deformation_method(dt, dx, dy, u, v, kx, ky, area) 计算次网格湍流扩散的水平扩散系数 in the Arakawa C grid。\n采用形变方法计算水平扩散系数，参考 Smagorinsky (1963)。 Smagorinsky J. General circulation experiments with the primitive equations:\nI. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164. 可以让评估更加平滑 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: dx (:,:) x方向长度: m real(kind=fp), intent(in) :: dy (:,:) y方向长度: m real(kind=fp), intent(in) :: u (:,:) 纬向风 u-stag: m/s real(kind=fp), intent(in) :: v (:,:) 经向风 v-stag: m/s real(kind=fp), intent(out) :: kx (:,:) horizontal diffusion coefficient in x: m&#94;2/s real(kind=fp), intent(out) :: ky (:,:) horizontal diffusion coefficient in y: m&#94;2/s real(kind=fp), intent(in), optional :: area 传入全局平均面积，避免重复计算: m&#94;2","tags":"","loc":"proc/cal_kh_by_deformation_method.html"},{"title":"hdiff1d_by_k_theory – diffusion","text":"public  subroutine hdiff1d_by_k_theory(dt, dx, kh, rho, c, dc, volume) 采用前向欧拉法求解 1D 次网格湍流扩散（K-theory），最外 1 圈为边界，不做更新。\n假设密度的扰动相对于平均密度很小，可以忽略。\n采用体积比计算扩散，确保污染物浓度分布与密度分布的一致性。 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: dx (:) 网格长度: m real(kind=fp), intent(in) :: kh (:) 扩散系数: m&#94;2/s real(kind=fp), intent(in) :: rho (:) 密度: kg/m&#94;3 real(kind=fp), intent(inout) :: c (:) 浓度: ug/m&#94;3 real(kind=fp), intent(out) :: dc (:) 浓度变化: ug/m&#94;3 real(kind=fp), intent(in), optional :: volume (:) 网格体积: m&#94;3","tags":"","loc":"proc/hdiff1d_by_k_theory.html"},{"title":"mod_tool – diffusion","text":"工具集 Used by module~~mod_tool~~UsedByGraph module~mod_tool mod_tool module~diffusion diffusion module~diffusion->module~mod_tool module~mod_hdiff mod_hdiff module~diffusion->module~mod_hdiff module~mod_vdiff mod_vdiff module~diffusion->module~mod_vdiff module~mod_hdiff->module~mod_tool module~mod_vdiff->module~mod_tool Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses. Variables Type Visibility Attributes Name Initial integer, public, parameter :: fp = 8 real(kind=fp), public, parameter :: eps = epsilon(1.0_fp) 极小值 Subroutines public  subroutine cal_cfl_time_step (dx, kx, dt, dt_, nt) 计算满足 CFL(Courant–Friedrichs–Lewy) 条件的时间积分步长 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dx (:) 网格分辨率: m real(kind=fp), intent(in) :: kx (:) 扩散系数: m2/s real(kind=fp), intent(in) :: dt 全局迭代步长: s real(kind=fp), intent(out) :: dt_ 局部迭代步长: s integer, intent(out) :: nt 迭代多少次","tags":"","loc":"module/mod_tool.html"},{"title":"diffusion – diffusion","text":"Uses mod_vdiff mod_hdiff mod_tool module~~diffusion~~UsesGraph module~diffusion diffusion module~mod_hdiff mod_hdiff module~diffusion->module~mod_hdiff module~mod_tool mod_tool module~diffusion->module~mod_tool module~mod_vdiff mod_vdiff module~diffusion->module~mod_vdiff module~mod_hdiff->module~mod_tool module~mod_vdiff->module~mod_tool Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses.","tags":"","loc":"module/diffusion.html"},{"title":"mod_vdiff – diffusion","text":"垂直扩散 Uses mod_tool module~~mod_vdiff~~UsesGraph module~mod_vdiff mod_vdiff module~mod_tool mod_tool module~mod_vdiff->module~mod_tool Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses. Used by module~~mod_vdiff~~UsedByGraph module~mod_vdiff mod_vdiff module~diffusion diffusion module~diffusion->module~mod_vdiff Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses. Subroutines public  subroutine vdiff_by_k_theory (dt, kz, dz, rho, conc) 采用后向欧拉法求解垂直方向的网格湍流扩散（K-theory）。\n用 Thomas 算法求解后向时间差分格式的垂直扩散方程。\n假设底层和顶层边界通量均为 0。 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: kz (:) 垂直扩散系数：m&#94;2/s real(kind=fp), intent(in) :: dz (:) 网格垂直层宽度: m real(kind=fp), intent(in) :: rho (:) 空气密度: kg/m&#94;3 real(kind=fp), intent(inout) :: conc (:) 浓度: ug/m&#94;3 public  subroutine thomas_solver (n, a, b, c, x) Thomas算法求解三对角矩阵线性方程组 Ax = d，\n其中，A = [a, b, c] Arguments Type Intent Optional Attributes Name integer, intent(in) :: n 矩阵规模 real(kind=fp), intent(in) :: a (n) 下对角线 real(kind=fp), intent(in) :: b (n) 主对角线 real(kind=fp), intent(in) :: c (n) 上对角线 real(kind=fp), intent(inout) :: x (n) 求解得到的解向量","tags":"","loc":"module/mod_vdiff.html"},{"title":"mod_hdiff – diffusion","text":"水平扩散 Uses mod_tool module~~mod_hdiff~~UsesGraph module~mod_hdiff mod_hdiff module~mod_tool mod_tool module~mod_hdiff->module~mod_tool Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses. Used by module~~mod_hdiff~~UsedByGraph module~mod_hdiff mod_hdiff module~diffusion diffusion module~diffusion->module~mod_hdiff Help Graph Key Nodes of different colours represent the following: Graph Key Module Module Submodule Submodule Subroutine Subroutine Function Function Program Program This Page's Entity This Page's Entity Solid arrows point from a submodule to the (sub)module which it is\ndescended from. Dashed arrows point from a module or program unit to \nmodules which it uses. Subroutines public  subroutine cal_kh_by_deformation_method (dt, dx, dy, u, v, kx, ky, area) 计算次网格湍流扩散的水平扩散系数 in the Arakawa C grid。\n采用形变方法计算水平扩散系数，参考 Smagorinsky (1963)。 Read more… Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: dx (:,:) x方向长度: m real(kind=fp), intent(in) :: dy (:,:) y方向长度: m real(kind=fp), intent(in) :: u (:,:) 纬向风 u-stag: m/s real(kind=fp), intent(in) :: v (:,:) 经向风 v-stag: m/s real(kind=fp), intent(out) :: kx (:,:) horizontal diffusion coefficient in x: m&#94;2/s real(kind=fp), intent(out) :: ky (:,:) horizontal diffusion coefficient in y: m&#94;2/s real(kind=fp), intent(in), optional :: area 传入全局平均面积，避免重复计算: m&#94;2 public  subroutine hdiff1d_by_k_theory (dt, dx, kh, rho, c, dc, volume) 采用前向欧拉法求解 1D 次网格湍流扩散（K-theory），最外 1 圈为边界，不做更新。\n假设密度的扰动相对于平均密度很小，可以忽略。\n采用体积比计算扩散，确保污染物浓度分布与密度分布的一致性。 Arguments Type Intent Optional Attributes Name real(kind=fp), intent(in) :: dt 积分时间: s real(kind=fp), intent(in) :: dx (:) 网格长度: m real(kind=fp), intent(in) :: kh (:) 扩散系数: m&#94;2/s real(kind=fp), intent(in) :: rho (:) 密度: kg/m&#94;3 real(kind=fp), intent(inout) :: c (:) 浓度: ug/m&#94;3 real(kind=fp), intent(out) :: dc (:) 浓度变化: ug/m&#94;3 real(kind=fp), intent(in), optional :: volume (:) 网格体积: m&#94;3","tags":"","loc":"module/mod_hdiff.html"},{"title":"module_tool.F90 – diffusion","text":"Files dependent on this one sourcefile~~module_tool.f90~~AfferentGraph sourcefile~module_tool.f90 module_tool.F90 sourcefile~diffusion.f90 diffusion.F90 sourcefile~diffusion.f90->sourcefile~module_tool.f90 sourcefile~module_hdiff.f90 module_hdiff.F90 sourcefile~diffusion.f90->sourcefile~module_hdiff.f90 sourcefile~module_vdiff.f90 module_vdiff.F90 sourcefile~diffusion.f90->sourcefile~module_vdiff.f90 sourcefile~module_hdiff.f90->sourcefile~module_tool.f90 sourcefile~module_vdiff.f90->sourcefile~module_tool.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Source Code module mod_tool !! 工具集 implicit none private public cal_cfl_time_step , fp , eps #ifdef REAL_KIND integer , parameter :: fp = REAL_KIND #else integer , parameter :: fp = 8 ! use, intrinsic :: iso_fortran_env, only: fp => real64 #endif real ( fp ), parameter :: eps = epsilon ( 1.0_fp ) !! 极小值 contains subroutine cal_cfl_time_step ( dx , kx , dt , dt_ , nt ) !! 计算满足 CFL(Courant–Friedrichs–Lewy) 条件的时间积分步长 real ( fp ), intent ( in ) :: dx (:) !! 网格分辨率: m real ( fp ), intent ( in ) :: kx (:) !! 扩散系数: m2/s real ( fp ), intent ( in ) :: dt !! 全局迭代步长: s real ( fp ), intent ( out ) :: dt_ !! 局部迭代步长: s integer , intent ( out ) :: nt !! 迭代多少次 ! local variable real ( fp ), parameter :: factor = 0.8_fp !! 满足CFL条件的程度 < 1.0 ! 计算 CFL 约束的最小时间步长 dt_ = factor * minval ( dx ** 2 / kx ) ! 确保 dt_ 能整除 dt，并计算迭代次数 if ( dt_ < dt ) then nt = int ( dt / dt_ ) + 1 dt_ = dt / nt else dt_ = dt nt = 1 end if end subroutine cal_cfl_time_step end module mod_tool","tags":"","loc":"sourcefile/module_tool.f90.html"},{"title":"diffusion.F90 – diffusion","text":"This file depends on sourcefile~~diffusion.f90~~EfferentGraph sourcefile~diffusion.f90 diffusion.F90 sourcefile~module_hdiff.f90 module_hdiff.F90 sourcefile~diffusion.f90->sourcefile~module_hdiff.f90 sourcefile~module_tool.f90 module_tool.F90 sourcefile~diffusion.f90->sourcefile~module_tool.f90 sourcefile~module_vdiff.f90 module_vdiff.F90 sourcefile~diffusion.f90->sourcefile~module_vdiff.f90 sourcefile~module_hdiff.f90->sourcefile~module_tool.f90 sourcefile~module_vdiff.f90->sourcefile~module_tool.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Source Code module diffusion use mod_tool , only : cal_cfl_time_step use mod_vdiff , only : vdiff_by_k_theory use mod_hdiff , only : cal_kh_by_deformation_method , hdiff1d_by_k_theory end module diffusion","tags":"","loc":"sourcefile/diffusion.f90.html"},{"title":"module_vdiff.F90 – diffusion","text":"This file depends on sourcefile~~module_vdiff.f90~~EfferentGraph sourcefile~module_vdiff.f90 module_vdiff.F90 sourcefile~module_tool.f90 module_tool.F90 sourcefile~module_vdiff.f90->sourcefile~module_tool.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Files dependent on this one sourcefile~~module_vdiff.f90~~AfferentGraph sourcefile~module_vdiff.f90 module_vdiff.F90 sourcefile~diffusion.f90 diffusion.F90 sourcefile~diffusion.f90->sourcefile~module_vdiff.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Source Code module mod_vdiff !! 垂直扩散 use mod_tool , only : fp , eps implicit none private public vdiff_by_k_theory , thomas_solver contains subroutine vdiff_by_k_theory ( dt , kz , dz , rho , conc ) !! 采用后向欧拉法求解垂直方向的网格湍流扩散（K-theory）。 !! 用 Thomas 算法求解后向时间差分格式的垂直扩散方程。 !! 假设底层和顶层边界通量均为 0。 real ( fp ), intent ( in ) :: dt !! 积分时间: s real ( fp ), intent ( in ) :: kz (:) !! 垂直扩散系数：m&#94;2/s real ( fp ), intent ( in ) :: dz (:) !! 网格垂直层宽度: m real ( fp ), intent ( in ) :: rho (:) !! 空气密度: kg/m&#94;3 real ( fp ), intent ( inout ) :: conc (:) !! 浓度: ug/m&#94;3 real ( fp ) :: aa ( size ( conc )) !! 下对角线 real ( fp ) :: bb ( size ( conc )) !! 主对角线 real ( fp ) :: cc ( size ( conc )) !! 上对角线 real ( fp ) :: kz_rho ( size ( conc )) !! 辅助变量 integer :: k , nz nz = size ( dz ) do k = 1 , nz - 1 kz_rho ( k ) = kz ( k ) * ( rho ( k ) + rho ( k + 1 )) / ( dz ( k ) + dz ( k + 1 )) end do ! 后向欧拉 ! Lower boundary condition aa ( 1 ) = 0. bb ( 1 ) = 1. + dt / dz ( 1 ) * kz_rho ( 1 ) / rho ( 1 ) ! Upper boundary condition bb ( nz ) = 1. + dt / dz ( nz ) * kz_rho ( nz - 1 ) / rho ( nz ) cc ( nz ) = 0. ! middle do k = 2 , nz aa ( k ) = - dt / dz ( k ) * kz_rho ( k - 1 ) / rho ( k - 1 ) end do do k = 2 , nz - 1 bb ( k ) = 1. + dt / dz ( k ) * ( kz_rho ( k - 1 ) / rho ( k ) + kz_rho ( k ) / rho ( k )) end do do k = 1 , nz - 1 cc ( k ) = - dt / dz ( k ) * kz_rho ( k ) / rho ( k + 1 ) end do ! A c_{t+1} = c_t call thomas_solver ( nz , aa , bb , cc , conc ) where ( conc < eps ) conc = 0. end subroutine vdiff_by_k_theory subroutine thomas_solver ( n , a , b , c , x ) !! Thomas算法求解三对角矩阵线性方程组 Ax = d， !! 其中，A = [a, b, c] implicit none integer , intent ( in ) :: n !! 矩阵规模 real ( fp ), intent ( in ) :: a ( n ) !! 下对角线 real ( fp ), intent ( in ) :: b ( n ) !! 主对角线 real ( fp ), intent ( in ) :: c ( n ) !! 上对角线 real ( fp ), intent ( inout ) :: x ( n ) !! 求解得到的解向量 real ( fp ) :: b_ ( n ) !! 临时数组用于存储修正后的对角线 integer :: i ! 前向消去, 消去下对角线(a) b_ ( 1 ) = b ( 1 ) x ( 1 ) = x ( 1 ) / b_ ( 1 ) ! 最开始的时候 d = x do i = 2 , n b_ ( i ) = b ( i ) - a ( i ) / b_ ( i - 1 ) * c ( i - 1 ) ! c(i-1) ==> b(i) x ( i ) = ( x ( i ) - a ( i ) * x ( i - 1 )) / b_ ( i ) ! x(i-1) ==> x(i) end do ! 回带求解 x ( n ) = x ( n ) do i = n - 1 , 1 , - 1 x ( i ) = x ( i ) - c ( i ) * x ( i + 1 ) / ( b_ ( i )) end do end subroutine thomas_solver end module mod_vdiff","tags":"","loc":"sourcefile/module_vdiff.f90.html"},{"title":"module_hdiff.F90 – diffusion","text":"This file depends on sourcefile~~module_hdiff.f90~~EfferentGraph sourcefile~module_hdiff.f90 module_hdiff.F90 sourcefile~module_tool.f90 module_tool.F90 sourcefile~module_hdiff.f90->sourcefile~module_tool.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Files dependent on this one sourcefile~~module_hdiff.f90~~AfferentGraph sourcefile~module_hdiff.f90 module_hdiff.F90 sourcefile~diffusion.f90 diffusion.F90 sourcefile~diffusion.f90->sourcefile~module_hdiff.f90 Help Graph Key Nodes of different colours represent the following: Graph Key Source File Source File This Page's Entity This Page's Entity Solid arrows point from a file to a file which it depends on. A file\nis dependent upon another if the latter must be compiled before the former\ncan be. Source Code module mod_hdiff !! 水平扩散 use mod_tool , only : fp , eps implicit none private public cal_kh_by_deformation_method , hdiff1d_by_k_theory contains subroutine cal_kh_by_deformation_method ( dt , dx , dy , u , v , kx , ky , area ) !! 计算次网格湍流扩散的水平扩散系数 in the Arakawa C grid。 !! 采用形变方法计算水平扩散系数，参考 Smagorinsky (1963)。 !! !! Smagorinsky J. General circulation experiments with the primitive equations: !! I. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164. !! real ( fp ), intent ( in ) :: dt !! 积分时间: s real ( fp ), intent ( in ) :: dx (:, :) !! x方向长度: m real ( fp ), intent ( in ) :: dy (:, :) !! y方向长度: m real ( fp ), intent ( in ) :: u (:, :) !! 纬向风 u-stag: m/s real ( fp ), intent ( in ) :: v (:, :) !! 经向风 v-stag: m/s real ( fp ), intent ( out ) :: kx (:, :) !! horizontal diffusion coefficient in x: m&#94;2/s real ( fp ), intent ( out ) :: ky (:, :) !! horizontal diffusion coefficient in y: m&#94;2/s real ( fp ), optional , intent ( in ) :: area !! 传入全局平均面积，避免重复计算: m&#94;2 ! local real ( fp ), parameter :: C_s = 1._fp / ( 4._fp * sqrt ( 2._fp )) !! Smagorinsky常数, 1/(4*sqrt(2)) = 0.18 real ( fp ), parameter :: kh_min = 1.0_fp !! 最小扩散系数: m&#94;2/s real ( fp ) :: kh_max !! 最大扩散系数: m&#94;2/s real ( fp ) :: kh0 !! 静态扩散系数 (Anthes and Warner, 1978)：m&#94;2/s real ( fp ) :: cell_area !! 网格面积：m&#94;2 real ( fp ) :: dudx , dudy , dvdx , dvdy !! 二维速度梯度张量（形变+旋转）的风量: s&#94;-1 integer :: i , j if ( present ( area )) then cell_area = area kh_max = max ( area / ( 3 2._fp * dt ), 1.e5_fp ) end if do j = 2 , size ( u , 2 ) - 1 do i = 2 , size ( u , 1 ) - 1 if (. not . present ( area )) then cell_area = dx ( i , j ) * dy ( i , j ) kh_max = max ( cell_area / ( 3 2._fp * dt ), 1.e5_fp ) end if kh0 = 0.003_fp * cell_area / dt !! 可以让评估更加平滑 ! Kx at u-stag cell dudx = ( u ( i + 1 , j ) - u ( i - 1 , j ) ) / ( 2 * dx ( i , j )) dudy = ( u ( i , j + 1 ) - u ( i , j - 1 ) ) / ( 2 * dy ( i , j )) dvdx = (( v ( i + 1 , j ) - v ( i , j )) + ( v ( i + 1 , j - 1 ) - v ( i , j - 1 ))) / ( 2 * dx ( i , j )) dvdy = (( v ( i , j ) - v ( i , j - 1 )) + ( v ( i + 1 , j ) - v ( i + 1 , j - 1 ))) / ( 2 * dy ( i , j )) kx ( i , j ) = kh0 + C_s * cell_area * sqrt (( dudy + dvdx ) ** 2 + ( dudx - dvdy ) ** 2 ) ! Ky at v-stag cell dudx = (( u ( i , j ) - u ( i - 1 , j )) + ( u ( i , j + 1 ) - u ( i - 1 , j + 1 ))) / ( 2 * dx ( i , j )) dudy = (( u ( i , j + 1 ) - u ( i , j )) + ( u ( i - 1 , j + 1 ) - u ( i - 1 , j ))) / ( 2 * dy ( i , j )) dvdx = ( v ( i + 1 , j ) - v ( i - 1 , j )) / ( 2 * dx ( i , j )) dvdy = ( v ( i , j + 1 ) - v ( i , j - 1 )) / ( 2 * dy ( i , j )) ky ( i , j ) = kh0 + C_s * cell_area * sqrt (( dudy + dvdx ) ** 2 + ( dudx - dvdy ) ** 2 ) ! 确保数值计算稳定性 kx ( i , j ) = max ( kh_min , min ( kx ( i , j ), kh_max )) ky ( i , j ) = max ( kh_min , min ( ky ( i , j ), kh_max )) end do end do end subroutine cal_kh_by_deformation_method subroutine hdiff1d_by_k_theory ( dt , dx , kh , rho , c , dc , volume ) !! 采用前向欧拉法求解 1D 次网格湍流扩散（K-theory），最外 1 圈为边界，不做更新。 !! 假设密度的扰动相对于平均密度很小，可以忽略。 !! 采用体积比计算扩散，确保污染物浓度分布与密度分布的一致性。 real ( fp ), intent ( in ) :: dt !! 积分时间: s real ( fp ), intent ( in ) :: dx (:) !! 网格长度: m real ( fp ), intent ( in ) :: kh (:) !! 扩散系数: m&#94;2/s real ( fp ), intent ( in ) :: rho (:) !! 密度: kg/m&#94;3 real ( fp ), intent ( inout ) :: c (:) !! 浓度: ug/m&#94;3 real ( fp ), intent ( out ) :: dc (:) !! 浓度变化: ug/m&#94;3 real ( fp ), optional , intent ( in ) :: volume (:) !! 网格体积: m&#94;3 ! local real ( fp ) :: u !! 虚拟速度，扩散等效为速度: m/s real ( fp ) :: fR , fL !! 某个网格右侧和左侧通量: ug/m&#94;2/s real ( fp ) :: flux ( size ( c )) !! 所有网格的右侧通量: ug/m&#94;2/s integer :: i dc = 0.0_fp ! 计算通量 do i = 1 , size ( c ) - 1 ! flux(i) > 为流出入 u = kh ( i ) * 2._fp / ( dx ( i ) + dx ( i + 1 )) flux ( i ) = ( rho ( i ) + rho ( i + 1 )) / 2.0 * u * ( c ( i + 1 ) / rho ( i + 1 ) - c ( i ) / rho ( i )) end do ! 时间积分, 前向欧拉 do i = 2 , size ( c ) - 1 fR = flux ( i ) fL = flux ( i - 1 ) if ( present ( volume )) then ! 进行校正体积校正，保证质量守恒 if ( fL < 0 ) fL = fL * volume ( i - 1 ) / volume ( i ) if ( fR > 0 ) fR = fR * volume ( i + 1 ) / volume ( i ) end if dc ( i ) = dt / dx ( i ) * ( fR - fL ) c ( i ) = c ( i ) + dc ( i ) if ( c ( i ) < eps ) c ( i ) = 0. ! 数值稳定性 end do end subroutine hdiff1d_by_k_theory end module mod_hdiff","tags":"","loc":"sourcefile/module_hdiff.f90.html"},{"title":"手册 – diffusion","text":"","tags":"","loc":"page/index.html"},{"title":"开发者手册 – diffusion","text":"混沌系统并不是随机系统。混沌研究的是秩序与随机性之间的过渡 --- 詹姆斯·格雷克 1.1 引言 大气的运动会引起污染物浓度在时空上的变化，在空气质量数值模式的语境中，\n污染物的湍流扩散是指由风引起，大尺度风场无法直接模拟的那部分污染物浓度变化过程。 风矢量是时空上的连续函数，但数值模式对大气方程进行了离散化，使得风矢量仅在网格点上计算。\n离散化后的风矢量（代表某一尺度的风）只能描述一个网格内、某段时间内（即一个时空网格）的平均状态或代表性风场。 从物理意义上讲，根据质量守恒原理(连续方程)，大气运动是连续发生的。而离散化的风矢量无法准确描述的、\n由次网格风引起的污染物浓度再分配，正是湍流扩散所对应的部分。 在风场的雷诺分解 (Reynolds Decomposition) 中， 代表扰动部分（次网格风），\n即小尺度、快速变化的风场成分。 通过对连续方程(一维)进行雷诺分解可得： 其中， 公式 1.4 右边的第二项 即为扩散项（湍流扩散项）的数学表达。\n该项描述了湍流对浓度分配的贡献，使得离散化后的连续方程在数学上得以闭合。 需要注意的是，雷诺分解后的时间平均项代表的是平滑的背景场（时空平滑），依然是一个连续函数。\n而数值离散化则可以理解为对雷诺分解后的时间平均项进行离散采样，\n即在有限的网格点上近似表示这一连续背景场。 1.2 K理论 由于 未知，因此需要采用参数化方案来估算扩散项的影响。\n在气象学和大气污染扩散研究中，\n通常采用一阶涡流粘度方法（k理论、一阶闭合、梯度扩散假设）来求解湍流项，\nK理论假设湍流扩散与分子扩散类似，引入了一个湍流扩散系数 来代替分子扩散系数 . 为什么能把湍流扩散类比为分子扩散？尽管它们的物理机制和作用范围不同（流体的馄钝性、分子的不规则运动），\n这两种扩散机制在物质传输过程中具有一定的相似性。 高浓度向低浓度扩散：大气中的湍流扩散和分子扩散的基本原理都是物质在不同浓度梯度下的流动，\n最终使得浓度达到均衡。 随机性：分子扩散是由分子的热运动引起的随机过程，\n而湍流扩散是由气流中的湍流湍动引起的随机过程（在雷诺分解表现为均值为0的白噪音）。 散的速率直观上与物质浓度差以及介质和环境条件（如气压、温度等）密切相关。根据菲克定律（可以看看 热传导方程 ），\n分子扩散的通量公式为： 为了描述湍流导致的浓度演变，引入湍流扩散系数 （ ）， 湍流扩散的通量公式为： 需要注意，等式右侧的浓度梯度是通过雷诺分解后的平均项获得的，在后续的公式中，\n在后续公式推导中，为了简化表示， 记作 。 为了解决一致性问题（气象模式和空气质量模式的离散方案不同），用混合比 来表示浓度， 其中 ，\n假设密度的扰动可以忽略（小尺度空气密度的变化应该远远小于污染物浓度的变化），\n也就是 。那么 因此 因此湍流扩散项为: 推广到三维情况 对于地球大气而言，水平方向的湍流主要由风速剪切引起，而垂直方向的湍流则主要受温度梯度的影响。\n因此，数值模式通常分别处理污染物的水平扩散和垂直扩散。 1.3 水平扩散 相较于水平输送，水平扩散的影响通常较小，在某些情况下甚至可以忽略不计。 根据 Smagorinsky (1963) 公式，可求解水平扩散系数 其中：\n- 是 Smagorinsky 常数（通常取 0.1 - 0.25 ，大气应用中常用 0.2 ）。\n- 是网格尺度（通常定义为网格间距的几何平均： ）。\n- 是水平速度的二阶张量的模量，定义为： 可以增加一个常数，让模拟的浓度变得更加光滑（Anthes and Warner, 1978）。 其中 的值为 为了便于通量计算，水平扩散系数 通常是在 Arakawa C 网格上定义，对于一个网格单元，其四条边上的每一条边均需计算一个 值，以准确描述水平扩散过程。 的计算涉及到速度的空间差分，采用中值差分方案求速度梯度。 在求解 v-stag 网格点上的 在 u-stag 网格上的梯度时，需要首先将 从 v-stag 网格插值到 mass 网格。为了简单，通常采用线性插值方法。 因此 对于扩散方程，采用前向欧拉进行时间差分。 而在计算浓度梯度时，采用中值差分方案计算浓度梯度（注意计算的是 u-stag 或者 v-stag 网格上的浓度梯度）。 因此，最后的积分表达式为 注意扩散是从高浓度到低浓度， 始终大于0，浓度梯度决定通量的方向。 1.4 垂直扩散 白天行星边界层中湍流作用较强，垂直扩散在污染物的垂直浓度分布中起着至关重要的作用。 垂直方向的湍流，受温度梯度的影响和风切边的功能影响，因此垂直扩散系数的求解更为复杂。\n该模块中不实现垂直扩散系数的计算。 在已知垂直扩散系数，用数值方法求解垂直扩散方程的难点在于垂直层网格的大小变化较为明显(从 10 ~ 1000m )。 采用显示的时间积分方案会要求 dt 非常小，会显著增加计算时间，因此通常采用隐式的时间积分方案。 与水平扩散一样，在计算浓度梯度时，采用中值差分方案计算浓度梯度（注意计算的是 w-stag 网格上的浓度梯度，做密度平均的时候，忽略了两层厚度之间的差异，同时任务两个格点中心之间的梯度处处相等）。 而时间差分方案选用后向欧拉。 带入浓度梯度公式 假设 则 将后向欧拉方程改写为矩阵形式 ，其中 为 的三对角元素表达式如下： 下对角线元素（第 行，第 列） ： 主对角线元素（第 行，第 列） ： 上对角线元素（第 行，第 列） ： 1.5 参考文献 Smagorinsky J. General circulation experiments with the primitive equations: I. The basic experiment[J]. Monthly weather review, 1963, 91(3): 99-164. Anthes R A, Warner T T. Development of hydrodynamic models suitable for air pollution and other mesometerological studies[J]. Monthly Weather Review, 1978, 106(8): 1045-1078.","tags":"","loc":"page/dev.html"},{"title":"使用者手册 – diffusion","text":"Warning 万丈高楼平地起 use diffusion , only : cal_hdiff_k , vdiff_by_k_theory , hdiff_by_k_theory call vdiff_by_k_theory ( dt , kz , dz , rho , conc )","tags":"","loc":"page/user.html"}]}