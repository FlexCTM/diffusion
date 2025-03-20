## diffusion

次网格湍流模块

## 安装
1. 使用MAKE安装
```
make lib
```
生成的库的路径为 ```build/libdiffusion.a```

2. 使用CMAKE安装
```
mkdir build; cd build

FC=ifort cmake ../ # 正常编译
FC=ifort cmake  -DCMAKE_BUILD_TYPE=Debug ../ # 调试

make
```

3. 使用包管理器 fpm 安装
```
fpm build --compiler 'ifort'
```

如果采用 fpm 管理依赖关系和编译过程，则只需要在项目配置文件中引入仓即可。
```
[dependencies]
diffusion.git = "git@github.com:FlexCTM/diffusion.git"
```

## 单元测试
```
fpm test --compiler 'ifort'
```

## 调试
```
COMPIFLE=gnu make gdb
```

### 案例代码

``` fortran
    use mod_tool, only: fp
    use diffusion, only: cal_hdiff_k, vdiff_by_k_theory, hdiff_by_k_theory

    implicit none

    integer, parameter :: nz = 5

    integer :: i
    real(fp) :: dt  !! 积分时间: s
    real(fp) :: kz(nz) !! 扩散系数: m2/s
    real(fp) :: dz(nz) !! 垂直网格宽度: m

    real(fp) :: rho(nz) !! 密度
    real(fp) :: conc(nz) !! 浓度

    dt = 10.
    ! 读取数据
    kz = 0.1
    dz = 1.0
    rho = 1.0
    conc = [10, 8, 5, 10, 3]

    ! 垂直扩散
    do i = 1, 20
        call vdiff_by_k_theory(dt, kz, dz, rho, conc)
    end do
```

## License
This software is distributed under a **dual-license model**:  

- **Open Source License (GPLv3)**:  
  - Academic institutions, researchers, and open-source developers may use, modify, and distribute the code under the terms of the **GNU General Public License v3 (GPLv3)**.  
  - Any modifications or derivative works **must also be open-sourced** under GPLv3.  

- **Commercial License**:  
  - Companies or organizations that wish to use **this software** in proprietary, closed-source, or commercial applications must obtain a separate **commercial license** from the author.  
  - For commercial licensing inquiries, please contact: **[xiaolh09@lzu.edu.cn]**  

© 2025 [xiaolh]. All rights reserved.
