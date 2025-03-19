## diffusion
Sub-grid turbulent diffusion

## 安装
1. 使用MAKE安装
```
make lib
```

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

## 单元测试
```
fpm test --compiler 'ifort'
```

## 调试
```
COMPIFLE=gnu make gdb
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
