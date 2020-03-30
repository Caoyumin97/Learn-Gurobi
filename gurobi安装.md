Gurobi+python安装

1. 获取license

2. 下载软件包

3. 运行快捷方式 --> 验证license --> 保存license文件（gurobi.lic）

4. 系统环境变量添加GRB_LICENSE_FILE --> license文件位置

5. /gurobi/win64/python36/lib/gurobipy --> copy to python36/lib/site-packages/

6. 运行/gurobi/win64/bin/pysetup.bat --> 输入python36文件夹位置

7. 安装完成，运行试验代码

   ```python
   from gurobipy import *
   
   model = Model("gurobi_test")
   ```

   