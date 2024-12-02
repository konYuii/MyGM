# 扫掠曲面实验
## 环境、编程语言及库
- Ubuntu 24.04
- C++
- 编译器：g++ 13.2.0
- GLFW等图形库由框架本身提供、添加了GSL库用于解微分方程
## 编译和运行
ubuntu系统下: 
进入project目录下执行以下脚本:
```sh
apt install libglfw3-dev # ubuntu
mkdir build
cd build
cmake ..
cmake --build .
./project
```
## 修改说明
### 对原有框架的修改
- src/GLProgram.cpp中将标架的绘制从单一颜色更改为了3个方向用3种颜色表示
- src/main.cpp中将sweptSurface、frame、profileCurve的生成改为调用相应的函数。将几种情况的输入改为预编译if语句，3种任务依次对应6个Case。
- 额外添加的功能代码均在wyx/目录文件下，其中myBasis.cpp为自定义的基函数库，用于计算基函数、曲线曲面点的坐标、导数等；sweptSurface.cpp为扫掠算法实现，作业中要求的3种扫掠算法均在该文件中，sweptSurface.h为头文件，包含了一些函数声明，并且包含形状控制矩阵的声明；frame.cpp中实现了3种标架的计算方法。
- 在CMakeLists.txt中添加wyx目录
### 数据结构
- 未增加新的数据结构
### 算法
#### 扫掠算法
**实现在sweptSurface.cpp中**
- ***SurfMeshParams***函数：计算曲面采样点对应的参数值
- ***GenerateKnotVector***函数：初始化knot vector
- ***GenerateSweptSurface0**函数：书上最简单情况的实现，作业中未要求，但实现了作为参考
- ***GenerateSweptSurface1**函数：作业要求的第一种实现，先用公式计算出采样点的三维坐标，再调用SurfMeshParams计算参数值后，用GenerateKnotVector初始化knot vector，最后用书上算法9.1实现插值计算NURBS曲面。
- ***GenerateSweptSurface2**函数：作业要求的第二种实现，利用论文中的公式直接计算出扫掠NURBS曲面的控制点，并沿用两条曲线的u、v等参数信息。
- ***GenerateSweptSurface3**函数：作业要求的第三种实现，用GenerateProfileCurve函数生成剖面线的一系列实例，然后用类似第一种实现中的插值方法插值处扫掠NURBS曲面的控制点。
- ***GenerateProfileCurve**函数：给出参数向量，生成剖面线。
#### 标架计算
**实现在frame.cpp中**
- ***getFrameMatrix***函数、***getInverseFrameMatrix***函数：将标架转为矩阵形式
- ***computeFrame1***函数：作业要求的第一种计算方法，计算Frenet标架
- ***computeFrame2***函数：作业要求的第二种计算方法，实现法线投影法
- ***computeFrame3***函数：作业要求的第三种计算方法，计算旋转最小量标架，参考Blom论文中给出的微分方程公式，调用GSL库求解。
- 3个***computeSingleFrame***函数：对应以上3种实现方法，分别计算单个标架。