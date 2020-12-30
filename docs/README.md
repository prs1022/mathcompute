# 数值分析 (Numerical Analysis)

## 第一章  数值分析与科学计算引论

计算机解决科学问题经历的几个过程

实际问题 --> 数学模型 --> 数值计算方法 --> 程序设计 --> 上机计算结果

数值分析也称为计算方法，是计算数学的一个主要部
分.

    计算数学是数学科学的一个分支，主要研究用计算机求解各种数学问题的数值计算方法及其理论与软件实现. 

数值分析的内容包括函数的数值逼近、数值微分与数
值积分、非线性方程数值解、数值线性代数、常微和偏微
数值解等.


## 第二章  插值法

- 拉格朗日插值公式及余项

定理: 已知 x0, x1 ,....,xn是n+1个互异的节点，则满足插值条件pn(xj)=yj=f(xj),j=0,1,...,n的n次多项式
![](http://latex.codecogs.com/gif.latex?p_n(x)=a_0+a_1x+a_2x^2+...+a_nx^n) 是**存在且唯一**的

- 牛顿插值
均差、插值公式及余项

- 分段插值
分段线形、分段三次埃尔米特、三次样条

## 第三章  函数逼近与曲线拟合

- 曲线拟合的最小二乘法
原理、性质、计算过程

- 正交多项式的构造
勒让德、埃尔米特、切比雪夫

- 最优一致逼近
用 切比雪夫多项式和泰勒级数展开作最优一致逼近的近似计算

## 第四章  数值积分与数值微分
- 等距节点的求积分公式、
- 牛顿型数值积分公式、
- 梯形公式、
- 辛普森公式、
- 复化梯形公式、
- 复化辛普森公式
- 龙贝格积分法、
- 不等距节点的高斯积分公式、
- 误差计算

## 第五章 迭代法
- 简单迭代法
- 牛顿迭代法及其变形
- 线性方程组的迭代法和误差估计
- 迭代松弛迭代法
- 非线性方程组的迭代法

## 第六章 矩阵的特征值和特征向量
- 乘幂法和反幂法计算方阵的特征值
- 雅克比
- 矩阵的householder变换
- QR分解

## 第七章 常微分方程初值问题的数值解法
- 欧拉法及改进欧拉法、向前、向后、中点
- 误差估计、绝对稳定性、隐式计算格式的使用、
- 泰勒级数法、龙格-库塔法 及其绝对稳定性分析
- 多步法（哈明公式、Adams公式）
- 外推法
- 微分方程组的常用计算格式

## 课后习题答案下载
[文档：PPT_数值分析_李庆扬_第五版.rar](http://note.youdao.com/noteshare?id=4766b25f16aa9530f269aae6cc628162)

[文档：数值分析课程第五版课后习题答案(李庆扬)](http://note.youdao.com/noteshare?id=1df76fadfd275d7b93447f5ce0932a52)
 

# Matlab

## 插值与拟合

### 拉格朗日插值
```
function y0 = Lagrange(X,Y, x0 )
%X,Y是已知的插值点坐标点
%x0是插值点
%y0是Lagrange插值点处的函数值
m=length(X);
%N为权系数
N=zeros(m,1);
y0=0;
for i=1:m
    N(i)  =1;
    for j=1:m
        if j~=i
            N(i)=N(i)*(x0-X(j))/(X(i)-X(j));
        end
    end
    y0=y0+Y(i)*N(i);
end

测试
clear;clc
%拉格朗日插值
x=0:0.1:0.4;
y=[0.5 0.5398 0.5793 0.6179 0.7554];
%求解f(0.13)
x1=0.13;
y1=Lagrange(x,y,x1)
disp('精确解->')
fprintf('%f\n',normcdf(x1,0,1))
disp('与精确解的误差->')
fprintf('%f\n',abs(y1-normcdf(x1,0,1)))

```
### 牛顿插值
```
function s = niudun(x,y,yt)
%NIUDUN Summary of this function goes here
%   Detailed explanation goes here
%定义符号变量p
syms p;
s=y(1);
n=length(x);
xishu=0;
dxs=1;
%求出牛顿插值系数
for i=1:n-1
    for j=i+1:n
        xishu(j)=(y(j)-y(i))/(x(j)-x(i));
    end
    temp1(i)=xishu(i+1);
    dxs=dxs*(p-x(i));
    s=s+temp1(i)*dxs;
    y=xishu;
end
simplify(s);
%若插值点函数缺省，则输出多项式
if(nargin==2)
    s=subs(s,'p','x');
    s=collect(s);
    s=vpa(s,4);
else
    %读取要插值点向量长度
    m=length(yt);
    for i=1:m
        temp2(i)=subs(s,'p',yt(i));
    end
    %得到系列插值点的插值结果
s=temp2;
end
```
### 最小二乘法拟合
```
clear;
clc;
x = -1:1:3;
y = [0 2 2.8 3.6 4.8];
x=x';
y=y';
%c为要求的变量
%形成线性方程的系数矩阵
a=[ fun(x(1))
    fun(x(2))
    fun(x(3))
    fun(x(4))
    fun(x(5))
     ]
b=y;
%法方程
A=a'*a;
B=a'*b;
c=A\B;
%c就是拟合基函数系数
%作出原离散数据，与拟合函数图像
x_n=-3:1:3;
y_n=c(1)*1+c(2)*x_n+c(3)*x_n.^2;
fprintf('拟合二次函数%d*x^2+%d*x+%d\n',c(3),c(2),c(1))
figure(2);
plot(x,y,'*',x_n,y_n)
grid
function f=fun(x)
    f(1) = 1;
    f(2) = x;
    f(3) = x^2;
end

```

## 积分与微分

### 复化梯形

### 梯形公式外推
```
function S=inte_tra_rec(f,a,b,t)
%梯形公式的外推法
%f为被积函数，a为积分上限，b为积分下限，t为外推次数
h=b-a;
S=h*(f(a)+f(b))/2;
if t>0
for i=1:t
   h=h/2;
    S=S/2+h*sum(f(linspace(a+h,b-h,2^(i-1))));
end
end
```

### 辛普森
```
function s = xps(fun, a,b,n )
%XPS Summary of this function goes here
%   Detailed explanation goes here
%a,b为积分区间
%n为划分的子区间个数
%如果n缺省，则默认划分100个区间
if nargin<3
    n=100;
end
%子区间长度
h=(b-a)/(2*n);
%设初值
s1=0;
s2=0;
for k=1:n
    x=a+h*(2*k-1);
    s1=s1+fun(x);
end
for k=1:(n-1)
    x=a+h*2*k;
    s2=s2+fun(x);
end
%结果想加
eps=0.00001;
s=h*(fun(a+eps)+fun(b)+4*s1+2*s2)/3;
end

```
### 复化辛普森
```
function S = FSimpson( f,a,b,N )
%FSIMPSON Summary of this function goes here
%   Detailed explanation goes here
% a,b表示区间端点
%N表示区间个数
%S是复化Simpson公式求得的积分值
h=(b-a)/N;
fa=feval(f,a);
fb=feval(f,b);
S=fa+fb;
x=a;
for i=1:N
    x=x+h/2;
    fx=feval(f,x);
    S=S+4*fx;
    x=x+h/2;
    fx=feval(f,x);
    S=S+2*fx;
end
S=h*S/6;
end
```

### 龙贝格
```
function s= rombg( a,b,TOL)
%ROMBG Summary of this function goes here
%   Detailed explanation goes here
n=1;
h=b-a;
%设置设计误差初值
delt=1;
x=a;
k=0;
R=zeros(4,4);
eps=0.00001;
R(1,1)=h*(fun(a+eps)+fun(b))/2;
while delt>TOL
    %如果两次计算的差值大于给定误差则进入循环
    k=k+1;
    h=h/2;
    s=0;
    for j=1:n
        x=a+h*(2*j-1);
        s=s+fun(x);
    end
    R(k+1,1)=R(k,1)/2+h*s;
    n=n*2;
    for i=1:k
        R(k+1,i+1)=((4^i)*R(k+1,i)-R(k,i))/(4^i-1);
    end
    %前后两次值的差
    delt=abs(R(k+1,k)-R(k+1,k+1));
end
s=R(k+1,k+1);
end
```

## 解线性方程组

1111


## 迭代

### 雅克比
```
function y = yacobi(A,b,k)
%--雅可比迭代法
% A为原始矩阵，b为列向量，k为迭代次数
n=size(A,2);
D=diag(diag(A));
L=-tril(A)+D;
U=-triu(A)+D;
M=D;
I=eye(n);
B=I-inv(D)*A;
f=inv(D)*b;
x=zeros(n,1);
%迭代过程
for kk=1:k
   for i=1:n
    sum=0;
    for j=1:n
        if j~=i
          sum=sum+A(i,j)*x(j);
        end
    end
    x(i)=(b(i)-sum)/A(i,i);
   end
end
y=x;
end

```

## 矩阵特征值和向量



## 常微分的初值问题

1111

## 扩展内容

1111

# 文档编写
[熟悉Markdown语法](https://www.zybuluo.com/mdeditor) </br>
[熟悉LaTeX语法](https://www.cnblogs.com/ywl925/p/3671827.html) </br>
这里推荐一款markdown的编辑软件[Typora](https://www.typora.io/)，里面集成了LaTex，非常方便美观


# 下载地址
[文档：Numerical Mathematics and Computing ...](http://note.youdao.com/noteshare?id=7a6deffb034a22b0ae9f056c1ccabdaf)
