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
定理一: 已知 x0, x1 ,....,xn是n+1个互异的节点，则满足插值条件pn(xj)=yj=f(xj),j=0,1,...,n的n次多项式
![](http://latex.codecogs.com/gif.latex?p_n(x)=a_0+a_1x+a_2x^2+...+a_nx^n) 是**存在且唯一**的

## 第三章  函数逼近与曲线拟合
111111

## 第四章  数值积分与数值微分
111111

## 第五章

111111


## 课后习题答案下载
[文档：PPT_数值分析_李庆扬_第五版.rar](http://note.youdao.com/noteshare?id=4766b25f16aa9530f269aae6cc628162)

[文档：数值分析课程第五版课后习题答案(李庆..](http://note.youdao.com/noteshare?id=1df76fadfd275d7b93447f5ce0932a52)
 

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

### 梯形
### 辛普森
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

1111

## 矩阵特征值和向量

1111

## 常微分的初值问题

1111

## 扩展内容

1111

# 文档编写
[熟悉Markdown语法](https://www.zybuluo.com/mdeditor)

网站demo里支持Latex，而github不支持

所以还有另一种方案书写数学公式
通过 ‘http://latex.codecogs.com/gif.latex?\\frac{1}{1+sin(x)}’ 设置为图片地址

效果如下

![](http://latex.codecogs.com/gif.latex?\\frac{1}{1+sin(x)})


虽然是相当的不美观，但是看起来还是可以的


![](http://latex.codecogs.com/gif.latex?\\frac{x_1}{x_2+sin(x^2)})


# 下载地址
[文档：Numerical Mathematics and Computing ...](http://note.youdao.com/noteshare?id=7a6deffb034a22b0ae9f056c1ccabdaf)
