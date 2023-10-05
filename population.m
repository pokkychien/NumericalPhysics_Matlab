%{
    population.m
    ~~~~~~~~~~~~~
    dy/dt = ky - cy^2為具備外在競爭關係的的人口增長模型，
    本模組用來計算人口增長的數值解。

    This module contains a method for population growth estimation 
    modelled by an one-dimention differential equation, named logistic 
    equation of population growth.

    The differential equation is: dy/dt = ky - cy^2, and will be
    estimated by numerical methods. The initial condition is y(0) = 1.

    Chang Kai-Po @ Jian Lab 2023/03/13
%}
clc;clear;close;
xlow = 0; xhigh = 50; h = 0.1;
k = 0.1; c = 0.01;
level = round((xhigh-xlow)/h);
x = linspace (xlow, xhigh, level);
y1 = ranged_finit_diff(1, xlow, xhigh, k, c, h);
y2 = ranged_exact_sol(1, xlow, xhigh, k, c, h);
plot(x, y1, x, y2);

function output = finit_diff (y, k, c, h)
    %{
    根據定義計算此微分方程的有限差分。
    若dy/dt = ky - cy^2，經過推導可得到
    f(x+h) = (1+hk)f(x) - hcf(x)^2
    %}
    output = (1+h*k)*y - h*c*y^2;
end 

function output = ranged_finit_diff (y0, xlow, xhigh, k, c, h)
    %{
    透過有限差分計算此微分方程在在xlow與xhigh此一範圍之間的估計值。
    %}
    output = [y0];
    level = round((xhigh-xlow)/h);
    x = linspace (xlow, xhigh, level);
    for i = 1:length(x)-1
        output = [output finit_diff(output(i), k, c, h)];
    end     
    %disp(output);
    %plot(x, output);
end 

function output = ranged_exact_sol(y0, xlow, xhigh, k, c, h)
    %{
    傳回此微分方程精確解在xlow與xhigh此一範圍之間的結果。
    詳解請見 https://ch-hsieh.blogspot.com/2016/03/blog-post_10.html
    %}    
    level = round((xhigh-xlow)/h);
    x = linspace (xlow, xhigh, level);
    a = k*y0;
    b = k-(c*y0);
    d = c*y0;
    output = a ./ (b*exp(-k*x) + d);
    %disp(output);
    %plot(x, output);
end 
