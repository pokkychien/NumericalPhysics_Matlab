%{    
    simpson_int.m
    ~~~~~~~~~~~
    辛普森數值積分 (Simpson's rule) 的簡單範例。
    設一函數f(x)，上下界a與b，切割成n分，
    simpson(f,a,b,n) 將會對f(x)在a與b之間進行數值積分。
    
    此積分公式為：
    integral(f(x),a,b) = h/3 * (f(a)+f(b)+4*sum(f(a+2i*h))+2*sum(f(a+(2i+1)*h)))
    其中h=(b-a)/n，i=1,2,...,n-1。
    
    同時，exact_integral(a,b)為f(x)在a與b之間的解析解，將數值解與解析解相比較。
    例如說，f(x)=sin(x)在a與b之間積分的解析解為-cos(b)+cos(a)。

    Created by Chang Kai-Po @ Jian Lab, 2023/3/19
%}
clc;clear;close;
a = 0.0;
b = 5;
n = 10;
fprintf("將 sin(x) 從 %g 到 %g 做積分\n", a, b);
fprintf("解析解: %g\n", exact_integral(a,b));
fprintf("數值解: %g\n", simpson(@f,a,b,n));

function y = f(x)
    % 目標函數
    y = 1/(x);
end

function s = simpson(f,a,b,n)
    % 辛普森積分的計算
    h = (b-a)/n;
    s = f(a)+f(b); % 邊界條件
    for i = 1:2:(n-1) % 奇數項
        s = s + 4*f(a+i*h);
    end
    for i = 2:2:(n-2) % 偶數項
        s = s + 2*f(a+i*h);
    end
    s = s*h/3;
end

function y = exact_integral(a,b)
    % 解析解
    y = -cos(b) + cos(a);
end



