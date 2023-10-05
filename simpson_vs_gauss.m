function simpson_vs_gauss
%{ 
    simpson_vs_gauss.m 
    ~~~~~~~~~~~~~~~~~~~
    比較辛普森積分和高斯積分的誤差。
    範例函數為 x^2*exp(x)。
    比較的範圍為 n= 3 到 10。
    高斯積分程式碼取自 https://www.mathworks.com/matlabcentral/fileexchange/4540-legendre-gauss-quadrature-weights-and-nodes
    
    Chang Kai-Po @ Jian Lab, NCTU, Taiwan, 2023/03/20

%} 

    % 目標函數
    f = @(x) x.^2 .* exp(x);

    % 解析解
    exact = @(a,b) (b.^2 - 2*b + 2) .* exp(b) - (a.^2 - 2*a + 2) .* exp(a);

    % 高斯積分
    function integral = gauss(f, a, b, n)
        [x, w] = leggauss(n, -1, 1);
        integral = (b - a) / 2 * sum(w .* f((b - a) / 2 * x + (b + a) / 2));
    end

    % 辛普森積分
    function integral = simpson(f, a, b, n)
        h = (b - a) / n;
        x = linspace(a, b, n + 1);
        y = f(x);
        integral = h / 3 * (sum(y(1:2:end-1)) + 4 * sum(y(2:2:end)) + sum(y(3:2:end)));
    end

    % 計算在 n = 5 時的精確積分值，以及辛普森積分和高斯積分的值
    exact_integral = exact(0, 10);
    simpson_integral = simpson(f, 0, 10, 10);
    gauss_integral = gauss(f, 0, 10, 5);

    % 列印結果
    simpson_error = abs(exact_integral - simpson_integral) / exact_integral;
    gauss_error = abs(exact_integral - gauss_integral) / exact_integral;

    % Print the results
    fprintf('Value for exact integral of x^2*exp(x) from 0 to 10 at n = 5: %g\n', exact_integral);
    fprintf('Value for Simpson integral of x^2*exp(x) from 0 to 10 at n = 5: %g\n', simpson_integral);
    fprintf('Value for Gauss integral of x^2*exp(x) from 0 to 10 at n = 5: %g\n', gauss_integral);
    fprintf("Simpson's rule error: %g\n", simpson_error);
    fprintf('Gaussian quadrature error: %g\n', gauss_error);

    % 計算n=3~10時的誤差
    simpson_errors = zeros(1, 7);
    gauss_errors = zeros(1, 7);
    for i = 3:9
        simpson_errors(i-2) = abs(exact_integral - simpson(f, 0, 10, 2*(i+1))) / exact_integral;
        gauss_errors(i-2) = abs(exact_integral - gauss(f, 0, 10, i+1)) / exact_integral;
    end

    % 對誤差作圖
    x = 3:9;
    plot(x, simpson_errors, 'o-', 'DisplayName', 'Simpson rule');
    hold on;
    plot(x, gauss_errors, 'o-', 'DisplayName', 'Gaussian quadrature');    
    xlabel('n');
    ylabel('Error');
    title('Errors of Simpson rule and Gaussian quadrature');
    legend('Location', 'northeast');
    
end 

function [x,w]=leggauss(N,a,b)

% lgwt.m
%
% This script is for computing definite integrals using Legendre-Gauss 
% Quadrature. Computes the Legendre-Gauss nodes and weights  on an interval
% [a,b] with truncation order N
%
% Suppose you have a continuous function f(x) which is defined on [a,b]
% which you can evaluate at any x in [a,b]. Simply evaluate it at all of
% the values contained in the x vector to obtain a vector f. Then compute
% the definite integral using sum(f.*w);
%
% Written by Greg von Winckel - 02/25/2004
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';

% Initial guess
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);

% Legendre-Gauss Vandermonde Matrix
L=zeros(N1,N2);

% Derivative of LGVM
Lp=zeros(N1,N2);

% Compute the zeros of the N+1 Legendre Polynomial
% using the recursion relation and the Newton-Raphson method

y0=2;

% Iterate until new points are uniformly within epsilon of old points
while max(abs(y-y0))>eps
    
    
    L(:,1)=1;
    Lp(:,1)=0;
    
    L(:,2)=y;
    Lp(:,2)=1;
    
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k;
    end
 
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    
    y0=y;
    y=y0-L(:,N2)./Lp;
    
end

% Linear map from[-1,1] to [a,b]
x=(a*(1-y)+b*(1+y))/2;      

% Compute the weights
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end 
