%{
    single_slit.m
    ~~~~~~~~~~~~~~
    使用高斯積分進行單狹縫繞射計算。    
    Created by Chang Kai-Po @ Jian Lab, NCTU, Taiwan, on 2023/3/26.    
%}

clc;clear;close;
y = linspace(-0.1, 0.1, 1000);
z=zeros(1000, 1);
for i = 1:1000
    z(i) = slit(2, y(i), 630, 0.1);
end
plot(y, z, 'LineWidth', 1.5)
xlabel('y')
ylabel('Intensity')
title('Intensity of light passing through a single slit')
grid on

function slit_intensity = slit(L, y, wavelength, d)
% Using gaussian quadrature to calculate the intensity of light passing through a single slit.

% 在泰勒展開近似後，單狹縫公式在某一個點的值為：
% I(y, x) = exp(i*2*pi*a*(y-x)^2)
% 其中 a = k0 / 2L = 2*pi / (wavelength * 2L) = pi / (wavelength * L)

% 波長的單位為nm
a = pi / (wavelength * 10^-9 * L);
% 等等要對x做積分，所以先將y帶入公式
f = @(x) exp(1j*a*((y-x).^2));
% 積分範圍為[-d/2, d/2]
slit_intensity = abs(integral(f, -d/2*10^-3, d/2*10^-3));

end

