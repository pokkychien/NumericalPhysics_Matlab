% potential_well.m
% ~~~~~~~~~~~~~~~~
% 使用有限差分法估計一個電子在一階位能井中的波函數機率分布，
% 以及其對應的能量。
% 位能井的長度為 1，電子的質量為electron_mass，位能井的深度為無限。

% Created by Chang Kai-Po @ Jian Lab, 2023/03/18

clc;clear;close;
electron_mass = 9.10938356e-31;
hbar = 1.0545718e-34;

[eigenvalues, eigenvectors] = coeff_matrix(500);
epsilon = 0.002;
disp("First ten eigenvalues:")
disp(arrayfun(@(e) sprintf('%.7f', real(e)), eigenvalues(1:10), 'UniformOutput', false))

figure
for i = 1:3
    % Energy = eigenvalue *0.5* hbar^2 / epsilon^2 / m
    energy = eigenvalues(i) * 0.5 * hbar^2 / epsilon^2 / electron_mass;
    subplot(3, 1, i)
    plot(eigenvectors(:, i), 'DisplayName', sprintf('n=%d E=%f', i, real(energy)))
    ylabel("$\psi(x)$", 'Interpreter', 'latex')
    legend('Location', 'northeast')
end 
xlabel('$x$', 'Interpreter', 'latex')


% 有限差分法的係數矩陣。
% 方程式為-(hbar^2 / 2m) d^2(psi(x) / dx^2) = E psi(x), 0 < x < L
% 邊界條件為psi(0) = psi(L) = 0
% 運算後，有限差分法的係數矩陣為
% A = [2 -1 0 0 0 ... 0]
%     [-1 2 -1 0 0 ... 0]
%     [0 -1 2 -1 0 ... 0]
%     [.................]       
%     [0 0 ... -1 2 -1]
%     [0 0 ... ... -1 2]
% 因此為一特徵值問題，將其轉換為求解A的特徵值和特徵向量。
% 其中特徵值為E*2*epsilon^2*m/hbar^2，特徵向量為psi(x)。

function [sorted_eigenvalues, sorted_eigenvectors] = coeff_matrix(N)
    A = zeros(N, N);
    A(1, 1) = 2;
    A(N, N) = 2;
    A(1, 2) = -1;
    A(N, N-1) = -1;
    for i = 2:N-1
        A(i, i-1) = -1;
        A(i, i) = 2;
        A(i, i+1) = -1;
    end
    A = sparse(A);
    [eigenvectors, eigenvalues] = eigs(A, N-2, 'sm');
    [sorted_eigenvalues, sorted_indices] = sort(diag(eigenvalues));
    sorted_eigenvectors = eigenvectors(:, sorted_indices);
end

