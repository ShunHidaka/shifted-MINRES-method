%
% sample program for shifted Krylov solver
% First update : 2024/12/12
% Last update  : 2024/12/17
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
%

% Prepare matrix $A, b, sigma$
% https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html
% http://www.elses.jp/matrix/
[A, rows, cols, entries] = mmread("ELSES_MATRIX_CLIQ6912std_A.mtx");
%[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT900h_A.mtx");
N = rows;
% Prepare shits $sigma^{(m)}$
M = 3;
sigma = zeros(M, 1);
for m = 1:1:M
    sigma(m) = 0.001*m + 0.01i;
end
b = ones(N, 1);

% Solve by shifted MINRES method
max_itr = 100000;
threshold = 1e-13;
[x, flag, rres, itrs] = shifted_minres(A, b, N, sigma, M, max_itr, threshold);

% verification of results
true_res = zeros(M,1);
for m = 1:1:M
    r = b - (A*x(:,m) + sigma(m)*x(:,m));
    true_res(m) = norm(r)/norm(b);
end
