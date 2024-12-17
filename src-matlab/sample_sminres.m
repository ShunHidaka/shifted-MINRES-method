%
% sample program for shifted Krylov solver
% First update : 2024/12/12
% Last update  : 2024/12/17
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
%

% Prepare matrix $A, b, sigma$
% https://math.nist.gov/MatrixMarket/mmio/matlab/mmiomatlab.html
% http://www.damp.tottori-u.ac.jp/~hoshi/elses_matrix/ELSES_MATRIX_VCNT4000std_20130517.tgz
[A, rows, cols, entries] = mmread("ELSES_MATRIX_VCNT4000std_A.mtx");
N = rows;
% Prepare shits $sigma^{(m)}$
M = 3;
sigma = zeros(M, 1);
for m = 1:1:M
    sigma(m) = 0.001*m;
end
b = ones(N, 1);

max_itr = 100000;
threshold = 1e-13;

% Solve by shifted MINRES method
[x, flag, rres, itrs] = shifted_MINRES(A, b, N, sigma, M, max_itr, threshold);

% verification of results
true_res = zeros(M,1);
for m = 1:1:M
    r = b - (A*x(:,m) + sigma(m)*x(:,m));
    true_res(m) = norm(r)/norm(b);
end
