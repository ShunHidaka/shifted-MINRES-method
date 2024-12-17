%
% shifted MINRES method by MATLAB
% First update : 2024/12/12
% Last update  : 2024/12/17
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
% 実数向けに作成、複素数用ではない
% Solve shifted linear systems
% (A + \sigma^{(m)}I_N)\vb{x}^{(m)} = \vb{b},  (m=1,\dots, M).
% argument：
%   "Matrix" A, "Right-Hand-Side-vector" rhs, "Matrix-size" N,
% 　"Shifts" sigma, "Number-of-shits" M,
%   "Maximum-iteration" max_itr, "relative-residual-norm-threshold" threshold
% Return：
%   "Approximate-solution" x (M \times N)
%   "Status" flag
%   "Relative-residual-norm" rres (M)
%   "Converged iterations" itrs (M)
%

function [x, flag, rres, itrs] = shifted_MINRES(A, rhs, N, sigma, M, max_itr, threshold)
    x    = zeros(N, M); % Approximate solutions
    flag = 1;           % if flag=0 then all Eqs converged
    rres = zeros(M, 1); % relative residual norm of each Eqs
    itrs = zeros(M, 1); % converged iterations
    % Declaer variables
    q       = zeros(N, 3);  % Lanczos basis, q_{j-1}, q_{j}, q_{j+1}
    alpha   = zeros(1, 1);  % Lanczos constant
    beta    = zeros(2, 1);  % Lanczos constant, beta_{j-1}, beta_{j}
    r       = zeros(3, M);  % Tridiagonal matrix
    G1      = zeros(2,2,M); % Givens rotation matrix, G_{j-2}
    G2      = zeros(2,2,M); % Givens rotation matrix, G_{j-1}
    G3      = zeros(2,2,M); % Givens rotation matrix, G_{j}
    p       = zeros(N,3,M); % auxiliary vec, p(:,1,m)=p_{j-2}^{(m)}, p(:,2,m)=p_{j-1}^{(m)}, p(:,3,m)=p_{j}^{(m)}
    f       = zeros(M);     % temporary variable for updating approximates
    h       = zeros(M);     % residual norm on Algorithm
    is_conv = zeros(M, 1);  % status of each Eqs
    % Initialize variables
    r0nrm = norm(rhs);
    q(:,2) = rhs / r0nrm;
    for m = 1:1:M
        itrs(m) = N;
        f(m) = 1;
        h(m) = r0nrm;
        is_conv(m) = 0;
    end
    % Main loop
    conv_num = 0;
    for j = 1:1:max_itr
        % Lanczos process
        q(:,3) = A*q(:,2) - beta(1)*q(:,1);
        alpha(1) = dot(q(:,3), q(:,2));
        q(:,3) = q(:,3) - alpha(1)*q(:,2);
        beta(2) = norm(q(:,3));
        % update each shift Eq
        for m = 1:1:M
            if is_conv(m) == 1
                continue;
            end
            r(1,m)=0; r(2,m)=beta(1); r(3,m)=alpha(1)+sigma(m);
            if j >= 3
                r(1:2,m) = G1(:,:,m)*r(1:2,m);
            end
            if j >= 2
                r(2:3,m) = G2(:,:,m)*r(2:3,m);
            end
            G3(:,:,m) = planerot([r(3,m); beta(2)]);
            r(3,m) = G3(1,1,m)*r(3,m) + G3(1,2,m)*beta(2);
            p(:,1,m) = p(:,2,m);
            p(:,2,m) = p(:,3,m);
            p(:,3,m) = ( q(:,2) - r(1,m)*p(:,1,m) - r(2,m)*p(:,2,m) ) / r(3,m);
            x(:,m) = x(:,m) + r0nrm*G3(1,1,m)*f(m)*p(:,3,m);
            f(m) = G3(2,1,m)*f(m);
            h(m) = abs(G3(2,1,m))*h(m);
            G1(:,:,m) = G2(:,:,m);
            G2(:,:,m) = G3(:,:,m);
            % Determine convergence for each shift Eq
            rres(m) = h(m)/r0nrm;
            if rres(m) <= threshold
                is_conv(m) = 1;
                itrs(m) = j;
                conv_num = conv_num + 1;
            end
        end
        q(:,3) = q(:,3)/beta(2);
        q(:,1) = q(:,2);
        q(:,2) = q(:,3);
        beta(1) = beta(2);
        % Determine convergenc for all Eqs
        if conv_num == M
            flag = 0;
            break;
        end
    end
