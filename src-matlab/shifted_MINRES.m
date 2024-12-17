%
% shifted MINRES method by MATLAB
% First update : 2024/12/12
% Last update  : 2024/12/13
% Created by "ShunHidaka (https://github.com/ShunHidaka)"
% 実数向けに作成、複素数用ではない
% 引数：
%   行列A, 右辺ベクトルrhs, 行列サイズN,
% 　シフトsigma, シフト数M,
%   最大反復回数max_itr, 閾値threshold
% 返り値：
%   近似解x(N行M列)
%   フラグflag
%   相対残差rres（M行）
%   収束したときの反復回数itrs（M行）
%

function [x, flag, rres, itrs] = shifted_MINRES(A, rhs, N, sigma, M, max_itr, threshold)
    x    = zeros(N, M);  % 近似解
    flag = 1;           % 正常に収束したかどうかを示すフラグ、0なら正常に実行
    rres = zeros(M, 1); % 各シフト方程式の相対残差
    itrs = zeros(M, 1); % 収束したときの反復回数
    % 使用する変数の宣言
    q       = zeros(N, 3);  % Lanczos基底 q_{j-1}, q_{j}, q_{j+1}
    alpha   = zeros(1, 1);  % Lanczos過程で得られる定数
    beta    = zeros(2, 1);  % Lanczos過程で得られる定数 beta_{j-1}, beta_{j}
    r       = zeros(3, M);  % 三重対角行列用のベクトル
    G1      = zeros(2,2,M); % Givens回転行列を格納する配列, G_{j-2}
    G2      = zeros(2,2,M); % Givens回転行列を格納する配列, G_{j-1}
    G3      = zeros(2,2,M); % Givens回転行列を格納する配列, G_{j}
    p       = zeros(N,3,M); % 補助ベクトル, p(:,1,m)=p_{j-2}^{(m)}, p(:,2,m)=p_{j-1}^{(m)}, p(:,3,m)=p_{j}^{(m)}
    f       = zeros(M);     % 解の更新に用いる一時変数
    h       = zeros(M);     % 残差のノルム
    is_conv = zeros(M, 1);  % 収束した方程式を管理するベクトル
    % 変数の初期化
    r0nrm = norm(rhs);
    q(:,2) = rhs / r0nrm;
    for m = 1:1:M
        f(m) = 1;
        h(m) = r0nrm;
        is_conv(m) = 0;
    end
    % メインループ
    conv_num = 0;
    for j = 1:1:max_itr
        % Lanczos過程
        q(:,3) = A*q(:,2) - beta(1)*q(:,1);
        alpha(1) = dot(q(:,3), q(:,2));
        q(:,3) = q(:,3) - alpha(1)*q(:,2);
        beta(2) = norm(q(:,3));
        % 各シフト方程式の更新
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
            % シフト方程式の収束判定
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
        % 全ての方程式が収束したかの判定
        if conv_num == M
            flag = 0;
            break;
        end
    end
