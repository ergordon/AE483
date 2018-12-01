%% 3.1.1
    clear;clc
M = [4 -5 7 10; 3 -7 1 -6; 4 -9 2 2; -6 6 1 5];
N = mat2str((M+M')/2)
    %% 3.1.2
    clear;clc
M1 = [-1.0 7.0 6.5; 7.0 2.0 3.5; 6.5 3.5 3.0];
M2 = [-1.0 8.0 6.5; 8.0 4.0 5.0; 6.5 5.0 -1.0];
M3 = [7.0 6.5 4.0; 6.5 9.0 5.5; 4.0 5.5 6.0];
M4 = [7.0 1.5 1.0; 1.5 4.0 7.0; 1.0 7.0 1.0];
M5 = [2.0 7.0 4.0; 7.0 3.0 2.0; 4.0 2.0 5.0];
    eig(M1)
    eig(M2)
    eig(M3)
    eig(M4)
    eig(M5)
    %% 3.1.3
    clear;clc
F1 = [3.0 -0.5; -0.5 3.0];
g1 = [1.5; -3.5];
h1 = [1.0];
F2 = [4.0 9.5; 9.5 0.0];
g2 = [-4.0; 3.5];
h2 = [-0.5];
    Z1 = -1*(F1\eye(size(F1)))*g1;
    Z2 = -1*(F2\eye(size(F2)))*g2;
    x = Z1; F = F1; g = g1; h = h1;     % test minimizer
    % x = Z2; F = F2; g = g2; h = h2;   % comment out the other one and remove
                                        % the comment on this one to see which
                                        % value is smaller, thats the answer
    disp(mat2str(x))
    disp(mat2str((1/2)*x'*F*x+g'*x+h))      
    %% 3.1.4
    clear;clc
Ad = [4.4 -1.0 0.8; 2.8 -1.1 -2.1; -0.6 -1.0 -1.7];
Bd = [-1.4 -4.8; -1.5 3.3; 3.8 -0.3];
Cd = [3.8 3.1 -2.6];
    F = Bd'*Cd'*Cd*Bd+eye(2);
    K = (Cd*Bd*(F\eye(size(F))))';
    P = 1-K'*F*K;
    
    mat2str(K)
    mat2str(P)
    %% 3.1.5
    clear;clc
Ad = [-4.6 4.8 -3.7; -3.9 4.0 -4.3; -3.6 -1.1 -2.4];
Bd = [-3.3 -1.2; 1.0 -4.9; 2.5 3.2];
Rd = [8.0 5.0; 5.0 8.0];
    F = (Bd'*Bd+Rd)';
    K = Bd*(F\eye(size(F)));
    mat2str(K')
    P = mat2str(eye(size(K,1))-K*F*K')
    %% 3.1.6
    clear;clc
Ad = [-4.1 0.1 -4.5; 1.2 -1.4 1.0; -0.9 -4.9 4.9];
Bd = [1.6 -0.6; -0.9 4.9; -2.2 4.0];
    F = Bd'*Bd+eye(size(Bd'*Bd));
    g = 2*(Bd'*Ad-Bd');
    c = 2;
    disp(mat2str(c*F\g))
    K = c*F\g;
    disp(mat2str((Ad'-K'*Bd'-eye(size(Ad)))*(Ad-Bd*K-eye(size(Ad)))+K'*K))
    %% 3.1.8
    clear;clc
Ad = [0.0 4.9 0.3; -1.4 4.8 1.5; -2.7 -4.2 -3.5];
Bd = [-1.5 -3.9; -0.2 -1.3; 0.5 -4.0];
Qd = [10.0 4.0 4.0; 4.0 4.0 3.0; 4.0 3.0 3.0];
Rd = [10.0 4.5; 4.5 8.0];
P_iplusone = [10.0 3.0 4.5; 3.0 1.0 2.0; 4.5 2.0 10.0];
    F = Bd'*P_iplusone*Bd+Rd;
    g = 2*Bd'*P_iplusone*Ad;
    h = Ad'*P_iplusone*Ad+Qd;
    mat2str(2*F\g)
    K = 2*F\g;
    % K'*F*K
    % mat2str(h-1/4*g'*F*g)
    P = mat2str(Ad'*P_iplusone*(Ad-Bd*K)-K'*Bd'*P_iplusone*(Ad-Bd*K)+Qd+K'*Rd*K)
    %% 3.2.1
    clear;clc
    % c'
    %% 3.2.2
    clear;clc
    M = [2 3; 3 6];
    x = [-5;-2];
    mat2str(2*x'*M)
    %% 3.2.3
    clear;clc
C = [0 -1; -3 3];
    disp(mat2str(C'*eye(2)*C))
    %% 3.2.4
    clear;clc
    a = 18;
    b = 6;
    c = 16;
    m11 = a;
    m12 = b/2;
    m21 = m12;
    m22 = c;
    disp(mat2str([m11 m12; m21 m22]))
    %% 3.2.5
    clear;clc
Ad = [1 -6 -4 8; 4 -1 -7 -5; 5 6 -2 -3; 3 9 -7 -3];
Bd = [5 -9; -2 -8; -1 3; 0 -3];
Qf = [8 0 0 0; 0 7 0 0; 0 0 8 0; 0 0 0 6];
Qd = [4 0 0 0; 0 2 0 0; 0 0 4 0; 0 0 0 5];
Rd = [6 0; 0 8];
z = [-5; 9; 5; 5];
    h = z'*(Ad'*Qf*Ad+Qd)*z;
    g = (2*Bd'*Qf*Ad)*z;
    F = Bd'*Qf*Bd+Rd;
    1/2*z'*Qf*z
    %% 3.2.6
    clear;clc
Ad = [8 -1 5 -7; -6 -8 -3 -7; -1 -6 -9 -6; 7 8 3 0];
Bd = [6 -3; -3 3; -2 0; -5 0];
Qf = [7 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 6];
Qd = [8 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1];
Rd = [9 0; 0 5];
z = [-9; -8; 5; 3];
    h = z'*(Ad'*Qf*Ad+Qd)*z
    g = (2*Bd'*Qf*Ad)*z
    F = Bd'*Qf*Bd+Rd
    vd = (h-1/4*g'*(F\eye(size(F)))*g)/2;
    ud = -1/2*(F\eye(size(F)))*g;
    vpa(vd,8)
    vpa(ud,8)