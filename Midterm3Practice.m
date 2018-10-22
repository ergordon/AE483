%% 3.1.1
clear;clc

M = [-1 2 -7 -6; -10 5 -5 -2; 5 -3 -9 -5; -10 3 0 -6];

N = mat2str((M+M')/2)

%% 3.1.2
clear;clc

M1 = [-1.0 10.0; 10.0 0.0];
M2 = [4.0 2.5; 2.5 4.0];
M3 = [6.0 8.5; 8.5 9.0];
M4 = [5.0 0.0; 0.0 4.0];
M5 = [9.0 5.5; 5.5 8.0];

eig(M1)
eig(M2)
eig(M3)
eig(M4)
eig(M5)

%% 3.1.3
clear;clc

F1 = [6.0 5.0 4.5; 5.0 5.0 2.5; 4.5 2.5 9.0];
g1 = [4.0; -4.5; -4.0];
h1 = [-3.0];
F2 = [8.0 7.0 8.5; 7.0 1.0 6.5; 8.5 6.5 1.0];
g2 = [-2.0; 1.0; 2.0];
h2 = [5.0];

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

Ad = [1.2 -4.9 0.3; 4.9 -3.8 4.2; 1.7 1.8 -0.2];
Bd = [-0.8 4.7; -3.6 -4.0; -2.5 2.2];
Cd = [4.1 4.3 0.0];

F = Bd'*Cd'*Cd*Bd+eye(2);
K = (Cd*Bd*(F\eye(size(F))))'
P = 1-K'*F*K

%% 3.1.5
clear;clc

Ad = [0.6 1.3 4.6; -2.3 2.5 -1.5; -4.1 -3.1 1.0];
Bd = [4.0 0.9; 2.6 2.8; -2.9 -0.4];
Rd = [6.0 1.5; 1.5 3.0];

F = (Bd'*Bd+Rd)';
K = Bd*(F\eye(size(F)));
mat2str(K')
P = mat2str(eye(size(K,1))-K*F*K')

%% 3.1.6
clear;clc

Ad = [-4.6 -4.8 4.0; -0.8 1.5 -2.2; -4.2 2.2 2.0];
Bd = [-1.8 -2.4; -1.6 -2.7; 0.0 0.1];

F = Bd'*Bd+eye(size(Bd'*Bd));
g = 2*(Bd'*Ad-Bd');

c = 2;
disp(mat2str(c*F\g))

K = c*F\g;
disp(mat2str((Ad'-K'*Bd'-eye(size(Ad)))*(Ad-Bd*K-eye(size(Ad)))+K'*K))

%% 3.1.8
clear;clc

Ad = [-1.1 4.5 0.3; 3.6 0.6 0.9; -2.6 -0.3 -2.5];
Bd = [3.8 1.1; -4.0 4.0; 0.4 3.3];
Qd = [4.0 2.5 2.5; 2.5 6.0 0.0; 2.5 0.0 8.0];
Rd = [4.0 3.5; 3.5 4.0];
P_iplusone = [6.0 0.5 2.5; 0.5 2.0 2.0; 2.5 2.0 4.0];

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

M = [10 9; 9 14];
x = [-4; 6];

mat2str(2*x'*M)

%% 3.2.3
clear;clc

C = [-2 0; -3 3];

disp(mat2str(C'*eye(2)*C))

%% 3.2.4
clear;clc

a = 12;
b = 22;
c = 16;

m11 = a;
m12 = b/2;
m21 = m12;
m22 = c;

disp(mat2str([m11 m12; m21 m22]))

%% 3.2.5
clear;clc

Ad = [-3 9 5 -1; 3 -1 -6 9; 7 6 1 -3; -1 -4 -5 -6];
Bd = [-1 0; 7 0; 4 3; 4 -3];
Qf = [5 0 0 0; 0 9 0 0; 0 0 4 0; 0 0 0 8];
Qd = [5 0 0 0; 0 4 0 0; 0 0 2 0; 0 0 0 9];
Rd = [9 0; 0 4];
z = [6; -8; 0; 6];

h = z'*(Ad'*Qf*Ad+Qd)*z;
g = (2*Bd'*Qf*Ad)*z;
F = Bd'*Qf*Bd+Rd;

1/2*z'*Qf*z

%% 3.2.6
clear;clc

Ad = [7 -4 -5; 2 -7 -5; 9 4 4];
Bd = [9; 9; -4];
Qf = [4 0 0; 0 9 0; 0 0 8];
Qd = [5 0 0; 0 6 0; 0 0 3];
Rd = [3];
z = [8; -6; -2];

h = z'*(Ad'*Qf*Ad+Qd)*z
g = (2*Bd'*Qf*Ad)*z
F = Bd'*Qf*Bd+Rd

vd = (h-1/4*g'*(F\eye(size(F)))*g)/2;
ud = -1/2*(F\eye(size(F)))*g;

vpa(vd,8)
vpa(ud,8)

