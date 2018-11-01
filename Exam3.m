function Exam3
clc;clear;
%% Control Instance 3

P3
end

function p1
Ad = [-2 8 -9; 6 -2 8; 9 1 8];
Bd = [1 -6; 7 -6; 7 4];
Qf = [9 0 0; 0 3 0; 0 0 2];
Qd = [6 0 0; 0 3 0; 0 0 6];
Rd = [5 0; 0 9];
Piplus = [2 4 5; 4 20 4; 5 4 20];
n=2

N=1
A=Ad; B=Bd; R=Rd; Q=Qd;

P = Qf;
F = zeros(size(Rd));
g = zeros(size(2*B'*Qf*A));
h = zeros(size(A'*P*A+Q));
K = zeros(size(g));

for i = 1:n-N
    F = B'*P*B+R;
    g = 2*B'*P*A;
    h = A'*P*A+Q;
    K = 2*F\g;
    P = h+g'/2*(-K)+(-K')*g/2+K'*F*K;
end

disp(mat2str(K))
disp(mat2str(P))
end

%% Working Code for Problem 1
function P1a

Ad = [-2 8 -9; 6 -2 8; 9 1 8];
Bd = [1 -6; 7 -6; 7 4];
Qf = [9 0 0; 0 3 0; 0 0 2];
Qd = [6 0 0; 0 3 0; 0 0 6];
Rd = [5 0; 0 9];
Piplus = [2 4 5; 4 20 4; 5 4 20];

F = Bd'*Piplus*Bd+Rd;
g = 2*Bd'*Piplus*Ad;
h = Ad'*Piplus*Ad+Qd;


K = 2*F\g;
P = Ad'*Piplus*(Ad-Bd*K)-K'*Bd'*Piplus*(Ad-Bd*K)+Qd+K'*Rd*K

mat2str(P)
mat2str(K)
end

function P2
Ad = [-3 -3; -5 0];
Bd = [-4; -7];
Qf = [8 0; 0 9];
Qd = [4 0; 0 6];
Rd = [7];

mat2str(Qf)
end

function P3
o_desired = [-0.22; 0.27; -0.63];

m = 0.47375000;
g = 9.81000000;
J = [0.00453000 0.00000000 0.00000000; 0.00000000 0.00450000 0.00000000; 0.00000000 0.00000000 0.00626000];
sigmamax = 1000.00000000;
kF = 0.00001000;
kM = 0.00000010;
L = 0.80000000;
dt = 0.05000000;

% Data in MATLAB format
Qd = diag([10, 1, 1, 10, 1, 100, 1, 100, 100, 1, 1, 1]);
Rd = diag([1, 10, 10, 10]);

x0 = [0.00; -0.10; -1.50; -0.60; 0.70; -0.40; -0.50; 0.10; 0.60; 0.30; -0.10; 0.70];
t0 = 0.00;
t1 = 0.35;

% Control policy
xe = zeros(12, 1);
xe(1:3)=o_desired;
ue = zeros(4, 1);
ue(4)=m*g;

mat2str(xe)
mat2str(ue)
    
    
end

function P4
Ad = [-2 -9 -8; 5 -5 6; -5 -7 -7];
Bd = [5 1; 1 -4; -7 -2];
Qf = [3 0 0; 0 3 0; 0 0 9];
Qd = [5 0 0; 0 2 0; 0 0 8];
Rd = [3 0; 0 5];
n = 4;

N=1;
A=Ad; B=Bd; R=Rd; Q=Qd;

P = Qf;
F = zeros(size(Rd));
g = zeros(size(2*B'*Qf*A));
h = zeros(size(A'*P*A+Q));
K = zeros(size(g));

for i = 1:n-N
    F = B'*P*B+R;
    g = 2*B'*P*A;
    h = A'*P*A+Q;
    K = 2*F\g;
    P = h+g'/2*(-K)+(-K')*g/2+K'*F*K;
end

disp(mat2str(P))
disp(mat2str(K))
end


function P5
F1 = [7.0 4.5 4.0 -0.5; 4.5 6.0 0.5 0.0; 4.0 0.5 7.0 4.5; -0.5 0.0 4.5 9.0];
g1 = [-2.5; -2.0; 4.5; -1.0];
h1 = [3.0];
F2 = [0.0 0.0 4.5 6.5; 0.0 0.0 7.0 9.0; 4.5 7.0 2.0 7.0; 6.5 9.0 7.0 2.0];
g2 = [5.0; -4.5; 3.5; -4.0];
h2 = [-4.5];

    Z1 = -1*(F1\eye(size(F1)))*g1;
    Z2 = -1*(F2\eye(size(F2)))*g2;
    x = Z1; F = F1; g = g1; h = h1;     % test minimizer
    % x = Z2; F = F2; g = g2; h = h2;   % comment out the other one and remove
                                        % the comment on this one to see which
                                        % value is smaller, thats the answer
    disp(mat2str(x))
    disp(mat2str((1/2)*x'*F*x+g'*x+h))      
end