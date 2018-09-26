function HW21
format short g
clc;
P211
%P212
%P213
%P214
end

function x = P211

syms L k_F k_M sig1 sig2 sig3 sig4 sig5

P11 = [0; -(L/4); 0];
F11 = [0; 0; k_F];
M11 = [0; 0; -k_M];

tau11 = skew(P11)*F11 + M11;

%% Point 2
P22 = [L; -(L/4); 0];
F22 = [0; 0; k_F];
M22 = [0; 0; -k_M];

tau22 = skew(P22)*F22 + M22;

%% Point 3
P33 = [0; (L-(L/4)); 0];
F33 = [0; 0; k_F];
M33 = [0; 0; -k_M];

tau33 = skew(P33)*F33 + M33;

%% Point 4
P44 = [-L; -(L/4); 0];
F44 = [0; 0; k_F];
M44 = [0; 0; k_M];

tau44 = skew(P44)*F44 + M44;

%%
tauRotor = [tau11 tau22 tau33 tau44];
fRotor = [F11 F22 F33 F44];

W = [tauRotor; fRotor(3,:)]

equationsToMatrix(W,[k_M k_F L])
%{
syms L k_F k_M sig1 sig2 sig3 sig4 sig5

P11 = [-(L/5); -(2*L/5); 0];
F11 = [0; 0; k_F];
M11 = [0; 0; -k_M];

tau11 = skew(P11)*F11 + M11;

%% Point 2
P22 = [-(L+(L/5)); -(2*L/5); 0];
F22 = [0; 0; k_F];
M22 = [0; 0; k_M];

tau22 = skew(P22)*F22 + M22;

%% Point 3
P33 = [(L/2)+((L/2)-(L/5)); -(2*L/5); 0];
F33 = [0; 0; k_F];
M33 = [0; 0; -k_M];

tau33 = skew(P33)*F33 + M33;

%% Point 4
P44 = [(L/2)+((L/2)-(L/5)); (L/2)+(L/2 - ((2*L)/5));0];
F44 = [0; 0; k_F];
M44 = [0; 0; k_M];

tau44 = skew(P44)*F44 + M44;

%% Point 5

P55 = [-(L/5);(L/2)+(L/2 - ((2*L)/5));0];
F55 = [0; 0; k_F];
M55 = [0; 0; -k_M];

tau55 = skew(P55)*F55 + M55;

%%
tauRotor = [tau11 tau22 tau33 tau44 tau55];
fRotor = [F11 F22 F33 F44 F55];

W = [tauRotor; fRotor(3,:)]
%}
 end

function x = P212
clc; clear;
sigmamax = 1000.00000000;
k_F = 0.00001000;
k_M = 0.00000010;
L = 0.80000000;
m = 0.40000000;
g = 9.81000000;
u_desired = [-0.42947206; 0.16817971; -0.01002652; 4.88679239];

P11 = [0; L; 0];
F11 = [0; 0; k_F];
M11 = [0; 0; k_M];

tau11 = skew(P11)*F11 + M11;

%% Point 2
P22 = [0; -L; 0];
F22 = [0; 0; k_F];
M22 = [0; 0; k_M];

tau22 = skew(P22)*F22 + M22;

%% Point 3
P33 = [-L; 0; 0];
F33 = [0; 0; k_F];
M33 = [0; 0; -k_M];

tau33 = skew(P33)*F33 + M33;

%% Point 4
P44 = [L; 0; 0];
F44 = [0; 0; k_F];
M44 = [0; 0; -k_M];

tau44 = skew(P44)*F44 + M44;

%%
tauRotor = [tau11 tau22 tau33 tau44];
fRotor = [F11 F22 F33 F44];

W = [tauRotor; fRotor(3,:)]

s_desired = (W\eye(4))*u_desired
s = zeros(4,1);
for i = 1:4
    if s_desired(i)< 0
        s(i) = 0;
    elseif s_desired(i)> sigmamax^2
        s(i) = sigmamax;
    else
        s(i) = s_desired(i)
    end    
end

mat2str(W*s)
end

function x = P213
sigmamax = 1000.00000000;
k_F = 0.00001000;
k_M = 0.00000010;
L = 0.80000000;
m = 2.40000000;
g = 9.81000000;

%% Upward
((sigmamax^2*4*k_F)/(m))-g

%% Downward
9.81
end

function x = P214
sigmamax = 1000.00000000;
k_F = 0.00001000;
k_M = 0.00000010;
L = 0.80000000;
m = 2.20000000;
g = 9.81000000;
theta1 = 0.73000000;
theta2 = -0.06000000;
theta3 = 0.31000000;
u = [-0.05330167; 0.76841014; -0.04818420; 4.95167374];

Rz = Rz(cos(theta1),sin(theta1));
Ry = Ry(cos(theta2),sin(theta2));
Rx = Rx(cos(theta3),sin(theta3));

R = Rz*Ry*Rx

fo = mat2str(-(R*[0;0;u(4)])+[0; 0; m*g])
tau = mat2str([u(1); u(2); u(3)])
end

%Make Skew Matrix
function S = skew(w)
S = [0 -w(3) w(2); w(3) 0 -w(1); -w(2) w(1) 0];
end

%% Rotation Matrix Maker
function X = Rx(c,s)
    X = [1 0 0; 0 c -s; 0 s c];
end

function Y = Ry(c,s)
    Y= [c 0 s; 0 1 0; -s 0 c];
end

function Z = Rz(c,s)
    Z = [c -s 0; s c 0; 0 0 1];
end