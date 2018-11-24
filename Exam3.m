function Exam3a
clc;clear;
%% Control Instance 4

P5
end

function P1
Ad = [9 4 4; 1 -8 8; 3 -8 7];
Bd = [3; 2; -3];
Qf = [3 0 0; 0 2 0; 0 0 9];
Qd = [9 0 0; 0 2 0; 0 0 6];
Rd = [6];
Piplus = [12 5 4; 5 8 10; 4 10 14];

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
o_desired = [0.23; 0.34; 0.07];

m = 0.75677000;
g = 9.81000000;
J = [0.00402000 0.00000000 0.00000000; 0.00000000 0.00472000 0.00000000; 0.00000000 0.00000000 0.00664000];
sigmamax = 1000.00000000;
kF = 0.00001000;
kM = 0.00000010;
L = 0.80000000;
dt = 0.10000000;

% Data in MATLAB format
Qd = diag([100, 100, 100, 1, 100, 100, 10, 1, 1, 100, 10, 10]);
Rd = diag([1, 100, 10, 10]);

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
 

syms o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3 u1 u2 u3 u4 real
xc = [o1 o2 o3 th1 th2 th3 v1 v2 v3 w1 w2 w3]';
uc = [u1 u2 u3 u4]';

Ac = double(subs(subs(jacobian(getRates(J,xc,uc,m,g),xc),xc,xe),uc,ue))
Bc = double(subs(subs(jacobian(getRates(J,xc,uc,m,g),uc),xc,xe),uc,ue))


[Ad, Bd] = ssdata(c2d(ss(Ac,Bc,[],[]),dt));
K = dlqr(Ad,Bd,Qd,Rd);
mat2str(K)

x0 = [0.30; 0.60; -1.00; 0.80; 0.70; 0.50; -0.50; 0.20; 0.30; -0.70; -0.90; 0.50];
t0 = 0.00;
t1 = 0.60;

tf = 0.4;
x = x0;
t = 0;


while t<tf
    u = ue - K*(x-xe);
    xdot=matlabFunction(dx);
    [t,x] = ode45(@(t,x) xdot(x(1),x(2),double(u)),[t t+dt],x);
    x = x(end,:)';
    t = t(end,:)';
end
x
end

function P4
Ad = [8 -7; 0 7];
Bd = [5; -5];
Qf = [2 0; 0 6];
Qd = [7 0; 0 2];
Rd = [4];
n = 4;


A=Ad; B=Bd; R=Rd; Q=Qd;

P = Qf;

for i = 1:n
    F = B'*P*B+R;
    g = 2*B'*P*A;
    h = A'*P*A+Q;
    K = 2*F\g;
    P = h+g'*(-K)+K'*F*K;
end

disp(mat2str(P))
disp(mat2str(K))
end


function P5
F1 = [8.0 3.5 6.0 2.5; 3.5 6.0 6.5 6.5; 6.0 6.5 2.0 2.0; 2.5 6.5 2.0 10.0];
g1 = [-4.0; -5.0; 1.0; -2.5];
h1 = [-4.0];
F2 = [10.0 5.0 5.0 7.0; 5.0 8.0 5.5 4.0; 5.0 5.5 6.0 5.5; 7.0 4.0 5.5 9.0];
g2 = [-2.5; 3.5; -4.0; 1.5];
h2 = [-4.0];

    Z1 = -1*(F1\eye(size(F1)))*g1
    -inv(F1')*g1
    Z2 = -1*(F2\eye(size(F2)))*g2;
    %x = Z1; F = F1; g = g1; h = h1;     % test minimizer
    x = Z2; F = F2; g = g2; h = h2;   % comment out the other one and remove
                                        % the comment on this one to see which
                                        % value is smaller, thats the answer
    disp(mat2str(x))
    disp(mat2str((1/2)*x'*F*x+g'*x+h))  
    
    
end

function p1again
A = [-1 7 1; -7 9 5; -4 5 -9];
B = [-3 -5; 6 -5; -1 -7];
Q = [6 0 0; 0 9 0; 0 0 8];
Qd = [1 0 0; 0 1 0; 0 0 3];
R = [3 0; 0 7];
Piplus = [20 7 9; 7 10 -1; 9 -1 8];

    F = B'*Piplus*B+R;
    g = 2*B'*Piplus*A;
    h = A'*Piplus*A+Qd;
    K = 2*F\g;
    P = h+g'*(-K)+K'*F*K;
    mat2str(K)
    mat2str(P)
end

function xdot = getRates(J,x,u,m,g)
q= x(1:3);
th= x(4:6);
v= x(7:9);
w = x(10:12);

R = GetR_ZYX(th);
f = [0;0;m*g]-R*[0;0;u(4)];
tau = [u(1); u(2); u(3)];

dq = v;
dth = GetThetaDot_ZYX(th, w);
dv = f/m;
dw = inv(J)*(tau - skew(w)*J*w);

xdot = [dq; dth; dv; dw];
end


function A = skew(a)

    A = [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    
end

function R = GetR_ZYX(theta)

    c1 = cos(theta(1));
    s1 = sin(theta(1));
    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    R = [c1*c2 c1*s2*s3-s1*c3 c1*s2*c3+s1*s3;
         s1*c2 s1*s2*s3+c1*c3 s1*s2*c3-c1*s3;
         -s2 c2*s3 c2*c3];

end

function thetadot = GetThetaDot_ZYX(theta, w)

    c2 = cos(theta(2));
    s2 = sin(theta(2));
    c3 = cos(theta(3));
    s3 = sin(theta(3));
    thetadot = (1/c2)*[0 s3 c3; 0 c2*c3 -c2*s3; c2 s2*s3 s2*c3]*w;
    
end